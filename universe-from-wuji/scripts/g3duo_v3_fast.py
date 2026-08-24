#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
g3duo_v3_fast.py
虫洞模型v2 的 DIM 维通用版 · 性能优化版
- 用NumPy数组替代Python对象
- 预计算候选邻接表，O(N_space)替代O(N_space²)
- BFS绑定MAX_HOPS，用层内邻接表加速
"""
import sys
import numpy as np
import os
import time
import csv as csv_module
from collections import deque
import warnings
warnings.filterwarnings("ignore")

# 配置
DIM = 3
N = 8
Nt = 32
MAX_HOPS = 16
P_RETURN = 0.15
WORMHOLE_PROB = 0.02

if len(sys.argv) >= 2: DIM = int(sys.argv[1])
if len(sys.argv) >= 3: N = int(sys.argv[2])
if len(sys.argv) >= 4: Nt = int(sys.argv[3])
if len(sys.argv) >= 5: MAX_HOPS = int(sys.argv[4])
if len(sys.argv) >= 6: P_RETURN = float(sys.argv[5])
if len(sys.argv) >= 7: WORMHOLE_PROB = float(sys.argv[6])
if len(sys.argv) >= 8: SEED = int(sys.argv[7])
else: SEED = 42

out_name = f"v2_D{DIM}_N{N}_Nt{Nt}_h{MAX_HOPS}_w{WORMHOLE_PROB}"
save_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(save_dir)

print("="*60)
print(f"WORMHOLE v2 · {DIM}+1D · N={N}, Nt={Nt}, h={MAX_HOPS}, w={WORMHOLE_PROB}, seed={SEED}")
print("="*60)

np.random.seed(SEED)

N_space = N ** DIM

# ========== 1. 预计算空间坐标与邻居 ==========
# 用整数索引直接算邻居，避免unravel_index
def make_coords(N, DIM):
    # 提前生成所有坐标，用numpy
    coords = np.zeros((N**DIM, DIM), dtype=np.int16)
    for i in range(N**DIM):
        x = i
        for d in range(DIM-1, -1, -1):
            coords[i, d] = x % N
            x //= N
    return coords

coords = make_coords(N, DIM)

# 邻居偏移：每个轴±1
offsets = np.zeros((2*DIM, DIM), dtype=np.int8)
for d in range(DIM):
    offsets[2*d, d] = 1
    offsets[2*d+1, d] = -1

# 预计算每个点的邻居索引（周期边界）
neighbors = np.zeros((N_space, 2*DIM), dtype=np.int32)
for d in range(DIM):
    # 正方向
    idx = np.arange(N_space)
    c = coords[:, d]
    c_pos = (c + 1) % N
    # 计算新索引
    temp = coords.copy()
    temp[:, d] = c_pos
    # 转换回索引（按行优先）
    new_idx = np.ravel_multi_index(temp.T, [N]*DIM)
    neighbors[:, 2*d] = new_idx
    # 负方向
    c_neg = (c - 1) % N
    temp[:, d] = c_neg
    new_idx = np.ravel_multi_index(temp.T, [N]*DIM)
    neighbors[:, 2*d+1] = new_idx

print(f"空间点: {N_space}, 邻居/点: {2*DIM}")

# ========== 2. 构建网络（用邻接表CSR风格） ==========
# events[t][sid] = 事件ID
total_events = N_space * Nt
event_t = np.zeros(total_events, dtype=np.int16)
event_sid = np.zeros(total_events, dtype=np.int32)

# 邻接表：用list of lists（内存可控）
out_lists = [[] for _ in range(total_events)]
in_lists = [[] for _ in range(total_events)]
closed = np.zeros(total_events, dtype=bool)

t_start = time.time()

# Seed layer (t=0): 完全图
seed_start = 0
seed_end = N_space
for i in range(seed_start, seed_end):
    event_t[i] = 0
    event_sid[i] = i
    # 种子层完全连接
    all_ids = list(range(seed_start, seed_end))
    out_lists[i] = [j for j in all_ids if j != i]
    closed[i] = True
for i in range(seed_start, seed_end):
    all_ids = list(range(seed_start, seed_end))
    in_lists[i] = [j for j in all_ids if j != i]

# 生长循环
eid = N_space  # 下一个事件ID
# 预计算每层的局部邻域关系
layer_size = N_space

for t in range(1, Nt):
    layer_start = t * layer_size
    layer_end = (t + 1) * layer_size
    prev_start = (t - 1) * layer_size
    prev_end = t * layer_size

    # 为当前层每个空间点预计算候选
    for sid in range(N_space):
        new_id = layer_start + sid
        event_t[new_id] = t
        event_sid[new_id] = sid

        # 邻域：自身+空间邻居
        nb_sids = np.concatenate([[sid], neighbors[sid]])

        # 上一层的候选：space邻居对应的事件
        candidates = []
        for nb_sid in nb_sids:
            prev_id = prev_start + nb_sid
            candidates.append(prev_id)

        if not candidates:
            closed[new_id] = False
            continue

        # 采样出度
        n_out = np.random.poisson(0.6) + 1
        n_out = min(n_out, len(candidates))

        # 权重
        weights = np.ones(len(candidates))
        for idx_c, c_sid in enumerate(nb_sids):
            if c_sid == sid:
                weights[idx_c] = 3.0
            else:
                weights[idx_c] = 1.5
        weights /= weights.sum()

        chosen = np.random.choice(len(candidates), size=n_out, replace=False, p=weights)

        for idx_t in chosen:
            p_id = candidates[idx_t]
            # e_new -> p
            out_lists[new_id].append(p_id)
            in_lists[p_id].append(new_id)
            # p -> e_new（本地回返）
            out_lists[p_id].append(new_id)
            in_lists[new_id].append(p_id)
            closed[new_id] = True
            closed[p_id] = True

    # 虫洞层：非邻域低概率连接
    for sid in range(N_space):
        new_id = layer_start + sid
        nb_sids = set(list(nb_sids)) if 'nb_sids' in dir() else set([sid] + spatial_neighbors[sid])
        for prev_sid in range(N_space):
            if prev_sid not in nb_sids:
                if np.random.rand() < WORMHOLE_PROB:
                    p_id = prev_start + prev_sid
                    out_lists[p_id].append(new_id)
                    in_lists[new_id].append(p_id)
                    closed[new_id] = True
                    closed[p_id] = True

    if t % 5 == 0:
        print(f"  t={t}/{Nt}, events={layer_end}")

print(f"网络构建完成: {total_events} 事件, {time.time()-t_start:.1f}s")

# ========== 3. 计算R值（向量化） ==========
def compute_R_fast(t, sid, win_t=3):
    t0 = t
    sid0 = sid
    t_start = max(0, t0 - win_t)
    t_end = min(Nt, t0 + win_t + 1)
    local_sids = np.concatenate([[sid0], neighbors[sid0]])
    R_sum = 0.0
    local_events = []
    for tt in range(t_start, t_end):
        for ss in local_sids:
            ev_id = tt * N_space + ss
            if ev_id < total_events:
                local_events.append(ev_id)
    N_len = len(local_events)
    if N_len <= 1:
        return 0.0
    max_edges = N_len * (N_len - 1)
    actual = sum(len(out_lists[e]) for e in local_events)
    R = 1.0 - actual / max_edges
    return max(0.0, min(1.0, R))

R_by_t = np.zeros(Nt)
for t in range(Nt):
    R_sum = 0.0
    count = 0
    for sid in range(N_space):
        R_sum += compute_R_fast(t, sid)
        count += 1
    R_by_t[t] = R_sum / count

print(f"R计算完成: mean={R_by_t.mean():.4f}")
print(f"  Seed R={R_by_t[0]:.3f} -> Late R={R_by_t[-1]:.3f}")

# ========== 4. 采样BFS（绑定） ==========
def periodic_dist(sid1, sid2):
    c1 = coords[sid1]
    c2 = coords[sid2]
    dc = np.abs(c1 - c2)
    dc = np.minimum(dc, N - dc)
    return np.sqrt(np.sum(dc**2))

def bfs_bounded(start, target, max_d):
    if start == target:
        return 0
    vis = {start: 0}
    q = deque([start])
    while q:
        c = q.popleft()
        if vis[c] >= max_d:
            return -1
        for nid in out_lists[c]:
            if nid == target:
                return vis[c] + 1
            if nid not in vis:
                vis[nid] = vis[c] + 1
                q.append(nid)
    return -1

n_samples = 1000
pairs = []
for _ in range(n_samples):
    i, j = np.random.choice(total_events, 2, replace=False)
    dt = abs(event_t[i] - event_t[j])
    d_space = periodic_dist(event_sid[i], event_sid[j])
    d_coord = np.sqrt(dt**2 + d_space**2)
    d_net = bfs_bounded(i, j, MAX_HOPS)
    if d_net >= 0:
        pairs.append((d_coord, d_net))

pairs = np.array(pairs)
print(f"有效路径: {len(pairs)}/{n_samples}")

slope = 0.0
if len(pairs) > 30:
    dc = pairs[:, 0]
    dn = pairs[:, 1]
    mask = (dc > 0) & (dn > 0)
    if mask.sum() > 20:
        lc, ln = np.log(dc[mask]), np.log(dn[mask])
        slope = np.cov(lc, ln)[0, 1] / np.var(lc)

print(f"斜率β = {slope:.4f}")
print(f"判定: {'STABLE' if 0.80 < slope < 1.20 else 'UNSTABLE'}")

# ========== 5. 保存 ==========
suspended = int((~closed).sum())
eff_path_rate = len(pairs) / n_samples
dn_vals = pairs[:, 1] if len(pairs) > 0 else np.array([])
if len(dn_vals) > 1:
    hist, _ = np.histogram(dn_vals, bins=20)
    p = hist / hist.sum()
    p = p[p > 0]
    entropy = -np.sum(p * np.log2(p))
else:
    entropy = 0.0

# CSV
csv_path = f"scan_{out_name}_seed{SEED}.csv"
with open(csv_path, 'w', newline='') as f:
    w = csv_module.writer(f)
    w.writerow(['dim','N','N_t','h','w','Rm','beta','seed','eff_path_rate','entropy','suspended'])
    w.writerow([DIM, N, Nt, MAX_HOPS, WORMHOLE_PROB, R_by_t[4:].mean(), slope, SEED, eff_path_rate, entropy, suspended])

# 时间序列
ts_path = f"timeseries_{out_name}_seed{SEED}.csv"
with open(ts_path, 'w', newline='') as f:
    w = csv_module.writer(f)
    w.writerow(['t', 'R'])
    for t in range(Nt):
        w.writerow([t, R_by_t[t]])

print(f"[已保存] {csv_path}")
print(f"[已保存] {ts_path}")
print("="*60)
