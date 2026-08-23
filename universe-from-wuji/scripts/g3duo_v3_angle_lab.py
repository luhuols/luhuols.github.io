#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
g3duo_v3_angle_lab.py
虫洞模型v2 · 角度约束实验版
支持命令行参数控制邻居类型、角度约束模式、容差等
用法:
  python g3duo_v3_angle_lab.py 3 8 32 16 0.15 0.02 42 axial none 10
  python g3duo_v3_angle_lab.py 3 8 32 16 0.15 0.02 42 moore golden 10
"""
import sys
import numpy as np
import os
import time
import csv as csv_module
from collections import deque
import warnings
warnings.filterwarnings("ignore")

# ==================== 配置（全部可由命令行覆盖）====================
DIM = 3
N = 8
Nt = 32
MAX_HOPS = 16
P_RETURN = 0.15
WORMHOLE_PROB = 0.02
SEED = 42

# 邻居类型: axial | moore
NEIGHBOR_TYPE = "axial"
# 角度约束: none | golden | reject_golden | only_90 | broad
ANGLE_MODE = "none"
# 角度容差（度）
ANGLE_TOL = 10.0

# 解析命令行
if len(sys.argv) >= 2: DIM = int(sys.argv[1])
if len(sys.argv) >= 3: N = int(sys.argv[2])
if len(sys.argv) >= 4: Nt = int(sys.argv[3])
if len(sys.argv) >= 5: MAX_HOPS = int(sys.argv[4])
if len(sys.argv) >= 6: P_RETURN = float(sys.argv[5])
if len(sys.argv) >= 7: WORMHOLE_PROB = float(sys.argv[6])
if len(sys.argv) >= 8: SEED = int(sys.argv[7])
if len(sys.argv) >= 9: NEIGHBOR_TYPE = sys.argv[8]
if len(sys.argv) >= 10: ANGLE_MODE = sys.argv[9]
if len(sys.argv) >= 11: ANGLE_TOL = float(sys.argv[10])

out_name = f"v2_D{DIM}_N{N}_Nt{Nt}_h{MAX_HOPS}_w{WORMHOLE_PROB}_nbr{NEIGHBOR_TYPE}_ang{ANGLE_MODE}_tol{ANGLE_TOL}"
save_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(save_dir)

print("="*70)
print(f"WORMHOLE v2 · ANGLE LAB · {DIM}+1D")
print(f"N={N}, Nt={Nt}, h={MAX_HOPS}, w={WORMHOLE_PROB}, seed={SEED}")
print(f"邻居={NEIGHBOR_TYPE}, 角度约束={ANGLE_MODE}, 容差={ANGLE_TOL}°")
print("="*70)

np.random.seed(SEED)
N_space = N ** DIM

# ========== 1. 坐标与邻居 ==========
def make_coords(N, DIM):
    coords = np.zeros((N**DIM, DIM), dtype=np.int16)
    for i in range(N**DIM):
        x = i
        for d in range(DIM-1, -1, -1):
            coords[i, d] = x % N
            x //= N
    return coords

coords = make_coords(N, DIM)

# 生成邻居偏移
if NEIGHBOR_TYPE == "moore" and DIM == 3:
    moore_offsets = []
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dz in [-1, 0, 1]:
                if dx == 0 and dy == 0 and dz == 0:
                    continue
                moore_offsets.append([dx, dy, dz])
    offsets = np.array(moore_offsets, dtype=np.int8)
    N_NEIGH = len(offsets)
    print(f"[Moore邻域] 邻居数={N_NEIGH}")
elif NEIGHBOR_TYPE == "broad" and DIM == 3:
    # 更大范围：2步Moore邻域（125-1=124邻居）
    broad_offsets = []
    for dx in [-2,-1,0,1,2]:
        for dy in [-2,-1,0,1,2]:
            for dz in [-2,-1,0,1,2]:
                if dx==0 and dy==0 and dz==0:
                    continue
                broad_offsets.append([dx,dy,dz])
    offsets = np.array(broad_offsets, dtype=np.int8)
    N_NEIGH = len(offsets)
    print(f"[Broad邻域] 邻居数={N_NEIGH}")
else:
    offsets = np.zeros((2*DIM, DIM), dtype=np.int8)
    for d in range(DIM):
        offsets[2*d, d] = 1
        offsets[2*d+1, d] = -1
    N_NEIGH = 2*DIM
    print(f"[Axial邻域] 邻居数={N_NEIGH}")

# 预计算邻居索引
neighbors = np.zeros((N_space, N_NEIGH), dtype=np.int32)
for k in range(N_NEIGH):
    off = offsets[k]
    temp = coords.copy()
    for d in range(DIM):
        temp[:, d] = (coords[:, d] + off[d]) % N
    new_idx = np.ravel_multi_index(temp.T, [N]*DIM)
    neighbors[:, k] = new_idx

# ========== 2. 角度约束函数 ==========
def shortest_dir_vec(sid1, sid2):
    c1 = coords[sid1].astype(np.float32)
    c2 = coords[sid2].astype(np.float32)
    dc = c2 - c1
    for d in range(DIM):
        if dc[d] > N // 2:
            dc[d] -= N
        elif dc[d] < -N // 2:
            dc[d] += N
    return dc

def angle_with_ref(sid1, sid2, ref_vec=None):
    """计算sid1->sid2的方向与参考向量的夹角（度）。ref_vec=None时用z轴。"""
    v = shortest_dir_vec(sid1, sid2)
    nv = np.linalg.norm(v)
    if nv < 1e-6:
        return 0.0
    if ref_vec is None:
        ref_vec = np.array([0.0, 0.0, 1.0], dtype=np.float32)
    nr = np.linalg.norm(ref_vec)
    if nr < 1e-6:
        return 0.0
    cosang = np.clip(np.dot(v, ref_vec)/(nv*nr), -1.0, 1.0)
    return np.degrees(np.arccos(cosang))

def check_angle_constraint(sid1, sid2, mode, tol):
    """返回True表示允许连接"""
    if mode == "none":
        return True
    # 计算方向与随机参考方向的夹角（为了各向同性，参考方向用sid1的坐标编码）
    # 更公平的做法：计算两个空间邻居之间的相对夹角，但这需要上下文
    # 简化：我们计算与"局部生长方向"的夹角——用sid1的坐标作为伪随机种子生成参考轴
    np.random.seed((sid1 * 73856093) % 2**31)
    ref = np.random.randn(3).astype(np.float32)
    ref /= np.linalg.norm(ref) + 1e-10
    theta = angle_with_ref(sid1, sid2, ref)

    if mode == "golden":
        target = 137.5
        return abs(theta - target) < tol
    elif mode == "reject_golden":
        target = 137.5
        return abs(theta - target) >= tol
    elif mode == "only_90":
        return abs(theta - 90.0) < tol
    elif mode == "broad_golden":
        # 允许黄金角±tol，以及其补角
        t1 = abs(theta - 137.5) < tol
        t2 = abs(theta - 42.5) < tol  # 180-137.5
        return t1 or t2
    return True

# ========== 3. 构建网络 ==========
total_events = N_space * Nt
event_t = np.zeros(total_events, dtype=np.int16)
event_sid = np.zeros(total_events, dtype=np.int32)
out_lists = [[] for _ in range(total_events)]
in_lists = [[] for _ in range(total_events)]
closed = np.zeros(total_events, dtype=bool)

t_start = time.time()

# Seed layer (t=0): 完全图
for i in range(N_space):
    event_t[i] = 0
    event_sid[i] = i
    out_lists[i] = [j for j in range(N_space) if j != i]
    closed[i] = True
for i in range(N_space):
    in_lists[i] = [j for j in range(N_space) if j != i]

# 生长循环
layer_size = N_space
for t in range(1, Nt):
    layer_start = t * layer_size
    prev_start = (t - 1) * layer_size

    for sid in range(N_space):
        new_id = layer_start + sid
        event_t[new_id] = t
        event_sid[new_id] = sid

        # 收集候选邻居（空间邻居+自身）
        nb_sids = np.concatenate([[sid], neighbors[sid]])

        # 【角度约束过滤】
        if ANGLE_MODE != "none":
            filtered = []
            for nb_sid in nb_sids:
                if check_angle_constraint(sid, nb_sid, ANGLE_MODE, ANGLE_TOL):
                    filtered.append(nb_sid)
            if len(filtered) == 0:
                # 如果全被过滤，保留自身
                filtered = [sid]
            nb_sids = np.array(filtered, dtype=np.int32)

        candidates = [prev_start + nb_sid for nb_sid in nb_sids]
        if not candidates:
            closed[new_id] = False
            continue

        n_out = np.random.poisson(0.6) + 1
        n_out = min(n_out, len(candidates))

        weights = np.ones(len(candidates))
        for idx_c, c_sid in enumerate(nb_sids):
            weights[idx_c] = 3.0 if c_sid == sid else 1.5
        weights /= weights.sum()

        chosen = np.random.choice(len(candidates), size=n_out, replace=False, p=weights)

        for idx_t in chosen:
            p_id = candidates[idx_t]
            out_lists[new_id].append(p_id)
            in_lists[p_id].append(new_id)
            out_lists[p_id].append(new_id)
            in_lists[new_id].append(p_id)
            closed[new_id] = True
            closed[p_id] = True

    # 虫洞层
    for sid in range(N_space):
        new_id = layer_start + sid
        nb_sids_set = set([sid] + neighbors[sid].tolist())
        for prev_sid in range(N_space):
            if prev_sid not in nb_sids_set:
                if np.random.rand() < WORMHOLE_PROB:
                    p_id = prev_start + prev_sid
                    out_lists[p_id].append(new_id)
                    in_lists[new_id].append(p_id)
                    closed[new_id] = True
                    closed[p_id] = True

    if t % 5 == 0:
        print(f"  t={t}/{Nt}, events={(t+1)*N_space}")

print(f"网络构建完成: {total_events} 事件, {time.time()-t_start:.1f}s")

# ========== 4. 计算R值 ==========
def compute_R_fast(t, sid, win_t=3):
    t0 = t
    sid0 = sid
    t_start_win = max(0, t0 - win_t)
    t_end = min(Nt, t0 + win_t + 1)
    local_sids = np.concatenate([[sid0], neighbors[sid0]])
    local_events = []
    for tt in range(t_start_win, t_end):
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
    for sid in range(N_space):
        R_sum += compute_R_fast(t, sid)
    R_by_t[t] = R_sum / N_space

print(f"R计算完成: mean={R_by_t.mean():.4f}, Seed={R_by_t[0]:.3f} -> Late={R_by_t[-1]:.3f}")

# ========== 5. 角度统计 ==========
angle_results = {}
if DIM == 3:
    vectors = []
    for eid in range(total_events):
        sid = event_sid[eid]
        for nid in out_lists[eid]:
            if event_t[eid] == event_t[nid]:
                sid_n = event_sid[nid]
                dc = shortest_dir_vec(sid, sid_n)
                if np.any(dc != 0):
                    vectors.append(dc)
    if len(vectors) > 0:
        vectors = np.array(vectors, dtype=np.float32)
        norms = np.linalg.norm(vectors, axis=1, keepdims=True)
        vectors_u = vectors / (norms + 1e-10)
        n_sample = min(5000, len(vectors_u))
        idx = np.random.choice(len(vectors_u), n_sample, replace=False)
        v_sample = vectors_u[idx]
        dot = np.clip(v_sample @ v_sample.T, -1.0, 1.0)
        triu = np.triu_indices(n_sample, k=1)
        angles_deg = np.degrees(np.arccos(np.clip(dot[triu], -1.0, 1.0)))

        special = {
            '45': ((40,50), 0),
            '54.7': ((50,58), 0),
            '90': ((85,95), 0),
            '109.5': ((105,114), 0),
            '137.5_golden': ((132.5,142.5), 0),
        }
        for name, ((lo,hi), _) in special.items():
            special[name] = ((lo,hi), ((angles_deg > lo) & (angles_deg < hi)).sum() / len(angles_deg))

        angle_results = {
            'mean': float(angles_deg.mean()),
            'std': float(angles_deg.std()),
            'n_vec': len(vectors),
            'special': special
        }
        print(f"角度统计: mean={angles_deg.mean():.1f}°, std={angles_deg.std():.1f}°")
        for name, (_, frac) in special.items():
            print(f"  {name}: {frac*100:.2f}%")

# ========== 6. BFS采样 ==========
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

print(f"β = {slope:.4f} → {'STABLE' if 0.80 < slope < 1.20 else 'UNSTABLE'}")

# ========== 7. 保存 ==========
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

csv_path = f"scan_{out_name}_seed{SEED}.csv"
with open(csv_path, 'w', newline='') as f:
    w = csv_module.writer(f)
    w.writerow(['dim','N','Nt','h','w','seed','nbr','ang_mode','ang_tol',
                'Rm','beta','eff_path_rate','entropy','suspended',
                'ang_mean','ang_std','ang_nvec'])
    row = [DIM, N, Nt, MAX_HOPS, WORMHOLE_PROB, SEED, NEIGHBOR_TYPE, ANGLE_MODE, ANGLE_TOL,
           R_by_t[4:].mean(), slope, eff_path_rate, entropy, suspended]
    if angle_results:
        row.extend([angle_results['mean'], angle_results['std'], angle_results['n_vec']])
    else:
        row.extend(['', '', ''])
    w.writerow(row)

ts_path = f"ts_{out_name}_seed{SEED}.csv"
with open(ts_path, 'w', newline='') as f:
    w = csv_module.writer(f)
    w.writerow(['t', 'R'])
    for t in range(Nt):
        w.writerow([t, R_by_t[t]])

print(f"[已保存] {csv_path}")
print(f"[已保存] {ts_path}")
print("="*70)
