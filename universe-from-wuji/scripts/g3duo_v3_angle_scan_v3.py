#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
g3duo_v3_angle_scan_v3.py
虫洞模型v2 · 角度扫描v3 · 「稀疏连接+链接间夹角」版
核心改动：
  - Moore 26方向为候选池
  - 强制每个节点只保留 k 条连接（k=3,6,12,26）
  - 在这 k 条中，优先选择链接间夹角 ≈ target 的组合
  - 分离「密度效应」和「角度效应"

用法:
  python g3duo_v3_angle_scan_v3.py 3 8 32 16 0.15 0.02 42 moore 0 180 10 3
  # 最后四个参数: 起始角, 终止角, 步长, k(连接数限制)
"""
import sys
import numpy as np
import os
import time
import csv as csv_module
from collections import deque
import warnings
warnings.filterwarnings("ignore")

# ==================== 配置 ====================
DIM = 3
N = 8
Nt = 32
MAX_HOPS = 16
P_RETURN = 0.15
WORMHOLE_PROB = 0.02
SEED = 42
NEIGHBOR_TYPE = "moore"
ANGLE_START = 0.0
ANGLE_END = 180.0
ANGLE_STEP = 10.0
K_MAX = 3  # 每个节点的最大连接数（v3核心参数）

if len(sys.argv) >= 2: DIM = int(sys.argv[1])
if len(sys.argv) >= 3: N = int(sys.argv[2])
if len(sys.argv) >= 4: Nt = int(sys.argv[3])
if len(sys.argv) >= 5: MAX_HOPS = int(sys.argv[4])
if len(sys.argv) >= 6: P_RETURN = float(sys.argv[5])
if len(sys.argv) >= 7: WORMHOLE_PROB = float(sys.argv[6])
if len(sys.argv) >= 8: SEED = int(sys.argv[7])
if len(sys.argv) >= 9: NEIGHBOR_TYPE = sys.argv[8]
if len(sys.argv) >= 10: ANGLE_START = float(sys.argv[9])
if len(sys.argv) >= 11: ANGLE_END = float(sys.argv[10])
if len(sys.argv) >= 12: ANGLE_STEP = float(sys.argv[11])
if len(sys.argv) >= 13: K_MAX = int(sys.argv[12])

save_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(save_dir)

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

# Moore邻域（26方向）
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
else:
    offsets = np.zeros((2*DIM, DIM), dtype=np.int8)
    for d in range(DIM):
        offsets[2*d, d] = 1
        offsets[2*d+1, d] = -1
    N_NEIGH = 2*DIM

neighbors = np.zeros((N_space, N_NEIGH), dtype=np.int32)
for k in range(N_NEIGH):
    off = offsets[k]
    temp = coords.copy()
    for d in range(DIM):
        temp[:, d] = (coords[:, d] + off[d]) % N
    new_idx = np.ravel_multi_index(temp.T, [N]*DIM)
    neighbors[:, k] = new_idx

# ========== 2. 方向向量工具 ==========
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

def angle_between(sid_center, sid_a, sid_b):
    da = shortest_dir_vec(sid_center, sid_a)
    db = shortest_dir_vec(sid_center, sid_b)
    na, nb = np.linalg.norm(da), np.linalg.norm(db)
    if na < 1e-6 or nb < 1e-6:
        return 0.0
    cosang = np.clip(np.dot(da, db) / (na * nb), -1.0, 1.0)
    return np.degrees(np.arccos(cosang))

# ========== 3. 【v3核心】稀疏+角度约束选k个连接 ==========
def select_k_neighbors(sid, nb_sids, target, tol, k_max):
    """
    从候选邻居中选出最多k_max个，优先保留链接间夹角≈target的组合。
    策略：
      1. 保留自身（sid）
      2. 从非自身候选中，找一对夹角最接近target的作为"锚点对"
      3. 继续添加与已有集合中链接夹角最接近target的候选
      4. 直到k_max或候选耗尽
    """
    non_self = [int(s) for s in nb_sids if int(s) != sid]
    if len(non_self) == 0:
        return np.array([sid], dtype=np.int32)

    selected = [sid]  # 自身必须保留

    # 如果没有角度约束（target<0），随机选k_max-1个
    if target < 0:
        rng = np.random.RandomState((sid * 73856093) % 2**31)
        n_pick = min(k_max - 1, len(non_self))
        chosen = rng.choice(non_self, size=n_pick, replace=False)
        selected.extend(chosen.tolist())
        return np.array(selected, dtype=np.int32)

    # 有角度约束：贪心选择
    # 先找一对夹角最接近target的
    best_pair = None
    best_err = 999.0
    for i in range(len(non_self)):
        for j in range(i+1, len(non_self)):
            ang = angle_between(sid, non_self[i], non_self[j])
            err = min(abs(ang - target), abs(ang - (180.0 - target)))
            if err < best_err:
                best_err = err
                best_pair = (non_self[i], non_self[j])

    if best_pair is None:
        # 只剩一个候选
        if len(non_self) > 0:
            selected.append(non_self[0])
        return np.array(selected, dtype=np.int32)

    selected.extend(list(best_pair))
    remaining = [s for s in non_self if s not in best_pair]

    # 贪心添加：每次选与已有集合中链接夹角最接近target的
    while len(selected) - 1 < k_max and len(remaining) > 0:
        best_cand = None
        best_err = 999.0
        for cand in remaining:
            # 计算cand与所有已选非自身节点的平均夹角误差
            errs = []
            for s in selected:
                if s == sid:
                    continue
                ang = angle_between(sid, cand, s)
                err = min(abs(ang - target), abs(ang - (180.0 - target)))
                errs.append(err)
            avg_err = np.mean(errs) if errs else 999.0
            if avg_err < best_err:
                best_err = avg_err
                best_cand = cand

        if best_cand is not None:
            selected.append(best_cand)
            remaining.remove(best_cand)
        else:
            break

    return np.array(selected, dtype=np.int32)

# ========== 4. 网络构建 ==========
def build_network(target_angle, tol, k_max):
    total_events = N_space * Nt
    event_t = np.zeros(total_events, dtype=np.int16)
    event_sid = np.zeros(total_events, dtype=np.int32)
    out_lists = [[] for _ in range(total_events)]
    in_lists = [[] for _ in range(total_events)]
    closed = np.zeros(total_events, dtype=bool)

    # Seed layer (t=0): 完全图
    for i in range(N_space):
        event_t[i] = 0
        event_sid[i] = i
        out_lists[i] = [j for j in range(N_space) if j != i]
        closed[i] = True
    for i in range(N_space):
        in_lists[i] = [j for j in range(N_space) if j != i]

    layer_size = N_space
    for t in range(1, Nt):
        layer_start = t * layer_size
        prev_start = (t - 1) * layer_size

        for sid in range(N_space):
            new_id = layer_start + sid
            event_t[new_id] = t
            event_sid[new_id] = sid

            # 候选邻居 = 自身 + 空间邻居
            nb_sids = np.concatenate([[sid], neighbors[sid]])

            # 【v3核心】稀疏+角度约束选k个
            nb_sids = select_k_neighbors(sid, nb_sids, target_angle, tol, k_max)

            candidates = [prev_start + int(nb_sid) for nb_sid in nb_sids]
            if not candidates:
                closed[new_id] = False
                continue

            n_out = np.random.poisson(0.6) + 1
            n_out = min(n_out, len(candidates))

            weights = np.ones(len(candidates))
            for idx_c, c_sid in enumerate(nb_sids):
                weights[idx_c] = 3.0 if int(c_sid) == sid else 1.5
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

    return out_lists, in_lists, closed, event_t, event_sid

# ========== 5. R值计算 ==========
def compute_R(t, sid, out_lists, neighbors, Nt, N_space):
    win_t = 3
    t0 = t
    sid0 = sid
    t_start_win = max(0, t0 - win_t)
    t_end = min(Nt, t0 + win_t + 1)
    local_sids = np.concatenate([[sid0], neighbors[sid0]])
    local_events = []
    for tt in range(t_start_win, t_end):
        for ss in local_sids:
            ev_id = tt * N_space + ss
            if ev_id < Nt * N_space:
                local_events.append(ev_id)
    N_len = len(local_events)
    if N_len <= 1:
        return 0.0
    max_edges = N_len * (N_len - 1)
    actual = sum(len(out_lists[e]) for e in local_events)
    R = 1.0 - actual / max_edges
    return max(0.0, min(1.0, R))

# ========== 6. BFS与Beta ==========
def periodic_dist(sid1, sid2):
    c1 = coords[sid1]
    c2 = coords[sid2]
    dc = np.abs(c1 - c2)
    dc = np.minimum(dc, N - dc)
    return np.sqrt(np.sum(dc**2))

def bfs_bounded(start, target, max_d, out_lists, event_t, event_sid):
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

def measure_beta(out_lists, in_lists, closed, event_t, event_sid, n_samples=800):
    total_events = Nt * N_space
    pairs = []
    for _ in range(n_samples):
        i, j = np.random.choice(total_events, 2, replace=False)
        dt = abs(event_t[i] - event_t[j])
        d_space = periodic_dist(event_sid[i], event_sid[j])
        d_coord = np.sqrt(dt**2 + d_space**2)
        d_net = bfs_bounded(i, j, MAX_HOPS, out_lists, event_t, event_sid)
        if d_net >= 0:
            pairs.append((d_coord, d_net))
    pairs = np.array(pairs)
    slope = 0.0
    eff_rate = len(pairs) / n_samples
    if len(pairs) > 30:
        dc = pairs[:, 0]
        dn = pairs[:, 1]
        mask = (dc > 0) & (dn > 0)
        if mask.sum() > 20:
            lc, ln = np.log(dc[mask]), np.log(dn[mask])
            slope = np.cov(lc, ln)[0, 1] / np.var(lc)
    return slope, eff_rate

# ========== 7. 主扫描循环 ==========
print("="*70)
print(f"ANGLE SCAN v3 · SPARSE+INTER-ANGLE · {DIM}+1D · N={N}, Nt={Nt}")
print(f"扫描范围: {ANGLE_START}° ~ {ANGLE_END}°, 步长={ANGLE_STEP}°")
print(f"邻居类型: {NEIGHBOR_TYPE}, 虫洞概率: {WORMHOLE_PROB}")
print(f"【v3核心】每个节点最多 k={K_MAX} 条连接（稀疏约束）")
print(f"角度约束 = 链接间夹角（对应SPUM球堆积角）")
print("="*70)

results = []
tol = ANGLE_STEP / 2.0

angles = np.arange(ANGLE_START, ANGLE_END + 0.1, ANGLE_STEP)
if 137.5 not in angles and ANGLE_START <= 137.5 <= ANGLE_END:
    angles = np.sort(np.append(angles, 137.5))
if 109.5 not in angles and ANGLE_START <= 109.5 <= ANGLE_END:
    angles = np.sort(np.append(angles, 109.5))
if 70.0 not in angles and ANGLE_START <= 70.0 <= ANGLE_END:
    angles = np.sort(np.append(angles, 70.0))

for target in angles:
    t0 = time.time()
    out_lists, in_lists, closed, evt_t, evt_sid = build_network(target, tol, K_MAX)

    R_sum = 0.0
    for sid in range(N_space):
        R_sum += compute_R(Nt//2, sid, out_lists, neighbors, Nt, N_space)
    R_mean = R_sum / N_space

    beta, eff_rate = measure_beta(out_lists, in_lists, closed, evt_t, evt_sid)
    stable = 1 if 0.80 < beta < 1.20 else 0
    dt = time.time() - t0

    print(f"θ={target:6.1f}° | β={beta:.4f} | R={R_mean:.4f} | 稳定={stable} | 路径率={eff_rate:.2f} | {dt:.1f}s")

    results.append({
        'target_angle': target,
        'beta': beta,
        'R_mean': R_mean,
        'stable': stable,
        'eff_rate': eff_rate,
        'time_sec': dt
    })

# ========== 8. 保存与分析 ==========
out_name = f"scan_v3_k{K_MAX}_{DIM}D_N{N}_Nt{Nt}_h{MAX_HOPS}_w{WORMHOLE_PROB}_nbr{NEIGHBOR_TYPE}"
csv_path = f"{out_name}_seed{SEED}.csv"

with open(csv_path, 'w', newline='') as f:
    w = csv_module.writer(f)
    w.writerow(['target_angle','beta','R_mean','stable','eff_rate','time_sec'])
    for r in results:
        w.writerow([r['target_angle'], r['beta'], r['R_mean'], r['stable'], r['eff_rate'], r['time_sec']])

# 自动分析
best = max(results, key=lambda x: x['beta'] if x['stable'] else -1)
print("\n" + "="*70)
print(f"扫描完成! 结果保存: {csv_path}")
print(f"最优角度: θ={best['target_angle']:.1f}° | β={best['beta']:.4f} | R={best['R_mean']:.4f}")

special = {
    '黄金角137.5°': 137.5,
    '四面体角109.5°': 109.5,
    'v1峰值70°': 70.0,
    '直角90°': 90.0,
}
print("\n特殊角度表现:")
for name, ang in special.items():
    closest = min(results, key=lambda x: abs(x['target_angle']-ang))
    print(f"  {name}: θ={closest['target_angle']:.1f}° β={closest['beta']:.4f} 稳定={closest['stable']}")

print("="*70)
