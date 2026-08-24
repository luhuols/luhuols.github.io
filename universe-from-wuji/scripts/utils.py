#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
互含框架 · 工具函数库
utils.py — 共享辅助函数

关系先于对象。
"""

import os
import re
import json
import numpy as np
from datetime import datetime


# ==================== 常量定义 ====================

# 斩杀线阈值
KILL_LINE = 1.0 - 1.0 / np.e  # ≈ 0.6321205588

# 稳定性阈值（度规涌现最低标准）
BETA_STABLE_THRESHOLD = 0.85
BETA_EXCELLENT_THRESHOLD = 0.95

# 标准几何角（度）
SPECIAL_ANGLES = {
    'tetrahedral': 109.4712,   # 四面体角 arccos(-1/3)
    'golden': 137.5078,         # 黄金角 360°/φ²
    'octahedral': 90.0,         # 八面体角
    'icosahedral': 63.4349,     # 二十面体角
    'right': 90.0,
    'straight': 180.0,
}

# 文件名解析正则
FILENAME_PATTERN = re.compile(
    r'scan_v3_k(?P<k>\d+)_3D_N(?P<N>\d+)_Nt(?P<Nt>\d+)_h(?P<h>\d+)_w(?P<w>[0-9.]+)_nbr(?P<nbr>\w+)_seed(?P<seed>\d+)(?:_(?P<tag>.+))?\.csv'
)


# ==================== 文件与路径 ====================

def ensure_dir(path):
    """确保目录存在"""
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def parse_filename(filepath):
    """
    从 CSV 文件名解析实验参数

    返回 dict: {k, N, Nt, h, w, nbr, seed, tag, basename, fullpath}
    """
    basename = os.path.basename(filepath)
    m = FILENAME_PATTERN.match(basename)
    if not m:
        return None

    return {
        'k': int(m.group('k')),
        'N': int(m.group('N')),
        'Nt': int(m.group('Nt')),
        'h': int(m.group('h')),
        'w': float(m.group('w')),
        'nbr': m.group('nbr'),
        'seed': int(m.group('seed')),
        'tag': m.group('tag') or '',
        'basename': basename,
        'fullpath': filepath,
    }


def find_scan_csvs(directory='.', pattern='scan_v3_*.csv'):
    """查找目录下所有扫描结果 CSV"""
    import glob
    files = glob.glob(os.path.join(directory, pattern))
    files.sort()
    return files


# ==================== 互含度与度规计算 ====================

def compute_R_mean(adj_matrix):
    """
    计算网络平均互含度 R_mean

    参数:
        adj_matrix: np.ndarray, 邻接矩阵 (N_nodes, N_nodes)

    返回:
        float: R_mean ∈ [0, 1]
    """
    N = adj_matrix.shape[0]
    if N == 0:
        return 0.0

    # 互含度 = 本地回返率 + 长程虫洞贡献的归一化度量
    # 简化实现：基于邻接矩阵的局部聚类系数加权
    degrees = adj_matrix.sum(axis=1)
    nonzero = degrees > 0
    if not nonzero.any():
        return 0.0

    # 局部互含度近似：邻居间的实际连接数 / 可能连接数
    local_R = np.zeros(N)
    for i in range(N):
        if degrees[i] < 2:
            local_R[i] = 1.0  # 度为0或1时，视为完全回返
            continue
        neighbors = np.where(adj_matrix[i])[0]
        if len(neighbors) < 2:
            local_R[i] = 1.0
            continue
        # 邻居子图的边数
        sub = adj_matrix[np.ix_(neighbors, neighbors)]
        actual = sub.sum() / 2  # 无向图，每条边算两次
        possible = len(neighbors) * (len(neighbors) - 1) / 2
        if possible > 0:
            local_R[i] = actual / possible
        else:
            local_R[i] = 1.0

    return float(local_R.mean())


def compute_beta_metric(distances, N_nodes, box_size):
    """
    计算度规指数 β

    参数:
        distances: np.ndarray, 网络最短距离矩阵
        N_nodes: int, 节点总数
        box_size: float, 盒子尺度（网络特征长度）

    返回:
        float: β 指数，~1.0 表示三维欧氏度规
    """
    # 提取有限距离（排除 inf 和 0）
    mask = (distances > 0) & np.isfinite(distances)
    if not mask.any():
        return 0.0

    d_vals = distances[mask]

    # 用盒计数法估算分形维数
    # 在三维网格中，距离 r 内的节点数应 ~ r³
    # 取中位数距离作为特征尺度
    r_med = np.median(d_vals)
    if r_med <= 0:
        r_med = 1.0

    # 统计 r_med 范围内的节点数比例
    count_within = (d_vals <= r_med).sum()
    total_pairs = len(d_vals)

    if total_pairs == 0:
        return 0.0

    # 理想三维：N(r) ~ r³，在 r = r_med 时，N/N_total ~ (r_med/L)³
    # 但实际网络有离散性，用对数斜率估算
    # 简化：用距离分布的均值/标准差来估算"紧致度"
    d_mean = d_vals.mean()
    d_std = d_vals.std()

    if d_mean <= 0:
        return 0.0

    # β 指数：衡量距离分布与理想三维网格的吻合度
    # β ≈ 1.0 表示完美三维欧氏度规
    # β < 0.85 表示不稳定（UNSTABLE）
    # β > 1.1 表示过紧致（可能维度塌陷）

    cv = d_std / d_mean  # 变异系数
    # 理想泊松分布 CV ≈ 1/sqrt(n)，三维网格 CV 较小
    # 用经验公式映射到 β
    beta = 1.0 / (1.0 + cv)

    # 修正：基于网络尺寸归一化
    expected_mean = box_size * (N_nodes ** (1/3)) / 3.0
    if expected_mean > 0:
        ratio = d_mean / expected_mean
        beta = beta * (2.0 - ratio) if ratio < 2.0 else beta * (2.0 / ratio)

    return float(np.clip(beta, 0.0, 2.0))


def stability_status(beta):
    """
    根据 β 指数判断稳定性

    返回: 'UNSTABLE' | 'STABLE' | 'EXCELLENT'
    """
    if beta < BETA_STABLE_THRESHOLD:
        return 'UNSTABLE'
    elif beta >= BETA_EXCELLENT_THRESHOLD:
        return 'EXCELLENT'
    else:
        return 'STABLE'


# ==================== 数据统计 ====================

def angle_stats(df_csv):
    """
    对角度扫描 CSV 做统计分析

    参数:
        df_csv: pandas DataFrame, 含 angle, beta, R_mean, stable 等列

    返回:
        dict: 统计摘要
    """
    import pandas as pd

    if isinstance(df_csv, str):
        df_csv = pd.read_csv(df_csv)

    result = {
        'n_angles': len(df_csv),
        'beta_mean': float(df_csv['beta'].mean()),
        'beta_std': float(df_csv['beta'].std()),
        'beta_min': float(df_csv['beta'].min()),
        'beta_max': float(df_csv['beta'].max()),
        'R_mean_avg': float(df_csv['R_mean'].mean()),
        'R_mean_std': float(df_csv['R_mean'].std()),
        'stable_ratio': float(df_csv['stable'].mean()),
        'optimal_angle': float(df_csv.loc[df_csv['beta'].idxmax(), 'angle']),
        'optimal_beta': float(df_csv['beta'].max()),
    }

    # 特殊角度表现
    for name, angle_val in SPECIAL_ANGLES.items():
        # 找最接近的角度
        idx = (df_csv['angle'] - angle_val).abs().idxmin()
        result[f'{name}_angle'] = float(df_csv.loc[idx, 'angle'])
        result[f'{name}_beta'] = float(df_csv.loc[idx, 'beta'])
        result[f'{name}_stable'] = int(df_csv.loc[idx, 'stable'])

    return result


def compare_phases(dfs_dict):
    """
    对比多相数据（k=3/6/26...）

    参数:
        dfs_dict: dict, {k_value: DataFrame}

    返回:
        pandas DataFrame: 对比表
    """
    import pandas as pd

    rows = []
    for k, df in sorted(dfs_dict.items()):
        stats = angle_stats(df)
        rows.append({
            'k': k,
            'optimal_angle': stats['optimal_angle'],
            'optimal_beta': stats['optimal_beta'],
            'beta_mean': stats['beta_mean'],
            'beta_std': stats['beta_std'],
            'R_mean': stats['R_mean_avg'],
            'stable_ratio': stats['stable_ratio'],
            'golden_beta': stats.get('golden_beta', np.nan),
            'tetrahedral_beta': stats.get('tetrahedral_beta', np.nan),
        })

    return pd.DataFrame(rows)


# ==================== 日志与时间 ====================

def timestamp_str():
    return datetime.now().strftime('%Y%m%d_%H%M%S')


def log_print(msg, level='INFO'):
    """带时间戳的打印"""
    ts = datetime.now().strftime('%H:%M:%S')
    print(f"[{ts}] [{level}] {msg}")


# ==================== 互含网络构建辅助 ====================

def build_grid_3d(N):
    """构建 N³ 三维网格坐标"""
    x = np.arange(N)
    y = np.arange(N)
    z = np.arange(N)
    coords = np.array(np.meshgrid(x, y, z, indexing='ij')).reshape(3, -1).T
    return coords  # shape: (N³, 3)


def grid_distance(a, b, N):
    """三维网格周期边界距离"""
    d = np.abs(a - b)
    d = np.minimum(d, N - d)  # 周期边界
    return np.sqrt((d ** 2).sum())


def moore_neighbors(idx, N):
    """获取 Moore 邻居索引（26邻域）"""
    x, y, z = idx % N, (idx // N) % N, idx // (N * N)
    neighbors = []
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dz in [-1, 0, 1]:
                if dx == 0 and dy == 0 and dz == 0:
                    continue
                nx, ny, nz = (x + dx) % N, (y + dy) % N, (z + dz) % N
                neighbors.append(nx + ny * N + nz * N * N)
    return neighbors


def axial_neighbors(idx, N):
    """获取轴向邻居索引（6邻域）"""
    x, y, z = idx % N, (idx // N) % N, idx // (N * N)
    neighbors = []
    for dx, dy, dz in [(-1,0,0),(1,0,0),(0,-1,0),(0,1,0),(0,0,-1),(0,0,1)]:
        nx, ny, nz = (x + dx) % N, (y + dy) % N, (z + dz) % N
        neighbors.append(nx + ny * N + nz * N * N)
    return neighbors


if __name__ == '__main__':
    print("互含框架 · 工具函数库")
    print(f"斩杀线: {KILL_LINE:.6f}")
    print(f"稳定性阈值: {BETA_STABLE_THRESHOLD}")
    print("关系先于对象。")
