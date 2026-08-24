#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
互含框架 · 网络拓扑可视化
plot_network.py — 生成互含网络图、虫洞路径高亮、拓扑分析

伸出去必然弯回来。
"""

import os
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from utils import (
    build_grid_3d, moore_neighbors, axial_neighbors,
    compute_R_mean, ensure_dir, log_print, timestamp_str
)


def generate_network(N, Nt, hops, w, lambda_coupling, k, seed, nbr_type='moore'):
    """
    生成互含网络（简化版，与主程序 g3duo_v3.py 兼容）

    返回:
        coords: (N³, 3) 节点坐标
        adj: (N³, N³) 邻接矩阵（稀疏约束 k）
        wormhole_pairs: list of (i, j) 虫洞连接对
    """
    np.random.seed(seed)
    N_nodes = N ** 3
    coords = build_grid_3d(N)

    # 初始化邻接矩阵
    adj = np.zeros((N_nodes, N_nodes), dtype=np.int8)
    wormhole_pairs = []

    # 本地邻居连接
    nbr_func = moore_neighbors if nbr_type == 'moore' else axial_neighbors

    for i in range(N_nodes):
        neighbors = nbr_func(i, N)
        # 稀疏约束：最多 k 条连接
        n_links = min(k, len(neighbors))
        if n_links > 0:
            chosen = np.random.choice(neighbors, size=n_links, replace=False)
            for j in chosen:
                adj[i, j] = 1
                adj[j, i] = 1

    # 虫洞连接（概率 w）
    n_wormholes = int(w * N_nodes)
    for _ in range(n_wormholes):
        a = np.random.randint(0, N_nodes)
        b = np.random.randint(0, N_nodes)
        if a != b and adj[a, b] == 0:
            adj[a, b] = 1
            adj[b, a] = 1
            wormhole_pairs.append((a, b))

    return coords, adj, wormhole_pairs


def plot_3d_network(coords, adj, wormhole_pairs=None, highlight_nodes=None,
                    title="互含网络拓扑", save_path=None, figsize=(12, 10)):
    """
    绘制 3D 网络拓扑图

    参数:
        coords: (N, 3) 坐标
        adj: (N, N) 邻接矩阵
        wormhole_pairs: 虫洞连接列表
        highlight_nodes: 高亮节点索引列表
        title: 图标题
        save_path: 保存路径
    """
    fig = plt.figure(figsize=figsize, dpi=150)
    ax = fig.add_subplot(111, projection='3d')

    N_nodes = len(coords)

    # 绘制节点
    node_colors = ['#74b9ff'] * N_nodes
    node_sizes = [20] * N_nodes

    if highlight_nodes:
        for idx in highlight_nodes:
            node_colors[idx] = '#ff6b6b'
            node_sizes[idx] = 80

    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
               c=node_colors, s=node_sizes, alpha=0.8, edgecolors='none')

    # 绘制本地连接（淡色）
    drawn = set()
    for i in range(N_nodes):
        for j in range(i + 1, N_nodes):
            if adj[i, j] == 1:
                # 判断是否为虫洞（距离远的连接）
                dist = np.linalg.norm(coords[i] - coords[j])
                is_wormhole = dist > 2.0  # 简单阈值

                if is_wormhole:
                    continue  # 虫洞单独画

                if (i, j) not in drawn:
                    ax.plot([coords[i, 0], coords[j, 0]],
                           [coords[i, 1], coords[j, 1]],
                           [coords[i, 2], coords[j, 2]],
                           color='#2a2d37', alpha=0.3, linewidth=0.5)
                    drawn.add((i, j))

    # 绘制虫洞（高亮）
    if wormhole_pairs:
        for a, b in wormhole_pairs[:50]:  # 最多画 50 条，避免太乱
            ax.plot([coords[a, 0], coords[b, 0]],
                   [coords[a, 1], coords[b, 1]],
                   [coords[a, 2], coords[b, 2]],
                   color='#f5a623', alpha=0.8, linewidth=1.5)

    # 设置
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title, fontsize=13, color='#e8e6e3')

    # 深色主题
    fig.patch.set_facecolor('#0a0c12')
    ax.set_facecolor('#0a0c12')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('#2a2d37')
    ax.yaxis.pane.set_edgecolor('#2a2d37')
    ax.zaxis.pane.set_edgecolor('#2a2d37')
    ax.tick_params(colors='#888')
    ax.xaxis.label.set_color('#888')
    ax.yaxis.label.set_color('#888')
    ax.zaxis.label.set_color('#888')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, facecolor='#0a0c12', edgecolor='none')
        log_print(f"网络图已保存: {save_path}")
    else:
        plt.show()

    plt.close()


def plot_degree_distribution(adj, save_path=None):
    """绘制度分布"""
    degrees = adj.sum(axis=1)

    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)

    bins = np.arange(degrees.min(), degrees.max() + 2) - 0.5
    ax.hist(degrees, bins=bins, color='#74b9ff', alpha=0.7, edgecolor='#0a0c12')
    ax.axvline(x=degrees.mean(), color='#f5a623', linestyle='--', 
               label=f'均值: {degrees.mean():.1f}')

    ax.set_xlabel('节点度', fontsize=12)
    ax.set_ylabel('节点数', fontsize=12)
    ax.set_title('互含网络 · 度分布', fontsize=13, color='#e8e6e3')
    ax.legend()
    ax.grid(True, alpha=0.2)

    fig.patch.set_facecolor('#0a0c12')
    ax.set_facecolor('#13161f')
    ax.tick_params(colors='#e8e6e3')
    ax.xaxis.label.set_color('#e8e6e3')
    ax.yaxis.label.set_color('#e8e6e3')
    ax.title.set_color('#e8e6e3')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, facecolor='#0a0c12', edgecolor='none')
        log_print(f"度分布图已保存: {save_path}")
    else:
        plt.show()

    plt.close()


def plot_adjacency_heatmap(adj, save_path=None):
    """绘制邻接矩阵热图"""
    fig, ax = plt.subplots(figsize=(8, 8), dpi=150)

    im = ax.imshow(adj, cmap='YlOrBr', interpolation='nearest', vmin=0, vmax=1)
    ax.set_xlabel('节点索引', fontsize=11)
    ax.set_ylabel('节点索引', fontsize=11)
    ax.set_title('互含网络 · 邻接矩阵', fontsize=13, color='#e8e6e3')
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='连接')

    fig.patch.set_facecolor('#0a0c12')
    ax.set_facecolor('#0a0c12')
    ax.tick_params(colors='#888')
    ax.xaxis.label.set_color('#888')
    ax.yaxis.label.set_color('#888')
    ax.title.set_color('#e8e6e3')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, facecolor='#0a0c12', edgecolor='none')
        log_print(f"邻接热图已保存: {save_path}")
    else:
        plt.show()

    plt.close()


def plot_wormhole_paths(coords, wormhole_pairs, save_path=None):
    """单独绘制虫洞路径图"""
    fig = plt.figure(figsize=(10, 10), dpi=150)
    ax = fig.add_subplot(111, projection='3d')

    # 只画虫洞
    for a, b in wormhole_pairs:
        ax.plot([coords[a, 0], coords[b, 0]],
               [coords[a, 1], coords[b, 1]],
               [coords[a, 2], coords[b, 2]],
               color='#f5a623', alpha=0.6, linewidth=1.0)

    # 节点
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
               c='#74b9ff', s=10, alpha=0.5)

    ax.set_title(f'虫洞路径图 (n={len(wormhole_pairs)})', fontsize=13, color='#e8e6e3')

    fig.patch.set_facecolor('#0a0c12')
    ax.set_facecolor('#0a0c12')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.tick_params(colors='#888')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, facecolor='#0a0c12', edgecolor='none')
        log_print(f"虫洞路径图已保存: {save_path}")
    else:
        plt.show()

    plt.close()


def main():
    parser = argparse.ArgumentParser(description='互含框架网络可视化工具')
    parser.add_argument('-N', type=int, default=8, help='网格尺寸 N (默认: 8)')
    parser.add_argument('--Nt', type=int, default=32, help='时间步 Nt (默认: 32)')
    parser.add_argument('--hops', type=int, default=16, help='最大 hops (默认: 16)')
    parser.add_argument('-w', type=float, default=0.02, help='虫洞概率 (默认: 0.02)')
    parser.add_argument('--lambda', dest='lambda_c', type=float, default=0.15, help='耦合强度 (默认: 0.15)')
    parser.add_argument('-k', type=int, default=6, help='稀疏约束 k (默认: 6)')
    parser.add_argument('--seed', type=int, default=42, help='随机种子 (默认: 42)')
    parser.add_argument('--nbr', default='moore', choices=['moore', 'axial'], help='邻居类型')
    parser.add_argument('-o', '--output', default='./plots', help='输出目录 (默认: ./plots)')
    parser.add_argument('--all', action='store_true', help='生成所有图')

    args = parser.parse_args()

    ensure_dir(args.output)

    log_print(f"生成网络: N={args.N}, k={args.k}, w={args.w}, seed={args.seed}")

    # 生成网络
    coords, adj, wormholes = generate_network(
        args.N, args.Nt, args.hops, args.w, args.lambda_c, args.k, args.seed, args.nbr
    )

    N_nodes = len(coords)
    R_mean = compute_R_mean(adj)

    log_print(f"节点数: {N_nodes}, 边数: {adj.sum()//2}, 虫洞数: {len(wormholes)}, R_mean: {R_mean:.4f}")

    ts = timestamp_str()
    prefix = f"net_N{args.N}_k{args.k}_w{args.w}_seed{args.seed}_{ts}"

    # 生成图
    if args.all:
        plot_3d_network(coords, adj, wormholes,
                       title=f"互含网络 3D | N={args.N}, k={args.k}, R={R_mean:.3f}",
                       save_path=os.path.join(args.output, f'{prefix}_3d.png'))

        plot_degree_distribution(adj,
                                save_path=os.path.join(args.output, f'{prefix}_degree.png'))

        if len(wormholes) > 0:
            plot_wormhole_paths(coords, wormholes,
                               save_path=os.path.join(args.output, f'{prefix}_wormholes.png'))

        # 小网络才画热图，太大内存不够
        if N_nodes <= 512:
            plot_adjacency_heatmap(adj,
                                  save_path=os.path.join(args.output, f'{prefix}_adj.png'))
    else:
        # 默认只画 3D 网络图
        plot_3d_network(coords, adj, wormholes,
                       title=f"互含网络 3D | N={args.N}, k={args.k}, R={R_mean:.3f}",
                       save_path=os.path.join(args.output, f'{prefix}_3d.png'))

    log_print("可视化完成")


if __name__ == '__main__':
    main()
