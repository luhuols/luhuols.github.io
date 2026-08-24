#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
互含框架 · β 指数批量分析与稳定性判据
analyze_beta.py — 对角度扫描结果做统计分析与报告生成

伸出去必然弯回来。
"""

import os
import sys
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from utils import (
    parse_filename, find_scan_csvs, angle_stats, compare_phases,
    stability_status, BETA_STABLE_THRESHOLD, BETA_EXCELLENT_THRESHOLD,
    SPECIAL_ANGLES, log_print, ensure_dir, timestamp_str
)


def analyze_single(csv_path, output_dir='./analysis'):
    """
    分析单个 CSV 文件

    生成:
        - 统计摘要 JSON
        - β-θ 曲线图
        - 稳定性报告
    """
    ensure_dir(output_dir)

    # 解析参数
    params = parse_filename(csv_path)
    if params is None:
        log_print(f"无法解析文件名: {csv_path}", 'ERROR')
        return None

    # 读取数据
    df = pd.read_csv(csv_path)
    if 'angle' not in df.columns or 'beta' not in df.columns:
        log_print(f"CSV 格式错误，缺少 angle/beta 列: {csv_path}", 'ERROR')
        return None

    # 统计分析
    stats = angle_stats(df)

    # 添加文件参数
    stats.update({
        'file': params['basename'],
        'k': params['k'],
        'N': params['N'],
        'Nt': params['Nt'],
        'h': params['h'],
        'w': params['w'],
        'nbr': params['nbr'],
        'seed': params['seed'],
        'tag': params['tag'],
    })

    # 稳定性评估
    if stats['stable_ratio'] >= 0.95:
        stats['overall_status'] = 'EXCELLENT'
    elif stats['stable_ratio'] >= 0.85:
        stats['overall_status'] = 'STABLE'
    else:
        stats['overall_status'] = 'UNSTABLE'

    # 保存 JSON
    json_name = params['basename'].replace('.csv', '_analysis.json')
    json_path = os.path.join(output_dir, json_name)
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(stats, f, indent=2, ensure_ascii=False, default=float)
    log_print(f"统计摘要已保存: {json_path}")

    # 绘图
    fig_name = params['basename'].replace('.csv', '_beta.png')
    fig_path = os.path.join(output_dir, fig_name)
    plot_beta_curve(df, stats, fig_path)
    log_print(f"β 曲线图已保存: {fig_path}")

    # 打印摘要
    print_summary(stats)

    return stats


def plot_beta_curve(df, stats, save_path):
    """绘制 β-θ 曲线"""
    fig, ax1 = plt.subplots(figsize=(10, 6), dpi=150)

    # 主坐标轴: beta
    color_beta = '#f5a623'
    ax1.plot(df['angle'], df['beta'], 'o-', color=color_beta, linewidth=2, markersize=5, label='β 指数')
    ax1.axhline(y=BETA_STABLE_THRESHOLD, color='#c0a060', linestyle='--', alpha=0.7, label=f'稳定阈值 β={BETA_STABLE_THRESHOLD}')
    ax1.axhline(y=BETA_EXCELLENT_THRESHOLD, color='#2ea44f', linestyle='--', alpha=0.7, label=f'优秀阈值 β={BETA_EXCELLENT_THRESHOLD}')

    # 标记最优角
    opt_angle = stats['optimal_angle']
    opt_beta = stats['optimal_beta']
    ax1.plot(opt_angle, opt_beta, 'r*', markersize=15, label=f'最优: θ={opt_angle}°, β={opt_beta:.4f}')

    # 标记特殊角
    special_colors = {'golden': '#ff6b6b', 'tetrahedral': '#4ecdc4', 'right': '#95e1d3'}
    for name, color in special_colors.items():
        key = f'{name}_angle'
        if key in stats:
            sa = stats[key]
            sb = stats.get(f'{name}_beta', 0)
            ax1.axvline(x=sa, color=color, linestyle=':', alpha=0.5)
            ax1.annotate(f'{name}\n{sa}°', xy=(sa, sb), fontsize=8, color=color,
                        ha='center', va='bottom')

    ax1.set_xlabel('扫描角度 θ (°)', fontsize=12)
    ax1.set_ylabel('度规指数 β', color=color_beta, fontsize=12)
    ax1.tick_params(axis='y', labelcolor=color_beta)
    ax1.set_ylim(0.7, 1.05)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='lower left', fontsize=9)

    # 次坐标轴: R_mean
    ax2 = ax1.twinx()
    color_R = '#74b9ff'
    ax2.plot(df['angle'], df['R_mean'], 's--', color=color_R, linewidth=1.5, markersize=4, alpha=0.7, label='R_mean')
    ax2.set_ylabel('平均互含度 R', color=color_R, fontsize=12)
    ax2.tick_params(axis='y', labelcolor=color_R)
    ax2.set_ylim(0.8, 1.0)

    # 标题
    k = stats.get('k', '?')
    N = stats.get('N', '?')
    title = f"互含框架 · β-θ 曲线 | k={k}, N={N}, seed={stats.get('seed', '?')}"
    plt.title(title, fontsize=13, color='#e8e6e3')

    fig.patch.set_facecolor('#0a0c12')
    ax1.set_facecolor('#13161f')
    ax1.tick_params(colors='#e8e6e3')
    ax1.xaxis.label.set_color('#e8e6e3')
    ax1.yaxis.label.set_color('#e8e6e3')
    ax2.tick_params(colors='#e8e6e3')
    ax2.yaxis.label.set_color('#e8e6e3')

    plt.tight_layout()
    plt.savefig(save_path, facecolor='#0a0c12', edgecolor='none')
    plt.close()


def print_summary(stats):
    """打印统计摘要"""
    print("\n" + "=" * 60)
    print("  互含框架 · β 指数分析报告")
    print("=" * 60)
    print(f"  文件: {stats['file']}")
    print(f"  参数: k={stats['k']}, N={stats['N']}, Nt={stats['Nt']}, seed={stats['seed']}")
    print("-" * 60)
    print(f"  扫描角度数: {stats['n_angles']}")
    print(f"  β 均值:     {stats['beta_mean']:.4f} ± {stats['beta_std']:.4f}")
    print(f"  β 范围:     [{stats['beta_min']:.4f}, {stats['beta_max']:.4f}]")
    print(f"  R_mean:     {stats['R_mean_avg']:.4f} ± {stats['R_mean_std']:.4f}")
    print(f"  稳定比例:   {stats['stable_ratio']*100:.1f}%")
    print(f"  综合状态:   {stats['overall_status']}")
    print("-" * 60)
    print(f"  最优角度:   θ = {stats['optimal_angle']:.1f}°")
    print(f"  最优 β:     {stats['optimal_beta']:.4f}")
    print("-" * 60)

    # 特殊角
    special_names = {
        'golden': '黄金角',
        'tetrahedral': '四面体角',
        'right': '直角',
    }
    for key, cn in special_names.items():
        a = stats.get(f'{key}_angle')
        b = stats.get(f'{key}_beta')
        s = stats.get(f'{key}_stable')
        if a is not None:
            status = '稳定' if s == 1 else '不稳定'
            print(f"  {cn}: θ={a:.1f}°, β={b:.4f} [{status}]")

    print("=" * 60 + "\n")


def batch_compare(csv_dir='./data', output_dir='./analysis'):
    """
    批量对比目录下所有 CSV

    生成对比表和综合图
    """
    ensure_dir(output_dir)

    csv_files = find_scan_csvs(csv_dir)
    if not csv_files:
        log_print(f"未找到 CSV 文件: {csv_dir}", 'ERROR')
        return

    log_print(f"找到 {len(csv_files)} 个数据文件")

    # 逐个分析
    all_stats = []
    dfs = {}
    for fpath in csv_files:
        params = parse_filename(fpath)
        if params is None:
            continue
        stats = analyze_single(fpath, output_dir)
        if stats:
            all_stats.append(stats)
            dfs[params['k']] = pd.read_csv(fpath)

    if len(all_stats) < 2:
        log_print("数据文件不足 2 个，跳过对比分析")
        return

    # 生成对比表
    compare_df = compare_phases(dfs)
    csv_out = os.path.join(output_dir, f'phase_comparison_{timestamp_str()}.csv')
    compare_df.to_csv(csv_out, index=False)
    log_print(f"对比表已保存: {csv_out}")

    # 打印对比表
    print("\n" + "=" * 80)
    print("  三相/多相对比")
    print("=" * 80)
    print(compare_df.to_string(index=False))
    print("=" * 80 + "\n")

    # 绘制对比图
    fig_path = os.path.join(output_dir, f'phase_comparison_{timestamp_str()}.png')
    plot_phase_comparison(dfs, fig_path)
    log_print(f"对比图已保存: {fig_path}")


def plot_phase_comparison(dfs_dict, save_path):
    """绘制多相 β 曲线对比图"""
    fig, ax = plt.subplots(figsize=(12, 7), dpi=150)

    colors = ['#ff6b6b', '#4ecdc4', '#f5a623', '#95e1d3', '#c7ceea']

    for i, (k, df) in enumerate(sorted(dfs_dict.items())):
        color = colors[i % len(colors)]
        label = f'k={k} (N={df.get("N", "?")})'
        ax.plot(df['angle'], df['beta'], 'o-', color=color, linewidth=2, 
                markersize=4, label=label, alpha=0.85)

    ax.axhline(y=BETA_STABLE_THRESHOLD, color='#c0a060', linestyle='--', alpha=0.5)
    ax.axhline(y=BETA_EXCELLENT_THRESHOLD, color='#2ea44f', linestyle='--', alpha=0.5)

    # 标记特殊角
    for name, angle in [('四面体角', 109.5), ('黄金角', 137.5)]:
        ax.axvline(x=angle, color='white', linestyle=':', alpha=0.3)
        ax.text(angle, 0.72, name, rotation=90, color='white', fontsize=9, ha='center')

    ax.set_xlabel('扫描角度 θ (°)', fontsize=12)
    ax.set_ylabel('度规指数 β', fontsize=12)
    ax.set_title('互含框架 · 多相 β-θ 对比', fontsize=14, color='#e8e6e3')
    ax.set_ylim(0.7, 1.05)
    ax.grid(True, alpha=0.2)
    ax.legend(loc='best', fontsize=10)

    fig.patch.set_facecolor('#0a0c12')
    ax.set_facecolor('#13161f')
    ax.tick_params(colors='#e8e6e3')
    ax.xaxis.label.set_color('#e8e6e3')
    ax.yaxis.label.set_color('#e8e6e3')
    ax.title.set_color('#e8e6e3')

    plt.tight_layout()
    plt.savefig(save_path, facecolor='#0a0c12', edgecolor='none')
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='互含框架 β 指数分析工具')
    parser.add_argument('input', nargs='?', default='./data', 
                        help='输入 CSV 文件或目录 (默认: ./data)')
    parser.add_argument('-o', '--output', default='./analysis', 
                        help='输出目录 (默认: ./analysis)')
    parser.add_argument('--batch', action='store_true', 
                        help='批量对比模式')

    args = parser.parse_args()

    if os.path.isfile(args.input):
        # 单文件分析
        analyze_single(args.input, args.output)
    elif os.path.isdir(args.input):
        if args.batch:
            batch_compare(args.input, args.output)
        else:
            # 分析目录下所有文件
            csv_files = find_scan_csvs(args.input)
            for f in csv_files:
                analyze_single(f, args.output)
    else:
        log_print(f"输入路径不存在: {args.input}", 'ERROR')
        sys.exit(1)


if __name__ == '__main__':
    main()
