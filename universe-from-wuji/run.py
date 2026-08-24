#!/usr/bin/env python3
"""
互含网络模拟 - 傻瓜式运行器 v2
自动循环多种子，自动读取CSV结果，零手动输入
用法: python run.py          # 看任务列表
      python run.py WH-001   # 跑指定任务
"""

import sys
import json
import os
import re
import time
import glob
import csv
from datetime import datetime

TASKS_FILE = "params/tasks.json"
RESULTS_DIR = "results"
SRC_FILE = "scripts/g3duo_v3_fast.py"

def load_tasks():
    with open(TASKS_FILE, 'r', encoding='utf-8') as f:
        return json.load(f)

def list_tasks():
    data = load_tasks()
    print("\n" + "="*70)
    print("📋 任务列表（选一个没有 ✅ 的跑）")
    print("="*70)
    for t in data["tasks"]:
        done = "✅" if t.get("done") else "🟢"
        mem = t['expected']['memory_gb']
        n_seeds = t['params'].get('seeds', 1)
        mins = t['expected']['runtime_minutes'] * n_seeds
        hrs = mins // 60
        rem = mins % 60
        time_str = f"{hrs}h{rem}m" if hrs else f"{rem}m"
        print(f"  {done} {t['task_id']:8} | N={t['params']['N']:2} Nt={t['params']['Nt']:2} "
              f"hops={t['params']['hops']:2} ×{n_seeds}seeds | {time_str} | {mem}GB | {t['description'][:22]}")
    print("="*70)
    print(f"\n👉 运行命令: python run.py <任务ID>")
    print(f"👉 示例:     python run.py WH-001\n")

def inject_params(task, seed):
    """生成临时脚本，注入参数和种子"""
    p = task['params']
    
    with open(SRC_FILE, 'r', encoding='utf-8') as f:
        code = f.read()
    
    def replace_var(code, name, value):
        pattern = rf'^(\s*{name}\s*=\s*)([^#\n]+)(.*)$'
        replacement = rf'\g<1>{value}\g<3>'
        code_new = re.sub(pattern, replacement, code, flags=re.MULTILINE)
        if code_new == code:
            print(f"⚠️  警告: 没找到变量 {name}，请检查 {SRC_FILE}")
        return code_new
    
    code = replace_var(code, 'N', p['N'])
    code = replace_var(code, 'Nt', p['Nt'])
    code = replace_var(code, 'MAX_HOPS', p['hops'])
    code = replace_var(code, 'SEED', seed)
    
    # 禁用命令行参数覆盖，防止注入的值被sys.argv覆盖
    disable_argv = 'import sys\nsys.argv = [sys.argv[0]]\n'
    code = disable_argv + code
    
    tmp_file = f"tmp_{task['task_id']}_seed{seed}.py"
    with open(tmp_file, 'w', encoding='utf-8') as f:
        f.write(code)
    
    return tmp_file

def read_result_from_csv(task, seed):
    """从生成的CSV文件中读取结果"""
    p = task['params']
    pattern = f"scan_v2_D3_N{p['N']}_Nt{p['Nt']}_h{p['hops']}_w0.02_seed{seed}.csv"
    files = glob.glob(pattern)
    if not files:
        return None
    
    with open(files[0], 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        row = next(reader)
        return {
            'beta': float(row['beta']),
            'Rm': float(row['Rm']),
            'eff_path_rate': float(row['eff_path_rate']),
            'entropy': float(row['entropy']),
            'suspended': int(row['suspended'])
        }

def run_task(task_id):
    data = load_tasks()
    task = next((t for t in data['tasks'] if t['task_id'] == task_id), None)
    if not task:
        print(f"❌ 找不到任务: {task_id}")
        list_tasks()
        return
    
    n_seeds = task['params'].get('seeds', 1)
    single_mins = task['expected']['runtime_minutes']
    total_mins = single_mins * n_seeds
    
    print(f"\n🚀 任务: {task_id} | {task['description']}")
    print(f"   参数: N={task['params']['N']}, Nt={task['params']['Nt']}, hops={task['params']['hops']}")
    print(f"   种子数: {n_seeds} | 单个种子约{single_mins}分钟 | 预计总时间: {total_mins}分钟")
    print("-" * 60)
    
    # 内存检查
    try:
        import psutil
        avail_gb = psutil.virtual_memory().available / (1024**3)
        need_gb = task['expected']['memory_gb']
        if avail_gb < need_gb * 0.7:
            print(f"⚠️  警告: 可用内存仅 {avail_gb:.1f}GB，建议 {need_gb}GB")
            ans = input("   仍要继续吗？(y/n): ").strip().lower()
            if ans != 'y':
                return
    except ImportError:
        pass
    
    results = []
    for seed in range(n_seeds):
        print(f"\n🌱 种子 {seed+1}/{n_seeds} (SEED={seed})...")
        
        tmp_file = inject_params(task, seed)
        start = time.time()
        
        import subprocess
        try:
            proc = subprocess.run(
                [sys.executable, tmp_file],
                capture_output=True, text=True, timeout=single_mins*120
            )
        except subprocess.TimeoutExpired:
            print(f"   ⏱️  种子 {seed} 超时")
            if os.path.exists(tmp_file):
                os.remove(tmp_file)
            continue
        finally:
            if os.path.exists(tmp_file):
                os.remove(tmp_file)
        
        elapsed = time.time() - start
        
        if proc.returncode != 0:
            print(f"   ❌ 失败")
            if proc.stderr:
                print(f"   错误: {proc.stderr[-300:]}")
            continue
        
        # 自动读取CSV
        res = read_result_from_csv(task, seed)
        if res:
            res['seed'] = seed
            res['runtime_seconds'] = round(elapsed, 1)
            results.append(res)
            status = "STABLE" if 0.80 < res['beta'] < 1.20 else "UNSTABLE"
            print(f"   ✅ beta={res['beta']:.4f}, Rm={res['Rm']:.4f}, {status}, {elapsed/60:.1f}分钟")
        else:
            print(f"   ⚠️  没找到CSV结果文件")
    
    if not results:
        print("\n❌ 所有种子都失败了，请检查脚本或参数")
        return
    
    # 汇总统计
    betas = [r['beta'] for r in results]
    rms = [r['Rm'] for r in results]
    beta_mean = sum(betas) / len(betas)
    beta_std = (sum((b - beta_mean)**2 for b in betas) / len(betas))**0.5 if len(betas) > 1 else 0
    rm_mean = sum(rms) / len(rms)
    
    summary = {
        'task_id': task_id,
        'timestamp': datetime.now().isoformat(),
        'contributor': {
            'name': input("\n你的昵称（用于致谢）: ").strip() or 'anonymous',
            'github': input('GitHub ID（可选）: ').strip() or ''
        },
        'environment': {
            'python': sys.version.split()[0],
            'platform': sys.platform
        },
        'params': task['params'],
        'results': {
            'beta_mean': round(beta_mean, 4),
            'beta_std': round(beta_std, 4),
            'Rm_mean': round(rm_mean, 4),
            'n_success': len(results),
            'n_total': n_seeds,
            'total_runtime_minutes': round(sum(r['runtime_seconds'] for r in results) / 60, 1),
            'stability': 'STABLE' if 0.80 < beta_mean < 1.20 else 'UNSTABLE'
        },
        'per_seed': results
    }
    
    os.makedirs(RESULTS_DIR, exist_ok=True)
    filename = f"{task_id}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    filepath = os.path.join(RESULTS_DIR, filename)
    with open(filepath, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    
    print(f"\n{'='*60}")
    print(f"🎉 任务完成！结果已保存: {filepath}")
    print(f"   beta = {beta_mean:.4f} ± {beta_std:.4f}")
    print(f"   Rm   = {rm_mean:.4f}")
    print(f"   判定 = {summary['results']['stability']}")
    print(f"   成功种子: {len(results)}/{n_seeds}")
    print(f"{'='*60}")
    print(f"\n👉 下一步: 把 {RESULTS_DIR}/{filename} 提交到仓库或发到群里")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        list_tasks()
    else:
        run_task(sys.argv[1])


