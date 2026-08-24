#!/usr/bin/env python3
"""
互含网络模拟 - 傻瓜式运行器
用法: python run.py          # 看任务列表
      python run.py WH-001   # 跑指定任务
"""

import sys
import json
import os
import re
import time
import shutil
from datetime import datetime

TASKS_FILE = "params/tasks.json"
RESULTS_DIR = "results"
SRC_FILE = "scripts/g3duo_v3_fast.py"

def load_tasks():
    with open(TASKS_FILE, 'r', encoding='utf-8') as f:
        return json.load(f)

def list_tasks():
    data = load_tasks()
    print("\n" + "="*60)
    print("📋 任务列表（选一个没有 ✅ 的跑）")
    print("="*60)
    for t in data["tasks"]:
        done = "✅" if t.get("done") else "🟢"
        mem = t['expected']['memory_gb']
        hrs = t['expected']['runtime_minutes'] // 60
        mins = t['expected']['runtime_minutes'] % 60
        time_str = f"{hrs}h{mins}m" if hrs else f"{mins}m"
        print(f"  {done} {t['task_id']:8} | N={t['params']['N']:2} Nt={t['params']['Nt']:2} "
              f"hops={t['params']['hops']:2} | {time_str} | {mem}GB | {t['description'][:20]}")
    print("="*60)
    print(f"\n👉 运行命令: python run.py <任务ID>")
    print(f"👉 示例:     python run.py WH-001\n")

def inject_params(task):
    """把参数注入到源码里，生成临时文件"""
    p = task['params']
    
    with open(SRC_FILE, 'r', encoding='utf-8') as f:
        code = f.read()
    
    # 替换配置区（只改大写开头的全局变量）
    def replace_var(code, name, value):
        # 匹配形如: N = 10 或 N=10 或 N = 10  # 注释
        pattern = rf'^(\s*{name}\s*=\s*)([^#\n]+)(.*)$'
        replacement = rf'\g<1>{value}\g<3>'
        code_new = re.sub(pattern, replacement, code, flags=re.MULTILINE)
        if code_new == code:
            print(f"⚠️  警告: 没找到变量 {name}，请检查 src/g3duo_v3.py 中是否有 {name} = xxx")
        return code_new
    
    code = replace_var(code, 'N', p['N'])
    code = replace_var(code, 'Nt', p['Nt'])
    code = replace_var(code, 'MAX_HOPS', p['hops'])
    code = replace_var(code, 'SEEDS', p['seeds'])
    
    # 如果有 MAX_HOPS 的别名 HOPS 也替换
    code = replace_var(code, 'HOPS', p['hops'])
    
    tmp_file = f"tmp_{task['task_id']}.py"
    with open(tmp_file, 'w', encoding='utf-8') as f:
        f.write(code)
    
    return tmp_file

def run_task(task_id):
    data = load_tasks()
    task = None
    for t in data["tasks"]:
        if t["task_id"] == task_id:
            task = t
            break
    
    if not task:
        print(f"❌ 找不到任务: {task_id}")
        list_tasks()
        return
    
    if task.get("done"):
        print(f"⚠️  任务 {task_id} 已经有结果了，确定要重跑吗？")
        ans = input("   输入 y 继续: ").strip().lower()
        if ans != 'y':
            return
    
    print(f"\n🚀 准备运行: {task_id}")
    print(f"   描述: {task['description']}")
    print(f"   参数: N={task['params']['N']}, Nt={task['params']['Nt']}, "
          f"k={task['params']['k']}, hops={task['params']['hops']}, seeds={task['params']['seeds']}")
    print(f"   预计: {task['expected']['runtime_minutes']}分钟, {task['expected']['memory_gb']}GB内存")
    print("-" * 50)
    
    # 检查内存
    try:
        import psutil
        avail_gb = psutil.virtual_memory().available / (1024**3)
        need_gb = task['expected']['memory_gb']
        if avail_gb < need_gb * 0.8:
            print(f"⚠️  警告: 可用内存 {avail_gb:.1f}GB，建议 {need_gb}GB，可能会卡")
            ans = input("   仍要继续吗？(y/n): ").strip().lower()
            if ans != 'y':
                return
    except ImportError:
        pass
    
    # 注入参数
    tmp_file = inject_params(task)
    print(f"✅ 已生成临时脚本: {tmp_file}")
    
    # 运行
    start = time.time()
    print(f"\n⏳ 开始计算...（别关窗口，预计 {task['expected']['runtime_minutes']} 分钟）\n")
    
    import subprocess
    try:
        proc = subprocess.Popen(
            [sys.executable, tmp_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )
        
        output_lines = []
        for line in proc.stdout:
            print(line, end='')
            output_lines.append(line)
        
        proc.wait()
        elapsed = time.time() - start
        output = ''.join(output_lines)
        
        if proc.returncode != 0:
            print(f"\n❌ 脚本异常退出，返回码: {proc.returncode}")
            return
            
    except KeyboardInterrupt:
        print("\n\n⛔ 用户中断")
        proc.kill()
        return
    finally:
        if os.path.exists(tmp_file):
            os.remove(tmp_file)
            print(f"\n🧹 已清理临时脚本")
    
    print(f"\n{'='*60}")
    print(f"✅ 计算完成！耗时: {elapsed/60:.1f} 分钟")
    print(f"{'='*60}")
    
    # 交互式收集结果
    print("\n📊 请从上方输出中复制关键数值（找不到就填 0）:")
    beta = input("   beta (β) 值: ").strip()
    r_mean = input("   R_mean 值: ").strip()
    candidate = input("   candidate_ratio (候选比例): ").strip() or "0"
    
    # 保存结果
    os.makedirs(RESULTS_DIR, exist_ok=True)
    result = {
        "task_id": task_id,
        "timestamp": datetime.now().isoformat(),
        "contributor": {
            "name": input("\n   你的昵称（用于致谢）: ").strip() or "anonymous",
            "github": input("   你的 GitHub ID（可选）: ").strip() or ""
        },
        "environment": {
            "python": sys.version.split()[0],
            "platform": sys.platform
        },
        "params": task['params'],
        "results": {
            "beta": float(beta) if beta else 0.0,
            "R_mean": float(r_mean) if r_mean else 0.0,
            "candidate_ratio": float(candidate),
            "runtime_seconds": round(elapsed, 1)
        },
        "raw_output_tail": output[-3000:]  # 保存最后3000字以便核查
    }
    
    filename = f"{task_id}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    filepath = os.path.join(RESULTS_DIR, filename)
    with open(filepath, 'w', encoding='utf-8') as f:
        json.dump(result, f, indent=2, ensure_ascii=False)
    
    print(f"\n🎉 结果已保存: {filepath}")
    print(f"\n👉 下一步: 把 {RESULTS_DIR}/ 目录下的这个文件")
    print(f"   上传到 GitHub（发 Pull Request）或发到群里")
    print(f"   文件名: {filename}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        list_tasks()
    else:
        run_task(sys.argv[1])
