# 互含框架 · 虫洞模型数值实验
**从无极涌现出三维宇宙：开源创世模拟器**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21817655.svg)](https://doi.org/10.5281/zenodo.21817655)

> **核心命题**：关系先于对象。三维时空不是被放进去的，是从互含网络的双向回返中**涌现**出来的。

## 什么是互含框架？

互含范畴（Inter-Containment Framework）是一套从纯粹关系网络出发重构物理学基础的理论尝试：

- **无极（Wuji）**：前宇宙基底，纯粹关系网络，具五性（对称/圆满/自洽/全息/无穷大）
- **弯（Wan）**：自发不对称，R 值度量对称破缺强度
- **双向回返**：伸出去必然弯回来，是度规涌现的必要条件
- **斩杀线**：R₀ = 1 − 1/e ≈ 0.632，文明最低生存线
- **三维锁定**：D=3 是唯一稳定的宏观维度

本仓库包含支撑论文 A/B/C 的全部数值实验脚本与原始数据。

## 仓库结构

```
universe-from-wuji/
├── scripts/          # Python 数值实验脚本
│   ├── g3duo_v3.py              # 主程序：度规涌现模拟
│   ├── g3duo_v3_angle_scan_v3.py   # 角度扫描：能量预算预测
│   ├── analyze_beta.py          # β 指数分析
│   ├── plot_network.py          # 网络可视化
│   └── utils.py                 # 工具函数
├── data/             # 实验输出数据（.csv）
├── docs/             # 实验备忘录与说明文档
├── README.md
├── ROADMAP.md
├── CONTRIBUTING.md
└── LICENSE
```

## 快速开始

### 环境要求
- Python ≥ 3.9
- NumPy, SciPy, Matplotlib, NetworkX

### 安装
```bash
git clone https://github.com/luhuols/luhuols.github.io.git
cd luhuols.github.io/universe-from-wuji
pip install numpy scipy matplotlib networkx
```

### 运行核心实验
```bash
cd scripts

# 标准三维度规涌现实验（N=8³, Nt=32, hops=16, k=6）
python g3duo_v3_angle_scan_v3.py 3 8 32 16 0.15 0.02 42 moore 0 180 10 6

# 精细扫描验证 109.5°（100°–140°，步长 2.5°）
python g3duo_v3_angle_scan_v3.py 3 8 32 16 0.15 0.02 42 moore 100 140 2.5 6

# N=16 大网格验证（需 32GB 内存）
python g3duo_v3_angle_scan_v3.py 3 16 40 10 0.15 0.02 42 moore 0 180 10 6
```

**参数顺序**：维度 N Nt hops lambda w seed nbr 起始角 结束角 步长 k

## 核心脚本说明

| 脚本 | 功能 | 对应论文 |
|------|------|----------|
| `g3duo_v3.py` | 互含网络生成、BFS 深度统计、β 指数计算 | A（度规涌现） |
| `g3duo_v3_angle_scan_v3.py` | 角度空间扫描、R_mean(θ) 与 β(θ) 映射 | C（斩杀线与能量预算） |
| `analyze_beta.py` | 批量 β 分析、稳定性判据、相图绘制 | A/C |
| `plot_network.py` | 网络拓扑可视化、虫洞路径高亮 | A |

## 关键结果速览

### 三相角度漂移（已验证）

| 相区 | 稀疏约束 k | 最优回返角 | 度规指数 β | 物理对应 |
|------|-----------|-----------|-----------|---------|
| 稀疏相 | k=3 | 180° | 0.9921 | 暗物质极限 |
| **临界相** | **k=6** | **109.5° / 130°** | **0.9786 / 0.9798** | **重子物质窗口** |
| 混沌相 | k=26 | 无显著峰值 | 0.93 ± 0.01 | 早期宇宙 |

- **维度锁定**：D=1,2,4,5 全 UNSTABLE，仅 D=3 STABLE
- **度规涌现**：N=5³–12³ 范围内 β = 0.85–1.11 ± 0.02
- **斩杀线**：R₀ = 1 − 1/e ≈ 0.632
- **角度预测**：109.5° 四面体角在 k=6 临界相中成为最优回返角

## 实验数据

三相角度扫描原始数据（CSV）已上传至 `data/` 目录：

- `scan_v3_k3_3D_N8_Nt32_h16_w0.02_nbrmoore_seed42.csv`
- `scan_v3_k6_3D_N8_Nt32_h16_w0.02_nbrmoore_seed42.csv`
- `scan_v3_k26_3D_N8_Nt32_h16_w0.02_nbrmoore_seed42.csv`

## 学术引用

本代码支撑以下论文（Zenodo DOI）：

- **论文A**：《基于分层互含因果网络的三维度规涌现》
- **论文C**：《互含度演化与宇宙学斩杀线》
- 完整 1–31 篇 DOI 索引见 [互含范畴论文列表](https://doi.org/10.5281/zenodo.21817655)

如果你使用了本代码，请引用：
```
刘巍. 互含框架虫洞模型数值实验 [Software]. Zenodo. https://doi.org/10.5281/zenodo.21817655
```

## 参与贡献

欢迎通过 Pull Request / Issue / 邮件提交角度扫描数据，共同验证宇宙能量预算预测：

**Ω_b : Ω_c : Ω_Λ ≈ 5% : 27% : 68%**

详见 [CONTRIBUTING.md](CONTRIBUTING.md)。

## 联系方式

- 作者：刘巍
- 邮箱：13320205668@163.com
- 项目主页：[https://www.luhuo.online/universe-from-wuji/](https://www.luhuo.online/universe-from-wuji/)

## 许可证

本项目采用 [MIT License](LICENSE) 开源。
