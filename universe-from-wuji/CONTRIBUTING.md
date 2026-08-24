# 贡献指南 · 互含框架虫洞模型

感谢你愿意伸手。关系先于对象，数据先于完美。

---

## 📋 文件名规范

脚本会自动生成文件名，**请勿手动修改**，直接上传即可。

**格式**：
```
scan_v3_k{最大连接数}_3D_N{网格尺寸}_Nt{时间步}_h{hops}_w{虫洞概率}_nbr{邻居类型}_seed{种子}.csv
```

**示例**：
```
scan_v3_k26_3D_N8_Nt32_h16_w0.02_nbrmoore_seed42.csv
scan_v3_k6_3D_N8_Nt32_h16_w0.02_nbrmoore_seed42.csv
scan_v3_k3_3D_N8_Nt32_h16_w0.02_nbrmoore_seed42.csv
```

| 字段 | 含义 | 示例 |
|------|------|------|
| `k26` | 每个节点最大连接数（稀疏约束） | k3, k6, k26 |
| `3D` | 空间维度 | 3D |
| `N8` | 网格尺寸 N³ | N8, N12 |
| `Nt32` | 时间步数 | Nt32, Nt40 |
| `h16` | 最大 hops | h16, h20 |
| `w0.02` | 虫洞概率 | w0.02, w0.15 |
| `nbrmoore` | 邻居类型 | nbrmoore, nbraxial |
| `seed42` | 随机种子 | seed42, seed123 |

> 💡 命令行第 5 个参数 `0.15` 是 **lambda（耦合强度）**，目前不在文件名中。如需区分不同 lambda，请在提交说明里备注。

---

## 🚀 三种提交方式（任选其一）

### 方式一：Pull Request（推荐，有 Git 基础）

**第 1 步：Fork 仓库**
1. 打开 [https://github.com/luhuols/luhuols.github.io](https://github.com/luhuols/luhuols.github.io)
2. 点击右上角 **Fork** → 创建你的个人副本

**第 2 步：上传数据**
1. 进入你 Fork 后的仓库 → `universe-from-wuji/data/`
2. 点击 **Add file** → **Upload files**
3. 拖入你的 `scan_v3_*.csv` 文件
4. Commit message 写：
   ```
   [数据提交] N=8_Nt=32_h=16_w=0.02_k=6_seed=42_你的ID
   ```
5. 点击 **Commit changes**

**第 3 步：发 Pull Request**
1. 回到你 Fork 的仓库首页
2. 点击 **Contribute** → **Open pull request**
3. PR 标题：
   ```
   [数据提交] N=8_Nt=32_h=16_w=0.02_k=6_seed=42_你的ID
   ```
4. PR 正文模板：
   ```markdown
   - 网格尺寸 N³：8³
   - 时间步 Nt：32
   - 最大 hops：16
   - 虫洞概率 w：0.02
   - 稀疏约束 k：6
   - 邻居类型：moore
   - 扫描范围：0°–180°，步长 10°
   - 随机种子：42
   - 运行环境：本地 / Google Colab / 服务器
   - 耗时：约 X 分钟
   - 备注：（如有异常数据点请说明）
   ```
5. 点击 **Create pull request**
6. 等待维护者审核合并

---

### 方式二：GitHub Issue（不会 Git，有 GitHub 账号）

1. 打开 [Issues 页面](https://github.com/luhuols/luhuols.github.io/issues)
2. 点击 **New issue**
3. 标题格式：
   ```
   [数据提交] N=8_Nt=32_h=16_w=0.02_k=6_seed=42_你的ID
   ```
4. 正文贴运行参数（见上方模板）
5. 直接拖拽上传 CSV 文件（单文件 ≤ 25MB）
6. 点击 **Submit new issue**

---

### 方式三：邮件（最简，无 GitHub 账号）

- **收件人**：`13320205668@163.com`
- **邮件标题**：
  ```
  [互含框架角度扫描] N=8_Nt=32_h=16_w=0.02_k=6_seed=42_你的ID
  ```
- **正文**：贴运行参数和简要说明
- **附件**：`scan_v3_*.csv`

> 收到后维护者会手动上传并标注贡献来源。

---

## 🖥️ 运行命令参考

```bash
# 标准角度扫描（0°–180°，步长 10°，k=6）
python g3duo_v3_angle_scan_v3.py 3 8 32 16 0.15 0.02 42 moore 0 180 10 6

# 参数说明（共 12 个）：
#  3       = 维度（固定为3）
#  8       = 网格尺寸 N（N³ 节点）
#  32      = 时间步 Nt
#  16      = 最大 hops
#  0.15    = 耦合强度 lambda
#  0.02    = 虫洞概率 w
#  42      = 随机种子
#  moore   = 邻居类型（moore / axial）
#  0       = 起始角度（°）
#  180     = 结束角度（°）
#  10      = 角度步长（°）
#  6       = 稀疏约束 k（每个节点最大连接数）

# 不同 k 值对比示例：
python g3duo_v3_angle_scan_v3.py 3 8 32 16 0.15 0.02 42 moore 0 180 10 3   # 极稀疏
python g3duo_v3_angle_scan_v3.py 3 8 32 16 0.15 0.02 42 moore 0 180 10 26  # 较密
```

> 💻 **内存要求**：N=8 仅需约 8GB 内存，普通笔记本即可跑。N=12 建议 16GB 以上，N=16 建议 32GB。

---

## ✅ 数据验证清单

提交前请确认 CSV 包含以下列：

| 列名 | 含义 | 示例值 |
|------|------|--------|
| `angle` | 扫描角度（°） | 0, 10, 20... |
| `beta` | 该角度下度规指数 | 0.9798 |
| `R_mean` | 该角度下平均互含度 | 0.9318 |
| `stable` | 稳定性标志（1=稳定） | 1 |
| `path_rate` | 路径率 | 0.74 |
| `time_sec` | 该角度耗时（秒） | 704.1 |

- 文件大小参考：N=8, 步长 10° 约 **5–15 KB**
- 若超过 **50MB**，请压缩为 `.zip` 后提交
- 若出现 `stable=0` 或 `beta=0` 的异常行，请在提交说明中标注

---

## 🏆 贡献者墙

| 贡献者 | 数据文件 | 关键发现 | 状态 |
|--------|---------|---------|------|
| 刘巍 | `scan_v3_k26_3D_N8_Nt32_h16_w0.02_nbrmoore_seed42.csv` | 最优 30° β=0.9915，k 较密时峰值前移 | ✅ 已合并 |
| 刘巍 | `scan_v3_k6_3D_N8_Nt32_h16_w0.02_nbrmoore_seed42.csv` | 最优 130° β=0.9798，黄金角窗口激活 | ✅ 已合并 |
| 刘巍 | `scan_v3_k3_3D_N8_Nt32_h16_w0.02_nbrmoore_seed42.csv` | 最优 180° β=0.9921，极稀疏下边界峰值 | ✅ 已合并 |
| 待填写 | — | — | — |

> 你的名字将出现在这里。伸出去必然弯回来，数据是另一种伸手。

---

## 📬 联系

- 作者：刘巍
- 邮箱：13320205668@163.com
- 项目主页：[https://www.luhuo.online/universe-from-wuji/](https://www.luhuo.online/universe-from-wuji/)
