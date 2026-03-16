当然可以。下面是**与英文版内容一一对应、可直接并列放在仓库中的中文 `README_zh.md`**。
术语保持与脚本/方法一致，便于你后续写文档或论文 Methods。

---

如果投递到slurm系统出错：

`sbatch --wrap="which R; R -q -e 'print(\"hi\")'; which Rscript; Rscript --version; Rscript /pmaster/xutingfeng/test.R"`

检查

请激活env:`conda activate /pmaster/xutingfeng/conda/envs/gwaslab`

# PRS 分析流程（基于 LDpred2）

本仓库提供一套 **基于 LDpred2 的稳健 PRS 分析流程**，包含权重生成、PRS 计算、tuning set 评估以及最优模型权重提取。

该流程使用 **Bash 脚本实现**，并具备 **逐步缓存（step-level caching）机制**：
**如果某一步的输出文件已存在且非空，该步骤将自动跳过**，非常适合长时间计算或中断后重跑。

---

## 1. 参考基因组版本

* **LD 参考（LDCache）**：**GRCh38**, snpID是chr:pos:ref:alt
* **GWAS 汇总统计（summary statistics）**：**必须为 GRCh38**, snpID是chr:pos:ref:alt
* **目标基因型数据（plink2 `--pfile`）**, ID是rsID，因为中间会自动做rsID转换；如果不是rsID的话，也可以全部chr:pos:ref:alt

> ⚠️ 若基因组版本不一致，将导致 LD 匹配错误，PRS 结果无效。

最后基因型计算是基于GRCh38得版本（取决于传入得geno）
---

## 2. GWAS 汇总统计输入格式

GWAS summary 文件需为 **制表符分隔**（`.tsv` 或 `.tsv.gz`），并包含以下**逻辑字段**：

| 字段名       | 含义            |
| --------- | ------------- |
| `chr`     | 染色体编号（1–22）   |
| `pos`     | GRCh38 上的碱基位置 |
| `rsid`    | 变异 ID（可为空或占位） |
| `a1`      | 效应等位基因        |
| `a0`      | 非效应等位基因       |
| `beta`    | 效应值           |
| `beta_se` | beta 的标准误     |
| `n_eff`   | 有效样本量         |

> 列顺序不做强制要求，但字段**语义必须正确**。

---

## 3. rsID 分配策略（HapMap3）

为了保证与 LDpred2 和 HapMap3 LD 参考的兼容性：

* 使用 **HapMap3 rsID 映射表（GRCh38）** 进行 rsID 分配
* 通过以下命令完成：

  ```bash
  python KeyMapReplacer.py -p EUR_38.rsid
  ```
* 该步骤会：

  * 基于 `(chr, pos)` 补全或替换 `rsid`
  * 将 SNP 限制在 HapMap3 集合内

📌 在运行 `plink2 --score` 之前，该步骤是 **必须的**。

---

## 4. 分析流程概览

整个 pipeline 包含以下 5 个主要步骤：

### 1️⃣ LDpred2 权重生成

* 通过 `LDpred2AndLassosum2LDCaches.R` 运行
* 输出：

  ```
  ldpred2_grid.tsv.gz
  ```

### 2️⃣ rsID 分配

* 基于 HapMap3 rsID 映射
* 输出：

  ```
  ldpred2_grid_rsid.tsv.gz
  ```

### 3️⃣ PRS 计算

* 使用 `plink2 --score`
* 并行计算多个 LDpred2 参数组合
* 标志输出文件：

  ```
  LDPred2_PRS.sscore
  ```

### 4️⃣ 模型评估

* 在 tuning set 中评估 PRS 表现
* 评估指标：**AUC**
* 输出：

  ```
  LDPred2_PRS.auc.tsv
  ```

### 5️⃣ 最优模型权重提取

* 自动选择 AUC 最优的 PRS 模型
* 提取对应的 LDpred2 SNP 权重
* 输出：

  ```
  ldpred2_grid.best.tsv.gz
  ```

---

## 5. 必须的命令行参数

以下 **三个参数是强制要求的**：

| 参数   | 说明                          |
| ---- | --------------------------- |
| `-g` | GWAS 汇总统计文件                 |
| `-f` | tuning set 表型文件（regenie 格式） |
| `-o` | 输出目录                        |

示例：

```bash
bash run_prs.sh \
  -g Input/MVP/GWAS/xxx.tsv.gz \
  -f Input/MVP/TuningSet/myocardial_infarction.regenie \
  -t 3 \
  -o Output/MVP/myocardial_infarction
```

---

## 6. 可选超参数（Hyperparameters）

| 参数                  | 默认值             | 说明                |
| ------------------- | --------------- | ----------------- |
| `-t`                | `3`             | tuning set 中表型所在列 |
| `-L`                | `LDCache/chr#.` | LD 参考路径           |
| `-p`                | UKB pfile 路径    | plink2 基因型前缀      |
| `--score-col-nums`  | `4-101`         | LDpred2 权重列范围     |
| `--score-col-start` | `5`             | AUC 计算起始列         |
| `--dry-run`         | 关闭              | 仅打印命令，不实际运行       |

---

## 7. 压缩策略

* 所有中间文件与最终结果均使用 **标准 `gzip`**
* 与以下工具完全兼容：

  * `zcat`
  * `plink2`
  * `csvkit`

> 本流程 **不需要** `bgzip` 或 `tabix` 索引。

---

## 8. 可复现性与断点续跑

* 每一步均会检查：

  * 输入文件是否存在且非空
  * 输出文件是否已存在
* 若输出存在 → **自动跳过**
* 适用于：

  * 计算中断后的重跑
  * 多参数实验
  * HPC 批量作业

---

## 9. 主要输出文件说明

| 文件                         | 含义                   |
| -------------------------- | -------------------- |
| `ldpred2_grid.best.tsv.gz` | 最优 LDpred2 模型 SNP 权重 |
| `LDPred2_PRS.sscore`       | 所有模型的 PRS 结果         |
| `LDPred2_PRS.auc.tsv`      | PRS 性能评估（AUC）        |
| `time.log`                 | 运行时间记录               |

---

## 10. 使用建议

* 强烈建议在分析前确认所有数据的 **基因组版本一致**
* 在正式使用最优模型前，检查 `LDPred2_PRS.auc.tsv`
* 大规模实验建议：

  * 先使用 `--dry-run` 验证流程
  * 结合 SLURM / LSF job array 使用

