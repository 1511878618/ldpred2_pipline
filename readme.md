下面是一份**更新后的、完整且可直接放到仓库里的 `README.md`**，内容已经与你现在的脚本实现（CLI 参数、gzip、skip 机制、LDpred2 → PRS → AUC → best weights）完全对齐，并在原 README 基础上做了结构化和补充说明。

---

# PRS Pipeline (LDpred2-based)

This repository provides a **robust, restartable PRS pipeline** based on **LDpred2**, including weight generation, PRS calculation, tuning-set evaluation, and best-model extraction.

The pipeline is implemented as a **Bash script with step-level caching**:
if an output file already exists and is non-empty, the corresponding step will be **automatically skipped**.

---

## 1. Reference & Genome Build

* **LD reference (LDCache)**: **GRCh38**, snpid assume to be chr:pos:ref:alt
* **GWAS summary statistics**: **must be GRCh38**, snpid assume to be chr:pos:ref:alt
* **Target genotype (plink2 `--pfile`)**: assumed to be rsID

> ⚠️ Genome build mismatch will lead to incorrect LD matching and invalid PRS results.

---

## 2. Input GWAS Summary Statistics

The GWAS summary file must be **tab-delimited** (`.tsv` or `.tsv.gz`) and contain the following **logical columns**:

| Column name | Description                                |
| ----------- | ------------------------------------------ |
| `chr`       | Chromosome (1–22)                          |
| `pos`       | Base-pair position (GRCh38)                |
| `rsid`      | Variant ID (may be missing or placeholder) |
| `a1`        | Effect allele                              |
| `a0`        | Non-effect allele                          |
| `beta`      | Effect size                                |
| `beta_se`   | Standard error of beta                     |
| `n_eff`     | Effective sample size                      |

> Column order is flexible, but column **semantics must match**.

---

## 3. rsID Assignment Strategy

To ensure compatibility with downstream tools and HapMap3-based LD references:

* rsID assignment is performed using **HapMap3 rsID mapping**
* Implemented via:

  ```bash
  python KeyMapReplacer.py -p EUR_38.rsid
  ```
* This step:

  * Replaces or fills `rsid` using `(chr, pos)`
  * Restricts variants to HapMap3-compatible SNPs

📌 This step is **mandatory** before `plink2 --score`.

---

## 4. Pipeline Overview

The pipeline consists of the following steps:

1. **LDpred2 weight generation**

   * Run via `LDpred2AndLassosum2LDCaches.R`
   * Output:

     ```
     ldpred2_grid.tsv.gz
     ```

2. **rsID assignment**

   * HapMap3-based mapping
   * Output:

     ```
     ldpred2_grid_rsid.tsv.gz
     ```

3. **PRS calculation**

   * Using `plink2 --score`
   * Multiple LDpred2 models are evaluated in parallel
   * Output (marker file):

     ```
     LDPred2_PRS.sscore
     ```

4. **Model evaluation**

   * Evaluate PRS performance on tuning set
   * Metric: AUC
   * Output:

     ```
     LDPred2_PRS.auc.tsv
     ```

5. **Best model extraction**

   * Automatically selects the best-performing score
   * Extracts the corresponding LDpred2 weights
   * Output:

     ```
     ldpred2_grid.best.tsv.gz
     ```

---

## 5. Required Command-Line Arguments

The pipeline enforces **three mandatory parameters**:

| Parameter | Description                                |
| --------- | ------------------------------------------ |
| `-g`      | GWAS summary statistics file               |
| `-f`      | Tuning set phenotype file (regenie format) |
| `-o`      | Output directory                           |

Example:

```bash
bash run_prs.sh \
  -g Input/MVP/GWAS/xxx.tsv.gz \
  -f Input/MVP/TuningSet/myocardial_infarction.regenie \
  -t 3 \
  -o Output/MVP/myocardial_infarction
```

---

## 6. Optional Hyperparameters

| Option              | Default          | Description                          |
| ------------------- | ---------------- | ------------------------------------ |
| `-t`                | `3`              | Phenotype column index in tuning set |
| `-L`                | `LDCache/chr#.`  | LDCache path pattern                 |
| `-p`                | (UKB pfile path) | plink2 genotype prefix               |
| `--score-col-nums`  | `4-101`          | LDpred2 score columns                |
| `--score-col-start` | `5`              | Score column start index for AUC     |
| `--dry-run`         | off              | Print commands without executing     |

---

## 7. Compression Policy

* All intermediate and final outputs use **standard `gzip`**
* Compatible with:

  * `zcat`
  * `plink2`
  * `csvkit`

> `bgzip` / `tabix` is **not required**, since no random-access indexing is used.

---

## 8. Reproducibility & Restartability

* Each step checks:

  * required input files exist and are non-empty
  * output file existence
* If output exists → **step is skipped**
* Safe for:

  * interrupted runs
  * parameter sweeps
  * HPC / batch execution

---

## 9. Expected Final Outputs

| File                       | Description                      |
| -------------------------- | -------------------------------- |
| `ldpred2_grid.best.tsv.gz` | Best LDpred2 SNP weights         |
| `LDPred2_PRS.sscore`       | PRS values for all tested models |
| `LDPred2_PRS.auc.tsv`      | AUC evaluation summary           |
| `time.log`                 | Runtime log                      |

---

## 10. Notes & Best Practices

* Always verify genome build consistency
* Inspect `LDPred2_PRS.auc.tsv` before trusting the selected model
* For large-scale experiments:

  * use `--dry-run` for validation
  * wrap the script in job arrays (SLURM / LSF)
