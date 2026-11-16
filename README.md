# Testing-Process Profiles and Mixed-Effects Modeling

### A Robust Multivariate Analysis of COVID-19 Diagnostic Data

**Advanced Multivariate Statistics · Academic Year 2025–2026**  
**Corinne D'Elia and Martina Filice · Prof. Chiara Masci**

## Project Overview

This project investigates the `covid_testing` dataset from the R `medicaldata` package using robust multivariate methods, model-based clustering, and logistic mixed-effects regression. It identifies demographic and testing-process profiles, quantifies heterogeneity across clinics or specimen-collection units, and compares retrospective predictive performance on held-out patients.

Each observation represents a diagnostic test. Some patients contributed multiple tests, so the analysis uses a patient-aware split, a patient-level random-intercept sensitivity model, and a patient-level bootstrap. The turnaround-time variables describe stages of the diagnostic workflow and can become available after specimen collection or result verification, which defines the retrospective scope of the study.

## Research Questions

1. **RQ1 — Robust latent testing-process profiles.** What demographic and testing-process profiles emerge, and how do conventional Gaussian mixture modeling and robust trimmed clustering (TCLUST) compare?
2. **RQ2 — Logistic mixed-effects modeling.** How are age and TCLUST profiles associated with positivity, accounting for collection-unit heterogeneity and repeated tests from the same patient?
3. **RQ3 — Out-of-sample predictive validation.** How does a clinic-level logistic GLMM compare with standard logistic regression on held-out patients, and how uncertain are the differences?

## Dataset and Patient-Aware Split

The raw dataset contains **15,524 diagnostic tests**, **12,344 distinct patient identifiers**, and **17 variables**. After removing invalid outcomes and keeping records complete in the core analysis fields, **15,223 tests** remain. The binary response `result_bin` encodes positive tests as `1` and negative tests as `0`.

| Variable | Description | Role |
|---|---|---|
| `age` | Patient age | Clustering feature and fixed effect |
| `pan_day` | Day on the recorded pandemic timeline | Clustering feature |
| `col_rec_tat` | Collection-to-receipt turnaround time in hours | Clustering feature |
| `rec_ver_tat` | Receipt-to-verification turnaround time in hours | Clustering feature |
| `clinic_name` | Clinic or specimen-collection unit | Clinic random intercept |
| `subject_id` | Patient identifier | Split, patient random intercept, and bootstrap |
| `result_bin` | Recorded binary test outcome | Response in the logistic models |

The analysis excludes the PCR cycle threshold `ct_result` to prevent direct target leakage. It also drops `payor_group` and `patient_class` because of extensive missingness.

The notebook randomly splits unique `subject_id` values 80/20 with `random_state=42` and keeps every patient's tests in the same partition. The split uses patient identifiers without outcome stratification.

| Partition | Tests |
|---|---:|
| Training | 12,196 |
| Held-out test | 3,027 |
| Total | 15,223 |

The held-out set contains **187 positive** and **2,840 negative** tests, for a positive class prevalence of **6.18%**.

## Methodology

### 1. Exploratory Analysis and Preprocessing

The notebook examines outcome imbalance, missingness, clinic sizes, repeated tests, and the distributions of the continuous features. It assesses redundancy using Spearman correlations and variance inflation factors (VIF).

The notebook computes robust squared Mahalanobis distances using the Minimum Covariance Determinant (MCD) estimator. The 99th percentile of a chi-square distribution with four degrees of freedom serves as an exploratory cutoff. The notebook retains flagged records for subsequent analysis. It fits Yeo–Johnson transformations and standardization on the training features and applies the fitted transformations to the held-out set.

### 2. Gaussian Mixtures and TCLUST (RQ1)

The R package `mclust` compares Gaussian mixture models (GMMs) using the Bayesian Information Criterion (BIC). BIC selects **eight components** with the **VVV** covariance structure, in which volume, shape, and orientation vary across components. A separate **three-component GMM** fitted with `scikit-learn` provides a directly comparable partition for TCLUST. A Ward-linkage dendrogram and the silhouette score for the three-component GMM provide complementary views of the training data.

The retained TCLUST solution uses `k = 3`, `alpha = 0.05`, and `restr.fact = 50`. The `ctlcurves` diagnostics compare choices of cluster count and trimming. The notebook uses three profiles for interpretability and comparability. TCLUST labels the trimmed training tests **Profile 0** and the three retained groups **Profiles 1–3**. The notebook compares the partitions on non-trimmed training tests using the adjusted Rand index (ARI) and visualizes the profiles on the first two principal components.

### 3. Logistic Mixed-Effects Modeling (RQ2)

The primary binomial GLMM uses `lme4::glmer` with age and categorical TCLUST profile as fixed effects and a random intercept for `clinic_name`. Profile 0 is the reference category. An extended sensitivity model also includes a random intercept for `subject_id`. The analysis compares their training fit using the Akaike Information Criterion (AIC), summarizes collection-unit heterogeneity using the latent-scale variance partition coefficient (VPC), and applies simulation-based residual diagnostics with `DHARMa`.

The notebook stores the odds-ratio forest plot in `forest_plot_OR.png`. It generates `marginal_effects.png` from the GLMM's fixed-effect estimates and their covariance matrix. The latter figure shows fitted probabilities with the clinic random intercept set to zero and **pointwise 95% Wald confidence intervals**, within the observed training age range of each profile.

### 4. Held-Out Evaluation and Bootstrap (RQ3)

The notebook maps each held-out test to the profile of its nearest non-trimmed training test in the transformed feature space. The cutoff is the 95th percentile of nearest-neighbor distances among non-trimmed training tests, excluding self-matches. Held-out tests beyond that cutoff receive Profile 0.

The predictive comparison uses the clinic-level GLMM and standard logistic regression with the same fixed predictors. Both models are assessed on the same held-out tests using average precision (AP), ROC-AUC, and Brier score. The notebook selects the GLMM classification threshold by maximizing F1 on training predictions and applies it to the held-out tests.

The paired bootstrap uses **500 resamples of held-out patients**, keeping each sampled patient's tests together. The percentile intervals quantify held-out sample variation in the metric differences conditional on the fitted models. Wald intervals quantify uncertainty in the GLMM fixed effects.

## Main Results

### RQ1 · Testing-Process Profiles

| Analysis | Result |
|---|---|
| MCD diagnostic flags | 5,402 training tests (44.29%) |
| BIC-selected GMM | 8 components, VVV covariance structure |
| Three-component comparison GMM | Silhouette score of 0.2089 |
| TCLUST | 3 retained profiles and 610 trimmed tests |
| GMM–TCLUST agreement on non-trimmed tests | ARI of 0.0474 |

The low ARI indicates limited agreement between the partitions. The profiles describe differences in demographics and testing workflows. Profile 0, for example, has longer average turnaround times than each retained profile.

### RQ2 · Association With Positivity

| Fixed effect | Odds ratio | 95% Wald CI |
|---|---:|---|
| Age, per unit | 1.0045 | [0.9983, 1.0107] |
| Profile 1 vs Profile 0 | 0.4096 | [0.2974, 0.5641] |
| Profile 2 vs Profile 0 | 0.6273 | [0.4506, 0.8734] |
| Profile 3 vs Profile 0 | 0.3333 | [0.2246, 0.4946] |

The age interval includes 1. At the same age and clinic, the fitted odds are lower for Profiles 1–3 than for the trimmed reference group. These estimates describe conditional associations.

The clinic random-intercept variance is **2.2074**, corresponding to a **40.15% latent-scale VPC**. Adding a patient random intercept lowers the training AIC from **4,982.92** to **4,925.03**. The predictive comparison evaluates the clinic-only GLMM.

### RQ3 · Held-Out Performance

| Model | Average precision ↑ | ROC-AUC ↑ | Brier score ↓ |
|---|---:|---:|---:|
| Standard logistic regression | 0.0705 | 0.5643 | 0.0582 |
| Clinic-level logistic GLMM | **0.0962** | **0.5989** | **0.0576** |

| Paired difference | Bootstrap mean | 95% percentile CI |
|---|---:|---|
| AP (GLMM minus baseline) | +0.0265 | [0.0074, 0.0510] |
| ROC-AUC (GLMM minus baseline) | +0.0341 | [−0.0019, 0.0722] |
| Brier score (baseline minus GLMM) | +0.0006 | [−0.0002, 0.0014] |

Only the AP interval excludes zero. At the training-selected threshold of **0.0765**, the GLMM detects **71 of 187** positive tests and yields **708 false positives**, with approximately **9% precision**, **38% recall**, and an **F1-score of 0.15** for the positive class.

## Interpretation and Limitations

- **Cluster interpretation.** A patient's repeated tests can fall into different profiles. Skewed turnaround times and the chi-square reference cutoff affect the MCD flagging rate. Held-out profile labels depend on the nearest-neighbor approximation.
- **Validation scope.** Some held-out patients can belong to clinics represented in training. The results measure performance for held-out patients within the study's healthcare setting.
- **Predictive scope.** The GLMM's absolute discrimination remains modest. The false-positive burden at the selected threshold limits the practical value of the predictive gain.

## Repository Structure

| File | Contents |
|---|---|
| [`Advanced_multivariate_statistics.ipynb`](Advanced_multivariate_statistics.ipynb) | Main Colab/Jupyter notebook with Python and R code and saved results |
| [`advanced_multivariate_statistics.py`](advanced_multivariate_statistics.py) | Colab export of the notebook |
| [`report.pdf`](report.pdf) | Final report with methodology, results, figures, and discussion |
| [`Presentazione Advanced.pptx`](Presentazione%20Advanced.pptx) | Final project presentation |
| `report.tex`, `figures/` | LaTeX source and external figures, where included |
| `README.md` | Project documentation |

## Technologies and Libraries

| Environment | Packages |
|---|---|
| Python | `pandas`, `numpy`, `scipy`, `scikit-learn`, `statsmodels`, `pyreadr`, `matplotlib`, `seaborn`, `rpy2` |
| R | `mclust`, `tclust`, `lme4`, `DHARMa`, `sjPlot`, `ggplot2`, `scales` |

`rpy2` connects the Python notebook to R. The final visualization cell uses `sjPlot` for the forest plot and `ggplot2` for the fitted probability curves.

## How to Run the Project

1. Open `Advanced_multivariate_statistics.ipynb` in Google Colab or in a local Jupyter Python kernel.
2. Make the `R` executable available on the kernel's `PATH` and install `rpy2` before running the first code cell. For local runs, install the listed Python packages in the kernel's environment and use an `rpy2` version that supports `pandas2ri.activate()`, which the notebook calls during setup.
3. Run the cells from top to bottom. The notebook installs `pyreadr` and the required R packages through its setup and analysis cells, and it downloads `covid_testing.rda` from the public `medicaldata` repository. Internet access is required for these steps.

The notebook is the executable entry point. The Python export retains IPython shell commands such as `!pip` and `!R`. The saved notebook outputs document the results reported here. Package versions and the runtime can affect exact values on re-execution.

## Data Source

The `covid_testing` dataset is distributed with the R `medicaldata` package. [Package documentation](https://higgi13425.github.io/medicaldata/) · [Source repository](https://github.com/higgi13425/medicaldata)
