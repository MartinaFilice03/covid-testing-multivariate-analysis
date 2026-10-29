
# Testing-Process Profiles and Mixed-Effects Modeling

### A Robust Multivariate Analysis of COVID-19 Diagnostic Data

**Advanced Multivariate Statistics — Academic Year 2025–2026**

**Authors:** Corinne D'Elia and Martina Filice

---

## Project Overview

This project investigates the `covid_testing` dataset from the R [`medicaldata`](https://higgi13425.github.io/medicaldata/) package using robust multivariate statistical methods, model-based clustering, and logistic mixed-effects regression.

The main objective is to identify interpretable demographic and diagnostic testing-process profiles, investigate heterogeneity across clinics and specimen-collection units, and evaluate the retrospective predictive performance of a hierarchical logistic model.

The analysis addresses several statistical challenges commonly encountered in real-world healthcare data:

- Multivariate outliers and highly skewed distributions.
- Severe class imbalance between positive and negative COVID-19 test results.
- Repeated diagnostic tests belonging to the same patient.
- Potential dependence among observations collected within the same clinic.
- The need for robust clustering and reliable out-of-sample predictive evaluation.

The project combines Python and R through `rpy2` to integrate statistical methods available in both environments.

**Important:** Each observation represents a diagnostic test rather than necessarily a distinct patient. Furthermore, some testing-process variables may only become available after the diagnostic procedure. Therefore, this study focuses on retrospective statistical analysis and does not establish a prospective clinical early-warning system.

---

## Research Questions

The investigation is structured around three main Research Questions (RQs).

### RQ1 — Robust Latent Testing-Process Profiles

What demographic and testing-process profiles can be identified through multivariate clustering, and how do conventional Gaussian Mixture Models (GMM) compare with robust trimmed clustering (TCLUST)?

The analysis includes:

- Multivariate outlier detection using Minimum Covariance Determinant (MCD).
- Gaussian Mixture Models and model selection through BIC.
- Robust clustering through TCLUST.
- Comparison of clustering partitions using the Adjusted Rand Index (ARI).
- Principal Component Analysis (PCA) for visual interpretation.

### RQ2 — Logistic Mixed-Effects Modeling

How are patient age and testing-process profiles associated with COVID-19 test positivity, and how much heterogeneity is associated with clinics and repeated patient measurements?

A logistic Generalized Linear Mixed Model (GLMM) is estimated using:

- **Fixed effects:** patient age and TCLUST cluster membership.
- **Random effect:** clinic or specimen-collection unit (`clinic_name`).
- **Sensitivity analysis:** an additional patient-level random intercept (`subject_id`).

The Variance Partition Coefficient (VPC) is used to describe clinic-level heterogeneity on the latent logistic scale.

### RQ3 — Out-of-Sample Predictive Validation

Does the clinic-level logistic GLMM provide better retrospective predictive performance than standard logistic regression on previously unseen patients?

The models are compared using:

- Average Precision (reported as PR-AUC).
- Receiver Operating Characteristic Area Under the Curve (ROC-AUC).
- Brier Score.
- Precision, recall, and F1-score.
- Confusion matrix.
- Patient-level bootstrap confidence intervals for predictive performance differences.

---

## Dataset Description

The analysis uses the publicly available `covid_testing` dataset from the R `medicaldata` package.

The original dataset contains **15,524 diagnostic-test records and 17 variables**.

After removing invalid test outcomes and selecting complete observations for the analysis, **15,223 valid diagnostic tests** remain.

The binary target variable is `result_bin`:

- `1`: Positive COVID-19 test.
- `0`: Negative COVID-19 test.

### Selected Variables

| Variable | Description | Analytical Role |
|---|---|---|
| `age` | Patient age | Clustering feature and GLMM fixed effect |
| `pan_day` | Day within the recorded pandemic timeline | Clustering feature |
| `col_rec_tat` | Collection-to-receipt turnaround time | Clustering feature |
| `rec_ver_tat` | Receipt-to-verification turnaround time | Clustering feature |
| `clinic_name` | Clinic or specimen-collection unit | GLMM random-effect grouping factor |
| `subject_id` | Patient identifier | Train/test splitting and sensitivity analysis |
| `result_bin` | Positive or negative diagnostic outcome | Binary response variable |

The PCR cycle-threshold variable `ct_result` is excluded from predictive features to prevent direct target leakage.

The variables `payor_group` and `patient_class` are also excluded because of extensive missingness.

### Patient-Aware Train/Test Split

Since individual patients may undergo multiple diagnostic tests, the dataset is divided using unique patient identifiers rather than randomly splitting individual records.

An 80/20 split is performed with `random_state = 42`.

| Partition | Number of Tests |
|---|---:|
| Training set | 12,196 |
| Test set | 3,027 |
| Total | 15,223 |

The held-out test set contains:

- 187 positive tests.
- 2,840 negative tests.
- Approximately 6.18% positive outcomes.

This approach prevents the same patient from appearing in both training and testing partitions.

The split is patient-aware but is not stratified by the outcome variable.

---

## Methodology

### 1. Exploratory Data Analysis and Preprocessing

The preliminary analysis investigates the characteristics and quality of the dataset.

It includes:

- Distribution of positive, negative, and invalid test outcomes.
- Missing-value analysis.
- Distribution of observations across clinics.
- Identification of repeated tests belonging to the same patient.
- Spearman correlation analysis.
- Variance Inflation Factor (VIF) assessment.
- Descriptive analysis of continuous variables.

To identify potentially anomalous multivariate observations, robust Mahalanobis distances are computed using the Minimum Covariance Determinant (MCD) estimator.

A chi-square reference threshold at the 99th percentile is used for diagnostic flagging.

The continuous clustering features are then transformed using Yeo–Johnson transformations and standardized using `StandardScaler`.

Both preprocessing procedures are fitted exclusively on the training data to prevent information leakage into the evaluation stage.

### 2. Model-Based Clustering (GMM)

Gaussian Mixture Models are applied to identify latent structure within the transformed and standardized data.

The R package `mclust` is used to compare alternative mixture specifications and select a model using the Bayesian Information Criterion (BIC).

The BIC-selected model contains **eight Gaussian components** with a VVV covariance structure.

A separate conventional Gaussian Mixture Model with three components is fitted using `scikit-learn` to facilitate comparison with TCLUST.

Ward hierarchical clustering and the Silhouette Score are used as complementary exploratory diagnostics.

### 3. Robust Clustering (TCLUST)

TCLUST is applied to reduce the influence of extreme observations during cluster estimation.

The selected configuration is:

- Number of clusters: `k = 3`.
- Trimming proportion: `alpha = 0.05`.
- Eigenvalue restriction factor: `restr.fact = 50`.

The trimmed observations are assigned to **Cluster 0**, while the three retained clusters are labeled 1, 2, and 3.

TCLUST diagnostic curves (`ctlcurves`) are examined across different numbers of clusters and trimming proportions.

The resulting profiles are summarized using their demographic and testing-process characteristics.

The Adjusted Rand Index (ARI) is used to compare TCLUST assignments with the conventional three-component GMM, excluding trimmed observations.

Principal Component Analysis (PCA) provides a two-dimensional visualization of the cluster structure.

### 4. Logistic Mixed-Effects Modeling (GLMM)

A binomial logistic Generalized Linear Mixed Model is fitted using the R package `lme4`.

The primary model includes:

- Patient age as a fixed effect.
- Categorical TCLUST profile membership as a fixed effect.
- Clinic-level random intercepts through `(1 | clinic_name)`.

The analysis reports:

- Fixed-effect regression coefficients.
- Odds Ratios (OR).
- Wald 95% confidence intervals.
- Statistical significance of fixed effects.
- Clinic-level random-effect variance.
- Variance Partition Coefficient (VPC).

An extended model additionally includes `(1 | subject_id)` to investigate dependence between repeated tests belonging to the same patient.

The primary and extended models are compared using the Akaike Information Criterion (AIC).

Simulation-based residual diagnostics are performed using the R package `DHARMa`.

Odds-ratio forest plots and model-based predicted probabilities are generated using `sjPlot` and `ggeffects`.

### 5. Out-of-Sample Predictive Evaluation

The trained models are evaluated on a held-out test set containing patients not present in the training data.

Test observations are transformed using the preprocessing parameters learned exclusively from the training set.

Since the notebook does not implement a native TCLUST prediction procedure, test-set cluster membership is approximated using a one-nearest-neighbor classifier fitted on non-trimmed training observations.

A training-derived nearest-neighbor distance threshold is used to assign sufficiently distant test observations to Cluster 0.

The clinic-level GLMM is compared with a standard logistic regression containing the same fixed-effect predictors.

Predictive performance is assessed using Average Precision, ROC-AUC, Brier Score, and threshold-based classification metrics.

The classification threshold is selected by maximizing F1-score on in-sample training predictions and subsequently applied unchanged to the test set.

### 6. Bootstrap and Statistical Uncertainty

A nonparametric bootstrap procedure is implemented to assess uncertainty in the differences between the predictive performance of the GLMM and standard logistic regression.

The analysis uses:

- 500 bootstrap resamples.
- Resampling at the patient level.
- Paired comparison of the two models.
- 95% percentile confidence intervals.

All test observations associated with a selected patient are retained together during resampling.

The bootstrap evaluates differences between predictions from already-fitted models; the complete modeling pipeline is not refitted in each iteration.

Uncertainty in the GLMM fixed-effect coefficients is evaluated separately using Wald confidence intervals.

---

## Main Results

### Robust Outlier Detection and Clustering

| Analysis | Result |
|---|---|
| MCD-flagged observations | 5,402 training tests (44.29%) |
| BIC-selected GMM | 8 components, VVV covariance structure |
| Conventional comparison GMM | 3 components |
| GMM Silhouette Score | 0.2089 |
| TCLUST configuration | 3 clusters, 5% trimming |
| TCLUST trimmed observations | 610 tests |
| GMM–TCLUST Adjusted Rand Index | 0.0474 |

The low ARI indicates limited agreement between the conventional GMM and TCLUST solutions.

The identified clusters primarily characterize demographic and diagnostic-process patterns rather than independently validated clinical risk groups.

### Mixed-Effects Modeling

The primary GLMM identifies substantial heterogeneity across clinic or specimen-collection units.

| Metric | Result |
|---|---:|
| Clinic random-effect variance | Approximately 2.207 |
| Latent-scale VPC | 40.15% |
| AIC — clinic-level GLMM | 4,982.92 |
| AIC — clinic and patient GLMM | 4,925.03 |

The extended model achieves a lower AIC, suggesting improved in-sample fit after accounting for repeated patient measurements.

However, a lower AIC does not automatically imply superior out-of-sample predictive performance.

### Predictive Performance

The following results are obtained on the held-out test set.

| Model | Average Precision / PR-AUC | ROC-AUC | Brier Score |
|---|---:|---:|---:|
| Standard Logistic Regression | 0.0705 | 0.5643 | 0.0582 |
| Clinic-Level Logistic GLMM | **0.0962** | **0.5989** | **0.0576** |

The GLMM achieves better performance across the three reported metrics, although the absolute discrimination remains limited.

### Bootstrap Confidence Intervals

| Performance Difference | Mean Difference | 95% Bootstrap CI |
|---|---:|---|
| PR-AUC improvement | 0.0265 | [0.0074, 0.0510] |
| ROC-AUC improvement | 0.0341 | [-0.0019, 0.0722] |
| Brier Score improvement | 0.0006 | [-0.0002, 0.0014] |

The confidence interval for the PR-AUC improvement excludes zero.

The corresponding intervals for ROC-AUC and Brier Score include zero, indicating that the evidence for improvement is less conclusive for these metrics.

Overall, the clinic-level GLMM demonstrates a modest relative predictive advantage, but the results do not support clinical deployment.

---

## Repository Structure

```text
├── Advanced_multivariate_statistics.ipynb  # Main Jupyter Notebook (Python + R)
├── advanced_multivariate_statistics.py     # Python export from Google Colab
├── report.tex                              # LaTeX source of the final report
├── report.pdf                              # Compiled academic report
├── figures/                                # Figures referenced by the LaTeX report, if external
└── README.md                               # Project documentation
```

The notebook is the primary executable analysis.

The `.py` file is an export of the Google Colab notebook and contains notebook-specific shell commands.

The `figures/` directory should be included if the LaTeX source references external image files.

---

## Technologies and Libraries

The project combines Python and R through the `rpy2` interface.

### Python

- `pandas` — data manipulation and cleaning.
- `numpy` — numerical operations and bootstrap resampling.
- `scipy` — statistical distributions and hierarchical clustering.
- `scikit-learn` — preprocessing, GMM, PCA, logistic regression, and predictive metrics.
- `statsmodels` — VIF calculations and standard logistic regression.
- `pyreadr` — importing the original R dataset.
- `matplotlib` and `seaborn` — statistical visualization.
- `rpy2` — integration between Python and R.

### R

- `mclust` — Gaussian model-based clustering.
- `tclust` — robust trimmed clustering.
- `lme4` — logistic generalized linear mixed models.
- `DHARMa` — simulation-based residual diagnostics.
- `sjPlot` — odds-ratio visualization.
- `ggeffects` — marginal predictions.
- `ggplot2` — graphical representations.

---

## How to Run the Project

### Recommended Environment: Google Colab

The project was developed using Google Colab with Python and R integration.

1. Open `Advanced_multivariate_statistics.ipynb` in Google Colab.
2. Ensure that R and the required Python libraries are available.
3. Execute the notebook cells sequentially.
4. Allow the notebook to install any missing dependencies.
5. Run the complete analysis to reproduce the tables, statistical results, and visualizations.

The notebook automatically downloads the `covid_testing.rda` dataset from the upstream GitHub repository.

An internet connection is required for the initial dataset download and installation of dependencies.

### Local Execution

A local environment requires Python, Jupyter Notebook, R, `rpy2`, and the packages listed above.

The exported `advanced_multivariate_statistics.py` contains Google Colab/Jupyter-specific shell commands such as `!pip` and `!R`.

Consequently, the Python export should be treated primarily as a readable representation of the notebook and is not guaranteed to execute unchanged using:

```bash
python advanced_multivariate_statistics.py
```

The Jupyter notebook is the recommended executable version.

Random seeds are used for the principal stochastic procedures, but exact results may depend on package versions and the execution environment.

---

## Interpretation and Limitations

Several limitations must be considered when interpreting the findings.

**Testing-process rather than clinical profiles:** The clustering variables describe demographic information and diagnostic workflow timing. The resulting clusters should not be interpreted as confirmed clinical severity or disease-risk categories.

**Sensitivity to clustering choices:** The GMM and TCLUST solutions exhibit low agreement. This does not independently establish the superiority of one method.

**Multivariate outliers:** The large proportion of MCD-flagged observations reflects distributional characteristics and the reference cutoff rather than proving widespread recording errors.

**Approximate test-set cluster assignment:** Held-out observations are assigned to TCLUST profiles using a nearest-neighbor approximation instead of a native TCLUST prediction procedure.

**Repeated patient observations:** Some patients contribute multiple tests. The analysis addresses this through patient-aware splitting, sensitivity modeling, and patient-level bootstrap resampling.

**Limited generalizability:** Held-out patients can belong to clinics already represented during training. The results do not demonstrate transportability to entirely new healthcare institutions.

**Retrospective evaluation:** Some operational variables may only be available after particular stages of diagnostic processing. Therefore, the model is not validated for prospective COVID-19 screening or early clinical prediction.

**Bootstrap scope:** Bootstrap confidence intervals quantify variation in predictive performance differences conditional on the fitted models. They do not include uncertainty from refitting the complete analysis pipeline.

These limitations do not invalidate the statistical analyses but restrict the conclusions that can be drawn from them.

---

## Conclusions

This project demonstrates an integrated statistical workflow combining robust multivariate diagnostics, Gaussian model-based clustering, TCLUST, and logistic mixed-effects modeling.

The clustering analyses identify heterogeneous demographic and testing-process profiles, while the mixed-effects models capture meaningful clinic-level heterogeneity and examine the role of repeated patient observations.

The clinic-level GLMM achieves a modest improvement in out-of-sample predictive performance compared with standard logistic regression, particularly in Average Precision.

However, absolute predictive performance remains limited, and the model cannot be considered a validated clinical decision-making tool.

The main contribution of the project is methodological: it illustrates the importance of robust exploratory analysis, appropriate grouping structures, patient-aware validation, and careful interpretation of predictive uncertainty in real-world diagnostic data.

---

## Report and Data Source

The complete methodology, statistical results, figures, and discussion are presented in the final report:

- [Final Report (PDF)](report.pdf)
- [LaTeX Source](report.tex)

**Dataset:** `covid_testing`, R `medicaldata` package.

**Documentation:** https://higgi13425.github.io/medicaldata/

**Source repository:** https://github.com/higgi13425/medicaldata

---

## Authors

**Corinne D'Elia**  
**Martina Filice**

Advanced Multivariate Statistics  
Academic Year 2025-2026  
Prof. Chiara Masci
