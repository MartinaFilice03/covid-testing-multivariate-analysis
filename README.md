# Testing-Process Profiles and Mixed-Effects Modeling: A Robust Multivariate Analysis of COVID-19 Diagnostic Data

Repository containing the code, methodology, and final report for the **Advanced Multivariate Statistics** project (Academic Year 2025–2026).

---

## Project Overview
This project analyzes the `covid_testing` medical dataset (sourced from the R `medicaldata` package) to extract operational and predictive insights while rigorously managing the complexities inherent to healthcare data—specifically multivariate outliers, severe class imbalance, and hierarchical facility-level structures.

The investigation is structured around three primary Research Questions (RQs):
1. **RQ1 (Robust Latent Profiles):** Identification of latent testing-process profiles using robust multivariate clustering, comparing standard Gaussian Mixture Models (GMM) with robust trimmed clustering (TCLUST).
2. **RQ2 (Hierarchical Mixed Effects):** Evaluation of patient-level characteristics on diagnostic outcomes using a Generalized Linear Mixed Model (GLMM) with clinic-level random intercepts (`1 | clinic_name`) to isolate facility-level heterogeneity.
3. **RQ3 (Out-of-Sample Predictive Validation):** Assessment of out-of-sample generalizability under severe class imbalance using PR-AUC, ROC-AUC, Brier scores, and rigorous data-leakage safeguards.

---

## Repository Structure
```text
├── Advanced_multivariate_statistics.ipynb   # Main Jupyter Notebook (Python/R integration via rpy2)
├── report.tex                               # Source code of the formal LaTeX report
├── report.pdf                               # Compiled final PDF report
└── README.md                                # Project documentation
