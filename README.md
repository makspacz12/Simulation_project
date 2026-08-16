# ST3247 Simulation — Adaptive-Network SIR Inference

**National University of Singapore** · Department of Statistics and Data Science  
**Module:** [ST3247 Simulation](https://nusmods.com/courses/ST3247/simulation) (4 Units)  
**Instructor:** [Assoc. Prof. Alexandre Thiéry](https://alexxthiery.github.io/)  
**Project brief:** [Simulation-Based Inference for an Adaptive-Network Epidemic Model](https://alexxthiery.github.io/teaching/SBI_infection/SBI-infection.html)

Group project by **Maksymilian Paczyński** and **Ruofei Fang**.

---

## Overview

This repository is the **graded group project** for [ST3247 Simulation](https://nusmods.com/courses/ST3247/simulation) (AY25/26 Semester 2) under **Assoc. Prof. Alexandre Thiéry**.

The assessment was to run a full **simulation-based inference** pipeline on a stochastic SIR epidemic on an **adaptive contact network** (Gross et al., 2006), where the likelihood is intractable. In practice that meant implementing and stress-testing the same ideas examined in the module — especially **Monte Carlo simulation**, **Bayesian inference**, and **Approximate Bayesian Computation (ABC)** — then pushing beyond basic rejection ABC with advanced methods.

**Primary deliverables:** written **report** (`report.pdf`) + reproducible **code / notebook** (`simulator.ipynb`). Project marking (per brief): **70% report** · **30% code**.

| Parameter | Meaning | Prior |
|---|---|---|
| β | Infection probability per S–I edge per step | Uniform(0.05, 0.50) |
| γ | Recovery probability per infected node per step | Uniform(0.02, 0.20) |
| ρ | Rewiring probability per S–I edge per step | Uniform(0.0, 0.8) |

**Final estimates** (regression-adjusted ABC, ε = 10%):

| Parameter | Median | 95% credible interval | Shrinkage |
|---|---:|---|---:|
| β | 0.176 | [0.108, 0.284] | 0.664 |
| γ | 0.089 | [0.066, 0.119] | 0.747 |
| ρ | 0.316 | [0.244, 0.393] | 0.834 |

Informed 8-statistic design resolves the β–ρ confound: posterior correlation \(r(\beta,\rho) = 0.071\) at ε = 5% (vs 0.818 for naive temporal means).

---

## Course assessment & syllabus (what this project tested)

Under Prof. Thiéry the module was heavily **implementation-first** (Python / NumPy). Typical component weights for AY25/26 Sem 2:

| Component | Weight |
|---|---:|
| Canvas quizzes | 10% |
| Midterm | 20% |
| **Group project (this repo)** | **30%** |
| Final exam | 40% |

**Topics covered in the course** (and exercised by the project / exams):

1. Monte Carlo estimation  
2. Inverse transform sampling and rejection sampling  
3. Importance sampling (including self-normalised IS)  
4. Bayesian inference  
5. **Approximate Bayesian Computation (ABC)** — core of this project; also relevant to finals preparation  
6. Markov chain basics  
7. Markov Chain Monte Carlo (lecture coverage varied by semester; we still implemented **ABC-MCMC** as an advanced project method)

The project brief required: (i) basic **rejection ABC**, (ii) careful **summary-statistic design** (β and ρ are mechanistically confounded), and (iii) at least one **advanced SBI method** (we used regression adjustment + ABC-MCMC), plus validation.

Midterm performance: **A** · list **A (Honours)** — context for the course standard under Prof. Thiéry.

---

## Repository layout

```text
.
├── report.pdf              # Final project report (primary deliverable)
├── report.tex              # LaTeX source
├── references.bib
├── simulator.ipynb         # Full analysis pipeline (exploration → ABC → reg-adj → MCMC)
├── data/                   # Observed trajectories (40 replicates)
│   ├── infected_timeseries.csv
│   ├── rewiring_timeseries.csv
│   └── final_degree_histograms.csv
└── graphics/               # Figures embedded in the report
```

Large simulation caches (`*.npz`) are gitignored and regenerated on first run.

---

## What we did

1. **Rejection ABC** — 10,000 prior simulations; baseline vs informed summary statistics; identifiability via \(|r(\beta,\rho)| < 0.3\).
2. **Summary-statistic design** — eight mechanistic probes exploiting Infection → Recovery → Rewiring phase order (pure-β first step, early rewire ratio for ρ, post-peak decay for γ).
3. **Regression adjustment** (Beaumont, Zhang & Balding, 2002) — primary reported posterior.
4. **ABC-MCMC** (Marjoram et al., 2003) — independent cross-check at a tighter tolerance.
5. **Validation** — synthetic-truth recovery and posterior predictive checks.

---

## How to reproduce

```bash
# 1. Clone
git clone https://github.com/makspacz12/st3247-adaptive-sir-abc.git
cd st3247-adaptive-sir-abc

# 2. Open and run the notebook (Python 3 + numpy, pandas, matplotlib, scikit-learn)
jupyter notebook simulator.ipynb
```

On first run the notebook builds a prior bank of 10,000 simulations (~several minutes) and caches it under `data/`. ABC-MCMC is the slowest step (~1–2 hours with the notebook settings).

Rebuild the PDF (optional):

```bash
# requires a TeX distribution (e.g. TeX Live / MiKTeX / tectonic)
pdflatex report.tex && bibtex report && pdflatex report.tex && pdflatex report.tex
```

---

## Authors

- Maksymilian Paczyński  
- Ruofei Fang  

Group project for **ST3247 Simulation**, NUS, under **Assoc. Prof. Alexandre Thiéry**.

---

## References (selected)

- Gross, T., D’Lima, C. J. D., & Blasius, B. (2006). Epidemic dynamics on an adaptive network. *PRL*.  
- Beaumont, M. A., Zhang, W., & Balding, D. J. (2002). Approximate Bayesian computation in population genetics. *Genetics*.  
- Marjoram, P. et al. (2003). Markov chain Monte Carlo without likelihoods. *PNAS*.  
- Course page: [alexxthiery.github.io/teaching/SBI_infection](https://alexxthiery.github.io/teaching/SBI_infection/SBI-infection.html)

Full bibliography: `references.bib`.
