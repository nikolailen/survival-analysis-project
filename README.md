# Survival Analysis of Heart Transplant Patients

🌐 [Project page](https://nikolailen.github.io/survival-analysis-project/)

👤 Project contact: [Nikolai Len](https://www.linkedin.com/in/niklen/)

## Overview

This project analyzes survival outcomes for heart transplant candidates using the `survival` package in R. It evaluates how transplant status, surgery, and age groups affect patient survival over time.

## Dataset and Methods

- Dataset: JASA/heart transplant data from the `survival` package
- Core methods:
  - Kaplan-Meier survival estimation
  - Log-rank tests for group comparison
  - Nelson-Aalen cumulative hazard analysis
  - Cox proportional hazards modeling

## Repository Structure

- `Finalized research/Survival_Analysis.Rmd`: analysis source notebook
- `Finalized research/Survival_Analysis.pdf`: original report export
- `index.md`: GitHub Pages entrypoint generated from the R Markdown report
- `index_files/`: generated figure assets used by `index.md`

## Reproducibility

Render the project page from the R Markdown source:

```powershell
Rscript -e "rmarkdown::render('Finalized research/Survival_Analysis.Rmd', output_format='github_document', output_file='index.md', output_dir='.')"
```

## License

This project is distributed under the terms of the MIT License. See `LICENSE`.
