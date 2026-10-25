# Project context

Master's Project 1: "Statistical Analysis of SPT-N Values and Shear Strength Parameters in Dhaka".

- Data: `data/Dhaka qu, spt N, depth.csv`, columns Depth (m), SPT N, Unconfined compressive strength qu (kPa).
- 134 samples, Dhaka only. Field SPT-N only (no N60 correction, no PI).
- Focus is qu vs SPT-N, not friction angle.
- Sampling areas: Kawlar, Ashkona, Purbachal, 300 ft road, and MRT Line 6 (Uttara to Motijheel).

# Rules for every script

- Write numbered scripts in `R/` (01_..., 02_..., ...).
- Save figures as 300-dpi PNG to `outputs/figures`, named `Fig4_x_...`.
- Save tables as CSV to `outputs/tables`, named `Table4_x_...`.
- Use ggplot2 with theme_bw and Times New Roman.
- Label axes with units: SPT-N (blows/300 mm), qu (kPa), Depth (m).
- Use set.seed(123).
- Print a short plain-English summary of each result to the console.

# First task

Create folders `R/`, `outputs/figures`, `outputs/tables`, and a script `R/00_packages.R` that installs (if missing) and loads:
tidyverse, janitor, psych, moments, car, lmtest, corrplot, GGally, caret, rsample, boot, broom, knitr, flextable, officer.
