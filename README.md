# MAT2007 Programming Project: E. coli Antibiotic Resistance Analysis

This project explores a public dataset of antibiotic resistance in *Escherichia coli* (E. coli), obtained from BV-BRC and published on Kaggle.

## Research Question
**Which antibiotics have the highest resistance rates in E. coli?**

## Dataset
- Source: [E. coli Resistance Dataset (Kaggle)](https://www.kaggle.com/datasets/valeriamaciel/e-coli-resistance-dataset)
- File used: `BVBRC_E.coli_Dataset.csv`
- The dataset contains 195,795 observations and 21 variables.

## What the script does
- Loads and explores the dataset
- Cleans the data by removing observations with missing antibiotic or resistance phenotype information
- Calculates resistance rates and standard errors for each antibiotic
- Excludes antibiotics with fewer than 100 observations
- Identifies the top 10 antibiotics with the highest resistance rates
- Calculates the resistance phenotype distributions for the top 10 antibiotics
- Calculates the Pearson correlation between sample size and resistance rate
- Visualizes:
  - Top 10 antibiotics by resistance rate, including standard error bars
  - Relationship between sample size and resistance rate

## Output
The script generates:
- A table of the top 10 antibiotics by resistance rate
- A table of resistance phenotype distributions for the top 10 antibiotics
- A horizontal bar chart showing the top 10 resistance rates
- A scatter plot showing the relationship between sample size and resistance rate

## How to run
Requires Python with `pandas`, `numpy`, and `matplotlib` installed.

Place `BVBRC_E.coli_Dataset.csv` in the same folder as the Python script.

```bash
pip install pandas numpy matplotlib
python antibiotic_resistance.py
```
