# Analysing and Forecasting Economic Time Series

## Swiss government bond yields

This example introduces descriptive analysis, correlograms, autoregressive models, residual checks, and out-of-sample forecast evaluation using monthly Swiss Confederation bond yields. Only the 10-year yield is forecast; the other maturities are available for student exercises.

### Files

- `swiss_bond_yields_lecture.R` — commented lecture script.
- `data/snb_swiss_government_bond_yields_monthly_clean.csv` — monthly 2-, 5-, 10-, 20-, and 30-year yields and the 10-year minus 2-year spread.
- `output/` — charts and forecast-error tables created when the script runs; not included in the repository.

### Data source

The data are end-of-month CHF Swiss Confederation zero-coupon spot yields, in percent per year, from January 2015 through October 2026. Source: [Swiss National Bank, *Spot interest rates on Swiss Confederation bonds and other bond issues – Month*, cube `rendeiduebm`](https://data.snb.ch/en/topics/ziredev/cube/rendeiduebm). The [SNB CSV API](https://data.snb.ch/api/cube/rendeiduebm/data/csv/en?dimSel=D0(CHF),D1(2J,5J,10J,20J,30J)&fromDate=2015-01) supplied the published data. The CSV here has one row per month and adds the 10-year minus 2-year spread.

### Run

Install R and the single add-on package once:

```r
install.packages("forecast")
```

From this folder, run:

```sh
Rscript swiss_bond_yields_lecture.R
```

The script selects AR(1)–AR(6) by BIC using data through 2021, fits an AR model with a mean and no differencing or seasonal terms, and forecasts the 10-year yield for 12 months. It recursively evaluates 1-, 6-, and 12-month forecasts from 2022 onward against an expanding historical-mean benchmark. Inspect the residual checks before interpreting the forecasts.
