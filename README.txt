Data and code in NJ-serology repository at https://github.com/freedmanari
- All code by Ari Freedman
- Scripts (not included) for acquiring data from New Jersey Department of Health all written by Justin Sheen

--------------------------

Data

From NJ Department of Health (all of these except delay_dists.csv are available only upon request):
- q1_pcr.csv: weekly first-time COVID-19 PCR positives by NJ county in study period
- q1_sero.csv: weekly first-time COVID-19 serology positives by NJ county in study period
- q2_pcr.csv: weekly first-time PCR positives with past PCR positive, by NJ county
- q2_sero.csv: weekly first-time serology positives with past PCR positive, by NJ county
- q3_pcr.csv: weekly first-time PCR positives without past PCR positive, by NJ county
- q3_sero.csv: weekly first-time serology positives without past PCR positive, by NJ county
- waiting_times.csv: for each week and county, lists the delays in weeks from positive PCR to positive serology for all positive serology results from that week with past PCR positive
- sero_tests_by_day.csv: daily number of serology tests given in NJ, broken up by positives and negatives
- sero_tests.csv: all COVID-19 serology tests from New Jersey in study period, with titer value and first PCR positive date given for each test when applicable
- model_sero.csv: just the serology tests applicable for the serology model (derived from sero_tests.csv)
- vaccinations_by_age.csv: weekly NJ COVID-19 vaccinations for different age groups
- delay_dists.csv: matrix of backwards delay distributions for serology tests taken in each test week, i.e. for each week the distribution of delays from past positive PCR test to serology tests taken in that week (derived from sero_tests.csv)

Data from other sources:
- co-est2021-alldata.csv: 2021 estimated NJ census data by county, from https://www.census.gov/data/tables/time-series/demo/popest/2020s-counties-total.html
- sc-est2021-agesex-civ.csv: 2021 estimated NJ census data by age and sex, from https://www.census.gov/data/tables/time-series/demo/popest/2020s-state-detail.html
- owid-covid-data.csv: weekly CFR estimates from US, from https://ourworldindata.org/mortality-risk-covid
- NJ_incidence_by_age.csv: weekly COVID-19 incidence in NJ per 100,000 by age group, from https://covid.cdc.gov/covid-data-tracker/#demographicsovertime
- Provisional_COVID-19_Death_Counts_by_Week_Ending_Date_and_State.csv: weekly COVID-19 death counts by state from https://data.cdc.gov/NCHS/Provisional-COVID-19-Death-Counts-by-Week-Ending-D/r8kw-7aab

--------------------------

Code

R files (to be run in this order):
- thetas.R: initializes most data sources, simulates COVID-19 true incidence curves, and simulates the log-odds ratios (thetas)
- thetas_plots.R: plots the simulated incidence curves and log-odds ratios
- summary_plots.R: plots various summary figures relating to NJ COVID-19 testing volume, positivity, and behavior for both PCR and serology tests, comparing different regions of NJ
- sero_model.R: initializes the rest of the data sources and runs the serology model for both the main text and the supplementary sensitivity analyses
- sero_model_plots.R: plots some of the simulated data fed into the serology model as well as the results of the serology model and its fitted parameters, both for the main text model and the supplementary models

Stan files:
- sero_model_prev_PCR.stan: encodes the serology submodel for just serology tests with past PCR positive
- sero_model.stan: encodes the full serology submodel covering all serology tests from model_sero.csv

