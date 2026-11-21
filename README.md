# Testing-the-robustness-of-mobility-as-a-mechanism-for-the-sustainability-of-demand-sharing
This is the code for 'Testing the robustness of mobility as a mechanism for the sustainability of demand sharing' by Erskine and Dyble

Code and data for **Erskine, S. & Dyble, M.** *Testing the robustness of mobility as a mechanism for the sustainability of demand sharing.* [citation / DOI to be added on publication]

This repository contains an R reproduction of the agent-based model of Lewis *et al.* (2014), together with the parameter sensitivity analysis and model amendments reported in the paper.

> Lewis HM, Vinicius L, Strods J, Mace R, Migliano AB. 2014. High mobility explains demand sharing and enforced cooperation in egalitarian hunter-gatherers. *Nature Communications* **5**, 5789. [doi:10.1038/ncomms6789](https://doi.org/10.1038/ncomms6789)

---

## Overview

Lewis *et al.* (2014) proposed that high residential mobility ("walking away" from free-riders) can sustain a demand-sharing hunter-gatherer population that would otherwise collapse through the proliferation of non-hunting free-riders. No code was released with the original model, so we rebuilt it from the published description and then asked how robust its conclusions are.

The analysis has three stages:

1. **Reproduction.** Our implementation closely matches the original results (VEcv = 99.2%).
2. **Global parameter sensitivity analysis.** 20,000 simulations in which all eight tested parameters are sampled at random from the ranges in the table below.
3. **Theoretical evaluation through model amendment.** Core assumptions about mobility and hunting behaviour are relaxed one at a time.

**Headline result:** the effect of mobility on population survival is robust to extensive parameter variation and to relaxing most key assumptions. It is sensitive to two assumptions: the starting proportion of hunters (the population must start with at least c.25% hunters to reliably survive) and the use of discrete alleles (if the propensity to hunt is a continuous heritable trait, populations eventually collapse, with none surviving to 15,000 years).

---

## Repository contents

| File | Description | Corresponds to |
| --- | --- | --- |
| `Base model re creation.R` | Reproduction of the Lewis *et al.* (2014) model. Parameters can be fixed at base values, or sampled randomly for sensitivity testing (see [Running the model](#running-the-model)) | Step 1 (Reproduction) and Step 2 (Parameter sensitivity) |
| `FN 20k sims full results.xlsx` | Output of the 20,000 global sensitivity simulations ("fix none": all parameters varied) | Step 2, Figure 2 |
| `Variable threshold amendment.R` | The migration threshold `thresh` is an individual-level, heritable trait subject to mutation | Mobility amendment 1, Figure 3 |
| `Groups in lattice variable MC.R` | The 20 locations are arranged on a 5 x 4 lattice and migration cost scales with Euclidean distance | Mobility amendment 2 |
| `Starting proportion of hunters.R` | The starting proportion of hunters is sampled between 0 and 1 | Hunting amendment 1, Figure 4A |
| `Binary w 0.05 0.875.R` | Hunting propensity is binary: every agent has the mean propensity of its phenotype (hunters 0.875, free-riders 0.05) | Hunting amendment 2 |
| `Continuous propensity.R` | Hunting propensity `w` is a continuous, heritable trait with no H/F alleles; run for extended durations | Hunting amendment 2, Figure 4B |

---

## Requirements

- **R** (4.0 or later recommended)
- R packages:
  - `future`
  - `future.apply`

Install the dependencies with:

```r
install.packages(c("future", "future.apply"))
```

The scripts run simulations in parallel using `future::plan(multisession, workers = parallel::detectCores() - 1)`, so they will use all but one of your CPU cores. A single simulation covers 1,500 years (547,500 daily timesteps), so large batches are computationally expensive; the 20,000-simulation sensitivity run is best done on a multi-core machine or cluster.

---

## Running the model

Each script is self-contained. Open it in R or RStudio and run it from top to bottom; results are returned as a data frame (`data`) with one row per simulation.

Two settings at the top of each script control the run:

```r
sims  <- 200    # number of simulations
years <- 1500   # simulated duration in years (365 days per year)
```

### Base values vs. random parameter sampling

In `Base model re creation.R`, each parameter has a fixed base-value line followed by a commented-out line that samples it at random. To run the **base model**, leave the file as it is. To run the **global sensitivity analysis**, comment out the fixed lines and uncomment the `runif()` lines, for example:

```r
# dens <- 0.9
dens = runif(1, min = 0, max = 1)
```

For the sedentary comparison, set `thresh = -20` (so agents never migrate).

### Model parameters

| Parameter | Description | Base value | Range tested in sensitivity analysis |
| --- | --- | --- | --- |
| `dens` | Density effect on hunting success | 0.9 | 0 to 1 |
| `HC` | Hunting energy cost | 0.4 | 0 to 2 |
| `MC` | Movement (migration) energy cost | 0.3 | 0 to 2 |
| `thresh` | Movement threshold | 0.5 | -5 to 5 |
| `peak` | Age at peak hunting ability | 35 | 25 to 45 |
| `RepT` | Reproduction threshold | 13.5 | 5 to 15 |
| `k` | Value of a kill | 4 | 0.5 to 10 |
| `age` | Age of independence | 15 | 13 to 18 |

Parameters held fixed throughout, as in the original model: daily energy cost `DC = 0.5`, maximum energy `f = 15`, mutation rate `mut = 0.005`, Gompertz-Makeham mortality parameters `a = 2.7e-7` and `b = 2.4e-4`, 20 locations, and a starting population of 100 agents.

---

## Output

Each script returns a data frame with one row per simulation. The columns in `Base model re creation.R` (and in `FN 20k sims full results.xlsx`) are:

| Column | Description |
| --- | --- |
| `dens_val`, `MC_val`, `HC_val`, `thresh_val`, `peak_val`, `RepT_val`, `k_val`, `age_val` | Parameter values used in the simulation |
| `time_to_collapse` | Day on which the population reached zero. Populations that survive the full run are recorded as 547,500 days (1,500 years) |
| `final_population_size` | Population size at the end of the simulation (0 if collapsed) |
| `final_prop_hunters`, `final_prop_free_riders` | Proportion of agents with the H and F alleles at the end of the simulation |
| `total_migrations` | Total number of migrations |
| `avg_days_between_migrations` | Mean number of days between migrations for an individual |
| `final_mean_energy` | Mean agent energy at the end of the simulation |
| `final_group_size_variance`, `final_mean_occupied_group_size` | Group size statistics at the end of the simulation |
| `final_mean_thresh`, `final_sd_thresh` | Mean and standard deviation of the migration threshold at the end |
| `mean_population_size`, `mean_prop_hunters`, `mean_prop_free_riders`, `mean_energy` | Means of yearly values across the simulation |
| `mean_group_size_variance`, `mean_occupied_group_size`, `mean_total_migrations` | Means of yearly group and migration statistics |
| `mean_thresh_over_time`, `mean_sd_thresh_over_time` | Mean migration threshold statistics over time |
| `mean_age_all_deaths`, `mean_age_energy_death`, `mean_age_natural_death` | Mean age at death, overall and by cause |
| `mean_prop_reproducing` | Mean proportion of the population reproducing each year |
| `mean_w`, `sd_w` | Mean and standard deviation of hunting propensity at the end |

The amendment scripts record additional outputs relevant to each amendment (for example, the migration threshold over time in `Variable threshold amendment.R`).

### Reproducing Figure 2

Figure 2 shows histograms of each parameter's values for the simulations in which the population survived the full 1,500 years. To reproduce it from `FN 20k sims full results.xlsx`:

1. Define a binary outcome: `survived = 1` if `time_to_collapse > 547499` (population > 0 at the end of the simulation), otherwise `0`.
2. Plot histograms of each `*_val` column for the surviving simulations only (954 of the 20,000).
3. Report the point-biserial correlation *r* between each parameter and the binary survival outcome, computed across **all** 20,000 simulations. (Computing it on the survivors alone would give no variation in the outcome.)

---

## Summary of amendments

| Original model | Amended model | Result |
| --- | --- | --- |
| Migration threshold `thresh` is fixed across all agents and time | Threshold is variable and heritable | No impact on population survival |
| All agents pay a fixed cost `MC` when migrating | Locations on a lattice; `MC` varies with distance | No impact on population survival |
| Hunters have propensity 0.75-1, free-riders 0-0.1 | Hunters hunt (*w* = 1), free-riders do not (*w* = 0) | No impact on population survival |
| As above | Hunters 0.875, free-riders 0.05 | No impact on population survival |
| As above | Propensity to hunt is continuous and heritable | Minimal impact within 1,500 years; populations do not survive beyond c.15,000 years |
| Simulation begins with all hunters | Starting proportion of hunters varies from 0 to 1 | Starting population must be c.25% hunters to reliably survive |
| Hunters' offspring can mutate to free-riders; free-rider offspring are always free-riders | Mutation occurs equally in both directions | No impact on population survival |

---

## Notes on reproducibility

- Simulations are stochastic. `future_lapply(..., future.seed = TRUE)` is used for parallel-safe random number generation, but no global seed is set, so individual runs will differ. Aggregate results (for example, survival rates across hundreds or thousands of runs) should be consistent.
- The model uses a day as its timestep; results are recorded yearly.

---

## Citation

If you use this code or data, please cite the paper:

```
Erskine, S. & Dyble, M. [Year]. Testing the robustness of mobility as a mechanism
for the sustainability of demand sharing. [Journal]. [DOI]
```

and the original model:

```
Lewis, H.M., Vinicius, L., Strods, J., Mace, R. & Migliano, A.B. (2014).
High mobility explains demand sharing and enforced cooperation in egalitarian
hunter-gatherers. Nature Communications 5, 5789. doi:10.1038/ncomms6789
```

## Funding

S.E. acknowledges funding support from the Cambridge ESRC: Cam-DPT and the Harding Distinguished Postgraduate Scholars programme.

## Contact

Sarah Erskine, [@SCErskine](https://github.com/SCErskine)
sce32@cam.ac.uk

## Licence

[Add a licence here, e.g. MIT for the code. See https://choosealicense.com]
