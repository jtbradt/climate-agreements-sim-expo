# Climate Agreements Simulation

An interactive, browser-based tool for simulating international climate negotiations. Designed for classroom exercises where student teams represent different country groups and negotiate emission reduction commitments and financial transfers.

**Live tool:** Hosted via GitHub Pages at [jtbradt.github.io/climate-agreements-sim-expo](https://jtbradt.github.io/climate-agreements-sim-expo/)

## Overview

Students set three policy parameters for each of six country groups:

1. **Emission reduction start year** (2015--2100)
2. **Annual emission reduction rate** (0--20%)
3. **Net financial transfer** ($B/yr; positive = paying, negative = receiving)

The tool instantly computes and visualizes:
- GHG emission pathways under business-as-usual and negotiated policies
- Global temperature change relative to pre-industrial levels
- Country-specific climate damages (% change in per capita GDP)
- Cumulative abatement costs by country

## Country Groups

| Group | Representative Countries | Pop (M, 2010) | GDP/cap (2010 $) | Baseline Temp (°C) |
|---|---|---|---|---|
| **United States** | USA | 309 | 43,952 | 13.6 |
| **European Union** | DEU, FRA, GBR, ITA, ESP | 316 | 33,168 | 11.2 |
| **China** | CHN | 1,338 | 2,870 | 14.2 |
| **India** | IND | 1,206 | 1,032 | 25.7 |
| **Sub-Saharan Africa** | NGA, ETH, KEN, TZA, ZAF | 384 | 1,338 | 22.7 |
| **Oil States** | SAU, ARE, KWT, IRQ, IRN | 144 | 7,185 | 19.8 |

Multi-country groups use population-weighted averages for all parameters.

## Data Sources

### Burke, Hsiang & Miguel (2015) Replication Package

The primary data source is the replication package for:

> Burke, M., Hsiang, S. M., & Miguel, E. (2015). Global non-linear effect of temperature on economic production. *Nature*, 527(7577), 235--239.

From this package we extract:

- **Baseline GDP per capita** — Country-level GDP per capita from the World Development Indicators (WDI), averaged over 1980--2010. Stored in `GrowthClimateDataset.csv`.
- **Population** — Country populations circa 2010, from the same dataset.
- **Population-weighted mean temperature** — Annual average temperature weighted by sub-national population distribution, from the University of Delaware gridded temperature dataset. Variable `UDel_temp_popweight` in `GrowthClimateDataset.csv`, averaged over 1980--2010.
- **Historical GDP growth rates** — Mean annual per capita GDP growth from WDI (`growthWDI`), averaged 1980--2010.
- **Temperature conversion factors (Tconv)** — From `CountryTempChange_RCP85.csv`. These factors translate global mean temperature change into country-specific local temperature change, derived from the CMIP5 RCP8.5 multi-model ensemble mean. Computed as the ratio of population-weighted country-level warming to global mean warming. For example, Tconv = 1.42 for Oil States means they warm 42% faster than the global average.
- **Damage function coefficients** — The quadratic relationship between temperature and GDP growth from Burke et al.'s pooled regression: `growth_effect(T) = 0.0127*T - 0.0005*T²`.

### Emissions Data

Baseline GHG emissions (Gt CO2e/year) and business-as-usual growth rates are drawn from:

- **EDGAR 2025 Report** (Emissions Database for Global Atmospheric Research, European Commission Joint Research Centre) — 2024 country-level total GHG emissions.
- **IEA Global Energy Review 2025** — Regional emissions trends and growth rates.
- **Climate Action Tracker** — Country-level current-policy emissions projections.

| Group | Emissions (Gt CO2e, 2024) | BAU Growth Rate | Source Notes |
|---|---|---|---|
| United States | 5.9 | -0.5%/yr | Market-driven decarbonization, partially offset by policy uncertainty |
| European Union | 3.2 | -2.0%/yr | EU ETS + Fit for 55 regulations driving continued decline |
| China | 15.5 | +0.5%/yr | Near-peak; massive renewable deployment offsetting coal |
| India | 4.4 | +3.5%/yr | Strong GDP growth and industrialization driving demand |
| Sub-Saharan Africa | 2.4 | +3.0%/yr | Low base, rapid population growth and urbanization |
| Oil States | 2.8 | +2.0%/yr | Domestic consumption growth, petrochemical expansion |

Sub-Saharan Africa and Oil States aggregates are summed from individual country EDGAR entries.

### Emissions-to-Concentration Model

The relationship between global GHG emissions and atmospheric CO2 concentration is estimated via OLS regression on historical data (2000--2024):

```
ΔCO2_concentration = 0.9893 + 0.0245 × total_emissions_lag
```

This simple linear model is applied iteratively to accumulate atmospheric concentrations from a base of 390.42 ppm in 2010.

## Methodology

### Emission Reduction Pathways

For each country, business-as-usual emissions grow exponentially from the base year:

```
emissions_BAU(t) = emissions_0 × (1 + g)^(t - 2010)
```

When a reduction policy is active (starting at `reduction_year` with annual rate `r`):

```
emissions(t) = emissions_BAU(reduction_year) × (1 - r)^(t - reduction_year + 1)
```

### Temperature

Global atmospheric CO2 concentration is accumulated from total emissions across all six groups. Global mean temperature anomaly (relative to pre-industrial) is computed using:

```
T_global = λ × log₂(CO2 / CO2_preindustrial)
```

where `λ = 2.367` is the climate sensitivity parameter and `CO2_preindustrial = 297.06 ppm`.

### Country-Specific Warming

Each country experiences warming proportional to their temperature conversion factor:

```
T_local(t) = T_base + (T_global(t) - T_global(2010)) × Tconv
```

This reflects that high-latitude regions warm faster than the tropics, and continental interiors warm faster than coastal areas.

### Climate Damages (Burke et al. 2015)

The Burke et al. damage function estimates the effect of temperature on economic growth:

```
damage(T) = 0.0127 × T - 0.0005 × T²
```

This function peaks at approximately 13°C, meaning countries near this temperature (US, EU, China) experience smaller marginal damages, while hotter countries (India, SSA, Oil States) face larger growth penalties. Local temperature is capped at 30°C to avoid out-of-sample extrapolation, following Burke et al.

Annual GDP per capita growth under climate change is:

```
g_CC(t) = g_base + [damage(T_local(t)) - damage(T_base)]
```

Climate damages are reported as the percentage difference in per capita GDP between the climate-change and no-climate-change scenarios.

### Abatement Costs

Marginal abatement costs follow a quadratic schedule:

```
MAC(t) = c₁ × A(t) + c₂ × A(t)²
```

where `A(t) = 1 - (1-r)^(years_since_start)` is the cumulative abatement fraction. The parameters `c₁` and `c₂` vary by country to reflect differences in decarbonization costs:

| Group | c₁ | c₂ | Rationale |
|---|---|---|---|
| United States | 200 | 50 | High per-capita emissions, gas infrastructure, sprawl |
| European Union | 200 | 50 | Expensive labor, legacy infrastructure, diminishing returns |
| China | 100 | 25 | Scale economies, manufacturing capacity, but large coal fleet |
| India | 80 | 20 | Low labor costs, greenfield potential, excellent solar |
| Sub-Saharan Africa | 80 | 15 | Good renewable resources, but financing constraints |
| Oil States | 250 | 75 | Structural fossil fuel dependence on both supply and demand side |

### Financial Transfers

Each country group can set a net annual transfer in $B/yr. Positive values indicate payments into a climate fund; negative values indicate receipts. The tool displays whether transfers are balanced (sum to zero) across all groups.

Transfers are tracked as a policy instrument for negotiation purposes. In the current version, they represent commitments that students must negotiate as part of the overall agreement.

## Technical Details

The simulation runs entirely client-side as a single HTML file with no server dependencies. The only external resource is Chart.js (loaded from CDN). All data is embedded directly in JavaScript.

- **Framework:** Vanilla HTML/CSS/JavaScript
- **Charting:** [Chart.js](https://www.chartjs.org/) v4.4.7
- **Deployment:** GitHub Pages (static hosting from `docs/` directory)
- **Performance:** All computations complete in <1ms; chart updates are instantaneous on slider drag

## Repository Structure

```
├── docs/
│   └── index.html          # Complete simulation app (single file)
├── data/
│   ├── BurkeHsiangMiguel2015_Replication.zip   # BHM replication package
│   ├── co2_conc.csv         # Historical CO2 concentrations
│   ├── ghg_emissions_baseline.csv  # Legacy baseline emissions (original 3-group model)
│   ├── gdp_baseline.csv     # Legacy baseline GDP (original 3-group model)
│   └── pop_baseline.csv     # Legacy baseline population (original 3-group model)
└── README.md
```

The `data/` CSV files are from the original 3-group version of the simulation and are retained as reference. The current 6-group model derives its parameters from the BHM replication package and external emissions sources as documented above.

## Author

Jacob Bradt (jacob.bradt@mccombs.utexas.edu)

## References

- Burke, M., Hsiang, S. M., & Miguel, E. (2015). Global non-linear effect of temperature on economic production. *Nature*, 527(7577), 235--239.
- EDGAR — Emissions Database for Global Atmospheric Research (2025). *GHG Emissions of All World Countries*. European Commission, Joint Research Centre.
- IEA (2025). *Global Energy Review: CO2 Emissions*.
- Climate Action Tracker (2025). Country assessments. https://climateactiontracker.org/
