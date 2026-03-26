# Climate Agreements Simulation

An interactive, browser-based tool for simulating international climate negotiations. Designed for classroom exercises where student teams represent different country groups and negotiate emission reduction commitments and financial transfers.

**Live tool:** Hosted via GitHub Pages at [jtbradt.github.io/climate-agreements-sim-expo](https://jtbradt.github.io/climate-agreements-sim-expo/)

## Overview

Students set three policy parameters for each of six country groups:

1. **Emission reduction start year** (2025--2100)
2. **Annual emission reduction rate** (0--20%)
3. **Net financial transfer** ($B/yr; positive = paying, negative = receiving)

The tool instantly computes and visualizes:
- GHG emission pathways under business-as-usual and negotiated policies
- Global temperature change relative to pre-industrial levels
- Country-specific climate damages (% change in per capita GDP)
- Cumulative abatement costs by country

## Country Groups

All economic values are in **2024 current (nominal) US dollars** unless otherwise noted. Base year for all projections is **2024**.

| Group | Representative Countries | Pop (M) | GDP/cap ($) | Baseline Temp (°C) | Tconv |
|---|---|---|---|---|---|
| **United States** | USA | 340 | 84,534 | 13.6 | 1.32 |
| **European Union** | EU-27 | 450 | 43,305 | 11.2 | 1.14 |
| **China** | CHN | 1,409 | 13,303 | 14.2 | 1.29 |
| **India** | IND | 1,451 | 2,695 | 25.7 | 1.14 |
| **Sub-Saharan Africa** | World Bank SSF aggregate | 1,291 | 1,533 | 22.7 | 1.12 |
| **Oil States** | SAU, ARE, KWT, IRQ, IRN | 189 | 14,339 | 19.8 | 1.42 |

Oil States GDP per capita is a population-weighted average of the five constituent countries.

## Data Sources

### GDP and Population

**Source:** World Bank World Development Indicators (WDI), 2024 data (last updated 2025-02-24).

- **GDP per capita:** Indicator `NY.GDP.PCAP.CD` — GDP per capita in current US dollars, 2024. The EU uses the World Bank's EU-27 aggregate (code `EUU`). Sub-Saharan Africa uses the World Bank regional aggregate (code `SSF`). Oil States is the population-weighted average of Saudi Arabia ($35,122), UAE ($50,274), Kuwait ($32,718), Iraq ($6,074), and Iran ($5,190).
- **Population:** Indicator `SP.POP.TOTL`, 2024.

### GDP Growth Rates

Long-run BAU growth rates are calibrated assumptions informed by IMF World Economic Outlook projections and OECD long-term baseline scenarios. These represent average annual per capita GDP growth rates projected forward from 2024, absent climate damages:

| Group | Growth Rate | Rationale |
|---|---|---|
| United States | 1.5%/yr | Mature economy, historical trend ~1.5-2% |
| European Union | 1.3%/yr | Mature economy, aging population |
| China | 3.5%/yr | Decelerating from ~5% current; convergence dynamics |
| India | 5.0%/yr | Rapid industrialization, demographic dividend |
| Sub-Saharan Africa | 2.5%/yr | Young population, low base, institutional constraints |
| Oil States | 1.5%/yr | Diversification efforts, volatile oil revenues |

### Baseline Temperature and Temperature Conversion Factors

**Source:** Burke, Hsiang & Miguel (2015) replication package.

> Burke, M., Hsiang, S. M., & Miguel, E. (2015). Global non-linear effect of temperature on economic production. *Nature*, 527(7577), 235--239.

- **Population-weighted mean temperature** — Annual average temperature weighted by sub-national population distribution, from the University of Delaware gridded temperature dataset (`UDel_temp_popweight` in `GrowthClimateDataset.csv`), averaged over 1980--2010. Multi-country groups use population-weighted averages.
- **Temperature conversion factors (Tconv)** — From `CountryTempChange_RCP85.csv`. These translate global mean temperature change into country-specific local temperature change, derived from the CMIP5 RCP8.5 multi-model ensemble mean. Computed as the ratio of population-weighted country-level warming to global mean warming. For example, Tconv = 1.42 for Oil States means they warm 42% faster than the global average.
- **Damage function coefficients** — The quadratic relationship between temperature and GDP growth from Burke et al.'s pooled regression (no lags): `growth_effect(T) = 0.0127*T - 0.0005*T²`.

### Emissions Data

**Source:** EDGAR 2025 Report (Emissions Database for Global Atmospheric Research, European Commission Joint Research Centre), IEA Global Energy Review 2025, and Climate Action Tracker country assessments.

All emissions are total GHG emissions in **Gt CO2e/year**, 2024 data.

| Group | Emissions (Gt CO2e) | BAU Growth Rate | Source Notes |
|---|---|---|---|
| United States | 5.9 | -0.5%/yr | Market-driven decarbonization, partially offset by policy uncertainty |
| European Union | 3.2 | -2.0%/yr | EU ETS + Fit for 55 regulations driving continued decline |
| China | 15.5 | +0.5%/yr | Near-peak; massive renewable deployment offsetting coal |
| India | 4.4 | +2.5%/yr | Strong GDP growth and industrialization driving demand |
| Sub-Saharan Africa | 2.4 | +2.5%/yr | Low base, rapid population growth and urbanization |
| Oil States | 2.8 | +1.5%/yr | Domestic consumption growth, petrochemical expansion |

Sub-Saharan Africa and Oil States aggregates are summed from individual country EDGAR entries. Emissions growth rates are calibrated BAU assumptions (no new policies beyond current implementation) informed by IEA and Climate Action Tracker projections.

### Atmospheric CO2 Concentration

- **Base concentration (2024):** 427.0 ppm, from NOAA Global Monitoring Laboratory.
- **Pre-industrial concentration:** 297.06 ppm, from Burke et al. (2015).
- **Emissions-to-concentration model:** OLS regression on historical data (2000--2024): `ΔCO2 = 0.9893 + 0.0245 × total_emissions_lag`. Applied iteratively from the 2024 base.

## Methodology

### Emission Reduction Pathways

For each country, business-as-usual emissions grow exponentially from the base year:

```
emissions_BAU(t) = emissions_0 × (1 + g)^(t - 2024)
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
T_local(t) = T_base + (T_global(t) - T_global(2024)) × Tconv
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

These parameters are pedagogical calibrations, not sourced from empirical marginal abatement cost curves.

### Discounting

All net present value (NPV) calculations in the Financial Summary use a **3% annual discount rate**, consistent with the US government's central rate for regulatory cost-benefit analysis. This applies to:

- **Cumulative abatement costs (NPV):** Annual abatement spending is discounted back to the 2024 base year before summing.
- **Avoided damages (NPV):** Cumulative dollar damages (per capita GDP loss × population, summed over all years) are discounted back to 2024, computed under both BAU and policy scenarios. Avoided damages = BAU damages − policy damages.

The abatement cost **chart** shows undiscounted cumulative costs for visual clarity. The **Financial Summary table** reports discounted NPV values.

Note: The Burke et al. damage function applies to GDP *growth rates*, which causes damages to compound over time. Combined with population growth, undiscounted cumulative damages can be very large. Discounting substantially reduces these figures but the magnitudes remain significant, particularly for hot, populous countries (India, Sub-Saharan Africa).

### Financial Transfers

Each country group can set a net annual transfer in $B/yr. Positive values indicate payments into a climate fund; negative values indicate receipts. The tool displays whether transfers are balanced (sum to zero) across all groups.

Transfers are tracked as a policy instrument for negotiation purposes. In the current version, they represent commitments that students must negotiate as part of the overall agreement.

### Financial Summary

The Financial Summary tab provides two views:

1. **Default view** (visible during negotiation): Shows NPV abatement costs, cumulative transfers, and climate damages as % GDP per capita in 2100. The unit mismatch between dollar costs and percentage damages is intentional — it prevents students from trivially computing net benefits, preserving negotiation tension.

2. **Detailed Breakdown** (hidden behind an instructor toggle): Reveals NPV avoided damages in $B and net benefit (avoided damages − abatement cost − transfers). Designed for use during the post-exercise debrief to illustrate the full cost-benefit picture and discuss why collective action was (or wasn't) worth it.

## Technical Details

The simulation runs entirely client-side as a single HTML file with no server dependencies. The only external resource is Chart.js (loaded from CDN). All data is embedded directly in JavaScript.

- **Framework:** Vanilla HTML/CSS/JavaScript
- **Charting:** [Chart.js](https://www.chartjs.org/) v4.4.7
- **Deployment:** GitHub Pages (static hosting from `docs/` directory)
- **Performance:** All computations complete in <1ms; chart updates are instantaneous on slider drag

## Repository Structure

```
├── docs/
│   ├── index.html          # Complete simulation app (single file)
│   └── slides.html         # Briefing slides for classroom projection (2 slides)
├── briefings/              # Confidential country briefing sheets (1 per team)
│   ├── united-states.md
│   ├── european-union.md
│   ├── china.md
│   ├── india.md
│   ├── sub-saharan-africa.md
│   └── oil-states.md
├── data/
│   ├── BurkeHsiangMiguel2015_Replication.zip   # BHM replication package
│   ├── co2_conc.csv         # Historical CO2 concentrations
│   ├── ghg_emissions_baseline.csv  # Legacy baseline emissions (original 3-group model)
│   ├── gdp_baseline.csv     # Legacy baseline GDP (original 3-group model)
│   └── pop_baseline.csv     # Legacy baseline population (original 3-group model)
├── FACILITATION_GUIDE.md    # Instructor guide for classroom exercise
└── README.md
```

The `data/` CSV files are from the original 3-group version of the simulation and are retained as reference. The current 6-group model derives its parameters from the sources documented above.

## Author

Jacob Bradt (jacob.bradt@mccombs.utexas.edu)

## References

- Burke, M., Hsiang, S. M., & Miguel, E. (2015). Global non-linear effect of temperature on economic production. *Nature*, 527(7577), 235--239.
- EDGAR — Emissions Database for Global Atmospheric Research (2025). *GHG Emissions of All World Countries*. European Commission, Joint Research Centre.
- IEA (2025). *Global Energy Review: CO2 Emissions*.
- Climate Action Tracker (2025). Country assessments. https://climateactiontracker.org/
- World Bank (2025). World Development Indicators. https://data.worldbank.org/
- NOAA Global Monitoring Laboratory (2025). Trends in Atmospheric Carbon Dioxide. https://gml.noaa.gov/ccgg/trends/
