############################## FRONT MATTER ##############################
#   SCRIPT:     climate_agreements.R
#   AUTHOR:     Jacob Bradt (jacob.bradt@mccombs.utexas.edu)

# Load dependencies:
pacman::p_load(data.table, tidyverse, shiny, fixest, bslib, shinylive, here)

# Set project root:
i_am("app.R")

# Set project root:
local_dir <- "C:/Users/jtb3862/Documents/Teaching/Global Net Zero/Spring 2025/Simulations/Climate Agreements/"

# STEP 1: IMPORT BASELINE DATA -----

# Import CO2 concentrations under baseline:
conc <- fread(here("data/co2_conc.csv"))

# CO2 emissions:
emissions <- fread(here("data/ghg_emissions_baseline.csv"))

# Get relationship between total emissions and concentrations:
emissions[, total := Developed + `Developing A` + `Developing B`]
conc <- merge(emissions, conc)
conc[, delta_conc := Baseline - shift(Baseline)]
conc[, emissions_lag := shift(total)]
conc_fit <- feols(delta_conc ~ emissions_lag, data = conc)
rm(conc, emissions)

# CO2 emissions for pathways:
emissions <- fread(here("data/ghg_emissions_baseline.csv"))
emissions <- melt(emissions, id.vars = "Year", variable.name = "country", value.name = "emissions")

# Concentrations for pathways:
conc <- fread(here("data/co2_conc.csv"))
conc_base <- conc[Year == 2000]$Baseline
conc_preindustrial <- 297.05671
climate_sensitivity <- 2.36701601
conc_fit$coefficients
rm(conc)

# STEP 2: FUNCTIONS TO CHANGE EMISSIONS PATHWAYS -----

# Emissions pathways function:
emissions_pathway <- function(selected_country, reductions_year, reductions_rate){
  
  # Subset to target country:
  emissions_temp <- emissions[country == selected_country]
  
  # Adjust peak emissions between peak_year, reductions_year:
  peak_emissions <- emissions_temp[Year == reductions_year]$emissions
  
  # Adjust reductions between reductions_year, 2100 based on input reductions_rate
  for (y in reductions_year:2100) {
    emissions_temp[Year == y]$emissions <- peak_emissions * (1-reductions_rate)^(y-reductions_year + 1)
  }
  
  # Return new data:
  return(emissions_temp)
  
}

# Full emissions pathways function:
all_emissions_pathways <- function(countries, reduction_years, reductions_rates){
  
  new_emissions <- lapply(1:length(countries), function(c){
    emissions_pathway(countries[c], reduction_years[c], reductions_rates[c])
  }) %>%
    rbindlist(.)
  
  return(new_emissions)
  
}

# STEP 3: FUNCTIONS TO CALCULATE TEMPERATURE CHANGE -----

# Function to translate total emissions into delta T:
emissions_to_temp <- function(total_emissions){
  
  # Initiate atmos_conc field:
  total_emissions$atmos_conc <- 0
  
  # Calculate atmospheric concentrations:
  for (y in 1:nrow(total_emissions)){
    if (y == 1) {
      total_emissions[y]$atmos_conc <- conc_base +  total_emissions[y]$total_emissions * conc_fit$coefficients[2] + conc_fit$coefficients[1]
    } else{
      total_emissions[y]$atmos_conc <- total_emissions[y-1]$atmos_conc +  total_emissions[y]$total_emissions * conc_fit$coefficients[2]+ conc_fit$coefficients[1]
      
    }
  }
  
  # Calculate temperature change
  total_emissions$temp_change <- climate_sensitivity * log(total_emissions$atmos_conc/conc_preindustrial, base = 2)
  
  return(total_emissions)
}

# Calculate baseline temperature change:
temp_base <- emissions_to_temp(emissions[, .(total_emissions = sum(emissions)), by = Year])

# STEP 4: FUNCTIONS TO CALCULATE DAMAGES/ADAPTATION COSTS -----

# Import GDP:
gdp <- fread(here("data/gdp_baseline.csv"))
gdp <- melt(gdp, id.vars = "Year", variable.name = "country", value.name = "gdp")

# Burke et al. (2018) damage function:
damage_func <- function(t){
  0.0127 * t - 0.0005 * (t^2)
}

# Baseline temps (2000):
temp_base_developed <- 15
temp_base_developingA <- 24.8
temp_base_developingB <- 28.0

# Baseline growth rates:
gdp_growth_developed <- 0.015
gdp_growth_developingA <- 0.04
gdp_growth_developingB <- 0.01

# Function to calculate damages for a given country, temperature pathway
temp_to_damages <- function(selected_country, gdp, temp_change) {
  
  # Subset to selected country:
  gdp_temp <- gdp[country == selected_country]
  
  # Set baseline temp:
  if (selected_country == "Developed") {
    temp_base_selected <- temp_base_developed
  }
  if (selected_country == "Developing A"){
    temp_base_selected <- temp_base_developingA
  } 
  if (selected_country == "Developing B") {
    temp_base_selected <- temp_base_developingB
  }
  
  # Set no-climate effect growth rate:
  if (selected_country == "Developed") {
    gdp_growth_selected <- gdp_growth_developed
  }
  if (selected_country == "Developing A"){
    gdp_growth_selected <- gdp_growth_developingA
  }
  if (selected_country == "Developing B") {
    gdp_growth_selected <- gdp_growth_developingB
  }
  
  # Normalize temperature changes to base year 2000:
  temp_change_2000 <- temp_change$temp_change - temp_change[Year == 2000]$temp_change
  
  # Calculate annual climate damages:
  temp_damages <- damage_func(temp_base_selected + temp_change_2000) - damage_func(temp_base_selected)
  
  # Calculate GDP trajectory:
  gdp_output <- data.table(Year = 2000:2100, country = selected_country, gdp = 0)
  
  for (y in 1:nrow(gdp_output)) {
    if (y == 1) {
      gdp_output[y]$gdp <- gdp_temp[Year == 2000]$gdp
    } else {
      gdp_output[y]$gdp <- gdp_temp[y-1]$gdp * (1 + gdp_growth_selected + temp_damages[y])
    }
    
  }
  
  return(gdp_output)
  
}

# Function to calculate abatement costs:
abatement_costs <- function(selected_country, reduction_year, reduction_rate, cost_param1, cost_param2) {
  
  # Construct abatement cost DT:
  abatement_cost <- data.table(Year = 2000:2100, country = selected_country, cumulative_abatement = 0.0)
  
  # Set cumulative abatement:
  abatement_cost[ Year >= reduction_year]$cumulative_abatement <- sapply(1:(2101 - reduction_year), function(x){1 - (1 - reduction_rate)^x})
  
  # Calculate marginal abatement costs over time:
  abatement_cost[, mac := cost_param1 * cumulative_abatement + cost_param2 * cumulative_abatement^2 ]
  
  # Calculate actual abatement:
  emissions_BAU <- emissions_pathway(selected_country, 2100, 0.0)
  emissions_policy <- emissions_pathway(selected_country, reduction_year, reduction_rate)
  abatement_cost[, abated := c(0.0, diff(emissions_BAU$emissions - emissions_policy$emissions))]
  
  # Calculate abatement costs ($B):
  abatement_cost[, abatement_cost := mac * abated]
  
  return(abatement_cost)
  
}

# Import population data:
pop <- fread(here("data/pop_baseline.csv"))
pop <- melt(pop, id.vars = "Year", variable.name = "country", value.name = "pop")

# STEP 4: SHINY STUFF -----

# Define UI for application
ui <- page_fillable(
  
  
  layout_columns(
    card(card_header("GHG Emissions Pathways"),
         plotOutput("emissionsPath")),
    card(card_header("Temperature Change"),
         plotOutput("tempPath"))
  ),
  layout_columns(
    card(card_header("Policy Inputs"),
         layout_columns(
           card(
             numericInput(inputId = "developed_reduc",
                          label = "Developed Reduction Year",
                          value = 2100,
                          min = 2025,
                          max = 2100,
                          step = 1),
             numericInput(inputId = "developingA_reduc",
                          label = "Developing A Reduction Year",
                          value = 2100,
                          min = 2025,
                          max = 2100,
                          step = 1),
             numericInput(inputId = "developingB_reduc",
                          label = "Developing B Reduction Year",
                          value = 2100,
                          min = 2025,
                          max = 2100,
                          step = 1)
           ),
           card(
             numericInput(inputId = "developed_reduc_rate",
                          label = "Developed Reduction Rate",
                          value = 0,
                          min = 0,
                          max = 100,
                          step = 1),
             numericInput(inputId = "developingA_reduc_rate",
                          label = "Developing A Reduction Rate",
                          value = 0,
                          min = 0,
                          max = 100,
                          step = 1),
             numericInput(inputId = "developingB_reduc_rate",
                          label = "Developing B Reduction Rate",
                          value = 0,
                          min = 0,
                          max = 100,
                          step = 1)
           ), col_widths = c(6,6))),
    card(layout_columns(
      card(card_header("Climate Damages from Post-2000 Warming"), 
           plotOutput("gdpDamages")),
      card(card_header("Cumulative Abatement Costs"), 
           plotOutput("gdpMitigation"))
    ),
    col_widths = c(4,4)
    ), col_widths = c(6,6) )
  
)


# Define server logic required to draw lines:
server <- function(input, output) {
  
  # Function to plot new emissions pathways data:
  draw_emissions_paths <- function(list_reduc, list_reduc_rate){
    
    # New emissions pathway:
    emissions_new <- all_emissions_pathways(
      c("Developed", "Developing A", "Developing B"),
      list_reduc,
      list_reduc_rate
    )
    
    # Combine emissions data:
    emissions_new[, scenario := "New Policies"]
    emissions_old <- emissions
    emissions_old[, scenario := "Business-as-Usual"]
    emissions_full <- rbindlist(list(emissions_new, emissions_old))
    
    # Plot:
    ggplot(data = emissions_full) +
      geom_line(aes(x = Year, y = emissions, color = country, linetype = scenario), linewidth = 2) +
      theme_classic() +
      xlab("Year") +
      ylab("Gigatons CO2 Equivalent/year")
  }
  
  # Function to plot new temperature change path:
  draw_temperature_path <- function(list_reduc, list_reduc_rate){
    
    # New emissions pathway:
    emissions_new <- all_emissions_pathways(
      c("Developed", "Developing A", "Developing B"),
      list_reduc,
      list_reduc_rate
    )
    
    # Convert to temperature change:
    temp_new <- emissions_to_temp(emissions_new[, .(total_emissions = sum(emissions)), by = Year])
    
    # Combine temperature data:
    temp_new[, scenario := "New Policies"]
    temp_old <- temp_base
    temp_old[, scenario := "Business-as-Usual"]
    temp_full <- rbindlist(list(temp_new, temp_old))
    
    # Plot:
    ggplot(data = temp_full) +
      geom_hline(yintercept = 1.5, linetype = "dotted", color = "black") +
      geom_hline(yintercept = 2, linetype = "dotted", color = "black") +
      geom_line(aes(x = Year, y = temp_change, color = scenario, linetype = scenario), linewidth = 2) +
      scale_color_manual(values = c("black", "blue")) + 
      theme_classic() +
      xlab("Year") +
      ylab("Degrees Celsius Above Pre-industrial Levels")
  }
  
  # Function to plot GDP impacts:
  draw_climate_damages <- function(list_reduc, list_reduc_rate){
    
    # New emissions pathway:
    emissions_new <- all_emissions_pathways(
      c("Developed", "Developing A", "Developing B"),
      list_reduc,
      list_reduc_rate
    )
    
    # Convert to temperature change:
    temp_new <- emissions_to_temp(emissions_new[, .(total_emissions = sum(emissions)), by = Year])
    
    # Construct damages under no further climate change:
    temp_nocc <- data.table(Year = 2000:2100, temp_change = 0.0)
    climate_damages_nocc <- rbindlist(list(temp_to_damages("Developed", gdp, temp_nocc), temp_to_damages("Developing A", gdp, temp_nocc), temp_to_damages("Developing B", gdp, temp_nocc)))
    
    # Construct damages under BAU:
    climate_damages_BAU <- rbindlist(list(temp_to_damages("Developed", gdp, temp_base), temp_to_damages("Developing A", gdp, temp_base), temp_to_damages("Developing B", gdp, temp_base)))
    climate_damages_BAU[, gdp_delta := (climate_damages_BAU$gdp - climate_damages_nocc$gdp) / climate_damages_nocc$gdp]
    climate_damages_BAU[, scenario := "Business-as-Usual"]
    
    # Construct damages under policies:
    climate_damages_policies <- rbindlist(list(temp_to_damages("Developed", gdp, temp_new), temp_to_damages("Developing A", gdp, temp_new), temp_to_damages("Developing B", gdp, temp_new)))
    climate_damages_policies[, gdp_delta := (climate_damages_policies$gdp - climate_damages_nocc$gdp) / climate_damages_nocc$gdp]
    climate_damages_policies[, scenario := "New Policies"]
    
    # Combine BAU, policies:
    climate_damages <- rbindlist(list(climate_damages_BAU, climate_damages_policies))
    
    # Plot:
    ggplot(data = climate_damages) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
      geom_line(aes(x = Year, y = gdp_delta, color = country, linetype = scenario), linewidth = 2) +
      theme_classic() +
      scale_y_continuous(labels = scales::percent) + 
      xlab("Year") +
      ylab("% Change in per Capita GDP")
    
  }
  
  # Function to plot abatement costs curves:
  draw_abatement_costs <- function(list_reduc, list_reduc_rate){
    
    # Calculate abatement costs:
    all_abatement_costs <- rbindlist(list(
      abatement_costs("Developed", list_reduc[1], list_reduc_rate[1], 200, 50),
      abatement_costs("Developing A", list_reduc[2], list_reduc_rate[2], 100, 25),
      abatement_costs("Developing B", list_reduc[3], list_reduc_rate[3], 100, 10)
    ))
    
    # Calculate cumulative abatement spending:
    all_abatement_costs[, abatement_cost_cum := cumsum(abatement_cost), by = "country"]
    
    # Plot:
    ggplot(data = all_abatement_costs) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
      geom_line(aes(x = Year, y = abatement_cost_cum, color = country), linewidth = 2) +
      theme_classic() +
      scale_y_continuous(labels = scales::dollar) + 
      xlab("Year") +
      ylab("Cumulative Abatement Costs ($B)")
  }
  
  
  # Make plot of bill depth vs. bill length
  output$emissionsPath <- renderPlot({
    draw_emissions_paths(c(input$developed_reduc, input$developingA_reduc, input$developingB_reduc),
                         c(input$developed_reduc_rate/100, input$developingA_reduc_rate/100, input$developingB_reduc_rate/100))
  })
  output$tempPath <- renderPlot({
    draw_temperature_path(c(input$developed_reduc, input$developingA_reduc, input$developingB_reduc),
                          c(input$developed_reduc_rate/100, input$developingA_reduc_rate/100, input$developingB_reduc_rate/100))
  })
  output$gdpDamages <- renderPlot({
    draw_climate_damages(c(input$developed_reduc, input$developingA_reduc, input$developingB_reduc),
                         c(input$developed_reduc_rate/100, input$developingA_reduc_rate/100, input$developingB_reduc_rate/100))
  })
  output$gdpMitigation <- renderPlot({
    draw_abatement_costs(c(input$developed_reduc, input$developingA_reduc, input$developingB_reduc),
                         c(input$developed_reduc_rate/100, input$developingA_reduc_rate/100, input$developingB_reduc_rate/100))
  })
  
}


shinyApp(ui = ui, server = server)

