## Load libraries, define common parameters ----
library(tidyverse)
library(patchwork)
library(ggtext)
library(MASS)
library(janitor)
library(purrr)
library(furrr)
library(gganimate)
library(transformr)
library(gifski)

select <- dplyr::select

abs_mortality_post_theme <- theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # text             = element_text(size = 24),
    axis.line        = element_line()
  ) 

plot_dim <- list(width = 1684, height = 1684 * 9/16)

## Function definitions ----

predict_year_rlm <- function(
    year_to_predict, 
    full_time_series, 
    date_col = NULL,
    values_col = NULL,
    num_years_baseline,
    formula = 'observed ~ week_number + 
    I(week_number^2) + 
    I(sin(2 * pi * week_number / 52.18)) + 
    I(cos(2 * pi * week_number / 52.18))',
    # rlm_params = list(
    #   method = NULL, scale.est = NULL, psi = NULL, c = NULL
    # ),
    verbose = TRUE,
    ...
) {
  
  # Parameter checks
  if (!is.numeric(year_to_predict)) stop('Please provide a numeric year value')
  if (is.null(date_col) | !is.character(date_col)) stop('Please provide the name of the date column in full_time_series')
  if (is.null(values_col) | !is.character(values_col)) stop('Please provide the name of the values column in full_time_series')
  
  # # Base RLM parameter list - add any missing params to rlm_params
  # base_rlm_params = list(
  #   method = NULL, scale.est = NULL, psi = NULL, c = NULL
  # )
  # rlm_params <- modifyList(base_rlm_params, rlm_params)
  
  # Initiate results list
  results <- list()
  
  # Get full list of years present in full_time_series
  all_years <- full_time_series %>% select(!!as.symbol(date_col)) %>% pull() %>% year() %>% unique()
  
  # Grab num_years_baseline years before the year to predict
  eligible_training_years <- all_years %>% Filter(function(x) x < year_to_predict, .)
  
  if (length(eligible_training_years) < num_years_baseline) 
    stop('Not enough years in historical data to train baseline of length num_years_baseline')
  
  training_years <- eligible_training_years[
    (length(eligible_training_years) - num_years_baseline + 1):length(eligible_training_years)
  ]
  
  # Split full series into train and test based on provided parameters
  train <- full_time_series %>% 
    filter(year(!!as.symbol(date_col)) %in% training_years)
  
  test <- full_time_series %>% 
    filter(year(!!as.symbol(date_col)) %in% year_to_predict)
  
  # Train model 
  # if (verbose) cat(str_interp(
  #   'Fitting RLM with params 
  #   | method    = ${if (is.null(rlm_params$method)) "M" else rlm_params$method}
  #   | scale.est = ${if (is.null(rlm_params$scale.est)) "MAD" else rlm_params$scale.est}
  #   | psi       = ${if (is.null(rlm_params$psi)) "psi.bisquare" else rlm_params$psi}
  #   | c         = ${if (is.null(rlm_params$c)) 4.685 else rlm_params$c}
  #   \n'
  # ))
  
  results$fit_rlm <- MASS::rlm(
    as.formula(formula), 
    data = train,
    ...
    # method    = if (is.null(rlm_params$method)) 'M' else rlm_params$method,  
    # scale.est = if (is.null(rlm_params$scale.est)) 'MAD' else rlm_params$scale.est,
    # psi       = if (is.null(rlm_params$psi)) 'psi.bisquare' else rlm_params$psi,
    # c         = if (is.null(rlm_params$c)) 4.685 else rlm_params$c
  )

  results$converged <- ifelse('converged' %in% names(results$fit_rlm), results$fit_rlm$converged, FALSE)
  
  # Predict over the test set and return
  if (verbose) cat('Fitting complete, returning outputs\n')
  results$results <- bind_rows(
    train %>% 
      bind_cols(predicted = predict(results$fit_rlm, train)) %>% 
      mutate(series = 'Training'),
    test %>% 
      bind_cols(predicted = predict(results$fit_rlm, test)) %>% 
      mutate(series = 'Testing')
  )
  
  return(results)
  
}
