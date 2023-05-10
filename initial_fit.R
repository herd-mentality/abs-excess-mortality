## Load libraries ----
library(tidyverse)
library(MASS)
library(janitor)
library(purrr)
library(gganimate)
library(transformr)
library(gifski)

select <- dplyr::select

abs_mortality_post_theme <- theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(),
    text             = element_text(size = 24)
  ) 

plot_dim <- list(width = 1684, height = 678)

## Ingest data ----
# Data sourced from: https://www.abs.gov.au/articles/measuring-australias-excess-mortality-during-covid-19-pandemic-doctor-certified-deaths
suppressWarnings(
  abs_deaths <- file.path(
    getwd(), 
    'input', 
    'Comparison of all cause baseline and COVID-19 period deaths against regression, January 2016 - February 2022 (a)(b)(c)(d)(e)(f).csv'
  ) %>% 
    read_csv(skip = 1, show_col_types = FALSE) %>% 
    separate_wider_delim(`95% bounds`, '|', names = c('95%_lower', '95%_higher')) %>% 
    mutate(
      across(starts_with('95%'), as.numeric),
      `Week starting date` = dmy(`Week starting date`),
      week_number = 1:n()
    ) %>% 
    filter(!is.na(`Week starting date`)) %>% 
    clean_names()
)

## Initial plot ----
abs_deaths %>%
  dplyr::select(week_starting_date, expected, observed) %>% 
  pivot_longer(expected:observed, names_to = 'series', values_to = 'value') %>%
  mutate(series = str_to_title(series)) %>% 
  ggplot(aes(x = week_starting_date, y = value, colour = series)) +
  geom_line() + 
  scale_y_continuous(labels = scales::comma) +
  annotate(
    'rect', 
    xmin = ymd('2020-01-01'), 
    xmax = cutoff_date, 
    ymin = -Inf, ymax = Inf, 
    alpha = 0.2, fill = '#aaaaaa'
  ) +
  annotate(
    'rect', 
    xmin = cutoff_date, 
    xmax = max(abs_deaths$week_starting_date), 
    ymin = -Inf, ymax = Inf, 
    alpha = 0.2, fill = '#FA9F42'
  ) +
  labs(
    x = 'Week starting date', 
    y = 'Mortality counts', 
    colour = element_blank(),
    title = 'ABS observed and expected mortality counts',
    caption = 'Forecast year 2021 highlighted'
  ) +
  abs_mortality_post_theme +
  theme(legend.position = 'bottom')

## Modelling parameters ----
cutoff_date <- ymd('2021-01-01')
train_set   <- abs_deaths %>% filter(week_starting_date >= cutoff_date %m-% years(5) & week_starting_date < cutoff_date)
test_set    <- abs_deaths %>% filter(week_starting_date >= cutoff_date & week_starting_date < cutoff_date %m+% years(1))

## Tests ----
# Fit robust regression using SAS default values
rlm_fit <- MASS::rlm(
  observed ~ week_number + 
    I(week_number^2) +
    I(sin(2 * pi * week_number / 52.18)) + 
    I(cos(2 * pi * week_number / 52.18)), 
  data = abs_deaths %>% mutate(week_number = 1:n()),
  method = 'M',  
  scale.est = 'MAD',
  psi = 'psi.bisquare',
  c = 4.685
)

# Plot weights
abs_deaths %>% 
  mutate(bisquare_weighting = rlm_fit$w) %>% 
  ggplot(aes(x = week_starting_date, y = observed, alpha = bisquare_weighting)) +
  geom_point()

## Predicting 2020 and 2021 ----

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
  rlm_params = list(
    method = NULL, scale.est = NULL, psi = NULL, c = NULL
  ),
  verbose = TRUE
) {
  
  # Parameter checks
  if (!is.numeric(year_to_predict)) stop('Please provide a numeric year value')
  if (is.null(date_col) | !is.character(date_col)) stop('Please provide the name of the date column in full_time_series')
  if (is.null(values_col) | !is.character(values_col)) stop('Please provide the name of the values column in full_time_series')
  
  # Base RLM parameter list - add any missing params to rlm_params
  base_rlm_params = list(
    method = NULL, scale.est = NULL, psi = NULL, c = NULL
  )
  rlm_params <- modifyList(base_rlm_params, rlm_params)

  # Initiate results list
  results <- list()
  
  # Get full list of years present in full_time_series
  all_years <- full_time_series %>% select(!!as.symbol(date_col)) %>% pull() %>% year() %>% unique()
  
  # Grab num_years_baseline years before the year to predict
  eligible_training_years <- all_years %>% Filter(function(x) x < year_to_predict, .)
  
  if (length(eligible_training_years) < num_years_baseline) 
    stop('Not enough years in historical data to train baseline of length num_years_baseline')
  
  training_years <- eligible_training_years[(length(eligible_training_years) - num_years_baseline + 1):length(eligible_training_years)]
  
  # Split full series into train and test based on provided parameters
  train <- full_time_series %>% 
    filter(year(!!as.symbol(date_col)) %in% training_years)
  
  test <- full_time_series %>% 
    filter(year(!!as.symbol(date_col)) %in% year_to_predict)
  
  # Train model 
  if (verbose) cat(str_interp(
    'Fitting RLM with params 
    | method    = ${if (is.null(rlm_params$method)) "M" else rlm_params$method}
    | scale.est = ${if (is.null(rlm_params$scale.est)) "MAD" else rlm_params$scale.est}
    | psi       = ${if (is.null(rlm_params$psi)) "psi.bisquare" else rlm_params$psi}
    | c         = ${if (is.null(rlm_params$c)) 4.685 else rlm_params$c}
    \n'
  ))
  
  results$fit_rlm <- MASS::rlm(
    as.formula(formula), 
    data      = train,
    method    = if (is.null(rlm_params$method)) 'M' else rlm_params$method,  
    scale.est = if (is.null(rlm_params$scale.est)) 'MAD' else rlm_params$scale.est,
    psi       = if (is.null(rlm_params$psi)) 'psi.bisquare' else rlm_params$psi,
    c         = if (is.null(rlm_params$c)) 4.685 else rlm_params$c
  )
  
  
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

# Usage
predict_year_rlm(
  year_to_predict    = 2021, 
  full_time_series   = abs_deaths, 
  'week_starting_date', 'observed', 
  num_years_baseline = 5,
  # We get closer to the ABS' values if we don't have the squared term for some reason
  formula = 'observed ~ week_number + 
    I(week_number^2) + 
    I(sin(2 * pi * week_number / 52.18)) + 
    I(cos(2 * pi * week_number / 52.18))'
)$results %>% 
  ggplot(aes(x = week_starting_date)) +
  geom_line(aes(y = observed, colour = 'Observed')) +
  geom_line(aes(y = expected, colour = 'ABS Predicted')) +
  geom_line(aes(y = predicted, colour = 'Herd Mentality Predicted')) +
  annotate(
    'rect', 
    xmin = cutoff_date, 
    xmax = max(abs_deaths$week_starting_date), 
    ymin = -Inf, ymax = Inf, 
    alpha = 0.2, fill = '#FA9F42'
  ) +
  scale_colour_manual(values = c(
    'Observed' = 'grey', 
    'ABS Predicted' = '#EC4899', 
    'Herd Mentality Predicted' = '#14B8A6'
  )) +
  labs(
    colour = element_blank()
  )

## Find ABS's parameters ----

# Grid search over values of psi, scale.est, and c to find the closest match
# to ABS's parameters

# It looks like removal of the squared time term brings us closer
possible_parameters <- bind_rows(
  crossing(
    c = seq(4.685 * 0.5, 4.685 * 1.5, 0.1),
    scale.est = c("MAD", "Huber"),
    psi = 'psi.bisquare'
  ),
  crossing(
    k = seq(1.345 * 0.5, 1.345 * 1.5, 0.1),
    scale.est = c("MAD", "Huber"),
    psi = 'psi.huber'
  ),
  crossing(
    a = seq(2 * 0.5, 2 * 1.5, 0.1),
    b = seq(4 * 0.5, 4 * 1.5, 0.1),
    c = seq(8 * 0.5, 8 * 1.5, 0.1),
    scale.est = c("MAD", "Huber"),
    psi = 'psi.hampel'
  )
)

suppressWarnings(
  possible_parameters_and_sse <- possible_parameters %>% 
    pmap_dfr(
      function(scale.est, psi, a, b, c, k) {
        # todo: update func to take in other tunng params as well
        prediction <- if (psi == 'psi.bisquare') {
          
          predict_year_rlm(
            year_to_predict    = 2021, 
            full_time_series   = abs_deaths, 
            'week_starting_date', 'observed', 
            num_years_baseline = 5,
            # We get closer to the ABS' values if we don't have the squared term for some reason
            formula = 'observed ~ week_number + 
              I(sin(2 * pi * week_number / 52.18)) + 
              I(cos(2 * pi * week_number / 52.18))',
            rlm_params = list(
              c         = c,
              scale.est = scale.est,
              psi       = psi
            ),
            verbose = FALSE
            )$results 
          
        } else {}
          
          
        sse <- prediction %>% 
          mutate(diff = (predicted - observed)^2) %>% 
          .$diff %>% 
          sum()
        
        tibble(
          c         = c,
          scale.est = scale.est,
          psi       = psi,
          sse       = sse
        )
        
      }
    )
)

possible_parameters_and_sse %>%
  group_by(psi) %>% 
  filter(sse == min(sse))



## Look at effect of varying c ----
varying_c <- seq(2.5, 15.5, 1) %>% 
  purrr::map_dfr(
    function(x) {
      
      model <- MASS::rlm(
        observed ~ week_number + 
          I(week_number^2) +
          I(sin(2 * pi * week_number / 52.18)) + 
          I(cos(2 * pi * week_number / 52.18)), 
        data = train_set,
        method = 'M',  
        scale.est = 'MAD',
        psi = 'psi.huber',
        c = x
      )
      
      train_set %>% 
        mutate(
          expected = predict(model, train_set),
          series = 'Historical'
        ) %>% 
        bind_cols(
          tibble(
            c = x,
            weights = model$w
          )
        ) %>% 
        bind_rows(
          test_set %>% 
            mutate(
              expected = predict(model, test_set), 
              series = 'Forecast',
              c = x, 
              weights = 1
            )
        )
      
    }
  )

animation_gif <- varying_c %>%
  ggplot(aes(x = week_starting_date, y = observed)) +
  geom_point(aes(alpha = weights)) +
  geom_line(aes(y = expected, colour = series), show.legend = FALSE) +
  scale_y_continuous(labels = scales::comma) +
  scale_colour_manual(values = c('Historical' = '#14B8A6', 'Forecast' = '#EC4899')) +
  annotate(
    'rect', 
    xmin = cutoff_date, 
    xmax = max(abs_deaths$week_starting_date), 
    ymin = -Inf, ymax = Inf, 
    alpha = 0.2, fill = '#FA9F42'
  ) +
  abs_mortality_post_theme +
  labs(
    x        = 'Week starting date',
    y        = 'Mortality counts',
    title    = 'Robust Regression',
    subtitle = 'Varying the outlier threshold, c',
    caption  = 'Forecast year 2021 highlighted, c = {formatC(frame_time, format = "f", digits = 2)}',
    alpha    = 'Weights'
  ) +
  transition_time(c, rev(range(varying_c$c))) +
  enter_fade() +
  exit_fade()

anim_save(
  file.path(getwd(), "robust_regression_varying_c_huber.gif"), 
  animate(
    animation_gif,
    width = plot_dim$width, height = plot_dim$height,
    # Duration of GIF is then nframes / fps
    fps = 10, nframes = 10 * 5
  )
)

## Find best RMSE/MAE/MAPE ----
seq(2.1, 5, 0.1) %>% 
  sapply(
    function(x) {
      
      model <- MASS::rlm(
        observed ~ week_number + 
          I(week_number^2) +
          I(sin(2 * pi * week_number / 52.18)) + 
          I(cos(2 * pi * week_number / 52.18)), 
        data = abs_deaths,
        method = 'M',  
        scale.est = 'MAD',
        psi = 'psi.bisquare',
        c = x
      )
      
      residuals <- abs_deaths$observed - predict(model, abs_deaths, type = 'response')

      sqrt(mean(residuals^2))
      
    }
  )

# AS c increases, RMSE over training set decreases. This is because a higher c means that the outlying points will be included when optimising the least squares loss - leads to a better RMSE (since RMSE is calculated over all points)
# TODO: 
#  - We want to look at forecasting over 2020, but also forecasting 2021 including 2020 in baseline period since that's topical
#  - Look at grid searching over all perms of c/scaleest/psi to find best RMSE using CV - then apply to 2021?



model <- MASS::rlm(
  observed ~ week_number +                 # ABS-modified Serfling formula
    I(week_number^2) + 
    I(sin(2 * pi * week_number / 52.18)) + 
    I(cos(2 * pi * week_number / 52.18)), 
  data      = train_set,
  method    = 'M',                         # M estimation
  scale.est = 'MAD',                       # Scale function, median absolute deviation
  psi       = 'psi.bisquare',              # Weight function, Tukey's bisquare   
  c         = 4.685                        # Outlier threshold, SAS default for bisquare
)


# TODO: plot CIs/PIs as well

train_set %>% 
  mutate(expected = predict(model, train_set, type = 'response'), series = 'train') %>% 
  bind_rows(
    test_set %>% 
      mutate(expected = predict(model, test_set, type = 'response'), series = 'test')
  ) %>% 
  dplyr::select(week_starting_date, series, expected, observed) %>% 
  pivot_longer(expected:observed, names_to = 'type', values_to = 'value') %>% 
  ggplot(aes(x = week_starting_date, y = value, colour = type)) +
  geom_point() +
  geom_vline(xintercept = cutoff_date)



seq(1.5, 5, 0.1) %>% 
  purrr::map_dfr(
    function(x) {
      
      model <- MASS::rlm(
        observed ~ week_number, 
        data = train_set,
        method = 'M',  
        scale.est = 'MAD',
        psi = 'psi.bisquare',
        c = x
      )
      
      residuals <- test_set$observed - predict(model, test_set, type = 'response')
      
      tibble(c = x, rmse = sqrt(mean(residuals^2)), mae = mean(abs(residuals)))
      
    }
  ) %>% 
  pivot_longer(rmse:mae, names_to = 'measure', values_to = 'value') %>% 
  ggplot(aes(x = c, y = value, colour = measure)) +
  geom_point()

## Animation ----


## Plot forecast ----


## Compare different scale estimators and psi functions ----
c('psi.huber', 'psi.hampel', 'psi.bisquare') %>% 
  map_dfr(
    function(x) {
      
      abs_deaths %>% 
        bind_cols(
          tibble(
            psi = x,
            weights = MASS::rlm(
              observed ~ week_number, 
              data = abs_deaths,
              method = 'M',  
              scale.est = 'MAD',
              psi = x
              # We dont set any tuning parameters here since the names change between the psis
            )$w
          )
        )
      
    }
  ) %>% 
  ggplot(aes(x = week_starting_date, y = observed, alpha = weights)) +
  geom_point(colour = '#14B8A6') +
  facet_grid(rows = vars(psi)) +
  scale_y_continuous(labels = scales::comma)

## Compare different scale estimators and psi functions ----
# TODO: Use gridExtra to place all two/three plots together 

c("MAD", "Huber", "proposal 2") %>% 
  map_dfr(
    function(x) {
      
      abs_deaths %>% 
        bind_cols(
          tibble(
            scaleest = x,
            weights = MASS::rlm(
              observed ~ week_number, 
              data = abs_deaths,
              method = 'M',  
              scale.est = x,
              psi = 'psi.bisquare',
              c = 1.5
            )$w
          )
        )
      
    }
  ) %>% 
  ggplot(aes(x = week_starting_date, y = observed, alpha = weights)) +
  geom_point(colour = '#14B8A6') +
  facet_grid(rows = vars(scaleest)) +
  scale_y_continuous(labels = scales::comma)

## Do grid search over all c, scaleest, psi - best forecast rmse is second model, plot against ABS' SAS defaults