# TODO:
#  - Create/format graphs using the defined theme in libs-utils where possible

## Load libs and function definitions ----
source(file.path(getwd(), 'libs-utils.R'))


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
    caption = 'Target years 2020 & 2021 highlighted'
  ) +
  abs_mortality_post_theme +
  theme(legend.position = 'bottom')

ggsave(
  file.path(getwd(), 'plots', "abs_initial_plot.jpg"),
  device = 'jpg',
  width = plot_dim$width,
  height = plot_dim$height,
  units = 'px'
)

# Excess mortality
abs_deaths %>% 
  filter(year(week_starting_date) %in% 2020:2021) %>% 
  mutate(expected_mortality = observed - expected) %>% 
  group_by(year = year(week_starting_date)) %>% 
  summarise(expected_mortality = sum(expected_mortality), .groups = 'drop')

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

plot_predicted <- function(df, legend = 'none', ...) {
  
  df %>% 
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
    scale_y_continuous(labels = scales::comma) +
    labs(
      colour = element_blank(),
      ...
    ) +
    theme(legend.position = legend)
  
}

with_t2 <- predict_year_rlm(
  year_to_predict    = 2021, 
  full_time_series   = abs_deaths, 
  date_col           = 'week_starting_date', 
  values_col         = 'observed', 
  num_years_baseline = 5,
  # We get closer to the ABS' values if we don't have the squared term for some reason
  formula = paste0(
    'observed ~ week_number + ',
    'I(week_number^2) + ',
    'I(sin(2 * pi * week_number / 52.18)) + ',
    'I(cos(2 * pi * week_number / 52.18))'
  ),
  method             = 'M',
  scale.est          = 'MAD',
  psi                = 'psi.bisquare',
  c                  = 4.685
)$results %>% 
  plot_predicted(
    legend = 'none', 
    y = 'Observed', 
    x = element_blank(),
    title = 'Herd Mentality model with/out t^2 term'
  )

no_t2 <- predict_year_rlm(
  year_to_predict    = 2021, 
  full_time_series   = abs_deaths, 
  date_col           = 'week_starting_date', 
  values_col         = 'observed', 
  num_years_baseline = 5,
  # We get closer to the ABS' values if we don't have the squared term for some reason
  formula = paste0(
    'observed ~ week_number + ',
    # 'I(week_number^2) + ',
    'I(sin(2 * pi * week_number / 52.18)) + ',
    'I(cos(2 * pi * week_number / 52.18))'
  ),
  method             = 'M',
  scale.est          = 'MAD',
  psi                = 'psi.bisquare',
  c                  = 4.685
)$results %>% 
plot_predicted(legend = 'bottom', x = 'Week', y = 'Observed') 

with_t2 + no_t2 + plot_layout(ncol = 1)

ggsave(
  file.path(getwd(), 'plots', "abs_plot_w_wo_t2.jpg"),
  device = 'jpg',
  width = plot_dim$width,
  height = plot_dim$height * 2,
  units = 'px'
)

## Find ABS's parameters ----

# Grid search over values of psi, scale.est, and c to find the closest match
# to ABS's parameters by optimising SSE against their outputs.

# It looks like removal of the squared time term brings us closer so we'll try that.

step_size   <- 0.05
param_range <- 0.7
scale.ests  <- c("MAD", "Huber")
formulae    <- c(
  paste0(
    'observed ~ week_number + ',
    'I(sin(2 * pi * week_number / 52.18)) + ',
    'I(cos(2 * pi * week_number / 52.18))'
  ),
  paste0(
    'observed ~ week_number + ',
    'I(week_number^2) + ',
    'I(sin(2 * pi * week_number / 52.18)) + ',
    'I(cos(2 * pi * week_number / 52.18))'
  )
)

possible_parameters <- bind_rows(
  # Bisquare
  crossing(
    c = seq(4.685 * (1 - param_range), 4.685 * (1 + param_range), step_size),
    scale.est = scale.ests,
    formula = formulae,
    psi = 'psi.bisquare'
  ),
  # Huber
  crossing(
    k = seq(1.345 * (1 - param_range), 1.345 * (1 + param_range), step_size),
    scale.est = scale.ests,
    formula = formulae,
    psi = 'psi.huber'
  )#,
  # # Hampel
  # crossing(
  #   a = seq(2 * 0.5, 2 * 1.5, 0.1),
  #   b = seq(4 * 0.5, 4 * 1.5, 0.1),
  #   c = seq(8 * 0.5, 8 * 1.5, 0.1),
  #   scale.est = c("MAD", "Huber"),
  #   psi = 'psi.hampel'
  # )
)

future::plan(multisession, workers = 8)

possible_parameters_and_sse <- possible_parameters %>% 
  future_pmap_dfr(
    function(formula, psi, scale.est, c, k) {
      
      model_rlm <- if (psi == 'psi.bisquare') {
        
        predict_year_rlm(
          year_to_predict    = 2021, 
          full_time_series   = abs_deaths, 
          date_col           = 'week_starting_date', 
          values_col         = 'observed', 
          num_years_baseline = 5,
          formula            = formula,
          verbose            = FALSE,
          method             = 'M',
          scale.est          = scale.est,
          psi                = 'psi.bisquare',
          c                  = c
        )
        
      } else if (psi == 'psi.huber') {
        
        predict_year_rlm(
          year_to_predict    = 2021, 
          full_time_series   = abs_deaths, 
          date_col           = 'week_starting_date', 
          values_col         = 'observed', 
          num_years_baseline = 5,
          formula            = formula,
          verbose            = FALSE,
          method             = 'M',
          scale.est          = scale.est,
          psi                = 'psi.huber',
          k                  = k
        )
        
      }
      
      delta <- model_rlm$results %>% mutate(
        delta = (observed - predicted)^2
      )
      
      tibble(
        scale.est = scale.est,
        psi       = psi,
        formula   = formula,
        c         = c,
        k         = k,
        rmse      = sqrt(mean(delta$delta)),
        converged = if ('converged' %in% names(model_rlm)) model_rlm$converged else FALSE
      )
      
    }
  )

future::plan(sequential)

# How does this look compared to the ABS's model?
abs_closest_candidate <- possible_parameters_and_sse %>% 
  filter(psi == 'psi.bisquare') %>% 
  filter(rmse == min(rmse), converged) %>% 
  filter(row_number() == 1)

abs_closest_candidate_results <- predict_year_rlm(
  year_to_predict    = 2021, 
  full_time_series   = abs_deaths, 
  date_col           = 'week_starting_date', 
  values_col         = 'observed', 
  num_years_baseline = 5,
  formula            = abs_closest_candidate$formula,
  verbose            = FALSE,
  method             = 'M',
  scale.est          = abs_closest_candidate$scale.est,
  psi                = abs_closest_candidate$psi,
  c                  = abs_closest_candidate$c
)

abs_closest_candidate_results$results %>% 
  ggplot(aes(x = week_starting_date)) +
  geom_line(aes(y = observed, colour = 'Observed')) +
  geom_line(aes(y = expected, colour = 'ABS Predicted')) +
  geom_line(aes(y = predicted, colour = 'Herd Mentality Predicted')) +
  scale_colour_manual(values = c(
    'Observed' = 'grey', 
    'ABS Predicted' = '#EC4899', 
    'Herd Mentality Predicted' = '#14B8A6'
  )) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = 'Week', y = 'Deaths',
    colour  = element_blank(),
    title   = 'Closest parameters to ABS model',
    caption = paste0(
      'Closest candidate model -\n',
      str_interp('  Scale estimator: ${abs_closest_candidate$scale.est}\n'),
      str_interp('  Weighting function: ${abs_closest_candidate$psi}\n'),
      str_interp('  Outlier threshold: ${ifelse(is.na(abs_closest_candidate$k), abs_closest_candidate$c, abs_closest_candidate$k)}')
    )
  ) +
  theme(
    plot.caption = element_text(hjust = 0), 
    legend.position = 'bottom'
  )

abs_closest_candidate_results$results %>% 
  filter(year(week_starting_date) == 2021) %>% 
  mutate(
    abs_em = observed - expected,
    hm_em  = observed - predicted
  ) %>% 
  summarise(across(ends_with('_em'), sum))


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
  file.path(getwd(), 'plots', "robust_regression_varying_c_huber.gif"), 
  animate(
    animation_gif,
    width = plot_dim$width, height = plot_dim$height,
    # Duration of GIF is then nframes / fps
    fps = 10, nframes = 10 * 5
  )
)


# ## Grid search for optimal RMSE ----
# 
# # Use the leftover parameter data-frame for parameter combinations.
# # Perform leave-one-year-out-CV.
# 
# train_set_yearly_lst <- train_set %>% 
#   group_by(year = year(week_starting_date)) %>% 
#   group_split()
# 
# future::plan(multisession, workers = 8)
# 
# possible_parameters_and_sse_cv <- possible_parameters %>% 
#   future_pmap_dfr(
#     function(formula, psi, scale.est, c, k) {
#       
#       results <- list()
#       
#       for (fold in seq_along(train_set_yearly_lst)) {
#         
#         train_set_fold <- train_set_yearly_lst[setdiff(seq_along(train_set_yearly_lst), fold)] %>% bind_rows()
#         test_set_fold  <- train_set_yearly_lst[[fold]]
#         
#         model_rlm <- if (psi == 'psi.bisquare') {
#           
#           rlm(
#             formula   = as.formula(formula),
#             data      = train_set_fold,
#             method    = 'M',
#             scale.est = scale.est,
#             psi       = 'psi.bisquare',
#             c         = c,
#             maxit     = 200
#           )
#           
#         } else if (psi == 'psi.huber') {
#           
#           rlm(
#             formula   = as.formula(formula),
#             data      = train_set_fold,
#             method    = 'M',
#             scale.est = scale.est,
#             psi       = 'psi.huber',
#             k         = k,
#             maxit     = 200
#           )
#           
#         }
#         
#         rmse <- test_set_fold %>% 
#           mutate(
#             predicted = predict(model_rlm, test_set_fold),
#             se        = (observed - predicted)^2
#           ) %>% 
#           .$se %>% 
#           mean() %>% 
#           sqrt()
#         
#         results[[fold]] <- tibble(
#           rmse = rmse,
#           converged = if ('converged' %in% names(model_rlm)) model_rlm$converged else FALSE
#         )
#         
#       }
#       
#       tibble(
#         formula = formula,
#         psi = psi,
#         scale.est = scale.est, 
#         c = c, 
#         k = k
#       ) %>% 
#         bind_cols(
#           results %>% 
#             bind_rows() %>% 
#             summarise(
#               mean_rmse = mean(rmse),
#               converged = sum(converged)
#             )
#         )
# 
#     },
#     .id = 'iteration'
#   )
# 
# future::plan(sequential)
# 
# best_candidate_params <- possible_parameters_and_sse_cv %>% 
#   filter(mean_rmse == min(mean_rmse)) %>% 
#   filter(converged == length(train_set_yearly_lst)) %>% 
#   filter(row_number() == 1)
# 
# # TODO: Not much has changed; but do all the stuff - graph, calculate excess mortality, etc.
# best_candidate_results <- predict_year_rlm(
#   year_to_predict    = 2021, 
#   full_time_series   = abs_deaths, 
#   date_col           = 'week_starting_date', 
#   values_col         = 'observed', 
#   num_years_baseline = 5,
#   formula            = best_candidate_params$formula,
#   verbose            = FALSE,
#   method             = 'M',
#   scale.est          = best_candidate_params$scale.est,
#   psi                = best_candidate_params$psi,
#   c                  = best_candidate_params$c,
#   maxit              = 200
# )
# 
# best_candidate_results$results %>% 
#   ggplot(aes(x = week_starting_date)) +
#   geom_line(aes(y = observed, colour = 'Observed')) +
#   geom_line(aes(y = expected, colour = 'ABS Predicted')) +
#   geom_line(aes(y = predicted, colour = 'Herd Mentality Predicted')) +
#   scale_colour_manual(values = c(
#     'Observed' = 'grey', 
#     'ABS Predicted' = '#EC4899', 
#     'Herd Mentality Predicted' = '#14B8A6'
#   )) +
#   scale_y_continuous(labels = scales::comma) +
#   labs(
#     x = 'Week', y = 'Deaths',
#     colour  = element_blank(),
#     title   = 'Closest parameters to ABS model',
#     caption = paste0(
#       'Closest candidate model -\n',
#       str_interp('  Scale estimator: ${abs_closest_candidate$scale.est}\n'),
#       str_interp('  Weighting function: ${abs_closest_candidate$psi}\n'),
#       str_interp('  Outlier threshold: ${ifelse(is.na(abs_closest_candidate$k), abs_closest_candidate$c, abs_closest_candidate$k)}')
#     )
#   ) +
#   theme(
#     plot.caption = element_text(hjust = 0), 
#     legend.position = 'bottom'
#   )
# 
# abs_closest_candidate_results$results %>% 
#   filter(year(week_starting_date) == 2021) %>% 
#   mutate(
#     abs_em = observed - expected,
#     hm_em  = observed - predicted
#   ) %>% 
#   summarise(across(ends_with('_em'), sum))

## Archive ----
if (FALSE) {
  
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
  
  c("MAD", "Huber") %>% 
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
  
}
