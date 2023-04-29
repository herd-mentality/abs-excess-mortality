library(tidyverse)
library(MASS)
library(janitor)

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
      `Week starting date` = dmy(`Week starting date`)
    ) %>% 
    filter(!is.na(`Week starting date`)) %>% 
    clean_names()
)

abs_deaths %>%
  dplyr::select(week_starting_date, expected, observed) %>% 
  pivot_longer(expected:observed, names_to = 'series', values_to = 'value') %>% 
  ggplot(aes(x = week_starting_date, y = value, colour = series)) +
  geom_line() +
  scale_y_continuous(labels = scales::comma)

