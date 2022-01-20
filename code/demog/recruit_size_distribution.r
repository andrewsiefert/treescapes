library(tidyverse)
library(brms)

options(mc.cores = parallel::detectCores())

# read in recruit size data
d <- read_csv("data/recruit_size_data.csv") %>%
  mutate(z1 = (dia-2.54)/remper)

# fit model
m <- brm(bf(z1 ~ 1,
            hu ~ 1),
    data = d, family = hurdle_gamma())

samples <- posterior_samples(m, pars = c("b", "shape"))
write_csv(samples, "results/demog_pars/recr_size_pars.csv")
