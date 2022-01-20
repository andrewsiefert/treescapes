library(tidyverse)
library(rethinking)
library(rstan)
library(bayesplot)

source("code/transformer.r")
tf <- readRDS("data/transformers.rds")


mat_levs <- 10 %>% transform(tf$env$mat)
size_levs <- 20 %>% transform(tf$tree$prevdia)
logsize_levs <- backtransform(size_levs, tf$tree$prevdia) %>% transform(tf$tree$log_prevdia, log = T)
crowd_seq <- seq(-3.5, 2.5, length.out = 100)

# Plot setup----

## Labels ----

ba_lab <- bquote(Neighbor~basal~area~(m^2))
growth_lab <- expression(atop(Diameter~growth~rate, (mm~yr^-1)))
s_surv_lab <- expression(atop(Sapling~survival, probability~(5~yr)))
c_surv_lab <- expression(atop(Canopy~tree~survival, probability~(5~yr)))
recruit_lab <- expression(atop(Recruitment, (recruits~ind^-1~yr^-1)))

title_size <- 10
label_size <- 8

## Theme----

theme_top <- theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = title_size),
        axis.text.y = element_text(size = label_size),
        strip.text.x = element_text(size = title_size, margin = margin(t=2, b=2, unit="pt")), 
        legend.justification = "left",
        legend.title = element_text(size = title_size),
        legend.text = element_text(size = label_size), 
        legend.key.height = unit(0.4, "cm"))

theme_middle <- theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = title_size),
        axis.text.y = element_text(size = label_size),
        strip.background = element_blank(),
        strip.text.x = element_blank(),         
        legend.justification = "left",
        legend.title = element_text(size = title_size),
        legend.text = element_text(size = label_size), 
        legend.key.height = unit(0.4, "cm"))

theme_bottom <- theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = title_size),
        axis.text.x = element_text(size = label_size),
        axis.title.y = element_text(size = title_size),
        axis.text.y = element_text(size = label_size),
        strip.background = element_blank(),
        strip.text.x = element_blank(),    
        legend.justification = "left",
        legend.title = element_text(size = title_size),
        legend.text = element_text(size = label_size), 
        legend.key.height = unit(0.4, "cm"))

# Growth ------------------------------------------------------------------

## Prepare data ----
data <- read_csv("data/growth_data_train.csv")

d <- data %>% 
  filter(ann_dia_growth < 40) %>%
  mutate(wd = transform(wood_density, tf$trait$wood_density_log, log = T),
         sla = transform(sla, tf$trait$sla), 
         hmax = transform(hmax, tf$trait$hmax_log, log = T),
         mat = transform(mat, tf$env$mat),
         size = transform(prevdia, tf$tree$prevdia),
         log_size = transform(prevdia, tf$tree$log_prevdia, log = T),
         crowd = transform(canopy_nbr_ba, tf$tree$log_canopy_nbr_ba, log = T), 
         wd_sq = wd^2,
         sla_sq = sla^2,
         hmax_sq = hmax^2,
         sp =  as.numeric(factor(sp)),
         site = as.numeric(factor(site)))

f <- lm(ann_dia_growth ~ 1 + mat + log_size + size + crowd + 
          wd + sla + hmax + 
          wd_sq + sla_sq + hmax_sq + 
          wd:sla + wd:hmax + sla:hmax + 
          wd:log_size + sla:log_size + hmax:log_size +
          wd:size + sla:size + hmax:size +
          mat:wd + mat:sla + mat:hmax + 
          mat:wd_sq + mat:sla_sq + mat:hmax_sq +
          mat:wd:sla + mat:wd:hmax + mat:sla:hmax,
        data = d)

X <- model.matrix(f)


## Load model ----
fit <- readRDS("results/models/growth_model.rds")

samples <- rstan::extract(fit)
colnames(samples$beta_pop) <- colnames(X)

post <- as_tibble(samples$beta_pop) %>% rename(int = 1) %>% janitor::clean_names()


## Average crowding response ----

### Calculate growth ----
df <- expand.grid(int = 1, 
                  mat = mat_levs, 
                  size = size_levs,
                  log_size = logsize_levs,
                  crowd = crowd_seq) %>%
  as_tibble() 

beta <- post[,names(df)] %>% as.matrix()

yo <- sapply(1:nrow(beta), function(i) exp(as.matrix(df) %*% beta[i,]))

df <- df %>% 
  mutate(Growth = rowMeans(yo),
         lo = apply(yo, 1, quantile, 0.05),
         hi = apply(yo, 1, quantile, 0.95),
         crowding = backtransform(crowd, tf$tree$log_canopy_nbr_ba, log = T))

### Plot average growth ----
g1 <- ggplot(df, aes(x = crowding, y = Growth, ymin = lo, ymax = hi)) +
  geom_path() +
  geom_ribbon(alpha = 0.2) +
  labs(x = ba_lab, y = growth_lab) +
  theme_top


## Species responses ----

### Calculate growth ----
df <- expand.grid(int = 1, 
                  #mat = mat_levs, 
                  #size = size_levs,
                  #log_size = logsize_levs,
                  crowd = crowd_seq) %>%
  as_tibble() 

beta <- post[,names(df)] %>% colMeans()
sp_int <- colMeans(samples$sp_eff[,,1])
sp_crowd <- colMeans(samples$sp_eff[,,5])

yo <- matrix(NA, nrow = 100, ncol = dim(samples$sp_eff)[2])

for(i in 1:100) {
  for(s in 1:97) {
    yo[i,s] <- exp(as.matrix(df[i,]) %*% beta + sp_int[s] + df$crowd[i] * sp_crowd[s])
  }
}

df2 <- yo %>% 
  as_tibble() %>% 
  mutate(crowding = df$crowd %>% backtransform(tf$tree$log_canopy_nbr_ba, log = T)) %>% 
  gather("sp", "growth", -crowding)


### Plot species growth ----
g2 <- ggplot(df2, aes(x = crowding, y = growth, group = sp)) +
  geom_path(alpha = 0.5) +
  labs(x = ba_lab, y = growth_lab) +
  theme_top



# Sapling survival ------------------------------------------------------------------

## Prepare data ----
data <- read_csv("data/survival_sapling_data_train.csv")

d <- data %>% 
  mutate(wd = transform(wood_density, tf$trait$wood_density_log, log = T),
         sla = transform(sla, tf$trait$sla), 
         hmax = transform(hmax, tf$trait$hmax_log, log = T),
         mat = transform(mat, tf$env$mat),
         logsize = transform(prevdia, tf$sapling$log_prevdia, log = T),
         crowd = transform(canopy_nbr_ba, tf$tree$log_canopy_nbr_ba, log = T), 
         wd_sq = wd^2,
         sla_sq = sla^2,
         hmax_sq = hmax^2)

f <- lm(surv ~ 1 + mat + logsize + crowd + 
          wd + sla + hmax + 
          wd_sq + sla_sq + hmax_sq + 
          wd:sla + wd:hmax + sla:hmax + 
          wd:logsize + sla:logsize + hmax:logsize +
          mat:wd + mat:sla + mat:hmax + 
          mat:wd_sq + mat:sla_sq + mat:hmax_sq +
          mat:wd:sla + mat:wd:hmax + mat:sla:hmax,
        data = d)

X <- model.matrix(f)


## Load model ----
fit <- readRDS("results/models/survival_sapling_model.rds")

samples <- rstan::extract(fit)
colnames(samples$beta_pop) <- colnames(X)

post <- as_tibble(samples$beta_pop) %>% rename(int = 1) %>% janitor::clean_names()


## Average crowding response ----

### Calculate growth ----
df <- expand.grid(int = 1, 
                  crowd = crowd_seq) %>%
  as_tibble() 

beta <- post[,names(df)] %>% as.matrix()

yo <- sapply(1:nrow(beta), function(i) plogis(as.matrix(df) %*% beta[i,])^5)

df <- df %>% 
  mutate(surv = rowMeans(yo),
         lo = apply(yo, 1, quantile, 0.05),
         hi = apply(yo, 1, quantile, 0.95),
         crowding = backtransform(crowd, tf$tree$log_canopy_nbr_ba, log = T))

### Plot growth ----
ss1 <- ggplot(df, aes(x = crowding, y = surv, ymin = lo, ymax = hi)) +
  geom_path() +
  geom_ribbon(alpha = 0.2) +
  labs(x = ba_lab, y = s_surv_lab) + 
  theme_middle


## Species responses ----

### Calculate growth ----
df <- expand.grid(int = 1, 
                  crowd = crowd_seq) %>%
  as_tibble() 

n <- nrow(df)
n_sp <- dim(samples$sp_eff)[2]

beta <- post[,names(df)] %>% colMeans()
sp_int <- colMeans(samples$sp_eff[,,1])
sp_crowd <- colMeans(samples$sp_eff[,,4])

yo <- matrix(NA, nrow = n, ncol = n_sp)

for(i in 1:n) {
  for(s in 1:n_sp) {
    yo[i,s] <- plogis(as.matrix(df[i,]) %*% beta + sp_int[s] + df$crowd[i] * sp_crowd[s])^5
  }
}

df2 <- yo %>% 
  as_tibble() %>% 
  mutate(crowding = df$crowd %>% backtransform(tf$tree$log_canopy_nbr_ba, log = T)) %>% 
  gather("sp", "surv", -crowding)


### Plot growth ----
ss2 <- ggplot(df2, aes(x = crowding, y = surv, group = sp)) +
  geom_path(alpha = 0.5) +
  labs(x = ba_lab, y = s_surv_lab) + 
  theme_middle



# Canopy survival ------------------------------------------------------------------

## Prepare data ----
data <- read_csv("data/survival_canopy_data_train.csv")

d <- data %>% 
  mutate(wd = transform(wood_density, tf$trait$wood_density_log, log = T),
         sla = transform(sla, tf$trait$sla), 
         hmax = transform(hmax, tf$trait$hmax_log, log = T),
         mat = transform(mat, tf$env$mat),
         size = transform(prevdia, tf$canopy$prevdia),
         log_size = transform(prevdia, tf$canopy$log_prevdia, log = T),
         crowd = transform(canopy_nbr_ba, tf$tree$log_canopy_nbr_ba, log = T), 
         wd_sq = wd^2,
         sla_sq = sla^2,
         hmax_sq = hmax^2)

f <- lm(surv ~ 1 + mat + size + log_size + crowd + 
          wd + sla + hmax + 
          wd_sq + sla_sq + hmax_sq + 
          wd:sla + wd:hmax + sla:hmax + 
          wd:size + sla:size + hmax:size +
          wd:log_size + sla:log_size + hmax:log_size +
          mat:wd + mat:sla + mat:hmax + 
          mat:wd_sq + mat:sla_sq + mat:hmax_sq +
          mat:wd:sla + mat:wd:hmax + mat:sla:hmax,
        data = d)

X <- model.matrix(f)


## Load model ----
fit <- readRDS("results/models/survival_canopy_model.rds")

samples <- rstan::extract(fit)
colnames(samples$beta_pop) <- colnames(X)

post <- as_tibble(samples$beta_pop) %>% rename(int = 1) %>% janitor::clean_names()


## Average crowding response ----

### Calculate growth ----
df <- expand.grid(int = 1, 
                  crowd = crowd_seq) %>%
  as_tibble() 

beta <- post[,names(df)] %>% as.matrix()

yo <- sapply(1:nrow(beta), function(i) plogis(as.matrix(df) %*% beta[i,])^5)

df <- df %>% 
  mutate(surv = rowMeans(yo),
         lo = apply(yo, 1, quantile, 0.05),
         hi = apply(yo, 1, quantile, 0.95),
         crowding = backtransform(crowd, tf$tree$log_canopy_nbr_ba, log = T))

### Plot growth ----
sc1 <- ggplot(df, aes(x = crowding, y = surv, ymin = lo, ymax = hi)) +
  geom_path() +
  geom_ribbon(alpha = 0.2) +
  labs(x = ba_lab, y = c_surv_lab) + 
  theme_middle


## Species responses ----

### Calculate growth ----
df <- expand.grid(int = 1, 
                  crowd = crowd_seq) %>%
  as_tibble() 

n <- nrow(df)
n_sp <- dim(samples$sp_eff)[2]

beta <- post[,names(df)] %>% colMeans()
sp_int <- colMeans(samples$sp_eff[,,1])
sp_crowd <- colMeans(samples$sp_eff[,,5])

yo <- matrix(NA, nrow = n, ncol = n_sp)

for(i in 1:n) {
  for(s in 1:n_sp) {
    yo[i,s] <- plogis(as.matrix(df[i,]) %*% beta + sp_int[s] + df$crowd[i] * sp_crowd[s])^5
  }
}

df2 <- yo %>% 
  as_tibble() %>% 
  mutate(crowding = df$crowd %>% backtransform(tf$tree$log_canopy_nbr_ba, log = T)) %>% 
  gather("sp", "surv", -crowding)


### Plot growth ----
sc2 <- ggplot(df2, aes(x = crowding, y = surv, group = sp)) +
  geom_path(alpha = 0.5) +
  labs(x = ba_lab, y = c_surv_lab) + 
  theme_middle



# Recruitment ------------------------------------------------------------------

## Prepare data ----
data <- readRDS("data/recruitment_data_train.rds")

d <- data$df %>% 
  mutate(wd = transform(wood_density, tf$trait$wood_density_log, log = T),
         sla = transform(sla, tf$trait$sla), 
         hmax = transform(hmax, tf$trait$hmax_log, log = T),
         mat = transform(mat, tf$env$mat),
         crowd_c = transform(canopy_nbr_ba, tf$recr$canopy_nbr_ba_log, log = T), 
         crowd_s = transform(sapling_nbr_ba, tf$recr$sapling_nbr_ba),
         wd_sq = wd^2,
         sla_sq = sla^2,
         hmax_sq = hmax^2)

f <- lm(ingrowth ~ mat + crowd_c + crowd_s + 
          wd + sla + hmax + 
          wd_sq + sla_sq + hmax_sq + 
          wd:sla + wd:hmax + sla:hmax + 
          mat:wd + mat:sla + mat:hmax +
          mat:wd_sq + mat:sla_sq + mat:hmax_sq +
          mat:wd:sla + mat:wd:hmax + mat:sla:hmax,
        data = d)

X <- model.matrix(f)[,-1]

## Load model ----
fit <- readRDS("results/models/recruitment_model.rds")

samples <- rstan::extract(fit)
colnames(samples$beta_pop) <- colnames(X)


## Average crowding response ----

### Calculate recruitment ----
df <- expand.grid(crowd_c = modelr::seq_range(d$crowd_c, 100)) %>%
  as_tibble() 

n <- nrow(df)
n_sample <- nrow(samples$beta_pop)

yo <- matrix(NA, nrow = n, ncol = n_sample)
for(s in 1:n_sample) {
  for(i in 1:nrow(df)) {
    str <- exp(samples$str_pop[s])
    filter <- exp(samples$beta_pop[s,'crowd_c']*df$crowd_c[i])
    yo[i,s] <- str * filter
  }
}

df <- df %>% 
  mutate(recr = rowMeans(yo),
         lo = apply(yo, 1, quantile, 0.05),
         hi = apply(yo, 1, quantile, 0.95),
         crowding = backtransform(crowd_c, tf$recr$canopy_nbr_ba_log, log = T))

### Plot growth ----
r1 <- ggplot(df, aes(x = crowding, y = recr, ymin = lo, ymax = hi)) +
  geom_path() +
  geom_ribbon(alpha = 0.2) +
  labs(x = ba_lab, y = recruit_lab) +
  theme_bottom


## Species responses ----

### Calculate growth ----
df <- expand.grid(crowd_c = modelr::seq_range(d$crowd_c, 100)) %>%
  as_tibble() 

n <- nrow(df)
n_sp <- dim(samples$str)[2]


str_sp <- colMeans(samples$str)
crowd_sp <- colMeans(samples$sp_eff[,,6]) + mean(samples$beta_pop[s,'crowd_c'])

yo <- matrix(NA, nrow = n, ncol = n_sp)
for(s in 1:n_sp) {
  for(i in 1:nrow(df)) {
    filter <- exp(crowd_sp[s] * df$crowd_c[i])
    yo[i,s] <- str_sp[s] * filter
  }
}

df2 <- yo %>% 
  as_tibble() %>% 
  mutate(crowding = df$crowd_c %>% backtransform(tf$recr$canopy_nbr_ba_log, log = T)) %>% 
  gather("sp", "recr", -crowding)


### Plot growth ----
r2 <- ggplot(df2, aes(x = crowding, y = recr, group = sp)) +
  geom_path(alpha = 0.5) +
  labs(x = ba_lab, y = recruit_lab) +
  theme_bottom

# Plot ------------------------------------------------------------------------

pdf("results/figures/ed_fig9.pdf", 6.5, 8)
gridExtra::grid.arrange(rbind(cbind(ggplotGrob(ss1), ggplotGrob(ss2)),
                              cbind(ggplotGrob(sc1), ggplotGrob(sc2)),
                              cbind(ggplotGrob(g1), ggplotGrob(g2)),
                              cbind(ggplotGrob(r1), ggplotGrob(r2))))
dev.off()



