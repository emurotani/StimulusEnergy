library(tidyverse)

set.seed(1234)

# in HFA, background luminance = 10 (cd)
dbTocntr <- function(db, bg = 10) {
  OPI::dbTocd(db)/bg
}

# michaelis-menten equation
mm <- function(cntr, saturation_spike = 100, cntr_50 = 1.00) {
  saturation_spike * cntr / (cntr_50 + cntr)
}

# simulate M-RGC spikes
# fig.1
tibble(dB = 0:40) %>% 
  mutate(contrast = dbTocntr(dB),
         f_rate = mm(contrast)) %>% 
  # filter(contrast < 2.0) %>%
  ggplot(aes(x = dB, y = f_rate)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  labs(x = "Threshold intensity (dB)",
       y = "Firing rates (spike/second)")

# calculate local M-RGC density
# threshold intensity at normal region was assumed to be 32dB.
tibble(dB = 0:40) %>% 
  mutate(contrast = dbTocntr(dB),
         f_rate = mm(contrast)) %>% 
  mutate(ref_spike = mm(dbTocntr(32))) %>% 
  mutate(relative_density = ref_spike / f_rate) -> df_local_dens

# add "<0dB" row
df_local_dens %>% 
  add_row(tibble(
    dB = -1L,
    contrast = NA,
    f_rate = NA,
    ref_spike = NA,
    relative_density = 0
  ),
  .before = 1) -> df_local_dens

df_local_dens

# for table in ARVO poster
df_local_dens %>% 
  filter(dB %in% c(32, 27, 22, 0)) %>% 
  arrange(desc(dB))

# fig.2
df_local_dens %>% 
  filter(between(dB, 0, 32)) %>% 
  ggplot(aes(x = dB, y = relative_density)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(breaks = seq(0, 1, 0.2), 
                     limits = c(0, 1.0),
                     labels = scales::percent) +
  labs(x = "Threshold intensity (dB)",
       y = "Relative M-RGC density")

# load parameters of empirical FOS curves
df_par_raw <- read_csv("estimated_parameters.csv")

df_par_raw
# sensitivity (full threshold) + 1 = sens (SITA equivalent)

df_par_raw %>% 
  transmute(sens = sensitivity + 1,
            par = parameter,
            value = mean) %>% 
  pivot_wider(names_from = "par",
              values_from = "value") -> df_par

df_par

# calculate weighted average of local M-RGC density
mean_dens <- function(mu, sigma, lamda) {
  tibble(threshold = 0:40, 
         weight = dnorm(0:40, mu, sigma)) %>% 
    mutate(weight = (1 - lamda) * weight) %>% 
    add_row(tibble(threshold = -1, weight = lamda), .before = 1) %>% 
    left_join(df_local_dens, by = c("threshold" = "dB")) %>% 
    mutate(dens = relative_density * weight) %>% 
    summarize(density = sum(dens)) %>% 
    as.numeric(density)
}


# example in ARVO poster
df_par %>% 
  filter(sens == 15)

mean_dens(19.7, 2.79, 0.503)

tibble(threshold = 0:40, 
       weight = dnorm(0:40, 19.7, 2.79)) %>% 
  mutate(weight = (1 - 0.503) * weight) %>% 
  add_row(tibble(threshold = -1, weight = 0.503), .before = 1) %>% 
  left_join(df_local_dens, by = c("threshold" = "dB")) %>% 
  mutate(dens = relative_density * weight) %>% 
  filter(threshold %in% c(-1, 19, 20, 21, 32)) %>% 
  select(1, 2, 3, 6, 7) %>% 
  arrange(desc(threshold))


# calculate mean M-RGC density in regions with glaucomatous damage
df_par %>% 
  mutate(density = pmap_dbl(.l = list(mu, sigma, lamda), 
                            .f = mean_dens)) -> df_mean_dens

# fig.3
df_mean_dens %>% 
  filter(sens <= 32) %>% 
  ggplot(aes(x = sens, y = density)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(breaks = seq(0, 1, 0.2), 
                     limits = c(0, 1.0),
                     labels = scales::percent) +
  labs(x = "Retinal sensitivity (dB)",
       y = "Predicted mean M-RGC density")


# for comparison
# fig 3A in Harwerth et al. IOVS 1999;40:2242-2250 
df_mean_dens %>% 
  mutate(loss_sens = 32 - sens,
         loss_density = 1 - density) %>% 
  filter(sens <= 32) %>% 
  ggplot(aes(x = loss_density, y = loss_sens)) +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 1, 0.2),
                     labels = scales::percent) +
  labs(x = "M-RGC loss (%)",
       y = "Retinal sensitivity loss (dB)")

