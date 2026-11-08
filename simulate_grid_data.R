# =============================================================================
# Simulated 100 m grid data for a fictional city
#
# Produces a file with the same columns as the real register grid
# (year, municipality, sb, fb_e, fb_ne, tot, n, e), so that
# divergence_multiscale.R runs on it unchanged.
#
# No real cell, location or count is used. The only inputs taken from the
# real data are a handful of city-wide aggregates (see "Targets" below):
# number of populated cells, quantiles of cell population, group shares and
# the DI by radius. The geography is invented: a round city with an
# irregular border, a dense centre, satellite centres, and a few
# neighbourhoods where people born outside Europe are concentrated.
#
# Groups: sb = born in Sweden, fb_e = born in Europe, fb_ne = born outside Europe
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)

set.seed(2022)
dir.create("data", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Targets (city-wide aggregates of the real data, Stockholm 2022)
# -----------------------------------------------------------------------------
# populated cells:         city ~10,700, buffer ~27,600
# persons per cell, city:  median 49, mean 92
# persons per cell, buffer: median 18, mean 35
# group shares, city:      sb 0.738, fb_e 0.103, fb_ne 0.159
# DI by radius (city):     150 m 0.079, 300 m 0.070, 600 m 0.060, 1000 m 0.053

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
r_city   <- 9000      # city radius (m); the border wobbles around it
r_region <- 19000     # city + 10 km buffer

# Smooth random surface: a sum of Gaussian bumps with random centres,
# widths and weights. Used for density, local noise and segregation.
bumps <- function(x, y, k, sd_min, sd_max, weight, r_min = 0, r_max = r_region) {
  rad <- sqrt(runif(k, (r_min / r_region)^2, (r_max / r_region)^2)) * r_region
  ang <- runif(k, 0, 2 * pi)
  x0  <- rad * cos(ang); y0 <- rad * sin(ang)
  s   <- runif(k, sd_min, sd_max)
  w   <- weight(k)
  z   <- numeric(length(x))
  for (i in seq_len(k)) z <- z + w[i] * exp(-((x - x0[i])^2 + (y - y0[i])^2) / (2 * s[i]^2))
  z
}
rescale01 <- function(z) (z - min(z)) / (max(z) - min(z))

# -----------------------------------------------------------------------------
# 1. Grid: 100 m cells inside the region; city = inside a wobbly circle
# -----------------------------------------------------------------------------
grid <- expand_grid(e = seq(-r_region, r_region - 100, 100),
                    n = seq(-r_region, r_region - 100, 100)) %>%
  mutate(cx = e + 50, cy = n + 50,                       # cell centre
         dist = sqrt(cx^2 + cy^2),
         ang  = atan2(cy, cx)) %>%
  filter(dist <= r_region) %>%
  mutate(border    = r_city * (1 + 0.12 * sin(3 * ang + 1) + 0.08 * sin(5 * ang)),
         in_target = dist <= border,
         # municipality: 180 = the city; buffer split into 8 fictional municipalities
         municipality = if_else(in_target, 180L, 901L + as.integer(floor((ang + pi) / (2 * pi) * 8)) %% 8L))

# -----------------------------------------------------------------------------
# 2. Where people live: density surface and populated cells
# -----------------------------------------------------------------------------
grid <- grid %>%
  mutate(
    # dense centre + satellite centres
    dens = rescale01(exp(-dist^2 / (2 * 7000^2)) +
                     bumps(cx, cy, 40, 600, 2000, function(k) runif(k, 0.2, 0.8), 2000, r_region)),
    # patchiness: parks, water, industrial areas (small-scale smooth noise)
    patch = bumps(cx, cy, 2500, 150, 400, function(k) rnorm(k)),
    patch = (patch - mean(patch)) / sd(patch),
    # a cell is populated with a probability rising with density
    p_pop = plogis(-1.6 + 3 * dens + 1.0 * patch),
    populated = runif(n()) < p_pop
  ) %>%
  filter(populated)

# -----------------------------------------------------------------------------
# 3. Who lives there: a segregation surface, persons per cell, group counts
# -----------------------------------------------------------------------------
# seg  = neighbourhoods (0.6-1.3 km across) where people born outside Europe
#        are concentrated: 6 in the outer city, 24 in the buffer. They are large
#        estates, so they also have more persons per cell.
# wide = milder broad variation across the region (2-7 km scale).
# The coefficients here and in group_probs() were tuned by hand until the DI
# by radius and the spread of local DI were close to the targets.
grid <- grid %>%
  mutate(seg  = pmin(1, bumps(cx, cy, 6, 600, 1300, function(k) runif(k, 0.6, 1), 3000, 8500) +    # in the city
                         bumps(cx, cy, 24, 600, 1300, function(k) runif(k, 0.6, 1), 9500, 17000)),   # in the buffer
         wide = bumps(cx, cy, 40, 1500, 3500, function(k) rnorm(k)),
         wide = (wide - mean(wide)) / sd(wide),
         # persons per cell: log-normal, larger where density is higher and in the estates
         pop  = pmax(1L, as.integer(round(exp(2.0 + 3.8 * dens + 1.2 * seg + 0.2 * patch + rnorm(n(), 0, 1.05))))),
         u_e  = rnorm(n(), 0, 0.5),                    # cell-level variation
         u_ne = rnorm(n(), 0, 0.9))

# Group probabilities (multinomial logit, sb = reference). The intercepts are
# tuned so the city's shares match the targets.
group_probs <- function(d, a_e, a_ne) {
  d %>% mutate(eta_e  = a_e  + 1.4 * seg + 0.4 * wide + u_e,
               eta_ne = a_ne + 5.0 * seg + 1.0 * wide + u_ne,
               den    = 1 + exp(eta_e) + exp(eta_ne),
               p_sb = 1 / den, p_e = exp(eta_e) / den, p_ne = exp(eta_ne) / den)
}
a_e <- -2; a_ne <- -2
for (i in 1:30) {                                       # simple fixed-point tuning
  s <- group_probs(grid, a_e, a_ne) %>% filter(in_target) %>%
    summarise(e = weighted.mean(p_e, pop), ne = weighted.mean(p_ne, pop))
  a_e  <- a_e  + log(0.103 / s$e)
  a_ne <- a_ne + log(0.159 / s$ne)
}

# Draw the counts: each person is assigned a group with the cell's probabilities
sim <- group_probs(grid, a_e, a_ne) %>%
  mutate(sb    = rbinom(n(), pop, p_sb),
         fb_e  = rbinom(n(), pop - sb, p_e / (p_e + p_ne)),
         fb_ne = pop - sb - fb_e) %>%
  transmute(year = 2022L, municipality, sb, fb_e, fb_ne, tot = pop, n, e)

write_csv(sim, "data/simulated_grid_2022.csv")

# -----------------------------------------------------------------------------
# 4. Check against the targets
# -----------------------------------------------------------------------------
sim %>%
  mutate(area = if_else(municipality == 180, "city", "buffer")) %>%
  group_by(area) %>%
  summarise(cells = n(), median_pop = median(tot), mean_pop = round(mean(tot), 1),
            sb = round(sum(sb) / sum(tot), 3), fb_e = round(sum(fb_e) / sum(tot), 3),
            fb_ne = round(sum(fb_ne) / sum(tot), 3)) %>%
  print()
