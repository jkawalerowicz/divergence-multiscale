# =============================================================================
# Multiscale Divergence Index (DI), Stockholm 2022
#
# Method:
#   - for every 100 m cell in Stockholm municipality (180), sum the group counts
#     of all cells whose centre lies within the radius (the cell itself included;
#     cells in the 10 km buffer count as neighbours)
#   - local DI = sum over groups of p * log(p / r), where p is the local share
#     and r the share in Stockholm municipality
#   - radii: 150, 300, 600, 1000 m
#
# Groups: sb = born in Sweden, fb_e = born in Europe, fb_ne = born outside Europe
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

dir.create("results", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Data: one row per 100 m cell (e, n = lower-left corner, SWEREF99 TM)
# -----------------------------------------------------------------------------
cnt <- read_csv("data/simulated_grid_2022.csv", show_col_types = FALSE,
                col_select = c(year, municipality, e, n, sb, fb_e, fb_ne)) %>%
  filter(year == 2022) %>%
  transmute(x = e, y = n, municipality, sb, fb_e, fb_ne,
            in_target = municipality == 180)

# Reference: composition of Stockholm municipality
ref <- cnt %>%
  filter(in_target) %>%
  summarise(pop = sum(sb) + sum(fb_e) + sum(fb_ne),
            r_sb = sum(sb) / pop, r_fb_e = sum(fb_e) / pop, r_fb_ne = sum(fb_ne) / pop)

# Target cells: Stockholm municipality only
tgt <- cnt %>% filter(in_target) %>% select(x, y)

# -----------------------------------------------------------------------------
# Local sums within each radius and the local DI
# -----------------------------------------------------------------------------
focal <- list()
for (b in c(150, 300, 600, 1000)) {
  # offsets (in metres) of the cells within radius b
  k   <- floor(b / 100) * 100
  off <- expand_grid(dx = seq(-k, k, 100), dy = seq(-k, k, 100)) %>% filter(dx^2 + dy^2 <= b^2)

  # every target cell paired with every neighbour position; neighbour counts summed per target cell
  focal[[as.character(b)]] <- cross_join(tgt, off) %>%
    mutate(nx = x + dx, ny = y + dy) %>%
    inner_join(cnt %>% select(nx = x, ny = y, sb, fb_e, fb_ne), by = c("nx", "ny")) %>%
    group_by(x, y) %>%
    summarise(across(c(sb, fb_e, fb_ne), sum), .groups = "drop") %>%
    cross_join(ref %>% select(-pop)) %>%
    mutate(pop = sb + fb_e + fb_ne,
           DI  = if_else(pop > 0,
                         if_else(sb    > 0, (sb    / pop) * log((sb    / pop) / r_sb),    0) +
                         if_else(fb_e  > 0, (fb_e  / pop) * log((fb_e  / pop) / r_fb_e),  0) +
                         if_else(fb_ne > 0, (fb_ne / pop) * log((fb_ne / pop) / r_fb_ne), 0),
                         NA_real_),
           bw = b) %>%
    select(x, y, bw, local_pop = pop, DI)

  cat(sprintf("local sums at %4dm: %d cells in the radius\n", b, nrow(off)))
}

# Local DI for every cell and radius, with the cell's own population
local <- bind_rows(focal) %>%
  left_join(cnt %>% transmute(x, y, own_pop = sb + fb_e + fb_ne), by = c("x", "y"))

# -----------------------------------------------------------------------------
# Summary: DI at each radius
# -----------------------------------------------------------------------------
# DI = mean of local DI weighted by each cell's own population
summary_bw <- local %>%
  group_by(bw) %>%
  summarise(DI_local_med = median(DI, na.rm = TRUE),
            DI_local_p90 = quantile(DI, 0.9, na.rm = TRUE),
            DI           = weighted.mean(DI, own_pop, na.rm = TRUE),   # last: it overwrites DI
            cells        = n(),
            pop          = sum(own_pop))
print(summary_bw)

write_csv(summary_bw, "results/DI_by_radius.csv")
write_csv(local,      "results/DI_local_cells.csv")   # cell level: keep out of GitHub

# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------
radius_labels <- c("150m", "300m", "600m", "1000m")
theme_set(theme_bw(10))

ggsave("results/fig_DI_by_radius.png", width = 6, height = 4, dpi = 200,
       ggplot(summary_bw, aes(bw, DI)) +
         geom_line() +
         geom_point(size = 2) +
         scale_x_continuous(breaks = c(150, 300, 600, 1000)) +
         expand_limits(y = 0) +
         labs(x = "Radius (m)", y = "Divergence Index", title = "Segregation by radius"))

ggsave("results/fig_map_local_DI.png", width = 9, height = 8, dpi = 200,
       ggplot(local %>% mutate(radius = factor(paste0(bw, "m"), radius_labels)),
              aes(x + 50, y + 50, fill = DI)) +
         geom_tile(width = 100, height = 100) +
         facet_wrap(~ radius) +
         scale_fill_viridis_c(option = "magma", direction = -1,
                              limits = c(0, quantile(local$DI, 0.99, na.rm = TRUE)), oob = scales::squish) +
         coord_equal() +
         labs(x = NULL, y = NULL, fill = "Local DI", title = "Local Divergence Index") +
         theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()))
