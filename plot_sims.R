suppressPackageStartupMessages({
  require(ggplot2)
  require(cowplot)
})

.args <- commandArgs(trailingOnly = TRUE)

dt <- readRDS(.args[1])

p <- ggplot(dt) + aes(time, cumcases, group=n) +
  geom_step(alpha = 0.1) + theme_minimal() +
  scale_x_continuous("day") +
  scale_y_continuous("cumulative cases")

save_plot(tail(.args, 1), p)