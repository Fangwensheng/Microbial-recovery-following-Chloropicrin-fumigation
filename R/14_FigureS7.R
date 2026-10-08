# ============================================================
# Figure S7

library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(scales)

dat <- read_csv("FigureS7_plot_data.csv", show_col_types = FALSE)
st  <- read_csv("FigureS7_statistics.csv", show_col_types = FALSE)

dat <- dat %>%
  mutate(
    Context = factor(Context, levels = c("Unplanted", "Tomato-planted")),
    Panel = factor(Panel, levels = c("FigS7a","FigS7b","FigS7c","FigS7d"))
  )


context_cols <- c(
  "Unplanted" = "#2A7F9E",
  "Tomato-planted" = "#D97943"
)

theme_fig <- theme_classic(base_size = 8.5) +
  theme(
    text = element_text(family = "Arial", colour = "black"),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8, colour = "black"),
    axis.line = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    legend.position = "none",
    plot.margin = margin(4, 6, 4, 4)
  )

fmt_q <- function(x) {
  ifelse(
    x < 0.001,
    format(x, scientific = TRUE, digits = 2),
    format(round(x, 3), nsmall = 3)
  )
}

get_stat <- function(panel) {
  st %>% filter(Panel == panel) %>% slice(1)
}

make_panel <- function(panel_id, xlab, ylab) {
  dd <- dat %>% filter(Panel == panel_id)
  ss <- get_stat(panel_id)
  context_name <- as.character(unique(dd$Context))

  stat_lab <- paste0(
    "partial \u03C1 = ", sprintf("%.3f", as.numeric(ss$Partial_rho)),
    "\nq = ", fmt_q(as.numeric(ss$Partial_FDR))
  )

  ggplot(
    dd,
    aes(
      x = Adjusted_X_z,
      y = Adjusted_Y_z,
      colour = Context
    )
  ) +
    geom_point(
      size = 1.15,
      alpha = 0.28,
      stroke = 0
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.7
    ) +
    scale_colour_manual(values = context_cols) +
    labs(
      x = xlab,
      y = ylab,
      subtitle = context_name
    ) +
    annotate(
      "label",
      x = Inf, y = Inf,
      label = stat_lab,
      hjust = 1.03,
      vjust = 1.03,
      size = 2.5,
      label.size = 0,
      fill = alpha("white", 0.85),
      label.padding = unit(0.10, "lines")
    ) +
    theme_fig +
    theme(
      plot.subtitle = element_text(
        size = 8.5,
        face = "plain",
        hjust = 0
      )
    )
}

p_a <- make_panel(
  "FigS7a",
  "Adjusted early response diversity",
  "Adjusted late compensation efficiency"
)

p_b <- make_panel(
  "FigS7b",
  "Adjusted early response diversity",
  "Adjusted late compensation efficiency"
)

p_c <- make_panel(
  "FigS7c",
  "Adjusted integrated compensation efficiency",
  "Adjusted cumulative residual"
)

p_d <- make_panel(
  "FigS7d",
  "Adjusted integrated compensation efficiency",
  "Adjusted cumulative residual"
)

figS8 <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(
      family = "Arial",
      size = 11,
      face = "bold"
    )
  )

ggsave(
  "FigureS7_R.pdf",
  figS8,
  width = 174,
  height = 145,
  units = "mm",
  device = cairo_pdf
)
