library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(grid)

dat <- read_csv("Figure7_V3_plot_data.csv", show_col_types = FALSE)
st  <- read_csv("Figure7_V3_statistics.csv", show_col_types = FALSE)

dat <- dat %>%
  mutate(
    Context = factor(Context, levels = c("Unplanted", "Tomato-planted")),
    Recovery_class = factor(
      Recovery_class,
      levels = c("Rapid","Overshoot/reversal","Gradual/incomplete","Persistent/strengthened")
    )
  )

# Replace these with the exact colors used in Fig.1–6 if necessary.
context_cols <- c("Unplanted"="#2A7F9E", "Tomato-planted"="#D97943")
context_shapes <- c("Unplanted"=16, "Tomato-planted"=17)
context_ltypes <- c("Unplanted"="solid", "Tomato-planted"="dashed")

theme_fig <- theme_classic(base_size = 8.5) +
  theme(
    text = element_text(family="Arial", colour="black"),
    axis.title = element_text(size=9),
    axis.text = element_text(size=8, colour="black"),
    axis.line = element_line(linewidth=0.45),
    axis.ticks = element_line(linewidth=0.4),
    legend.title = element_blank(),
    legend.text = element_text(size=7.5),
    plot.margin = margin(4,6,4,4)
  )

fmt_q <- function(x){
  ifelse(x < 0.001,
         format(x, scientific=TRUE, digits=2),
         format(round(x,3), nsmall=3))
}

get_stat <- function(panel, context){
  st %>% filter(Panel==panel, Context==context) %>% slice(1)
}

# ---------- Fig. 7a ----------
sa_u <- get_stat("Fig7a","Unplanted")
sa_t <- get_stat("Fig7a","Tomato-planted")
lab_a <- paste0(
  "Unplanted: \u03C1=",sprintf("%.3f",sa_u$Spearman_rho),", q=",fmt_q(sa_u$FDR),
  "\nTomato-planted: \u03C1=",sprintf("%.3f",sa_t$Spearman_rho),", q=",fmt_q(sa_t$FDR)
)
yr_a <- range(dat$Early_RD_7_14, na.rm=TRUE)
ypad_a <- diff(yr_a)*0.12

p_a <- ggplot(dat, aes(Provider_richness, Early_RD_7_14,
                       colour=Context, shape=Context, linetype=Context)) +
  geom_point(size=1.15, alpha=0.24, stroke=0) +
  geom_smooth(method="lm", formula=y~x, se=FALSE, linewidth=0.7) +
  scale_colour_manual(values=context_cols) +
  scale_shape_manual(values=context_shapes) +
  scale_linetype_manual(values=context_ltypes) +
  coord_cartesian(
    xlim=c(0.5, NA),
    ylim=c(yr_a[1]-ypad_a, yr_a[2]+ypad_a),
    clip="off"
  ) +
  labs(x="Provider richness", y="Early response diversity") +
  annotate(
    "label",
    x=Inf, y=Inf, label=lab_a,
    hjust=1.03, vjust=1.03, size=2.45,
    label.size=0, fill="white"
  ) +
  theme_fig +
  theme(legend.position="bottom")

# ---------- Fig. 7b ----------
sb_u <- get_stat("Fig7b","Unplanted")
sb_t <- get_stat("Fig7b","Tomato-planted")
lab_b <- paste0(
  "Unplanted: \u03C1=",sprintf("%.3f",sb_u$Spearman_rho),", q=",fmt_q(sb_u$FDR),
  "\nTomato-planted: \u03C1=",sprintf("%.3f",sb_t$Spearman_rho),", q=",fmt_q(sb_t$FDR)
)
yr_b <- range(dat$Late_CE_28_70, na.rm=TRUE)
ylim_b <- c(max(0,yr_b[1]-0.08), min(1.02,yr_b[2]+0.08))

p_b <- ggplot(dat, aes(Early_RD_7_14, Late_CE_28_70,
                       colour=Context, shape=Context, linetype=Context)) +
  geom_point(size=1.15, alpha=0.24, stroke=0) +
  geom_smooth(method="lm", formula=y~x, se=FALSE, linewidth=0.7) +
  scale_colour_manual(values=context_cols) +
  scale_shape_manual(values=context_shapes) +
  scale_linetype_manual(values=context_ltypes) +
  coord_cartesian(ylim=ylim_b, clip="off") +
  labs(x="Early response diversity", y="Late compensation efficiency") +
  annotate(
    "label",
    x=Inf, y=Inf, label=lab_b,
    hjust=1.03, vjust=1.03, size=2.45,
    label.size=0, fill="white"
  ) +
  theme_fig +
  theme(legend.position="bottom")

# ---------- Fig. 7c ----------
# ========== 1. 参数预设（原配色、线型完全保留） ==========
context_cols <- c("#1F78B4", "#E67E22")  # 蓝、橙，和原图色调一致
context_ltypes <- c("solid", "dashed")     # 实线、虚线分组

# ========== 2. 纵轴范围（可按需微调边距） ==========
yr_c <- range(dat$Late_CE_28_70, na.rm = TRUE)
ylim_c <- c(max(-0.2, yr_c[1] - 0.5), min(2, yr_c[2] + 0.2))

# ========== 3. 绘图主体 ==========
p_c <- ggplot(dat, aes(x = Recovery_class, 
                       y = Late_CE_28_70,
                       colour = Context, 
                       linetype = Context)) +
  
  # --- 新增：原始散点图层（放在箱线下方，对齐分组） ---
  geom_jitter(
    position = position_dodge(width = 0.6),  # 和箱体位置严格对齐
    size = 1.2,     # 散点大小
    alpha = 0.7,    # 半透明避免重叠遮挡
    shape = 16      # 实心圆点样式
  ) +
  
  # --- 箱线图：缩窄箱体 + 确保误差棒两端显示 ---
  geom_boxplot(
    fill = NA,                     # 空心箱体，和原图风格一致
    position = position_dodge(width = 0.6),
    width = 0.45,                  # 箱体宽度缩窄，解决“箱体太长/太宽”
    outlier.shape = NA,            # 不重复显示异常值（已叠加散点）
    linewidth = 0.6,
    fatten = 2.5,                  # 中位数线加粗，更清晰
    outlier.colour = NA
  ) +
  stat_summary(
    fun = mean,
    fun.min = function(x) mean(x) - sd(x)/sqrt(length(x)),
    fun.max = function(x) mean(x) + sd(x)/sqrt(length(x)),
    geom = "errorbar",
    width = 0.15,        # 误差棒端帽宽度
    linewidth = 0.5,
    position = position_dodge(width = 0.6)
  ) +
  
  # 颜色、线型手动映射（完全沿用原逻辑）
  scale_colour_manual(values = context_cols) +
  scale_linetype_manual(values = context_ltypes) +
  
  # 纵轴范围
  coord_cartesian(ylim = ylim_c, clip = "off") +
  
  # x轴标签换行
  scale_x_discrete(
    labels = c(
      "Rapid" = "Rapid",
      "Overshoot/reversal" = "Overshoot/\nreversal",
      "Gradual/incomplete" = "Gradual/\nincomplete",
      "Persistent/strengthened" = "Persistent/\nstrengthened"
    )
  ) +
  
  labs(x = NULL, y = "Late compensation efficiency") +
  
  # 科研风格主题
  theme_classic() +
  theme(
    axis.text = element_text(size = 9, colour = "black"),
    axis.text.x = element_text(size = 9),
    axis.title.y = element_text(size = 10, margin = margin(r = 8)),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.width = unit(1.2, "cm"),
    axis.line = element_line(colour = "black", linewidth = 0.4),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )

p_c

# ---------- Fig. 7d ----------
sd_u <- get_stat("Fig7d","Unplanted")
sd_t <- get_stat("Fig7d","Tomato-planted")
lab_d <- paste0(
  "Unplanted: \u03C1=",sprintf("%.3f",sd_u$Spearman_rho),", q=",fmt_q(sd_u$FDR),
  "\nTomato-planted: \u03C1=",sprintf("%.3f",sd_t$Spearman_rho),", q=",fmt_q(sd_t$FDR)
)
xr_d <- range(dat$Integrated_CE_7_70, na.rm=TRUE)
yr_d <- range(dat$Cumulative_residual_7_70, na.rm=TRUE)
xlim_d <- c(max(0,xr_d[1]-0.05), min(1.02,xr_d[2]+0.05))
ylim_d <- c(max(0,yr_d[1]-0.05), yr_d[2]+0.15)

p_d <- ggplot(dat, aes(Integrated_CE_7_70, Cumulative_residual_7_70,
                       colour=Context, shape=Context, linetype=Context)) +
  geom_point(size=1.15, alpha=0.24, stroke=0) +
  geom_smooth(method="lm", formula=y~x, se=FALSE, linewidth=0.7) +
  scale_colour_manual(values=context_cols) +
  scale_shape_manual(values=context_shapes) +
  scale_linetype_manual(values=context_ltypes) +
  coord_cartesian(xlim=xlim_d, ylim=ylim_d, clip="off") +
  labs(
    x="Integrated compensation efficiency",
    y="Time-integrated absolute residual"
  ) +
  annotate(
    "label",
    x=Inf, y=Inf, label=lab_d,
    hjust=1.03, vjust=1.03, size=2.45,
    label.size=0, fill="white"
  ) +
  theme_fig +
  theme(legend.position="bottom")

# ---------- Assemble ----------
fig7 <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(tag_levels="a") &
  theme(plot.tag=element_text(family="Arial", size=11, face="bold"))
fig7

ggsave("Figure7_V3_R.pdf", fig7, width=174, height=150, units="mm", device=cairo_pdf)

