# 06_Figure6.R
source("R/00_setup.R")
long_file <- file.path(D_REC,"provider_replacement_long_Rrerun.csv")
if(file.exists(long_file)){
  prlong <- read_csv(long_file,show_col_types=FALSE)
  a0 <- prlong %>% filter(Day==70)
  b0 <- prlong %>%
    mutate(Unplanted_defined_class=Unplanted_class) %>%
    group_by(Unplanted_defined_class,Context,Day) %>%
    summarise(Median_abs_function_deviation=median(Absolute_function_deviation,na.rm=TRUE),
              Median_provider_replacement_BC=median(Provider_replacement_BC,na.rm=TRUE),.groups="drop")
  c0 <- prlong %>% group_by(Context,Day) %>%
    summarise(Median_treatment_to_reference_BC=median(Provider_replacement_BC,na.rm=TRUE),
              Median_reference_temporal_drift_BC=median(Reference_drift_from_7d_BC,na.rm=TRUE),.groups="drop")
} else {
  a0 <- read_csv(file.path(D_PANEL,"Main_Figures","Figure6",
                           "Fig6a_70d_function_deviation_provider_replacement.csv"),show_col_types=FALSE)
  b0 <- read_csv(file.path(D_PANEL,"Main_Figures","Figure6",
                           "Fig6b_class_function_vs_provider_trajectory.csv"),show_col_types=FALSE)
  c0 <- read_csv(file.path(D_PANEL,"Main_Figures","Figure6",
                           "Fig6c_provider_replacement_vs_reference_drift.csv"),show_col_types=FALSE)
}
d_file <- file.path(D_REC,"K00360_nasB_provider_composition_Rrerun.csv")
if(file.exists(d_file)){
  d0 <- read_csv(d_file,show_col_types=FALSE)
} else {
  d0 <- read_csv(file.path(D_PANEL,"Main_Figures","Figure6",
                           "Fig6d_K00360_nasB_provider_composition.csv"),show_col_types=FALSE)
}

# a: near-reference threshold = 0.25; strong provider replacement = 0.50.
a <- ggplot(a0,aes(Absolute_function_deviation,Provider_replacement_BC))+
  geom_vline(xintercept=.25,linetype=2,colour="grey60")+
  geom_hline(yintercept=.50,linetype=2,colour="grey60")+
  geom_point(alpha=.25,size=.8)+facet_wrap(~Context)+
  labs(x=expression("|log"[2]*"FC| at 70 d"),
       y="Provider replacement (Bray-Curtis)")+theme_pub()

# b: class-specific time dynamics. Use two aligned rows rather than a dual axis.
bf <- b0 %>% select(Unplanted_defined_class,Context,Day,
                    Median_abs_function_deviation,Median_provider_replacement_BC) %>%
  pivot_longer(starts_with("Median_"),names_to="Metric",values_to="Value") %>%
  mutate(Metric=recode(Metric,
    Median_abs_function_deviation="Function deviation",
    Median_provider_replacement_BC="Provider replacement"))
b <- ggplot(bf,aes(Day,Value,colour=Context,linetype=Context))+
  geom_line(linewidth=.6)+geom_point(size=1.3)+
  facet_grid(Metric~Unplanted_defined_class,scales="free_y")+
  scale_colour_manual(values=c(Unplanted="#4C6F8A",`Tomato-planted`="#2A9D8F"))+
  scale_x_continuous(breaks=c(7,14,28,42,56,70))+
  labs(x="Days after fumigation",y=NULL)+theme_pub()

# c: treatment->reference replacement versus reference temporal drift
cl <- c0 %>%
  select(Context,Day,Median_treatment_to_reference_BC,Median_reference_temporal_drift_BC) %>%
  pivot_longer(starts_with("Median_"),names_to="Comparison",values_to="BC") %>%
  mutate(Comparison=recode(Comparison,
    Median_treatment_to_reference_BC="Fumigated -> time-matched reference",
    Median_reference_temporal_drift_BC="Reference temporal drift"))
c <- ggplot(cl,aes(Day,BC,colour=Comparison,linetype=Comparison))+
  geom_line(linewidth=.65)+geom_point(size=1.4)+facet_wrap(~Context)+
  scale_x_continuous(breaks=c(7,14,28,42,56,70))+
  labs(x="Days after fumigation",y="Median provider Bray-Curtis")+theme_pub()

# d: K00360/nasB, Unplanted exemplar only.
# Keep the largest contributors individually and collapse the tail for editable plotting.
dd <- d0 %>% filter(Context=="Unplanted") %>%
  group_by(MAG_provider) %>% summarise(mx=max(Provider_contribution_fraction),.groups="drop") %>%
  arrange(desc(mx)) %>% slice_head(n=8) %>% pull(MAG_provider)
d <- d0 %>% filter(Context=="Unplanted") %>%
  mutate(Provider=if_else(MAG_provider %in% dd,MAG_provider,"Other")) %>%
  group_by(Day,Group,Provider) %>% summarise(Fraction=sum(Provider_contribution_fraction),.groups="drop") %>%
  ggplot(aes(factor(Day),Fraction,fill=Provider))+
  geom_col(width=.8)+facet_wrap(~Group,ncol=1)+
  labs(x="Days after fumigation",y="Provider contribution fraction")+theme_pub()

p <- a/b/(c|d)+plot_annotation(tag_levels="a")
save_pdf(p,"Figure6_R.pdf",width=7.2,height=8.8)
