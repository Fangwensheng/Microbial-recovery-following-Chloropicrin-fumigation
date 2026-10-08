# 02_Figure2.R
# Important annotation rule:
#   panels a-d = UNIVARIATE three-factor model;
#   panels e-f = COMPOSITION PERMANOVA on the full CLR/Aitchison composition matrix.
# Thus the P values in e-f are NOT tests of the plotted distance curves.

source("R/00_setup_Fig2_S2_v2.R")

raw <- read_csv(file.path(D,"Fig2a-d_abundance_and_Shannon_raw_FIXED.csv"),show_col_types=FALSE) %>%
  mutate(Treatment=factor(Treatment,levels=c("CK","CP","TCK","TCP")))
dist <- read_csv(file.path(D,"Fig2e-f_Aitchison_distance_to_time_matched_reference.csv"),show_col_types=FALSE)
st <- read_csv(file.path(D,"Fig2_displayed_statistics_v2.csv"),show_col_types=FALSE)

sm <- raw %>%
  group_by(Panel,Metric,Treatment,Day) %>%
  summarise(mean=mean(Plot_value,na.rm=TRUE),
            se=sd(Plot_value,na.rm=TRUE)/sqrt(sum(is.finite(Plot_value))),
            .groups="drop")

stat_label <- function(panel){
  z <- st %>% filter(Panel==panel)
  paste0(z$Analysis_type,"\n",z$Term,", ",z$Display_P)
}

make_ts <- function(panel,ylab,metric_lab,ylim=NULL){
  z <- sm %>% filter(Panel==panel)
  p <- ggplot(z,aes(Day,mean,colour=Treatment,shape=Treatment,
                    linetype=Treatment,group=Treatment))+
    geom_line(linewidth=0.62)+
    geom_point(size=1.7)+
    geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=1.6,linewidth=0.34)+
    scale_colour_manual(values=treatment_cols,drop=FALSE)+
    scale_shape_manual(values=treatment_shapes,drop=FALSE)+
    scale_linetype_manual(values=treatment_ltypes,drop=FALSE)+
    scale_x_continuous(breaks=c(7,14,28,42,56,70))+
    labs(x="Days after fumigation",y=ylab)+
    theme_nc()+
    annotate("text",x=7,y=Inf,label=metric_lab,hjust=-0.02,vjust=1.35,size=2.7)+
    annotate("label",x=7,y=-Inf,label=stat_label(panel),
             hjust=-0.02,vjust=-0.25,size=2.2,label.size=NA,
             fill=scales::alpha("white",0.85),lineheight=0.95)
  if(!is.null(ylim)) p <- p + coord_cartesian(ylim=ylim)
  p
}

a <- make_ts("Fig2a",expression(log[10]*" gene copies g"^{-1}*" dry soil"),"16S rRNA",c(7.74,8.33))
b <- make_ts("Fig2b",expression(log[10]*" gene copies g"^{-1}*" dry soil"),"ITS",c(5.15,6.10))
c <- make_ts("Fig2c","Shannon index","Bacteria",c(7.00,8.22))
d <- make_ts("Fig2d","Shannon index","Fungi",c(2.75,3.55))

# e/f: plotted curves are DESCRIPTIVE distances to contemporaneous reference centroids.
# The displayed P value comes from factorial PERMANOVA of the FULL CLR/Aitchison
# composition matrix, not from a model fitted to these two distance curves.
dsm <- dist %>%
  filter(Distance_type=="Treatment_to_reference") %>%
  mutate(Context=factor(Context,levels=c("Unplanted","Tomato-planted"))) %>%
  group_by(Profile,Context,Day) %>%
  summarise(mean=mean(Aitchison_distance),
            se=sd(Aitchison_distance)/sqrt(n()),
            .groups="drop")

make_dist <- function(profile,panel,ylim){
  z <- dsm %>% filter(Profile==profile)
  ggplot(z,aes(Day,mean,group=Context,linetype=Context,shape=Context))+
    geom_line(colour=COL_FUM,linewidth=0.66)+
    geom_point(colour=COL_FUM,size=1.7)+
    geom_errorbar(aes(ymin=mean-se,ymax=mean+se),colour=COL_FUM,width=1.6,linewidth=0.34)+
    scale_linetype_manual(values=c(Unplanted="solid",`Tomato-planted`="dashed"))+
    scale_shape_manual(values=c(Unplanted=16,`Tomato-planted`=15))+
    scale_x_continuous(breaks=c(7,14,28,42,56,70))+
    coord_cartesian(ylim=ylim)+
    labs(x="Days after fumigation",
         y="Aitchison distance to\ntime-matched reference")+
    theme_nc()+
    annotate("text",x=7,y=Inf,label=profile,hjust=-0.02,vjust=1.35,size=2.7)+
    annotate("label",x=7,y=-Inf,label=stat_label(panel),
             hjust=-0.02,vjust=-0.25,size=2.2,label.size=NA,
             fill=scales::alpha("white",0.85),lineheight=0.95)
}

e <- make_dist("Bacteria","Fig2e",c(237,268)) + theme(legend.position="none")
f <- make_dist("Fungi","Fig2f",c(69,82)) + theme(legend.position="none")

p <- ((a|b)/(c|d)/(e|f)) +
  plot_annotation(tag_levels="a",
                  theme=theme(plot.tag=element_text(face="bold",size=9))) +
  plot_layout(guides="collect") &
  theme(legend.position="top") &
  guides(linetype=guide_legend(nrow=1),
         shape=guide_legend(nrow=1),
         colour=guide_legend(nrow=1))

save_both(p,"Figure2",width=7.2,height=8.2)
