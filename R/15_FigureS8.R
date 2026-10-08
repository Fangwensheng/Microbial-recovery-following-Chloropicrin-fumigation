# 13_FigureS8.R
source("R/00_setup.R")
all <- read_csv(file.path(D_PANEL,"Supplementary","FigureS8",
                          "FigS8_all_contrasts_and_descriptive_stats.csv"),show_col_types=FALSE)
plant <- read_csv(file.path(D_PANEL,"Supplementary","FigureS8",
                            "FigS8_plant_raw_with_individual_grades.csv"),show_col_types=FALSE)
disease <- read_csv(file.path(D_PANEL,"Supplementary","FigureS8",
                              "FigS8f_wilt_disease_index_descriptive.csv"),show_col_types=FALSE)

# a functional/pathogen qPCR heatmap, Tomato-planted context
qa <- all %>% filter(Section=="Functional/pathogen qPCR",Context=="Tomato-planted") %>%
  mutate(Metric=factor(Metric,levels=c("nifH","AOB","AOA","nirS","nirK","nosZ",
                                      "phoD","pomA","mcrA","fungcbl","FS","RS")))
a <- ggplot(qa,aes(factor(Day),Metric,fill=log2_mean_ratio))+
  geom_tile(colour="white",linewidth=.15)+
  scale_fill_gradient2(low="#3B6FB6",mid="white",high="#C94C4C",midpoint=0)+
  theme_pub()+labs(x="Days after fumigation",y=NULL,fill=expression(log[2]*"(TCP/TCK)"))

# b selected soil/process heatmap
selproc <- c("铵态氮","硝态氮","表观硝化速率","反硝化速率","N2O","CO2","MBC-2","MBN","MBP","有效磷")
qb <- all %>% filter(Section=="Soil/process",Context=="Tomato-planted",Metric %in% selproc) %>%
  mutate(Metric=factor(Metric,levels=selproc))
b <- ggplot(qb,aes(factor(Day),Metric,fill=log2_mean_ratio))+
  geom_tile(colour="white",linewidth=.15)+
  scale_fill_gradient2(low="#3B6FB6",mid="white",high="#C94C4C",midpoint=0)+
  theme_pub()+labs(x="Days after fumigation",y=NULL,fill=expression(log[2]*"(TCP/TCK)"))

# c N-cycle gene contrasts
c0 <- qa %>% filter(as.character(Metric) %in% c("AOB","AOA","nirS","nirK","nosZ"))
c <- ggplot(c0,aes(Day,log2_mean_ratio,colour=Metric))+
  geom_hline(yintercept=0,linetype=2,colour="grey60")+
  geom_line(linewidth=.6)+geom_point(size=1.3)+theme_pub()+
  labs(x="Days after fumigation",y=expression(log[2]*"(TCP/TCK)"))

# d N pools / process
d0 <- qb %>% filter(as.character(Metric) %in% c("铵态氮","硝态氮","表观硝化速率","反硝化速率","N2O"))
d <- ggplot(d0,aes(Day,log2_mean_ratio,colour=Metric))+
  geom_hline(yintercept=0,linetype=2,colour="grey60")+
  geom_line(linewidth=.6)+geom_point(size=1.3)+theme_pub()+
  labs(x="Days after fumigation",y=expression(log[2]*"(TCP/TCK)"))

# e pathogens
e0 <- qa %>% filter(as.character(Metric) %in% c("FS","RS"))
e <- ggplot(e0,aes(Day,log2_mean_ratio,colour=Metric))+
  geom_hline(yintercept=0,linetype=2,colour="grey60")+
  geom_line(linewidth=.65)+geom_point(size=1.4)+theme_pub()+
  labs(x="Days after fumigation",y=expression(log[2]*"(TCP/TCK)"))

# f wilt disease index - descriptive treatment-level index only; no inferential test.
f <- ggplot(disease,aes(Day,Wilt_disease_index,colour=Treatment,shape=Treatment,linetype=Treatment))+
  geom_line(linewidth=.65)+geom_point(size=1.5)+
  scale_colour_manual(values=treatment_cols)+scale_shape_manual(values=treatment_shapes)+
  scale_linetype_manual(values=treatment_ltypes)+theme_pub()+
  labs(x="Days after fumigation",y="Wilt disease index")

# g shoot dry weight
g0 <- plant %>% group_by(Treatment,Day) %>%
  summarise(mean=mean(`植株干重g`,na.rm=TRUE),se=sd(`植株干重g`,na.rm=TRUE)/sqrt(sum(!is.na(`植株干重g`))),.groups="drop")
g <- ggplot(g0,aes(Day,mean,colour=Treatment,shape=Treatment,linetype=Treatment))+
  geom_line(linewidth=.65)+geom_point(size=1.4)+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=2,linewidth=.3)+
  scale_colour_manual(values=treatment_cols)+scale_shape_manual(values=treatment_shapes)+
  scale_linetype_manual(values=treatment_ltypes)+theme_pub()+
  labs(x="Days after fumigation",y="Shoot dry weight (g)")

# h microbial biomass responses
h0 <- qb %>% filter(as.character(Metric) %in% c("MBC-2","MBN","MBP"))
h <- ggplot(h0,aes(Day,log2_mean_ratio,colour=Metric))+
  geom_hline(yintercept=0,linetype=2,colour="grey60")+
  geom_line(linewidth=.65)+geom_point(size=1.3)+theme_pub()+
  labs(x="Days after fumigation",y=expression(log[2]*"(TCP/TCK)"))

save_pdf((a|b)/(c|d)/(e|f)/(g|h)+plot_annotation(tag_levels="a"),
         "FigureS8_R.pdf",width=7.2,height=10.2)
