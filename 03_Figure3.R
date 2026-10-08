# 03_Figure3.R
source("R/00_setup.R")
alpha <- read_csv(file.path(D_PANEL,"Main_Figures","Figure3",
                            "Fig3a-b_KO_COG_alpha_raw.csv"),show_col_types=FALSE)
dist_file <- file.path(D_REC,"Aitchison_distance_all_profiles_Rrerun.csv")
rel_file <- file.path(D_REC,"relative_distance_7d_equals_1_Rrerun.csv")
if(file.exists(dist_file)){
  dist <- read_csv(dist_file,show_col_types=FALSE) %>% filter(Profile %in% c("KO","COG"))
} else {
  dist <- read_csv(file.path(D_PANEL,"Main_Figures","Figure3",
                             "Fig3c-d_KO_COG_Aitchison_distance.csv"),show_col_types=FALSE)
}
if(file.exists(rel_file)){
  rel <- read_csv(rel_file,show_col_types=FALSE)
} else {
  rel <- read_csv(file.path(D_PANEL,"Main_Figures","Figure3",
                            "Fig3e-f_relative_distance_7d_equals_1.csv"),show_col_types=FALSE)
}

sa <- alpha %>% group_by(Panel,Metric,Treatment,Day) %>%
  summarise(mean=mean(Value),se=sd(Value)/sqrt(n()),.groups="drop")
make_alpha <- function(panel,y){
  ggplot(sa %>% filter(Panel==panel),
         aes(Day,mean,colour=Treatment,shape=Treatment,linetype=Treatment))+
    geom_line(linewidth=.6)+geom_point(size=1.6)+
    geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=2,linewidth=.35)+
    scale_colour_manual(values=treatment_cols)+scale_shape_manual(values=treatment_shapes)+
    scale_linetype_manual(values=treatment_ltypes)+
    scale_x_continuous(breaks=c(7,14,28,42,56,70))+
    labs(x="Days after fumigation",y=y)+theme_pub()
}
a<-make_alpha("Fig3a","KO richness"); b<-make_alpha("Fig3b","COG richness")

sdist <- dist %>% group_by(Profile,Context,Day,Distance_type) %>%
  summarise(mean=mean(Aitchison_distance),se=sd(Aitchison_distance)/sqrt(n()),.groups="drop")
make_dist <- function(profile){
  z<-sdist %>% filter(Profile==profile)
  ggplot(z %>% filter(Distance_type=="Treatment_to_reference"),
         aes(Day,mean,colour=Context,linetype=Context))+
    geom_line(linewidth=.7)+geom_point(size=1.6)+
    geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=2,linewidth=.35)+
    scale_colour_manual(values=c(Unplanted="#4C6F8A",`Tomato-planted`="#2A9D8F"))+
    scale_linetype_manual(values=c(Unplanted="solid",`Tomato-planted`="dashed"))+
    scale_x_continuous(breaks=c(7,14,28,42,56,70))+
    labs(x="Days after fumigation",y="Aitchison distance\nto time-matched reference")+theme_pub()
}
c<-make_dist("KO"); d<-make_dist("COG")

# e-f: relative distance is NOT an absolute recovery index.
# It is each profile's mean treatment-to-reference Aitchison distance divided by
# that profile's own 7-d mean. Thus 1 = the 7-d disturbance-state distance.
make_rel <- function(context){
  ggplot(rel %>% filter(Context==context),
         aes(Day,Relative_distance_7d_eq_1,colour=Profile,shape=Profile))+
    geom_hline(yintercept=1,linetype=2,colour="grey60")+
    geom_line(linewidth=.7)+geom_point(size=1.7)+
    scale_colour_manual(values=profile_cols)+
    scale_x_continuous(breaks=c(7,14,28,42,56,70))+
    labs(x="Days after fumigation",
         y="Relative distance to time-matched reference\n(7 d = 1)")+
    theme_pub()
}
e<-make_rel("Unplanted"); f<-make_rel("Tomato-planted")
p<-(a|b)/(c|d)/(e|f)+plot_annotation(tag_levels="a")
save_pdf(p,"Figure3.pdf",width=7.2,height=8.0)
