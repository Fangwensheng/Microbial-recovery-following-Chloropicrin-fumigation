# 00e_recompute_KO_trajectories_and_classes.R
source("R/00_setup.R")
if(!requireNamespace("DESeq2",quietly=TRUE)) stop("Install Bioconductor DESeq2.")
suppressPackageStartupMessages(library(DESeq2))

ko <- read_csv(file.path(D_RAW,"KO_all_counts.csv"),show_col_types=FALSE)
meta <- read_csv(file.path(D_RAW,"metadata_metagenome.csv"),show_col_types=FALSE)
day0 <- read_csv(file.path(D_RAW,"DESeq2_day0_all_KO.csv"),show_col_types=FALSE)
disturbed <- day0 %>% filter(!is.na(padj),padj<.10,abs(log2FoldChange)>=.5) %>%
  select(KO,Day0_log2FC=log2FoldChange)

count <- as.data.frame(ko); rownames(count)<-count[[1]]; count[[1]]<-NULL
count <- as.matrix(count); storage.mode(count)<-"integer"
count <- count[disturbed$KO,,drop=FALSE]

run_contrast <- function(trt,ref,wk,mode){
  mm <- meta %>% filter(Week==wk,Treatment %in% c(ref,trt)) %>%
    mutate(condition=factor(Treatment,levels=c(ref,trt)))
  dds <- DESeqDataSetFromMatrix(countData=count[,mm$SampleID,drop=FALSE],
    colData=as.data.frame(mm %>% column_to_rownames("SampleID")),design=~condition)
  dds <- DESeq(dds,quiet=TRUE)
  rr <- results(dds,contrast=c("condition",trt,ref),independentFiltering=FALSE)
  tibble(KO=rownames(rr),Mode=mode,Week=wk,Day=wk*7,
         later_log2FC=rr$log2FoldChange)
}
later <- bind_rows(
  map_dfr(c(1,2,4,6,8,10),~run_contrast("CP","CK",.x,"Natural")),
  map_dfr(c(1,2,4,6,8,10),~run_contrast("TCP","TCK",.x,"Plant"))
) %>% left_join(disturbed,by="KO") %>%
  mutate(Oriented_relative_effect=later_log2FC/Day0_log2FC)
write_csv(later,file.path(D_REC,"KO_7to70d_oriented_trajectories_Rrerun.csv"))

trapz_mean_abs <- function(day,effect){
  o<-order(day); day<-day[o]; effect<-abs(effect[o])
  sum(diff(day)*(head(effect,-1)+tail(effect,-1))/2)/(max(day)-min(day))
}
desc <- later %>% group_by(KO,Mode) %>%
  summarise(Time_integrated_absolute_residual_7to70d=trapz_mean_abs(Day,Oriented_relative_effect),
            Absolute_residual_at_70d=abs(Oriented_relative_effect[Day==70]),.groups="drop")
write_csv(desc,file.path(D_REC,"KO_recovery_descriptors_Rrerun.csv"))

# Unplanted k=3 core clustering; k selection shown separately by silhouette.
wideN <- later %>% filter(Mode=="Natural") %>% select(KO,Day,Oriented_relative_effect) %>%
  pivot_wider(names_from=Day,values_from=Oriented_relative_effect) %>% arrange(KO)
X <- as.matrix(wideN[,-1])
set.seed(20260817)
km <- kmeans(X,centers=3,nstart=100,iter.max=1000)
cent <- km$centers
core_order <- order(rowMeans(abs(cent)))
lab_core <- setNames(c("Rapid","Gradual/incomplete","Persistent/strengthened"),core_order)
natural_base <- unname(lab_core[as.character(km$cluster)])
natural_class <- ifelse(natural_base=="Rapid" & apply(X,1,min,na.rm=TRUE)<=-.5,
                        "Overshoot/reversal",natural_base)

# Tomato-planted assigned to nearest Unplanted core centroid, then same overshoot rule.
wideP <- later %>% filter(Mode=="Plant") %>% select(KO,Day,Oriented_relative_effect) %>%
  pivot_wider(names_from=Day,values_from=Oriented_relative_effect) %>%
  arrange(match(KO,wideN$KO))
XP <- as.matrix(wideP[,-1])
nearest <- apply(XP,1,function(v) which.min(apply(cent,1,function(cc)sum((v-cc)^2))))
plant_base <- unname(lab_core[as.character(nearest)])
plant_class <- ifelse(plant_base=="Rapid" & apply(XP,1,min,na.rm=TRUE)<=-.5,
                      "Overshoot/reversal",plant_base)

classes <- tibble(KO=wideN$KO,Natural_base_class=natural_base,Natural_class=natural_class,
                  Plant_base_class=plant_base,Plant_class=plant_class) %>%
  left_join(disturbed,by="KO") %>%
  mutate(Day0_direction=if_else(Day0_log2FC>0,"Enriched","Suppressed"))
write_csv(classes,file.path(D_REC,"KO_trajectory_classes_Rrerun.csv"))

# k silhouette diagnostic
sil <- map_dfr(2:7,function(k){
  set.seed(20260817)
  zz<-kmeans(X,centers=k,nstart=100)
  tibble(k=k,Silhouette=mean(cluster::silhouette(zz$cluster,dist(X))[,3]))
})
write_csv(sil,file.path(D_REC,"KO_cluster_silhouette_Rrerun.csv"))
