# 00d_recompute_composition_metrics.R
source("R/00_setup.R")

amp_meta <- read_csv(file.path(D_RAW,"metadata_amplicon.csv"),show_col_types=FALSE) %>%
  mutate(Day=Group*7)
mg_meta <- read_csv(file.path(D_RAW,"metadata_metagenome.csv"),show_col_types=FALSE) %>%
  mutate(Day=Week*7)

read_feature_csv <- function(path){
  z <- read_csv(path,show_col_types=FALSE)
  rn <- z[[1]]
  x <- as.matrix(z[,-1]); rownames(x)<-rn
  storage.mode(x)<-"double"; x
}
bac <- read_feature_csv(file.path(D_RAW,"otutab_Bacteria_Non_Genji.csv"))
fun <- read_feature_csv(file.path(D_RAW,"otutab_Fungi_Non_Genji.csv"))
ko  <- read_feature_csv(file.path(D_RAW,"KO_all_counts.csv"))
cog <- read_feature_csv(file.path(D_RAW,"COG_all_counts.csv"))

distance_long <- function(mat,meta,profile){
  x <- relative_abundance(mat)
  xclr <- clr_matrix(x,pseudocount=1e-6)
  rows <- list(); ii<-1
  for(ctx in c("Unplanted","Tomato-planted")){
    trt <- if(ctx=="Unplanted") "CP" else "TCP"
    ref <- if(ctx=="Unplanted") "CK" else "TCK"
    for(dd in c(7,14,28,42,56,70)){
      trt_ids <- meta %>% filter(Day==dd,Treatment==trt) %>% pull(SampleID) %>% intersect(colnames(xclr))
      ref_ids <- meta %>% filter(Day==dd,Treatment==ref) %>% pull(SampleID) %>% intersect(colnames(xclr))
      cent <- rowMeans(xclr[,ref_ids,drop=FALSE])
      for(s in trt_ids){
        rows[[ii]]<-tibble(Profile=profile,Context=ctx,Day=dd,SampleID=s,
                           Treatment=trt,Reference=ref,
                           Distance_type="Treatment_to_reference",
                           Aitchison_distance=sqrt(sum((xclr[,s]-cent)^2)));ii<-ii+1
      }
      for(s in ref_ids){
        rows[[ii]]<-tibble(Profile=profile,Context=ctx,Day=dd,SampleID=s,
                           Treatment=ref,Reference=ref,
                           Distance_type="Control_internal",
                           Aitchison_distance=sqrt(sum((xclr[,s]-cent)^2)));ii<-ii+1
      }
    }
  }
  bind_rows(rows)
}

dist_all <- bind_rows(
  distance_long(bac,amp_meta,"Bacteria"),
  distance_long(fun,amp_meta,"Fungi"),
  distance_long(ko,mg_meta,"KO"),
  distance_long(cog,mg_meta,"COG")
)
write_csv(dist_all,file.path(D_REC,"Aitchison_distance_all_profiles_Rrerun.csv"))

rel <- dist_all %>% filter(Distance_type=="Treatment_to_reference") %>%
  group_by(Profile,Context,Day) %>% summarise(Mean_Aitchison_distance=mean(Aitchison_distance),.groups="drop") %>%
  group_by(Profile,Context) %>%
  mutate(Relative_distance_7d_eq_1=Mean_Aitchison_distance/
           Mean_Aitchison_distance[Day==7]) %>% ungroup()
write_csv(rel,file.path(D_REC,"relative_distance_7d_equals_1_Rrerun.csv"))

# Day-0 PCoA: raw counts + 0.5 pseudocount -> CLR -> Euclidean.
pcoa_day0 <- function(mat,meta,profile){
  ids <- meta %>% filter(Day==0,Treatment %in% c("CK","CP")) %>% pull(SampleID) %>% intersect(colnames(mat))
  x <- clr_matrix(mat[,ids,drop=FALSE],pseudocount=.5)
  dd <- dist(t(x),method="euclidean")
  pc <- cmdscale(dd,k=2,eig=TRUE,add=FALSE)
  eig <- pc$eig[pc$eig>0]
  pct <- 100*eig/sum(eig)
  tibble(Profile=profile,SampleID=rownames(pc$points),
         PCoA1=pc$points[,1],PCoA2=pc$points[,2],
         PCoA1_percent=pct[1],PCoA2_percent=pct[2]) %>%
    left_join(meta %>% select(SampleID,Treatment,Day),by="SampleID")
}
pc0 <- bind_rows(
  pcoa_day0(bac,amp_meta,"Bacteria"),
  pcoa_day0(fun,amp_meta,"Fungi"),
  pcoa_day0(ko,mg_meta,"KO"),
  pcoa_day0(cog,mg_meta,"COG")
)
write_csv(pc0,file.path(D_REC,"day0_PCoA_coordinates_Rrerun.csv"))
