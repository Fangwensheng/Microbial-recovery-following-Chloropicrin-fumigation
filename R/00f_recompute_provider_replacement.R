# 00f_recompute_provider_replacement.R
source("R/00_setup.R")
ab <- read_tsv(file.path(D_MAG,"NR_MAG_abundance_78samples.txt"),show_col_types=FALSE)
copy <- read_tsv(file.path(D_MAG,"NR_MAG_KO_copy_number_representative.txt"),
                 show_col_types=FALSE,progress=FALSE)
meta <- read_csv(file.path(D_RAW,"metadata_metagenome.csv"),show_col_types=FALSE)
traj_file <- file.path(D_REC,"KO_7to70d_oriented_trajectories_Rrerun.csv")
if(file.exists(traj_file)){
  traj <- read_csv(traj_file,show_col_types=FALSE)
} else {
  traj <- read_csv(file.path(D_RAW,"DESeq2_disturbed_KO_7to70d_trajectories.csv"),show_col_types=FALSE)
}
class_file <- file.path(D_REC,"KO_trajectory_classes_Rrerun.csv")
if(file.exists(class_file)){
  classes <- read_csv(class_file,show_col_types=FALSE)
} else {
  classes <- read_csv(file.path(D_RAW,"DESeq2_KO_trajectory_classes.csv"),show_col_types=FALSE)
}
driver_file <- file.path(D_REC,"Fig5_driver_master_reclassified_Rrerun.csv")
if(!file.exists(driver_file)) driver_file <- file.path(D_REC,"Fig5_driver_master_reclassified_v2.csv")
driver <- read_csv(driver_file,show_col_types=FALSE)
kos <- intersect(driver$KO,names(copy))

# matrix: representative MAG x samples
A <- as.matrix(ab %>% select(-ANI_cluster,-Representative_MAG))
rownames(A)<-ab$Representative_MAG
# matrix: representative MAG x KO
C <- as.matrix(copy %>% select(all_of(kos)))
rownames(C)<-copy$Representative_MAG
storage.mode(A)<-"double";storage.mode(C)<-"double"

mean_ab <- function(treatment,wk){
  ids<-meta %>% filter(Treatment==treatment,Week==wk) %>% pull(SampleID)
  rowMeans(A[,ids,drop=FALSE])
}
bray <- function(p,q){
  if(sum(p)==0 || sum(q)==0) return(NA_real_)
  p<-p/sum(p);q<-q/sum(q)
  sum(abs(p-q))/sum(p+q)
}
rows<-list(); ii<-1
for(ctx in c("Unplanted","Tomato-planted")){
  trt<-if(ctx=="Unplanted")"CP" else "TCP"
  ref<-if(ctx=="Unplanted")"CK" else "TCK"
  ref7<-mean_ab(ref,1)
  for(wk in c(1,2,4,6,8,10)){
    at<-mean_ab(trt,wk); ar<-mean_ab(ref,wk)
    for(ko in kos){
      cc<-C[,ko]; prov<-which(cc>0)
      if(length(prov)<2) next
      qt<-at[prov]*cc[prov]; qr<-ar[prov]*cc[prov]; q7<-ref7[prov]*cc[prov]
      if(sum(qt)==0 || sum(qr)==0 || sum(q7)==0) next
      bc<-bray(qt,qr); drift<-bray(qr,q7)
      lf<-traj %>% filter(KO==ko,Mode==if(ctx=="Unplanted")"Natural" else "Plant",Week==wk) %>%
        pull(later_log2FC)
      rows[[ii]]<-tibble(KO=ko,Context=ctx,Day=wk*7,
                         Provider_richness=length(prov),
                         Function_log2FC=lf,
                         Absolute_function_deviation=abs(lf),
                         Provider_replacement_BC=bc,
                         Reference_drift_from_7d_BC=drift);ii<-ii+1
    }
  }
}
out<-bind_rows(rows) %>% left_join(classes %>% select(KO,Natural_class,Plant_class),by="KO") %>%
  rename(Unplanted_class=Natural_class,Tomato_planted_class=Plant_class)
write_csv(out,file.path(D_REC,"provider_replacement_long_Rrerun.csv"))

# ---- 70-d threshold sensitivity and reference-drift summaries ----
thresholds <- c(0.10,0.25,0.50)
thr <- purrr::map_dfr(c("Unplanted","Tomato-planted"), function(ctx){
  z <- out %>% filter(Context==ctx,Day==70)
  purrr::map_dfr(thresholds,function(th){
    ev <- z %>% filter(Absolute_function_deviation<=th)
    tibble(Context=ctx,Near_reference_threshold=th,
           Strong_replacement_threshold=0.50,
           N_evaluable=nrow(ev),
           N_strong=sum(ev$Provider_replacement_BC>=0.50,na.rm=TRUE),
           Percent_strong=100*N_strong/N_evaluable)
  })
})
write_csv(thr,file.path(D_REC,"provider_replacement_threshold_sensitivity_Rrerun.csv"))

refsum <- out %>% group_by(Context,Day) %>%
  summarise(Median_treatment_to_reference_BC=median(Provider_replacement_BC,na.rm=TRUE),
            Median_reference_temporal_drift_BC=median(Reference_drift_from_7d_BC,na.rm=TRUE),
            Percent_treatment_gt_reference_drift=100*mean(Provider_replacement_BC>Reference_drift_from_7d_BC,na.rm=TRUE),
            .groups="drop")
write_csv(refsum,file.path(D_REC,"provider_replacement_reference_drift_summary_Rrerun.csv"))

# ---- K00360/nasB provider contribution example, Unplanted only ----
example_ko <- "K00360"
if(example_ko %in% colnames(C)){
  cc <- C[,example_ko]
  prov <- which(cc>0)
  ex <- list(); jj <- 1
  for(wk in c(1,2,4,6,8,10)){
    for(grp in c("CK","CP")){
      aa <- mean_ab(grp,wk)
      raw_contrib <- aa[prov]*cc[prov]
      frac <- raw_contrib/sum(raw_contrib)
      ex[[jj]] <- tibble(KO=example_ko,Context="Unplanted",Day=wk*7,Group=grp,
                         MAG_provider=rownames(C)[prov],KO_copy_number=cc[prov],
                         Mean_MAG_abundance=aa[prov],Provider_contribution_fraction=frac)
      jj <- jj+1
    }
  }
  write_csv(bind_rows(ex),file.path(D_REC,"K00360_nasB_provider_composition_Rrerun.csv"))
}
