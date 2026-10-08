# 00_setup.R
# Run all scripts from the package root:
# setwd(".../Chloropicrin_R_Figure_Reproduction_V1")

required <- c(
  "ggplot2","dplyr","tidyr","readr","readxl","purrr","stringr","tibble",
  "vegan","ape","patchwork","scales","ggrepel","pheatmap","ggalluvial","cluster"
)
miss <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(miss)>0){
  stop("Please install missing CRAN packages: ", paste(miss, collapse=", "))
}
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(readxl)
  library(purrr); library(stringr); library(tibble); library(vegan); library(ape)
  library(patchwork); library(scales); library(ggrepel); library(cluster)
})

ROOT <- normalizePath(".", winslash="/", mustWork=TRUE)
D_PANEL <- file.path(ROOT,"data","panel")
D_RAW   <- file.path(ROOT,"data","raw")
D_MAG   <- file.path(ROOT,"data","mag")
D_REC   <- file.path(ROOT,"data","recomputed")
D_ECO   <- file.path(ROOT,"data","ecology")
OUT     <- file.path(ROOT,"outputs")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

COL_REF <- "#55AEBB"
COL_FUM <- "#E76F51"
COL_BAC <- "#4C78A8"
COL_FUNGI <- "#8C6BB1"
COL_KO <- "#E76F51"
COL_COG <- "#E9A23B"

treatment_cols <- c(CK=COL_REF, CP=COL_FUM, TCK=COL_REF, TCP=COL_FUM)
treatment_shapes <- c(CK=16, CP=16, TCK=15, TCP=15)
treatment_ltypes <- c(CK="solid", CP="solid", TCK="dashed", TCP="dashed")
profile_cols <- c(Bacteria=COL_BAC, Fungi=COL_FUNGI, KO=COL_KO, COG=COL_COG)

theme_pub <- function(base_size=8){
  theme_classic(base_size=base_size, base_family="Arial") +
    theme(
      axis.title=element_text(size=base_size),
      axis.text=element_text(size=base_size-1, colour="black"),
      legend.title=element_blank(),
      legend.text=element_text(size=base_size-1),
      strip.background=element_blank(),
      strip.text=element_text(face="bold", size=base_size),
      plot.title=element_blank(),
      plot.margin=margin(4,4,4,4)
    )
}
theme_set(theme_pub())

mean_se <- function(x){
  x <- x[is.finite(x)]
  tibble(mean=mean(x), se=sd(x)/sqrt(length(x)), n=length(x))
}

clr_matrix <- function(mat, pseudocount=0.5){
  # input: features x samples, non-negative.
  x <- as.matrix(mat)
  storage.mode(x) <- "double"
  x <- x + pseudocount
  lx <- log(x)
  sweep(lx, 2, colMeans(lx), "-")
}

relative_abundance <- function(mat){
  x <- as.matrix(mat)
  storage.mode(x) <- "double"
  sweep(x,2,colSums(x),"/")
}

aitchison_to_time_matched_reference <- function(count_mat, meta, treatment, reference,
                                                pseudocount=1e-6,
                                                convert_to_relative=TRUE){
  # Returns per-treatment-sample distance to the CLR centroid of the
  # contemporaneous reference group at each day.
  x <- as.matrix(count_mat)
  if(convert_to_relative) x <- relative_abundance(x)
  xclr <- clr_matrix(x, pseudocount=pseudocount)
  days <- sort(unique(meta$Day))
  out <- list()
  z <- 1
  for(dd in days){
    ref_ids <- meta %>% filter(Day==dd, Treatment==reference) %>% pull(SampleID)
    trt_ids <- meta %>% filter(Day==dd, Treatment==treatment) %>% pull(SampleID)
    ref_ids <- intersect(ref_ids,colnames(xclr)); trt_ids <- intersect(trt_ids,colnames(xclr))
    if(length(ref_ids)<2 || length(trt_ids)<1) next
    centroid <- rowMeans(xclr[,ref_ids,drop=FALSE])
    for(s in trt_ids){
      out[[z]] <- tibble(Day=dd, SampleID=s,
                         Aitchison_distance=sqrt(sum((xclr[,s]-centroid)^2)))
      z <- z+1
    }
  }
  bind_rows(out)
}

panel_letter <- function(p, lab){
  p + annotate("text", x=-Inf, y=Inf, label=lab, hjust=-0.4, vjust=1.2,
               fontface="bold", size=3)
}

save_pdf <- function(plot, filename, width=7.2, height=5.5){
  ggsave(file.path(OUT,filename), plot, width=width, height=height,
         units="in", device=cairo_pdf)
}
