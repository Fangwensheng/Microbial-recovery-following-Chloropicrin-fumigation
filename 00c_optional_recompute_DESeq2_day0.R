# 00c_optional_recompute_DESeq2_day0.R
# Optional full regeneration of the 0-d disturbed-KO screen.
source("R/00_setup.R")
if(!requireNamespace("DESeq2",quietly=TRUE)){
  stop("Install Bioconductor DESeq2 before running this script.")
}
suppressPackageStartupMessages(library(DESeq2))

ko <- read_csv(file.path(D_RAW,"KO_all_counts.csv"),show_col_types=FALSE)
meta <- read_csv(file.path(D_RAW,"metadata_metagenome.csv"),show_col_types=FALSE)

count <- as.data.frame(ko)
rownames(count) <- count[[1]]
count[[1]] <- NULL
count <- as.matrix(count)
storage.mode(count) <- "integer"

m0 <- meta %>% filter(Week==0,Treatment %in% c("CK","CP")) %>%
  mutate(condition=factor(Treatment,levels=c("CK","CP")))
count0 <- count[,m0$SampleID,drop=FALSE]

keep <- rowSums(count0)>=10 & rowSums(count0>=2)>=2
count0f <- count0[keep,,drop=FALSE]
stopifnot(nrow(count0f)==6099)

dds <- DESeqDataSetFromMatrix(countData=count0f,
                              colData=as.data.frame(m0 %>% column_to_rownames("SampleID")),
                              design=~condition)
dds <- DESeq(dds)
res <- results(dds,contrast=c("condition","CP","CK"),alpha=0.10)
out <- as.data.frame(res) %>% rownames_to_column("KO") %>%
  mutate(
    Disturbed_primary=!is.na(padj) & padj<0.10 & abs(log2FoldChange)>=0.5,
    High_confidence=!is.na(padj) & padj<0.05 & abs(log2FoldChange)>=1
  )
write_csv(out,file.path(D_REC,"DESeq2_day0_all_KO_Rrerun.csv"))
write_csv(out %>% filter(Disturbed_primary),
          file.path(D_REC,"DESeq2_disturbed_KO_primary_Rrerun.csv"))
message("Prefilter KOs: ",nrow(out),"; primary disturbed: ",sum(out$Disturbed_primary))
