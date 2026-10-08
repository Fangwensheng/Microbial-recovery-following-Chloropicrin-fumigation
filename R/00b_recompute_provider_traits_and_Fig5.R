# 00b_recompute_provider_traits_and_Fig5.R
source("R/00_setup.R")

# This script replaces the previously undocumented five provider functional-trait
# categories with a transparent, gene-level, annotation-based classification.
# These are operational trait proxies, not universal KEGG categories.
# Classification uses KO_description + Level3 only (not Swiss-Prot free text),
# making the rules easier to audit and reproduce.
# A gene may belong to >1 trait if it independently satisfies >1 rule.

gene <- read_tsv(file.path(D_MAG,"MAG-KO-基因.txt"), show_col_types=FALSE, progress=FALSE)
rep <- read_tsv(file.path(D_MAG,"NR_MAG_representatives.txt"), show_col_types=FALSE)
traits0 <- read_tsv(file.path(D_MAG,"NR_MAG_traits.txt"), show_col_types=FALSE)
rrn <- read_tsv(file.path(D_MAG,"NR_MAG_rrnDB_copy_number.txt"), show_col_types=FALSE)
ab <- read_tsv(file.path(D_MAG,"NR_MAG_abundance_78samples.txt"), show_col_types=FALSE)
copy <- read_tsv(file.path(D_MAG,"NR_MAG_KO_copy_number_representative.txt"),
                 show_col_types=FALSE, progress=FALSE)
meta <- read_csv(file.path(D_RAW,"metadata_metagenome.csv"),show_col_types=FALSE)

rules <- tribble(
  ~Trait_category, ~regex,
  "Transporter", "abc transporters|phosphotransferase system|\\btransporter\\b|\\btransport system\\b|\\bpermease\\b|\\bsymporter\\b|\\bantiporter\\b|\\befflux\\b|\\buptake system\\b",
  "OxStress", "oxidative stress|superoxide dismutase|catalase|peroxiredoxin|glutathione peroxidase|glutathione reductase|alkyl hydroperoxide reductase|organic hydroperoxide resistance|rubrerythrin|methionine sulfoxide reductase",
  "DNArepair", "base excision repair|nucleotide excision repair|mismatch repair|homologous recombination|non-homologous end-joining|dna repair|excision repair|mismatch repair|recombination protein|photolyase|dna ligase|endonuclease iii|endonuclease iv|uracil-dna glycosylase|apurinic/apyrimidinic|dna damage|alkyltransferase",
  "Motility", "bacterial chemotaxis|flagellar assembly|chemotaxis|chemotactic|flagell|motility|twitching motility|swarming",
  "Dormancy", "sporulation|spore|germination|dormancy|small acid-soluble spore|dipicolinate|germinant"
)

g <- gene %>%
  filter(Genome %in% rep$Representative_MAG) %>%
  mutate(KO_description=replace_na(KO_description,"-"),
         Level3=replace_na(Level3,"-"),
         class_text=str_to_lower(paste(KO_description,Level3,sep=" | ")))

membership <- map_dfr(seq_len(nrow(rules)), function(i){
  g %>%
    filter(str_detect(class_text, regex(rules$regex[i], ignore_case=TRUE))) %>%
    transmute(`Gene ID`, Genome, KO, KEGG_name, KO_description, Level3,
              Trait_category=rules$Trait_category[i])
})
write_csv(membership,file.path(D_REC,"provider_trait_gene_membership_Rrerun.csv"))

counts <- membership %>%
  count(Genome,Trait_category,name="N") %>%
  complete(Genome=rep$Representative_MAG,
           Trait_category=rules$Trait_category, fill=list(N=0)) %>%
  pivot_wider(names_from=Trait_category,values_from=N)

traits <- traits0 %>%
  left_join(rrn %>% select(Representative_MAG,rrn_copy_number,Match_level),
            by="Representative_MAG") %>%
  left_join(counts,by=c("Representative_MAG"="Genome")) %>%
  mutate(
    Transporter_density=1000*Transporter/Gene_number,
    OxStress_density=1000*OxStress/Gene_number,
    DNArepair_density=1000*DNArepair/Gene_number,
    Motility_density=1000*Motility/Gene_number,
    Dormancy_density=1000*Dormancy/Gene_number,
    CAZyme_density=1000*CAZyme_total/Gene_number,
    Genome_size_Mb=Genome_size_bp/1e6,
    tRNA_per_Mb=tRNA_number/Genome_size_Mb
  )
write_csv(traits,file.path(D_REC,"NR_MAG_traits_reclassified_Rrerun.csv"))

# CK0 provider weights p_ik = mean(A_i,CK0) * C_ik / sum_i(...)
ck0 <- meta %>% filter(Treatment=="CK",Week==0) %>% pull(SampleID)
ab_ck0 <- ab %>% mutate(CK0_mean=rowMeans(across(all_of(ck0)))) %>%
  select(Representative_MAG,CK0_mean)

driver_old <- read_csv(file.path(D_PANEL,"Main_Figures","Figure5",
                                 "Fig5b_KO_level_driver_raw.csv"),
                       show_col_types=FALSE)
kos <- driver_old$KO

# reshape only the 1,177 disturbed KOs used in Fig.5
copy_long <- copy %>%
  select(Representative_MAG,all_of(kos)) %>%
  pivot_longer(-Representative_MAG,names_to="KO",values_to="KO_copy") %>%
  filter(KO_copy>0) %>%
  left_join(ab_ck0,by="Representative_MAG") %>%
  mutate(w=CK0_mean*KO_copy) %>%
  group_by(KO) %>%
  mutate(p=w/sum(w)) %>%
  ungroup() %>%
  left_join(traits %>% select(Representative_MAG,
                              Transporter_density,OxStress_density,DNArepair_density,
                              Motility_density,Dormancy_density),
            by="Representative_MAG")

wtraits <- copy_long %>% group_by(KO) %>%
  summarise(
    Weighted_Transporter_density=sum(p*Transporter_density),
    Weighted_OxStress_density=sum(p*OxStress_density),
    Weighted_DNArepair_density=sum(p*DNArepair_density),
    Weighted_Motility_density=sum(p*Motility_density),
    Weighted_Dormancy_density=sum(p*Dormancy_density),
    .groups="drop"
  )

driver <- driver_old %>%
  select(-Weighted_Transporter_density,-Weighted_OxStress_density,
         -Weighted_DNArepair_density,-Weighted_Motility_density,
         -Weighted_Dormancy_density) %>%
  left_join(wtraits,by="KO")
write_csv(driver,file.path(D_REC,"Fig5_driver_master_reclassified_Rrerun.csv"))

# Univariate Spearman
pred <- tribble(
~Predictor,~Variable,~Block,
"Initial disturbance magnitude","Initial_disturbance_magnitude","Initial state",
"Initial direction","Initial_direction_binary","Initial state",
"Baseline KO abundance","CK0_baseline_KO_abundance","Initial state",
"Provider richness","Provider_richness","Functional architecture",
"Effective provider diversity","Effective_provider_diversity","Functional architecture",
"Provider phylogenetic breadth","Provider_phylogenetic_breadth","Functional architecture",
"KO genomic copy number","Weighted_KO_copy","Functional architecture",
"Provider GC content","Weighted_GC_percent","Provider traits",
"Provider genome size","Genome_size_Mb","Provider traits",
"Provider gene density","Weighted_Gene_density","Provider traits",
"Provider tRNA density","Weighted_tRNA_per_Mb","Provider traits",
"rrn copy number","Weighted_rrn_copy","Provider traits",
"CAZyme density","Weighted_CAZyme_density","Provider traits",
"Transporter density","Weighted_Transporter_density","Provider traits",
"Oxidative-stress density","Weighted_OxStress_density","Provider traits",
"DNA-repair density","Weighted_DNArepair_density","Provider traits",
"Motility/chemotaxis","Weighted_Motility_density","Provider traits",
"Dormancy/sporulation","Weighted_Dormancy_density","Provider traits"
)

univ_one <- function(yvar,context){
  z <- pmap_dfr(pred,function(Predictor,Variable,Block){
    ok <- complete.cases(driver[,c(Variable,yvar)])
    ct <- cor.test(driver[[Variable]][ok],driver[[yvar]][ok],
                   method="spearman",exact=FALSE)
    tibble(Context=context,Predictor=Predictor,Variable=Variable,Block=Block,
           N=sum(ok),rho=unname(ct$estimate),P=ct$p.value)
  })
  z %>% mutate(FDR=p.adjust(P,"BH"))
}
univ <- bind_rows(univ_one("Unplanted_residual","Unplanted"),
                  univ_one("Tomato_planted_residual","Tomato-planted"))
write_csv(univ,file.path(D_REC,"Fig5_univariate_reclassified_Rrerun.csv"))

# partial Spearman via rank residualization
partial_spearman <- function(x,y,covs,df){
  z <- df %>% select(all_of(c(x,y,covs))) %>% drop_na() %>%
    mutate(across(everything(),rank,ties.method="average"))
  rx <- resid(lm(reformulate(covs,response=x),data=z))
  ry <- resid(lm(reformulate(covs,response=y),data=z))
  ct <- cor.test(rx,ry,method="pearson")
  tibble(N=nrow(z),partial_rho=unname(ct$estimate),P=ct$p.value)
}
covars <- c("Initial_disturbance_magnitude","CK0_baseline_KO_abundance",
            "Initial_direction_binary","Weighted_completeness","Weighted_contamination")
partial <- crossing(
  Context=c("Unplanted","Tomato-planted"),
  Predictor=c("Provider richness","Effective provider diversity",
              "Provider phylogenetic breadth")
) %>%
  mutate(
    x=recode(Predictor,
      "Provider richness"="Provider_richness",
      "Effective provider diversity"="Effective_provider_diversity",
      "Provider phylogenetic breadth"="Provider_phylogenetic_breadth"),
    y=if_else(Context=="Unplanted","Unplanted_residual","Tomato_planted_residual")
  ) %>%
  pmap_dfr(function(Context,Predictor,x,y){
    partial_spearman(x,y,covars,driver) %>%
      mutate(Context=Context,Predictor=Predictor,.before=1)
  })
write_csv(partial,file.path(D_REC,"Fig5_partial_spearman_Rrerun.csv"))

# ---------------- Nested CV ridge ----------------
# Exact outer/inner fold assignments are bundled so the analysis is reproducible
# independent of R's random-number generator.
outer <- read_csv(file.path(D_REC,"CV_outer_fold_assignments.csv"),show_col_types=FALSE)
inner <- read_csv(file.path(D_REC,"CV_inner_fold_assignments.csv"),show_col_types=FALSE)
lambda_grid <- 10^seq(-4,4,length.out=81)

impute_scale <- function(Xtrain,Xtest){
  med <- apply(Xtrain,2,function(x) median(x,na.rm=TRUE))
  for(j in seq_len(ncol(Xtrain))){
    Xtrain[is.na(Xtrain[,j]),j] <- med[j]
    Xtest[is.na(Xtest[,j]),j] <- med[j]
  }
  mu <- colMeans(Xtrain)
  s <- apply(Xtrain,2,sd); s[!is.finite(s)|s==0] <- 1
  list(train=sweep(sweep(Xtrain,2,mu,"-"),2,s,"/"),
       test=sweep(sweep(Xtest,2,mu,"-"),2,s,"/"))
}
ridge_predict <- function(Xtrain,ytrain,Xtest,lambda){
  zz <- impute_scale(Xtrain,Xtest)
  ybar <- mean(ytrain); yc <- ytrain-ybar
  beta <- solve(crossprod(zz$train)+lambda*diag(ncol(zz$train)),
                crossprod(zz$train,yc))
  as.numeric(ybar + zz$test %*% beta)
}
blocks <- list(
  `Initial state`=c("Initial_disturbance_magnitude","Initial_direction_binary",
                    "CK0_baseline_KO_abundance"),
  `Functional architecture`=c("Provider_richness","Effective_provider_diversity",
                              "Provider_phylogenetic_breadth","Weighted_KO_copy"),
  `Provider traits`=c("Weighted_GC_percent","Genome_size_Mb","Weighted_Gene_density",
                      "Weighted_tRNA_per_Mb","Weighted_rrn_copy",
                      "Weighted_CAZyme_density","Weighted_Transporter_density",
                      "Weighted_OxStress_density","Weighted_DNArepair_density",
                      "Weighted_Motility_density","Weighted_Dormancy_density")
)
blocks$Combined <- c(blocks[[1]],blocks[[2]],blocks[[3]])

nested_cv <- function(yvar,context){
  allres <- list(); ff <- list()
  for(bn in names(blocks)){
    cols <- blocks[[bn]]
    X <- as.matrix(driver[,cols]); storage.mode(X) <- "double"
    y <- driver[[yvar]]
    fold_r2 <- numeric(10); selected <- numeric(10)
    for(of in 1:10){
      train <- which(outer$OuterFold!=of); test <- which(outer$OuterFold==of)
      inmap <- inner %>% filter(OuterFold==of) %>% select(RowIndex0,InnerFold)
      inf <- setNames(inmap$InnerFold,inmap$RowIndex0+1)
      mse <- numeric(length(lambda_grid))
      for(li in seq_along(lambda_grid)){
        sse <- 0; nn <- 0
        for(k in 1:10){
          val <- train[inf[as.character(train)]==k]
          tr  <- train[inf[as.character(train)]!=k]
          pr <- ridge_predict(X[tr,,drop=FALSE],y[tr],X[val,,drop=FALSE],lambda_grid[li])
          sse <- sse + sum((y[val]-pr)^2); nn <- nn+length(val)
        }
        mse[li] <- sse/nn
      }
      lam <- lambda_grid[which.min(mse)]
      pr <- ridge_predict(X[train,,drop=FALSE],y[train],X[test,,drop=FALSE],lam)
      fold_r2[of] <- 1-sum((y[test]-pr)^2)/sum((y[test]-mean(y[train]))^2)
      selected[of] <- lam
    }
    allres[[bn]] <- tibble(Context=context,Block=bn,N=nrow(driver),
                           Mean_heldout_R2=mean(fold_r2),
                           SD_outer_fold_R2=sd(fold_r2),
                           Median_selected_lambda=median(selected))
    ff[[bn]] <- tibble(Context=context,Block=bn,OuterFold=1:10,
                       Heldout_R2=fold_r2,Selected_lambda=selected)
  }
  list(summary=bind_rows(allres),folds=bind_rows(ff))
}
cvU <- nested_cv("Unplanted_residual","Unplanted")
cvT <- nested_cv("Tomato_planted_residual","Tomato-planted")
write_csv(bind_rows(cvU$summary,cvT$summary),
          file.path(D_REC,"Fig5_nestedCV_reclassified_Rrerun.csv"))
write_csv(bind_rows(cvU$folds,cvT$folds),
          file.path(D_REC,"Fig5_nestedCV_outerfolds_Rrerun.csv"))

print(univ %>% filter(Predictor %in% c("Transporter density","Oxidative-stress density",
                                      "DNA-repair density","Motility/chemotaxis",
                                      "Dormancy/sporulation")))
print(bind_rows(cvU$summary,cvT$summary))
