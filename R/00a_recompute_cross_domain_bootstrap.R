# 00a_recompute_cross_domain_bootstrap.R
source("R/00_setup.R")
# Hierarchical bootstrap accounts for both biological-replicate uncertainty and
# the fact that each domain is represented by two profiles:
# taxonomy = Bacteria + Fungi; function = KO + COG.
# In each bootstrap iteration:
# 1) sample the two profiles within each domain with replacement;
# 2) for each selected profile x day, resample biological replicate distances;
# 3) compute the replicate mean at each day and normalize by its own 7-d mean;
# 4) fit relative_distance ~ Day/10;
# 5) contrast = mean(function slopes) - mean(taxonomy slopes).
#
# Main point estimate is therefore interpretable as change in normalized
# treatment-to-contemporaneous-reference distance per 10 d.

distance_file <- file.path(D_REC,"Aitchison_distance_all_profiles_Rrerun.csv")
if(file.exists(distance_file)){
  dat <- read_csv(distance_file,show_col_types=FALSE) %>%
    filter(Distance_type=="Treatment_to_reference")
} else {
  tax <- read_csv(file.path(D_PANEL,"Main_Figures","Figure2",
                            "Fig2e-f_Aitchison_distance_to_time_matched_reference.csv"),
                  show_col_types=FALSE)
  fun <- read_csv(file.path(D_PANEL,"Main_Figures","Figure3",
                            "Fig3c-d_KO_COG_Aitchison_distance.csv"),
                  show_col_types=FALSE)
  dat <- bind_rows(tax,fun) %>% filter(Distance_type=="Treatment_to_reference")
}

boot_context <- function(context, B=50000, seed=20260817){
  set.seed(seed)
  days <- c(7,14,28,42,56,70)
  tax_pool <- c("Bacteria","Fungi")
  fun_pool <- c("KO","COG")
  contrast <- numeric(B)
  tax_slope <- numeric(B)
  fun_slope <- numeric(B)

  one_profile_slope <- function(profile){
    m <- sapply(days,function(dd){
      v <- dat %>% filter(Context==context, Profile==profile, Day==dd) %>% pull(Aitchison_distance)
      mean(sample(v, length(v), replace=TRUE))
    })
    rel <- m/m[1]
    coef(lm(rel ~ I(days/10)))[2]
  }

  for(b in seq_len(B)){
    tx <- sample(tax_pool,2,replace=TRUE)
    fn <- sample(fun_pool,2,replace=TRUE)
    tax_slope[b] <- mean(vapply(tx, one_profile_slope, numeric(1)))
    fun_slope[b] <- mean(vapply(fn, one_profile_slope, numeric(1)))
    contrast[b] <- fun_slope[b]-tax_slope[b]
  }

  p_emp <- 2*min((sum(contrast>=0)+1)/(B+1),
                 (sum(contrast<=0)+1)/(B+1))
  p_emp <- min(p_emp,1)
  summ <- tibble(
    Context=context, B=B, Seed=seed,
    Taxonomy_slope_mean=mean(tax_slope),
    Function_slope_mean=mean(fun_slope),
    Function_minus_taxonomy=mean(contrast),
    CI_low=quantile(contrast,0.025),
    Median=median(contrast),
    CI_high=quantile(contrast,0.975),
    Empirical_two_sided_P=p_emp
  )
  draws <- tibble(Context=context, Iteration=seq_len(B),
                  Slope_contrast_function_minus_taxonomy=contrast)
  list(summary=summ, draws=draws)
}

b1 <- boot_context("Unplanted")
b2 <- boot_context("Tomato-planted")
summary_out <- bind_rows(b1$summary,b2$summary)
draws_out <- bind_rows(b1$draws,b2$draws)

write_csv(summary_out,file.path(D_REC,"cross_domain_bootstrap_summary_Rrerun.csv"))
write_csv(draws_out,file.path(D_REC,"cross_domain_bootstrap_draws_Rrerun.csv"))
print(summary_out)
