library(cmdstanr)
set_cmdstan_path("~/cmdstan")

Sys.setenv(OMP_NUM_THREADS=1)
Sys.setenv(STAN_NUM_THREADS=1)
Sys.setenv(MKL_NUM_THREADS=1)
Sys.setenv(OPENBLAS_NUM_THREADS=1)

mod <- cmdstan_model("migration.stan")

fit <- mod$sample(
  data = "mig.json",
  chains = 3,
  parallel_chains = 1,
  iter_warmup = 1000,
  iter_sampling = 1000
)

n_pop <- floor(sqrt(length(fit$summary("Migration")$mean)) + 0.5)
cat(fit$summary("Migration")$mean, file = "est_mig.txt", sep = " ")

summary_all <- fit$summary(variables = NULL,
            posterior::default_summary_measures()[1:4],
            quantiles = ~ quantile(., probs = c(0.025, 0.05, 0.95, 0.975)),
            posterior::default_convergence_measures()
)
write.csv(summary_all, "summary_mig.csv", row.names = FALSE)
