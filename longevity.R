library(longevity)
library(ggplot2)
library(dplyr, warn.conflicts = FALSE)

par(pch = 20, bty = "l")
data(japanese2, package = "longevity")
data(dutch, package = "longevity")

# Table 1
female_japanese <- japanese2 |>
  dplyr::filter(gender == "female") |>
  dplyr::select(!gender)
knitr::kable(
  tidyr::pivot_wider(
    female_japanese,
    names_from = bcohort,
    values_from = count
  ),
  booktabs = TRUE,
  longtable = FALSE,
  centering = TRUE,
  row.names = FALSE,
  linesep = "",
  caption = "Death count by birth cohort and age band for female Japanese."
)


# Extract sampling frame from attributes of data set
yr_samp <- year(attr(x = dutch, which = "sampling_frame"))
# Preprocess data for analysis
dutch1 <- dutch |>
  subset(!is.na(ndays)) |>
  # Remove interval censored data for the time being
  mutate(
    time = ndays / 365.25, # age at death
    time2 = time,
    # min/max age to be included in sampling frame
    ltrunc = ltrunc / 365.25,
    rtrunc = rtrunc / 365.25,
    event = 1
  ) |> # observed failure time (event=1)
  subset(time > 98) |>
  select(time, time2, ltrunc, rtrunc, event, gender, byear)
# Subset all interval censored, interval truncated records
dutch2 <- dutch |>
  subset(is.na(ndays)) |>
  mutate(
    time2 = ceiling_date(
      dmy(paste("01-", dmonth, "-", dyear)),
      unit = "month"
    ) -
      1 -
      dmy(paste("01-01-", byear)),
    time = dmy(paste("01-", dmonth, "-", dyear)) -
      dmy(paste("31-12-", byear)),
    ltrunc = dmy(paste("01-01-1986")) - dmy(paste("31-12-", byear)),
    rtrunc = dmy(paste("31-12-2015")) - dmy(paste("01-01-", byear))
  ) |>
  select(time, time2, ltrunc, rtrunc, gender, byear) |>
  # Transform data from days to years for interpretability
  mutate(
    time = as.numeric(time) / 365.25, # lower censoring limit
    time2 = as.numeric(time2) / 365.25, # upper censoring limit
    ltrunc = as.numeric(ltrunc) / 365.25, # min age to be included
    rtrunc = as.numeric(rtrunc) / 365.25, # max age to be included
    event = 3
  ) |> # interval censoring (event=3)
  subset(time > 98) # subset exceedances above 98 years
# Combine databases
dutch_data <- rbind(dutch1, dutch2)


## Dutch preview
set.seed(2025)
# Sample some observations
icens_d <- sample(which(dutch_data$event == 1), 3)
obs_d <- sample(which(dutch_data$event == 3), 2)

## Table 2
knitr::kable(
  dutch_data[c(obs_d, icens_d), ],
  digits = 2,
  booktabs = TRUE,
  longtable = FALSE,
  centering = TRUE,
  row.names = FALSE,
  linesep = "",
  caption = "Sample of five Dutch records, formatted so that the inputs match the function arguments used by the package. Columns give the age in years at death (or plausible interval), lower and upper truncation bounds giving minimum and maximum age for inclusion, an integer indicating the type of censoring, gender and birth year."
)


# Keep only non-empty cells
japanese2 <- japanese2[japanese2$count > 0, ]
# Define arguments that are recycled
japanese2$rtrunc <- 2020 -
  as.integer(substr(japanese2$bcohort, 1, 4))
# The line above extracts the earliest year of the birth cohort
# Create a list with all arguments common to package functions
args_japan <- with(
  japanese2,
  list(
    time = age, # lower censoring bound
    time2 = age + 1L, # upper censoring bound
    event = 3, # define interval censoring
    type = "interval2",
    rtrunc = rtrunc, # right truncation limit
    weights = count
  )
) # counts as weights


## Model-comparison
thresh <- 108
model0 <- fit_elife(arguments = args_japan, thresh = thresh, family = "exp")

(model1 <- fit_elife(arguments = args_japan, thresh = thresh, family = "gomp"))


anova(model1, model0)
# Information criteria
c("exponential" = BIC(model0), "Gompertz" = BIC(model1))


nEW <- 179
library(lubridate)
set.seed(2023)
# First observation from 1856, so maximum age for truncation is 55 years
ub <- pgamma(q = 55 * 365.25, scale = 9.945 * 365.25, shape = 1.615)
# Sample right truncated record
bdate_EW <- lubridate::ymd("1910-12-31") -
  qgamma(p = runif(n = nEW) * ub, scale = 9.945 * 365.25, shape = 1.615)
# Obtain truncation bounds given sampling frame
ltrunc_EW <- pmax(0, (ymd("1968-01-01") - bdate_EW) / 365.25 - 110)
rtrunc_EW <- as.numeric(ymd("2020-12-31") - bdate_EW) / 365.25 - 110
sim_EW <- longevity::samp_elife(
  n = nEW,
  scale = 1.2709, # parameters obtained from fitting IDL data for E&W
  shape = -0.0234,
  lower = ltrunc_EW, # smallest age for left truncation limit
  upper = rtrunc_EW, # maximum age
  family = "gp", # generalized Pareto
  type2 = "ltrt"
) # left truncated right truncated
ewsim4 <- data.frame(
  time = sim_EW,
  ltrunc = ltrunc_EW,
  rtrunc = rtrunc_EW
)


## Bootstrap comparison
set.seed(2022)
# Count of unique right truncation limit
db_rtrunc <- aggregate(
  count ~ rtrunc,
  FUN = "sum",
  data = japanese2,
  subset = age >= thresh
)
B <- 1000L # Number of bootstrap replications
boot_anova <- numeric(length = B)
boot_gof <- numeric(length = B)
for (b in seq_len(B - 1L)) {
  boot_samp <- # Generate bootstrap sample
    do.call(
      rbind, #merge data frames
      apply(db_rtrunc, 1, function(x) {
        # for each rtrunc and count
        count <- table(
          #tabulate count
          floor(
            #round down
            samp_elife(
              # sample right truncated exponential
              n = x["count"],
              scale = model0$par,
              family = "exp", #null model
              upper = x["rtrunc"] - thresh,
              type2 = "ltrt"
            )
          )
        )
        data.frame(
          # return data frame
          count = as.integer(count),
          rtrunc = as.numeric(x["rtrunc"]) - thresh,
          eage = as.integer(names(count))
        )
      })
    )
  boot_mod0 <- # Fit null model to bootstrap sample
    with(
      boot_samp,
      fit_elife(
        time = eage,
        time2 = eage + 1L,
        rtrunc = rtrunc,
        type = "interval",
        event = 3,
        family = "exp",
        weights = count
      )
    )
  boot_mod1 <- # Fit alternative model to bootstrap sample
    with(
      boot_samp,
      fit_elife(
        time = eage,
        time2 = eage + 1L,
        rtrunc = rtrunc,
        type = "interval",
        event = 3,
        family = "gomp",
        weights = count
      )
    )
  boot_anova[b] <- deviance(boot_mod0) -
    deviance(boot_mod1)
}
# Add original statistic
boot_anova[B] <- deviance(model1) - deviance(model0)
# Bootstrap p-value
(pval <- rank(boot_anova)[B] / B)


## Threshold diagnostic tools
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
# Threshold sequence
u <- 100:110
# Threshold stability plot
tstab(
  arguments = args_japan,
  family = "gp",
  method = "profile",
  which.plot = "shape",
  thresh = u
)
# Northrop-Coleman diagnostic based on score tests
nu <- length(u) - 1L
nc_score <- nc_test(arguments = c(args_japan, list(thresh = u)))
score_plot <- plot(nc_score)
graphics.off()


## Probability-probability and quantile-quantile plots
## for generalized Pareto model fitted above age 105
## to Dutch data.
fit_dutch <- fit_elife(
  arguments = dutch_data,
  event = 3,
  type = "interval2",
  family = "gp",
  thresh = 105,
  export = TRUE
)
par(mfrow = c(1, 2))
plot(fit_dutch, which.plot = c("pp", "qq"))


## Probability-probability (left) and generalized Pareto quantile-quantile (right) plots for the simulated England and Wales supercentenarian data.
fit_EW <- with(
  ewsim4,
  longevity::fit_elife(
    time = time,
    ltrunc = ltrunc,
    rtrunc = rtrunc,
    family = "gp",
    export = TRUE
  )
)
plots <- plot(
  fit_EW,
  which.plot = c("pp", "qq"),
  plot.type = "ggplot",
  plot = FALSE
)
library(patchwork)
plots[[1]] + plots[[2]]

## Bootstrap goodness-of-fit
# Create contingency table with observations
# grouping all above ubound
get_observed_table <-
  function(data, ubound) {
    # Generate all combinations
    df_combo <-
      expand.grid(
        eage = 0:ubound,
        rtrunc = unique(data$rtrunc)
      )
    df_combo$count <- 0
    # Merge data frame
    # (ensures some empty category appear)
    df_count <- data |>
      dplyr::select(eage, rtrunc, count) |>
      dplyr::full_join(df_combo, by = c("eage", "rtrunc", "count")) |>
      dplyr::mutate(eage = pmin(eage, ubound)) |>
      # Regroup observations above ubound
      dplyr::count(eage, rtrunc, wt = count, name = "count")
  }
# Compute expected counts, conditioning
# on number per right truncation limits
get_expected_count <-
  function(model, data, ubound) {
    data$prob <-
      (pelife(
        q = ifelse(data$eage == ubound, data$rtrunc, data$eage + 1),
        family = model$family,
        scale = model$par
      ) -
        pelife(q = data$eage, family = model$family, scale = model$par)) /
      pelife(q = data$rtrunc, family = model$family, scale = model$par)
    count_rtrunc <- data |>
      dplyr::group_by(rtrunc) |>
      dplyr::summarize(tcount = sum(count), .groups = "keep")
    merge(data, count_rtrunc, by = "rtrunc") |>
      dplyr::transmute(observed = count, expected = tcount * prob)
  }
# Compute chi-square statistic
chisquare_stat <- function(data) {
  with(data, sum((observed - expected)^2 / expected))
}
# Upper bound for pooling
ubound <- 5L
boot_gof <- numeric(length = B)
for (b in seq_len(B - 1L)) {
  # Generate bootstrap sample
  boot_samp <-
    do.call(
      rbind, #merge data frames
      apply(db_rtrunc, 1, function(x) {
        # for each rtrunc and count
        count <- table(
          #tabulate count
          floor(
            #round down
            # sample right truncated exponential
            samp_elife(
              n = x["count"],
              scale = model0$par,
              family = "exp", #null model
              upper = x["rtrunc"] - thresh,
              type2 = "ltrt"
            )
          )
        )
        # return data frame
        data.frame(
          count = as.integer(count),
          rtrunc = as.numeric(x["rtrunc"]) - thresh,
          eage = as.integer(names(count))
        )
      })
    )
  # Fit null model to bootstrap sample
  boot_mod0 <-
    with(
      boot_samp,
      fit_elife(
        time = eage,
        time2 = eage + 1L,
        rtrunc = rtrunc,
        type = "interval",
        event = 3,
        family = "exp",
        weights = count
      )
    )
  ctab <- get_observed_table(boot_samp, ubound = ubound)
  boot_gof[b] <-
    chisquare_stat(
      data = get_expected_count(model = boot_mod0, data = ctab, ubound = ubound)
    )
}
# Add original statistic
db_origin <-
  aggregate(
    count ~ rtrunc + age,
    FUN = "sum",
    data = japanese2,
    subset = age >= thresh
  ) |>
  dplyr::mutate(eage = age - thresh)
ctab <- get_observed_table(db_origin, ubound = ubound)
boot_gof[B] <-
  chisquare_stat(
    data = get_expected_count(model = model0, data = ctab, ubound = ubound)
  )
# Bootstrap p-value
boot_gof_pval <- rank(boot_gof)[B] / B


## Covariate test
print(
  test_elife(
    arguments = args_japan,
    thresh = 110,
    family = "gp",
    covariate = japanese2$gender
  )
)


## MLE and 95% confidence intervals for endpoint
# Create grid of threshold values
thresholds <- 105:110
# Grid of values at which to evaluate profile
psi <- seq(120, 200, length.out = 101)
# Calculate the profile for the endpoint
# of the generalized Pareto at each threshold
endpt_tstab <- do.call(
  endpoint.tstab,
  args = c(
    args_japan,
    list(psi = psi, thresh = thresholds, plot = FALSE)
  )
)
# Compute corresponding confidence intervals
profile <- endpoint.profile(
  arguments = c(args_japan, list(thresh = 110, psi = psi))
)
# Plot point estimates and confidence intervals
g1 <- autoplot(endpt_tstab, plot = FALSE, ylab = "lifespan (in years)")
# Plot the profile curve with cutoffs for conf. int. for 110
g2 <- autoplot(profile, plot = FALSE)
patchwork::wrap_plots(g1, g2)


## Hazard and posterior distribution
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
threshold <- 108
# Note that we cannot have an argument 'arguments' in 'ru'
post_samp <- rust::ru(
  logf = lpost_elife,
  weights = args_japan$weights,
  rtrunc = args_japan$rtrunc,
  event = 3,
  type = "interval2",
  time = args_japan$time,
  time2 = args_japan$time2,
  thresh = threshold,
  family = "gp",
  trans = "BC",
  n = 1000,
  d = 2,
  init = c(1.67, -0.08),
  lower = c(0, -1)
)

plot(post_samp, xlab = "scale", ylab = "shape", bty = "l")
age <- threshold + seq(0.1, 10, length.out = 101)
haz_samp <- apply(X = post_samp$sim_vals, MARG = 1, FUN = function(par) {
  helife(x = age - threshold, scale = par[1], shape = par[2], family = "gp")
})
# Functional boxplots
fbox <- fda::fbplot(
  fit = haz_samp,
  x = age,
  xlim = range(age),
  ylim = c(0, 2.8),
  xaxs = "i",
  yaxs = "i",
  xlab = "age",
  ylab = "hazard (in years)",
  outliercol = "gray90",
  color = "gray60",
  barcol = "black",
  bty = "l"
)

## Empirical CDF
ecdf <- np_elife(arguments = args_japan, thresh = 108)
# Summary statistics, accounting for censoring
round(summary(ecdf), digits = 2)
# Plots of fitted parametric model and nonparametric CDF
model_gp <- fit_elife(
  arguments = args_japan,
  thresh = 108,
  family = "gp",
  export = TRUE
)
# ggplot2 plots, wrapped to display side by side
patchwork::wrap_plots(
  autoplot(
    model_gp, # fitted model
    plot = FALSE, # return list of ggplots
    which.plot = c("dens", "cdf"),
    breaks = seq(0L, 8L, by = 1L) # set bins for histogram
  )
)
