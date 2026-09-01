# ============================================================
# A PRIORI MONTE CARLO POWER ANALYSIS
#
# X = Attachment avoidance
# M = Responsive caregiving
# Y = Separation-related behaviour
#
# Assumed standardized paths:
# a  = -0.15
# b  = -0.15
# c' =  0.05
#
# Indirect effect = a*b = 0.0225
#
# Target power = .80
# Alpha = .05
# ============================================================


if (!requireNamespace("lavaan", quietly = TRUE)) {
  install.packages("lavaan")
}

if (!requireNamespace("MASS", quietly = TRUE)) {
  install.packages("MASS")
}

library(lavaan)
library(MASS)


set.seed(2026)

target_power <- 0.80
alpha <- 0.05

n_sim <- 1000

# Total sample size = China + Austria combined
candidate_N <- seq(500, 800, by = 50)


a_pop  <- -0.15   # avoidance -> responsive caregiving
b_pop  <- -0.15   # responsive caregiving -> separation problems
cp_pop <-  0.05   # direct effect avoidance -> separation problems

indirect_pop <- a_pop * b_pop

cat("Population indirect effect =", indirect_pop, "\n")


# Correlation between attachment avoidance and anxiety
r_avoid_anx <- 0.30

# For the power analysis, these are treated as nuisance
# variables and assumed to have no population effect.
anx_to_M <- 0
anx_to_Y <- 0

country_to_M <- 0
country_to_Y <- 0





analysis_model <- '

  # Mediator equation
  caregiving ~ a*avoidance +
               anxiety +
               country

  # Outcome equation
  separation ~ cprime*avoidance +
               b*caregiving +
               anxiety +
               country

  # Indirect effect
  indirect := a*b

  # Total effect
  total := cprime + a*b

'


simulate_dataset <- function(N) {

  # ----------------------------------
  # Avoidance and anxiety
  # ----------------------------------

  Sigma_X <- matrix(
    c(
      1, r_avoid_anx,
      r_avoid_anx, 1
    ),
    nrow = 2,
    byrow = TRUE
  )

  X <- MASS::mvrnorm(
    n = N,
    mu = c(0, 0),
    Sigma = Sigma_X
  )

  avoidance <- X[, 1]
  anxiety   <- X[, 2]


  # ----------------------------------
  # Country
  #
  # -1 = China
  # +1 = Austria
  # ----------------------------------

  country <- sample(
    c(-1, 1),
    size = N,
    replace = TRUE
  )


  # ----------------------------------
  # Responsive caregiving
  # M = aX + error
  # ----------------------------------

  var_pred_M <-
    a_pop^2 +
    anx_to_M^2 +
    country_to_M^2 +
    2 * a_pop * anx_to_M * r_avoid_anx

  resid_var_M <- 1 - var_pred_M

  caregiving <-
    a_pop * avoidance +
    anx_to_M * anxiety +
    country_to_M * country +
    rnorm(
      N,
      mean = 0,
      sd = sqrt(resid_var_M)
    )


  # ----------------------------------
  # Calculate covariance structure
  # needed for outcome residual
  # variance
  # ----------------------------------

  cov_avoid_M <-
    a_pop +
    anx_to_M * r_avoid_anx

  cov_anx_M <-
    a_pop * r_avoid_anx +
    anx_to_M

  cov_country_M <- country_to_M


  Sigma_predictors <- matrix(
    c(
      1,
      cov_avoid_M,
      r_avoid_anx,
      0,

      cov_avoid_M,
      1,
      cov_anx_M,
      cov_country_M,

      r_avoid_anx,
      cov_anx_M,
      1,
      0,

      0,
      cov_country_M,
      0,
      1
    ),
    nrow = 4,
    byrow = TRUE
  )


  beta_Y <- c(
    cp_pop,
    b_pop,
    anx_to_Y,
    country_to_Y
  )


  var_pred_Y <- as.numeric(
    t(beta_Y) %*%
      Sigma_predictors %*%
      beta_Y
  )

  resid_var_Y <- 1 - var_pred_Y


  separation <-
    cp_pop * avoidance +
    b_pop * caregiving +
    anx_to_Y * anxiety +
    country_to_Y * country +
    rnorm(
      N,
      mean = 0,
      sd = sqrt(resid_var_Y)
    )


  data.frame(
    avoidance  = avoidance,
    anxiety    = anxiety,
    country    = country,
    caregiving = caregiving,
    separation = separation
  )
}


analyse_one <- function(N) {

  dat <- simulate_dataset(N)


  fit <- try(
    lavaan::sem(
      analysis_model,
      data = dat
    ),
    silent = TRUE
  )


  if (inherits(fit, "try-error")) {
    return(NA)
  }


  if (!lavInspect(fit, "converged")) {
    return(NA)
  }


  pars <- parameterEstimates(fit)


  indirect_row <- pars[
    pars$lhs == "indirect" &
      pars$op == ":=",
  ]


  if (nrow(indirect_row) != 1) {
    return(NA)
  }


  # 1 = indirect effect significant
  # 0 = indirect effect non-significant
  as.numeric(
    indirect_row$pvalue < alpha
  )
}


estimate_power <- function(N) {

  results <- replicate(
    n_sim,
    analyse_one(N)
  )


  valid_results <- results[
    !is.na(results)
  ]


  data.frame(
    N = N,
    power = mean(valid_results),
    valid_replications = length(valid_results),
    failed_replications =
      n_sim - length(valid_results)
  )
}


power_results <- do.call(
  rbind,
  lapply(
    candidate_N,
    function(N) {

      result <- estimate_power(N)

      cat(
        "N =",
        N,
        "| Power =",
        round(result$power, 3),
        "\n"
      )

      result
    }
  )
)


print(power_results)


adequate <- power_results[
  power_results$power >= target_power,
]


if (nrow(adequate) > 0) {

  minimum_N <- adequate$N[1]
  minimum_power <- adequate$power[1]

  cat(
    "\n=====================================\n",
    "Minimum total N reaching 80% power:",
    minimum_N,
    "\nEstimated power:",
    round(minimum_power, 3),
    "\n=====================================\n"
  )

} else {

  cat(
    "\nNo sample size in the tested range",
    "reached 80% power.\n"
  )
}


plot(
  power_results$N,
  power_results$power,
  type = "b",
  pch = 19,
  xlab = "Total sample size (China + Austria)",
  ylab = "Power for indirect effect",
  ylim = c(0, 1)
)

abline(
  h = 0.80,
  lty = 2
)


