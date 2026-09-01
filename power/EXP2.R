# ============================================================
# Power simulation:
# Responsive vs Unresponsive physiology coordination
# ============================================================

set.seed(123)

# -----------------------------
# Parameters
# -----------------------------

n_sim <- 1000          # number of simulated experiments

n_points <- 12         # number of time points per condition
window_size <- 5       # window size
lag <- 1               # lag

alpha <- 0.05          # significance level

# Assumed true physiological coordination
rho_responsive <- 0.20
rho_unresponsive <- 0.05


# ============================================================
# 1. Simulate physiology data for one dyad in one condition
# ============================================================

simulate_dyad <- function(n, lag, rho) {

  # Caregiver physiology
  caregiver <- rnorm(n)

  # Dog physiology
  dog <- numeric(n)

  for (t in (lag + 1):n) {

    dog[t] <- rho * caregiver[t - lag] +
      sqrt(1 - rho^2) * rnorm(1)

  }

  data.frame(
    caregiver = caregiver,
    dog = dog
  )
}


# ============================================================
# 2. Calculate windowed cross-lagged correlation
# ============================================================

windowed_clc <- function(df, window_size, lag) {

  # Number of windows
  n_windows <- nrow(df) - window_size + 1

  corrs <- numeric(n_windows)

  for (i in 1:n_windows) {

    start <- i
    end <- start + window_size - 1

    caregiver_window <-
      df$caregiver[start:(end - lag)]

    dog_window <-
      df$dog[(start + lag):end]

    corrs[i] <- cor(
      caregiver_window,
      dog_window
    )
  }

  # Average coordination for this dyad
  mean(corrs, na.rm = TRUE)
}


# ============================================================
# 3. Run ONE simulated experiment
# ============================================================

simulate_one_experiment <- function(n_dyads) {

  responsive_corrs <- numeric(n_dyads)
  unresponsive_corrs <- numeric(n_dyads)

  for (d in 1:n_dyads) {

    # -------------------------
    # Responsive condition
    # -------------------------

    responsive_data <- simulate_dyad(
      n = n_points,
      lag = lag,
      rho = rho_responsive
    )

    responsive_corrs[d] <- windowed_clc(
      responsive_data,
      window_size = window_size,
      lag = lag
    )


    # -------------------------
    # Unresponsive condition
    # -------------------------

    unresponsive_data <- simulate_dyad(
      n = n_points,
      lag = lag,
      rho = rho_unresponsive
    )

    unresponsive_corrs[d] <- windowed_clc(
      unresponsive_data,
      window_size = window_size,
      lag = lag
    )
  }


  # -------------------------
  # Compare conditions
  # paired t-test because
  # -------------------------

  test <- t.test(
    responsive_corrs,
    unresponsive_corrs,
    paired = TRUE
  )

  test$p.value < alpha
}


# ============================================================
# 4. Estimate power for ONE sample size
# ============================================================

estimate_power <- function(n_dyads) {

  results <- replicate(
    n_sim,
    simulate_one_experiment(n_dyads)
  )

  mean(results)
}


# ============================================================
# 5. Try different sample sizes
# ============================================================

sample_sizes <- c(
  20,
  25,
  30,
  35,
  40,
  45,
  50,
  60
)

power_results <- sapply(
  sample_sizes,
  estimate_power
)


# ============================================================
# 6. Results
# ============================================================

results_table <- data.frame(
  n_dyads = sample_sizes,
  power = round(power_results, 3)
)

print(results_table)


# ============================================================
# 7. Plot
# ============================================================

plot(
  sample_sizes,
  power_results,
  type = "b",
  ylim = c(0, 1),
  xlab = "Number of dyads",
  ylab = "Power",
  main = "Power: Responsive vs Unresponsive"
)

abline(
  h = 0.80,
  lty = 2
)
