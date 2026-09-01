# ============================================================
# Power simulation: Responsive vs Unresponsive
# Windowed cross-lagged correlation + paired t-test
# ============================================================

set.seed(123)

# Parameters
n_sim <- 1000
n_dyads <- 75

n_points <- 12
window_size <- 5
lag <- 1

rho_responsive <- 0.25
rho_unresponsive <- 0.05

alpha <- 0.05


# ------------------------------------------------------------
# Simulate one dyad in one condition
# ------------------------------------------------------------

simulate_dyad <- function(n, lag, rho) {

  caregiver <- rnorm(n)
  dog <- numeric(n)

  for (t in (lag + 1):n) {

    dog[t] <-
      rho * caregiver[t - lag] +
      sqrt(1 - rho^2) * rnorm(1)
  }

  data.frame(
    caregiver = caregiver,
    dog = dog
  )
}


# ------------------------------------------------------------
# Windowed cross-lagged correlation
# ------------------------------------------------------------

windowed_clc <- function(df, window_size, lag) {

  n_windows <- floor(
    (nrow(df) - window_size) / lag
  )

  corrs <- numeric(n_windows)

  for (i in 1:n_windows) {

    start <- (i - 1) * lag + 1
    end <- start + window_size - 1

    caregiver_win <-
      df$caregiver[start:(end - lag)]

    dog_win <-
      df$dog[(start + lag):end]

    corrs[i] <- cor(
      caregiver_win,
      dog_win
    )
  }

  mean(corrs, na.rm = TRUE)
}


# ------------------------------------------------------------
# Simulate one complete experiment
# ------------------------------------------------------------

simulate_one <- function() {

  responsive_scores <- numeric(n_dyads)
  unresponsive_scores <- numeric(n_dyads)

  for (d in 1:n_dyads) {

    # Responsive condition
    responsive_data <- simulate_dyad(
      n = n_points,
      lag = lag,
      rho = rho_responsive
    )

    responsive_scores[d] <- windowed_clc(
      responsive_data,
      window_size,
      lag
    )


    # Unresponsive condition
    unresponsive_data <- simulate_dyad(
      n = n_points,
      lag = lag,
      rho = rho_unresponsive
    )

    unresponsive_scores[d] <- windowed_clc(
      unresponsive_data,
      window_size,
      lag
    )
  }


  # Paired t-test
  test <- t.test(
    responsive_scores,
    unresponsive_scores,
    paired = TRUE
  )

  test$p.value < alpha
}


# ------------------------------------------------------------
# Estimate power
# ------------------------------------------------------------

results <- replicate(
  n_sim,
  simulate_one()
)

power <- mean(results)

mcse <- sqrt(
  power * (1 - power) / n_sim
)


cat(
  "Sample size:", n_dyads,
  "\nEstimated power:", round(power, 3),
  "\nMonte Carlo SE:", round(mcse, 3),
  "\n"
)
