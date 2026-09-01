set.seed(123)

n_sim <- 5000
n_dyads <- 90
dz <- 0.30
alpha <- 0.05

simulate_one <- function() {

  # Simulate within-dyad Responsive - Unresponsive differences
  difference <- rnorm(
    n_dyads,
    mean = dz,
    sd = 1
  )

  test <- t.test(
    difference,
    mu = 0
  )

  test$p.value < alpha
}

results <- replicate(
  n_sim,
  simulate_one()
)

power <- mean(results)

cat(
  "Sample size =", n_dyads,
  "\nEstimated power =", round(power, 3),
  "\n"
)
