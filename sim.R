library(dplyr)
library(survival)

set.seed(428361)

n_experiments = 1000

n_high_risk <- 5001
n_low_risk <- 5000
x_diff_mean <- 0.1
x_sd <- 1.0
lambda0 <- 1.0

n <- n_low_risk + n_high_risk

coeff_cohort <- rep(NA, n_experiments)
coeff_subcohort <- rep(NA, n_experiments)
coeff_wsubcohort <- rep(NA, n_experiments)

for (i in 1:n_experiments) {
    print(i)
    cohort <- tibble(.rows=n)
    cohort$risk_cat <- c(rep(0, n_low_risk), rep(1, n_high_risk))
    cohort$x <- rnorm(n, x_diff_mean*cohort$risk_cat, x_sd)
    cohort$lambda <- lambda0 * exp(cohort$x)

    cohort$time_to_event <- rexp(n, cohort$lambda)
    cohort$censor_time <- runif(n, .1, .3)
    cohort$time <- pmin(cohort$censor_time, cohort$time_to_event)
    cohort$event <- cohort$censor_time >= cohort$time_to_event

    # Cohort
    m_cohort <- coxph(Surv(time, event) ~ x, cohort)
    coeff_cohort[i] <- m_cohort$coefficients[1]

    # Subcohort
    stopifnot(sum(cohort$event == 0) > sum(cohort$event == 1))
    case_idx <- which(cohort$event == 1)
    control_idx <- which(cohort$event == 0)
    subcohort_idx <- c(case_idx, control_idx[1:sum(cohort$event == 1)])
    subcohort <- cohort[subcohort_idx, ]
    m_subcohort <- coxph(Surv(time, event) ~ x, subcohort)
    coeff_subcohort[i] <- m_subcohort$coefficients[1]

    # Weighting
    subcohort$w <- ifelse(
        subcohort$event == 0,
        sum(cohort$event == 0) / sum(subcohort$event == 0),
        sum(cohort$event == 1) / sum(subcohort$event == 1))
    m_wsubcohort <- coxph(Surv(time, event) ~ x, subcohort, weights=w)
    coeff_wsubcohort[i] <- m_wsubcohort$coefficients[1]
}
boxplot(list(cohort=coeff_cohort, subcohort=coeff_subcohort, weighted_subcohort=coeff_wsubcohort))
