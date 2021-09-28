library(dplyr)
library(MatchIt)
requireNamespace('multipleNCC')
library(survival)

set.seed(428361)

n_experiments = 100

n_high_risk <- 500
n_low_risk <- 500
x_diff_mean <- 0.1
x_sd <- 1.0
lambda0 <- 1.0

n <- n_low_risk + n_high_risk


methods <- list(
    cohort = function(subcohort, cohort, subcohort_idx) {
        coxph(Surv(time, event) ~ x, cohort)
    },
    subcohort = function(subcohort, cohort, subcohort_idx) {
        coxph(Surv(time, event) ~ x, subcohort)
    },
    naive = function(subcohort, cohort, subcohort_idx) {
        subcohort$w <- ifelse(
            subcohort$event == 0,
            sum(cohort$event == 0) / sum(subcohort$event == 0),
            sum(cohort$event == 1) / sum(subcohort$event == 1))
        coxph(Surv(time, event) ~ x, subcohort, weights=w)
    },
    kaplan_meier = function(subcohort, cohort, subcohort_idx) {
        cohort$sampling_p <- multipleNCC::KMprob(
            survtime = cohort$follow_up_time,
            samplestat = case_when(
                !cohort$sampled ~ 0,
                cohort$sampled & cohort$event == 0 ~ 1,
                cohort$sampled & cohort$event == 1 ~ 2),
            m=1)
        subcohort$w <- cohort$sampling_p[subcohort_idx]
        subcohort <- subcohort[subcohort$w > 0,]
        coxph(Surv(time, event) ~ x, subcohort, weights=w)
    },
    km_exact = function(subcohort, cohort, subcohort_idx) {
        Sur <- summary(survfit(Surv(time, event) ~ 1, cohort),
                       times=subcohort$time)$surv
        devKM <- function(w) {
            subcohort$w <- w / sum(w)
            SurW <- summary(survfit(Surv(time, event) ~ 1, subcohort, weights = subcohort$w),
                            times = subcohort$time)$surv
            mean((SurW-Sur)^2)
        }
        w <- ifelse(
            subcohort$event,
            sum(cohort$event) / sum(subcohort$event),
            sum(1-cohort$event) / sum(1-subcohort$event))
        w <- w / sum(w)
        res <- optim(
            par=w,
            fn=devKM,
            gr = NULL,
            method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent")[4],
            lower = 0,
            upper = 1,
            control = list(trace=3),
            hessian = FALSE)
        subcohort$w <- res$par / sum(res$par)
        subcohort <- subcohort[subcohort$w > 0,]
        coxph(Surv(time, event) ~ x, subcohort, weights=w)
    }

)


coeff = list()
for (method in names(methods)) {
    coeff[[method]] <- rep(NA, n_experiments)
}

for (i in 1:n_experiments) {
    message("Experiment ", i)

    # Cohort
    cohort <- tibble(.rows=n)
    cohort$risk_cat <- c(rep(0, n_low_risk), rep(1, n_high_risk))
    cohort$x <- rnorm(n, x_diff_mean*cohort$risk_cat, x_sd)
    cohort$lambda <- lambda0 * exp(cohort$x)

    cohort$time_to_event <- rexp(n, cohort$lambda)
    cohort$follow_up_time <- runif(n, .1, .3)
    cohort$time <- pmin(cohort$follow_up_time, cohort$time_to_event)
    cohort$event <- cohort$follow_up_time >= cohort$time_to_event

    # Case-control
    cohort$rounded_follow_up_time <- round(cohort$follow_up_time, 2)
    matching <- matchit(
        event ~ rounded_follow_up_time,
        data=cohort,
        exact='rounded_follow_up_time',
        ratio=1)
    subcohort_idx <- which(!is.na(matching$subclass))
    subcohort <- cohort[subcohort_idx, ]
    cohort$sampled <- F
    cohort$sampled[subcohort_idx] <- T

    print(table(cohort$event))
    print(table(subcohort$event))

    saveRDS(cohort, paste0('experiment_', i, '.Rds'))

    # Compute coefficients with different methods
    for (method in names(methods)) {
        model <- methods[[method]](subcohort, cohort, subcohort_idx)
        coeff[[method]][i] <- model$coefficients[1]
    }
}

pdf('coeff_estimates.pdf')
boxplot(coeff)
abline(h=1.0, lty=2)
boxplot(lapply(coeff, exp))
abline(h=exp(1.0), lty=2)
dev.off()
