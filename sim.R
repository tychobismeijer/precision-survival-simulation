library(dplyr)
requireNamespace('MatchIt')
requireNamespace('multipleNCC')
library(purrr)
library(survival)

set.seed(428361)

n_experiments = 10
risk_group_proportion = 0.5
coeff_risk = 1.0

wb_intercept = 7.03 # estimated from Dutch nation-wide cohort
wb_shape=0.86 # estimated from Dutch nation-wide cohort

x_diff_mean <- 1.0
x_sd <- 1.0

months_per_year=12
fu_up_to = 2011

# From NKR
incidence <- tribble(
    ~year, ~n,
    1989,  325,
    1990,  338,
    1991,  496,
    1992,  632,
    1993,  700,
    1994,  757,
    1995,  818,
    1996,  792,
    1997,  927,
    1998,  980,
    1999, 1020,
    2000, 1104,
    2001, 1120,
    2002, 1038,
    2003, 1142,
    2004, 1300)


n <- sum(incidence$n)

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
    }
)


coeff = list()
for (method in names(methods)) {
    coeff[[method]] <- rep(NA, n_experiments)
}

for (i in 1:n_experiments) {
    message("Experiment ", i)

    # Cohort
    cohort <- pmap_dfr(incidence, function(year, n) {
        n_high_risk = floor(n * risk_group_proportion)
        n_low_risk = n - n_high_risk
        tibble(
            follow_up_time = (fu_up_to - year - runif(n)) * months_per_year,
            year=year,
            risk_cat = rep(c(0, 1), times=c(n_low_risk, n_high_risk)) 
    )})
    cohort$scale <- exp(wb_intercept + coeff_risk * scale(cohort$risk_cat))
    cohort$x <- rnorm(nrow(cohort), x_diff_mean*cohort$risk_cat, x_sd)

    cohort$time_to_event <- rweibull(nrow(cohort), scale=cohort$scale, shape=wb_shape)
    cohort$time <- pmin(cohort$follow_up_time, cohort$time_to_event)
    cohort$event <- cohort$follow_up_time >= cohort$time_to_event

    # Case-control
    matching <- MatchIt::matchit(
        event ~ year,
        data=cohort,
        exact='year',
        ratio=1)
    subcohort_idx <- which(!is.na(matching$subclass))
    subcohort <- cohort[subcohort_idx, ]
    cohort$sampled <- F
    cohort$sampled[subcohort_idx] <- T

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
boxplot(cohort$x ~ cohort$event)
dev.off()
