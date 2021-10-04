library(dplyr)
library(doMC)
library(foreach)
library(ggplot2)
library(MatchIt)
requireNamespace('multipleNCC')
library(purrr)
library(survival)

registerDoMC(8)

n_experiments = 10
risk_group_proportion = 0.5
x_diff_mean <- 1.0
x_sd <- 0.5
lambda0 <- 0.001175099 # per month
hazard_ratio <- 5.0


set.seed(428361)
seeds <- map_int(seq(n_experiments), ~ sample(.Machine$integer.max/2, 1))

generate_cohort <- function() {
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

    # Cohort
    cohort <- pmap_dfr(incidence, function(year, n) {
        n_high_risk = floor(n * risk_group_proportion)
        n_low_risk = n - n_high_risk
        tibble(
            follow_up_time = (fu_up_to - year - runif(n)) * months_per_year,
            year=year,
            risk_cat = rep(c(0, 1), times=c(n_low_risk, n_high_risk)) 
    )})
    cohort$x <- rnorm(nrow(cohort), x_diff_mean*cohort$risk_cat, x_sd)
    cohort$lambda <- lambda0 * exp(log(hazard_ratio) * scale(cohort$risk_cat, scale=F))
    cohort$time_to_event <- rexp(nrow(cohort), cohort$lambda)
    cohort$time <- pmin(cohort$follow_up_time, cohort$time_to_event)
    cohort$event <- cohort$follow_up_time >= cohort$time_to_event
    cohort$sample <- seq(nrow(cohort))

    cohort
}

select_case_control <- function(cohort) {
    matching <- matchit(
        event ~ year,
        data=cohort,
        exact='year',
        ratio=1)
    subcohort_idx <- which(!is.na(matching$subclass))
    subcohort <- cohort[subcohort_idx, ]
    subcohort
}

make_ncc_method <- function(f, ...) {
    function(subcohort, cohort) {
        cohort$sampled <- cohort$sample %in% subcohort$sample
        cohort$sampling_p <- f(
            survtime = cohort$follow_up_time,
            samplestat = case_when(
                !cohort$sampled ~ 0,
                cohort$sampled & cohort$event == 0 ~ 1,
                cohort$sampled & cohort$event == 1 ~ 2),
            ...)
        subcohort_idx <- match(subcohort$sample, cohort$sample)
        subcohort$w <- cohort$sampling_p[subcohort_idx]
        subcohort <- subcohort[subcohort$w > 0,]
        coxph(Surv(time, event) ~ x, subcohort, weights=w)
    }
}

methods <- list(
    cohort = function(subcohort, cohort) {
        coxph(Surv(time, event) ~ x, cohort)
    },
    subcohort = function(subcohort, cohort) {
        coxph(Surv(time, event) ~ x, subcohort)
    },
    ncc_km = make_ncc_method(multipleNCC::KMprob, m=1),
    ncc_gam = make_ncc_method(multipleNCC::GAMprob),
    ncc_glm = make_ncc_method(multipleNCC::GLMprob)
)


coeff_l <- foreach(i=1:n_experiments) %dopar% {
    message("Experiment ", i)
    set.seed(seeds[i])

    cohort <- generate_cohort()
    subcohort <- select_case_control(cohort)

    # Compute coefficients with different methods
    coeff <- list()
    for (method in names(methods)) {
        model <- methods[[method]](subcohort, cohort)
        coeff[[method]] <- model$coefficients[1]
    }
    coeff
}

coeff = list()
for (method in names(methods)) {
    coeff[[method]] <- rep(NA, n_experiments)
}
for (i in seq_along(coeff_l)) {
    for (method in names(methods)) {
        coeff[[method]][i] <- coeff_l[[i]][[method]]
    }
}

pdf('coeff_estimates.pdf')
boxplot(coeff)
abline(h=0.0, lty=2)
abline(h=mean(coeff$cohort), lty=3)
boxplot(lapply(coeff, exp))
abline(h=1.0, lty=2)
abline(h=exp(mean(coeff$cohort)), lty=3)
plot_df <- map_dfr(coeff, ~ tibble(coeff=.), .id='model') %>%
    group_by(model) %>%
    summarize(sd=sd(coeff), coeff=mean(coeff), se=sd/n(), lower=coeff-se, upper=coeff+se)
ggplot(plot_df, aes(x=model, y=exp(coeff), ymin=exp(lower), ymax=exp(upper))) +
    scale_x_discrete(limits=names(coeff)) +
    geom_point() +
    geom_errorbar() +
    scale_y_continuous("Coefficient")
ggplot(plot_df, aes(x=model, y=coeff, ymin=lower, ymax=upper)) +
    scale_x_discrete(limits=names(coeff)) +
    geom_point() +
    geom_errorbar() +
    scale_y_continuous("Hazard Ratio")
dev.off()
