# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

performance_test1 <- function(input) {
    .Call(`_UnifProb_performance_test1`, input)
}

performance_test2 <- function(input) {
    .Call(`_UnifProb_performance_test2`, input)
}

performance_test3 <- function(input) {
    .Call(`_UnifProb_performance_test3`, input)
}

performance_test4 <- function(n) {
    invisible(.Call(`_UnifProb_performance_test4`, n))
}

compute_prob_fft <- function(m, gt, ht, diff_t, debug = FALSE) {
    .Call(`_UnifProb_compute_prob_fft`, m, gt, ht, diff_t, debug)
}

set_plan_flag <- function(flag) {
    invisible(.Call(`_UnifProb_set_plan_flag`, flag))
}

set_fft_rounding <- function(x) {
    invisible(.Call(`_UnifProb_set_fft_rounding`, x))
}

set_fft_min_size <- function(x) {
    invisible(.Call(`_UnifProb_set_fft_min_size`, x))
}

