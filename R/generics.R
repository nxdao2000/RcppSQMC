## SMC (particle filter) and use for fixed lag smoothing
setGeneric("pfilter2",function(object,...)standardGeneric("pfilter2"))

## iterated smoothing
setGeneric('is2',function(object,...)standardGeneric("is2"))

HilbertResampler <- function(quasi, parPhi, x, dx, N, W, xh) {
    invisible(.Call('is2_HilbertResampler', PACKAGE = 'RcppSQMC', quasi, parPhi, x, dx, N, W, xh))
}

SQMC_SV <- function(y, dy, dx, T, theta, seed, ns, N, qmc, src, computeExp, parPsi, lik, expx, expx2) {
    invisible(.Call('is2_SQMC_SV', PACKAGE = 'RcppSQMC', y, dy, dx, T, theta, seed, ns, N, qmc, src, computeExp, parPsi, lik, expx, expx2))
}

SQMCBack_SV <- function(y, dy, dx, T, theta, seed, ns, N, Nb, qmc, qmcB, Marg, parPsi, lik, expx, expx2) {
    .Call('is2_SQMCBack_SV', PACKAGE = 'RcppSQMC', y, dy, dx, T, theta, seed, ns, N, Nb, qmc, qmcB, Marg, parPsi, lik, expx, expx2)
}

SQMC2F_SV <- function(y, dy, dx, tstar, T, theta, theta2F, seed, ns, N, Nb, qmc, parPsi, lik, expx, expx2) {
    .Call('is2_SQMC2F_SV', PACKAGE = 'RcppSQMC', y, dy, dx, tstar, T, theta, theta2F, seed, ns, N, Nb, qmc, parPsi, lik, expx, expx2)
}

SQMC_Univ <- function(y, dy, dx, T, theta, seed, ns, N, qmc, src, computeExp, lik, expx, expx2) {
    .Call('is2_SQMC_Univ', PACKAGE = 'RcppSQMC', y, dy, dx, T, theta, seed, ns, N, qmc, src, computeExp, lik, expx, expx2)
}

SQMCBack_Univ <- function(y, dy, dx, T, theta, seed, ns, N, Nb, qmc, qmcb, M, lik, expx, expx2) {
    .Call('is2_SQMCBack_Univ', PACKAGE = 'RcppSQMC', y, dy, dx, T, theta, seed, ns, N, Nb, qmc, qmcb, M, lik, expx, expx2)
}

SQMC_Neuro <- function(y, dy, dx, T, theta, seed, ns, N, qmc, src, computeExp, parPsi, lik, expx, expx2) {
    .Call('is2_SQMC_Neuro', PACKAGE = 'RcppSQMC', y, dy, dx, T, theta, seed, ns, N, qmc, src, computeExp, parPsi, lik, expx, expx2)
}

