install.packages('reticulate')
library(reticulate)

rm(list = ls())
cat("\014")
source("plot_theo_sim.R")
source_python("alpha_estimation.py")

# PARAMETERS TO ADJUST
##########################################################
alpha = 12
beta = 4
n_min = 20
n_max = 2000
n_step = 20
alpha_start = 1.0 # a start point in root-finding Newton-Raphson method
stop_threshold = 5e-3 # the value of the real valued function in Newton-Raphson Method, for which the algorithm
# will stop
derivative_step = 1e-6 # represents a small coefficient close to zero, which will be used to calculate derivatives
integration_step = 1e-2 # the distance between samples, based on which integrates will be calculated
end_integration_step = 1e-2 # In theory, the integration in gamma function is calculated from zero up to infinity,
# but for larger arguments the integrated function tends asymptotically to zero. Thus, the numerical integration
# will stop if a value of the integrated function for a new sample will be less than `end_integration_step` times
# already integrated area.
#########################################################

print(sprintf('Estimating Gamma Distribution parameters with alpha = %f, beta = %f for probes from %d to %d by %d', alpha, beta, n_min, n_max, n_step))

n_vec = seq(n_min, n_max, n_step) # sizes of a probe
n_len = length(n_vec)
est_alpha = numeric(n_len)
est_beta = numeric(n_len)

idx = 1
for (n in n_vec)
{
  probe = array(rgamma(n, shape = alpha, scale = beta))
  est_alpha[idx] = alpha_est_newton_raphson_method(probe,
                                                   alpha_start,
                                                   stop_threshold,
                                                   derivative_step,
                                                   integration_step,
                                                   end_integration_step)
  moment_1 = mean(probe)
  est_beta[idx] = moment_1 / est_alpha[idx]
  print(sprintf('%d / %d', n, n_max))
  print(sprintf('alpha = %f, beta = %f', est_alpha[idx], est_beta[idx]))
  idx = idx + 1
}

plot_theo_sim(n_vec,
              rep(alpha, n_len),
              est_alpha,
              "alpha's estimation",
              "number of samples",
              "alpha")
plot_theo_sim(n_vec,
              rep(beta, n_len),
              est_beta,
              "beta's estimation",
              "number of samples",
              "beta")