from pymbar import timeseries
# determine equilibrated region
[t, g, Neff_max] = timeseries.detectEquilibration(A_t)
# extract equilibrated region
A_t_equilibrated = A_t[t:]


