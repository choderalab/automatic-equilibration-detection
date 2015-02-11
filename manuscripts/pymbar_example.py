from pymbar import timeseries
[t, g, Neff_max] = timeseries.detectEquilibration(A_t) # determine equilibrated region
A_t_equilibrated = A_t[t:] # extract equilibrated part of the data for analysis

