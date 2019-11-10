
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
from scipy.stats import genpareto
from statsmodels.distributions.empirical_distribution import ECDF


import sys


data = xr.open_dataset('/Users/albacid/Projects/SERDP/twl_POT.nc')

variable = data['area'].values


# Plot histograms
nbins = 10

fig, ax1 = plt.subplots()
hist_data = ax1.hist(variable, bins=nbins)
ax1.set_xlabel('area (m*h)')
ax1.set_ylabel('nº events')


# fig, ax = plt.subplots(3)
#
# ax[0].hist(data['twl_exceedances'].values, bins=nbins)
# ax[0].set_xlabel('exceedances (m)')
# ax[0].set_ylabel('nº events')
#
# ax[1].hist(data['area'].values, bins=nbins)
# ax[1].set_xlabel('area (m*h)')
# ax[1].set_ylabel('nº events')
#
# ax[2].hist(data['duration'].values, bins=nbins)
# ax[2].set_xlabel('duration (h)')
# ax[2].set_ylabel('nº events')
#
# plt.tight_layout()
# plt.show()
# sys.exit()

# plot empirical PDF
temp = np.histogram(variable, bins=nbins)
twl_pdf = temp[0]/len(variable)

twl_data = []
for i in range(len(temp[1])-1):

    twl_data.append(temp[1][i] + ((temp[1][i+1] - temp[1][i]))/2.0)

# ecdf = ECDF(variable)
# cdf = ecdf(temp[1])
# print(cdf)

ax2 = ax1.twinx()
ax2.plot(twl_data, twl_pdf, 'r', label='empirical pdf')
ax2.set_ylabel('probability')
ax2.set_ylim([0, 0.7])
#plt.show()


# fit a generalized pareto and get params (MLE)
shape, location, scale = genpareto.fit(variable, -1, loc=0, scale=3)
print('shape=', round(shape,4), 'loc=', round(location, 4), 'scale=', round(scale, 4))

# get generalized pareto CDF
ejex = np.arange(min(variable), max(variable), 0.05)
twl_pareto = genpareto.pdf(ejex, shape, location, scale)
ax2.plot(ejex, twl_pareto, 'k', label='pareto pdf')

ejex = twl_data
twl_pareto = genpareto.pdf(ejex, shape, location, scale)
ax2.plot(ejex, twl_pareto, '-.k')

plt.legend()
plt.show()


kk = st.kstest(variable, genpareto.cdf, [shape, location, scale])
print(kk) # pvalue debería ser mayor a 0.05 para que fuesen de la misma distribucion??

# If the p-value is less than 0.05, we reject the null hypothesis that there's no difference between
# the means and conclude that a significant difference does exist. If the p-value is larger than 0.05,
# we cannot conclude that a significant difference exists.



