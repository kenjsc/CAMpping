import math
import numpy as np
from scipy.stats import lognorm, norm
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt


def log_normal(log_vals):
	"""Function accepts a list of values to be fitted with a log normal distribution and
	returns the mean mu and standard deviation of the normal distribution (natural log of the
	lognormal values)
	"""
	ex = np.mean(log_vals)
	var = np.var(log_vals)
	
	std = (np.log((var / ex ** 2) + 1)) ** 2
	mu = np.log(ex) - 0.5 * (std ** 2)
	#Return the mean and standard deviation of the list values
	return mu, std
	
def power_normal(vals, power):
        """Function accepts a list of values to be fitted with a power normal distribution
        and returns the mean mu and standard deviation of the normal distirbution (power root
        of the power normal values). The power value is the power to which the normal
        distribution values were rased (the inverse of which will be used)
        """
        vals = np.array(vals)
        vals = vals ** (1.0/power)
        return np.mean(vals), np.std(vals)
        
	
def plot(vals, name, bins=100, power=1, nmu=None, nstd=None):
	"""Function accepts a dictionary of orf length - frequency pairs, and the length of the
	ORF of interest, also optional plot flag for if data histogram should be plotted. Fits a
	lognormal distribution to the orf length data. Returns dictionaries of the probabilites of
	each length, the Predicted p-value for each length based on log normal, the observed
	p-value for each length, and the p-value for the length of interest.
	"""	
#	mu, std = log_normal(vals)
        loged = np.log(vals)
        if not nmu and not nstd:
            nmu, nstd = np.mean(loged), np.std(loged)
	lgn = lognorm([nstd], scale=math.e ** nmu)
	
	#Plot the fit line
	f, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
                
        ax1.scatter(range(len(vals)), sorted(vals))
        hi, pos, patch = ax2.hist(vals, bins, normed=True)
        x = np.linspace(min(vals), max(vals), 500)
        y1 = lgn.pdf(x)
#        y2 = n.pdf(x) ** power
        ax2.plot(x, y1, color='green', lw=3.0)
#        ax2.plot(x, y2, color='red', lw=3.0)
        plt.text(0.4, max(hi) * 3 / 4, "E(X) = {:4f}\nVAR(X) = {:4f}\nSTD(X) = {:4f}".format(np.mean(vals), np.var(vals), np.std(vals)), fontsize=8)
        plt.savefig(name)
        plt.close()
        
def plot_any(vals, distribution, bins=50):
        dist = distribution(*distribution.fit(vals))
        plt.hist(vals, bins=bins, normed=True)
        x = np.linspace(min(vals), max(vals), 100)
        plt.plot(x, dist.pdf(x))
        plt.show()