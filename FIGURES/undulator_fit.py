import numpy
import matplotlib.pylab as plt

from silx.math.fit.functions import sum_gauss, sum_agauss
from silx.math.fit import fittheories
from silx.math.fit.fitmanager import FitManager


x = numpy.linspace(-4,4,100)
Gamma= numpy.pi / 2 * x**2
y = (numpy.sin(Gamma)/Gamma)**2


# Fitting
p = [y.max(),0,2]
fit = FitManager()
fit.setdata(x=x, y=y)
fit.loadtheories(fittheories)
fit.settheory('Gaussians')
# fit.settheory('Area Gaussians')
# fit.settheory('Slit')
fit.estimate()
fit.runfit()

print("Searched parameters = %s" % p)
print("Obtained parameters : ")
dummy_list = []
for param in fit.fit_results:
    # print(param)
    print(param['name'], ' = ', param['fitresult'])
    if param['name'] == "FWHM1":
        print("sigma1 = ",param['fitresult']/(2*numpy.sqrt(2*numpy.log(2))))
    dummy_list.append(param['fitresult'])
print("chisq = ", fit.chisq)

y1 = sum_gauss(x, *dummy_list)
# y1 = sum_agauss(x, *dummy_list)

print("Area original: ",y.sum()*(x[1]-x[0]))
print("Area fit: ",y1.sum()*(x[1]-x[0]))

plt.rc('text', usetex=True)
plt.plot(x,y)
plt.xlabel(r'$\theta_r=\theta\sqrt{L/\lambda}$') # =\sigma_{\mu}/\sigma_I')
plt.ylabel(r'$(sin(\Gamma)/\Gamma)^2$')
plt.plot(x,y1)
plt.savefig('undulator_fit.eps')
plt.show()


