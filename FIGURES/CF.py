import numpy
import scipy.constants as codata


def get_sigmas(name):

    if name == "ESRF-HB":
        return {'sigmaX':387.8e-6,"sigmaZ":3.5e-6,"sigmaX'":10.3e-6,"sigmaZ'":1.2e-6}
        # return {'sigmaX':414.97e-6,"sigmaZ":3.43e-6,"sigmaX'":10.31e-6,"sigmaZ'":1.16e-6}
    elif name == "EBS":
        # return {'sigmaX':27.2e-6, "sigmaZ":3.4e-6,"sigmaX'":5.2e-6,"sigmaZ'":1.4e-6}
        return {'sigmaX':30.18e-6, "sigmaZ":3.64e-6,"sigmaX'":4.37e-6,"sigmaZ'":1.37e-6}
    elif name == "EBS_LIN":
        ex = 147e-12
        ez = 5e-12
        bx = 5.2
        bz = 2.4
        return {'sigmaX':numpy.sqrt(ex*bx), "sigmaZ":numpy.sqrt(ez*bz),"sigmaX'":numpy.sqrt(ex/bx),"sigmaZ'":numpy.sqrt(ez/bz)}
    elif name == "ELETTRA2":
        return {'sigmaX':8.071e-05,"sigmaZ":4.66e-06,"sigmaX'":8.968e-06,"sigmaZ'":1.553e-06 }
    else:
        raise Exception("Unknown storage ring name")

def get_sigmas_radiation(photon_energy,undulator_length,method='elleaume'):
    lambdan = 1e-10 * codata.h*codata.c/codata.e*1e10 / photon_energy # in m
    if method == 'elleaume':
        return 2.740/4/numpy.pi*numpy.sqrt(lambdan*undulator_length),0.69*numpy.sqrt(lambdan/undulator_length),
    elif method == 'elleaume_approx':
        return numpy.sqrt(lambdan*undulator_length/2/numpy.pi**2),numpy.sqrt(lambdan/2/undulator_length),
    elif method == 'kim':
        return 1/4/numpy.pi*numpy.sqrt(2*lambdan*undulator_length),numpy.sqrt(lambdan/2/undulator_length),
    elif method == 'lin':
        return 1/4/numpy.pi*numpy.sqrt(2*lambdan*undulator_length),numpy.sqrt(lambdan/2/undulator_length),

if __name__ == "__main__":

    from srxraylib.plot.gol import plot, plot_table 
    import matplotlib.pylab as plt

    undulator_length = 2.0
    photon_energies = numpy.logspace(numpy.log10(10.),numpy.log10(100000.),100)
    #photon_energies = [5000,30000]

   
    m = get_sigmas("EBS")


    methods = ["elleaume","elleaume_approx","kim","lin"]
    methods_labels = ["$CF^E$","$CF^{EA}$","$CF^K$","lin"]

    out = numpy.zeros((photon_energies.size,1+len(methods)))

    print("Electron beam sigmas: ",m)

    for imethod,method in enumerate(methods):
        for iphoton_energy,photon_energy in enumerate(photon_energies):
            lambdan = 1e-10 * codata.h*codata.c/codata.e*1e10 / photon_energy # in m
            #print("wavelength [m]",lambdan)
            sr,srp = get_sigmas_radiation(photon_energy,undulator_length,method=method)
    
            #print("photon size, angle: ",sr,srp)
            f2dot35 = 2*numpy.sqrt(2*numpy.log(2))
            sx,sz,sxp,szp = m['sigmaX'],m['sigmaZ'],m["sigmaX'"],m["sigmaZ'"]
    
            Sx2 = sx**2 + sr**2
            Sz2 = sz**2 + sr**2
            Sxp2 = sxp**2 + srp**2
            Szp2 = szp**2 + srp**2
    
            Sx = numpy.sqrt( Sx2 )
            Sz = numpy.sqrt( Sz2 )
            Sxp = numpy.sqrt( Sxp2 )
            Szp = numpy.sqrt( Szp2 )
    
    
            if method == "kim":
                coeff = 1.0 / 4 / numpy.pi
            else:
                coeff = 1.0 / 2 / numpy.pi
            
            CF = (lambdan*coeff)**2 / numpy.sqrt( Sx2 * Sz2 * Sxp2 * Szp2 )
    
            Wx = f2dot35 * Sx
            Wxp = f2dot35 * Sxp
            Wz = f2dot35 * Sz
            Wzp = f2dot35 * Szp
    
            print(">>> Energy: %8.1f eV, FWHM H: %5.2f urad, FWHM V: %5.2f urad, CF: %g"%(photon_energy,Wxp,Wzp,CF))
    
            out[iphoton_energy,0] = photon_energy
            out[iphoton_energy,1+imethod] = CF


    # plot(out[:,0],out[:,1],xlog=1,ylog=1)
    p = plot_table(out[:,0].T,out[:,1:4].T,xlog=1,ylog=1,legend=methods_labels,
        xtitle="photon enetgy [eV]",ytitle="Coherent Fraction",show=0)
    plt.savefig("CF.eps")
    p.show()
