from spectral_cube import SpectralCube
#Calculate various physical quantities from
#spectral cubes, spectra, and shell parameters.

#Equation 1 in Arce+ 2011, from Rohlfs and Wilson 1996
def Tex(Tpeak, thick=True):
	"""
	Find the excitation temperature given peak temperature
	of an optically thick line. 
	"""
	if thick:
    	Tex = 5.53 / np.log(1 + (5.53)/(Tpeak+0.82))
    	return Tex
    else:
    	raise("Line must be optically thick.")


def extract_shell(cube=None, ra_center=None, dec_center=None,
    radius=None, thickness=None, velocity_range=None):
	    


    if velocity_range:


    return shell_pixels


def momentum(mass, velocity):
    return mass * velocity
def energy(mass, velocity):
    return 0.5 * mass * velocity ** 2.