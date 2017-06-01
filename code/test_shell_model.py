import numpy as np
import shell_model
import os
class TestPPP2PPV:
	def setUp(self):
		pass
	def tearDown(self):
		pass
	def test_random_preserve_total_flux(self):
		n = 10
		ppp = np.random.random((n, n, n))
		vel = np.random.random((n, n, n))
		vcen = np.linspace(0.05, 0.95, 10)
		ppv = shell_model.ppp2ppv(ppp, vel, vcen)

		assert abs((np.sum(ppv) - np.sum(ppp)) / np.sum(ppp)) < 1e-3 


class TestModel:
	def setUp(self):
		pass
	def tearDown(self):
		os.remove("new.fits")
	def test_default_par_write(self):
		shell_model.ppv_model(outfile="new.fits")