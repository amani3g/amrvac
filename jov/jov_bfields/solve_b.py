import numpy as np
import csv
import scipy
from scipy.special import legendre
import pyshtools
import matplotlib.pyplot as plt

####################
# This code solves for JRM09
# Amani Garvin
# NuSTAR
# Summer 2019
####################

def coeffs(nmax):
	#initialize a 2 x nmax+1 x nmax array
	coeffs = np.zeros((2, nmax+1, nmax+1))
	with open("jrm09_coeffs.csv") as csv_file:
		csv_reader = csv.reader(csv_file, delimiter=' ')
		for row in csv_reader:
			n, m = int(row[3]), int(row[4])
			if n<=nmax:
				if row[2] == 'g':
					coeffs[0,n,m] = row[1]
				elif row[2] == 'h':
					coeffs[1,n,m] = row[1]

	return coeffs

def jrm09(r):

	#Solves JRM09 on a sphere of radius multrj*rj

	#jupiter radius
	rj = 71492000 #m
	nmax = 10

	#Get coefficients
	ghs = coeffs(nmax)

	
	#Initialize Coeffs, Creating an SHMagCOeffs Class Instance
	clm = pyshtools.SHMagCoeffs.from_array(ghs, rj)

	mag = clm.expand(lmax = 10, a=r, sampling = 1)
	print('mag calculated!')

	#SHMagGrid saves the Br, Btheta, Bphi, Btotal, Vpot as SHGrid Class Instances

	#Get arrays of raw data
	Br = mag.rad.to_array()
	Btheta = mag.theta.to_array()
	Bphi = mag.phi.to_array()
	Bmag = mag.total.to_array()
	Vp = mag.pot.to_array()


	#Get lat and lon vectors from Br class
	lat = mag.rad.lats()
	lon = mag.rad.lons()

	return mag, Br, Btheta, Bphi, Bmag, Vp, lat, lon

def sphere_cartesian(r, lat, lon):

	x = r * np.cos(lon) * np.sin(lat)
	y = r * np.sin(lon) * np.sin(lat)
	z = r * np.cos(lat)

	return x, y, z

def dipole(r):

	#Solves for a dipole field in spherical coordinates
	rj = 71492000 #m
	M = 4.170 #G -> need to convert to nT

	#set up lat and lon in a similar manner to pshy
	lat_ = np.linspace(-np.pi/2, np.pi/2, 22)
	lon_ = np.linspace(-np.pi, np.pi, 22)
	lat, lon = np.meshgrid(lat_,lon_)

	Br = (2*M*np.cos(lat))/(r**3)
	Btheta = (M*np.sin(lat))/(r**3)
	Bphi = np.zeros((22,22))
	Bmag = (Br**2 + Btheta**2 + Bphi**2)**.5

	
	return Br, Btheta, Bphi, Bmag, lat_, lon_

mag, Br, Btheta, Bphi, Bmag, Vp, lat, lon = jrm09(71492000)
mag.plot()
plt.show()
