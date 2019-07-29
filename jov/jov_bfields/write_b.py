import numpy as np
import scipy
import solve_b
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyshtools
import csv
import os

#############################################
# Plot and Save JRM09 B-Fields Solved for mult*rj
# Amani Garvin
# NuSTAR Summer 2019
#############################################
rj = 71492000 #m


def write_jrm09():

	############################
	# Write JRM09 to CSV grid in Spherical Coords
	#
	#
	#
	############################


	#Create a folder titled JRM09_GRID
	jrm09 = '/Users/amani3g/jupiter/amrvac/jov/JRM09_grid'
	os.mkdir(jrm09)

	#Solve for Vectors
	r_vec = np.linspace(rj, 30*rj, 22)
	theta_vec = solve_jrm09.jrm09(rj)[6] 
	phi_vec = solve_jrm09.jrm09(rj)[7]

	print('Writing csvs!')
	# 1: Create JRM09_Vectors.csv
	jrm09_vectors = jrm09 + '/JRM09_Vectors.csv'
	with open(jrm09_vectors, mode='w') as vector_file:
		vector_writer = csv.writer(vector_file, delimiter = ',')

		# Write coordinate vectors to JRM09_Vectors.csv as columns, eg:
		# r_vec, phi_vec, theta_vec
		# r1, phi1, theta1
		# r2, phi2, theta2
		# .
		# .
		# .
		# rN, phiN, thetaN
		print('\t Writing vectors to csv...')
		for i in range(r_vec.size):
			vector_writer.writerow([r_vec[i], theta_vec[i], phi_vec[i]])
		print('\t Done!')
	# 2: For each r in r_vec
	for r in r_vec:
		mult = int(r/rj)
		r_surface_dir = jrm09 + '/' + str(mult)+'_surface'
		# Create a folder r_surface
		os.mkdir(r_surface_dir)
		# Solve for Br, Btheta, Bphi, Bmag
		mag, Br, Btheta, Bphi, Bmag, Vp, lat, lon = solve_jrm09.jrm09(r)

		# Write Br, Btheta, Bphi, Bmag to csv

		# First write Br r_surface/Br.csv
		Br_values = r_surface_dir + '/Br.csv'
		with open(Br_values, mode = 'w') as Br_file:
			Br_writer = csv.writer(Br_file, delimiter=',')
			
			# The values iterate as Br[theta, phi], eg:
			# Br[theta1, phi1], Br[theta1, phi2], Br[theta1, phi3],...., Br[theta1, phiN]
			# Br[theta2, phi1], Br[theta2, phi2], Br[theta2, phi3],...., Br[theta2, phiN]
			# .
			# . 	 	 	 	 	 	Br[theta(i), phi(j)]
			# .														.
			# Br[thetaN, phi1], Br[thetaN, phi2], Br[thetaN, phi3],...., Br[thetaN, phiN]
			print('\t Writing radial csv...')
			for i in range(theta_vec.size):
				Br_writer.writerow(Br[i,])
			print('\t Done!')

		# Then write Btheta r_surface/Btheta.csv
		Btheta_vales = r_surface_dir + '/Btheta.csv'
		with open(Btheta_vales, mode = 'w') as Btheta_file:
			Btheta_writer = csv.writer(Btheta_file, delimiter = ',')

			# The values iterate as Btheta[theta, phi], eg:
			# Btheta[theta1, phi1], Btheta[theta1, phi2], Btheta[theta1, phi3],...., Btheta[theta1, phiN]
			# Btheta[theta2, phi1], Btheta[theta2, phi2], Btheta[theta2, phi3],...., Btheta[theta2, phiN]
			# .
			# . 	 	 	 	 	 	.
			# .														.
			# Btheta[thetaN, phi1], Btheta[thetaN, phi2], Btheta[thetaN, phi3],...., Btheta[thetaN, phiN]
			print('\t Writing theta csv...')
			for i in range(theta_vec.size):
				for j in range(phi_vec.size):
					Btheta_writer.writerow(Btheta[i,])
			print('\t Done!')

		# Last write Bphi r_surface/Bphi.csv
		Bphi_values = r_surface_dir + '/Bphi.csv'
		with open(Bphi_values, mode = 'w') as Bphi_file:
			Bphi_writer = csv.writer(Bphi_file, delimiter = ',')
			# The values iterate as Bphi[theta, phi], eg:
			# Bphi[theta1, phi1], Bphi[theta1, phi2], Bphi[theta1, phi3],...., Bphi[theta1, phiN]
			# Bphi[theta2, phi1], Bphi[theta2, phi2], Bphi[theta2, phi3],...., Bphi[theta2, phiN]
			# .
			# . 	 	 	 	 	 	.
			# .														.
			# Bphi[thetaN, phi1], Bphi[thetaN, phi2], Bphi[thetaN, phi3],...., Bphi[thetaN, phiN]
			print('\t Writing phi csv...')
			for i in range(theta_vec.size):
				Bphi_writer.writerow(Bphi[i,])
			print('\t Done!')

		print('All csvs are written!')

	return

def write_dipole():

	############################
	# Write Dipole to CSV grid in Spherical Coords
	#
	#
	#
	############################


	#Create a folder titled DIPOLE_GRID
	dipole = '/Users/amani3g/jupiter/amrvac/jov/dipole_grid'
	os.mkdir(dipole)

	#Solve for Vectors
	r_vec = np.linspace(rj, 30*rj, 22)
	theta_vec = solve_jrm09.dipole(rj)[4] 
	phi_vec = solve_jrm09.dipole(rj)[5]

	print('Writing csvs!')
	# 1: Create dipole_vectors.csv
	dipole_vectors = dipole + '/dipole_vectors.csv'
	with open(dipole_vectors, mode='w') as vector_file:
		vector_writer = csv.writer(vector_file, delimiter = ',')
		# Write coordinate vectors to dipole_Vectors.csv as columns, eg:
		# r_vec, phi_vec, theta_vec
		# r1, phi1, theta1
		# r2, phi2, theta2
		# .
		# .
		# .
		# rN, phiN, thetaN
		print('\t Writing vectors to csv...')
		for i in range(r_vec.size):
			vector_writer.writerow([r_vec[i], theta_vec[i], phi_vec[i]])
		print('\t Done!')

	# 2: For each r in r_vec
	for r in r_vec:
		mult = int(r/rj)
		r_surface_dir = dipole + '/' + str(mult)+'_surface'
		# Create a folder r_surface
		os.mkdir(r_surface_dir)
		# Solve for Br, Btheta, Bphi, Bmag
		Br, Btheta, Bphi, Bmag, lat, lon = solve_jrm09.dipole(r)

		# Write Br, Btheta, Bphi, Bmag to csv

		# First write Br r_surface/Br.csv
		Br_values = r_surface_dir + '/Br.csv'
		with open(Br_values, mode = 'w') as Br_file:
			Br_writer = csv.writer(Br_file, delimiter=',')
			
			# The values iterate as Br[theta, phi], eg:
			# Br[theta1, phi1], Br[theta1, phi2], Br[theta1, phi3],...., Br[theta1, phiN]
			# Br[theta2, phi1], Br[theta2, phi2], Br[theta2, phi3],...., Br[theta2, phiN]
			# .
			# . 	 	 	 	 	 	Br[theta(i), phi(j)]
			# .														.
			# Br[thetaN, phi1], Br[thetaN, phi2], Br[thetaN, phi3],...., Br[thetaN, phiN]
			print('\t Writing radial csv...')
			for i in range(theta_vec.size):
				Br_writer.writerow(Br[i,])
			print('\t Done!')

		# Then write Btheta r_surface/Btheta.csv
		Btheta_vales = r_surface_dir + '/Btheta.csv'
		with open(Btheta_vales, mode = 'w') as Btheta_file:
			Btheta_writer = csv.writer(Btheta_file, delimiter = ',')

			# The values iterate as Btheta[theta, phi], eg:
			# Btheta[theta1, phi1], Btheta[theta1, phi2], Btheta[theta1, phi3],...., Btheta[theta1, phiN]
			# Btheta[theta2, phi1], Btheta[theta2, phi2], Btheta[theta2, phi3],...., Btheta[theta2, phiN]
			# .
			# . 	 	 	 	 	 	.
			# .														.
			# Btheta[thetaN, phi1], Btheta[thetaN, phi2], Btheta[thetaN, phi3],...., Btheta[thetaN, phiN]
			print('\t Writing theta csv...')
			for i in range(theta_vec.size):
				for j in range(phi_vec.size):
					Btheta_writer.writerow(Btheta[i,])
			print('\t Done!')

		# Last write Bphi r_surface/Bphi.csv
		Bphi_values = r_surface_dir + '/Bphi.csv'
		with open(Bphi_values, mode = 'w') as Bphi_file:
			Bphi_writer = csv.writer(Bphi_file, delimiter = ',')
			# The values iterate as Bphi[theta, phi], eg:
			# Bphi[theta1, phi1], Bphi[theta1, phi2], Bphi[theta1, phi3],...., Bphi[theta1, phiN]
			# Bphi[theta2, phi1], Bphi[theta2, phi2], Bphi[theta2, phi3],...., Bphi[theta2, phiN]
			# .
			# . 	 	 	 	 	 	.
			# .														.
			# Bphi[thetaN, phi1], Bphi[thetaN, phi2], Bphi[thetaN, phi3],...., Bphi[thetaN, phiN]
			print('\t Writing phi csv...')
			for i in range(theta_vec.size):
				Bphi_writer.writerow(Bphi[i,])
			print('\t Done!')

		print('All csvs are written!')

	return

write_dipole()
