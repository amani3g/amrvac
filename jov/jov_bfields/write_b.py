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

def write_dipole_mult(n_mult):

	############################
	# Write Dipole to CSV grid in Spherical Coords
	# n_mult is the last multiple of rj the grid is calculated to
	#
	#
	############################


	#Create a folder titled DIPOLE_GRID
	dipole = '/Users/amani3g/jupiter/amrvac/jov/dipole_grid'
	os.mkdir(dipole)

	#Solve for Vectors
	mult_vec = np.linspace(1, n_mult, 22)
	r_vec = np.linspace(rj, n_mult*rj, 22)
	theta_vec = solve_b.dipole(rj)[4] 
	phi_vec = solve_b.dipole(rj)[5]

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
			vector_writer.writerow([mult_vec[i], r_vec[i], theta_vec[i], phi_vec[i]])
		print('\t Done!')

	# 2: For each r in r_vec
	for i in range(0, mult_vec.size, 1):
		mult = mult_vec[i]
		r_surface_dir = dipole + '/' + str(mult) +'_rj_surface'
		# Create a folder r_surface
		os.mkdir(r_surface_dir)
		# Solve for Br, Btheta, Bphi, Bmag
		Br, Btheta, Bphi, Bmag, lat, lon = solve_b.dipole(r_vec[i])

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

def write_jrm09_mult(n_mult):

	############################
	# Write JRM09 to CSV grid in Spherical Coords
	# 
	#
	############################


	#Create a folder titled JRM09_GRID
	jrm09 = '/Users/amani3g/jupiter/amrvac/jov/jrm09_grid'
	os.mkdir(jrm09)

	#Solve for Vectors
	mult_vec = np.linspace(1, n_mult, 22)
	r_vec = np.linspace(rj, n_mult*rj, 22)
	theta_vec = solve_b.jrm09(rj)[6] 
	phi_vec = solve_b.jrm09(rj)[7]

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
			vector_writer.writerow([mult_vec[i], r_vec[i], theta_vec[i], phi_vec[i]])
		print('\t Done!')
	# 2: For each r in r_vec
	for i in range(0, mult_vec.size, 1):
		mult = int(mult_vec[i])
		if mult < 10: 
			r_surface_dir = jrm09 + '/0' + str(mult) +'_rj_surface'
		else:
			r_surface_dir = jrm09 + '/' + str(mult) +'_rj_surface'
		# Create a folder r_surface
		os.mkdir(r_surface_dir)
		# Solve for Br, Btheta, Bphi, Bmag
		mag, Br, Btheta, Bphi, Bmag, Vp, lat, lon = solve_b.jrm09(r_vec[i])

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

def write_jrm09_morePrecise(N1rj, N2rj, numrj1, N3rj):

	############################
	# Write JRM09 to CSV grid in Spherical Coords
	# Spacing is delrj1 between N1rj and N2rj
	# Spacing is delrj2 between N2rj and N3rj
	# The number of elements in the r vector is 22 always 
	# This is due to the SHTOOLS module which only gives 
	# 22 elements to each angular vector
	############################

	#Folder name to store JRM09 grid
	jrm09 = '/Users/amani3g/jupiter/amrvac/jov/jrm09_grid_2'
	#If folder does not exist
	if os.path.exists(jrm09) == False:
		#Create a folder titled JRM09_GRID
		os.mkdir(jrm09)
	else: 
		print(jrm09 + ' Exists, overwrite files?')
		overwrite = input('Overwrite (yes or no)? ')
		if overwrite=='no':
			jrm09 = input('New folder name: ')
			os.mkdir(jrm09)
			print(jrm09 + 'Created')
		if overwrite=='yes':
			print('Overwriting...')

	#Solve for Vectors
	print('Solving for Vectors')

	#Solve for Index Vector
	numrj2 = 22 - numrj1
	index_vec1 = np.linspace(N1rj, N2rj, numrj1)
	index_vec2 = np.linspace(N2rj+1, N3rj, numrj2)
	index_vec = np.hstack((index_vec1, index_vec2))

	#Format String Index Vector
	#Change decimal vals as strings where 000.00 -> 00000
	index_vec_str = []
	for i in index_vec:
		i_= str(int(i*100))
		#1.25 -> 125
		#22 -> 2200
		if i < 10:
			#125->0125
			i_str = '0'+i_
		else:
			i_str = i_
		index_vec_str.append(i_str)
	index_vec_str = np.asarray(index_vec_str)

	#Solve for R-Vec
	r_vec = rj * index_vec
	#Solve for Theta-Vec
	theta_vec = solve_b.jrm09(rj)[6] 
	#Solve for Phi-Vec
	phi_vec = solve_b.jrm09(rj)[7]

	#Begin Writing CSVs
	print('Writing csvs!')

	# 1: Create JRM09_Vectors.csv
	jrm09_vectors = jrm09 + '/JRM09_Vectors.csv'
	with open(jrm09_vectors, mode='w') as vector_file:
		vector_writer = csv.writer(vector_file, delimiter = ',')

		# Write coordinate vectors to JRM09_Vectors.csv as columns, eg:
		# index_vec_str, r_vec, phi_vec, theta_vec
		# istr1, r1, phi1, theta1
		# istr2, r2, phi2, theta2
		# .
		# .
		# .
		# istrN, rN, phiN, thetaN
		print('\t Writing vectors to csv...')
		for i in range(r_vec.size):
			vector_writer.writerow([index_vec_str[i], r_vec[i], theta_vec[i], phi_vec[i]])
		print('\t Done!')
	# 2: For each r in r_vec
	for i in range(0, index_vec.size, 1):

		#Create nindex_rj_surface/ directory
		r_surface_dir = jrm09 + '/' + index_vec_str[i] + '_rj_surface'
		os.mkdir(r_surface_dir)

		# Solve for Br, Btheta, Bphi, Bmag at that rj
		mag, Br, Btheta, Bphi, Bmag, Vp, lat, lon = solve_b.jrm09(r_vec[i])

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

write_jrm09_morePrecise(1, 5, 10, 50)
