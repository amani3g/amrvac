import numpy as np
import matplotlib.pyplot as plt
import get_data

################################################
# Routines to plot particle outputs from AMRVAC
# Amani Garvin 
# NuSTAR Summer 2019
################################################

def jupiter_wireframe(ax):


	#plot jupiter wireframe
	theta, phi = np.mgrid[0:2*np.pi:20j,0:np.pi:10j]

	rj = 6991100000#[cm]
	xj = rj*np.cos(theta)*np.sin(phi)
	yj = rj*np.sin(theta)*np.sin(phi)
	zj = rj*np.cos(phi)
	ax.plot_wireframe(xj,yj,zj,color='C1') #the orange spot!

	return

def get_file_names(dir, file_base, num):

	file_names = []
	particle_name = 'particle_' #at most 6 trailing zeros
	#loop through each particle track
	for i in range(1,num+1):
		#format file name based on particle number
		if i < 10:
			particle_num = '00000' + str(i)
		elif i<100:
			particle_num = '0000' + str(i)
		elif i<1000:
			particle_num = '000' + str(i)
		elif i<10000:
			particle_num = '00' + str(i)
		elif i<100000:
			particle_num = '0' + str(i)
		elif i<1000000:
			particle_num = str(i)
		#get data from file
		file_name = file_base + particle_name + particle_num + '.csv'
		#print(file_name)
		file_names.append(file_name)
	return file_names

def plot_3Dtrajectory_raw(file_names, dir, ax):

	for file in file_names:
		xline1, yline1, zline1, tline1 = get_data.get_lists(dir, file)
		ax.plot(xline1, yline1, zline1)

	return

def particle_3Dtrajectories(file_names, dir, rmax, ax):


	for file in file_names:
		xline1, yline1, zline1, tline1 = get_data.get_lists_rmax(dir,file, rmax)
		ax.plot(xline1, yline1, zline1)

	return

def particle_3Dposition(file_names, dir, rmax, err, ax):



	for file in file_names:
		x, y, z, t = get_position_r(dir, file, rmax, err)
		ax.plot(x,y,z)

	return

def flat_particlepos_auroralreg(file_names, dir, rmax, err):

	#Set Mollweide Axes with Jupiter Radius



	#Get Particles final position
	for file in file_names:
		x, y, z, t = get_position_r(dir, file, rmax, err)
		#lat = 
		#lon = 
		#ax.plot(lat, lon)

	return

