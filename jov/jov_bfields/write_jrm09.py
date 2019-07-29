import numpy as np
import scipy
import solve_jrm09
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

#-----------------------------------------
# Testing Linear Interpolation in AMRVAC
# With Dipole Field Example
#-----------------------------------------
'''
mult = 1
Br, Btheta, Bphi, Bmag, lat_, lon_ = solve_jrm09.dipole(mult)

#-----------------------------------------
# Plot Dipole Field and Save
#-----------------------------------------
#set up plotting devices
fig = plt.figure(figsize=(8.0, 8.0))
ax = fig.add_subplot(111, projection='3d')

#plot jupiter wireframe
theta, phi = np.mgrid[0:2*np.pi:20j,0:np.pi:10j]

rj = 71492000#[m]
xj = rj*np.cos(theta)*np.sin(phi)
yj = rj*np.sin(theta)*np.sin(phi)
zj = rj*np.cos(phi)
ax.plot_wireframe(xj,yj,zj,color='C1') #the orange spot!

x = r*np.cos(lat_)*np.sin(lon_)
y = r*np.sin(lat_)*np.sin(lon_)
z = r*np.cos(lon_)
Bx = Br * np.cos(Btheta) * np.sin(Bphi)
By = Br * np.sin(Btheta) * np.sin(Bphi)
Bz = Br * np.cos(Bphi)

plt.title('Radial Component of Dipole Field')
l_arrow = 10**10
ax.quiver(x, y, z, Bx, By, Bz, length = l_arrow, normalize=True)
plt.show()

'''

# Understanding the Output Arrays
#---------------------------------------
# Br, Btheta, Bphi are 2D arrays
# solving for the field at r=const
# B_[lat_, lon_]
# size = 22x22

# lat and lon are vectors of length 22
# lat = [-81.81....90.0, .136 spacing] #in degrees
# lon = [0.0.....343.63, .272 spacing] #in degrees





#########################
# Writes JRM09 to CSV: Not efficient at all!
#########################
'''
mults = np.linspace(1, 10**10,22)
rs = np.linspace(rj, 10**11, 22)
with open('dipole_grid_sphere_2.csv', mode = 'w') as csvfile:
	jrm09writer = csv.writer(csvfile, delimiter = ',')
	for r in rs:
		#print(r)
		#mag, Br, Btheta, Bphi, Bmag, Vp, lat, lon = solve_jrm09.jrm09(r)
		Br, Btheta, Bphi, Bmag, lat, lon = solve_jrm09.dipole(r)
		for i in range(lat.size):
			#get lon
			theta = lat[i]
			for j in range(lon.size):
				#get br, btheta, bphi
				br = Br[i, j]
				btheta = Btheta[i, j]
				bphi = Bphi[i, j]

				#lat
				phi = lon[j]

				jrm09writer.writerow([r, theta, phi, br, btheta, bphi])
'''

############################
# Write JRM09 to CSV grids
############################


#Create a folder titled JRM09_GRID

# Create JRM09_Vectors.csv
# Write coordinate vectors to JRM09_Vectors.csv as columns, eg:
# r_vec, phi_vec, theta_vec
# r1, phi1, theta1
# r2, phi2, theta2
# rN, phiN, thetaN


# 2: For each r in r_vec
# Create a folder r_surface
# Go into the folder r_surface
# Solve for Br, Btheta, Bphi, Bmag 
# Write Br, Btheta, Bphi, Bmag to csvs, the values iterate as B[theta, phi]



