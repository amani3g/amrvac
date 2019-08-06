import numpy as np
import matplotlib.pyplot as plt
import get_data


#file name and directory                                                                                          
dir = '/Users/amani3g/jupiter/amrvac/jov/test_output/'
outdir = '/Users/amani3g/jupiter/amrvac/jov/output/plots/'


#set up plotting devices
fig = plt.figure(figsize=(8.0, 8.0))
ax = fig.add_subplot(111, projection='3d')

#plot jupiter wireframe
theta, phi = np.mgrid[0:2*np.pi:20j,0:np.pi:10j]

rj = 7149200000#[m]
xj = rj*np.cos(theta)*np.sin(phi)
yj = rj*np.sin(theta)*np.sin(phi)
zj = rj*np.cos(phi)
ax.plot_wireframe(xj,yj,zj,color='C1') #the orange spot!


#################################
# Plot multiple particle tracks 
##################################

#get data from file
test_name = 'iprob_4_x0_5rj_v0_100eV_qe_qm_'
particle_name = 'particle_0000'

#loop through each particle track
num_particles = 10 #set the number of particles
for i in range(1, num_particles+1):

	if i < 10:
		file_name = test_name + particle_name + str(0) + str(i) + '.csv'
		xline1, yline1, zline1, tline1 = get_data.get_lists(dir, file_name,rj)
		#plot data
		plt.plot(xline1, yline1, zline1)
	elif i >= 10: 
		file_name = test_name + particle_name + str(i) + '.csv'
		print(file_name)
		xline1, yline1, zline1, tline1 = get_data.get_lists(dir, file_name,rj)
		#plot data
		plt.plot(xline1, yline1, zline1)

	
	
	


########################################
# Plot bfield along particles trajectory
#########################################

'''
#get data from magnetic field file
#bcsv = "bfield_dipole_2.csv"
#header = [' x(1)', 'x(2)', 'x(3)', 'B(1)', 'B(2)', 'B(3)']

dir = "/Users/amani3g/jupiter/amrvac/jov/jrm09_output/"
bcsv = "bfield_testing_inter_dipole_4_jul26.csv"

xb, yb, zb, b1, b2, b3, bmag = get_data.get_blists(dir, bcsv)

#shorten all arrays
group = 1000 #for array of size 192000, 191976

xb = get_data.shorten_arr(xb,group)
yb = get_data.shorten_arr(yb,group)
zb = get_data.shorten_arr(zb,group)

b1 = get_data.shorten_arr(b1,group)
b2 = get_data.shorten_arr(b2,group)
b3 = get_data.shorten_arr(b3,group)


#plot quiver plot of bfield
#length of arrows should be proportional to 10nT
l_arrow = 10**11
ax.quiver(xb, yb, zb, b1, b2, b3, length = l_arrow, normalize=True)
plt.title('Dipole Field')
'''
plt.title('JRM09 Test Case 1')
plt.show()
#plt.savefig(outdir + test_name + str(num_particles) + '.png')


