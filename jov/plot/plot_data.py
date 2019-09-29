import numpy as np
import matplotlib.pyplot as plt
import get_data
import plot_routines as plr

#constants
rj = 6991100000 #cm
num = 100

#file name and directory                                                                                          
dir = '/Users/amani3g/jupiter/amrvac/jov/output_jrm09_2/5rj/100keV/'
outdir = '/Users/amani3g/jupiter/amrvac/jov/output_jrm09_2/5rj/100keV/'
file_base = 'run_onlyz_singlenormal_'

#set up plotting devices
fig = plt.figure(figsize=(8.0, 8.0))
ax = fig.add_subplot(111, projection='3d')

plr.jupiter_wireframe(ax)
file_names = plr.get_file_names(dir, file_base, num)
#print(file_names)
#index = [49,50,51]
#new_file_names = np.delete(file_names, index)
#print(new_file_names[49])
#plr.plot_3Dtrajectory_raw(file_names, dir, ax)
plr.particle_3Dtrajectories(file_names, dir, rj, ax)
#ax.margins(rj,rj, 1)




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
plt.title('JRM09 Electrons at 5Rj 100keV')
plt.show()
#plt.savefig(outdir + test_name + str(num_particles) + '.png')


