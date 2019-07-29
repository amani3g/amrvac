import csv
import numpy as np
import matplotlib.pyplot as plt
import yt
from mpl_toolkits import mplot3d


##########################################
# Script to plot AMRVAC output and Bfield
# Amani Garvin
# NuSTAR 2019
##########################################
   
#define the payloads that can be chosen in AMRVAC 
#payload_dict = {}


def get_lists(dir, file_name, payloads): #get_lists(dir, file_name, payloads)
    ##### really want to make full use of dict object in this function but it works so that's a very low priority

    #####################################
    #get payloads from payload dict
    #####################################


    csv_file = dir+file_name
    # Returns a dictionary containing the csv file
    with open(csv_file, newline = '') as csvfile:
        dict = csv.DictReader(csvfile, delimiter=',')

        ##############################################################   
        #initialize a n_payloadxm_rows array data_arr[n_payload,m_row]
        #for row in dict:
            #for payload in payloads:
                #data_arr[n_payload, ].append(row.get(key_payload))

        #data_arr = np.asarray(data_arr, dtype=np.float64)
        ##############################################################

        #for row in dict: print(row[' x1'])
        xline = []
        yline = []
        zline = []
        tline = []

        for row in dict:
            xline.append(row.get(' x1'))
            yline.append(row.get(' x2'))
            zline.append(row.get(' x3'))
            tline.append(row.get(' time'))
        #print(len(xline), len(yline), len(zline), len(tline))

        xline = np.asarray(xline, dtype=np.float64)
        yline = np.asarray(yline, dtype=np.float64)
        zline = np.asarray(zline, dtype=np.float64)
        tline = np.asarray(tline, dtype=np.float64)
        #print(xline,yline,zline,tline)

    return xline, yline, zline, tline

def get_lists_tmax(dir, file_name, tmax):
    
    csv_file = dir+file_name
    # Returns a dictionary containing the csv file
    with open(csv_file, newline = '') as csvfile:
        dict = csv.DictReader(csvfile, delimiter=',')
    
        xline = []
        yline = []
        zline = []
        tline = []

        for row in dict:
            trow = float(row.get(' time'))
            if trow < tmax:
                xline.append(row.get(' x1'))
                yline.append(row.get(' x2'))
                zline.append(row.get(' x3'))
                tline.append(row.get(' time'))

        xline = np.asarray(xline, dtype = np.float64)
        yline = np.asarray(yline, dtype = np.float64)
        zline = np.asarray(zline, dtype = np.float64)
        tline = np.asarray(tline, dtype = np.float64)

                           
    return xline, yline, zline, tline

def get_lists_rmax(dir, file_name, rlim):
    csv_file = dir+file_name
    # Returns a dictionary containing the csv file
    with open(csv_file, newline = '') as csvfile:
        dict = csv.DictReader(csvfile, delimiter=',')
    
        xline = []
        yline = []
        zline = []
        tline = []

        for row in dict:
            x = float(row.get(' x1'))
            y = float(row.get(' x2'))
            z = float(row.get(' x3'))
            
            #calculate r
            r = np.sqrt(x**2 + y**2 + z**2)

            #

            if(r > rlim):
                #print(r)
                xline.append(row.get(' x1'))
                yline.append(row.get(' x2'))
                zline.append(row.get(' x3'))
                tline.append(row.get(' time'))

        xline = np.asarray(xline, dtype = np.float64)
        yline = np.asarray(yline, dtype = np.float64)
        zline = np.asarray(zline, dtype = np.float64)
        tline = np.asarray(tline, dtype = np.float64)

    return xline, yline, zline, tline

def get_minmax(dir, filename):
    
    #get the min max of the particle's trajectory
    #useful for plotting the b-field in the area of the particle

    xline, yline, zline, tline = get_lists(dir, filename)

    xvals = [np.amin(xline), np.amax(xline)]
    yvals = [np.amin(yline), np.amax(yline)]
    zvals = [np.amin(zline), np.amax(zline)]

    

    return xvals, yvals, zvals

def get_blists(dir, file_name, header=None):

    file = dir+file_name
    if header==None:
        #If no header is provided, assume that the data is arranged
        #by x1, x2, x3, b1, b2, b3
        with open(file, newline = '') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')

            x1 = []
            x2 = []
            x3 = []
            b1 = []
            b2 = []
            b3 = []

            for row in reader:
                x1.append(row[0])
                x2.append(row[1])
                x3.append(row[2])
                b1.append(row[3])
                b2.append(row[4])
                b3.append(row[5])

            x1 = np.asarray(x1, dtype = np.float64)
            x2 = np.asarray(x2, dtype = np.float64)
            x3 = np.asarray(x3, dtype = np.float64)
            b1 = np.asarray(b1, dtype = np.float64)
            b2 = np.asarray(b2, dtype = np.float64)
            b3 = np.asarray(b3, dtype = np.float64)
            bmag = np.sqrt(b1**2 + b2**2 + b3**2)
    elif header !=None:
        #If a header is provided, it will be an array of keys
        #which wil be used for dictreader 
        with open(file, newline = '') as csvfile:
            dict = csv.DictReader(csvfile, delimiter=',')

            x1 = []
            x2 = []
            x3 = []
            b1 = []
            b2 = []
            b3 = []

            for row in dict:
                x1.append(row.get(header[0]))
                x2.append(row.get(header[1]))
                x3.append(row.get(header[2]))
                b1.append(row.get(header[3]))
                b2.append(row.get(header[4]))
                b3.append(row.get(header[5]))

            x1 = np.asarray(x1, dtype = np.float64)
            x2 = np.asarray(x2, dtype = np.float64)
            x3 = np.asarray(x3, dtype = np.float64)
            b1 = np.asarray(b1, dtype = np.float64)
            b2 = np.asarray(b2, dtype = np.float64)
            b3 = np.asarray(b3, dtype = np.float64)
            bmag = np.sqrt(b1**2 + b2**2 + b3**2)

    return x1, x2, x3, b1, b2, b3, bmag

#####################################
# These get BList functions specific to AMRVAC
# Need to change how Blists are output in AMRVAC
######################################

def get_amrvac_blists(dir, file_name):
    ##### really want to make full use of dict object in this function but it works so that's a very low priority

    csv_file = dir+file_name
    # Returns a dictionary containing the csv file
    with open(csv_file, newline = '') as csvfile:
        dict = csv.DictReader(csvfile, delimiter=',')
    
        #for row in dict: print(row[' x1'])
        xline = []
        yline = []
        zline = []
        b1 = []
        b2 = []
        b3 = []

        for row in dict:
            #these headers are specific to fields saved through AMRVAC
            xline.append(row.get(' x(1)'))
            yline.append(row.get('x(2)'))
            zline.append(row.get('x(3)'))
            b1.append(row.get('B(1)'))
            b2.append(row.get('B(2)'))
            b3.append(row.get('B(3)'))

        #print(len(xline), len(yline), len(zline), len(tline))

        xline = np.asarray(xline, dtype=np.float64)
        yline = np.asarray(yline, dtype=np.float64)
        zline = np.asarray(zline, dtype=np.float64)
        b1 = np.asarray(b1, dtype=np.float64)
        b2 = np.asarray(b2, dtype=np.float64)
        b3 = np.asarray(b3, dtype=np.float64)
        #print(xline,yline,zline,tline)
        bmag = np.sqrt(b1**2 + b2**2 + b3**2)

    return xline, yline, zline, b1, b2, b3, bmag

def get_blists_rlim(dir, file_name, rlim):

    csv_file = dir+file_name
    # Returns a dictionary containing the csv file
    with open(csv_file, newline = '') as csvfile:
        dict = csv.DictReader(csvfile, delimiter=',')
    
        #for row in dict: print(row[' x1'])
        xline = []
        yline = []
        zline = []
        b1 = []
        b2 = []
        b3 = []

        for row in dict:
            #headers specific to fields saved through AMRVAC
            x1 = row.get(' x(1)')
            x2 = row.get('x(2)')
            x3 = row.get('x(3)')
            r = np.sqrt(x1**2 + x2**2 + x3**2)
            if r<rlim or r==rlim:
                xline.append(row.get(' x(1)'))
                yline.append(row.get('x(2)'))
                zline.append(row.get('x(3)'))
                b1.append(row.get('B(1)'))
                b2.append(row.get('B(2)'))
                b3.append(row.get('B(3)'))

        #print(len(xline), len(yline), len(zline), len(tline))

        xline = np.asarray(xline, dtype=np.float64)
        yline = np.asarray(yline, dtype=np.float64)
        zline = np.asarray(zline, dtype=np.float64)
        b1 = np.asarray(b1, dtype=np.float64)
        b2 = np.asarray(b2, dtype=np.float64)
        b3 = np.asarray(b3, dtype=np.float64)
        #print(xline,yline,zline,tline)
        bmag = np.sqrt(b1**2 + b2**2 + b3**2)

    return xline, yline, zline, b1, b2, b3, bmag

def shorten_arr(arr, n_prime):
    n = int(arr.size) 
    group = int(n/n_prime)
    if(n%group == 0): #group has to be a multiple of n      
        arr_prime = arr.reshape(-1,group).mean(axis=1) #sliding window average over array
        arr_prime[0] = arr[0] #keep edges
        arr_prime[n_prime-1] = arr[n-1] #keep edges
    else:
        print("Enter new array size that is a multiple of the original array size: " + str(n))
        return arr 
    return arr_prime

