#Here is the code that plots carbons in their ideal position with spins replacing them given a ratio
import numpy as np
from random import randint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def N_by_three_function(N, density, average_radius):
    #N is number of cells (of 8 carbons each)
    #density is number of spins/ carbon
    #funciton returns the rounded to nearest cubable number (8, 27, 64...)
    #ask what this values should be set at
    average_radius=1;

    #this funnction will take in N and will generate a Nx3 ideal matrix
    #representing the x, y, z location of the carbon

    #generating Nx3 array with x, y, z coord locations (starting at 0,0,0)
    coord_array=np.zeros((N**2,3))


#making the ideal matrix! (as in all spins are right on the cube corners)
#a spin on the top right corner should have x position = cube root of N
#y position = cube root of N
#z position = cube root of N
#likewise the other positions are calculated with that idea
    ii=0;
    for i in range(0,round(N**(1./3))): #x
        for j in range(0, round(N**(1./3))): #y
            for z in range(0,round(N**(1./3))): #z
                coord_array[ii, 0]=i*average_radius
                coord_array[ii, 1]=j*average_radius
                coord_array[ii, 2]=z*average_radius
                ii+=1;
                #print(ii)

    coord_array=coord_array[:ii, ::]
    return coord_array

def plot_three_d_carbons(N_by_three_array, matrix_of_spin_loc, random_spin_index):
    [m,n]=np.shape(N_by_three_array)
    ax = plt.axes(projection='3d')
    zdata = []
    xdata = []
    ydata = []

    zdata_ideal = []
    xdata_ideal = []
    ydata_ideal = []

    for i in range(0,m):

        zdata=zdata+[N_by_three_array[i,2]]
        ydata=ydata+[N_by_three_array[i,1]]
        xdata=xdata+[N_by_three_array[i,0]]

    #removing the spin carbons
    for i in range(0, len(random_spin_index)):
        zdata=zdata[: random_spin_index[i]]+zdata[random_spin_index[i]+1:]
        xdata=xdata[: random_spin_index[i]]+xdata[random_spin_index[i]+1:]
        ydata=ydata[: random_spin_index[i]]+ydata[random_spin_index[i]+1:]

    ax.scatter3D(xdata, ydata, zdata,c='r')
    #now graphing spins
    [m,n]=np.shape(matrix_of_spin_loc)
    for i in range(m):
        ax.scatter3D(matrix_of_spin_loc[i][0], matrix_of_spin_loc[i][1], matrix_of_spin_loc[i][2],c='black')
    ax.set_title("Location of Spins")
    ax.set_xlabel('X Distance')
    ax.set_ylabel('Y Distance')
    ax.set_zlabel('Z Distance')

def randomly_select_carbons_to_be_spins(density, Ncarbons, N_by_three_array):
    #This function will first calculate the number of spins based on Number of carbons and density

    number_of_spins=round(density*Ncarbons)
    [m,n]=np.shape(N_by_three_array)

    #now this function needs to randomly select number_of_spins carbons to be spins

    random_spin_index=[]
    for i in range(0,number_of_spins):
        random_spin_index=random_spin_index+[randint(0,m)]

    returnmatrix=np.zeros((number_of_spins,3))


    for i in range(0, number_of_spins):
        row=np.matrix(N_by_three_array[random_spin_index[i]])
        returnmatrix[i]=row

    return returnmatrix, random_spin_index

#This cell will  #function maxes at around 30,000 carbons
#1. create the N by 3 location array for ideal carbon distribution
#2. make a function that takes in density of spins and randomly distrubtes the spins
#3. Graph those carbons in 3D
#4. modify plot of carbons to include those spins

#1. create the N by 3 location array for ideal carbon distrubtion:
average_radius=1
Number_of_carbons=1000
density=1/1000
N_by_three_array=N_by_three_function(Number_of_carbons, density, average_radius)

#2. make a function that takes in the desnity of spins and randomly distrubtes the spins
matrix_of_spin_loc, random_spin_index=randomly_select_carbons_to_be_spins(density, Number_of_carbons, N_by_three_array)
print("The spins are located at these distances: [x, y, z]")
print(matrix_of_spin_loc)

#3. Graph those carbons in 3D
plot_of_carbons_and_spins=plot_three_d_carbons(N_by_three_array, matrix_of_spin_loc, random_spin_index)
