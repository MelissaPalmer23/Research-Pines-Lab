import numpy as np
from numpy import linalg
import scipy
import math
import sys

# D = 3.2698e-19
# w0 = 1.9618e+11 Hz
cc = complex(0.0,1.0)
w0=1.9618e+11
D=3.2698e-19
Lx = 2; Ly = 2; Lz = 2
N=Lx*Ly*Lz
sz = np.matrix([[0.50, 0.0], [0.0, -0.50]])

sx = np.matrix([[0.0, 0.50], [0.50, 0.0]])

sy = np.matrix([[0.0, -cc*0.50], [cc*0.50, 0.0]])

I  = np.matrix([[1.0, 0.0], [0.0, 1.0]])



#####################generate position of P1 spins ##########################
# Position of sites
Pos = np.zeros((N,3))
site=0
for i in range(Lx):
	for j in range(Ly):
		for k in range(Lz):
			Pos[site,0]=i; Pos[site,1]=j; Pos[site,2]=k;
#			print site,Pos[site,:]
			site = site + 1
#############################################################################


#####################generate distance between P1 spins ##########################
Dist= np.zeros((N,N))
for i in range(N):
	for j in range(i+1,N):
		Dist[i,j]=r_unit*r*math.sqrt((Pos[i,0]-Pos[j,0])**2+(Pos[i,1]-Pos[j,1])**2+(Pos[i,2]-Pos[j,2])**2)
#		print i,j,Dist[i,j]


######################################################################################

##########generate cos(angle), angle between vector(i,j) and magnetic file ###########
Ang= np.zeros((N,N))
for i in range(N):
        for j in range(i+1,N):
                Ang[i,j]=np.dot([Pos[j,0]-Pos[i,0],Pos[j,1]-Pos[i,1],Pos[j,2]-Pos[i,2]],[0,0,1])/Dis[i,j]
######################################################################################


HA = np.zeros((2**N,2**N))

hz = 1
hy = 1
hx = 1


for i in range(N):
	for j in range(i+1,N):
		hz = 1 ; hy = 1; hx = 1;
		for k in range(N):
			#Sz X Sz
			if k == i:
				Z = sz
				X = sx
				Y = sy;
			elif k == j:
				Z = sz; X = sx; Y = sy;
			else:
				Z = I; X = I; Y = I;

		hz = np.kron(hz,Z);  hy =  np.kron(hy,Y); hx = np.kron(hx,X);
		HA = HA + (3*(Ang[i,j]**2.0)-1.0)/(Dis[i,j]**3.0)*(hz - 1.0/2.0*( hy + hx ) )


# w0 Sz
HSZ = np.zeros((2**(N), 2**(N)))
for i in range(N):
	hz = 1
	for k in range(N):
		if k == i:
			Z = sz;
		else:
			Z = I;
		hz = np.kron(hz,Z);
	HSZ= HSZ + hz
# w0 Sx

T1 = np.zeros((2**N, 2**N))
T2 = np.zeros((2**N, 2**N))
T3 = np.zeros((2**N, 2**N))
INFO = np.zeros(N)

for i in range(3**N):
	O = 1
	for j in range(N):
		INFO[j] = (i//(3**(N-1-j)))%3 - 1
	if sum(INFO) == 1:
		for j in range(N):
			S = (1-abs(INFO[j]))*I + abs(INFO[j])*(sx+INFO[j]*cc*sy)
			O =  np.kron(O,S)
		T1 = T1 +O
	elif sum(INFO) == 2:
		for j in range(N):
			S = (1-abs(INFO[j]))*I + abs(INFO[j])*(sx+INFO[j]*cc*sy)
			O =  np.kron(O,S)
		T1 = T1 +O
	elif sum(INFO) == 3:
		for j in range(N):
			S = (1-abs(INFO[j]))*I + abs(INFO[j])*(sx+INFO[j]*cc*sy)
			O =  np.kron(O,S)
		T3 = T3 +O





#print(HA)
HH = w0*HSZ  + D*HA
w, v = np.linalg.eig(HH)

w = w.real
print (2**N,np.count_nonzero(HH))
#for i in range(2**N):
#	for j in range(2**N):#
#		if HH[i,j] != 0:
			#print i+1,j+1,HH[i,j],0.0
for i in range(2**N):
	print (i,w[i])
for i in range(2**N):
	for j in range(2**N):
		p1=abs(np.matmul(np.matmul(np.transpose(np.conjugate(v[:,i])), T1),v[:,j]))
		p2=abs(np.matmul(np.matmul(np.transpose(np.conjugate(v[:,i])), T2),v[:,j]))
		p3=abs(np.matmul(np.matmul(np.transpose(np.conjugate(v[:,i])), T3),v[:,j]))
		print (w[i]-w[j],p1[0,0],p2[0,0],p3[0,0],p1[0,0]+p2[0,0]+p3[0,0])
##############################################################################

#Melissa code (NV_centers.py)
#all input variables in this file are in terms of floats

def NV_center(density, number_of_nv_desired):
    #This function takes in a density and number of nv centers desired and returns
    #the number carbon samples you will need to sample from in the ideal world

    #Where density is of the form nv/carbons (ie 1./10000)
    return (number_of_nv_desired*(1/density))

def NV_center_approximate(density, number_of_nv_desired, standard_deviation):
    #This function takes in a density and number of nv centers desired and returns
    #the approximate number of carbon samples you will need to sample from in the
    #realistic world

    #all input variables are in terms of floats

    import numpy as np
    #Again get the standard deviation you want to use:
    random_offset=np.random.normal(0, standard_deviation)
    return abs(NV_center(density, number_of_nv_desired+(random_offset)))

def approximate_cube_root(N):
    #this function will take the cube root of a number.
    #this will be used as the length of the cube of spins of the sample
    return N**(1./3)

def generate_nv_center_location_array(density, N, standard_deviation):
    import random
    import numpy as np
    #This function will generate a N by 3 location array for the NV centers
    #where columns 1,2,3 represent x,y,z location respectively in terms of r_average
    #(the average distance between carbons)
    number_of_nv_centers=NV_center_approximate(density, N, standard_deviation)
    approximate_length_of_sample=(approximate_cube_root(number_of_nv_centers))
    coord_array=np.zeros((int(N), 3))

    for i in range(0,int(N)):
        #assign each coordinate to a random integer between 1 and the length of the sample
        #this will be the distance in terms of Radius_average_between_spins
        #that the spin has with respect to (0,0,0);

        coord_array[i][0]=random.randint(1, int(approximate_length_of_sample));
        coord_array[i][1]=random.randint(1, int(approximate_length_of_sample));
        coord_array[i][2]=random.randint(1, int(approximate_length_of_sample));

    return coord_array


#Written last week but copied to this file...
def distance_between_spins(firstlocationarray, secondlocationarray):
    import numpy as np
    #this function takes in two 1x3 arrays representing the x,y,z coord of the spin
    #and will return the distnace between them (classic distance formula)

    x_start=float(firstlocationarray[0])
    y_start=firstlocationarray[1]
    z_start=firstlocationarray[2]

    x_end=secondlocationarray[0]
    y_end=secondlocationarray[1]
    z_end=secondlocationarray[2]
    #distance formula
    a= np.sqrt( (z_end-z_start)**2+ (y_end-y_start)**2+ (x_end-x_start)**2)
    return a

def upper_tri_function(location_array):
    import numpy as np

    #this function generates the NxN upper triangular matrix representing the spins
    #distance to each other
    [N, threee]=np.shape(location_array)
    ideal_upper_tri=np.zeros((N,N))
    #calculate the distance between every spin and set that number in the correct
#position in ideal upper tri matrix
    for i in range(0,N): #all the rows
        for j in range(0,N): #all the columns
            a=distance_between_spins(location_array[j], location_array[i])
            ideal_upper_tri[i][j]=float(a)

    #making the lower triangle all zeros:
    for i in range(0,N):
        for j in range(0,0+i):
            ideal_upper_tri[i][j]=0

    return ideal_upper_tri
    