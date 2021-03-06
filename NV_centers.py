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
    print("number of nv centers")
    print(number_of_nv_centers)

    approximate_length_of_sample=(approximate_cube_root(number_of_nv_centers))
    print("approximatelength of sampe")
    print(approximate_length_of_sample)

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

#just for formating purposes
    #np.set_printoptions(precision=2, suppress=True, linewidth=120)

    return ideal_upper_tri

def plot_histogram(upper_tri_matrix):
    import numpy as np
    #import matplotlib.pyplot as plt
    #import plotly.plotly as py
    #this function will plot a histogram of the distance between NV spins:

    [m,n]=np.shape(upper_tri_matrix)
    array_of_Dij=[]
    #getting every element
    for i in range(0, n): #for every column
        for j in range(0,m): #for every row
            array_of_Dij=array_of_Dij+[upper_tri_matrix[j][i]]

    clean_array=[]
    for i in range(len(array_of_Dij)):
        if (array_of_Dij[i]!=0):
            clean_array=clean_array+[array_of_Dij[i]]

    print(clean_array)


    #plt.hist(clean_array, normed=True, bins=30)
    plt.ylabel('Probability');
    plt.xlabel('D_ij Distance between spins')
    plt.show()


location_of_nv_spins_array=generate_nv_center_location_array(1./10000, 12., .2)
upper_tri=upper_tri_function(location_of_nv_spins_array)
print(upper_tri)
#c=plot_histogram(upper_tri)
