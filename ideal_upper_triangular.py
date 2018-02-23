
def Model_matrix(N):
    import numpy as np
    #this function takes in N (the number of spins) and will generate
    #a Nx3 matrix where the 0th, 1st, and 2nd column represent x, y ,z coor
    #of each spin respectively.

    #its based on generating random numbers pulled from a guassian distribtion
    #centered at zero (the idea is normally there will be no difference between
    #the ideal locations and the slightly off locations at the spin where variance
    #is controlled by the user... the bigger the variance the more difference
    #the model matrix will be from the ideal matrix.

    #basically this random offset gets added to every location from the ideal
    #position of spins


    #calling the ideal matrix of location of spins
    the_ideal_matrix= ideal_matrix(N)
    [m,n]=np.shape(the_ideal_matrix)

    #Things you need from Yihua: (*also change the average distance between spin
    #in later function
    standard_deviation=.2;
    average_distance_between_spints=1;

    #Randomly select a number from a guassian distribtion with variance selected from
    #above and mean zero;

    for i in range(0,m): #for every row
        for j in range(0,n): #for every column
            #this function selects a number from gaussian distribtion centered
            #at zero w some standard deviation
            random_offset=np.random.normal(0, standard_deviation)
            the_ideal_matrix[i,j]=the_ideal_matrix[i,j]+random_offset

    return the_ideal_matrix #which is no longer ideal becuase of the offset






def ideal_matrix(N):

        #ask what this values should be set at
    average_radius=1;

    #this funnction will take in N and will generate a Nx3 ideal matrix
    #representing the x, y, z location of the spin
    #N must be a cube root because this models the spins as a cube
    #ie 8, 27, 64, 125 (2^2, 3^3, 4^3, 5^3...
    #import matplotlib.pyplot as plt
    import numpy as np

    # If N does not have a cube root, print a warning to user and break out of function
    root= N**(1/3)
    if (root!=int(root)):
        print ("N needs to be a cubed number!")
        return

    root=int(root) #converting the type of root from float to int

    #generating Nx3 array with x, y, z coord locations (starting at 0,0,0)
    coord_array=np.zeros((N,3))
    x=0;
    y=0;
    z=0;

#making the ideal matrix! (as in all spins are right on the cube corners)
#a spin on the top right corner should have x position = cube root of N
#y position = cube root of N
#z position = cube root of N
#likewise the other positions are calculated with that idea
    ii=0;
    for i in range(root): #x
        for j in range(root): #y
            for z in range(root): #z
                coord_array[ii, 0]=i*average_radius
                coord_array[ii, 1]=j*average_radius
                coord_array[ii, 2]=z*average_radius

                ii+=1;

    return coord_array


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


def ideal_upper_triangular(N):
    #this function calls a function to generate all the locations (x,y,z) of the spins
    #then turns that nx3 matrix into an upper triangular matrix nxn that
    #holds the distance between each spin (ie location 1, 2 would be
    #the distance between the first and second spin)

    import numpy as np

    # If N does not have a cube root, print a warning to user and break out of function
    root= N**(1/3)
    if (root!=int(root)):
        print ("N needs to be a cubed number!")
        return

    #initilizing a zero matrix
    ideal_upper_tri=np.zeros((N,N))

    #calling the previous function to the location matrix Nx3 w x,y,z position
    #the_ideal_matrix=ideal_matrix(N) #if you want the true ideal use this
    the_ideal_matrix=Model_matrix(N) #if you want the offset included use this

#calculate the distance between every spin and set that number in the correct
#position in ideal upper tri matrix
    for i in range(0,N): #all the rows
        for j in range(0,N): #all the columns
            a=distance_between_spins(the_ideal_matrix[j], the_ideal_matrix[i])
            ideal_upper_tri[i][j]=float(a)

    #making the lower triangle all zeros:
    for i in range(0,N):
        for j in range(0,0+i):
            ideal_upper_tri[i][j]=0

#just for formating purposes
    np.set_printoptions(precision=2, suppress=True, linewidth=120)
    return ideal_upper_tri



#Here is running your function
y=ideal_upper_triangular(8)
print(y)
