def populate(dx,dy,d):
    import numpy as np
    """
    
    Assumption:
        The mass of each particle is 1 unit mass
    
    This function generates a 2D array of size dx x dy and randomly populates 
    it with particles of value 1 and empty space of 0. 
    
    The array can be populated to result in a specific density using 'd'. The 
    Governing equation is 
    
    dencity = (number of particles)x(mass of one particle)/(1D volume)


    Parameters
    ----------
    dx : horizontal unit length of the volume section
    dy : virtical unit length of the volume section
    d : expected dencity

    Returns
    -------
    array : an array populated with particles with a dencity of 'd'

    """
    count = 0
    limit =d*(dx*dy)
    nums = np.ones(dx*dy)
    nums[:(dx*dy)-int(d*(dx*dy))] = 0
    np.random.shuffle(nums)
    array = np.reshape(np.array(nums), (-1, dx))
                    
    return array


def animate(array,sample_size):
    """
    This function takes in the generated array and the desired sample size, which
    is the total number of particles + the empty space in the array, and itterates 
    through a wide range of horizontal length which is used to ditermine the sample 
    volume. At each ittiration the dencity is calculated assumeing a mass of a molocule 
    is assumed to be 1 unit mass.
    

    Parameters
    ----------
    array : generated array of the particles and free space that is in the 
            simulated container.
    sample_size : the total number of particles + the empty space in the container

    Returns
    -------
    d : density at each itteration of the horizontal length which is used to ditermine the sample 
        volume
    vv : the horizontal length which is used to ditermine the sample volume
    nn : Number of particle in the sample volume at any 'vv'

    """
    
    import numpy as np
    
    d = []
    vv = []
    nn = []
    for i in np.arange(1,sample_size,1):
        n = len(array[0:i,0:i][array[0:i,0:i]==1])
        v = i*i
        ds = n/v
        vv.append(i)
        nn.append(n)
        d.append(ds)
    return d,vv,nn
        

def time_run(array,dx,dl):
    """
    This function moves the shifts the sample volume area
    from left to right to simulate horizontaly moving particles.
    at every timestep  
    Parameters
    ----------
    array : generated array of the particles and free space that is in the 
            simulated container.
    dx : horizontal unit length of the volume section

    dl : 1D sample volume = dl*dl

    Returns
    -------
    vv : Sample Volume
    nn : Number of particles in Sample Volume
    d :Density of the particles in the sample volume
    t : Time step

    """
    import numpy as np
    
    d = []
    cd = []
    vv = []
    nn = []
    t=[]
    ds=0
    for i in np.arange(0,dx,1):     
        if i+dl < dx:
            n = len(array[:dl,i:i+dl][array[:dl,i:i+dl]==1])
            v = dl**2
            cd.append(n/v-ds)
            ds = n/v
            vv.append(np.sqrt(v))
            nn.append(n)
            d.append(ds)
            t.append(i)
    return vv,nn,d,t,cd

def find_nearest(array, value):    
    import numpy as np
    
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def plot(dx,dy,d,dl,diviation):
    """
    This function generates 3 containers with the same density of particles 
    but the particles are reshuffled. 

    Parameters
    ----------
    dx: horizontal unit length of the volume section
    by: vertical unit length of the volume section
    d: expected density
    dl : 1D sample volume = dl*dl
    
    deviation: adjusting the ylimits of the plots 

    Returns
    -------
    recomended_Dl: the recommended dl that can satisfy the continuum hypothesis

    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    sample_size = min(dx,dy)

    array1 =populate(dx,dy,d)
    d1,vv1,nn1 = animate(array1, sample_size)
    
    
    array2 =populate(dx,dy,d)
    d2,vv2,nn2 = animate(array2, sample_size)
    
    
    array3 =populate(dx,dy,d)
    d3,vv3,nn3 = animate(array3, sample_size)
    
    
    
    vvv,nnv,dv,tv,cdv = time_run (array1,dx,dl)
    vvv2,nnv2,dv2,tv2,cdv2 = time_run (array2,dx,dl)
    vvv3,nnv3,dv3,tv3,cdv3 = time_run (array3,dx,dl)
    
    fig = plt.figure(figsize=(15,7))   
    plt.figure(1)
    
    ax = plt.subplot(212)
    ax.imshow(array3)
    ax.plot((dl,dl),(0,dl),'k')
    ax.plot((0,dl),(dl,dl),'k')
    ax.set_title( r'Sample Size Example')
    
    ax1 = plt.subplot(221)
    ax1.plot(vv1,d1,label='Particle Shuffle 1')
    ax1.plot(vv2,d2,label='Particle Shuffle 2')
    ax1.plot(vv3,d3,label='Particle Shuffle 3')
    ax1.plot((0,max(vv2)),(d,d),'k',label='Dencity of the sample')
    ax1.set_xlim(min(vv1),max(vv1))
    ax1.set_ylim(d-d*diviation,d+d*diviation)
    ax1.set_xlabel('Length Scale')
    ax1.set_ylabel( r'$\rho$ ($n \times m/v$)')
    ax1.legend()
    
    
    ax1 = plt.subplot(222)
    ax1.plot(np.array(tv),np.array(dv),label='Particle Shuffle 1')
    ax1.plot(np.array(tv2),np.array(dv2),label='Particle Shuffle 2')
    ax1.plot(np.array(tv3),np.array(dv3),label='Particle Shuffle 3')
    ax1.plot((0,max(tv2)),(d,d),'k',label='Dencity of the sample')
    ax1.set_ylim(d-d*diviation,d+d*diviation)
    ax1.set_xlim(0,max(tv))
    ax1.set_xlabel('Time Scale')
    ax1.set_ylabel( r'$\rho$ ($n \times m/v$)')
    ax1.legend()
    
    plt.show()
    fig.tight_layout()
    
    perfect_d = np.ones_like(d1)*d
    
    d_matrix = np.abs(np.median(np.array([d1-perfect_d,d2-perfect_d,d3-perfect_d]),axis=0))
    
    recomended_Dl = np.where(d_matrix == find_nearest(d_matrix, 0.0003))[0][0]
    
    return recomended_Dl
