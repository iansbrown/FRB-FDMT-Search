import numpy as np
#Verbose = False

# replace with print statements, if necessary
#def PDB(message, var = ''):
    #if Verbose :
        #print message + str(var)
        
def FDMT_initialization(Image, f_min, f_max, maxDT, dataType):
    """
    Input: Image - power matrix I(f,t)
        f_min,f_max - are the base-band begin and end frequencies.
            The frequencies can be entered in both MHz and GHz, units are factored out in all uses.
        maxDT - the maximal delay (in time bins) of the maximal dispersion.
            Appears in the paper as N_{\Delta}
            A typical input is maxDT = N_f
        dataType - To naively use FFT, one must use floating point types.
            Due to casting, use either complex64 or complex128.
    Output: 3d array, with dimensions [N_f,N_d0,Nt]
            where N_d0 is the maximal number of bins the dispersion curve travels at one frequency bin
    
    For details, see algorithm 1 in Zackay & Ofek (2014)
    """
    # Data initialization is done prior to the first FDMT iteration
    # See Equations 17 and 19 in Zackay & Ofek (2014)

    [F,T] = Image.shape

    deltaF = (f_max - f_min)/float(F)
    deltaT = int(np.ceil((maxDT-1) *(1./f_min**2 - 1./(f_min + deltaF)**2) / (1./f_min**2 - 1./f_max**2)))

    Output = np.zeros([F,deltaT+1,T],dataType)
    Output[:,0,:] = Image
    
    for i_dT in xrange(1,deltaT+1):
        Output[:,i_dT,i_dT:] = Output[:,i_dT-1,i_dT:] + Image[:,:-i_dT]
    return Output

def FDMT_iteration(Input, maxDT, F, f_min, f_max, iteration_num, dataType): #, Verbose = False):
    """
        Input: 
            Input - 3d array, with dimensions [N_f,N_d0,Nt]
            f_min,f_max - are the base-band begin and end frequencies.
                The frequencies can be entered in both MHz and GHz, units are factored out in all uses.
            maxDT - the maximal delay (in time bins) of the maximal dispersion.
                Appears in the paper as N_{\Delta}
                A typical input is maxDT = N_f
            dataType - To naively use FFT, one must use floating point types.
                Due to casting, use either complex64 or complex128.
            iteration num - Algorithm works in log2(Nf) iterations, each iteration changes all the 
                sizes (like in FFT)
        Output: 
            3d array, with dimensions [N_f/2,N_d1,Nt]
        where N_d1 is the maximal number of bins the dispersion curve travels at one output frequency band
        
        For details, see algorithm 1 in Zackay & Ofek (2014)
    """

    input_dims = Input.shape
    output_dims = list(input_dims)
    
    deltaF = 2**(iteration_num) * (f_max - f_min)/float(F)
    dF = (f_max - f_min)/float(F)
    # the maximum deltaT needed to calculate at the i'th iteration
    deltaT = int(np.ceil((maxDT-1) *(1./f_min**2 - 1./(f_min + deltaF)**2) / (1./f_min**2 - 1./f_max**2)))
    #PDB("deltaT = ",deltaT)
    #PDB("N_f = ",F/2.**(iteration_num))
    #PDB('input_dims', input_dims)
    
    output_dims[0] = output_dims[0]/2
    
    
    output_dims[1] = deltaT + 1
    #PDB('output_dims', output_dims)
    Output = np.zeros(output_dims,dataType)
    
    # No negative D's are calculated => no shift is needed
    # If you want negative dispersions, this will have to change to 1+deltaT,1+deltaTOld
    # Might want to calculate negative dispersions when using coherent dedispersion, to reduce the number of 
    # trial dispersions by a factor of 2 (reducing the complexity of the coherent part of the hybrid)
    ShiftOutput = 0
    ShiftInput = 0
    T = output_dims[2]
    F_jumps = output_dims[0]
    
    # For some situations, it is beneficial to play with this correction.
    # When applied to real data, one should carefully analyze and understand the effect of 
    # this correction on the pulse he is looking for (especially if convolving with a specific pulse profile)
    if iteration_num>0:
        correction = dF/2.
    else:
        correction = 0
    for i_F in range(F_jumps):
        
        f_start = (f_max - f_min)/float(F_jumps) * (i_F) + f_min
        f_end = (f_max - f_min)/float(F_jumps) *(i_F+1) + f_min
        f_middle = (f_end - f_start)/2. + f_start - correction
        # it turned out in the end, that putting the correction +dF to f_middle_larger 
        # (or -dF/2 to f_middle, and +dF/2 to f_middle larger)
        # is less sensitive than doing nothing when dedispersing a coherently dispersed pulse.
        # The confusing part is that the hitting efficiency is better with the corrections (!?!).
        f_middle_larger = (f_end - f_start)/2 + f_start + correction
        deltaTLocal = int(np.ceil((maxDT-1) *(1./f_start**2 - 1./(f_end)**2) / (1./f_min**2 - 1./f_max**2)))
        
        for i_dT in range(deltaTLocal+1):
            dT_middle = round(i_dT * (1./f_middle**2 - 1./f_start**2)/(1./f_end**2 - 1./f_start**2))
            dT_middle_index = dT_middle + ShiftInput
            
            dT_middle_larger = round(i_dT * (1./f_middle_larger**2 - 1./f_start**2)\
                                     /(1./f_end**2 - 1./f_start**2))            
                     
            
            dT_rest = i_dT - dT_middle_larger
            dT_rest_index = dT_rest + ShiftInput
            
            i_T_min = 0
            
            i_T_max = dT_middle_larger	
            Output[i_F,i_dT + ShiftOutput,i_T_min:i_T_max] = Input[2*i_F, dT_middle_index,i_T_min:i_T_max]
            
            
            i_T_min = dT_middle_larger
            i_T_max = T
            
            
            
            Output[i_F,i_dT + ShiftOutput,i_T_min:i_T_max] = Input[2*i_F, dT_middle_index,i_T_min:i_T_max] \
            + Input[2*i_F+1, dT_rest_index,i_T_min - dT_middle_larger:i_T_max-dT_middle_larger]
    
    return Output

def FDMT(Image, f_min, f_max, maxDT , dataType): #, Verbose = True):
    """
    This function implements the  FDMT algorithm.
    Input: Input power matrix I(f,t)
           f_min,f_max are the base-band begin and end frequencies.
                   The frequencies should be entered in MHz 
           maxDT - the maximal delay (in time bins) of the maximal dispersion.
                   Appears in the paper as N_{\Delta}
                   A typical input is maxDT = N_f
           dataType - a valid numpy dtype.
                      reccomended: either int32, or int64.
    Output: The dispersion measure transform of the Input matrix.
            The output dimensions are [Input.shape[1],maxDT]
    
    For details, see algorithm 1 in Zackay & Ofek (2014)
    """
    F,T = Image.shape
    f = int(np.log2(F))
    if (F not in [2**i for i in range(1,30)]): # or (T not in [2**i for i in range(1,30)]) :
        raise NotImplementedError("Input dimensions must be a power of 2")

    State = FDMT_initialization(Image,f_min,f_max,maxDT,dataType)
    #PDB('initialization ended')
    
    for i_t in range(1,f+1):
        State = FDMT_iteration(State, maxDT, F, f_min, f_max, i_t, dataType) #, Verbose)
    #PDB('total_time:', time.time() - x)
    [F,dT,T] = State.shape
    DMT = np.reshape(State,[dT,T])
    return DMT
