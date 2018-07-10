#!/usr/bin/env python
## Written by Bruce Wu
## Modified by Ian Brown

# import packages
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
plt.ioff()
#import datetime as dtime
from FDMT_functions import FDMT
from astropy.io import fits
import os
import errno
import timeit
from optparse import OptionParser

parser = OptionParser() #create parser for program options

parser.add_option("-o", "--obsid",type='string',action="store", dest="obsid",help='obsid to be used')#obsid for hypercube
parser.add_option('-p','--path',type='string', action= 'store', dest='path',help='path to creat fits files in', default='FDMT')#path to store FDMT output
parser.add_option('-t','--test',type='string', action= 'store', dest='test',help='True to generate a test pulse', default='False')#set True to create a test pulse
options, args = parser.parse_args()

# define functions
def Test():
    '''generate a single pixel pulse to feed into FDMT for test purposes
    '''
    DM=550
    df=.08
    f0=138.594
    dt=2.0
    tmax=112

    f=f0+(np.arange(384)*df) 

    delay=(4.148808*((((f/1000)**-2)-(f.max()/1000)**-2)*DM))/1000
    D=np.zeros((int((tmax/dt)),len(f)))
    #for i in range(len(f)):
   	# D[:,i]=np.random.normal(0,0.1,224)	
    D[np.int32((np.floor((delay/dt))+10)),np.arange(len(f))]+=1
    nt,nf=np.shape(D)
    D.shape=(nt,nf,1,1)
    #plt.clf()
    #plt.imshow(D.T,origin='lower',aspect=1,vmin=0,vmax=1,extent=(0,tmax,f.min(),f.max()))
    return D
    

def Control_pulse(D,fbins,df=.08,f0=138.594,dt=2.0,DM=550.0, tp=384.0):
    '''generate a pulse and places it into a single spacial pixel in image to ensure non-detection is not a failure of program
    '''
    f=f0+(np.arange(fbins)*df)
    delay=(4.148808*((((f/1000)**-2)-(f.max()/1000)**-2)*DM))/1000
    D[np.arange(len(f)),np.int32((np.floor((delay/dt))))]+=(tp/len(f))
    #plt.clf()
    #plt.imshow(D.T,origin='lower',aspect=1,vmin=0,vmax=1,extent=(0,tmax,f.min(),f.max()))
    return D    

def t_pulse(t_2, f_2, f, DM):
    """
    Function for generating the time values of a burst given the frequency
    values, according to the dispersion equation.
    t_pulse, t_2 in s
    f_2, f in MHz
    DM in pc cm^-3
    """
    return t_2 + 4.148808*DM*((f/1000.)**(-2) - (f_2/1000.)**(-2))/1000.

def mod_FDMT(im,f_min,f_max,df):
    """
    Modular FDMT. Dispersion curves which leave the left boundary of the
    image come back through the right boundary. Achieved by duplicating
    the last N_f time bins of the image to the front, taking the FDMT, then
    cutting the duplicated time bins off. This makes the FDMT output in the
    region where the y index is greater than the x index more reliable.
    """
    # make modular image
    mod_im = np.empty((N_f, N_f+N_t)) # time axis extended by N_f
    mod_im[:,:N_f] = im[:,N_t-N_f:] # duplicate last N_f times in front
    mod_im[:,N_f:] = im # the rest is original image
    # take FDMT and truncate
    A = FDMT(mod_im, f_min-df/2., f_max+df/2., N_f, 'float64')
    return np.delete(A, np.arange(N_f), axis=1)

def calculate_SNR(row,tbins):
    """
    Calculates the peak SNR for a time series de-dispersed at one DM. The
    parameters may need tweaking.
    """
    	
    ind = np.argmax(row)
    box1=row[0:(ind-1)]
    box2=row[(ind+1):tbins]
    noise = np.concatenate((box1, box2))
    rms = np.sqrt(np.mean(noise**2)) # calculate rms	
    """
    ind = np.argmax(row) # time index of highest value
    room = ind-66
    if room >= 0:
        right1 = ind-5
        left1 = room
        right2 = min(N_t, ind+72)
        left2 = min(ind+11, right2)
        box1 = row[left1:right1]
        box2 = row[left2:right2]
        noise = np.concatenate((box1, box2))
    else:
        	right1 = max(0, ind-5)
        	left1 = 0
        	right2 = ind+72
        	left2 = ind+11
        	right3 = N_t
        	left3 = N_t+room-1
        	box1 = row[left1:right1]
        	box2 = row[left2:right2]
        	box3 = row[left3:right3] # wraps around
        	noise = np.concatenate((box3, box1, box2))
    rms = np.sqrt(np.mean(noise**2)) # calculate rms
    """
    return row[ind]/rms, t[ind], rms # SNR = peak # over rms. Also return when the peak occurs

#pos = [0, 256, 512, 768, 1024]
counter = 0 # total number of detections
path=options.path
obsid=options.obsid
start = timeit.default_timer()
test=options.test
'''
for p in xrange(16):
    print p, dtime.datetime.now()
    top = pos[p/4]
    bottom = pos[p/4 + 1]
    left = pos[p%4]
    right = pos[p%4 + 1]
    '''
# load data
#filename = '1113366704_y'+str(top)+'-'+str(bottom)+'_x'+str(left)+'-'+str(right) # fill this in
filepath = os.getcwd()#os.path.join(os.getcwd(),path,obsid) # fits file containing data with axes y,x,f,t

if test == 'True':
    sub_image=Test()
else:    
    with fits.open('%s/%sI.fits'%(filepath,obsid)) as hdulist:
        sub_image = hdulist[0].data[:,:,:,:] # 8:520, 8:583
        print("time to load file:%f")%(timeit.default_timer()-start)
N_t, nf, N_x, N_y   = np.shape(sub_image) # nf is not power of two
print (N_t,nf,N_x,N_y)
# data parameters
n=4
while 2**n<nf:
	n+=1
N_f = 2**n # padded up to nearest power of two
f_min = 138.594 # MHz
dt =2  # s
df = .08 # MHz
f_max = f_min + df*N_f
print (f_min,f_max)
t = np.arange(N_t)*dt
f = np.arange(N_f)*df + f_min
const = 4.148808*((f_min/1000.)**(-2) - (f_max/1000.)**(-2)) # for converting 
                                                    # delay to DM
DMs = np.arange(N_t)*dt*1000./const # maximum delay is maximum bins
                         # spanned by pulse as inputed into FDMT
                         # multiplied by the timestep, then turned into
                         # a DM
dDM = DMs[1] - DMs[0] # DM step


outdir = ('%s/%s/'%(filepath,path))
#create outdir if does not already exist.
if not os.path.exists(os.path.dirname(outdir)):
    try:
        os.makedirs(os.path.dirname(outdir))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise

# pad frequency axis to power of two with zeros
z=timeit.default_timer()
image = np.zeros((N_t, N_f, N_x, N_y), dtype='float64')
image[:,:nf,:,:] = sub_image # fill in nf < N_f bins with data
print("total time to pad using npzeros:%f")%(timeit.default_timer()-z)
sub_image = 0

# main program
sigma = 6.5  # S/N detection threshold
y_lst, x_lst, SNR_lst, DM_lst, t_lst, n_lst, DM_lst, TI_lst= [], [], [], [], [], [], [], [] # store detections
tm = N_t*dt # for modulus purposes
'''
outdir = ('%s/%s/'%(filepath,path))
#create outdir if does not already exist.
if not os.path.exists(os.path.dirname(outdir)):
    try:
        os.makedirs(os.path.dirname(outdir))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise'''
disp=[150,550,950,1350]
x=np.random.randint(0,N_x,size=(100))
y=np.random.randint(0,N_x,size=(100))
intensity=np.linspace (0,384,10)

for i in disp:
    for TI in intensity:
        for j in xrange(x):
            
            imtemp = image[:,:,x[j],y[j]]
            img=np.transpose(imtemp,(1,0))
            img = Control_pulse(img,nf,df,f_min,dt,i,TI)
            med_pow = np.median(img,axis=1) #calculate median power for each fine channel 
            im=img - (med_pow[:,None]) #subtract median power to remove constant signals
            A = FDMT(im, f_min, f_max, (N_t), 'float64')#mod_FDMT(im,f_min,f_max,df) # take modular FDMT
            a,b=np.where(A==np.max(A))
            
            for k in a:# loop through DMs where A is at max value
                SNR, t_max, noise = calculate_SNR(A[k,:],N_t) # calculate SNR
                y_lst.append(i)
                x_lst.append(j)
                SNR_lst.append(SNR)
                DM_lst.append(DMs[k])
                t_lst.append(t_max) # store results
                n_lst.append(noise)
                DM_lst.append(i)
                TI_lst.append(TI)
                counter+=1
                plt.figure(figsize=(18,8))
                plt.imshow(im[:nf,:], origin='lower', cmap='Greys_r', interpolation='nearest', \
                       extent=(t[0], t[-1]+dt, f_min-df/2., \
                               f_max+df/2.), aspect='auto')
                t_pmax = t_pulse(t_max, f_min, f, DMs[k]) # superimpose
                plt.plot((t_pmax + 5.*dt) % tm, f, 'r') # dispersion curve, shifted 5 bins to 
                # the right to not block data in case it is hard to see, also modded
                if t_max < N_f*dt: # need to show whole image to see modded part
                    plt.xlim(0., tm)
                else:
                    plt.xlim(max(t_max-65.,0.), min(t_max+50.,tm)) # zoom in around pulse
                plt.ylim(f_min-df/2., f_max+df/2.) # set limits or there will be white 
                # space
                plt.colorbar()
                plt.xlabel('Time (s)', size=28)
                plt.ylabel('Frequency (MHz)', size=28)
                plt.title('Pixel(%s,%s)  S/N=%.2f DM=%.2f'%(i,j,SNR,DMs[k]), size=28)
                plt.savefig(outdir+'event'+str(counter)+'_curve.png')
                plt.close()
            
                #A = mod_FDMT(im)
                plt.figure(figsize=(18,8))
                plt.imshow(A, origin='lower', cmap='hot', interpolation='nearest', \
                       extent=(t[0], t[-1]+dt, DMs[0]-dDM/2., DMs[-1]+dDM/2.), \
                       aspect='auto')
                plt.plot(t_max+0.5*dt, DMs[k], 'c*', markersize=8.)
                plt.xlim(max(t_max-65.,0.), min(t_max+50.,tm))
                plt.ylim(DMs[0]-dDM/2., DMs[-1]+dDM/2.)
                plt.colorbar()
                plt.xlabel('Time (s)', size=28)
                plt.ylabel('Dispersion Measure (pc cm^-3)', size=28)
                plt.title('Pixel (%s,%s)  S/N=%.2f  DM=%.2f'%(i,j,SNR,DMs[k]), size=28)
                plt.savefig(outdir+'event'+str(counter)+'_FDMT.png')
                plt.close()
	#print("total time for SNR step:%f")%(timeit.default_timer()-Stime)
	
file1 = open(outdir+'/testevents.txt','w')
for l in xrange(len(y_lst)):
    q=('%s,%s,%s,%s,%s,%s,%s,%s,%s\n'% ((l+1),SNR_lst[l],n_lst[l],DM_lst[l],t_lst[l],y_lst[l],x_lst[l],DM_lst[l],TI_lst[l]))
    file1.write(q) 
file1.flush()
file1.close() 
    
    
print (counter)
print("total time for program:%f")%(timeit.default_timer()-start)
