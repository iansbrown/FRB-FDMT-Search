## Written by Bruce Wu

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

# define functions
def t_pulse(t_2, f_2, f, DM):
    """
    Function for generating the time values of a burst given the frequency
    values, according to the dispersion equation.
    t_pulse, t_2 in s
    f_2, f in MHz
    DM in pc cm^-3
    """
    return t_2 + 4.148808*DM*((f/1000.)**(-2) - (f_2/1000.)**(-2))/1000.

def mod_FDMT(im):
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

def calculate_SNR(row):
    """
    Calculates the peak SNR for a time series de-dispersed at one DM. The
    parameters may need tweaking.
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
    return row[ind]/rms, t[ind] # SNR = peak 
            # over rms. Also return when the peak occurs

#pos = [0, 256, 512, 768, 1024]
counter = 0 # total number of detections
path="default"
obsid="1165925976"
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
filepath = os.path.join(os.getcwd(),path,obsid) # fits file containing data with axes y,x,f,t
with fits.open('%s/%sI.fits'%(filepath,obsid)) as hdulist:
    sub_image = hdulist[0].data[:,:,:,:] # 8:520, 8:583
N_t, nf, N_x, N_y   = np.shape(sub_image) # nf is not power of two

# data parameters
N_f = 128 # padded up to nearest power of two
f_min = 169.755 # MHz
f_max = 210.395
dt = 0.5 # s
df = 0.32 # MHz
t = np.arange(N_t)*dt
f = np.arange(N_f)*df + f_min
const = 4.148808*((f_min/1000.)**(-2) - (f_max/1000.)**(-2)) # for converting 
                                                    # delay to DM
DMs = np.arange(N_f)*dt*1000./const # maximum delay is maximum bins
                         # spanned by pulse as inputed into FDMT
                         # multiplied by the timestep, then turned into
                         # a DM
dDM = DMs[1] - DMs[0] # DM step

# pad frequency axis to power of two with zeros
image = np.zeros((N_t, N_f, N_x, N_y), dtype='float64')
image[:,:nf,:,:] = sub_image # fill in nf < N_f bins with data
sub_image = 0

# main program
sigma = 6.5 # S/N detection threshold
y_lst, x_lst, SNR_lst, DM_lst, t_lst = [], [], [], [], [] # store detections

for i in xrange(N_y):
    for j in xrange(N_x): # loop through positions on sky
        im = image[:,:,j,i]
        A = mod_FDMT(im) # take modular FDMT
        for k in xrange(31, np.size(DMs)): # loop through DMs > 300
            SNR, t_max = calculate_SNR(A[k,:]) # calculate SNR
            if SNR > sigma: # if SNR is greater than detection threshold
                y_lst.append(i)
                x_lst.append(j)
                SNR_lst.append(SNR)
                DM_lst.append(DMs[k])
                t_lst.append(t_max) # store results

outdir = ('%s'%(filepath)) # fill this in
num = len(y_lst) # number of detections

# save figures
tm = N_t*dt # for modulus purposes
print (num)
for i in xrange(num):
    im = image[y_lst[i],x_lst[i],:,:]
    plt.figure(figsize=(18,8))
    plt.imshow(im, origin='lower', cmap='Greys_r', interpolation='nearest', \
           extent=(t[0], t[-1]+dt, f_min-df/2., \
                   f_max+df/2.), aspect='auto')
    t_pmax = t_pulse(t_lst[i], f_min, f, DM_lst[i]) # superimpose
    plt.plot((t_pmax + 5.*dt) % tm, f, 'r') # dispersion curve, shifted 5 bins to 
    # the right to not block data in case it is hard to see, also modded
    if t_lst[i] < N_f*dt: # need to show whole image to see modded part
        plt.xlim(0., tm)
    else:
        plt.xlim(max(t_lst[i]-65.,0.), min(t_lst[i]+50.,tm)) # zoom in around pulse
    plt.ylim(f_min-df/2., f_max+df/2.) # set limits or there will be white 
    # space
    plt.colorbar()
    plt.xlabel('Time (s)', size=16)
    plt.ylabel('Frequency (MHz)', size=16)
    plt.title('%s,%s  S/N=%.2f'%((y_lst[i]),(x_lst[i]),SNR_lst[i]), size=18)
    plt.savefig(outdir+'event'+str(counter+i)+'_curve.png')
    plt.close()

    A = mod_FDMT(im)
    plt.figure(figsize=(18,8))
    plt.imshow(A, origin='lower', cmap='hot', interpolation='nearest', \
           extent=(t[0], t[-1]+dt, DMs[0]-dDM/2., DMs[-1]+dDM/2.), \
           aspect='auto')
    plt.plot(t_lst[i]+0.5*dt, DM_lst[i], 'c*', markersize=8.)
    plt.xlim(max(t_lst[i]-65.,0.), min(t_lst[i]+50.,tm))
    plt.ylim(DMs[0]-dDM/2., DMs[-1]+dDM/2.)
    plt.colorbar()
    plt.xlabel('Time (s)', size=16)
    plt.ylabel('Dispersion Measure (pc cm^-3)', size=16)
    plt.title('%s,%s  S/N=%.2f'%((y_lst[i]),(x_lst[i]),SNR_lst[i]), size=18)
    plt.savefig(outdir+'event'+str(counter+i)+'_FDMT.png')
    plt.close()

counter += num
image = 0

print (counter)
