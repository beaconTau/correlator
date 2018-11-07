from __future__ import print_function
import os
import matplotlib.pyplot as plt
from scipy import signal
from scipy.fftpack import rfft, irfft, fftfreq
from scipy.signal import iirfilter,lfilter,butter
from scipy import signal
from scipy.fftpack import fft
import numpy as np
import sys
import math
#import BeaconTau as bt
import beacon_data_reader
from Pysolar.solar import *
import datetime
import astropy.coordinates as coord
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
import time

class Antenna:
    def __init__(self,x0,y0,z0):
        self.x=x0
        self.y=y0
        self.z=z0

class Vector:
    def __init__(self,x0,y0,z0):
        self.x=x0
        self.y=y0
        self.z=z0



def setup(run_num, ent_num, pol): #grab waveforms and time of event
    #Ben's Analysis Code:
    #dd= bt.DataDirectory()
    #r=dd.run(run_num)
    #e = r.get_entry(ent_num)
    #if(pol == 1): #Hpol
    #    ADC0 = np.asarray(e.channel(0))
    #    ADC1 = np.asarray(e.channel(2))
    #    ADC2 = np.asarray(e.channel(4))
    #    ADC3 = np.asarray(e.channel(6))
    #else:
    #    ADC0 = np.asarray(e.channel(1))
    #    ADC1 = np.asarray(e.channel(3))
    #    ADC2 = np.asarray(e.channel(5))
    #    ADC3 = np.asarray(e.channel(7))
        
    #volt0 = (ADC0-sum(ADC0)/len(ADC0))
    #volt1 = (ADC1-sum(ADC1)/len(ADC1))
    #volt2 = (ADC2-sum(ADC2)/len(ADC2))
    #volt3 = (ADC3-sum(ADC3)/len(ADC3))
    #times = r.get_entry(ent_num).times()

    lon = -118.238
    lat = 37.589

    d=beacon_data_reader.Reader("/project2/avieregg/beacon/telem/root",run_num)
    d.setEntry(ent_num)
    times =d.t()

    if(pol == 1): #Hpol
        volt0 = d.wf(0)-sum(d.wf(0))/len(d.wf(0))
        volt1 = d.wf(2)-sum(d.wf(2))/len(d.wf(2))
        volt2 = d.wf(4)-sum(d.wf(4))/len(d.wf(4))
        volt3 = d.wf(6)-sum(d.wf(6))/len(d.wf(6))

    else: #Vpol
        volt0 = d.wf(1)-sum(d.wf(1))/len(d.wf(1))
        volt1 = d.wf(3)-sum(d.wf(3))/len(d.wf(3))
        volt2 = d.wf(5)-sum(d.wf(5))/len(d.wf(5))
        volt3 = d.wf(7)-sum(d.wf(7))/len(d.wf(7))

    #filter events:
    volt0 = filter(volt0,times)
    volt1 = filter(volt1,times)
    volt2 = filter(volt2,times)
    volt3 = filter(volt3,times)

    #grab time stamp of event and format for later calculations
    h = d.header()
    event_time = datetime.datetime.utcfromtimestamp(h.readout_time)
    hour=event_time.hour+event_time.minute/60.0+event_time.second/3600.0
    loc = coord.EarthLocation(lon=lon*u.deg,lat=lat*u.deg)
    time=Time(event_time,location=loc)#-offset
    alt,az=0,0
    
    return (volt0,volt1,volt2,volt3,times,alt,az,time,loc,hour)





def get_norm(theta, phi): #negatives because wave travels towards antennas, so normal vector should be opposite direction from angles
    x=math.cos(math.radians(phi))*math.sin(math.radians(theta))*-1
    y=math.sin(math.radians(phi))*math.sin(math.radians(theta))*-1
    z=math.cos(math.radians(theta))*-1
    return Vector(x, y, z)

#time delay calculation.
#postive time delay: hits vector "smaller valued" antenna first
#negative time delay: hits second "larger valued" antenna first
def time_delay(ant,vec):
    #print((ant.x*vec.x+ant.y*vec.y+ant.z*vec.z)/3e8*1e9)
    return (ant.x*vec.x+ant.y*vec.y+ant.z*vec.z)/3e8*1e9

def asub(ant2, ant1):
    return Antenna(ant2.x-ant1.x, ant2.y-ant1.y, ant2.z-ant1.z)


def cosd(angle):
    return math.cos(math.radians(angle))

def sind(angle):
    return math.sin(math.radians(angle))

#def rotation(x,y,z,alpha,beta):#alpha=alt, beta=az. No longer using this
#    xp=cosd(alpha)*cosd(beta)*x-cosd(alpha)*sind(beta)*y+sind(alpha)*z
#    yp=sind(beta)*x+cosd(beta)*y
#    zp=-sind(alpha)*cosd(beta)*x+sind(beta)*sind(alpha)*y+cosd(alpha)*z
#    return(xp,yp,zp)

#def rotate_antennas(ant0,ant1,ant2,ant3,alpha,beta):
#    antA=Antenna(0,0,0)
#    antB=Antenna(0,0,0)
#    antC=Antenna(0,0,0)
#    antD=Antenna(0,0,0)
    
#    antA.x, antA.y, antA.z = rotation(ant0.x,ant0.y, ant0.z,alpha,beta)
#    antB.x, antB.y, antB.z = rotation(ant1.x,ant1.y, ant1.z,alpha,beta)
#    antC.x, antC.y, antC.z = rotation(ant2.x,ant2.y, ant2.z,alpha,beta)
#    antD.x, antD.y, antD.z = rotation(ant3.x,ant3.y, ant3.z,alpha,beta)
#    return(antA,antB,antC,antD)



#Loop over thetas and phis:

def generate_time_arrays(A0,A1,A2,A3,alt,az):
    low = 91
    high = 181
    
    theta_vec = np.linspace(0,180,low).astype(int)
    phi_vec = np.linspace(-180,180,high).astype(int)
    #theta_vec = np.linspace(0,90,31).astype(int)
    #phi_vec = np.linspace(-90,90,61).astype(int)
    
    t1 = np.zeros([len(theta_vec),len(phi_vec)])
    t2 = np.zeros([len(theta_vec),len(phi_vec)])
    t3 = np.zeros([len(theta_vec),len(phi_vec)])
    t4 = np.zeros([len(theta_vec),len(phi_vec)])
    t5 = np.zeros([len(theta_vec),len(phi_vec)])
    t6 = np.zeros([len(theta_vec),len(phi_vec)])
    
    A0p,A1p,A2p,A3p=A0,A1,A2,A3

    
    #Calculate time delays for each correlation box:

    
    for i in range(len(theta_vec)):
        for j in range(len(phi_vec)):
            norm_vector = get_norm(theta_vec[i],phi_vec[j])#norm vector at A0=

            #Start with Antenna 0. Calculate dot product of norm with v1,v2,v3.
            t1[i,j]=-time_delay(asub(A1p,A0p),norm_vector) #time in nanoseconds between A0 and A1
            t2[i,j]=-time_delay(asub(A2p,A0p),norm_vector) #A0 and A2
            t3[i,j]=-time_delay(asub(A3p,A0p),norm_vector) #A0 and A3
            t4[i,j]=-time_delay(asub(A2p,A1p),norm_vector) #A1 and A2
            t5[i,j]=-time_delay(asub(A3p,A1p),norm_vector) #A1 and A3
            t6[i,j]=-time_delay(asub(A3p,A2p),norm_vector) #A2 and A3

    return(t1,t2,t3,t4,t5,t6)

def normalize(array):
    return(array/len(array))




#def correlator(volt0,volt1,volt2,volt3,t1,t2,t3,t4,t5,t6):
def correlator(run_num,ent_num,pol,old_corr_val,tracker,t1,t2,t3,t4,t5,t6,hour_vector):

    low = 91
    high = 181
    
    #print(ent_num)
    theta_vec = np.linspace(0,180,low).astype(int)
    phi_vec = np.linspace(-180,180,high).astype(int)
    
    volt0,volt1,volt2,volt3,times,alt,az,e_time,loc,hour=setup(run_num,ent_num,pol)

    hour_vector.append(hour)
    dt = times[1]-times[0]
    window_num=1.0
    window_size = int(round(len(volt0)/window_num))
    radec_total = np.zeros([low,high])
    tracker_total=np.zeros([low,high])

    
    for i in range(0,int(window_num)):
        low_bound = i*window_size
        high_bound = i*window_size+int(window_size)
        
        volt0_norm= volt0[low_bound:high_bound]/volt0[low_bound:high_bound].std()
        volt1_norm= volt1[low_bound:high_bound]/volt1[low_bound:high_bound].std()
        volt2_norm= volt2[low_bound:high_bound]/volt2[low_bound:high_bound].std()
        if(volt3[low_bound:high_bound].std()!=0.0):
            volt3_norm = volt3[low_bound:high_bound]/volt3[low_bound:high_bound].std()
        else:
            return(old_corr_val,tracker)

        
        center = len(volt0_norm)


        single=1
        if (single==0):
        
            cor1 = normalize(np.asarray(signal.correlate(volt0_norm,volt1_norm)))
            cor2 = normalize(np.asarray(signal.correlate(volt0_norm,volt2_norm)))
            cor3 = normalize(np.asarray(signal.correlate(volt0_norm,volt3_norm)))
            cor4 = normalize(np.asarray(signal.correlate(volt1_norm,volt2_norm)))
            cor5 = normalize(np.asarray(signal.correlate(volt1_norm,volt3_norm)))
            cor6 = normalize(np.asarray(signal.correlate(volt2_norm,volt3_norm)))

            d1=cor1[np.rint((t1/dt+center)).astype(int)]#0 and 1
            d2=cor2[np.rint((t2/dt+center)).astype(int)]#0 and 2
            d3=cor3[np.rint((t3/dt+center)).astype(int)]#0 and 3
            d4=cor4[np.rint((t4/dt+center)).astype(int)]#1 and 2
            d5=cor5[np.rint((t5/dt+center)).astype(int)]#1 and 3
            d6=cor6[np.rint((t6/dt+center)).astype(int)]#2 and 3

            corr_val_fun=np.mean(np.array([d1,d2,d3,d4,d5,d6]),axis=0)
            corr_val2 = corr_val_fun

        else:
            cor1 = normalize(np.asarray(signal.correlate(volt0_norm,volt3_norm)))
            d1=cor1[np.rint((t3/dt+center)).astype(int)]#0 and 1
            corr_val2 = d1
            
        #corr_val_fun=corr_val_fun/np.max(corr_val_fun)
        #corr_val_fun = d2/np.max(d2)

        
        #good:
        radec_map,tracker=radec_plotter2(corr_val2,e_time,loc,hour,old_corr_val,tracker)
        old_corr_val = radec_map
        #radec_plotter2(corr_val2,e_time,loc,hour,old_corr_val,tracker)
        
    return(radec_map,tracker,hour_vector)

def radec_plotter2(corr_val1,e_time,loc,hour,new_map,tracker):#converts corr_val1 (the correlation map) from alt/az to RA and DEC

    low = 91 #61
    high = 181 #121
    degree = 180.0/float((low-1))
    
    RA_vec = np.linspace(0,360,high).astype(int)
    DEC_vec = np.linspace(90,-90,low).astype(int)

    #radec_map = np.zeros([61,121])
    #tracker = np.zeros([61,121])
    radec_map = new_map
    lon = -118.238
    lat = 37.589

    d=e_time.mjd-51544.5
    LST=(100.46+(0.985647*d)+lon+15*hour)%360.0

    #print(len(DEC_vec),len(RA_vec))
    for i in range(0,len(DEC_vec)):
        for j in range(0,len(RA_vec)):
            DEC = DEC_vec[i]
            RA = RA_vec[j]
            HA = LST-RA
            if(HA<0):
                HA=HA+360

            alt = np.arcsin(sind(DEC)*sind(lat)+cosd(DEC)*cosd(lat)*cosd(HA))*180.0/np.pi
            a = np.arccos((sind(DEC)-sind(alt)*sind(lat))/(cosd(alt)*cosd(lat)))*180.0/np.pi
            if sind(HA)<0:
                az=a
            else:
                az=360-a
            #print(az,alt)
            if(alt>0 and az>0 and az<180):#is it visible by BEACON?240
                
                y = (((low-1)/2)-int(alt/degree))%low
                x = (int(az/degree)+(high-1)/4)%high
                x_low = (x-1)%high
                x_high = (x+1)%high
                y_low = (y-1)
                y_high = (y+1)
                if (y==high):
                    y_high = y
                if (y==low):
                    y_low = y
                #print(i,j)
                
                #if (x>=121):
                #    x = x-121
                
                #print(x,y)
                if(tracker[i,j]==0):
                    radec_map[i,j]=(corr_val1[y,x]+corr_val1[y_low,x]+corr_val1[y_high,x]+corr_val1[y,x_low]+corr_val1[y,x_high])/6.0
                    tracker[i,j]=1.0
                else:
                    radec_map[i,j]=radec_map[i,j]+(corr_val1[y,x]+corr_val1[y_low,x]+corr_val1[y_high,x]+corr_val1[y,x_low]+corr_val1[y,x_high])/6.0
                    #print(y,x)
                    tracker[i,j]=tracker[i,j]+1.0
    return(radec_map,tracker)
                
                
    

def radec_plotter(corr_val1,e_time,loc,hour,new_map,tracker):
    #only want top half of corr_val
    alt_vec = np.linspace(90,0,31).astype(int)
    phi_vec = np.linspace(0,180,61).astype(int)

    lon = -118.238
    lat = 37.589


    
    for i in range(0,len(alt_vec)):
        for j in range(0,len(phi_vec)):
            #print(i,j)
            alt=alt_vec[i]
            phi=phi_vec[j]

            d=e_time.mjd-51544.5
 
            DEC=asind(sind(alt)*sind(lat)+cosd(alt)*cosd(lat)*cosd(phi))
            #another attempt at H and DEC:
            z=90.0-alt
            az=phi
            H2 = np.arcsin(-sind(z)*sind(az)/cosd(DEC))*180.0/np.pi
            arg2 = (cosd(lat)*cosd(z)-sind(lat)*sind(z)*cosd(az))/cosd(DEC)
            if(arg2>1.0):# or arg2<-1.0):
                arg2 = 1.0
            if(arg2<-1.0):
                arg2 = -1.0
            H1 = np.arccos(arg2)*180.0/np.pi
            if(H1==H2):
                H=H1
            if(H1>=90.0 and H2>=0.0):
                H=H1
            if(H1>=90.0 and H2<=0.0):
                H=360.0-H1
            if(H1<=90.0 and H2<=0.0):
                H=360.0-H1
            #print(H)


            LST=(100.46+(0.985647*d)+lon+15*hour)%360.0

            RA=(LST-H)%360.0

            y = int(round((90.0-DEC)/3.0)) #should go from 0 to 60
            x = int(round(RA/3.0)) #should go from 0 to 120
            j_new = j+30

            if(new_map[y,x]==0):
                #print("hello!")
                new_map[y,x]=corr_val1[i,j_new]

            else:
                tracker[y,x]=tracker[y,x]+1.0
                new_map[y,x]=(new_map[y,x]+corr_val1[i,j_new])

    return(new_map,tracker)

def sind(angle):
    return(math.sin(math.radians(angle)))

def cosd(angle):
    return(math.cos(math.radians(angle)))

def asind(number):
    return(math.degrees(math.asin(number)))

def acosd(number):
    if number>1.0:
        number=1.0
    if number<-1.0:
        number=-1.0
    a = math.acos(number)
    return(math.degrees(a))

def filter(data,time):
    
    fs = 500000000.0  # Sample frequency (Hz)
    Q = 30.0  # Quality factor
    f0 = 55000000.0 #frequency to notch (55 MHz)
    w0= f0/(fs/2)
    c, d = signal.iirnotch(w0, Q)
    filtered_data1 = lfilter(c,d,data)
    cutoff = 80.0e6
    nyq = 0.5 * fs
    order = 5.0
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    filtered_data2 = lfilter(b, a, filtered_data1)
    
    return(filtered_data2)


def main():
    start = time.time()# for timing purposes

    #Dimensions of correlation array:
    low = 91
    high = 181
    
    #Set these variables before running:
    tracker = np.zeros([low,high])
    run_num = 240
    ent_num = 0
    pol = 1 #1 for Hpol, 0 for Vpol
    number_events = 400.0 #total number of events
    plot = 1 #do you want a plot? 1 = yes, 0 = no
    
    #Ben's Code: Uncomment This!
    #os.environ["BEACON_DATA_DIR"]="/project2/avieregg/beacon/telem/raw"
    
    #Antenna Positions
    A0 = Antenna(0,0,0)
    A1 = Antenna(-6.039,-1.618,2.275)#east, north, elevation compared to A0
    A2 = Antenna(-1.272,-10.362,1.282)
    A3 = Antenna(3.411,-11.897,-0.432)

    hour_vec = []
    t1,t2,t3,t4,t5,t6 = generate_time_arrays(A0,A1,A2,A3,0.0,0.0)

    corr_val = np.zeros([low,high])

    full_sky = 1 #Do you want to evenly sample across 24 hour period? 1= yes, 0= no

    if(full_sky ==1):
        split = int(400000.0/number_events)
        run_list = [271, 268, 272, 247, 273, 248, 274, 275, 249, 250, 189, 265, 269, 266, 270, 267]#these runs and theevents on the next line were specifically chosen so that they would evenly sample over 24 hours. 
        num_list = [[0,49707],[39590,40683],[0,48965],[8845,10300],[0,48780],[8130,10330],[0,50503],[0,14999],#up to 275
                    [25335,50875],[0,8999],[124100,161500],[0,50229],[9940,50724],[40660,41110],[0,49999],[40020,40900]]
        total = 0
        for j in range(0,len(run_list)):#loop over number of runs
            print('run ', run_list[j])
            
            minval = num_list[j][0]
            maxval = num_list[j][1]
            
            for i in range(minval,maxval):#loop over range of events
                
                if ((total +i-minval)%split==0 and run_list[j]!=189):
                    corr_val,tracker,hour_vec=correlator(run_list[j],i,pol,corr_val,tracker,t1,t2,t3,t4,t5,t6,hour_vec)
                if (run_list[j]==189 and ((total+i-minval)%(split*2048/768))==0):#run 189 has a different time length so you have to make sure it isn't over sampled
                    corr_val,tracker,hour_vec=correlator(run_list[j],i,pol,corr_val,tracker,t1,t2,t3,t4,t5,t6,hour_vec)
            total = total + (maxval-minval)
            
            #np.save('partial_corr_val_1.npy',corr_val)
            #np.save('partial_tracker_1.npy',tracker)
    else:
        for i in range(0,1):
            corr_val,tracker,hour_vec = correlator(run_num,i,pol,corr_val,tracker,t1,t2,t3,t4,t5,t6,hour_vec)

            
    a1,a2=np.unravel_index(tracker.argmax(),tracker.shape)

    tracker[tracker == 0] = 1.0 #get rid of 'nan' values
    corr_val=np.divide(corr_val,tracker)

    np.save('corr_val.npy',corr_val)# save to file

    finish = time.time()
    print(finish-start,'seconds to finish', len(hour_vec),'events')

    

    if(plot):

        fig1=plt.figure(1)
        ax1 = fig1.add_subplot(1,1,1)
        im1 = ax1.imshow(corr_val, interpolation='none',extent=[0,360,-90,90],cmap=plt.cm.coolwarm) #cmap=plt.cm.jet)
        fig1.colorbar(im1, ax=ax1)
        plt.plot([299.75,350.75,266.25],[40,58,-29],marker='o',markersize=5,color="green")
        plt.xlim([0,360])
        plt.ylim([-90,90])
        plt.ylabel('Declination (Degrees)')
        plt.xlabel('Right Ascension (Degrees)')
        plt.savefig('complete_plot2.jpg')
        #plt.show()
        fig2=plt.figure(2)
        ax2 = fig2.add_subplot(1,1,1)
        im2 = ax2.imshow(tracker, interpolation='none', extent=[0,360,-90,90],cmap=plt.cm.coolwarm) #cmap=plt.cm.jet)
        fig2.colorbar(im2, ax=ax2)
        plt.plot([299.75,350.75,266.25],[40,58,-29],marker='o',markersize=5,color="green")
        plt.xlim([0,360])
        plt.ylim([-90,90])
        plt.ylabel('Declination (Degrees)')
        plt.xlabel('Right Ascension (Degrees)')
        
        #plt.show()
        plt.figure(3)
        plt.plot(hour_vec,marker='o',markersize=5,color='blue')
        plt.show()
    #plt.show()

if __name__=="__main__":
    main()
