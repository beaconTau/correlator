from __future__ import print_function
import os
import matplotlib.pyplot as plt
from scipy import signal
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


#   To run this code, you must have Ben Strutt's BeaconTau python package installed!
#   If you just run "python generate_correlator.py" it will run the code in "main"
#   If you want to reference these functions in another code, include the following on the top of the code:
#   from generate_correlator import Antenna
#   from generate_correlator import correlator #(or whatever other function you want)



def setup(run_num, ent_num, pol):
    #dd= bt.DataDirectory()
    #r=dd.run(run_num)
    #e = r.get_entry(ent_num)
    #times = r.get_entry(ent_num).times()

    lon = -118.238
    lat = 37.589

    
    d=beacon_data_reader.Reader("/project2/avieregg/beacon/telem/root",run_num)
    d.setEntry(ent_num)
    times =d.t()

    if(pol == 1): #Hpol
        #ADC0 = np.asarray(e.channel(0))
        #ADC1 = np.asarray(e.channel(2))
        #ADC2 = np.asarray(e.channel(4))
        #ADC3 = np.asarray(e.channel(6))
        volt0 = d.wf(0)-sum(d.wf(0))/len(d.wf(0))
        volt1 = d.wf(2)-sum(d.wf(2))/len(d.wf(2))
        volt2 = d.wf(4)-sum(d.wf(4))/len(d.wf(4))
        volt3 = d.wf(6)-sum(d.wf(6))/len(d.wf(6))

    else: #Vpol
        #ADC0 = np.asarray(e.channel(1))
        #ADC1 = np.asarray(e.channel(3))
        #ADC2 = np.asarray(e.channel(5))
        #ADC3 = np.asarray(e.channel(7))
        volt0 = d.wf(1)-sum(d.wf(1))/len(d.wf(1))
        volt1 = d.wf(3)-sum(d.wf(3))/len(d.wf(3))
        volt2 = d.wf(5)-sum(d.wf(5))/len(d.wf(5))
        volt3 = d.wf(7)-sum(d.wf(7))/len(d.wf(7))
    #volt0 = (ADC0-sum(ADC0)/len(ADC0))
    #volt1 = (ADC1-sum(ADC1)/len(ADC1))
    #volt2 = (ADC2-sum(ADC2)/len(ADC2))
    #volt3 = (ADC3-sum(ADC3)/len(ADC3))

    h = d.header()
    event_time = datetime.datetime.utcfromtimestamp(h.readout_time)
    #e2 = datetime.datetime.utcfromtimestamp(1538447820)
    #alt2=90-GetAltitude(lat,lon,e2)
    #print(alt2)
    #print("Run number is", run_num)
    #print("event is:", ent_num)
    #print("time is:", event_time)
    alt=90-GetAltitude(lat,lon,event_time)
    az =GetAzimuth(lat,lon,event_time)+270
    #print(alt,az)

    #loc = coord.EarthLocation(lon=lon*u.deg,
                              #lat=lat*u.deg)
    #print(loc)
    #now = Time(event_time)
    #altaz = coord.AltAz(location=loc, obstime=now)
    #sun = coord.get_sun(now)

    #print(sun.transform_to(altaz).alt)
    #print(sun.transform_to(altaz).az)

    
    return (volt0,volt1,volt2,volt3,times,-alt,-az)



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

def rotation(x,y,z,alpha,beta):#alpha=alt, beta=az
    #xp=cosd(beta)*x-sind(beta)*y
    #yp=cosd(alpha)*sind(beta)*x+cosd(alpha)*cosd(beta)*y-sind(alpha)*z
    #zp=sind(alpha)*sind(beta)*x+sind(alpha)*cosd(beta)*y+cosd(alpha)*z
    xp=cosd(alpha)*cosd(beta)*x-cosd(alpha)*sind(beta)*y+sind(alpha)*z
    yp=sind(beta)*x+cosd(beta)*y
    zp=-sind(alpha)*cosd(beta)*x+sind(beta)*sind(alpha)*y+cosd(alpha)*z
    return(xp,yp,zp)

def rotate_antennas(ant0,ant1,ant2,ant3,alpha,beta):
    ant0.x, ant0.y, ant0.z = rotation(ant0.x,ant0.y, ant0.z,alpha,beta)
    ant1.x, ant1.y, ant1.z = rotation(ant1.x,ant1.y, ant1.z,alpha,beta)
    ant2.x, ant2.y, ant2.z = rotation(ant2.x,ant2.y, ant2.z,alpha,beta)
    ant3.x, ant3.y, ant3.z = rotation(ant3.x,ant3.y, ant3.z,alpha,beta)
    return(ant0,ant1,ant2,ant3)



#Loop over thetas and phis:

def generate_time_arrays(A0,A1,A2,A3,alt,az):

    theta_vec = np.linspace(0,180,60).astype(int)
    phi_vec = np.linspace(-180,180,120)
    
    t1 = np.zeros([len(theta_vec),len(phi_vec)])
    t2 = np.zeros([len(theta_vec),len(phi_vec)])
    t3 = np.zeros([len(theta_vec),len(phi_vec)])
    t4 = np.zeros([len(theta_vec),len(phi_vec)])
    t5 = np.zeros([len(theta_vec),len(phi_vec)])
    t6 = np.zeros([len(theta_vec),len(phi_vec)])

    
    A0,A1,A2,A3=rotate_antennas(A0,A1,A2,A3,alt,az)
    
    #Calculate time delays for each correlation box:

    
    for i in range(len(theta_vec)):
        for j in range(len(phi_vec)):
            norm_vector = get_norm(theta_vec[i],phi_vec[j])#norm vector at A0=

            #Start with Antenna 0. Calculate dot product of norm with v1,v2,v3.
            t1[i,j]=-time_delay(asub(A1,A0),norm_vector) #time in nanoseconds between A0 and A1
            t2[i,j]=-time_delay(asub(A2,A0),norm_vector) #A0 and A2
            t3[i,j]=-time_delay(asub(A3,A0),norm_vector) #A0 and A3
            t4[i,j]=-time_delay(asub(A2,A1),norm_vector) #A1 and A2
            t5[i,j]=-time_delay(asub(A3,A1),norm_vector) #A1 and A3
            t6[i,j]=-time_delay(asub(A3,A2),norm_vector) #A2 and A3

    return(t1,t2,t3,t4,t5,t6)

#def correlator(volt0,volt1,volt2,volt3,t1,t2,t3,t4,t5,t6):
def correlator(run_num,ent_num,pol,A0,A1,A2,A3):

    theta_vec = np.linspace(0,180,60).astype(int)
    phi_vec = np.linspace(-180,180,120)

    
    volt0,volt1,volt2,volt3,times,alt,az=setup(run_num,ent_num,pol)
    volt0 = volt0/volt0.std()
    volt1 = volt1/volt1.std()
    volt2 = volt2/volt2.std()
    volt3 = volt3/volt3.std()
    
    dt = times[1]-times[0]


    if(np.max(volt0)<20):
        print(run_num,ent_num)
        t1,t2,t3,t4,t5,t6 = generate_time_arrays(A0,A1,A2,A3,-alt,az)
        center = len(times)
    
        cor1 = np.asarray(signal.correlate(volt0,volt1))
        cor2 = np.asarray(signal.correlate(volt0,volt2))
        cor3 = np.asarray(signal.correlate(volt0,volt3))
        cor4 = np.asarray(signal.correlate(volt1,volt2))
        cor5 = np.asarray(signal.correlate(volt1,volt3))
        cor6 = np.asarray(signal.correlate(volt2,volt3))

        d1=cor1[np.rint((t1/dt+center)).astype(int)]#0 and 1
        d2=cor2[np.rint((t2/dt+center)).astype(int)]#0 and 2
        d3=cor3[np.rint((t3/dt+center)).astype(int)]#0 and 3
        d4=cor4[np.rint((t4/dt+center)).astype(int)]#1 and 2
        d5=cor5[np.rint((t5/dt+center)).astype(int)]#1 and 3
        d6=cor6[np.rint((t6/dt+center)).astype(int)]#2 and 3

        corr_val=np.mean(np.array([d1,d2,d3,d4,d5,d6]),axis=0)
        #print(corr_val)
        return(corr_val)
    else:
        return(np.zeros([len(theta_vec),len(phi_vec)]))



def main():

    #Set these variables before running:

    
    run_num = 191
    ent_num = 33
    pol = 0 #1 for Hpol, 0 for Vpol
    #os.environ["BEACON_DATA_DIR"]="/project2/avieregg/beacon/telem/raw"
    d = beacon_data_reader.Reader("/project2/avieregg/beacon/telem/root",run_num)

    
    #Antenna Positions
    A0 = Antenna(0,0,0)
    A1 = Antenna(-6.039,-1.618,2.275)#east, north, elevation
    A2 = Antenna(-1.272,-10.362,1.282)
    A3 = Antenna(3.411,-11.897,-0.432)

    
    
    #t1,t2,t3,t4,t5,t6 = generate_time_arrays(A0,A1,A2,A3)

    corr_val1 = np.zeros([60,120])
    corr_val2 = np.zeros([60,120])
    corr_val3 = np.zeros([60,120])
    corr_val4 = np.zeros([60,120])
    tot_corr = np.zeros([60,120])

    #Day runs with more than 100,000 events: 191 (17-20), 172(19-22), 174,177, 179,190
    #Day runs with more than 10,000 events: 173, 175, 178,192

    #Night runs with more than 100,000 events: 180,181,182, 183, 186, 187,188,189
    for i in range(0,500):
        val =i*200
        corr_val1=corr_val1+correlator(191,val,pol,A0,A1,A2,A3)
        corr_val2=corr_val2+correlator(172,val,pol,A0,A1,A2,A3)
        corr_val3=corr_val3+correlator(174,val,pol,A0,A1,A2,A3)
        corr_val4=corr_val4+correlator(177,val,pol,A0,A1,A2,A3)
        tot_corr = tot_corr+corr_val1+corr_val2+corr_val3+corr_val4

    fig=plt.figure(1)
    ax = fig.add_subplot(1,1,1)
    im = ax.imshow(corr_val1, interpolation='none', extent=[-180,180,180,0],cmap=plt.cm.coolwarm) #cmap=plt.cm.jet)
    fig.colorbar(im, ax=ax)
    plt.xlabel('Azimuth Angle (Degrees)')
    plt.ylabel('Zenith Angle (Degrees)')

    fig=plt.figure(2)
    ax = fig.add_subplot(1,1,1)
    im = ax.imshow(corr_val2, interpolation='none', extent=[-180,180,180,0],cmap=plt.cm.coolwarm) #cmap=plt.cm.jet)
    fig.colorbar(im, ax=ax)
    plt.xlabel('Azimuth Angle (Degrees)')
    plt.ylabel('Zenith Angle (Degrees)')

    fig=plt.figure(3)
    ax = fig.add_subplot(1,1,1)
    im = ax.imshow(corr_val3, interpolation='none', extent=[-180,180,180,0],cmap=plt.cm.coolwarm) #cmap=plt.cm.jet)
    fig.colorbar(im, ax=ax)
    plt.xlabel('Azimuth Angle (Degrees)')
    plt.ylabel('Zenith Angle (Degrees)')

    fig=plt.figure(4)
    ax = fig.add_subplot(1,1,1)
    im = ax.imshow(corr_val4, interpolation='none', extent=[-180,180,180,0],cmap=plt.cm.coolwarm) #cmap=plt.cm.jet)
    fig.colorbar(im, ax=ax)
    plt.xlabel('Azimuth Angle (Degrees)')
    plt.ylabel('Zenith Angle (Degrees)')

    fig=plt.figure(5)
    ax = fig.add_subplot(1,1,1)
    im = ax.imshow(tot_corr, interpolation='none', extent=[-180,180,180,0],cmap=plt.cm.coolwarm) #cmap=plt.cm.jet)
    fig.colorbar(im, ax=ax)
    plt.xlabel('Azimuth Angle (Degrees)')
    plt.ylabel('Zenith Angle (Degrees)')
    
    plt.show()

if __name__=="__main__":
    main()
