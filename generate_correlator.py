import os
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import sys
import math
import BeaconTau as bt

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
    dd= bt.DataDirectory()
    r=dd.run(run_num)
    e = r.get_entry(ent_num)
    times = r.get_entry(ent_num).times()


    if(pol == 1): #Hpol
        ADC0 = np.asarray(e.channel(0))
        ADC1 = np.asarray(e.channel(2))
        ADC2 = np.asarray(e.channel(4))
        ADC3 = np.asarray(e.channel(6))

    else: #Vpol
        ADC0 = np.asarray(e.channel(1))
        ADC1 = np.asarray(e.channel(3))
        ADC2 = np.asarray(e.channel(5))
        ADC3 = np.asarray(e.channel(7))
    
    volt0 = (ADC0-sum(ADC0)/len(ADC0))
    volt1 = (ADC1-sum(ADC1)/len(ADC1))
    volt2 = (ADC2-sum(ADC2)/len(ADC2))
    volt3 = (ADC3-sum(ADC3)/len(ADC3))
    return (volt0,volt1,volt2,volt3,times)
        
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


#Loop over thetas and phis:

def generate_time_arrays(A0,A1,A2,A3):

    theta_vec = np.linspace(0,180,60).astype(int)
    phi_vec = np.linspace(-180,180,120)
    
    t1 = np.zeros([len(theta_vec),len(phi_vec)])
    t2 = np.zeros([len(theta_vec),len(phi_vec)])
    t3 = np.zeros([len(theta_vec),len(phi_vec)])
    t4 = np.zeros([len(theta_vec),len(phi_vec)])
    t5 = np.zeros([len(theta_vec),len(phi_vec)])
    t6 = np.zeros([len(theta_vec),len(phi_vec)])

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

    
    volt0,volt1,volt2,volt3,times=setup(run_num,ent_num,pol)
    dt = times[1]-times[0]
    
    t1,t2,t3,t4,t5,t6 = generate_time_arrays(A0,A1,A2,A3)
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

    fig=plt.figure(1)
    ax = fig.add_subplot(1,1,1)
    im = ax.imshow(corr_val, interpolation='none', extent=[-180,180,180,0],cmap=plt.cm.coolwarm) #cmap=plt.cm.jet)
    fig.colorbar(im, ax=ax)
    plt.xlabel('Azimuth Angle (Degrees)')
    plt.ylabel('Zenith Angle (Degrees)')

    #print('maximum of corr_val is: ', np.unravel_index(corr_val.argmax(),corr_val.shape))
    a1,a2=np.unravel_index(corr_val.argmax(),corr_val.shape)
    print(corr_val[a1,a2])
    theta_best=theta_vec[a1]
    phi_best=phi_vec[a2]
    while(theta_best>125):
        corr_val[a1,a2]=0
        a1,a2=np.unravel_index(corr_val.argmax(),corr_val.shape)
        theta_best=theta_vec[a1]
        phi_best=phi_vec[a2]
    t1_best=t1[a1,a2]
    t2_best=t2[a1,a2]
    t3_best=t3[a1,a2]
    t4_best=t4[a1,a2]
    t5_best=t5[a1,a2]
    t6_best=t6[a1,a2]



    print("From the correlation plot:")
    print("Best zenith angle:",theta_best)
    print("Best azimuth angle:",phi_best)
    print('predicted time delays between A0 and A1:', t1_best)
    print('predicted time delays between A0 and A2:', t2_best)
    print('predicted time delays between A0 and A3:', t3_best)


    plt.figure(2)
    plt.plot(times, volt0, label="Antenna 0")
    plt.plot(times+t1_best, volt1, label="Antenna 1")
    plt.plot(times+t2_best, volt2, label="Antenna 2")
    plt.plot(times+t3_best, volt3, label="Antenna 3")
    plt.legend()
    plt.title("Aligned Pulses from Best Correlator Bin")


    d1_best=(np.argmax(cor1)-np.size(cor1)/2.)*dt
    d2_best=(np.argmax(cor2)-np.size(cor2)/2.)*dt
    d3_best=(np.argmax(cor3)-np.size(cor3)/2.)*dt

    print('')
    print('Compare to perfect time delays (not from map):')
    print('A0 and A1:',d1_best)
    print('A0 and A2:',d2_best)
    print('A0 and A3:',d3_best)

    #plt.figure(3)
    #plt.plot(times,volt0)
    #plt.plot(times+d1_best, volt1)
    #plt.plot(times+d2_best, volt2)
    #plt.plot(times+d3_best, volt3)

    plt.show()

def main():

    #Set these variables before running:
    run_num = 191
    ent_num = 33
    pol = 0 #1 for Hpol, 0 for Vpol
    os.environ["BEACON_DATA_DIR"]="/project2/avieregg/beacon/telem/raw"

    #Antenna Positions
    A0 = Antenna(0,0,0)
    A1 = Antenna(-6.039,-1.618,2.275)#east, north, elevation
    A2 = Antenna(-1.272,-10.362,1.282)
    A3 = Antenna(3.411,-11.897,-0.432)

    
    #t1,t2,t3,t4,t5,t6 = generate_time_arrays(A0,A1,A2,A3)
    correlator(run_num,ent_num,pol,A0,A1,A2,A3)

if __name__=="__main__":
    main()
