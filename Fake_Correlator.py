from __future__ import print_function
import os
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import sys
import math
#import BeaconTau as bt
import beacon_data_reader

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

font = {
        'weight' : 'bold',
        'size' : 16}

plt.rc('font',**font)

#   To run this code, you must have Ben Strutt's BeaconTau python package installed!
#   If you just run "python generate_correlator.py" it will run the code in "main"
#   If you want to reference these functions in another code, include the following on the top of the code:
#   from generate_correlator import Antenna
#   from generate_correlator import correlator #(or whatever other function you want)

rootdatadir = os.environ("BEACON_ROOT_DATA")

def fake_setup():
    volt = np.load("../BeaconTau/run191_pulse126_averaged.npy")
    

def setup(run_num, ent_num, pol):
    d = beacon_data_reader.Reader("/project2/avieregg/beacon/telem/root",run_num)
    times = d.t()
    d.setEntry(ent_num)
    
    #dd= bt.DataDirectory()
    #r=dd.run(run_num)
    #e = r.get_entry(ent_num)
    #times = r.get_entry(ent_num).times()


    if(pol == 1): #Hpol
        ADC0 = np.asarray(d.wf(0))
        ADC1 = np.asarray(d.wf(2))
        ADC2 = np.asarray(d.wf(4))
        ADC3 = np.asarray(d.wf(6))

    else: #Vpol
        ADC0 = np.asarray(d.wf(1))
        ADC1 = np.asarray(d.wf(3))
        ADC2 = np.asarray(d.wf(5))
        ADC3 = np.asarray(d.wf(7))
    
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

    theta_vec = np.linspace(0,180,61).astype(int)
    phi_vec = np.linspace(-180,180,121)
    
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


def reflection_killer(volt0,volt1,volt2,volt3,times):
    finish=0
    dt = times[1]-times[0]
    #Find width of first pulse:
    power=volt0*volt0
    for i in range(len(times)):
        if(power[i]>500):
            start=i-140/4
            finish=i+400/4
            break

        
    print("i am here!")
    print(times[start],times[finish])
        
    #plt.show()
    fmin =start
    fmax =finish
    volt0=volt0[fmin:fmax]
    volt1=volt1[fmin:fmax]
    volt2=volt2[fmin:fmax]
    volt3=volt3[fmin:fmax]
    times=times[fmin:fmax]
    
    return(volt0,volt1,volt2,volt3,times)

def normalize(array):
    return(2*array/len(array))

#def correlator(volt0,volt1,volt2,volt3,t1,t2,t3,t4,t5,t6):
def correlator(run_num,ent_num,pol,A0,A1,A2,A3):

    theta_vec = np.linspace(0,180,61).astype(int)
    phi_vec = np.linspace(-180,180,121).astype(int)

    
    #volt0,volt1,volt2,volt3,times=fake_setup(run_num,ent_num,pol)
    volt0,volt1,volt2,volt3,times=setup(run_num,ent_num,pol)

    
    dt = times[1]-times[0]

    #volt0,volt1,volt2,volt3,times=reflection_killer(volt0,volt1,volt2,volt3,times)
    volt0_norm=volt0/volt0.std()#max(abs(volt0))
    volt1_norm=volt1/volt1.std()#max(abs(volt1))
    volt2_norm=volt2/volt2.std()#max(abs(volt2))
    volt3_norm=volt3/volt3.std()#max(abs(volt3))
    
    t1,t2,t3,t4,t5,t6 = generate_time_arrays(A0,A1,A2,A3)


    
    center = len(times)
    
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


    
    #corr_val1=np.zeros([len(d1[:,1]),len(d1[1,:])])
    #corr_val2=np.zeros([len(d1[:,1]),len(d1[1,:])])
    corr_val=np.mean(np.array([d1,d2,d3,d4,d5,d6]),axis=0)
    #corr_val=corr_val/np.max(corr_val)

    fig=plt.figure(pol+1)
    ax = fig.add_subplot(1,1,1)
    im = ax.imshow(corr_val, interpolation='none', extent=[-180,180,180,0],cmap=plt.cm.coolwarm) #cmap=plt.cm.jet)
    fig.colorbar(im, ax=ax)
    plt.xlabel('Azimuth Angle (Degrees)')
    plt.ylabel('Zenith Angle (Degrees)')
    plt.title('Correlation Plot')

    
    #print('maximum of corr_val is: ', np.unravel_index(corr_val.argmax(),corr_val.shape))
    a1,a2=np.unravel_index(corr_val.argmax(),corr_val.shape)
    print(corr_val[a1,a2])
    theta_best=theta_vec[a1]
    phi_best=phi_vec[a2]
    b1,b2=a1,a2
    while(theta_best>125):
        corr_val[b1,b2]=0
        b1,b2=np.unravel_index(corr_val.argmax(),corr_val.shape)
        theta_best=theta_vec[b1]
        phi_best=phi_vec[b2]
    #    print('here!')
    t1_best=t1[b1,b2]
    t2_best=t2[b1,b2]
    t3_best=t3[b1,b2]
    t4_best=t4[b1,b2]
    t5_best=t5[b1,b2]
    t6_best=t6[b1,b2]

    print("From the correlation plot:")
    print("Best zenith angle:",theta_best)
    print("Best azimuth angle:",phi_best)
    print('predicted time delays between A0 and A1:', t1_best)
    print('predicted time delays between A0 and A2:', t2_best)
    print('predicted time delays between A0 and A3:', t3_best)

    #corr_val=corr_val1+corr_val2    


    plt.figure(pol+3)
    plt.plot(times, volt0, label="Antenna 0")
    plt.plot(times+t1_best, volt1, label="Antenna 1")
    plt.plot(times+t2_best, volt2, label="Antenna 2")
    plt.plot(times+t3_best, volt3, label="Antenna 3")
    plt.legend()
    plt.title("Aligned Pulses from Best Correlator Bin")

    return(corr_val,b1,b2)

def mini_corr(volt,step):
    
    if(step<0):
        volt0=volt
        volt1=np.append(np.linspace(0,0,abs(step)),volt[:(len(volt)-int(abs(step)))])
    else:
        #volt0=np.append(np.linspace(0,0,abs(step)),volt[:(len(volt)-int(abs(step)))])
        volt0=volt
        volt1=np.append(volt[(int(abs(step))):],np.linspace(0,0,abs(step)))

    cor=normalize(np.asarray(signal.correlate(volt0,volt1)))
    #print(len(volt1))
    #print(len(volt0))
    return(cor,volt1)

def fake_correlator(theta_pos,phi_pos,A0,A1,A2,A3,counter):
    #print('here!')
    theta_vec = np.linspace(0,180,61).astype(int)
    phi_vec = np.linspace(-180,180,121).astype(int)
    #print(theta_vec[0])
    #print(phi_vec[0])
    
    volt=np.load("../BeaconTau/run191_pulse12_averaged.npy")
    times=np.linspace(0,len(volt)*2,len(volt))

    
    center=len(volt)
    dt=2
    
    t1,t2,t3,t4,t5,t6 = generate_time_arrays(A0,A1,A2,A3)

    volt0=volt/volt.std()#max(abs(volt))
    #volt1=[np.linspace(0,0,t1[theta_pos,phi_pos]*2
    step1=round(t1[theta_pos,phi_pos]/2.)
    step2=round(t2[theta_pos,phi_pos]/2.)
    step3=round(t3[theta_pos,phi_pos]/2.)
    step4=round(t4[theta_pos,phi_pos]/2.)
    step5=round(t5[theta_pos,phi_pos]/2.)
    step6=round(t6[theta_pos,phi_pos]/2.)

    print(step1,step2,step3)

    cor1,volt1=mini_corr(volt0,step1)
    cor2,volt2=mini_corr(volt0,step2)
    cor3,volt3=mini_corr(volt0,step3)
    cor4,volt2=mini_corr(volt1,step4)
    cor5,volt3=mini_corr(volt1,step5)
    cor6,volt3=mini_corr(volt2,step6)

    plt.figure(13)
    plt.plot(volt0)
    plt.plot(volt1)
    plt.plot(volt2)
    plt.plot(volt3)

    d1=cor1[np.rint((t1/dt+center)).astype(int)]
    d2=cor2[np.rint((t2/dt+center)).astype(int)]
    d3=cor3[np.rint((t3/dt+center)).astype(int)]
    d4=cor4[np.rint((t4/dt+center)).astype(int)]
    d5=cor5[np.rint((t5/dt+center)).astype(int)]
    d6=cor6[np.rint((t6/dt+center)).astype(int)]

    #fig=plt.figure(16)
    #ax = fig.add_subplot(1,1,1)
    #im = ax.imshow(d1, interpolation='none', extent=[-180,180,180,0],cmap=plt.cm.coolwarm) #cmap=plt.cm.jet)
    #fig.colorbar(im, ax=ax)
    #plt.xlabel('Azimuth Angle (Degrees)')
    #plt.ylabel('Zenith Angle (Degrees)')


    
    corr_val=np.mean(np.array([d1,d2,d3,d4,d5,d6]),axis=0)
    
    #corr_val=corr_val/len(corr_val)

    fig=plt.figure(17)
    ax = fig.add_subplot(1,1,1)
    im = ax.imshow(corr_val, interpolation='none', extent=[-180,180,180,0],cmap=plt.cm.coolwarm) #cmap=plt.cm.jet)
    fig.colorbar(im, ax=ax)
    plt.xlabel('Azimuth Angle (Degrees)')
    plt.ylabel('Zenith Angle (Degrees)')
    plt.title('Correlation Map for Simulated Signal')

    
    #plt.savefig('gif_32/correlator_'+str(counter)+'.png')
    
    a1,a2=np.unravel_index(corr_val.argmax(),corr_val.shape)
    print(a1,a2)
    print(theta_pos,phi_pos)
    if(abs(a1-theta_pos)<=2 or abs(a2-phi_pos)<=2 ):
        print("this one works!")
        print(theta_vec[a1],phi_vec[a2])
        print(corr_val[a1,a2])
    theta_best=theta_vec[a1]
    phi_best=phi_vec[a2]
    while(theta_best>125):
        corr_val[a1,a2]=0
        a1,a2=np.unravel_index(corr_val.argmax(),corr_val.shape)
        theta_best=theta_vec[a1]
        phi_best=phi_vec[a2]
        #print('here!')
    t1_best=t1[a1,a2]
    t2_best=t2[a1,a2]
    t3_best=t3[a1,a2]
    t4_best=t4[a1,a2]
    t5_best=t5[a1,a2]
    t6_best=t6[a1,a2]

    print("From the Fake correlation plot:")
    print("Best zenith angle:",theta_best)
    print("Best azimuth angle:",phi_best)
    print('predicted time delays between A0 and A1:', t1_best)
    print('predicted time delays between A0 and A2:', t2_best)
    print('predicted time delays between A0 and A3:', t3_best)

    #corr_val=corr_val1+corr_val2    


    plt.figure(18)
    plt.plot(times, volt0, label="Antenna 0")
    plt.plot(times+t1_best, volt1, label="Antenna 1")
    plt.plot(times+t2_best, volt2, label="Antenna 2")
    plt.plot(times+t3_best, volt3, label="Antenna 3")
    plt.legend()
    plt.title("Aligned Pulses from Simulated Signal")



    
    
def main():
    theta_vec = np.linspace(0,180,61).astype(int)
    phi_vec = np.linspace(-180,180,121).astype(int)
    
    #Set these variables before running:
    run_num = 191
    ent_num = 43#23
    pol = 1 #1 for Hpol, 0 for Vpol
    #os.environ["BEACON_DATA_DIR"]="/project2/avieregg/beacon/telem/raw"
    
    #Antenna Positions
    A0 = Antenna(0,0,0)
    A1 = Antenna(-6.039,-1.618,2.275)#east, north, elevation
    A2 = Antenna(-1.272,-10.362,1.282)
    A3 = Antenna(3.411,-11.897,-0.432)


    #t1,t2,t3,t4,t5,t6 = generate_time_arrays(A0,A1,A2,A3)
    corr_val,theta_pos,phi_pos=correlator(run_num,ent_num,pol,A0,A1,A2,A3)

    counter=0
    ##for theta_pos in range(0,60):
    #for phi_pos in range(0,40):
    #    for theta_pos in range(0,20):
    #        fake_correlator(theta_pos*3,phi_pos*3,A0,A1,A2,A3,counter)
    #        counter=counter+1

    fake_correlator(theta_pos,phi_pos,A0,A1,A2,A3,counter)
    plt.show()
    #pol=1
    #corr_valH=correlator(run_num,ent_num,pol,A0,A1,A2,A3)


    #fig=plt.figure(9)
    #ax = fig.add_subplot(1,1,1)
    #im = ax.imshow(corr_valV+corr_valH, interpolation='none', extent=[-180,180,180,0],cmap=plt.cm.coolwarm) #cmap=plt.cm.jet)
    #fig.colorbar(im, ax=ax)
    #plt.xlabel('Azimuth Angle (Degrees)')
    #plt.ylabel('Zenith Angle (Degrees)')
    
    #plt.show()
    
if __name__=="__main__":
    main()
