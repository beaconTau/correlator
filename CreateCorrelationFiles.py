from __future__ import print_function
#import BeaconTau as bt
import os
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import signal, fftpack
import numpy as np
import pickle
import sys
import time
import beacon_data_reader



#os.environ["BEACON_DATA_DIR"]="/project2/avieregg/beacon/telem/raw" #Ben's analysis

start_time = time.time()

print(start_time)
run_num=191


if (len(sys.argv)>1):
    event_min = int(sys.argv[1])
    event_max = int(sys.argv[2])
    run=int(sys.argv[3])
else:
    event_min = 0
    event_max = 1000
    run=int('1')

print('Event Minimum is: ', event_min)
print('Event Maximum is: ', event_max)
print('run number is: ', run)


#Ben's analysis:
#dd = bt.DataDirectory()
#r = dd.run(141)

#Cosmin's analysis:
d = beacon_data_reader.Reader("/project2/avieregg/beacon/telem/root",run_num)
d.setEntry(45)
h = d.header())

#r.draw('trigger_thresholds')
h_delay =[]
v_delay =[]
event_num = []

#quit()

#loop over multiple entries

num_events = 131000
leng = 1535 #length of time array. THIS CAN CHANGE

c1 = np.zeros((event_max-event_min,leng))
c2 = np.zeros((event_max-event_min,leng))
c3 = np.zeros((event_max-event_min,leng))
c4 = np.zeros((event_max-event_min,leng))
c5 = np.zeros((event_max-event_min,leng))
c6 = np.zeros((event_max-event_min,leng))
c7 = np.zeros((event_max-event_min,leng))
c8 = np.zeros((event_max-event_min,leng))
c9 = np.zeros((event_max-event_min,leng))
c10 = np.zeros((event_max-event_min,leng))
c11 = np.zeros((event_max-event_min,leng))
c12 = np.zeros((event_max-event_min,leng))

d1 = np.zeros((event_max-event_min))
d2 = np.zeros((event_max-event_min))
d3 = np.zeros((event_max-event_min))
d4 = np.zeros((event_max-event_min))
d5 = np.zeros((event_max-event_min))
d6 = np.zeros((event_max-event_min))
d7 = np.zeros((event_max-event_min))
d8 = np.zeros((event_max-event_min))
d9 = np.zeros((event_max-event_min))
d10 = np.zeros((event_max-event_min))
d11 = np.zeros((event_max-event_min))
d12 = np.zeros((event_max-event_min))

for i in range(0,event_max):#141740
    if(i % 1000 ==0):
        num= 100.-float(event_max-i)/float(event_max-event_min)*100.
        print('Run is ', num, 'percent finished')

    t1 = time.time()
    #e = r.get_entry(i) #Ben's analysis package
    d.setEntry(i) #Cosmin's analysis package

    #ADC1 = np.asarray(e.channel(1))#Ben's analysis package
    #ADC2 = np.asarray(e.channel(2))
    #ADC3 = np.asarray(e.channel(3))
    #ADC4 = np.asarray(e.channel(4))
    #ADC5 = np.asarray(e.channel(5))
    #ADC6 = np.asarray(e.channel(6))
    #ADC7 = np.asarray(e.channel(7))
    #ADC0 = np.asarray(e.channel(0))

    volt0 = d.wf(0)-sum(d.wf(0))/len(d.wf(0))#Cosmin's analysis package
    volt1 = d.wf(1)-sum(d.wf(1))/len(d.wf(1))
    volt2 = d.wf(2)-sum(d.wf(2))/len(d.wf(2))
    volt3 = d.wf(3)-sum(d.wf(3))/len(d.wf(3))
    volt4 = d.wf(4)-sum(d.wf(4))/len(d.wf(4))
    volt5 = d.wf(5)-sum(d.wf(5))/len(d.wf(5))
    volt6 = d.wf(6)-sum(d.wf(6))/len(d.wf(6))
    volt7 = d.wf(7)-sum(d.wf(7))/len(d.wf(7))

    #print(ADC0)

 
    times = d.t()
    dt = times[1]-times[0]
    t2 = time.time()

    t1 = time.time()


    volt0 = volt0/volt0.std()
    volt1 = volt1/volt1.std()
    volt2 = volt2/volt2.std()
    volt3 = volt3/volt3.std()
    volt4 = volt4/volt4.std()
    volt5 = volt5/volt5.std()
    volt6 = volt6/volt6.std()
    volt7 = volt7/volt7.std()


    #Horizontal Cross Correlation
    c1[i-event_min,:]=np.asarray(signal.correlate(volt0,volt2))
    c2[i-event_min,:]=np.asarray(signal.correlate(volt0,volt4))
    c3[i-event_min,:]=np.asarray(signal.correlate(volt0,volt6))
    c4[i-event_min,:]=np.asarray(signal.correlate(volt2,volt4))
    c5[i-event_min,:]=np.asarray(signal.correlate(volt2,volt6))
    c6[i-event_min,:]=np.asarray(signal.correlate(volt4,volt6))


    
    #Vertical Cross Correlation
    c7[i-event_min,:]=np.asarray(signal.correlate(volt1,volt3))
    c8[i-event_min,:]=np.asarray(signal.correlate(volt1,volt5))
    c9[i-event_min,:]=np.asarray(signal.correlate(volt1,volt7))
    c10[i-event_min,:]=np.asarray(signal.correlate(volt3,volt5))
    c11[i-event_min,:]=np.asarray(signal.correlate(volt3,volt7))
    c12[i-event_min,:]=np.asarray(signal.correlate(volt5,volt7))


    #Record best fit delay in these arrays:
    d1[i-event_min]=(np.argmax(c1[i-event_min,:])-(np.size(c1[i-event_min,:]))/2.)*dt
    d2[i-event_min]=(np.argmax(c2[i-event_min,:])-np.size(c2[i-event_min,:])/2.)*dt
    d3[i-event_min]=(np.argmax(c3[i-event_min,:])-np.size(c3[i-event_min,:])/2.)*dt
    d4[i-event_min]=(np.argmax(c4[i-event_min,:])-np.size(c4[i-event_min,:])/2.)*dt
    d5[i-event_min]=(np.argmax(c5[i-event_min,:])-np.size(c5[i-event_min,:])/2.)*dt
    d6[i-event_min]=(np.argmax(c6[i-event_min,:])-np.size(c6[i-event_min,:])/2.)*dt

    
    d7[i-event_min]=(np.argmax(c7[i-event_min,:])-np.size(c7[i-event_min,:])/2.)*dt
    d8[i-event_min]=(np.argmax(c8[i-event_min,:])-np.size(c8[i-event_min,:])/2.)*dt
    d9[i-event_min]=(np.argmax(c9[i-event_min,:])-np.size(c9[i-event_min,:])/2.)*dt
    d10[i-event_min]=(np.argmax(c10[i-event_min,:])-np.size(c10[i-event_min,:])/2.)*dt
    d11[i-event_min]=(np.argmax(c11[i-event_min,:])-np.size(c11[i-event_min,:])/2.)*dt
    d12[i-event_min]=(np.argmax(c12[i-event_min,:])-np.size(c12[i-event_min,:])/2.)*dt


    t2 = time.time()

    
elapsed_time = time.time()-start_time

print("total time passed: ", elapsed_time, " seconds")

folder = 'run'+str(run_num)+'/'

np.save(folder+'c1_'+str(run)+'.npy',c1)
np.save(folder+'c2_'+str(run)+'.npy',c2)
np.save(folder+'c3_'+str(run)+'.npy',c3)
np.save(folder+'c4_'+str(run)+'.npy',c4)
np.save(folder+'c5_'+str(run)+'.npy',c5)
np.save(folder+'c6_'+str(run)+'.npy',c6)
np.save(folder+'c7_'+str(run)+'.npy',c7)
np.save(folder+'c8_'+str(run)+'.npy',c8)
np.save(folder+'c9_'+str(run)+'.npy',c9)
np.save(folder+'c10_'+str(run)+'.npy',c10)
np.save(folder+'c11_'+str(run)+'.npy',c11)
np.save(folder+'c12_'+str(run)+'.npy',c12)


np.save(folder+'d1_'+str(run)+'.npy',d1)
np.save(folder+'d2_'+str(run)+'.npy',d2)
np.save(folder+'d3_'+str(run)+'.npy',d3)
np.save(folder+'d4_'+str(run)+'.npy',d4)
np.save(folder+'d5_'+str(run)+'.npy',d5)
np.save(folder+'d6_'+str(run)+'.npy',d6)
np.save(folder+'d7_'+str(run)+'.npy',d7)
np.save(folder+'d8_'+str(run)+'.npy',d8)
np.save(folder+'d9_'+str(run)+'.npy',d9)
np.save(folder+'d10_'+str(run)+'.npy',d10)
np.save(folder+'d11_'+str(run)+'.npy',d11)
np.save(folder+'d12_'+str(run)+'.npy',d12)

