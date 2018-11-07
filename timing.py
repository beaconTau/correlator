from __future__ import print_function
import matplotlib.pyplot as plt
from scipy import signal, fftpack
import numpy as np
import pickle
import sys
import time
from matplotlib.colors import LogNorm
from collections import defaultdict
import beacon_data_reader

run_num =191
folder = 'run'+str(191)+'/'

#d = beacon_data_reader.Reader("/project2/avieregg/beacon/telem/root",run_num)
dt=2
#times=d.t()

def delay_finder(volt0,volt1):
    center = len(volt0)
    volt0 = volt0/volt0.std()
    volt1 = volt1/volt1.std()
    cor = np.asarray(signal.correlate(volt0,volt1))
    delay = (np.argmax(cor)-(np.size(cor)+1)/2.)*dt
    return(delay)

def setup(run_num, ent_num, pol):

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
    
    return (volt0,volt1,volt2,volt3,times)

#import timing delays, h and v pol
h1 = np.load(folder+'d1_1.npy')
h2 = np.load(folder+'d2_1.npy')
h3 = np.load(folder+'d3_1.npy')
h4 = np.load(folder+'d4_1.npy')
h5 = np.load(folder+'d5_1.npy')
h6 = np.load(folder+'d6_1.npy')
#v1 = np.load(folder+'d7_1.npy')
#v2 = np.load(folder+'d8_1.npy')
#v3 = np.load(folder+'d9_1.npy')
#v4 = np.load(folder+'d10_1.npy')
#v5 = np.load(folder+'d11_1.npy')
#v6 = np.load(folder+'d12_1.npy')

#length = np.size(h1)

max1 = 22.18+10
max2 = 35.07+10
max3 = 41.29+10
max4 = 24.73+10
max5 = 42.84+10
max6 = 22.90+10

counter = 0

h1p = []
h2p = []
h3p = []
h4p = []
h5p = []
h6p = []
totalh= []
pulse1 = []
a=defaultdict(list)
test=[]
dict_list = []
err = 2

counter=1

event_counter = []
cluster_counter = []

#cut out events that aren't possible based on geometry
for i in range(9000):

    #print(i)
    #volt0,volt1,volt2,volt3,times= setup(run_num,i,0)
    d1=h1[i]#delay_finder(volt0,volt1)
    d2=h2[i]#delay_finder(volt0,volt2)
    d3=h3[i]#delay_finder(volt0,volt3)
    d4=h4[i]#delay_finder(volt1,volt2)
    d5=h5[i]#delay_finder(volt1,volt3)
    d6=h6[i]#delay_finder(volt2,volt3)
    if (np.abs(d1)<max1 and np.abs(d2)<max2 and np.abs(d3)<max3 and np.abs(d4)<max4 and np.abs(d5)<max5 and np.abs(d6)<max6):
        h6p.append(d6)
        h3p.append(d3)
        h1p.append(d1)
        h2p.append(d2)
        h4p.append(d4)
        h5p.append(d5)
        totalh.append(d1+d2+d3+d4+d5+d6)
        for j in range(len(h1p)):
            success=0
            if(abs(h1p[j]-d1)<=err and abs(h2p[j]-d2)<=err and abs(h3p[j]-d3)<=err and abs(h4p[j]-d4)<=err and abs(h5p[j]-d5)<=err and abs(h6p[j]-d6)<=err):

                #print('Event',i,'matches earlier event',j)
                event_counter.append(i)
                cluster_counter.append(j)
                #print(i,j)
                a[j].append(i)
                #print(a[j])
                if j not in dict_list:
                    dict_list.append(j)
                success=1
                break
            else:
                alte=1
                #print('here!!')
        if(j==len(h1p)-1):
            #print('here!')
            event_counter.append(i)
            cluster_counter.append(i)
            
            #else:
                #print('this is the first pulse of its kind!')
print(dict_list)
num_pulses=[]
print('total number of good pulses out of 96000:', counter)
d.setEntry(1)
length_wf = len(d.wf(0))



print(len(dict_list))
print(length_wf)
print('test')
for i in dict_list:#each i is a different pulse
    #print(len(a[i]))
    #a[i].append(i)
    if len(a[i])>100:
        #print("these pulses are clustered: ", a[i])
        entry = a[i][5]
        volt0 = np.zeros([len(a[i]),length_wf])
        d.setEntry(entry)
        counter = 0
        for j in a[i]:
            #if (j == 33):
                #print('the pulse label is:', i)
            d.setEntry(j)
            volt0[counter,:]=np.asarray(d.wf(0)-sum(d.wf(0))/len(d.wf(0)))
            counter = counter+1
            #print(h1[j],h2[j],h3[j])
        #volt2 = d.wf(2)
        #volt4 = d.wf(4)
        #volt6 = d.wf(6)
        print('Pulse ',i,'is found in events:',a[i][0],a[i][1],a[i][2])
        added_v0=np.zeros([length_wf])
        plt.figure(i)
        #if (i==126):
            #print(a[i])
            #print('that was 126!')
        for j in range(0,len(a[i])):
            #print('test',j)
            corr = np.asarray(signal.correlate(volt0[0,:],volt0[j,:]))
            delay=((np.argmax(corr)-(np.size(corr)+1)/2.)*dt)
            offset = int((delay/dt))
            #print(offset)
            if(offset>0):
                added_v0 = added_v0[:(len(volt0[0,:])-offset)]+volt0[j,offset:]
                added_v0 = np.append(added_v0,np.linspace(0,0,offset))
            else:
                #print((len(volt0[0,:])-abs(offset)))
                added_v0 = added_v0[abs(offset):]+volt0[j,:(len(volt0[0,:])-abs(offset))]
                added_v0 = np.append(np.linspace(0,0,abs(offset)),added_v0)
            #volt0[0,:]+
            plt.plot(times+delay,volt0[j,:])
            #print(added_v0)
        #plt.plot(times,added_v0/j)
        np.save('run191_pulse'+str(i)+'_averaged.npy',added_v0)
                    
        #plt.plot(times,volt0[counter-1,:])
        #plt.plot(times+h1[entry],volt2)
        #plt.plot(times+h2[entry],volt4)
        #plt.plot(times+h3[entry],volt6)
    num_pulses.append(len(a[i]))
plt.show()

#print(num_pulses)
print(max(num_pulses))
#print(a[3][5])
#print(a.items())      
np.save('pulse1.npy',np.asarray(pulse1))
#print(counter)
#print(len(list_pulses))
#print(len(counter))
print(np.max(h6p))
h6p = np.asarray(h6p)
h3p = np.asarray(h3p)
#make plots
#H1, yedges, xedges = np.histogram2d(h1p,h2p, bins=(50,50));#first y, then x
#extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
#fig = plt.figure(1)
#ax = fig.add_subplot(1, 1, 1)
#im = ax.imshow(H1, cmap=plt.cm.jet, extent=extent, norm=LogNorm(), interpolation='none')
#fig.colorbar(im, ax=ax)

#H2, yedges, xedges = np.histogram2d(h1p,h3p, bins=(50,50));#first y, then x
#extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
#fig = plt.figure(2)
#ax = fig.add_subplot(1, 1, 1)
#im = ax.imshow(H2, cmap=plt.cm.jet, extent=extent, norm=LogNorm(), interpolation='none')
#fig.colorbar(im, ax=ax)

#H2, yedges, xedges = np.histogram2d(h1p,h4p, bins=(50,50));#first y, then x
#extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
#fig = plt.figure(3)
#ax = fig.add_subplot(1, 1, 1)
#im = ax.imshow(H2, cmap=plt.cm.jet, extent=extent, norm=LogNorm(), interpolation='none')
#fig.colorbar(im, ax=ax)


#H2, yedges, xedges = np.histogram2d(h5p,h4p, bins=(50,50));#first y, then x
#extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
#fig = plt.figure(4)
#ax = fig.add_subplot(1, 1, 1)
#im = ax.imshow(H2, cmap=plt.cm.jet, extent=extent, norm=LogNorm(), interpolation='none')
#fig.colorbar(im, ax=ax)


#H2, yedges, xedges = np.histogram2d(h2p,h4p, bins=(50,50));#first y, then x
#extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
#fig = plt.figure(5)
#ax = fig.add_subplot(1, 1, 1)
#im = ax.imshow(H2, cmap=plt.cm.jet, extent=extent, norm=LogNorm(), interpolation='none')
#fig.colorbar(im, ax=ax)

#H2, yedges, xedges = np.histogram2d(h2p,h6p, bins=(50,50));#first y, then x
#extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
#fig = plt.figure(6)
#ax = fig.add_subplot(1, 1, 1)
#im = ax.imshow(H2, cmap=plt.cm.jet, extent=extent, norm=LogNorm(), interpolation='none')
#fig.colorbar(im, ax=ax)


#plt.figure(7)
##a = np.histogram(totalh)
#plt.hist(totalh, bins=100)

#plt.show()
