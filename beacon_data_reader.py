#! /usr/bin/env python 

# Example for how to read beacon data with python This mostly defines a
# beacon_data_reader class, but if you run it as a script it  will try to plot
# some waveforms from run 141. Normally, you would do something like
#    
#   import beacon_data_reader 
#   d = beacon_data_reader.Reader("/path/to/data",runno)  
#   d.setEntry(someInterestingEntry) 
#   d.wf(0) # get a numpy array of channel 0 



import ROOT
import numpy 
import sys


## You need to load the libbeaconroot library
## For the following line to work, the shared library must be in 
## your dynamic loader path (LD_LIBRARY_PATH or DYLD_LIBRARY_PATH, accordingly) 

ROOT.gSystem.Load("libbeaconroot.so"); 


## A class to read in data from a run
##  
##   e.g. d = beacon_data_reader("/path/to/runs",runno) 
##
##   To select the 100th entry in the run
##     d.setEntry(65079); 
##   Or, you can select by event number instead 
##     d.setEvent(360000065079); 
##  
##   You can then obtain the corresponding header, event or status objects
##     d.header()
##     d.event() 
##     d.status()  (this gets the status with the closest readout time to the event) 
##    
##   For now, you can look at the c++ doc to see what is in each
##  
##  You can also numpy arrays of the waveform values and time using 
## 
##  d.wf(channel)  
##   
##  For convenience, d.t() makes a time array with the right number of samples (the sample rate is always 0.5 GSPS) 

class Reader:

  def __init__(self,base_dir, run):

    self.run = run; 
    self.base_dir = base_dir

    self.event_file = ROOT.TFile("%s/run%d/event.root" % (base_dir, run))
    self.event_tree = self.event_file.Get("event") 
    self.evt = ROOT.beacon.Event() 
    self.event_entry = -1; 
    self.event_tree.SetBranchAddress("event",self.evt) 

    self.head_file = ROOT.TFile("%s/run%d/header.root" % (base_dir, run))
    self.head_tree = self.head_file.Get("header") 
    self.head = ROOT.beacon.Header(); 
    self.head_entry = -1
    self.head_tree.SetBranchAddress("header",self.head) 
    self.head_tree.BuildIndex("header.event_number") 

    self.status_file = ROOT.TFile("%s/run%d/status.root" % (base_dir, run))
    self.status_tree = self.status_file.Get("status") 
    self.stat= ROOT.beacon.Status(); 
    self.status_tree.SetBranchAddress("status",self.stat) 
    self.status_tree.BuildIndex("status.readout_time","status.readout_time_ns"); 
    self.status_entry =-1; 

    self.current_entry = 0; 
    
  def setEntry(self,i): 
    if (i < 0 or i >= self.head_tree.GetEntries()):
      sys.stderr.write("Entry out of bounds!") 
    else: 
      self.current_entry = i; 

  def setEvent(self,i):
    setEntry(self.head_tree.GetEntryNumberWithIndex(i)) 


  def event(self,force_reload = False): 
    if (self.event_entry != self.current_entry or force_reload):
      self.event_tree.GetEntry(self.current_entry)
      self.event_entry = self.current_entry 
    return self.evt 


  def wf(self,ch = 0):  
    ## stupid hack because for some reason it doesn't always report the right buffer length 
    g = self.event().getGraph(ch % 8) 
    return numpy.frombuffer(g.GetY(), numpy.dtype('float64'), self.event().getBufferLength()) - 64 
    g = None #try to forget 

  def t(self):
    return numpy.linspace(0, self.event().getBufferLength()*2, self.event().getBufferLength()) 

  def header(self,force_reload = False): 
    if (self.head_entry != self.current_entry or force_reload): 
      self.head_tree.GetEntry(self.current_entry); 
      self.head_entry = self.current_entry 
    return self.head 

  def status(self,force_reload = False): 
    if (self.status_entry != self.current_entry or force_reload): 
      self.status_tree.GetEntry(self.status_tree.GetEntryNumberWithBestIndex(self.header().readout_time[0], self.header().readout_time_ns[0]))
      self.status_entry = self.current_entry

    return self.stat


  def N(self): 
    return self.head_tree.GetEntries() 


# Plot some waveforms is run as a binary 
if __name__=="__main__": 

## pyplot for plotting
  import matplotlib.pyplot as plt

# If your data is elsewhere, pass it as an argument
  datapath = sys.argv[1] if len(sys.argv) > 1 else "/data/beacon/root"

# look at run 141
  d = Reader(datapath,141) 
# this is a random event
  d.setEntry(53) 

## dump the headers and status, just to show they're there
  d.header().Dump(); 
  d.status().Dump(); 
  print d.N() 

# plot all waveforms
  for i in range(8): 
    plt.subplot(2,4,i+1); 
    plt.plot(d.t(), d.wf(i))
  
  plt.show() 




  








