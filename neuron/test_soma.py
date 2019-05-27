# Test multiple synaptic activation of cell models

from neuron import h
import numpy as np
import time as cookie
import pylab as plt
import pickle
np.random.seed(149)

h.load_file("stdrun.hoc")
h.load_file("cell.hoc")

##################
# Creating cells #
##################

cell_soma = h.mkcell()
syn = h.Exp2Syn(cell_soma.soma(0.5))
syn.tau1 = 0.5
syn.tau2 = 1.5
syn.e = 0
 
cell_soma.soma.insert("hh")
cell_soma.soma.ena = 50
cell_soma.soma.ek = -77

################################
# Create spike times for input #
################################
tstop = 10000 # unts: ms
frequency = 5 # units: Hz

vecstims = h.VecStim()
evecs = h.Vector()
vec = []

intervals = []
mu = 1000./frequency # Convert to ms
elapsed_time = 0
flag = 1
while flag:
    roll = np.random.uniform(0,1)
    interval = -mu*np.log(roll)
        
    intervals.append(interval)
    elapsed_time += interval
        
    if elapsed_time > tstop:
        flag = 0
        spikes = np.cumsum(intervals)[:-1]
        for spike in spikes:
            evecs.append(spike)
            vec.append(spike)
           
        vecstims.play(evecs)

#####################
# Connecting inputs #
#####################
nc = h.NetCon(vecstims, syn)
nc.weight[0] = 1.17
nc.delay = 0

################################
# Setting up vectors to record #
################################
v = h.Vector()
v.record(cell_soma.soma(0.5)._ref_v)

t = h.Vector()
t.record(h._ref_t)

#########################
# Setting up simulation #
#########################
h.v_init = -70
h.t = 0
h.dt = 0.001
h.celsius = 35.0
h("tstep = 0")
h("period = 2")
h.tstop = tstop
h("steps_per_ms = 10")
h.load_file('init.hoc')

##################
# Run simulation #
##################
print("Starting...!")
ST = cookie.time()
h.run()
ET = cookie.time()-ST
print("Finished in %f seconds" % ET)

########
# Plot #
########
_=plt.plot(t,v)
_=plt.xlabel('Time (ms)')
_=plt.ylabel('Somatic Voltage (mV)')
plt.show()
