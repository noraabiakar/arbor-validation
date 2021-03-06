# Test multiple synaptic activation of cell models

from neuron import h
import numpy as np
import time as cookie
import pylab as plt
import pickle
import json
import sys

with open(sys.argv[1]) as json_file:
    in_param = json.load(json_file)

np.random.seed(149)

h.load_file("stdrun.hoc")
h.load_file("cell.hoc")

##################
# Creating cells #
##################

cell = h.mkcell()
if in_param["syn_seg"] == 0 :
    syn = h.Exp2Syn(cell.soma(in_param["syn_loc"]))
else :
    syn = h.Exp2Syn(cell.dend(in_param["syn_loc"]))

syn.tau1 = in_param["tau1_syn"]
syn.tau2 = in_param["tau2_syn"]
syn.e = in_param["e_syn"]

if in_param["soma_hh"] :
    cell.soma.insert("hh")
    cell.soma.ena = in_param["hh_ena"]
    cell.soma.ek = in_param["hh_ek"]
    cell.soma.gnabar_hh = in_param["hh_gnabar"]
    cell.soma.gkbar_hh = in_param["hh_gkbar"]
    cell.soma.gl_hh = in_param["hh_gl"]
else :
    cell.soma.insert("pas")
    cell.soma.e_pas = in_param["pas_e"]
    cell.soma.g_pas = in_param["pas_g"]

if in_param["dend_hh"] :
    cell.dend.insert("hh")
    cell.dend.ena = in_param["hh_ena"]
    cell.dend.ek = in_param["hh_ek"]
    cell.dend.gnabar_hh = in_param["hh_gnabar"]
    cell.dend.gkbar_hh = in_param["hh_gkbar"]
    cell.dend.gl_hh = in_param["hh_gl"]
else :
    cell.dend.insert("pas")
    cell.dend.e_pas = in_param["pas_e"]
    cell.dend.g_pas = in_param["pas_g"]

################################
# Create spike times for input #
################################
tstop = 200 # unts: ms
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

for v in vec:
    print v

#####################
# Connecting inputs #
#####################
nc = h.NetCon(vecstims, syn)
nc.weight[0] = in_param["weight"]
nc.delay = 0

################################
# Setting up vectors to record #
################################
v = h.Vector()
v.record(cell.soma(0.5)._ref_v)

t = h.Vector()
t.record(h._ref_t)

#########################
# Setting up simulation #
#########################
h.v_init = in_param["vinit"]
h.t = 0
h.dt = in_param["dt_neuron"]
h.celsius = in_param["temp"]
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
