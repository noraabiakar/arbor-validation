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

cell.soma.insert('ichan2')
cell.soma.gnatbar_ichan2 = 0.120 * in_param["gnatbar_ichan2"]
cell.soma.gkfbar_ichan2  = 0.016 * in_param["gkfbar_ichan2"]
cell.soma.gksbar_ichan2  = 0.006 * in_param["gksbar_ichan2"]
cell.soma.gl_ichan2    = 0.00004 * in_param["gl_ichan2"]
cell.soma.el_ichan2     =          in_param["el_ichan2"]

cell.soma.insert('borgka')
cell.soma.gkabar_borgka  = 0.001 * in_param["gkabar_borgka"]

cell.soma.insert('nca')
cell.soma.gncabar_nca    = 0.001 * in_param["gncabar_nca"]

cell.soma.insert('lca')
cell.soma.glcabar_lca    = 0.005 * in_param["glcabar_lca"]

cell.soma.insert('cat')
cell.soma.gcatbar_cat = 0.000037 * in_param["gcatbar_cat"]

cell.soma.insert('gskch')
cell.soma.gskbar_gskch   = 0.001 * in_param["gskbar_gskch"]

cell.soma.insert('cagk')
cell.soma.gkbar_cagk    = 0.0006 * in_param["gkbar_cagk"]

cell.soma.insert('ccanl')
cell.soma.catau_ccanl    = 10     * in_param["catau_ccanl"]
cell.soma.caiinf_ccanl   = 5.0e-6 * in_param["caiinf_ccanl"]

cell.soma.cm                = 1.0 * in_param["cm_mult"]
cell.soma.Ra                =       in_param["ra"]

syn = h.Exp2Syn(cell.soma(in_param["syn_loc"]))
syn.tau1 = in_param["tau1_syn"]
syn.tau2 = in_param["tau2_syn"]
syn.e = in_param["e_syn"]

cell.soma.ek        = in_param["ek"]
cell.soma.cao       = in_param["cao"]

cell.soma.enca      = 0
cell.soma.elca      = in_param["elca"]
cell.soma.etca      = in_param["etca"]
cell.soma.esk       = in_param["esk"]
cell.soma.enat      = in_param["enat"]
cell.soma.ekf       = in_param["ekf"]
cell.soma.eks       = in_param["eks"]

################################
# Create spike times for input #
################################
tstop = 200 # unts: ms
frequency = 5 # units: Hz

vecstims = h.VecStim()
evecs = h.Vector()

intervals = []

for spike in in_param["spikes"]:
    if spike > tstop:
        break
    evecs.append(spike)
vecstims.play(evecs)

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
h.steps_per_ms = 1/h.dt
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
