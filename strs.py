import strengths as strn
import numpy as np
import strengths.plot as strnplt

system = strn.load_rdsystem("system.json")

# setting A and B quantities chemostated to 200 and 0 on
# one side of the system, and to 0 and 200 on the other.

w = system.space.w
h = system.space.h

for x in range(system.space.w) :
    system.set_state("A", (x, 0,   0),     200)
    system.set_state("B", (x, 0,   0),     0  )
    system.set_state("A", (x, h-1, 0),     0  )
    system.set_state("B", (x, h-1, 0),     200)

    system.set_chemostat("A", (x, 0,   0), True)
    system.set_chemostat("B", (x, 0,   0), True)
    system.set_chemostat("A", (x, h-1, 0), True)
    system.set_chemostat("B", (x, h-1, 0), True)

# now we perform the simulation and plot the results
# exactly as we have done in the previous example

output = strn.simulate(
    system = system,
    t_sample = strn.UnitArray([0, 100, 200, 400], "s"),
    time_step = 1,
    engine = strn.engine_collection.tauleap_engine(),
    )

for sample in range(output.nsamples()) :
    strnplt.plot_sample_state_2D(output, "A", sample)

