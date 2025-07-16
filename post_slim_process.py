import tskit
import pyslim
import msprime

ts = tskit.load("fixed_sweep.trees")

# SLiM apparently creates an older version of ts than what pyslim uses. So I will update it.
ts = pyslim.update(ts)

metadata = ts.metadata['SLiM']['user_metadata']

Ne = metadata['L1'][0] * metadata['L2'][0] * metadata['rho'][0] # There might be extra correction as in Wakeley's paper, but I don't remember the formula.
# Trace back until we find the MRCA, extending the tree sequence from SLiM simulation
rts = pyslim.recapitate(ts, ancestral_Ne=Ne, recombination_rate=metadata['r'][0], random_seed=6)

# Sprinkle neutral mutations on top of the recapitated tree sequence
next_id = pyslim.next_slim_mutation_id(rts)
rts = msprime.sim_mutations(
            rts, rate=1e-8, random_seed=7, keep=True,
            model=msprime.SLiMMutationModel(type=0, next_id=next_id)
)

# Todo : continue simulation, as in https://tskit.dev/pyslim/docs/stable/vignette_continuing.html

