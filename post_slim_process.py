import tskit
import pyslim
import msprime
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import sys

ts = tskit.load("fixed_sweep.trees")
print(f"ts loaded at {datetime.now()}")
# SLiM apparently creates an older version of ts than what pyslim uses. So I will update it.
ts = pyslim.update(ts)
print(f"ts updated at {datetime.now()}")

metadata = ts.metadata['SLiM']['user_metadata']
L1 = metadata['L1'][0]
L2 = metadata['L2'][0]
rho = metadata['rho'][0]
r = metadata['r'][0]
m = metadata['m'][0]

demography = msprime.Demography.from_tree_sequence(ts)
for i in range(L1):
    for j in range(L2):
        home = i + j * L1 + 1
        demography["p"+str(home)].initial_size=rho
        if i>0:
            demography.set_migration_rate(home, home-1, m/4)
        if i<L1-1:
            demography.set_migration_rate(home, home+1, m/4)
        if j>0:
            demography.set_migration_rate(home, home-L1, m/4)
        if j<L2-1:
            demography.set_migration_rate(home, home+L1, m/4)

print(f"start recapitating at {datetime.now()}")

rts = pyslim.recapitate(
        ts, demography=demography,
        recombination_rate=1e-8,
        random_seed=4
)
print(f"recapitation completed. Adding neutral mutations at {datetime.now()}")

# Sprinkle neutral mutations on top of the recapitated tree sequence
next_id = pyslim.next_slim_mutation_id(rts)
rts = msprime.sim_mutations(
            rts, rate=1e-8, random_seed=7, keep=True,
            model=msprime.SLiMMutationModel(type=0, next_id=next_id)
)
rts.dump("recapitated_ts.trees")
print(f"neutral mutations added. Start running forward-in-time sims post SLiM at {datetime.now()}")

# run forward-in-time sim after sweep for the number of generation provided in command line.
new_time = int(sys.argv[1])
num_demes = L1 * L2
samples = {}
for i in range(num_demes):
    samples["p"+str(i+1)] = rho


new_ts = msprime.sim_ancestry(
              samples=samples,
              demography=demography,
              end_time=new_time,
              sequence_length=rts.sequence_length,
              recombination_rate=r,
              random_seed=9)
new_ts.dump("post_slim_ts.trees")
print(f"msprime sim ancestry done. now adding mutations to the new ts at {datetime.now()}")

new_ts = msprime.sim_mutations(
                 new_ts, rate=1e-8, random_seed=10, keep=True,
                 model=msprime.SLiMMutationModel(type=0)
        )

new_tables = new_ts.tables

new_nodes = np.where(new_tables.nodes.time == new_time)[0]
print(f"There are {len(new_nodes)} nodes from the start of the new simulation. Current time: {datetime.now()}")

slim_nodes = rts.samples(time=0)

print(f"nodes sampled from the original ts at {datetime.now()}")

# randomly give new_nodes IDs in rts
node_map = np.repeat(tskit.NULL, new_tables.nodes.num_rows)
node_map[new_nodes] = np.random.choice(slim_nodes, len(new_nodes), replace=False)

# shift times: in nodes and mutations
# since tree sequences are not mutable, we do this in the tables directly
# also, unmark the nodes at the end of the SLiM simulation as samples
tables = rts.tables
tables.nodes.flags = tables.nodes.flags & ~np.uint32(tskit.NODE_IS_SAMPLE)
tables.nodes.time = tables.nodes.time + new_time
tables.mutations.time = tables.mutations.time + new_time

print(f"merging tables at {datetime.now()}")
# merge the two sets of tables
tables.union(new_tables, node_map,
            add_populations=False,
            check_shared_equality=False)

# get back the tree sequence
full_ts = tables.tree_sequence()
full_ts.dump("full_ts.trees")

fig, ax = plt.subplots(1, 1, figsize=(10, 10))
afs = ts.allele_frequency_spectrum(polarised=True, mode='branch')
ax.loglog(np.arange(1, ts.num_samples+1) / ts.num_samples, afs[1:])
ax.set_xlabel("f")
ax.set_ylabel("n(f)")
fig.savefig("afs.png")
