import tskit
import pyslim
import msprime
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import sys
import gc

ts_name = str(sys.argv[1])
post_sweep_time = int(sys.argv[2])

ts = tskit.load(f"{ts_name}.trees")
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
root_times = set([ts.node(n).time for t in ts.trees() for n in t.roots])
recap_time = root_times.pop()

total_pops = L1 * L2

# Add the final ancestral population
demography.add_population(
    name="ancestral",
    description="ancestral population simulated by msprime",
    initial_size=total_pops * rho
)

# Helper function to generate unique time offsets
def get_next_time(current_time, iteration):
    return np.nextafter(current_time, (iteration + 2) * current_time)

# Pre-set all initial population sizes (much faster than doing it in the loop)
print(f"Setting initial population sizes at {datetime.now()}")
for i in range(total_pops):
    pop_name = f"p{i+1}"
    if pop_name in demography:
        demography[pop_name].initial_size = rho

# Build all demographic events in a batch
print(f"Building demographic events at {datetime.now()}")
level = 0
time_offset = recap_time
num_current_pops = total_pops

# Collect all populations and splits to add at once
intermediate_pops = []
all_splits = []

while num_current_pops > 1:
    num_groups = (num_current_pops + 98) // 99

    for group_idx in range(num_groups):
        start_idx = group_idx * 99
        end_idx = min((group_idx + 1) * 99, num_current_pops)

        # Generate names for this group
        if level == 0:
            derived = [f"p{i+1}" for i in range(start_idx, end_idx)]
        else:
            derived = [f"p{i}_ancestral_level{level-1}" for i in range(start_idx, end_idx)]

        # Determine ancestral name
        if num_groups == 1:
            ancestral_name = "ancestral"
        else:
            ancestral_name = f"p{group_idx}_ancestral_level{level}"
            # Store intermediate population info to add later
            intermediate_pops.append((ancestral_name, len(derived) * rho))

        # Store split info
        split_time = get_next_time(time_offset, level * 1000 + group_idx)
        all_splits.append((split_time, derived, ancestral_name))

        del derived

    num_current_pops = num_groups
    level += 1
    time_offset = get_next_time(time_offset, level * 1000)

print(f"Adding {len(intermediate_pops)} intermediate populations at {datetime.now()}")
# Add all intermediate populations at once
for pop_name, initial_size in intermediate_pops:
    demography.add_population(name=pop_name, initial_size=initial_size)

print(f"Adding {len(all_splits)} population splits at {datetime.now()}")
# Add all splits at once
for split_time, derived, ancestral in all_splits:
    demography.add_population_split(split_time, derived=derived, ancestral=ancestral)

# Clean up
del intermediate_pops, all_splits
gc.collect()

print(f"Demography setup complete with {level} levels of merging at {datetime.now()}")


rts = pyslim.recapitate(
        ts, demography=demography,
        recombination_rate=r,
        random_seed=4
)
print(f"recapitation completed. Adding neutral mutations at {datetime.now()}")

# delete ts and demography to save working memory
del ts, demography
gc.collect()

# Sprinkle neutral mutations on top of the recapitated tree sequence
next_id = pyslim.next_slim_mutation_id(rts)
rts = msprime.sim_mutations(
            rts, rate=1e-8, random_seed=7, keep=True,
            model=msprime.SLiMMutationModel(type=0, next_id=next_id)
)
rts.dump(f"recapitated_{ts_name}.trees")
print(f"neutral mutations added. Start running forward-in-time sims post SLiM at {datetime.now()}")

# run forward-in-time sim after sweep for the number of generation provided in command line.

well_mixed_model = msprime.Demography()
well_mixed_model.add_population(initial_size=L1*L2*rho, name='p_all')
new_ts = msprime.sim_ancestry(
              samples={'p_all' : L1*L2*rho},
              demography=well_mixed_model,
              end_time=post_sweep_time,
              sequence_length=rts.sequence_length,
              recombination_rate=r,
              random_seed=9)

new_ts = msprime.sim_mutations(
                 new_ts, rate=1e-8, random_seed=10, keep=True,
                 model=msprime.SLiMMutationModel(type=0)
        )

new_ts.dump(f"post_slim_{ts_name}.trees")
print(f"msprime sim ancestry done. now adding mutations to the new ts at {datetime.now()}")
new_tables = new_ts.tables

del new_ts
gc.collect()

new_nodes = np.where(new_tables.nodes.time == post_sweep_time)[0]
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

del rts
gc.collect()

tables.nodes.flags = tables.nodes.flags & ~np.uint32(tskit.NODE_IS_SAMPLE)
tables.nodes.time = tables.nodes.time + post_sweep_time
tables.mutations.time = tables.mutations.time + post_sweep_time

print(f"merging tables at {datetime.now()}")
# merge the two sets of tables
tables.union(new_tables, node_map,
            add_populations=False,
            check_shared_equality=False)

# get back the tree sequence
full_ts = tables.tree_sequence()
full_ts.dump(f"full_{ts_name}.trees")


fig, ax = plt.subplots(1, 1, figsize=(10, 10))
afs = full_ts.allele_frequency_spectrum(polarised=True, mode='branch')
ax.loglog(np.arange(1, full_ts.num_samples+1) / full_ts.num_samples, afs[1:])
ax.set_xlabel("f")
ax.set_ylabel("n(f)")
fig.savefig(f"afs_{ts_name}.png")
