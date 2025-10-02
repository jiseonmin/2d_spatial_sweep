import tskit
import pyslim
import msprime
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import sys

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

# Keep track of all populations at each level
current_level_pops = [f"p{i+1}" for i in range(total_pops)]
level = 0
time_offset = recap_time

# Merge populations hierarchically until we have one final ancestral population
while len(current_level_pops) > 1:
    next_level_pops = []
    num_groups = (len(current_level_pops) + 98) // 99  # Ceiling division

    for i in range(num_groups):
        # Get up to 99 populations for this group
        start_idx = i * 99
        end_idx = min((i + 1) * 99, len(current_level_pops))
        derived = current_level_pops[start_idx:end_idx]

        # Determine ancestral population name
        if num_groups == 1:
            # This is the final merge into "ancestral"
            ancestral_name = "ancestral"
        else:
            # Create intermediate ancestral population
            ancestral_name = f"p{i}_ancestral_level{level}"
            next_level_pops.append(ancestral_name)
            demography.add_population(
                name=ancestral_name,
                initial_size=len(derived) * rho
            )

        # Set initial sizes for derived populations
        for pop_name in derived:
            if pop_name in demography:
                demography[pop_name].initial_size = rho

        # Add the split event
        split_time = get_next_time(time_offset, level)
        demography.add_population_split(
            split_time,
            derived=derived,
            ancestral=ancestral_name
        )

    current_level_pops = next_level_pops
    level += 1
    time_offset = get_next_time(time_offset, level)

print(f"Demography setup complete with {level} levels of merging at {datetime.now()}")

#for i in range(L1):
#    for j in range(L2):
#        home = i + j * L1 + 1
#        demography["p"+str(home)].initial_size=rho

       # if i>0:
        #    demography.set_migration_rate(home, home-1, m/4)
        #if i<L1-1:
        #    demography.set_migration_rate(home, home+1, m/4)
        #if j>0:
        #    demography.set_migration_rate(home, home-L1, m/4)
        #if j<L2-1:
        #    demography.set_migration_rate(home, home+L1, m/4)

#merge all temporary rows of subpopulation into a single population. This has to happen after the temporary merges occur, so we use nextafter again.

#demography.add_population_split(
#        np.nextafter(np.nextafter(recap_time, 2 * recap_time), 2 * recap_time),
#        derived=ancestral_names,
#        ancestral="ancestral")

print(f"start recapitating at {datetime.now()}")

rts = pyslim.recapitate(
        ts, demography=demography,
        recombination_rate=r,
        random_seed=4
)
print(f"recapitation completed. Adding neutral mutations at {datetime.now()}")

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
afs = ts.allele_frequency_spectrum(polarised=True, mode='branch')
ax.loglog(np.arange(1, ts.num_samples+1) / ts.num_samples, afs[1:])
ax.set_xlabel("f")
ax.set_ylabel("n(f)")
fig.savefig(f"afs_{ts_name}.png")
