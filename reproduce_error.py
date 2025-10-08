import msprime
L1=10
L2=10
rho=10

demography = msprime.Demography()

for k in range(L1 * L2 + 1):
    demography.add_population(name=f"p{k}", initial_size=rho)
demography.add_population_split(
                derived=[f"p{i}" for i in range(L1 *L2)],
                ancestral=f"p{L1*L2}",
                time=10,
        )

ts = msprime.sim_ancestry(
    samples=[
        msprime.SampleSet(1, population=1),
        msprime.SampleSet(2, population=0)],
    demography=demography,
    end_time=0)
