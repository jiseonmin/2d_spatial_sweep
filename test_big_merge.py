import tskit
import pyslim

ts = tskit.load("r=0.trees")
rts = pyslim.recapitate(ts, ancestral_Ne=10000, recombination_rate=0, random_seed=1)
