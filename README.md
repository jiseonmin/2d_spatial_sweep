SLiM simulation generalizing [Min et al.](https://academic.oup.com/genetics/article/222/3/iyac139/6696215)

- Set L1 or L2 to 1 to make it into 1D.
- Run simulation on terminal not SLiMgui if you want to run a bigger simulation. For instance, you would type `slim -d L1=500 -d L2=500 2d_spatial_sweep.slim` to run the simulation overwriting the default width and height of deme grid, increasing each of them by 2 orders of magnitude.
- You might need to increase the runtime (currently 10k gen) depending on how large the population is to make sure the final frequency is close to or equal to 1.
- Use [pyslim](https://tskit.dev/pyslim/docs/stable/tutorial.html) to recapitate and add neutral mutations, and [tskit](https://tskit.dev/tskit/docs/stable/introduction.html) for further analyses.

  
