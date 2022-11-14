[![DOI](https://zenodo.org/badge/538920643.svg)](https://zenodo.org/badge/latestdoi/538920643)


# Code supporting the ["Saturated PSF paper"](http://dx.doi.org/10.2139/ssrn.4267855) (under review) 


This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> SaturatedPSFpaper_code

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that simulated data are not included in
   in this repo and need to be downloaded independently from [here](https://doi.org/10.5281/zenodo.7096585). 
   Unpack the archive in `data/sims` folder.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

   This will install all necessary packages for you to be able to run the scripts and everything should work out of the box, including correctly finding local paths.

**Warning:** the packages need to be precompiled on the first run, this can take quite a long time.

## TOC
The scripts in the `scripts` folder produce the plots used in Section 5 of the article and the source files for the supplementary "result browsers". If lualatex is installed in the system, scripts also compile the pdf files from the /tex sources. The plots are save in the `plots` folder.

- `SimulatedData_full.jl` uses the simulated data to produce the plots for the whole parameter space, resulting in 3360-page result browser.  This takes about 2.5 hours.
- `sim_noisy_thded_successpolts.jl` and `simT_noisy_thded_successpolts.jl` use the simulated data (for the "low-order" and "turbulent" phases respectively) in the reduced parameter space, with 10 runs with random initial value for each parameter set and generate also the result tables shown in Section 5.
- `z31_F-ds2_multipleruns_browser_plots.jl` does the same on the experimental data.

