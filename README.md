# SaturatedPSFpaper_code

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

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.
