using DrWatson
@quickactivate "SaturatedPSFpaper_code"

# using PhaseBases
using PhaseRetrieval
using SampledDomains
using AlternatingProjections
using Images
using FFTW
using HDF5
using LaTeXStrings
using ColorSchemes
# Main subfolder whith all different experimental data 
psffolder = datadir("exp") 

include("core.jl")
include("phaseutils.jl")
include("plotting_functionsMakie.jl")