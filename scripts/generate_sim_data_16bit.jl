## initialisation
using DrWatson
@quickactivate "Phase"

using PhaseBases
using PhaseRetrieval
using SampledDomains
# using AlternatingProjections
using Images
# using Plots
# using LaTeXStrings
# using FFTW
# using CSV
using HDF5
using DelimitedFiles
# using LinearAlgebra
# using Statistics
using Noise
# default(aspect_ratio=:equal)

# include(srcdir("plotting_functions.jl"))


struct SimConfig
    name::String
    ims::PhaseRetrieval.ImagingSensor
    f::Float64
    λ::Float64
    d::Float64
    q::Int
    roi::CartesianDomain2D
    dualroi::CartesianDomain2D
    ap::Array{Float64,2}
    mask::Array{Float64,2}
    # diversity::Array{Float64,2} # not implemented
end

function SimConfig(name::String, ims::PhaseRetrieval.ImagingSensor, λ::Float64)
    q = PhaseRetrieval.upscaleFactor(ims, λ)
    f = ims.lens.focallength
    d = ims.lens.aperture
    roi = PhaseRetrieval.make_centered_domain2D(ims)
    dualroi = dualDomain(roi, q) * f * λ
    ap, mask = PhaseBases.aperture(dualroi, d)
    # diversity = PhaseRetrieval.constructSHPhaseDiversity(wfs, d, λ)
    return SimConfig(name, ims, f, λ, d, q, roi, dualroi, ap, mask)
end


"""
Generate test configurations, test phases for full-frame sensors
"""
function prepare_sim_configs()
    # First define common parameters for all simulations
    λ = 632.E-6

    # we will run our simulations for several configurations
    # Each config uses the same lens with different aperture and same type of camera
    cam1540 = PhaseRetrieval.CameraChip(pixelsize = 5.2E-3, imagesize = (1280, 1024), bitdepth = 8)
    lensF300A25 = PhaseRetrieval.ImagingLens(300, 25)


    # apertures = (
    dbig = 25.0
    dsmall = 5.0
    dmedium = 12.0
    # )




    # Now define three sensors
    configs = [PhaseRetrieval.ImagingSensor(
        cam = cam1540,
        lens = PhaseRetrieval.diaphragm(lensF300A25, d))
               for d in [dmedium, dbig, dsmall]]


    sims = [SimConfig("config$i", config, λ) for (i, config) in enumerate(configs)]


    return sims

end

function calculate_test_phase(sim::SimConfig, phasetype::String)
    # phasetypes =["loworder", "turbulent", "sparse","ref"]

    if phasetype == "ref"
        phase = float(sim.ap)
    elseif phasetype == "loworder"

        # Parameters for the low-order phase#generate random zernike phase
        # amp = 1pi Nice example of convergence to wrong solution
        amp = 0.75pi
        coefn = [0.0, 0.04487997792421934, 0.061213166549227316, 0.0372359177320852, -0.025807672973317094, -0.10907698627242028, -0.0025831843312457176, -1.2854500898834675, 0.5100406297505503, -0.02063276561330821, -0.14316991684722485, -0.025345357367834, 0.1393506816540738, -0.19999956756592352, 0.13278232836541945, 0.1475377973125301, 0.0015068840385565028, 0.3662886530633074, 0.18953907467066136, 0.33815214098154595, 0.04366857995779775]
        zord = 5
        zer = ZernikeBW(sim.dualroi, sim.d, zord)
        phlow = compose(zer, coefn * amp)
        phase = phlow

    else
        if phasetype == "turbulent"
            phtable = readdlm(datadir("exp_raw/turb512/", "wf200.csv"), ',', Float64)
        elseif phasetype == "sparse"
            phtable = PhaseRetrieval.rescale(readdlm(datadir("exp_raw/sparse/", "TUDelftSeal.csv"), '\t', Float64))
        end

        dims = sim.ims.cam.imagesize
        if size(phtable) != reverse(dims) .* sim.q
            phtable = imresize(phtable, reverse(dims) .* sim.q)
        end

        phase = phtable
    end


    return phase

end



sims = prepare_sim_configs()

phasedir = "saturatedPSF"

simconfigs = [1] #[1, 2, 3]
phasetypes = ["loworder", "turbulent", "sparse", "ref"]

## Ideal PSFs
# Temporary dataframe: configname, phasetype → ideal PSF

params_idealPSF = Dict(
    "config" => simconfigs,
    "phasetype" => phasetypes
)

dicts = dict_list(params_idealPSF)




function prepare_ideal_psfs(d::Dict)
    @unpack config, phasetype = d
    sim = sims[config]
    phase = calculate_test_phase(sim, phasetype)
    idealpsf = psf(PhaseRetrieval.field(sim.ap, phase))
    # wsave(datadir("sims",phasedir,"phases",savename("phase",d,"csv")),phase)
    # wsave(datadir("sims",phasedir,"phases",savename("psf",d,"csv")),idealpsf)
    fulld = copy(d)
    fulld["phase"] = phase
    fulld["psf"] = idealpsf
    wsave(datadir("sims", phasedir, "idealPSFs", savename("phasepsf", d, "hdf5")), fulld)
    return fulld
end


# for d in dicts
#     prepare_ideal_psfs(d)
# end

# can be used later as 
# phasepsf= collect_results(datadir("sims",phasedir,"phases"), valid_filetypes =[".hdf5"])
#  or for a single psf

function load_ideal_psf(config, phasetype)
    d = @dict(config, phasetype)
    filename = datadir("sims", phasedir, "idealPSFs", savename("phasepsf", d, "hdf5"))
    @unpack psf = load(filename)
    return psf
end


## 2a. configname, phasetype, crop → ground truth sampled phase, sampled aperture

crops = [128, 256,512]

function crop_SimConfig(sim::SimConfig, crop)
    ims = PhaseRetrieval.roi(sim.ims, crop)
    name = sim.name * "_crop=$crop"
    λ = sim.λ
    q = PhaseRetrieval.upscaleFactor(ims, λ)
    f = ims.lens.focallength
    d = ims.lens.aperture
    roi = PhaseRetrieval.make_centered_domain2D(ims)
    dualroi = dualDomain(roi, q) * f * λ
    ap, mask = PhaseBases.aperture(dualroi, d)
    # diversity = PhaseRetrieval.constructSHPhaseDiversity(wfs, d, λ)
    return SimConfig(name, ims, f, λ, d, q, roi, dualroi, ap, mask)
end


params_sampledphase = Dict(
    "config" => simconfigs,
    "phasetype" => phasetypes,
    "crop" => crops
)


dicts = dict_list(params_sampledphase)

function prepare_sampled_phase(d::Dict)
    @unpack config, phasetype, crop = d
    sim = crop_SimConfig(sims[config],crop)
    phase = calculate_test_phase(sim, phasetype)
    fulld = copy(d)
    fulld["phase"] = phase
    fulld["ap"] = sim.ap
    wsave(datadir("sims", phasedir, "phases", savename("phaseap", d, "hdf5")), fulld)
    wsave(datadir("sims", phasedir, "phases", savename("ap", @dict(config, crop), "png")), sim.ap)
    return fulld
end

# params_sampledap = Dict(
#     "config" => simconfigs,
#     "crop" => crops
# )

# dicts = dict_list(params_sampledphase)


for d in dicts
    prepare_sampled_phase(d)
end


## 2b. ideal PSF (→=configname, phasetype), crop, noise, saturation → noisy PSF
satlevels = [0.95, 2,4,8,16]
noiselevels = [0,2,4]

params_noisyPSF = Dict(
    "config" => simconfigs,
    "phasetype" => phasetypes,
    "crop" => crops,
    "sat" => satlevels,
    "N" => noiselevels
)

dicts = dict_list(params_noisyPSF)

function croppsf(image, crop)
    if length(crop) == 1
        crop = (crop, crop)
    end
    cpixel = ceil.(Int,(size(image) .+ 1) ./ 2)
    ccrop = ceil.(Int,(crop .+ 1) ./ 2)
    corner = cpixel .- ccrop 
    return image[corner[1]+1 : corner[1]+crop[1], corner[2]+1 : corner[2]+crop[2]]
end

grayformats = Dict(
    8 => Gray{N0f8},
    12 => Gray{N4f12},
    16 => Gray{N0f16}
)

function prepare_noisy_psf(d::Dict, bpp = 8)
    level_number = 2^bpp
    @unpack config, phasetype, crop, sat, N = d
    psfimage = croppsf(load_ideal_psf(config, phasetype),crop)

    psfimageS = psfimage ./ maximum(psfimage) .* sat
    psfimageS = add_gauss!(psfimageS, N/float(level_number), 0.0, clip =true)
    # psfimageS = quantization(psfimageS, level_number)
    fulld = copy(d)
    fulld["bpp"] = bpp
    wsave(datadir("sims", phasedir, "psfs", savename("PSF", fulld, "tif")),  grayformats[bpp].(psfimageS))
    fulld["psf"] = psfimageS
    return fulld
end

## Save the PSFs

for bpp in [8,12,16]
    for d in dicts
        prepare_noisy_psf(d,bpp)
    end
end


