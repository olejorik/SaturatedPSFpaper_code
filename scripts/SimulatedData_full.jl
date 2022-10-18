# test of the method on simulated data and build a PDF browser of the results

# load all requirements
include("../src/ini.jl")

# Read files and read/update hardware config
phasedir = "saturatedPSF"
experimentname = "simdataconf1"
supnumber = "2"

figdir = plotsdir(experimentname)

simconfigs = [1] 
phasetypes = ["loworder", "turbulent"]
crops = [128,256] 
satlevels = [0.95, 2, 4, 8, 16]
noiselevels = [0, 2, 4]

params_input_data = Dict(
    :config => simconfigs,
    :phasetype => phasetypes,
    :cropsim => crops,
    :sat => satlevels,
    :N => noiselevels,
)

input_dicts = dict_list(params_input_data)

# this is dictionary of selected keys, for which the plot browser will be made
# They ar actuall all keys that change value in simulation
sel_dict = Dict(
    #:config => simconfigs,
    :phasetype => phasetypes,
    :cropsim => crops,
    :sat => satlevels,
    :N => noiselevels,
)

hw = hwConfig("UI1540", 300, 632e-6, 12)


## Preprocessing parameters
# Preprocessing; for simulated data we don't want to change the sizes in order to be able to use the precomputed ground truth
preprocess_params_all = (
    # pad=((0, 0), (0, 0)),
    # crop=((0, 0), (130, 126)),
    # finalsize=256,
    threshold=[0, 2, 4, 8],
    sat_level=255 / 255,
    low_level=0 / 255,
    downsample=1,
    offset=0, 
    PRalg= ["GSSatF", "GSSatP"],
    beta = [0.0, 0.125, 0.25, 0.5, 0.75, 0.875, 1]
)

preprocess_dict = dict_list(ntuple2dict(preprocess_params_all))

merge!(sel_dict,Dict(
    :threshold => preprocess_params_all.threshold,
    :PRalg => preprocess_params_all.PRalg,
    :beta => preprocess_params_all.beta
))

alldicts = dict_list(merge(params_input_data, ntuple2dict(preprocess_params_all)))

## Main part

for d in alldicts
    psfname, psfimage, ap, gt = read_sim_data(d)
    mask = ap2mask(ap)

    if hw.cam.channelbitdepth == 16
        if hw.cam.bitdepth == 12
            psfimages = reinterpret.(Gray{N4f12}, psfimages)
        elseif hw.cam.bitdepth == 10
            psfimages = reinterpret.(Gray{N6f10}, psfimages)
        elseif hw.cam.bitdepth == 14
            psfimages = reinterpret.(Gray{N2f14}, psfimages)
        end
    end

    # display(psfimage)
    # display(highlightsaturation(showpsf(psfimage, 5)))


    input_params = @ntuple psfname psfimage hw
    preprocess_params = dict2ntuple(d)
    psfprocessed, _ = preprocesssingle(input_params, preprocess_params)

    # display(highlightsaturation(showpsf(psfprocessed, 5)))

    ## Calculations a-la spectral method
    y0 = copy(psfprocessed)
    # binarize!(y0, 0.95 * maximum(Float64.(y0)))
    hardthreshold!(y0, 0.95 * maximum(Float64.(y0)))
    x0 = fftshift(ap) .* ifft((fftshift(Float32.(y0))))
    # showpsf(abs2.(x0))


    ## Phase retrieval part
    @unpack PRalg, beta = d
    checkpoint = 250
    PRparams = (
        PRalg = PRalg,
        apalg="DRAP",
        iterparams=(iternum=checkpoint, 
        # snapshots=[-2, -1, 0] .+ checkpoint, 
        snapshots = [],
        keephistory=true, tol=1e-17),
        DRAPbeta = beta,
        # inimethod="data",
        xinit=x0,
        # inimethod = "spectral",
        inimethod="random",
        ϕscale=0.0π,
        appolish=true,
        polishlength=100,
        processblur=false,
        blurparam=BlurPhaseparam(blurwidth=5, maxit=1)
    )

    gtfun(x) = calculate_phi_err(x, gt, mask)[2]
    sol, ap, mask, pr = retrievePhasesingle(input_params, preprocess_params, PRparams, gtfun=gtfun)

    fig = plot_results(sol, ap, mask, gt, psfprocessed, d[:sat])
    Label(fig[0, :], makelabel(psfname, input_params, preprocess_params, PRparams), textsize=12)

    plotname = savename(experimentname, d,"pdf", ignores=(:downsample, :low_level, :offset, :sat_level))
    wsave(plotsdir(figdir,plotname), fig)
end

## Build plot browser

frametemplate(plotname, filename, navigator) = """
\\begin{frame}[plain]
    \\begin{columns}
        \\column{.65\\textwidth}  			
            \\hypertarget{$plotname}{\\includegraphics{$filename}}
        \\column{.35\\textwidth}  	
            $navigator
    \\end{columns}
\\end{frame}

"""


decoratehiglight(s) = "\\colorbox{SpringGreen!20}{$s}"
decorate(s) = "\\colorbox{Gray!20}{$s}"

openstring = open(plotsdir("browser_template_begin.txt")) do file
    read(file, String)
end

openstring = replace(openstring, "experimentname"=>experimentname)
openstring = replace(openstring, "supnumber"=>"$supnumber")

open(plotsdir("satPSFarticle","browser_$experimentname.tex"), "w") do f
    write(f, openstring)


    plotname(d) = savename(experimentname, d,"pdf", ignores=(:downsample, :low_level, :offset, :sat_level))

    for d in alldicts
        pname = plotname(d)
        filename = pname
        navigator = ""
        for k in keys(sel_dict)
            navigator *= decorate(k) * "\\\\ \n"
            for val in sel_dict[k]
                target = plotname(neighbourdict(d, (k,val)))
                if d[k] == val
                    navigator *= "\\hyperlink{$target}{$(decoratehiglight(val))}"
                else
                    navigator *= "\\hyperlink{$target}{$(decorate(val))}"
                end
            end
            navigator *= "\\\\[0.6\\baselineskip] \n"
        end
        write(f, frametemplate(pname,filename, navigator))
    end

    write(f, "\\end{document}")


end
