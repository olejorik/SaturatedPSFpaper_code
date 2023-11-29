#  Here, we use HDR-PSF without any saturation

# load all requirements
include("../src/ini.jl")

exppar = Dict(
    :num => [11,14], #[11,14],
    :av => [1,4,16,256],
    :DRAPlength => [350],
    :expnum => 13, 
    :folder => "HDR/2023-11-17",
    :bgmean => [(s1 =4.607, s2=4.795, s3 = 5.009,s4 = 5.43, s8 = 6.274, s16 = 7.408)],
)

function run_exp(d)
    @unpack num, av, DRAPlength, expnum, folder = d
    @info "num = $num, av = $av"

    # Read files and read/update hardware config
    phasedir = "$folder/$num"
    # experimentname = "data220121-z31-d1"
    experimentname = "12bit-$expnum-$num-$av-$DRAPlength"
    supnumber = "$expnum-$num-$av"

    figdir = plotsdir(experimentname)

    crops = [400]#,128] #, [128,256,512]
    satlevels = [1]#, 2, 3, 4, 8, 16]#, 32]
    ds = 1

    params_input_data = Dict(
        :cropsim => crops,
        :sat => satlevels,
        :psfdir => "$folder/$num",
        :psfnamefunc => (sat->"PSF_av=$(av).tif"),
        :hw => hwConfig("UI3860", 300, 633e-6, 15)
    )
    nruns = 10

    input_dicts = dict_list(params_input_data)

    # this is dictionary of selected keys, for which the plot browser will be made
    # They ar actuall all keys that change value in simulation
    sel_dict = Dict{Symbol, Any}(
        #:config => simconfigs,
        # :phasetype => phasetypes,
        # :cropsim => crops,
        :sat => satlevels,
        # :N => noiselevels,
        # :av => av
    )




    ## Preprocessing parameters
    # Preprocessing; for simulated data we don't wnat to change the sizes in order to be able to use the precomputed ground truth
    satoffset =0
    preprocess_params_all = (
        pad=((0, 0), (0, 0)),
        # crop=((0, 0), (130, 126)),
        finalsize=crops,
        threshold= [0],  # threshold is not needed for HDR image
        sat_level=(255-satoffset) / 255, # slighty less than maximum
        low_level=5/4095,
        downsample=ds,
        offset=[5], 
        bgmean = [d[:bgmean]],
        PRalg= ["GS"],# "GSSatP"], #"GSSatF
        beta = 
        # [ [fill(0.99,30); fill(0.95, 30); fill(0.8,30); fill(0.7,30); fill(0.5, 30); fill(0.25, 30)] ] # put it inside a vector of length 1 to prevent expansio
        #  [  [fill(0.99,100); collect(range(1,0.5,200))] ]
        0.95 
        # [ collect(range(1,0.5,2500)) ]
    )

    preprocess_dict = dict_list(ntuple2dict(preprocess_params_all))

    merge!(sel_dict,Dict(
        :threshold => preprocess_params_all.threshold,
        :downsample => preprocess_params_all.downsample,
        :finalsize => preprocess_params_all.finalsize,
        :PRalg => preprocess_params_all.PRalg,
        :beta => preprocess_params_all.beta,        
        # :low_level => preprocess_params_all.low_level
    ))

    alldicts = dict_list(merge(params_input_data, ntuple2dict(preprocess_params_all)))

    ## Main part
    # plotdicts = copy(alldicts)
    plotdicts = []

    for d in alldicts
        # d = alldicts[1]
        psfname, psfimage = read_exp_data(d, phasedir)
        @info "Read $psfname"

        @unpack hw = d
        
        # here we have HDR image, can just normalise
        psfimage .= psfimage ./ maximum(psfimage)
        
        # display(psfimage)
        display(highlightsaturation(showpsf(psfimage, 5)))
        
        
        input_params = @ntuple psfname psfimage hw
        preprocess_params = dict2ntuple(d)
        psfprocessed, ap, mask = preprocesssingle(input_params, preprocess_params)
        d[:inputpsf] = psfprocessed

        refphase = nothing
        ap2, mask2 = apmask(hw, reverse(size(psfprocessed)), downsampling=d[:downsample],trim=0.1)


        # display(highlightsaturation(showpsf(psfprocessed, 5)))

        ## Calculations a-la spectral method
        y0 = copy(psfprocessed)
        # binarize!(y0, 0.95 * maximum(Float64.(y0)))
        hardthreshold!(y0, 0.5 * maximum(Float64.(y0)))
        y0 .= sqrt.(y0)
        showpsf(abs2.(y0)) # |> display
        x0 = fftshift(ap) .* ifft((fftshift(Float32.(y0))))
        showpsf(abs2.(x0))  # |> display


        ## Phase retrieval part
        @unpack PRalg, beta = d
        checkpoint = DRAPlength
        PRparams = (
            PRalg = PRalg,
            # PRalg = "GS",
            # PRalg = "GSSup",
            # PRalg = "GSSupNormed",
            # PRalg = "GSSatSN",
            # PRalg="GSSat",
            # PRalg = "GSClipped",
            apalg="DRAP",
            iterparams=(iternum=checkpoint, 
            # snapshots=[-2, -1, 0] .+ checkpoint, 
            snapshots = [],
            keephistory=true, tol=1e-10),
            DRAPbeta = beta,
            xinit=x0,
            # inimethod = "spectral",
            inimethod="random",
            # inimethod="data",
            ϕscale=2π,
            appolish=true,
            polishlength=200,
            processblur=false,
            blurparam=BlurPhaseparam(blurwidth=5, maxit=1)
        )

        for nrun in 1:nruns
        # nrun =1 # debug
            
            sol, ap, mask, pr = retrievePhasesingle(input_params, preprocess_params, PRparams)

            # Make plot dictionary
            dplot=copy(d)
            dplot[:nrun] = nrun
            dplot[:solhist]= sol
            dplot[:resampn] = resnumber(sol, ap, resmask = mask2)
            dplot[:resn] = resphinumber(sol, ap, resmask = mask2)
            
            if isnothing(refphase)
                refphase =  ifftshift(angle.(solution(sol)))
            end
            
            restoredpsf = psf(ap, angle.(ifftshift(solution(sol))))
            dplot[:restoredpsf] = restoredpsf

            cm = satpaletteTwoBand(dplot[:sat]; offset = satoffset, nsteps = 256)
            fig = plot_resultsNoGT(sol, ap, mask, dplot[:inputpsf], dplot[:sat]; refphase =refphase, resmask = mask2, colormap = cm)
            Label(fig[0, :], makelabel(psfname, input_params, preprocess_params, PRparams), fontsize=12)
            fig # |> display


            plotname = savename(experimentname, dplot,"pdf", ignores=(:resn, :resampn, :cropsim, :psfdir, :low_level, :offset, :sat_level))
            wsave(plotsdir(figdir,plotname), fig)

            push!(plotdicts,dplot)
        end
        
        # save restored phase to be used as ground truth later
        # First select the phase with the minimum error
        lastdata = filter_dict_list(plotdicts, d)
        distances = map(x->lastdist(x[:solhist]), lastdata)
        bestsol = lastdata[argmin(distances)][:solhist]
        restPhase = angle.(solution(bestsol))  
        wsave(datadir("pro",phasedir,"restPhase_av=$av.tif"), Gray.(restPhase))
    end

    
end

for d in dict_list(exppar)
    # d= dict_list(exppar)[2]
    run_exp(d)
end