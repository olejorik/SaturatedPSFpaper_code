#  A template to process data for GSSat article

# load all requirements
include("../src/ini.jl")

# Read files and read/update hardware config
phasedir = "saturatedPSF"
experimentname = "simT"
supnumber = "4"

figdir = plotsdir(experimentname)

simconfigs = [1] #, 2, 3]
phasetypes = [ "turbulent"]#, "sparse"] # , "ref"] # we don't need reference phase
crops = [128]#,256] #, [128,256,512]
satlevels = [0.95, 2, 4, 8, 16]
noiselevels = [0, 2, 4]

params_input_data = Dict(
    :config => simconfigs,
    :phasetype => phasetypes,
    :cropsim => crops,
    :sat => satlevels,
    :N => noiselevels,
    :hw => hwConfig("UI1540", 300, 632e-6, 12),
)
nruns = 10

input_dicts = dict_list(params_input_data)

# this is dictionary of selected keys, for which the plot browser will be made
# They ar actuall all keys that change value in simulation
sel_dict = Dict{Symbol,Any}(
    #:config => simconfigs,
    # :phasetype => phasetypes,
    # :cropsim => crops,
    :sat => satlevels,
    :N => noiselevels,
)


## Preprocessing parameters
# Preprocessing; for simulated data we don't wnat to change the sizes in order to be able to use the precomputed ground truth
const satoffset = 0
preprocess_params_all = (
    # pad=((0, 0), (0, 0)),
    # crop=((0, 0), (130, 126)),
    finalsize=128,
    threshold=[0, 2, 4, 8, 16],
    sat_level=(255 - satoffset) / 255,
    low_level=0 / 255,
    downsample=1,
    offset=0,
    PRalg=["GSSatF"],# "GSSatP"],
    beta=[0.5],
)

preprocess_dict = dict_list(ntuple2dict(preprocess_params_all))

merge!(
    sel_dict,
    Dict(
        :threshold => preprocess_params_all.threshold,
        # :downsample => preprocess_params_all.downsample,
        # :finalsize => preprocess_params_all.finalsize,
        # :PRalg => preprocess_params_all.PRalg,
        # :beta => preprocess_params_all.beta,
    ),
)

alldicts = dict_list(merge(params_input_data, ntuple2dict(preprocess_params_all)))

## Main part
plotdicts = []

for d in alldicts
    # d = alldicts[22]
    psfname, psfimage, ap, gt = read_sim_data(d)
    @unpack hw = d
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
    d[:inputpsf] = psfprocessed

    refphase = gt
    ap2, mask2 = apmask(
        hw, reverse(size(psfprocessed)); downsampling=d[:downsample], trim=0.1
    )
    # display(highlightsaturation(showpsf(psfprocessed, 5)))

    ## Calculations a-la spectral method
    y0 = copy(psfprocessed)
    # binarize!(y0, 0.95 * maximum(Float64.(y0)))
    hardthreshold!(y0, 0.95 * maximum(Float64.(y0)))
    x0 = fftshift(ap) .* ifft((fftshift(Float32.(y0))))
    # showpsf(abs2.(x0))

    ## Phase retrieval part
    @unpack PRalg, beta = d
    checkpoint = 150
    PRparams = (
        PRalg=PRalg,
        # PRalg = "GS",
        # PRalg = "GSSup",
        # PRalg = "GSSupNormed",
        # PRalg = "GSSatSN",
        # PRalg="GSSat",
        # PRalg = "GSClipped",
        apalg="DRAP",
        iterparams=(
            iternum=checkpoint,
            # snapshots=[-2, -1, 0] .+ checkpoint, 
            snapshots=[],
            keephistory=true,
            tol=1e-17,
        ),
        DRAPbeta=beta,
        # inimethod="data",
        xinit=x0,
        # inimethod = "spectral",
        inimethod="random",
        ϕscale=2.0π,
        appolish=true,
        polishlength=50,
        processblur=false,
        blurparam=BlurPhaseparam(; blurwidth=5, maxit=1),
    )

    gtfun(x) = calculate_phi_err(x, gt, mask)[2]
    for nrun in 1:nruns
        # nrun=1
        dplot = copy(d)

        dplot[:nrun] = nrun

        sol, ap, mask, pr = retrievePhasesingle(
            input_params, preprocess_params, PRparams; gtfun=gtfun
        )

        dplot[:solhist] = sol

        dplot[:resampn] = resnumber(sol, ap; resmask=mask2)
        dplot[:resn] = resphinumber(sol, ap; resmask=mask2)

        if isnothing(refphase)
            refphase = ifftshift(angle.(solution(sol)))
        end

        restoredpsf = psf(ap, angle.(ifftshift(solution(sol))))
        dplot[:restoredpsf] = restoredpsf


        # fig = plot_results(sol, ap, mask, gt, psfprocessed, d[:sat])
        # Label(
        #     fig[0, :],
        #     makelabel(psfname, input_params, preprocess_params, PRparams);
        #     textsize=12,
        # )

        # plotname = savename(
        #     experimentname,
        #     dplot,
        #     "pdf";
        #     ignores=(:resn, :resampn, :cropsim, :psfdir, :low_level, :offset, :sat_level),
        # )
        cm = satpaletteTwoBand(dplot[:sat]; offset=satoffset, nsteps=256)
        fig = plot_resultsNoGT(
            sol,
            ap,
            mask,
            dplot[:inputpsf],
            dplot[:sat];
            refphase=refphase,
            resmask=mask2,
            colormap=cm
        )
        Label(
            fig[0, :],
            makelabel(psfname, input_params, preprocess_params, PRparams);
            textsize=12
        )

        plotname = savename(
            experimentname,
            dplot,
            "pdf";
            ignores=(:resn, :resampn, :cropsim, :psfdir, :low_level, :offset, :sat_level)
        )

        wsave(plotsdir(figdir, plotname), fig)
        push!(plotdicts, dplot)
    end
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

openstring = open(srcdir("browser_template_begin.txt")) do file
    read(file, String)
end

openstring = replace(openstring, "experimentname" => experimentname)
openstring = replace(openstring, "supnumber" => "$supnumber")

fname = plotsdir("browser_$experimentname.tex")
open(plotsdir(fname), "w") do f
    write(f, openstring)

    function plotname(d)
        return savename(
            experimentname, d, "pdf";
            ignores=(:resn, :resampn, :cropsim, :psfdir, :low_level, :offset, :sat_level)
        )
    end

    merge!(sel_dict, Dict(:nrun => collect(1:nruns)))

    for d in plotdicts
        pname = plotname(d)
        filename = pname
        navigator = ""
        for k in keys(sel_dict)
            navigator *= decorate(k) * "\\\\ \n"
            for val in sel_dict[k]
                target = plotname(neighbourdict(d, (k, val)))
                if d[k] == val
                    navigator *= "\\hyperlink{$target}{$(decoratehiglight(val))}"
                else
                    navigator *= "\\hyperlink{$target}{$(decorate(val))}"
                end
            end
            navigator *= "\\\\[0.6\\baselineskip] \n"
        end
        write(f, frametemplate(pname, filename, navigator))
    end

    write(f, "\\end{document}")
end

dn = dirname(fname)
bn = basename(fname)
curd = pwd()
cd(dn)
run(`lualatex.exe -synctex=1 -interaction=nonstopmode $bn`)
# run(`latexmk -c`)
cd(curd)

## Making table of input psfs for fixed noise level, changing sat and threshold
sel_dicts = dict_list(Dict(
    #     :phasetype => "loworder",
    #     :cropsim => 128,
    #     # :PRalg => preprocess_params_all[:PRalg],
    #     :beta => preprocess_params_all[:beta],
    :N => params_input_data[:N],
))
# sel_dicts = dict_list(delete!(copy(sel_dict),:nrun))

# saving number of succesive runs in a dictionary of tables
errdict = Dict()
distdict = Dict()
for sel_dict in sel_dicts
    errtable = zeros(length(preprocess_params_all[:threshold]), length(params_input_data[:sat]))
    disttable = zeros(
        length(preprocess_params_all[:threshold]), length(params_input_data[:sat])
    )

    axlist = []
    fig = Figure(; resolution=(1000, 1000))
    sel_plots = filter_dict_list(plotdicts, sel_dict)
    for (j, sat) in enumerate(params_input_data[:sat])
        for (i, threshold) in enumerate(preprocess_params_all[:threshold])
            dallruns = filter_dict_list(sel_plots, Dict(:sat => sat, :threshold => threshold))
            d = dallruns[1]
            ax, img = image(
                fig[i, j], rotr90(showpsf(d[:inputpsf], 5)); axis=(aspect=DataAspect(),)
            )
            hidedecorations!(ax)
            if j == 1
                Box(fig[i, 0]; color=:gray90)
                Label(
                    fig[i, 0],
                    L"t=%$threshold";
                    tellheight=false,
                    rotation=pi / 2,
                    padding=(5, 5, 5, 5)
                )
            end
            errtable[i, j] = length(filter(x -> x == 0, getindex.(dallruns, :resn))) #number of runs with zero residues
            disttable[i, j] = minimum(last.(map(x -> distgthist(x[:solhist]), dallruns))) #smallest distance
        end
        Box(fig[1, j, Top()]; color=:gray90)
        Label(fig[1, j, Top()], L"s=%$sat"; tellwidth=false, rotation=0, padding=(5, 5, 5, 5))
        colsize!(fig.layout, j, Aspect(1, 1.0))
    end
    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)
    N = sel_dict[:N]
    Label(fig[0, :], "Input PSF, noise level N = $N", textsize=30)
    resize_to_layout!(fig)
    display(fig)

    wsave(plotsdir(experimentname, "inputPSFtable-$N.pdf"), fig)

    errdict[N] = errtable
    distdict[N] = disttable
end

## Making similar table for the output PSF

for sel_dict in sel_dicts
    axlist = []
    fig = Figure(; resolution=(1000, 1000))
    sel_plots = filter_dict_list(plotdicts, sel_dict)
    for (j, sat) in enumerate(params_input_data[:sat])
        for (i, threshold) in enumerate(preprocess_params_all[:threshold])
            dallruns = filter_dict_list(sel_plots, Dict(:sat => sat, :threshold => threshold))
            d = sort(dallruns; by=x -> getindex(x, :resn))[1]
            d = sort(dallruns; by=x -> last(disthist(x[:solhist])))[1]
            ax, img = image(
                fig[i, j], rotr90(showpsf(d[:restoredpsf], 5)); axis=(aspect=DataAspect(),)
            )
            hidedecorations!(ax)
            if j == 1
                Box(fig[i, 0]; color=:gray90)
                Label(
                    fig[i, 0],
                    L"t=%$threshold";
                    tellheight=false,
                    rotation=pi / 2,
                    padding=(5, 5, 5, 5)
                )
            end
        end
        Box(fig[1, j, Top()]; color=:gray90)
        Label(fig[1, j, Top()], L"s=%$sat"; tellwidth=false, rotation=0, padding=(5, 5, 5, 5))
        colsize!(fig.layout, j, Aspect(1, 1.0))
    end
    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)
    N = sel_dict[:N]
    Label(fig[0, :], "Restored PSF, noise level N = $N", textsize=30)
    resize_to_layout!(fig)
    display(fig)

    wsave(plotsdir(experimentname, "restoredPSFtable-$N.pdf"), fig)
end

## Making similar table for the ouptut phase

refphase = angle.(ifftshift(solution(plotdicts[121][:solhist]))) # fix it later to chose automatically
let (mask, mask2) = (nothing, nothing)
    for sel_dict in sel_dicts
        axlist = []
        fig = Figure(; resolution=(1000, 1000))
        sel_plots = filter_dict_list(plotdicts, sel_dict)
        for (j, sat) in enumerate(params_input_data[:sat])
            for (i, threshold) in enumerate(preprocess_params_all[:threshold])
                dallruns = filter_dict_list(
                    sel_plots, Dict(:sat => sat, :threshold => threshold)
                )
                d = sort(dallruns; by=x -> getindex(x, :resn))[1]
                d = sort(dallruns; by=x -> last(disthist(x[:solhist])))[1]

                if isnothing(mask)
                    ap, mask = apmask(
                        params_input_data[:hw], d[:finalsize]; downsampling=d[:downsample]
                    )
                    _, mask2 = apmask(
                        params_input_data[:hw],
                        d[:finalsize];
                        downsampling=d[:downsample],
                        trim=2
                    )
                    mask2 .*= (-0.1 .< refphase .< 1)
                end
                phi =
                    select_twin(angle.(ifftshift(solution(d[:solhist]))), refphase, mask2) .*
                    mask
                # phidiff, phierr = calculate_phi_err(solution(sol), refphase, mask)
                ax, img = showphasetight(
                    bboxview(removepiston(removetiptilt(phi, mask2), mask2), mask, 2),
                    fig[i, j];
                    hidedec=false
                )

                hidedecorations!(ax)
                if j == 1
                    Box(fig[i, 0]; color=:gray90)
                    Label(
                        fig[i, 0],
                        L"t=%$threshold";
                        tellheight=false,
                        rotation=pi / 2,
                        padding=(5, 5, 5, 5)
                    )
                end
            end
            Box(fig[1, j, Top()]; color=:gray90)
            Label(
                fig[1, j, Top()], L"s=%$sat"; tellwidth=false, rotation=0, padding=(5, 5, 5, 5)
            )
            colsize!(fig.layout, j, Aspect(1, 1.0))
        end
        colgap!(fig.layout, 5)
        rowgap!(fig.layout, 5)
        N = sel_dict[:N]
        Label(fig[0, :], "Restored phase, noise level N = $N", textsize=30)
        resize_to_layout!(fig)
        display(fig)

        wsave(plotsdir(experimentname, "restoredphasetable-$N.pdf"), fig)
    end
end

## Making similar table for the gt error

for sel_dict in sel_dicts
    N = sel_dict[:N]

    # vals= errdict[N]
    # vals = Int.(errtable)
    vals = distdict[N]
    fig = Figure()
    cmin, cmax = extrema(vals)
    # cmin, cmax = (0, 0.75)
    # chigh = :red
    cscheme = ColorSchemes.isoluminant_cgo_70_c39_n256
    cscheme = reverse(ColorSchemes.diverging_rainbow_bgymr_45_85_c67_n256)
    cscheme = reverse(ColorSchemes.RdYlGn_10)

    for (j, sat) in enumerate(params_input_data[:sat])
        for (i, threshold) in enumerate(preprocess_params_all[:threshold])
            v = vals[i, j]
            c = v > cmax ? chigh : get.(Ref(cscheme), (v - cmin) / (cmax - cmin))
            Box(fig[i, j]; color=c, width=150, height=150)
            Label(fig[i, j], "$(round(v, digits=3))"; textsize=30)
            if j == 1
                Box(fig[i, 1, Left()]; color=:gray90, padding=(5, 5, 5, 5))
                Label(
                    fig[i, 1, Left()],
                    L"t=%$threshold";
                    tellheight=false,
                    rotation=pi / 2,
                    padding=(5, 5, 5, 5)
                )
            end
            if i == 1
                Box(fig[1, j, Top()]; color=:gray90)
                Label(
                    fig[1, j, Top()],
                    L"s=%$sat";
                    tellwidth=false,
                    rotation=0,
                    padding=(5, 5, 5, 5)
                )
            end
        end
    end
    Label(fig[-1, :], "Phase error, noise level N = $N", textsize=30)
    Colorbar(
        fig[end+1, :];
        limits=(cmin, cmax),
        colormap=cscheme,
        vertical=false,
        flipaxis=false
        # highclip = :red,
        # width = 150, height = 150
    )
    rowgap!(fig.layout, 5)
    colgap!(fig.layout, 5)
    resize_to_layout!(fig)
    display(fig)
    wsave(plotsdir(experimentname, "successruns-$N.pdf"), fig)
end

## Making similar table for the number of residuals

for sel_dict in sel_dicts
    N = sel_dict[:N]

    vals= Int.(errdict[N])
    # vals = Int.(errtable)
    # vals= disttable
    fig = Figure()
    cmin, cmax = extrema(vals)
    # cmin, cmax = (0, 0.75)
    # chigh = :red
    cscheme = ColorSchemes.isoluminant_cgo_70_c39_n256
    cscheme = reverse(ColorSchemes.diverging_rainbow_bgymr_45_85_c67_n256)
    cscheme = ColorSchemes.RdYlGn_10

    for (j, sat) in enumerate(params_input_data[:sat])
        for (i, threshold) in enumerate(preprocess_params_all[:threshold])
            v = vals[i, j]
            c = v > cmax ? chigh : get.(Ref(cscheme), (v - cmin) / (cmax - cmin))
            Box(fig[i, j]; color=c, width=150, height=150)
            # Label(fig[i,j,], "$(round(v, digits=3))", textsize =30) 
            Label(fig[i, j], "$v"; textsize=30)
            if j == 1
                Box(fig[i, 1, Left()]; color=:gray90, padding=(5, 5, 5, 5))
                Label(
                    fig[i, 1, Left()],
                    L"t=%$threshold";
                    tellheight=false,
                    rotation=pi / 2,
                    padding=(5, 5, 5, 5)
                )
            end
            if i == 1
                Box(fig[1, j, Top()]; color=:gray90)
                Label(
                    fig[1, j, Top()],
                    L"s=%$sat";
                    tellwidth=false,
                    rotation=0,
                    padding=(5, 5, 5, 5)
                )
            end
        end
    end
    Label(fig[-1, :], "Number of successful runs, noise level N = $N", textsize=30)
    Colorbar(
        fig[end+1, :];
        limits=(cmin, cmax),
        colormap=cscheme,
        vertical=false,
        flipaxis=false
        # highclip = :red,
        # width = 150, height = 150
    )
    rowgap!(fig.layout, 5)
    colgap!(fig.layout, 5)
    resize_to_layout!(fig)
    display(fig)
    wsave(plotsdir(experimentname, "ressuccessruns-$N.pdf"), fig)
end