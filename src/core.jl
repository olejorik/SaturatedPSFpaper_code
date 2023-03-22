# Fucntion for read/save and processing

# TODO Refactor to function readfiles, as it has notheing to do with PSFs and move to utils
function readpsfs(psfdir, psfnamemask)
    filenames = readdir(datadir(psffolder, psfdir))
    psfnames = filter(x -> occursin(psfnamemask, x), filenames)
    psfimages = [load(datadir(psffolder, psfdir, psfname)) for psfname in psfnames]

    return psfnames, psfimages
end

Base.@kwdef struct hwConfig
    cam::PhaseRetrieval.CameraChip
    f::Float64
    λ::Float64
    d::Float64
end

hwConfig(s::String, f, λ, d) = hwConfig(PhaseRetrieval.camerasdict[s], f, λ, d)

function preprocess(input_params, preprocess_params)
    preprocess_def = (
        pad=((0, 0), (0, 0)),
        crop=((0, 0), (0, 0)),
        threshold=2,
        sat_level=255,
        downsample=1,
        offset=0,
        finalsize=0,
    )
    preprocess_params = merge(preprocess_def, preprocess_params)

    @unpack psfnames, psfimages, hw = input_params

    processed = [float64.(psfimage) for psfimage in psfimages]

    for psfimage in processed
        psfimage .+= -(preprocess_params.offset / (2^hw.cam.bitdepth - 1))
        hardthreshold!(psfimage, preprocess_params.threshold / (2^hw.cam.bitdepth - 1))
    end

    if preprocess_params.finalsize == 0
        crop = preprocess_params.crop
        processed = [
            psfimage[
                (1 + crop[1][1]):(end - crop[1][2]), (1 + crop[2][1]):(end - crop[2][2])
            ] for psfimage in processed
        ]
        processed = [
            PhaseRetrieval.binning(
                padarray(psfimage, Fill(0, preprocess_params.pad...)),
                preprocess_params.downsample,
            ) for psfimage in processed
        ]
    else
        cropsize = preprocess_params.finalsize * preprocess_params.downsample
        c = round.(Int, centroid(hardthreshold(processed[1],0.2,true)))
        processed = [
            PhaseRetrieval.binning(
                PaddedView(
                    0,
                    psfimage,
                    (
                        (c[1] - cropsize ÷ 2):(c[1] - cropsize ÷ 2 + cropsize - 1),
                        (c[2] - cropsize ÷ 2):(c[2] - cropsize ÷ 2 + cropsize - 1),
                    ),
                ),
                preprocess_params.downsample,
            ) for psfimage in processed
        ]
    end

    crop = reverse(size(processed[1]))
    ap, mask = apmask(hw, crop; downsampling=preprocess_params.downsample)

    return processed, ap, mask
end

function preprocesssingle(input_params, preprocess_params)
    preprocess_def = (
        pad=((0, 0), (0, 0)),
        crop=((0, 0), (0, 0)),
        threshold=2,
        sat_level=255,
        downsample=1,
        offset=0,
        finalsize=0,
    )
    preprocess_params = merge(preprocess_def, preprocess_params)

    @unpack psfname, psfimage, hw = input_params

    processed = Float64.(deepcopy(psfimage))
    # @info typeof(processed)

    processed .-= (preprocess_params.offset / (2^hw.cam.bitdepth))
    hardthreshold!(processed, 0) #subtract offset and make all negative zero
    hardthreshold!(processed, preprocess_params.threshold / (2^hw.cam.bitdepth))

    if preprocess_params.finalsize == 0
        crop = preprocess_params.crop
        processed = processed[
            (1 + crop[1][1]):(end - crop[1][2]), (1 + crop[2][1]):(end - crop[2][2])
        ]
        processed = PhaseRetrieval.binning(
            padarray(processed, Fill(0, preprocess_params.pad...)),
            preprocess_params.downsample,
        )
    else
        cropsize = preprocess_params.finalsize * preprocess_params.downsample
        c = round.(Int, centroid(hardthreshold(processed,0.2,true)))
        processed = PhaseRetrieval.binning(
            PaddedView(
                0,
                processed,
                (
                    (c[1] - cropsize ÷ 2):(c[1] - cropsize ÷ 2 + cropsize - 1),
                    (c[2] - cropsize ÷ 2):(c[2] - cropsize ÷ 2 + cropsize - 1),
                ),
            ),
            preprocess_params.downsample,
        )
    end

    crop = reverse(size(processed))
    ap, mask = apmask(hw, crop; downsampling=preprocess_params.downsample)

    return processed, ap, mask
end
function hardthreshold!(a::AbstractArray, th)
    a[a .<= th] .= 0
    return a
end

function hardthreshold!(a::AbstractArray, th, relative::Bool)
    if relative
        a[a .<= th * maximum(a)] .= 0
        return a
    else
        return hardthreshold!(a,th)
    end
end

hardthreshold(a, args...) = hardthreshold!(copy(a), args...)

function centroid(a::AbstractArray)
    c = [0, 0]
    for i in CartesianIndices(a)
        c += [Tuple(i)...] * gray(a[i])
    end
    if sum(gray, a) == 0
        return size(a) ./ 2
    else
        return c / sum(gray, a)
    end
end

function apmask(hw::hwConfig, crop=missing; downsampling=1, trim=0)
    @unpack cam, f, λ, d = hw
    if !ismissing(crop)
        cam = PhaseRetrieval.roi(cam, crop)
    end
    q = PhaseRetrieval.upscaleFactor(cam.pixelsize * downsampling, f, d, λ)
    roi = make_centered_domain2D(cam.imagesize..., cam.pixelsize * downsampling)
    dualroi = dualDomain(roi, q) * f * λ
    return PhaseRetrieval.aperture(dualroi, d * (1 - trim)) # (ap, mask)
end

function makeprproblems(input_params, preprocess_params, PRparams)
    params_def = (sat_level=1.0, low_level=1 / 255, PRalg="GS")

    params = merge(params_def, input_params, preprocess_params, PRparams)
    @info params.sat_level
    psfprocessed, ap, mask = preprocess(input_params, preprocess_params)

    a = fftshift(sqrt.(Float64.(ap)))
    prs = []
    alg = params.PRalg

    for (psfimage, filename) in zip(psfprocessed, input_params.psfnames)
        @info "Preprocessing $filename"

        A = fftshift(sqrt.(collect(Float64, psfimage))) # need it for Fourier plan
        size(a) == size(A) || error("Array sizes do not match: $(size(a)) and $(size(A))")

        N = sqrt(sum(abs2, A))
        n = sqrt(sum(abs2, a))

        if alg == "GS"
            # normalise A to a
            A = A ./ N .* n
            pr = TwoSetsFP(
                ConstrainedByAmplitude(a), FourierTransformedSet(ConstrainedByShape(A))
            )
        elseif alg == "GSSup"
            A = A ./ N .* n
            pr = TwoSetsFP(
                ConstrainedBySupport(a), FourierTransformedSet(ConstrainedByShape(A))
            )
        elseif alg == "GSSupNormed"
            A = A ./ N .* n
            pr = TwoSetsFP(
                ConstrainedBySupportNormed(a), FourierTransformedSet(ConstrainedByShape(A))
            )
        elseif alg == "GSSatF"
            sat = A .< sqrt(params.sat_level)
            @info "making saturated set"
            @info "extremsa of A $(extrema(A))"
            @info "using sat level $(params.sat_level)"
            A = A ./ N .* n
            A ./= params.sat_level
            pr = TwoSetsFP(
                ConstrainedByAmplitude(a),
                FourierTransformedSet(ConstrainedByShapeSaturated(A, sat)),
            )
        elseif alg == "GSSatP"
            sat = A .< sqrt(params.sat_level)
            @info "making saturated set"
            @info "extremsa of A $(extrema(A))"
            @info "using sat level $(params.sat_level)"
            A = A ./ N .* n
            A ./= params.sat_level
            pr = TwoSetsFP(
                ConstrainedByShape(a),
                FourierTransformedSet(ConstrainedByShapeSaturated(A, sat)),
            )
        elseif alg == "GSClipped"
            A = A ./ N .* n
            @info "making clipped set"
            @info "extremsa of A $(extrema(A))"
            @info "using sat level $(params.sat_level)"
            @info "using low level $(params.low_level)"
            pr = TwoSetsFP(
                ConstrainedByAmplitude(a),
                FourierTransformedSet(
                    ConstrainedByShapeClipped(
                        A,
                        float(sqrt(params.sat_level)) ./ N .* n,
                        sqrt(params.low_level) ./ N .* n,
                    ),
                ),
            )
        elseif alg == "GSSatSN"
            sat = A .< sqrt(params.sat_level)
            A = A ./ N .* n
            @info "using sat level $(params.sat_level)"
            A ./= params.sat_level
            pr = TwoSetsFP(
                ConstrainedBySupportNormed(a),
                FourierTransformedSet(ConstrainedByShapeSaturated(A, sat)),
            )
        else
            error("Method $alg is unknow")
        end

        push!(prs, pr)
    end

    return prs, ap, mask
end

function makeprproblem(input_params, preprocess_params, PRparams)
    params_def = (sat_level=1.0, low_level=1 / 255, PRalg="GS")

    params = merge(params_def, input_params, preprocess_params, PRparams)
    @info params.sat_level
    psfprocessed, ap, mask = preprocesssingle(input_params, preprocess_params)

    a = fftshift(sqrt.(Float64.(ap)))
    prs = []
    alg = params.PRalg

    (psfimage, filename) = (psfprocessed, input_params.psfname)

    @info "Preprocessing $filename"
    # showpsf(psfimage, 5) |> display

    A = fftshift(sqrt.(collect(Float64, psfimage))) # need it for Fourier plan
    size(a) == size(A) || error("Array sizes do not match: $(size(a)) and $(size(A))")

    N = sqrt(sum(abs2, A))
    n = sqrt(sum(abs2, a))

    if alg == "GS"
        # normalise A to a
        A = A ./ N .* n
        # pr = TwoSetsFP(ConstrainedByAmplitude(a), FourierTransformedSet(ConstrainedByAmplitude(A)))
        # quick fix for zer-threshold
        pr = TwoSetsFP(
            ConstrainedByAmplitude(a), FourierTransformedSet(ConstrainedByShape(A))
        )
    elseif alg == "GSSup"
        A = A ./ N .* n
        # pr = TwoSetsFP(ConstrainedBySupport(a), FourierTransformedSet(ConstrainedByAmplitude(A)))
        pr = TwoSetsFP(
            ConstrainedBySupport(a), FourierTransformedSet(ConstrainedByShape(A))
        )
    elseif alg == "GSSupNormed"
        A = A ./ N .* n
        # pr = TwoSetsFP(ConstrainedBySupportNormed(a), FourierTransformedSet(ConstrainedByAmplitude(A)))
        pr = TwoSetsFP(
            ConstrainedBySupportNormed(a), FourierTransformedSet(ConstrainedByShape(A))
        )
    elseif alg == "GSSatF"
        sat = A .< sqrt(params.sat_level)
        @info "making saturated set"
        @info "extremsa of A $(extrema(A))"
        @info "using sat level $(params.sat_level)"
        A = A ./ N .* n
        A ./= params.sat_level
        pr = TwoSetsFP(
            ConstrainedByAmplitude(a),
            FourierTransformedSet(ConstrainedByShapeSaturated(A, sat)),
        )
    elseif alg == "GSSatP"
        sat = A .< sqrt(params.sat_level)
        @info "making saturated set"
        @info "extremsa of A $(extrema(A))"
        @info "using sat level $(params.sat_level)"
        A = A ./ N .* n
        A ./= params.sat_level
        pr = TwoSetsFP(
            ConstrainedByShape(a),
            FourierTransformedSet(ConstrainedByShapeSaturated(A, sat)),
        )
    elseif alg == "GSClipped"
        A = A ./ N .* n
        @info "making clipped set"
        @info "extremsa of A $(extrema(A))"
        @info "using sat level $(params.sat_level)"
        @info "using low level $(params.low_level)"
        pr = TwoSetsFP(
            ConstrainedByAmplitude(a),
            FourierTransformedSet(
                ConstrainedByShapeClipped(
                    A,
                    float(sqrt(params.sat_level)) ./ N .* n,
                    sqrt(params.low_level) ./ N .* n,
                ),
            ),
        )
    elseif alg == "GSSatSN"
        sat = A .< sqrt(params.sat_level)
        A = A ./ N .* n
        @info "using sat level $(params.sat_level)"
        A ./= params.sat_level
        pr = TwoSetsFP(
            ConstrainedBySupportNormed(a),
            FourierTransformedSet(ConstrainedByShapeSaturated(A, sat)),
        )
    else
        error("Method $alg is unknow")
    end

    return pr, ap, mask
end
function makesolvers(PRparams, ap)
    @unpack PRalg,
    apalg, iterparams, appolish, DRAPbeta, polishlength, ϕscale, processblur, blurparam,
    inimethod = PRparams

    # initial phase
    ϕ0 = rand(Float64, size(ap)) * ϕscale # Random phase
    xinitrand = fftshift(PhaseRetrieval.field(ap, ϕ0))
    if inimethod == "random"
        xinit = xinitrand
    elseif inimethod == "spectral"
        xinit = zeros(ComplexF64, size(ap))
        xinit[1, 1] = 1
        xinit = xinit .* xinitrand
    elseif inimethod == "data"
        @unpack xinit = PRparams
        xinit = xinit .* xinitrand
    elseif error("Initialization method $inimethod is not defined")
    end

    iternum = iterparams.iternum
    snapshots = iterparams.snapshots #collect(1:snapshots:(iternum + 1))
    keephistory = iterparams.keephistory
    tol = iterparams.tol

    method = Dict("AP" => APparam, "DR" => DRparam, "DRAP" => DRAPparam)[apalg]
    if apalg == "DRAP"
        solver1 = method(;
            x⁰=xinit,
            maxϵ=tol,
            maxit=iternum + 1,
            keephistory=keephistory,
            snapshots=snapshots,
            β=DRAPbeta,
        )
    else
        solver1 = method(;
            x⁰=xinit,
            maxϵ=tol,
            maxit=iternum + 1,
            keephistory=keephistory,
            snapshots=snapshots,
        )
    end

    if appolish
        solver2 = APparam(; x⁰=xinit, maxϵ=tol, maxit=polishlength, keephistory=keephistory)
        if processblur
            solverblur = blurparam
            return (solver1, solverblur, solver2)
        else
            return (solver1, solver2)
        end
    else
        if processblur
            solverblur = blurparam
            return (solver1, solverblur)
        else
            return (solver1)
        end
    end
end


function retrievePhasesingle(input_params, preprocess_params, PRparams; kwargs...)
    prproblem, ap, mask = makeprproblem(input_params, preprocess_params, PRparams)

    solvers = makesolvers(PRparams, ap)

    @unpack psfname = input_params
    solhist = (@info "\nProcessing $psfname"; solve(prproblem, solvers; kwargs...))

    return solhist, ap, mask, prproblem
end

function retrievePhaseRandomRuns(
    input_params, preprocess_params, PRparams, runs=3, ϕscale=2π; args...
)
    prproblems, ap, mask = makeprproblems(input_params, preprocess_params, PRparams)

    solvers = makesolvers(PRparams, ap)

    @unpack psfnames = input_params
    solutions = [
        begin
            @info "\n\nProcessing $psfname, run #$i"
            ϕ0 = rand(Float64, size(ap)) * ϕscale # Random phase
            xinitrand = fftshift(PhaseRetrieval.field(ap, ϕ0))
            solve(pr, solvers; x⁰=xinitrand, args...)
        end for i in 1:runs, (pr, psfname) in zip(prproblems, psfnames)
    ]

    return solutions, ap, mask, prproblems
end

function satpaletteTwoBand(s; offset=0, nsteps=256)
    cgray = reverse(colormap("Blues", nsteps - 1 - offset))

    addsteps = round(Int, s * nsteps - (nsteps - 1 - offset))
    if addsteps <= 1
        cred = [colorant"red"]
    else
        cred = range(colorant"red"; stop=colorant"white", length=addsteps)
    end
    return ColorScheme(vcat(cgray, cred))
end


function plot_resultsNoGT(
    sol,
    ap,
    mask,
    inputpsf,
    sat;
    refphase=nothing,
    calculategtpsf=false,
    usegtphase=false,
    lowthr=0.5,
    resmask=nothing,
    colormap=nothing,
)
    fig = Figure(; resolution=(1000, 900))
    phi = angle.(ifftshift(solution(sol))) .* mask
    restoredpsf = psf(ap, angle.(ifftshift(solution(sol))))
    if isnothing(refphase)
        refphase = zeros(size(phi))
    end
    if calculategtpsf
        gtpsf = psf(ap, refphase)
    else
        gtpsf = inputpsf
        yst = abs.(ifftshift(fft(lasty(sol))))
        valmask = 50 / 255 .< inputpsf .< 200 / 255
        svector = yst[valmask] .^ 2 ./ inputpsf[valmask]
        restoredsatpsf = yst .^ 2 ./ mean(svector)
    end
    phidiff, phierr = calculate_phi_err(solution(sol), refphase, mask)
    phirms = round(phasenorm(removepiston(removetiptilt(phidiff .+ refphase))); digits=3)
    phierr = round(phierr; digits=3)

    if isnothing(resmask)
        resmask = erode(erode(erode(mask)))
    end

    y = lasty(sol)
    ampy = abs.(ifftshift(y))
    # ampplot = ampy - ap
    ampplot2 = ampy .* resmask
    ap2 = mask2ap(resmask)
    norming = norm(filter(!isnan, ap .* resmask)) / norm(filter(!isnan, ampplot2))
    ampplot = ampy .* mask .* norming

    # resnumber = length(filter(x -> x < lowthr, ampplot2))
    resn = resnumber(sol, ap; resmask=resmask)
    # @show resnumber

    ga = fig[1, 1] = GridLayout()
    plotpad = 2
    phi = bboxview(
        PhaseRetrieval.removepiston(PhaseRetrieval.removetiptilt(phidiff .+ refphase)),
        mask,
        plotpad,
    )
    axphi, plotphi = showphasetight(phi, ga[1, 1]; hidedec=false)
    # axphierr, plotphierr = showphasetight(bboxview(phidiff .* mask |> PhaseRetrieval.removetiptilt |> PhaseRetrieval.removepiston), ga[2, 1], 
    # hidedec=false)

    ph = bboxview(
        resmask .*
        (PhaseRetrieval.removepiston(PhaseRetrieval.removetiptilt(phidiff .+ refphase))),
        mask,
        plotpad,
    )
    # ph = bboxview(angle.(ifftshift(solution(sol))) .* resmask, mask)  # count only residues inside resmask, but use indexes of mask
    # ph = bboxview(angle.(ifftshift(solution(sol))) .* mask)  # if the above doesn't work yet
    posx, posy, resmap = getresmapsparce(ph)
    # posx .+= plotpad
    # posy .+= plotpad
    respos = resmap .> 0
    resneg = resmap .< 0
    nphires = length(resmap)
    scatter!(axphi, posx[respos], posy[respos]; color="white", markersize=3)
    scatter!(axphi, posx[resneg], posy[resneg]; color="black", markersize=3)

    contourrange = map(x -> 0.5:1:x, size(bboxview(ampplot, mask, plotpad)))
    axy, hmy = heatmap(
        ga[2, 1],
        contourrange[2],
        contourrange[1],
        bboxview(ampplot, mask, plotpad);
        axis=(
            aspect=DataAspect(),
            title="Restored aperture",
            # title = L"Restored aperture, $|y^*[a > 0]|$"
        ),
        colormap=:YlGn_6,
        colorrange=(lowthr, 1.5),
        lowclip=:red,
    )
    contour!(
        axy,
        contourrange[2],
        contourrange[1],
        bboxview(mask2ap(resmask), mask, 2);
        color=:blue,
        linewidth=1,
        levels=[0, 1],
    )

    cb = Colorbar(
        ga[1, 1, Left()];
        # plotphi, 
        size=10,
        tellheight=false,
        vertical=true,
        flipaxis=false,
        limits=(-π, π),
        ticks=((-π):(π / 2):π, ["-π", "-π/2", "0", "π/2", "π"]),
        colormap=:cyclic_mygbm_30_95_c78_n256,
    )
    cb.halign = :left

    cb = Colorbar(
        ga[2, 1, Left()],
        hmy;
        size=10,
        tellheight=false,
        vertical=true,
        flipaxis=false,
        # limits=(-π, π),
        # ticks=(-π:π/2:π, ["-π", "-π/2", "0", "π/2", "π"]),
        # colormap=:cyclic_mygbm_30_95_c78_n256
    )
    cb.halign = :left

    axphi.title = "Restored phase"
    axphi.subtitle = "rms = $phirms, nres = $nphires"
    axy.subtitle = "$(resn) inner points below $lowthr"

    if isnothing(colormap)
        hl = highlightsaturation ∘ showpsf
    else
        hl(x) = get(colormap, float.(x), (0.0, float(sat)))
    end

    # gb = fig[1,2] = GridLayout()
    gtpsfim = image(
        ga[1, 2],
        rotr90(showpsf(gtpsf, 5));
        axis=(aspect=DataAspect(), title="Input PSF", subtitle="log scale"),
    )
    restoredpsfim = image(
        ga[2, 2],
        rotr90(showpsf(restoredpsf, 5));
        axis=(aspect=DataAspect(), title="Reconstructed PSF", subtitle="log scale"),
    )
    inputpsfim = image(
        ga[1, 3],
        rotr90(hl(inputpsf));
        axis=(aspect=DataAspect(), title="Input PSF", subtitle="linear scale"),
    )
    restoredsatpsfim = image(
        ga[2, 3],
        rotr90(hl(restoredsatpsf));
        axis=(
            aspect=DataAspect(),
            title="Reconstructed sat. pixels",
            subtitle="highlighted at levels [$(round(maximum(inputpsf), digits = 2)),$sat]",
        ),
    )

    if !isnothing(colormap)
        cb = Colorbar(
            ga[1:2, 3, Right()];
            colormap=colormap,
            size=10,
            tellheight=false,
            vertical=true,
            flipaxis=true,
            limits=(0, sat),
            # ticks=(-π:π/2:π, ["-π", "-π/2", "0", "π/2", "π"]),
            # colormap=:cyclic_mygbm_30_95_c78_n256
        )
        cb.halign = :right
    end

    if length(distgthist(sol)) > 0
        distances = distgthist(sol)[1:itersteps(sol)]
        distlabel = L"Phase error rms$"
        distlegend = "phase error"
    else
        distances = disthist(sol)[1:itersteps(sol)]
        distlabel = L"Distance $||x^{k+1}-y^{k}|| / \;\;\;  ||x^{k+1}||$"
        distlegend = "distance"
    end

    # ax1 = Makie.Axis(ga[3, :], xlabel="Iteration number", ylabel=L"Convergence $\tfrac{||x^{k}-x^{k-1}||}{||x^{k}||}$", yscale=log10, yticklabelspace=50.0)
    ax1 = Makie.Axis(
        ga[3, :];
        xlabel="Iteration number",
        ylabel=L"Convergence $||x^{k+1}-x^{k}|| / \;\;\;  ||x^{k}||$",
        yscale=log10,
        yticklabelspace=50.0,
    )
    ax2 = Makie.Axis(ga[3, :]; xlabel="x axis", ylabel=distlabel, yticklabelspace=35.0)
    convplot = lines!(ax1, errhist(sol)[1:itersteps(sol)]; label="convergence")
    distplot = lines!(ax2, distances; color=:orange, label="distance")
    ax2.yaxisposition = :right
    ax2.yticklabelalign = (:left, :center)
    ax2.xticklabelsvisible = false
    ax2.xlabelvisible = false
    linkxaxes!(ax1, ax2)

    axislegend(ax2, [convplot, distplot], ["convergence", distlegend])

    colsize!(ga, 1, Aspect(1, 1.0))
    colsize!(ga, 2, Aspect(1, 1.0))
    colsize!(ga, 3, Aspect(1, 1.0))
    return fig
end


function makelabel(name, input_params, preprocess_params, PRparams)
    lbl = (
        "$name, $(PRparams.PRalg), $((PRparams.apalg) 
 *  (PRparams.apalg == "DRAP" ? "($(PRparams.DRAPbeta))" : "")
 *  (PRparams.appolish ? "+AP" : "")) ds=$(preprocess_params.downsample)" *
        " t=$(preprocess_params.threshold)"
    )
    return lbl
end

"""
    BlurPhase(blurwidth, ntimes)

Construct an iterative algorithm which blurs the phase with given kernelDocument
"""
Base.@kwdef struct BlurPhaseparam <: IterativeAlgorithm
    blurwidth = missing
    x⁰ = missing
    maxϵ::Union{Float64,Missing} = missing
    maxit::Union{Missing,Int64} = missing
    keephistory::Bool = false
    snapshots::Array{Int64} = Int64[]
end

import AlternatingProjections.solve
function solve(
    p::Problem,
    alg::BlurPhaseparam,
    x⁰,
    maxϵ,
    maxit,
    keephistory::Bool,
    snapshots::Vector{Int64},
)

    # process the default parameters
    !ismissing(x⁰) || error("BlurPhase algorithm needs input array")
    !ismissing(maxϵ) || (maxϵ = 1e-15)
    !ismissing(maxit) || (maxit = 1)

    img = angle.(fftshift(x⁰))
    blurop = Opiterator(img, x -> imfilter(x, Kernel.gaussian(alg.blurwidth)))
    blurop = Iterators.take(blurop, maxit)

    state = loop(blurop)

    img .= ifftshift(state.xᵏ)

    # This only blurs the phase and does not update the intensity
    # return abs.(x⁰) .* exp.(im .* img), (lasty = x⁰, errhist = Vector{Float64}[], xhist = typeof(x⁰)[], disthist = Vector{Float64}[], k= 0)
    # this will use the known amplitude
    return AlternatingProjections.amp(p.A) .* exp.(im .* img),
    (
        lasty=x⁰,
        errhist=Vector{Float64}[],
        xhist=typeof(x⁰)[],
        disthist=Vector{Float64}[],
        k=0,
    )
end

function loop(iter)
    x = nothing
    for y in iter
        x = y
    end
    return x
end

function read_sim_data(d::Dict; ext = "png")
    @unpack config, phasetype, cropsim, sat, N, bpp = d
    crop = cropsim
    dload = @dict(config, phasetype, crop, sat, N, bpp)
    psfshortnname = savename("PSF", dload, ext)
    psfname = datadir("sims", phasedir, "psfs", psfshortnname)
    psfimage = load(psfname)

    d_phaseap = @dict(config, crop, phasetype)
    filename = datadir("sims", phasedir, "phases", savename("phaseap", d_phaseap, "hdf5"))
    @unpack ap, phase = load(filename)

    return psfshortnname, psfimage, ap, phase
end

function read_exp_data(d::Dict)
    @unpack sat, psfnamefunc = d

    psfshortnname = psfnamefunc(sat)
    psfname = datadir("exp", phasedir, psfshortnname)
    psfimage = load(psfname)

    return psfshortnname, psfimage
end



neighbourdict(d, (key, value)) = merge(d, Dict(key => value))

sameval(dict, sel) = all([dict[key] == sel[key] for key in keys(sel)])
filter_dict_list(dict, sel) = filter(d -> sameval(d, sel), dict)

safediv(x, y) = (x ≈ 0 && y ≈ 0) ? 1.0 : x / y