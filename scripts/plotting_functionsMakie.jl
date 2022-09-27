using CairoMakie
using LaTeXStrings

showpsf(psf) = Gray.(PhaseRetrieval.rescale(psf))
showpsf(psf, α) = Gray.(PhaseRetrieval.logrescale(float(psf), α))

function highlightsaturation(img)
    mask = img .>= Gray(1.0)
    img = RGB.(img)
    img[mask] .= RGB(1, 0, 0)
    return img
end

function highlightsaturation(img, level)
    mask = img .>= Gray(level)
    img = RGB.(img)
    img[mask] .+= RGB(0.5, 0, 0)
    return img
end

function showarray(arr)
    return image(rotr90(Gray.(PhaseRetrieval.rescale(arr))); axis=(
        # yflip=false, 
        aspect=DataAspect(),
        # xlims=0.5 .+ (0, size(arr, 2)),
        # ylims=0.5 .+ (0, size(arr, 1)),
        # frame = :box)
    ))
end

function showarray(arr, α)
    return plot(Gray.(PhaseRetrieval.logrescale(arr, α)); axis=(
        # yflip=false, 
        aspect=DataAspect(),
        # xlims=0.5 .+ (0, size(arr, 2)),
        # ylims=0.5 .+ (0, size(arr, 1)),
        # frame = :box
    ))
end

function showarray_max(arr)
    return plot(
        Gray.(PhaseRetrieval.rescale(arr));
        yflip=false,
        aspect_ratio=1,
        borders=:none,
        size=(512, 512),
        annotations=(
            size(arr, 2) / 2.0,
            25,
            text("maxvalue = $(round(maximum(arr), sigdigits=3))", 16, :middle, :white),
        ),
    )
end

function showarray_max(arr, α)
    return plot(
        Gray.(PhaseRetrieval.logrescale(arr, α));
        yflip=false,
        aspect_ratio=1,
        borders=:none,
        size=(512, 512),
        annotations=(
            size(arr, 2) / 2.0,
            25,
            text("maxvalue = $(round(maximum(arr), sigdigits=3))", 16, :middle, :white),
        ),
    )
end

function showphase(inarr, fig=Figure(), picsize=512, cm=:cyclic_mygbm_30_95_c78_n256)
    if max(size(inarr)...) > picsize
        arr = imresize(inarr, picsize)
    else
        arr = inarr
    end
    rows, cols = size(arr)
    # fig, ax, hm = heatmap(phwrap.(arr), show=true, colormap =cm, 
    # size=(picsize, picsize),
    # aspect_ratio=cols / rows,
    # xaxis=([(0.5, cols + 0.5), L"\sigma_x", , ([0.5, cols / 2 + 0.25, cols + 0.5], [L"-\frac{1}{2s}", "0", L"\frac{1}{2s}"])),
    # yaxis=(L"\sigma_y", (0.5, rows + 0.5), ([0.5, rows / 2 + 0.25, rows + 0.5], [L"-\frac{1}{2s}", "0", L"\frac{1}{2s}"])),
    # right_margin = 5Plots.mm # quick fix
    # )
    # plotgrid = fig[1,1] = GridLayout()
    ax = CairoMakie.Axis(
        fig[1, 1];
        aspect=1,
        #    autolimitaspect = 1,
        #    xlabel = L"\sigma_x",
        #    ylabel = L"\sigma_y",
        #    xticks = ([0.5, cols / 2 + 0.5, cols + 0.5], [L"-\frac{1}{2s}", L"0", L"\frac{1}{2s}"]),
        #    yticks = ([0.5, rows / 2 + 0.5, rows + 0.5], [L"-\frac{1}{2s}", L"0", L"\frac{1}{2s}"])
    )
    hm = heatmap!(ax, phwrap.(arr); colormap=cm)
    cb = Colorbar(fig[1, 2], hm; width=10, tellheight=true)
    return ax, cb
end

function showphasetight(
    inarr,
    fig=Figure();
    picsize=512,
    cm=:cyclic_mygbm_30_95_c78_n256,
    hidedec=true,
    kwarg...,
)
    if max(size(inarr)...) > picsize
        arr = imresize(inarr, picsize)
    else
        arr = inarr
    end
    rows, cols = size(arr)

    if typeof(fig) == GridPosition
        pos = fig
    else
        pos = fig[1, 1]
    end
    ax = CairoMakie.Axis(
        pos;
        aspect=AxisAspect(1),
        #    autolimitaspect = 1,
        #    xlabel = L"\sigma_x",
        #    ylabel = L"\sigma_y",
        #    xticks = ([0.5, cols / 2 + 0.5, cols + 0.5], [L"-\frac{1}{2s}", L"0", L"\frac{1}{2s}"]),
        #    yticks = ([0.5, rows / 2 + 0.5, rows + 0.5], [L"-\frac{1}{2s}", L"0", L"\frac{1}{2s}"])
    )
    # contourrange = map(x-> 0.5:1:x, size(arr))
    hm = heatmap!(ax, phwrap.(arr); colormap=cm, colorrange=(-π, π), kwarg...)
    # hm = heatmap!(ax, contourrange[2], contourrange[1], phwrap.(arr), colormap = cm, colorrange =(-π,π), kwarg...)
    if hidedec
        hidedecorations!(ax; grid=false)
    end
    # cb = Colorbar(fig[1,2], hm, width = 10, tellheight=true)
    # cb.ticks = (-π:π/2:π, ["-π","-π/2", "0","π/2","π"])
    return ax, hm
end

function ap2mask(ap)
    mask = zero(ap) .+ 1
    mask[ap .== 0] .= NaN
    return mask
end

function mask2ap(mask)
    ap = zero(mask) .+ 1
    ap[isnan.(mask)] .= 0
    return ap
end

showphase(phaseap::Tuple, args...) = showphase(phaseap[1] .* ap2mask(phaseap[2]), args...)

function complexcolor(z)
    return HSV(angle(z) / pi * 180, abs(z), 1)
end

function bboxview(arr)
    idx = findall(!isnan, arr)
    region = minimum(idx):maximum(idx)
    return @view arr[region]
end

"""
    bboxview(arr, mask, pad = 0)

Make a box corresponding to the mask not-NaN elements surrounded by pad.
"""
function bboxview(arr, mask, pad=0)
    idx = findall(!isnan, mask)
    cpad = CartesianIndex(pad, pad)
    region = (minimum(idx) - cpad):(maximum(idx) + cpad)
    return @view arr[region]
end