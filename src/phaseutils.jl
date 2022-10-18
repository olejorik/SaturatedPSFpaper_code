using Statistics

phasenorm(ϕ) = std(phwrap.(filter(!isnan, ϕ)); corrected=false)

function calculate_phi_err(x, phase, mask; remtiptilt=true)
    ϕ = angle.(ifftshift(x))
    ϕdiff1 = (ϕ .- phase) .* mask
    ϕdiff2 = (twinphase(ϕ) .- phase) .* mask # in case the even part has wrong sign
    if remtiptilt
        cleanphase = removepiston ∘ removetiptilt
    else
        cleanphase = removepiston
    end
    ϕdiff1 = cleanphase(ϕdiff1)
    ϕdiff2 = cleanphase(ϕdiff2)
    ϕerr1 = phasenorm(ϕdiff1)
    ϕerr2 = phasenorm(ϕdiff2)
    ϕdiff = ϕerr1 < ϕerr2 ? ϕdiff1 : ϕdiff2
    ϕerr = min(ϕerr1, ϕerr2)

    return ϕdiff, ϕerr
end

function select_twin(ϕ, refphase, mask=nothing; remtiptilt=true)
    if isnothing(mask)
        mask = zero(ϕ) .+ 1
    end
    ϕdiff1 = (ϕ .- refphase) .* mask
    ϕdiff2 = (twinphase(ϕ) .- refphase) .* mask # in case the even part has wrong sign
    if remtiptilt
        cleanphase = removepiston ∘ removetiptilt
    else
        cleanphase = removepiston
    end
    ϕdiff1 = cleanphase(ϕdiff1)
    ϕdiff2 = cleanphase(ϕdiff2)
    ϕerr1 = phasenorm(ϕdiff1)
    ϕerr2 = phasenorm(ϕdiff2)
    ret = ϕerr1 < ϕerr2 ? ϕ : twinphase(ϕ)

    return ret
end

function resnumber(sol, ap; resmask=nothing, lowthr=0.5)
    y = lasty(sol)
    ampy = abs.(ifftshift(y))
    if isnothing(resmask)
        resmask = erode(erode(erode(ap2mask(ap))))
    end

    ampplot2 = ampy .* resmask
    norming = norm(filter(!isnan, ampplot2)) / norm(filter(!isnan, ap .* resmask))

    return length(filter(x -> x < lowthr * norming, ampplot2))
end

function resphinumber(sol, ap; resmask=nothing, lowthr=0.5)
    ph = angle.(ifftshift(solution(sol)))
    if isnothing(resmask)
        ph = ph .* erode(erode(erode(ap2mask(ap))))
    else
        ph = ph .* resmask
    end
    _, _, resmap = getresmapsparce(bboxview(ph))

    return length(resmap)
end

function resblock(A)
    return sum(
        phwrap, [A[2, 1] - A[1, 1], A[2, 2] - A[2, 1], A[1, 2] - A[2, 2], A[1, 1] - A[1, 2]]
    )
end

function getresmap(phase, positions=true)
    I, J = size(phase)
    resmap = similar(phase, I - 1, J - 1)
    posx = similar(phase, I - 1, J - 1)
    posy = similar(phase, I - 1, J - 1)
    for i in 1:(I - 1), j in 1:(J - 1)
        resmap[i, j] = resblock(phase[i:(i + 1), j:(j + 1)]) ÷ (2π)
        posx[i, j] = i + 1 / 2
        posy[i, j] = j + 1 / 2
    end
    resmap[resmap .== 0] .= NaN
    return posy, posx, resmap
end

function getresmapsparce(phase)
    I, J = size(phase)
    resmap = Int[]
    posx = Float16[]
    posy = Float16[]
    for i in 1:(I - 1), j in 1:(J - 1)
        res = resblock(phase[i:(i + 1), j:(j + 1)]) ÷ (2π)
        if !isnan(res) && res != 0
            push!(resmap, res)
            push!(posx, i + 1 / 2)
            push!(posy, j + 1 / 2)
        end
    end
    return posx, posy, resmap
end
