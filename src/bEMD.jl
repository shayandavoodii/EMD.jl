function midpoints(s, t)
    _, _, exts = find_extrema(s)
    n   = length(exts)
    mpx = similar(s, n-1)
    mpy = similar(mpx)
    for i=1:n-1
        mpx[i] = (t[exts[i]] + t[exts[i+1]]) /2.0
        mpy[i] = (s[exts[i]] + s[exts[i+1]]) /2.0
    end
    mpx,mpy
end

function get_spline(s, t, sparam)
    mpt, mps = midpoints(s, t)
    spl      = Spline1D(mpt, mps, bc="extrapolate", s=sparam)
    return spl(t)
end

function bsift(s, t; sparam=0.0)
    count = 1
    sd    = 1.0
    hlast = copy(s)
    hnew  = copy(s)
    while sd > 0.25
        spls     = get_spline(hlast, t, sparam)
        @. hnew  = hlast - spls
        sd       = sum(abs2.(hlast-hnew))/sum(hlast.^2)
        @. hlast = hnew
        count   +=1
    end
    return hnew
end

function bEMD(s, t; maximfs=5, sparam=0.0)
    sig  = copy(s)
    N    = length(s)
    imfs = ElasticArray{Float64}(undef, N, 0)
    for i in 1:maximfs
        h = bsift(sig, t; sparam=sparam)
        append!(imfs, h)
        sig .-= h
    end
    return convert(Array, imfs)
end
