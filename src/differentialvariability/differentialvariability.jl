module DifferentialVariability

export dv_ftest, dv_levene, dv_brownforsythe, dv_splinefit

using HypothesisTests: VarianceFTest, UnequalVarianceTTest, pvalue
using MultipleTesting
using Statistics: mean, std, median, var, quantile
using Dierckx: Spline1D

# ── Existing tests ─────────────────────────────────────────────────────────────

"""
    dv_ftest(X, Y) -> (p, q)

Parametric F-test for equality of variance for each gene (row).
Assumes normally distributed expression. For heavy-tailed or zero-inflated
data (typical in scRNA-seq), prefer `dv_brownforsythe`.
Returns raw p-values and Benjamini-Hochberg adjusted q-values.
"""
function dv_ftest(X::AbstractMatrix, Y::AbstractMatrix)
    @assert size(X,1) == size(Y,1)
    n = size(X,1)
    p = zeros(n)
    for k in 1:n
        p[k] = pvalue(VarianceFTest(float.(X[k,:]), float.(Y[k,:])))
    end
    q = adjust(p, BenjaminiHochberg())
    return p, q
end

"""
    dv_levene(X, Y) -> (p, q)

Levene's test for equality of variance for each gene (row), centering on
group means. Returns raw p-values and BH-adjusted q-values.
"""
function dv_levene(X::AbstractMatrix, Y::AbstractMatrix)
    @assert size(X,1) == size(Y,1)
    _dv_absdev(X, Y, mean)
end

"""
    dv_brownforsythe(X, Y) -> (p, q)

Brown-Forsythe test for equality of variance for each gene (row).
Variant of Levene's test using group medians as the center — more robust
to the heavy-tailed, zero-inflated distributions typical in scRNA-seq data.
Returns raw p-values and BH-adjusted q-values.
"""
function dv_brownforsythe(X::AbstractMatrix, Y::AbstractMatrix)
    @assert size(X,1) == size(Y,1)
    _dv_absdev(X, Y, median)
end

# Levene/Brown-Forsythe: Welch t-test on absolute deviations from group center.
function _dv_absdev(X::AbstractMatrix, Y::AbstractMatrix, center)
    n = size(X,1)
    p = zeros(n)
    for k in 1:n
        xk = float.(X[k,:])
        yk = float.(Y[k,:])
        zx = abs.(xk .- center(xk))
        zy = abs.(yk .- center(yk))
        p[k] = pvalue(UnequalVarianceTTest(zx, zy))
    end
    q = adjust(p, BenjaminiHochberg())
    return p, q
end

# ── Splinefit method [PMID:40113778] ──────────────────────────────────────────

# Library-size normalization (local copy; avoids cross-module dependency).
function _libsize_norm(X::AbstractMatrix)
    libsz = vec(sum(X, dims=1))
    γ     = mean(libsz)
    return X ./ libsz' .* γ
end

# Per-gene statistics: log1p(mean), log1p(CV), dropout rate.
function _genestat(X::AbstractMatrix)
    u     = vec(mean(X, dims=2))
    σ     = vec(std(X,  dims=2))
    cv    = σ ./ max.(u, eps())
    lgu   = log1p.(u)
    lgcv  = log1p.(cv)
    dropr = vec(mean(iszero.(X), dims=2))
    valid = isfinite.(lgu) .& isfinite.(lgcv)
    return lgu, lgcv, dropr, valid
end

# Cumulative arc-length along the sorted (lgu, lgcv, dropr) path.
function _arc_length(lgu, lgcv, dropr)
    n = length(lgu)
    s = zeros(n)
    for i in 2:n
        s[i] = s[i-1] + sqrt((lgu[i]-lgu[i-1])^2 +
                              (lgcv[i]-lgcv[i-1])^2 +
                              (dropr[i]-dropr[i-1])^2)
    end
    return s
end

# Nearest point on the fitted spline curve for each data point.
# Uses a bounded window search — exact for smooth curves, O(n·window).
function _nearest_on_spline(xyz_data::Matrix, xyz_spline::Matrix; window::Int=150)
    n       = size(xyz_data, 1)
    nearidx = zeros(Int, n)
    d       = zeros(n)
    for i in 1:n
        lo     = max(1, i - window)
        hi     = min(n, i + window)
        best_d = Inf
        best_j = i
        for j in lo:hi
            dij = sqrt((xyz_data[i,1]-xyz_spline[j,1])^2 +
                       (xyz_data[i,2]-xyz_spline[j,2])^2 +
                       (xyz_data[i,3]-xyz_spline[j,3])^2)
            if dij < best_d
                best_d = dij
                best_j = j
            end
        end
        nearidx[i] = best_j
        d[i]       = best_d
    end
    return nearidx, d
end

# Fit spline to one group; returns per-gene stats in original gene order.
# Mirrors sc_splinefit.m: arc-length parameterisation → Dierckx smoothing
# spline (k=3, s=n*smooth, matching splinefit(s, xyz', 15, 0.75)) → nearest
# point on curve → distance adjustments → half-normal p-values.
function _splinefit_group(X::AbstractMatrix; smooth::Float64=0.75)
    lgu, lgcv, dropr, valid = _genestat(X)
    ng = length(lgu)
    vi = findall(valid)

    lgu_v, lgcv_v, dropr_v = lgu[vi], lgcv[vi], dropr[vi]

    # Sort genes by (lgu, lgcv, dropr) to order them along the feature curve.
    order   = sortperm(collect(zip(lgu_v, lgcv_v, dropr_v)))
    lgu_s   = lgu_v[order]
    lgcv_s  = lgcv_v[order]
    dropr_s = dropr_v[order]

    s  = _arc_length(lgu_s, lgcv_s, dropr_s)
    n  = length(s)
    ds = n * smooth   # Dierckx smoothing factor: larger → smoother curve

    # Fit one cubic smoothing spline per coordinate (matches MATLAB splinefit).
    spl_u  = Spline1D(s, lgu_s;   k=3, s=ds, bc="nearest")
    spl_cv = Spline1D(s, lgcv_s;  k=3, s=ds, bc="nearest")
    spl_dr = Spline1D(s, dropr_s; k=3, s=ds, bc="nearest")

    fit_lgu   = spl_u.(s)
    fit_lgcv  = spl_cv.(s)
    fit_dropr = spl_dr.(s)

    xyz_fit    = hcat(fit_lgu, fit_lgcv, fit_dropr)
    xyz_sorted = hcat(lgu_s, lgcv_s, dropr_s)

    # For each gene find the nearest point on the fitted curve.
    nearidx_s, d_s = _nearest_on_spline(xyz_sorted, xyz_fit)

    # Distance adjustments (mirrors sc_splinefit.m):
    #   genes beyond the high/low-mean range get de-emphasized;
    #   genes below the CV curve (less variable than fitted) get strongly
    #   de-emphasized — they should not appear as differentially variable.
    fitmean = fit_lgu
    d_s[lgu_s .> maximum(fitmean)] ./= 100
    d_s[lgu_s .< minimum(fitmean)] ./= 10
    d_s[(lgcv_s .- fit_lgcv) .< 0] ./= 100

    # P-values: fit a symmetric half-normal to distances ≤ 90th percentile,
    # then upper-tail probability via erfc (= 1 − normcdf, no Distributions.jl).
    q90 = quantile(d_s, 0.9)
    dx  = d_s[d_s .<= q90]
    σ   = max(std(vcat(-dx, dx)), eps())
    pval_s = 0.5 .* erfc.(d_s ./ (σ * sqrt(2)))

    nearpt_s = xyz_fit[nearidx_s, :]   # coordinates of nearest spline point

    # Map: sorted-valid → unsorted-valid → full ng-length arrays
    inv_order          = invperm(order)
    d_v                = d_s[inv_order]
    pval_v             = pval_s[inv_order]
    nearpt_v           = nearpt_s[inv_order, :]

    d_out              = fill(NaN, ng);     d_out[vi]         = d_v
    pval_out           = fill(NaN, ng);     pval_out[vi]      = pval_v
    nearpt_out         = fill(NaN, ng, 3);  nearpt_out[vi, :] = nearpt_v

    return (lgu=lgu, lgcv=lgcv, dropr=dropr, valid=valid,
            nearpt=nearpt_out, d=d_out, pval=pval_out)
end

"""
    dv_splinefit(X1, X2; g, smooth) -> NamedTuple

Differential variability analysis between two cell groups using the splinefit
method [PMID:40113778].

Each group's genes are characterised in 3-D feature space
(log1p-mean, log1p-CV, dropout-rate). A Dierckx cubic smoothing spline
(matching MATLAB's `splinefit`) is fit through the genes ordered by arc length.
The deviation of each gene from its group's spline is compared between groups.

# Arguments
- `X1`, `X2`  : gene-by-cell expression matrices (rows must match in number and order)
- `g`         : shared gene name vector (optional; indices used if omitted)
- `smooth`    : smoothing strength ∈ [0, 1] (default 0.75; higher = smoother curve)

# Returns
Named tuple sorted by `DiffDist` descending:
- `gene`      : gene names
- `lgu1/2`    : log1p(mean expression) per group
- `lgcv1/2`   : log1p(coefficient of variation) per group
- `dropr1/2`  : dropout rate per group
- `d1/2`      : deviation from group spline (adjusted)
- `DiffDist`  : ‖v₁ − v₂‖ where vᵢ = gene_pt − nearest_spline_pt in group i
- `DiffSign`  : +1 more variable in group 1, −1 more variable in group 2
- `pval`      : one-sided p-value from z-score of DiffDist
"""
function dv_splinefit(X1::AbstractMatrix, X2::AbstractMatrix;
                      g::Union{Vector{String}, Nothing}=nothing,
                      smooth::Float64=0.75)
    @assert size(X1,1) == size(X2,1) "X1 and X2 must have the same number of genes"
    ng    = size(X1,1)
    genes = isnothing(g) ? ["Gene$i" for i in 1:ng] : g

    X1n = _libsize_norm(float.(X1))
    X2n = _libsize_norm(float.(X2))

    r1 = _splinefit_group(X1n; smooth)
    r2 = _splinefit_group(X2n; smooth)

    valid = r1.valid .& r2.valid
    vi    = findall(valid)

    lgu1_v, lgcv1_v, dropr1_v = r1.lgu[vi], r1.lgcv[vi], r1.dropr[vi]
    lgu2_v, lgcv2_v, dropr2_v = r2.lgu[vi], r2.lgcv[vi], r2.dropr[vi]

    # Deviation vectors: gene point − nearest spline point
    v1 = hcat(lgu1_v, lgcv1_v, dropr1_v) .- r1.nearpt[vi, :]
    v2 = hcat(lgu2_v, lgcv2_v, dropr2_v) .- r2.nearpt[vi, :]

    DiffDist = vec(sqrt.(sum((v1 .- v2).^2, dims=2)))
    DiffSign = sign.(vec(sqrt.(sum(v1.^2, dims=2))) .-
                     vec(sqrt.(sum(v2.^2, dims=2))))

    # Z-score → upper-tail p-value (mirrors sc_dvg.m normcdf logic)
    μ, σ = mean(DiffDist), std(DiffDist)
    ddz  = (DiffDist .- μ) ./ max(σ, eps())
    pval = 0.5 .* erfc.(ddz ./ sqrt(2))

    srt = sortperm(DiffDist, rev=true)

    return (
        gene     = genes[vi][srt],
        lgu1     = lgu1_v[srt],    lgcv1  = lgcv1_v[srt],   dropr1 = dropr1_v[srt],
        d1       = r1.d[vi][srt],
        lgu2     = lgu2_v[srt],    lgcv2  = lgcv2_v[srt],   dropr2 = dropr2_v[srt],
        d2       = r2.d[vi][srt],
        DiffDist = DiffDist[srt],
        DiffSign = DiffSign[srt],
        pval     = pval[srt],
    )
end

end
