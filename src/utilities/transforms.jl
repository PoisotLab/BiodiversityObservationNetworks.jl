
# -----------------------------------------------------------------------------
# Begin ZScore 
function _zscore_transform_stack(z, stack)
    zstack = deepcopy(stack)
    ztransform = BONs.StatsBase.transform(z, features(stack)[2])
    for i in eachindex(stack)
        zstack[i].raster.grid[zstack[i].raster.indices] .= ztransform[i,:]
    end
    return zstack
end

function _zscore_baseline(baseline)
    _, mat = features(baseline)
    z = BONs.StatsBase.fit(BONs.StatsBase.ZScoreTransform, mat)     
    return z
end
# End  ZScore 
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Begin PCA 

function _pca_transform_stack(pca, stack)
    L = map(r->r.raster, stack)
    pcad_layers = RasterStack(BONs.MVStats.predict(pca, L))
    return pcad_layers
end

function _fit_baseline_pca(baseline)
    L = map(r->r.raster, baseline)
    pca = BONs.MVStats.fit(BONs.MVStats.PCA, L)
    return pca
end

# End PCA
# -----------------------------------------------------------------------------