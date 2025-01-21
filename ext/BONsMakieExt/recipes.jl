const MAX_CORNERPLOT_DIMS_BEFORE_PCA = 20

# NOTES:
# The easiest way to make this compatable with GeoMakie (if loaded) or 
# else default back to normal Axis is to have it work on a ::Makie.GridPosition,
# which 


# TODO: pass any transform. Annoying because its a kwarg and could be of type
# <:StatsBase.AbstractDataTransform (as are some things provided by StatsBase,
# e.g. ZScoreTransform, and some things provided by MVStats, e.g. Whitening), 
# but also of type MVStats.AbstractDimensionalityReductoin (e.g. PCA, PPCA)

function BiodiversityObservationNetworks.cornerplot(
    layers::RasterStack;
    pca_layers = false,
    sz = (1600,1600)
)

    _, mat = BiodiversityObservationNetworks.features(layers)
    num_layers = length(layers)
    if num_layers > MAX_CORNERPLOT_DIMS_BEFORE_PCA || pca_layers
        pca = BiodiversityObservationNetworks.MultivariateStats.fit(BiodiversityObservationNetworks.MultivariateStats.PCA, mat)
        @info length(pca.prinvars)
        num_layers = length(pca.prinvars)
        mat = BiodiversityObservationNetworks.MultivariateStats.transform(pca, mat)    
    end 

    f = Figure(size=sz)
    for i in 1:num_layers-1, j in 1:num_layers
        if j > i
            ax = Axis(
                f[j-1,i],
                xlabel = j == num_layers ? "$i" : "",
                ylabel = i == 1 ? "$j" : "",
                xticksvisible=false,
                yticksvisible=false,
                xticklabelsvisible=false,
                yticklabelsvisible=false,
            )
            hexbin!(ax, mat[i,:], mat[j,:], bins=40)
        end
    end 
    f
end


function BiodiversityObservationNetworks.bonplot(
    position::GridPosition, 
    bon::BiodiversityObservationNetwork;
    axistype = Makie.Axis
)
    ax = axistype(position)
    plot = scatter!(ax, [node[1] for node in bon], [node[2] for node in bon], color=(:red))
    Makie.AxisPlot(ax, plot)
end


function BiodiversityObservationNetworks.bonplot(
    position::GridPosition,
    bon::BiodiversityObservationNetwork,
    poly::Any;
    axistype = Makie.Axis
)
    bonspoly = BiodiversityObservationNetworks._convert_to_bons_polygon(poly)
    BiodiversityObservationNetworks.bonplot(position, bon, bonspoly; axistype=axistype)
end


function BiodiversityObservationNetworks.bonplot(
    position::GridPosition,
    bon::BiodiversityObservationNetwork,
    poly::Polygon;
    axistype=Makie.Axis
)
    ax = axistype(position)
    poly!(ax, poly.geometry, strokewidth=2, color=(:grey, 0.1))
    plot = scatter!(ax, [node[1] for node in bon], [node[2] for node in bon], color=(:red))
    Makie.AxisPlot(ax, plot)
end


function BiodiversityObservationNetworks.bonplot(
    position::GridPosition,
    bon::BiodiversityObservationNetwork,
    poly::Vector{<:Polygon};
    axistype=Makie.Axis
)
    ax = axistype(position)
    for p in poly
        poly!(ax, p.geometry, strokewidth=1, color=(:grey, 0.1))
    end
    plot = scatter!(ax, [node[1] for node in bon], [node[2] for node in bon], color=(:red))
    Makie.AxisPlot(ax, plot)
end




function BiodiversityObservationNetworks.bonplot(
    position::GridPosition,
    bon::BiodiversityObservationNetwork,
    raster::Raster;
    axistype=Makie.Axis
)
    ax = axistype(position)
    heatmap!(ax, raster.raster)
    plot = scatter!(ax, [node[1] for node in bon], [node[2] for node in bon], color=(:red))
    Makie.AxisPlot(ax, plot)
end