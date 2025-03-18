const MAX_CORNERPLOT_DIMS_BEFORE_PCA = 20
const BONs = BiodiversityObservationNetworks

# NOTES:
# The easiest way to make this compatable with GeoMakie (if loaded) or 
# else default back to normal Axis is to have it work on a ::Makie.GridPosition,
# which 


# TODO: pass any transform. Annoying because its a kwarg and could be of type
# <:StatsBase.AbstractDataTransform (as are some things provided by StatsBase,
# e.g. ZScoreTransform, and some things provided by MVStats, e.g. Whitening), 
# but also of type MVStats.AbstractDimensionalityReductoin (e.g. PCA, PPCA)

function BONs.cornerplot(
    layers::RasterStack;
    pca_layers = false,
    sz = (1600,1600)
)

    _, mat = BiodiversityObservationNetworks.features(layers)
    num_layers = length(layers)
    if num_layers > MAX_CORNERPLOT_DIMS_BEFORE_PCA || pca_layers
        pca = BONs.MultivariateStats.fit(BiodiversityObservationNetworks.MultivariateStats.PCA, mat)
        num_layers = length(pca.prinvars)
        mat = BONs.MultivariateStats.transform(pca, mat)    
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


function BONs.bonplot(
    position::GridPosition, 
    bon::BiodiversityObservationNetwork;
    axistype = Makie.Axis
)
    ax = axistype(position)
    plot = scatter!(ax, [node[1] for node in bon], [node[2] for node in bon], color=(:red))
    Makie.AxisPlot(ax, plot)
end

function BONs.bonplot(
    position::GridPosition,
    bon::BiodiversityObservationNetwork,
    geom::T;
    axistype = Makie.Axis
) where T
    GEOM_TYPE = BONs._what_did_you_pass(geom)
    isnothing(GEOM_TYPE) && throw(ArgumentError("$T cannot be coerced to a valid Geometry"))
    BONs.bonplot(position, bon, Base.convert(GEOM_TYPE, geom); axistype=axistype)
end


function BONs.bonplot(
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


function BONs.bonplot(
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


function BONs.bonplot(
    position::GridPosition,
    bon::BiodiversityObservationNetwork,
    raster::Raster;
    axistype=Makie.Axis
)
    ax = axistype(position)
    heatmap!(ax, raster.raster)
    plot = scatter!(ax, [node.coordinate for node in bon], color=(:red))
    Makie.AxisPlot(ax, plot)
end


function Makie.voronoiplot(
    position::GridPosition,
    bon::BiodiversityObservationNetwork, 
    geom::Polygon;
    axistype = Makie.Axis
)
    vor = voronoi(bon, geom)
    
    ax = axistype(position)
    hidedecorations!(ax)
    map(v->poly!(ax, v, strokewidth=1), vor)
    if size(bon) < 50
        scatter!(ax, [n.coordinate for n in bon], color=:red)
    end
    poly!(ax, geom, color=(:white, 0), strokewidth=1)
end 


# Makie poly overloads
Makie.poly(polygon::Polygon; kw...) = poly(polygon.geometry; kw...)
Makie.poly(polygons::Vector{Polygon}; kw...) = begin
    poly(first(polygons); kw...)
    map(p->poly!(p; kw...), polygons[2:end])
    current_figure()
end 
Makie.poly!(polygon::Polygon; kw...) = poly!(polygon.geometry; kw...)
Makie.poly!(polygons::Vector{Polygon}; kw...) = begin
    map(p->poly!(p; kw...), polygons)
    current_figure()
end 
Makie.poly!(ax, polygon::Polygon; kw...) = poly!(ax, polygon.geometry; kw...)
Makie.poly!(ax, polygons::Vector{Polygon}; kw...) = map(p->poly!(ax,p; kw...), polygons)



