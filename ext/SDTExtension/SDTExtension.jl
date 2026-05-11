module SDTExtension
    using SpeciesDistributionToolkit
    using BiodiversityObservationNetworks
    using TestItems
    using BiodiversityObservationNetworks.StatsBase

    _apply_sdm_mask!(::BitMatrix, ::SDMLayer, ::Missing) = nothing
    function _apply_sdm_mask!(valid::BitMatrix, layer::SDMLayer, mask::AbstractMatrix)
        size(mask) == size(layer) || throw(ArgumentError("Mask size $(size(mask)) must match layer size $(size(layer))"))
        valid .&= Bool.(mask)
        return
    end
    function _apply_sdm_mask!(valid::BitMatrix, layer::SDMLayer, mask::SDMLayer)
        SimpleSDMLayers._layers_are_compatible(layer, mask) 
        valid .&= mask.indices
        return
    end
    function _apply_sdm_mask!(valid::BitMatrix, layer::SDMLayer, mask::SpeciesDistributionToolkit.SimpleSDMPolygons.AbstractGeometry)
        masked_layer = SpeciesDistributionToolkit.mask(layer, mask)
        valid .&= masked_layer.indices
        return
    end

    function _sdm_coordinates(layer::SDMLayer, keys::Vector{CartesianIndex{2}})
        Es = SimpleSDMLayers.eastings(layer)
        Ns = SimpleSDMLayers.northings(layer)
        coords = hcat([[Es[k[2]], Ns[k[1]] ] for k in keys]...)
        return coords
    end

    function BiodiversityObservationNetworks._extract_and_process_inclusion(inclusion::SDMLayer, keys, n)
        BiodiversityObservationNetworks._process_inclusion([inclusion[k] for k in keys], n)
    end

    """
        CandidatePool(layer::SDMLayer; mask=missing, inclusion=missing)

    Convert an `SDMLayer` into a [`CandidatePool`](@ref). Valid cells (from
    `layer.indices`) become candidates with geographic (lon, lat) coordinates.

    An optional `mask` (a `BitMatrix` or another `SDMLayer`) restricts which
    cells are included. An optional `inclusion` provides per-cell weights.
    """
    function BiodiversityObservationNetworks.CandidatePool(
        layer::SDMLayer; 
        mask=missing, 
        inclusion=missing
    )
        valid = copy(layer.indices)
        _apply_sdm_mask!(valid, layer, mask)

        if !ismissing(inclusion)
            size(inclusion) == size(layer) || throw(ArgumentError(
            """
            Inclusion must have same resolution as input layer. Received a $(size(layer)) domain and a $(size(inclusion)) inclusion.
            
            If you are using a polygon domain, the default rasterization size is (100, 100), but can be changed by passing the resolution keyword argument to sample.
            """
            ))
        end

        keys = vec(findall(valid))
        n = length(keys)
        n > 0 || throw(ArgumentError("No valid candidates after masking"))

        coords = _sdm_coordinates(layer, keys)
        incl = BiodiversityObservationNetworks._extract_and_process_inclusion(inclusion, keys, n)

        features = Matrix(layer[keys]')
        return CandidatePool(n, keys, coords, features, incl)
    end

    """
        CandidatePool(layers::Vector{<:SDMLayer}; mask=missing, inclusion=missing)

    Convert a vector of same-extent `SDMLayer`s into a [`CandidatePool`](@ref)
    with features. The valid pool is the intersection of valid cells across all
    layers. Each layer contributes one feature row. Coordinates are geographic
    (lon, lat).
    """
    function BiodiversityObservationNetworks.CandidatePool(
        layers::Vector{<:SDMLayer}; 
        mask=missing, 
        inclusion=missing
    )
        # Check if layers have same size, extent and crs
        SimpleSDMLayers._layers_are_compatible(layers)

        first_layer = first(layers)
        valid = reduce((a, b) -> a .& b, [l.indices for l in layers])
        _apply_sdm_mask!(valid, first_layer, mask)

        keys = vec(findall(valid))
        n = length(keys)
        n > 0 || throw(ArgumentError("No valid candidates after masking"))

        coords = _sdm_coordinates(first_layer, keys)
        features = Matrix{Float64}(undef, length(layers), n)
        for (row, l) in enumerate(layers)
            for (col, k) in enumerate(keys)
                features[row, col] = Float64(l[k])
            end
        end

        incl = BiodiversityObservationNetworks._extract_and_process_inclusion(inclusion, keys, n)

        return CandidatePool(n, keys, coords, features, incl)
    end
    
    """
        CandidatePool(poly::SimpleSDMPolygons.AbstractGeometry; resolution=0.5, mask=missing, inclusion=missing)
    
    Rasterize a polygon geometry into a [`CandidatePool`](@ref). A regular grid at
    the given `resolution` (in decimal degrees) is created over the polygon's
    bounding box, cells outside the polygon are excluded, and the result is returned
    as a pool with geographic (lon, lat) coordinates.

    Accepts any `SimpleSDMPolygons.AbstractGeometry`: `Polygon`, `MultiPolygon`,
    `Feature`, or `FeatureCollection`.

    ## Arguments
    - `poly`: the polygon to rasterize
    - `resolution`: the size to construct a raster 
    - `mask`: optional additional `BitMatrix` or `SDMLayer` mask applied after polygon rasterization
    - `inclusion`: optional per-cell weight `Matrix` or `Vector`
    """
    function BiodiversityObservationNetworks.CandidatePool(
        poly::SimpleSDMPolygons.AbstractGeometry;
        resolution = (100, 100),
        mask = missing,
        inclusion = missing,
    )
        bbox = SimpleSDMPolygons.boundingbox(poly)
        layer = SDMLayer(
            ones(resolution);
            x = (Float64(bbox.left),  Float64(bbox.right)),
            y = (Float64(bbox.bottom), Float64(bbox.top)),
        )

        # turns off cells whose centres fall outside the polygon
        mask!(layer, poly)

        return CandidatePool(layer; mask, inclusion)
    end

    @testitem "CandidatePool from Polygon" setup=[TestModule] begin
        poly = SDT.Polygon(
            (-5.0, 45.0), (5.0, 45.0), (5.0, 55.0), (-5.0, 55.0),
        )
        cp = CandidatePool(poly; resolution = (50, 30))
        @test cp isa CandidatePool
        @test cp.n > 0
        @test sum(cp.inclusion) ≈ 1.0
    end

    @testitem "CandidatePool from SDMLayer" setup=[TestModule] begin
        layer = SDT.SDMLayer(rand(Float32, 10, 15))
        cp = CandidatePool(layer)
        @test cp isa CandidatePool{CartesianIndex{2}}
        @test cp.n == 150
        @test size(cp.coordinates) == (2, 150)
        @test sum(cp.inclusion) ≈ 1.0
    end

    @testitem "CandidatePool from SDMLayer with matrix mask" setup=[TestModule] begin
        layer = SDT.SDMLayer(rand(Float32, 10, 10))
        mask = falses(10, 10)
        mask[1:5, :] .= true
        cp = CandidatePool(layer; mask)
        @test cp.n == 50
        @test all(k -> k[1] <= 5, cp.keys)
    end

    @testitem "CandidatePool from SDMLayer with SDMLayer mask" setup=[TestModule] begin
        layer = SDT.SDMLayer(rand(Float32, 10, 10))
        mask_layer = SDT.SDMLayer(rand(Float32, 10, 10))
        mask_layer.indices[6:end, :] .= false
        cp = CandidatePool(layer; mask=mask_layer)
        @test cp.n == 50
        @test all(k -> k[1] <= 5, cp.keys)
    end

    @testitem "CandidatePool from SDMLayer vector" setup=[TestModule] begin
        layers = [SDT.SDMLayer(rand(Float32, 8, 6)) for _ in 1:4]
        cp = CandidatePool(layers)
        @test cp isa CandidatePool{CartesianIndex{2}}
        @test cp.n == 48
        @test size(cp.features) == (4, 48)
        @test sum(cp.inclusion) ≈ 1.0
    end


    # =============================================================================
    # Climate rarity
    # =============================================================================


    """
        BiodiversityObservationNetworks.rarity(metric::RarityMetric, bon::BiodiversityObservationNetwork, layers::Vector{<:SDMLayer}; kwargs...)

        Applies a [`RarityMetric`](@ref) to a vector of `SDMLayer`s by passing the first layer's valid indices as a mask. 
    """
    function BiodiversityObservationNetworks.evaluate(
        metric::RarityMetric, 
        layers::Vector{<:SDMLayer};
        kwargs...
    )
        rar = BiodiversityObservationNetworks.evaluate(
            metric,
            [l.grid for l in layers];
            mask = layers[1].indices,
            kwargs...
        )
        
        rar_layer = copy(BiodiversityObservationNetworks._FLOAT_TYPE.(first(layers)))
        rar_layer.grid .= rar
        return quantize(rar_layer)
    end
    
    """
        BiodiversityObservationNetworks.rarity(metric::RarityMetric, bon::BiodiversityObservationNetwork, layers::Vector{<:SDMLayer}; kwargs...)

    Applies a [`RarityMetric`](@ref) that requires a [`BiodiversityObservationNetwork`](@ref) to a vector of `SDMLayer`s by passing the first layer's valid indices as a mask. 
    """
    function BiodiversityObservationNetworks.evaluate(
        metric::RarityMetric, 
        layers::Vector{<:SDMLayer},
        bon::BiodiversityObservationNetwork;
        kwargs...
    )
        rar = BiodiversityObservationNetworks.evaluate(
            metric,
            [l.grid for l in layers],
            bon;
            mask = layers[1].indices,
            kwargs...
        )
        
        rar_layer = copy(BiodiversityObservationNetworks._FLOAT_TYPE.(first(layers)))
        rar_layer.grid .= rar
        return quantize(rar_layer)
    end


    # =============================================================================
    # Climate velocity
    # =============================================================================

    function _ols(x, y)
        X = hcat(ones(size(x, 1)), x) 
        XᵀX⁻¹ = inv(X' * X)
        α, β = XᵀX⁻¹ * X' * y
        return α, β
    end
    function temporal_gradient(years, timeseries)
        baseline = first(timeseries)
        temporal_grad = BiodiversityObservationNetworks._FLOAT_TYPE.(copy(baseline))
        for x in eachindex(baseline)
            y = [l[x] for l in timeseries]
            _, β = _ols(years, y)
            temporal_grad[x] = β
        end
        return temporal_grad
    end 

    function spatial_gradient(layer::SDMLayer)
        offset = CartesianIndices((-1:1, -1:1))
        Δx, Δy = -(SpeciesDistributionToolkit.eastings(layer)[[2,1]]...), -(SpeciesDistributionToolkit.northings(layer)[[2,1]]...)
        spatial_grad = BiodiversityObservationNetworks._FLOAT_TYPE.(copy(layer))

        for x in eachindex(layer)
            @inbounds l = layer.grid[x .+ offset]
            @inbounds inc = layer.indices[x .+ offset]

            l[.!(inc)] .= layer.grid[x]
            a,b,c,d,e,f,g,h,i = [l[j,i] for i in 1:3, j in 1:3]

            ∂x = ((c + 2f + i)-(a + 2d + g)) / 8Δx
            ∂y = ((g + 2h + i)-(a + 2b + c)) / 8Δy
            spatial_grad[x] = sqrt((∂x)^2 + (∂y)^2)
        end
        return spatial_grad
    end 

    function BiodiversityObservationNetworks.evaluate(::Loarie2009, years, timeseries::Vector{<:SDMLayer})
        sg = spatial_gradient(timeseries[1])
        tg = temporal_gradient(years, timeseries)
        vel = tg / sg
        return quantize(vel)
    end

    function BiodiversityObservationNetworks.evaluate(metric::Loarie2009, years, timeseries::Vector{<:Vector{<:SDMLayer}}; max_pca_dim = 5)
        # PCA then do velocity on each pair separately
        baseline = timeseries[begin]
        pca = fit(PCA, baseline; maxoutdim=max_pca_dim)
        pca_timeseries = [SpeciesDistributionToolkit.predict(pca, L) for L in timeseries]

        velos = [
            BiodiversityObservationNetworks.evaluate(
                metric,
                years,
                [pca_timeseries[i][j] for i in eachindex(pca_timeseries)]
            ) for j in eachindex(pca_timeseries[begin])
        ]
        quantize(sum(velos))
    end


    
    # =============================================================================
    # Evaluation
    # =============================================================================

    function BiodiversityObservationNetworks._get_inclusion_indicator(raster::SDMLayer, bon)
        mat = zeros(Bool, size(raster))
        mat[bon.sites] .= 1
        return [mat[idx] for idx in findall(raster.indices)]
    end

    """
        evaluate(::MoransI, domain::SDMLayer, bon::BiodiversityObservationNetwork)
    """
    function BiodiversityObservationNetworks.evaluate(
        ::MoransI, 
        domain::SDMLayer,
        bon::BiodiversityObservationNetworks.BiodiversityObservationNetwork
    )
        Es, Ns = eastings(domain), northings(domain)
        coords = hcat([[Es[i], Ns[j]] for i in eachindex(Es), j in eachindex(Ns) if domain.indices[i,j]]...)
        inclusion_indicator = BiodiversityObservationNetworks._get_inclusion_indicator(domain, bon) 

        BiodiversityObservationNetworks._morans_i(coords, inclusion_indicator, bon)
    end

    function BiodiversityObservationNetworks.evaluate(::VoronoiVariance, domain::SDMLayer, bon::BiodiversityObservationNetworks.BiodiversityObservationNetwork)
        bon_coordinates = BiodiversityObservationNetworks._FLOAT_TYPE.(bon.coordinates)
        Es, Ns = eastings(domain), northings(domain)
        domain_coordinates = hcat([[Es[i], Ns[j]] for i in eachindex(Es), j in eachindex(Ns) if domain.indices[i,j]]...)
        return BiodiversityObservationNetworks._voronoi_variance(bon_coordinates, domain_coordinates)    
    end

    function BiodiversityObservationNetworks.evaluate(js::JensenShannon, layers::Vector{<:SpeciesDistributionToolkit.SDMLayer}, bon::BiodiversityObservationNetworks.BiodiversityObservationNetwork)
        Xfull = hcat([[l[i] for l in layers] for i in findall(first(layers).indices)]...)
        Xbon = hcat([[l[i] for l in layers] for i in bon.sites]...)
    
        return BiodiversityObservationNetworks._jensen_shannon(Xfull, Xbon; nbins = js.nbins)
    end
end
