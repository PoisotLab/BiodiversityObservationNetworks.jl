module SDTExtension
    using SpeciesDistributionToolkit
    using BiodiversityObservationNetworks
    using TestItems


    _apply_sdm_mask!(::BitMatrix, ::SDMLayer, ::Missing) = nothing
    function _apply_sdm_mask!(valid::BitMatrix, layer::SDMLayer, mask::AbstractMatrix)
        size(mask) == size(layer) || throw(ArgumentError("Mask size $(size(mask)) must match layer size $(size(layer))"))
        valid .&= Bool.(mask)
        return
    end
    function _apply_sdm_mask!(valid::BitMatrix, layer::SDMLayer, mask::SDMLayer)
        SimpleSDMLayers._layers_are_compatible(layer, mask) || throw(ArgumentError("Layer and mask are not compatible. Ensure size, extent, and CRS match."))
        valid .&= mask.indices
        return
    end

    function _sdm_coordinates(layer::SDMLayer, keys::Vector{CartesianIndex{2}})
        Es = SimpleSDMLayers.eastings(layer)
        Ns = SimpleSDMLayers.northings(layer)
        coords = hcat([[Es[k[2]], Ns[k[1]] ] for k in keys]...)
        return coords
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
        first_layer = first(layers)
        
        all(l -> SimpleSDMLayers._layers_are_compatible(l, first_layer), layers) ||
            throw(ArgumentError("All SDMLayers must have the same size, extent, and CRS"))

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
    - `resolution`: grid cell size in decimal degrees (default `0.5°`)
    - `mask`: optional additional `BitMatrix` or `SDMLayer` mask applied after polygon rasterization
    - `inclusion`: optional per-cell weight `Matrix` or `Vector`
    """
    function BiodiversityObservationNetworks.CandidatePool(
        poly::SimpleSDMPolygons.AbstractGeometry;
        resolution = 0.5,
        mask = missing,
        inclusion = missing,
    )
        bbox = SimpleSDMPolygons.boundingbox(poly)
        nrows = max(1, round(Int, (bbox.top - bbox.bottom) / resolution))
        ncols = max(1, round(Int, (bbox.right - bbox.left) / resolution))

        layer = SDMLayer(
            ones(nrows, ncols);
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
        cp = CandidatePool(poly; resolution = 1.0)
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
end
