# BiodiversityObservationNetworks.jl

The purpose of this package is to provide a high-level, extensible, modular
interface to the selection of sampling point for biodiversity processes in
space. It is based around a collection of types representing point selection
algorithms, used to select the most informative sampling points based on raster
data. Specifically, the algorithms work from a layer indicating *entropy* at
each location.

!!! warning "This package is in development"
    The `BiodiversityObservationNetworks.jl` package is currently under development. The API is not expected to change a lot, but it may change in order to facilitate the integration of new features.

## High-level types

```@docs
BONSampler
BONSeeder
BONRefiner
```

## Seeder and refiner functions

```@docs
seed
seed!
```

```@docs
refine
refine!
```

## Seeder algorithms

```@docs
BalancedAcceptance
```

## Refiner algorithms

```@docs
AdaptiveSpatial
```