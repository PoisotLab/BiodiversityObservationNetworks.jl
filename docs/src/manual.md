# BiodiversityObservationNetworks.jl

The purpose of this package is to provide a high-level, extensible, modular interface to the selection of sampling point for biodiversity processes in space. It is based around a collection of types representing point selection algorithms, used to select the most informative sampling points based on raster data. 

!!! warning "This package is in development"
    The `BiodiversityObservationNetworks.jl` package is currently under development. At this point (`v0.4` onwards) the API is not expected to change a lot, but it may change in order to facilitate the integration of new features.

## Installation

### Installing `julia`

Julia can be installed [here](https://julialang.org/downloads/), and is best installed using `juliaup` to manage different versions. On Unix based systems, `juliaup` can be installed with 

```bash
curl -fsSL https://install.julialang.org | sh
```

and on Windows systems with

```sh
winget install julia -s msstore
```

### Installing `BiodiversityObservationNetworks.jl`

`BiodiversityObservationNetworks.jl` is published in the Julia general repository, and can be installed with:

```julia
import Pkg
Pkg.add("BiodiversityObservationNetworks") # [!code highlight]
```

## Package Structure

![Overview of the package](./assets/structure.drawio.svg)


## Manual

The manual follows the [Diataxis](https://diataxis.fr/) schema, where resources fall into four categories:

1. [Tutorials](./tutorials): pedagogy focused
2. [How Tos](./howto): goal focused
3. [Reference](./reference): comprehensive information focused
4. Explanation: focused on understanding the design of the software (primarily for those interested in
   contributing to the package)
