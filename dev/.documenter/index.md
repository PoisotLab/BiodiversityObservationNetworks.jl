
# BiodiversityObservationNetworks.jl

The purpose of this package is to provide a high-level, extensible, modular interface to the selection of sampling point for biodiversity processes in space. It is based around a collection of types representing point selection algorithms, used to select the most informative sampling points based on raster data. Specifically, many algorithms work from a layer indicating _entropy_ of a model based prediction at each location.

::: warning This package is in development

The `BiodiversityObservationNetworks.jl` package is currently under development. The API is not expected to change a lot, but it may change in order to facilitate the integration of new features.

:::

## High-level types {#High-level-types}
<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.BONSampler' href='#BiodiversityObservationNetworks.BONSampler'><span class="jlbinding">BiodiversityObservationNetworks.BONSampler</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
abstract type BONSampler end
```


A `BONSampler` is any algorithm for proposing a set of sampling locations.


[source](https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/6322582496bd794d91f70bf1ce5deae9ccf32e28/src/types.jl#L1-L5)

</details>


::: warning Missing docstring.

Missing docstring for `BONSeeder`. Check Documenter&#39;s build log for details.

:::

::: warning Missing docstring.

Missing docstring for `BONRefiner`. Check Documenter&#39;s build log for details.

:::

## Seeder and refiner functions {#Seeder-and-refiner-functions}

::: warning Missing docstring.

Missing docstring for `seed`. Check Documenter&#39;s build log for details.

:::

::: warning Missing docstring.

Missing docstring for `seed!`. Check Documenter&#39;s build log for details.

:::

::: warning Missing docstring.

Missing docstring for `refine`. Check Documenter&#39;s build log for details.

:::

::: warning Missing docstring.

Missing docstring for `refine!`. Check Documenter&#39;s build log for details.

:::

## Seeder algorithms {#Seeder-algorithms}
<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.BalancedAcceptance' href='#BiodiversityObservationNetworks.BalancedAcceptance'><span class="jlbinding">BiodiversityObservationNetworks.BalancedAcceptance</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
BalancedAcceptance
```


A `BONSeeder` that uses Balanced-Acceptance Sampling (Van-dem-Bates et al. 2017 https://doi.org/10.1111/2041-210X.13003)


[source](https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/6322582496bd794d91f70bf1ce5deae9ccf32e28/src/balancedacceptance.jl#L1-L6)

</details>


## Refiner algorithms {#Refiner-algorithms}

::: warning Missing docstring.

Missing docstring for `AdaptiveSpatial`. Check Documenter&#39;s build log for details.

:::
<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.Uniqueness' href='#BiodiversityObservationNetworks.Uniqueness'><span class="jlbinding">BiodiversityObservationNetworks.Uniqueness</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Uniqueness
```


A `BONSampler`


[source](https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/6322582496bd794d91f70bf1ce5deae9ccf32e28/src/uniqueness.jl#L1-L5)

</details>


## Helper functions {#Helper-functions}

::: warning Missing docstring.

Missing docstring for `squish`. Check Documenter&#39;s build log for details.

:::
<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.entropize!' href='#BiodiversityObservationNetworks.entropize!'><span class="jlbinding">BiodiversityObservationNetworks.entropize!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
entropize!(U::Matrix{AbstractFloat}, A::Matrix{Number})
```


This function turns a matrix `A` (storing measurement values) into pixel-wise entropy values, stored in a matrix `U` (that is previously allocated).

Pixel-wise entropy is determined by measuring the empirical probability of randomly picking a value in the matrix that is either lower or higher than the pixel value. The entropy of both these probabilities are calculated using the -p×log(2,p) formula. The entropy of the pixel is the _sum_ of the two entropies, so that it is close to 1 for values close to the median, and close to 0 for values close to the extreme of the distribution.


[source](https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/6322582496bd794d91f70bf1ce5deae9ccf32e28/src/entropize.jl#L1-L13)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.entropize' href='#BiodiversityObservationNetworks.entropize'><span class="jlbinding">BiodiversityObservationNetworks.entropize</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
entropize(A::Matrix{Number})
```


Allocation version of `entropize!`.


[source](https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/6322582496bd794d91f70bf1ce5deae9ccf32e28/src/entropize.jl#L30-L34)

</details>

