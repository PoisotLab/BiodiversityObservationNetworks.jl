const MAX_CORNERPLOT_DIMS_BEFORE_PCA = 20

function BiodiversityObservationNetworks.cornerplot(
    layers::RasterStack;
    pca_layers = false
)

    _, mat = BiodiversityObservationNetworks.features(layers)
    num_layers = length(layers)
    if num_layers > MAX_CORNERPLOT_DIMS_BEFORE_PCA || pca_layers
        pca = BiodiversityObservationNetworks.MultivariateStats.fit(BiodiversityObservationNetworks.MultivariateStats.PCA, mat)
        @info length(pca.prinvars)
        num_layers = length(pca.prinvars)
        mat = BiodiversityObservationNetworks.MultivariateStats.transform(pca, mat)    
    end 

    f = Figure(size=(900,900))
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