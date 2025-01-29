
# How to measure the distance between a BON and the whole env space in layer
# stack?
# Marginals instead of MvNormal bc everything gets fucky 
# But can penalize w/ weighted distance between Covariance matrices

function _js_thing(rs::RasterStack, bon::BiodiversityObservationNetwork)
    function _js(P,Q) 
        M = BONs.Distributions.MixtureModel([P,Q], [0.5,0.5])
        div = 0.5*BONs.Distributions.kldivergence(P,M) + 0.5*BONs.Distributions.kldivergence(Q,M)
        return sqrt(div / log(2))    
    end
    _, Xfull = BONs.features(rs)
    Xsampled = rs[bon]

    nlayers = length(rs)
    
    Σ₁, Σ₂ = zeros(nlayers, nlayers), zeros(nlayers, nlayers)

    𝓛_js = 0.
    for i in axes(Xfull,1)
        𝓝₁ = BONs.Distributions.fit(BONs.Distributions.Normal, Xfull[i,:])
        𝓝₂ = BONs.Distributions.fit(BONs.Distributions.Normal, Xsampled[i,:])
        𝓛_js += _js(𝓝₁, 𝓝₂) 
        for j in i+1:nlayers
            Σ₁[i,j] = BONs.Distributions.cov(Xfull[:,i],Xfull[:,j])
            Σ₁[j,i] = Σ₁[i,j]
            # How to measure the distance between a BON and the whole env space in layer
            # stack?
            # Marginals instead of MvNormal bc everything gets fucky 
            # But can penalize w/ weighted distance between Covariance matrices
            
            function _js_thing(rs::RasterStack, bon::BiodiversityObservationNetwork)
                function _js(P,Q) 
                    M = BONs.Distributions.MixtureModel([P,Q], [0.5,0.5])
                    div = 0.5*BONs.Distributions.kldivergence(P,M) + 0.5*BONs.Distributions.kldivergence(Q,M)
                    return sqrt(div / log(2))    
                end
                _, Xfull = BONs.features(rs)
                Xsampled = rs[bon]
            
                nlayers = length(rs)
                
                Σ₁, Σ₂ = zeros(nlayers, nlayers), zeros(nlayers, nlayers)
            
                𝓛_js = 0.
                for i in axes(Xfull,1)
                    𝓝₁ = BONs.Distributions.fit(BONs.Distributions.Normal, Xfull[i,:])
                    𝓝₂ = BONs.Distributions.fit(BONs.Distributions.Normal, Xsampled[i,:])
                    𝓛_js += _js(𝓝₁, 𝓝₂) 
                    for j in i+1:nlayers
                        Σ₁[i,j] = BONs.Distributions.cov(Xfull[:,i],Xfull[:,j])
                        Σ₁[j,i] = Σ₁[i,j]
            
                        Σ₂[i,j] = BONs.Distributions.cov(Xsampled[:,i],Xsampled[:,j])
                        Σ₂[j,i] = Σ₂[i,j]
                    end 
                end
                𝓛_covariance = sqrt(sum((Σ₁ .- Σ₂).^2))
            
                𝓛_js, 𝓛_covariance
            end 
            
            
            n_reps = 500
            n_nodes = [2^i for i in 4:9]
            begin 
                f = Figure()
                ax = Axis(f[1,1])
                for n in n_nodes
                    density!(ax, [_js_thing(bioclim, sample(SimpleRandom(n), bioclim))[1] for _ in 1:n_reps])
                end 
                f
            end 

            Σ₂[i,j] = BONs.Distributions.cov(Xsampled[:,i],Xsampled[:,j])
            Σ₂[j,i] = Σ₂[i,j]
        end 
    end
    𝓛_covariance = sqrt(sum((Σ₁ .- Σ₂).^2))

    𝓛_js, 𝓛_covariance
end 


n_reps = 500
n_nodes = [2^i for i in 4:9]
begin 
    f = Figure()
    ax = Axis(f[1,1])
    for n in n_nodes
        density!(ax, [_js_thing(bioclim, sample(SimpleRandom(n), bioclim))[1] for _ in 1:n_reps])
    end 
    f
end 