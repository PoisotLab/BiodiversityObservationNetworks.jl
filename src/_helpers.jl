function makeenv(;H = 0.8, dims=(100,100)) 
	sdm = rand(MidpointDisplacement(H), dims...)
end

function makesdm(;H = 0.8, dims=(100,100), cutoff= 0.75) 
	sdm = rand(MidpointDisplacement(H), dims...)
	sdm = sdm ./ sum(sdm)
	ind = findall(x->x==1, sdm .<= quantile(vec(sdm), cutoff))
	sdm[ind] .= 0
	sdm = sdm ./ sum(sdm)
	sdm
end

function makeoccurrence(sdm; npts=100)
    occmat = zeros(size(sdm))
    for _ in 1:npts
        ind = rand(Categorical(vec(sdm)))
        occmat[ind] += 1
    end
    occmat
    return findall(i->i>0, occmat)
end

function makebon(sdm, n=50, padding=0.15)
	x,y = size(sdm)
	pts = rand((n,2)...)

	pts[:,1] = floor.((pts[:,1] .* (1-padding)*x) .+ (0.5padding*x))
	pts[:,2] = floor.((pts[:,2] .* (1-padding)*y) .+ (0.5padding*y))
	convert.(Int32,pts)
end