using BiodiversityObservationNetworks
using NeutralLandscapes
using Plots

dims = (100, 100)
uncert = rand(MidpointDisplacement(0.5), dims)
uncertplt = heatmap(uncert, frame=:box, axis=:none, colorbartitle="SDM Uncertainty")


npts = 200

noweighting = rand(BalancedAcceptance(numpoints=npts, α=-100), uncert)
strongweighting = rand(BalancedAcceptance(numpoints=npts, α=100.), uncert)

lowplt = heatmap(uncert,titlefontsize=8, colorbar=:none, title="No weighting", size=(350,300), frame=:box, axis=:none, colorbartitle="SDM Uncertainty")
scatter!(lowplt, [x[1] for x in noweighting], [x[2] for x in noweighting], ms=2.5, mc=:white, label="")

hiplt = heatmap(uncert, colorbar=:none, titlefontsize=8, title="Weighted against low uncertainty", size=(350,300), frame=:box, axis=:none, colorbartitle="SDM Uncertainty")
scatter!(hiplt, [x[1] for x in strongweighting], [x[2] for x in strongweighting], ms=2.5, mc=:white, label="")


plot(uncertplt, lowplt, hiplt, size=(1000, 250), dpi=300, layout=(1,3) )

savefig("balanced_acceptance_with_weights.png")