using BiodiversityObservationNetworks
using NeutralLandscapes
using Plots

dims = (200, 200)
uncert = rand(MidpointDisplacement(0.95), dims)
uncertplt = heatmap(uncert', frame=:box, axis=:none, colorbartitle="SDM Uncertainty")


npts = 30

weight = 10
noweighting = rand(BalancedAcceptance(numpoints=npts, α=-weight), uncert)
strongweighting = rand(BalancedAcceptance(numpoints=npts, α=weight), uncert)

lowplt = heatmap(uncert',titlefontsize=8, colorbar=:none, title="Weighted toward low uncertainty", size=(350,300), frame=:box, axis=:none, colorbartitle="SDM Uncertainty")
scatter!(lowplt, [x[1] for x in noweighting], [x[2] for x in noweighting], ms=2.5, mc=:white, label="")

hiplt = heatmap(uncert', colorbar=:none, titlefontsize=8, title="Weighted toward high uncertainty", size=(350,300), frame=:box, axis=:none, colorbartitle="SDM Uncertainty")
scatter!(hiplt, [x[1] for x in strongweighting], [x[2] for x in strongweighting], ms=2.5, mc=:white, label="")

plot(uncertplt, lowplt, hiplt, size=(1000, 250), dpi=300, layout=(1,3) )

savefig("balanced_acceptance_with_weights.png")