using BiodiversityObservationNetworks
using NeutralLandscapes
using Plots

dims = (250,250)
uncert = rand(MidpointDisplacement(0.99), dims)
uncertplt = heatmap(uncert, frame=:box, axis=:none, colorbartitle="SDM Uncertainty")


npts = 100

noweighting = rand(BalancedAcceptance(numpoints=npts, α=0.0001), uncert)
strongweighting = rand(BalancedAcceptance(numpoints=npts, α=1000.), uncert)

lowplt = heatmap(uncert, size=(350,300), frame=:box, axis=:none, colorbartitle="SDM Uncertainty")
scatter!(lowplt, [x[1] for x in noweighting], [x[2] for x in noweighting], ms=2.5, mc=:white, label="")

hiplt = heatmap(uncert, size=(350,300), frame=:box, axis=:none, colorbartitle="SDM Uncertainty")
scatter!(hiplt, [x[1] for x in strongweighting], [x[2] for x in strongweighting], ms=2.5, mc=:white, label="")


plot(uncertplt, lowplt, hiplt, size=(1000, 250), dpi=300, layout=(1,3) )