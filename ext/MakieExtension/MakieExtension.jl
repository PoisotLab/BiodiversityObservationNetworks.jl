module MakieExtension
    import Makie
    using BiodiversityObservationNetworks

    function Makie.convert_arguments( P::Makie.PointBased, bon::BiodiversityObservationNetwork) 
        return Makie.convert_arguments(P, [Makie.Point2(n...) for n in eachcol(bon.coordinates)])
    end    


end