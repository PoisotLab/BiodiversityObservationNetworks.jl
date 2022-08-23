using Test
using BiodiversityObservationNetworks

tests = [
    "\033[1m\033[34mSEEDER\033[0m  Balance Acceptance Sampling" => "balancedacceptance.jl",
    "\033[1m\033[35mREFINER\033[0m Adaptive Spatial Sampling" => "adaptivespatialsampling.jl",
]
 
global anyerrors = false
for test in tests
    try
       include(test.second)
       println("\033[1m\033[32m✓\033[0m\t$(test.first)")
    catch e
       global anyerrors = true
       println("\033[1m\033[31m×\033[0m\t$(test.first)")
       println("\033[1m\033[38m→\033[0m\ttest/$(test.second)")
       showerror(stdout, e, backtrace())
       println()
       break
    end
 end
 
 if anyerrors
    throw("Tests failed")
 end

