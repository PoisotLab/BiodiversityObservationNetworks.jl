using TestItemRunner


@testmodule TestModule begin
    using SpeciesDistributionToolkit    
    const SDT = SpeciesDistributionToolkit

    export SDT
end

@run_package_tests
