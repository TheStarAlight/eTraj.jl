
"""
**SemiclassicalSFI.jl**
Implementation of classical/semiclassical trajectory methods in strong-field ionization of atoms and molecules.
"""
module SemiclassicalSFI

include("Lasers/Lasers.jl")
include("Targets/Targets.jl")
include("ElectronSamplers/ElectronSamplers.jl")
using .Lasers
using .Targets
using .ElectronSamplers
include("TrajectorySimulation.jl")

export perform_traj_simulation, Lasers, Targets

using Pkg
function get_SCSFI_version()
    for (k,v::Pkg.API.PackageInfo) in Pkg.dependencies()
        if v.name == "SemiclassicalSFI"
            return v.version
        end
    end
end

end
