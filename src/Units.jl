
"""
    module Units

The module `Units` provides some commonly-used units.
"""
module Units

using Reexport: @reexport
@reexport using Unitful: pm, Å, nm, μm
@reexport using Unitful: as, fs, ps, ns
@reexport using Unitful: eV
@reexport using Unitful: W, GW, TW, PW, cm, mm, m
@reexport using Unitful: rad, °
deg = °
export deg

end
