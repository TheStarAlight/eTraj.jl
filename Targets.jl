

"""
The Targets module provides information about the targeting atoms or molecules.
"""
module Targets

export Target, SAEAtom, HydrogenLikeAtom
export IonPotential, AsympNuclCharge, TargetPotential, TargetForce, ADKRateExp


abstract type Target end


"Represents an atom under single-active-electron (SAE) approximation."
abstract type SAEAtom <: Target end
# should implement TargetPotential, TargetForce, ADKRateExp.


begin :HydrogenLikeAtom
    "Represents a Hydrogen-like atom."
    struct HydrogenLikeAtom <: SAEAtom
        "Ionization potential of the atom."
        IonPotential;
        "Asymptotic charge of the inner nucleus."
        NuclCharge;
    end
    "Gets the ionization potential of the atom."
    IonPotential(t::HydrogenLikeAtom) = t.IonPotential
    "Gets the asymptotic nuclear charge of the atom."
    AsympNuclCharge(t::HydrogenLikeAtom) = t.NuclCharge
    "Gets the potential function of the atom."
    TargetPotential(t::HydrogenLikeAtom) = (x,y,z) -> -t.NuclCharge*(x^2+y^2+z^2+1.0)^(-0.5)
    "Gets the force exerted on the electron from the atom (which is the neg-grad of potential)."
    TargetForce(t::HydrogenLikeAtom) = (x,y,z) -> -t.NuclCharge*(x^2+y^2+z^2+1.0)^(-1.5) .* (x,y,z)
    """
    Gets the exponential term of ADK rate which depends on
    Field strength `F`,
    Azimuthal angle of field `φ`,
    momentum's transverse component `pd` (in xy plane),
    and propagation-direction (which is Z axis) component `pz`.
    """
    ADKRateExp(t::HydrogenLikeAtom) = (F,φ,pd,pz) -> exp(-2(pd^2+pz^2+2*t.IonPotential)^1.5/3F)
end

end