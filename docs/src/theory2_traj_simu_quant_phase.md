# [Theory: Trajectory Evolution and Quantum Phase](@id theory_traj_phase)

Given the initial conditions, the electrons released from the tunnel exit subsequently evolve classically in the combination of Coulomb and laser fields, following a classical trajectory, and the scheme is named the *Classical Trajectory Monte Carlo (CTMC)*.
In addition to the position and momentum, a quantum phase can be attributed to the evolving trajectories within frameworks such as the *Quantum Trajectory Monte Carlo (QTMC)* and *Semiclassical Two-Step (SCTS) Model*. These models significantly preserve quantum effects in the final PMD, in contrast to the purely classical CTMC.

In this section we review the scheme of trajectory evolution and introduce the quantum phase methods available in `eTraj`.

```@contents
Pages = ["theory2_traj_simu_quant_phase.md"]
Depth = 2
```

---------------------------

## [Classical Trajectory Monte Carlo (CTMC)](@id CTMC)


Within the CTMC framework, each electron, carrying an associated probability ``W``, evolves along a classical trajectory until reaching a final momentum ``\pp_\infty = \pp\rvert_{t=\infty}``, which constitutes the primary observable of interest.

The electrons, emerging from the tunnel with distinct ionization times, initial positions, and momenta, evolve according to Hamilton's equations:
```math
\begin{equation}
    \dot{\rr} = \bm{\nabla}_{\pp} H, \quad \dot{\pp} = - \bm{\nabla}_{\rr} H,
\end{equation}
```
and we use the Hamiltonian under the LG.

After the termination of the laser pulse, the electron interacts solely with the residual parent ion.
At a distance from the parent ion, the electron experiences the Coulomb tail of the potential, and its Runge-Lenz vector
```math
\begin{equation}
    \aa = \pp \times \LL - Z\rr/r
\end{equation}
```
can be considered asymptotically conserved.
By leveraging the conservation of ``\aa`` alongside the angular momentum and energy, we derive the expression for the final momentum [^ShvetsovShilovski_2012]:
```math
\begin{equation}
\begin{aligned}
    \pp_\infty      &= p_\infty \frac{p_\infty (\LL\times\aa)-\aa}{1+p_\infty^2 L^2}, \\
    p_\infty^2/2    &= p^2/2 - Z/r, \\
    \LL             &= \rr\times\pp, \\
    \aa             &= \pp\times\LL - Z\rr/r,
\end{aligned}
\end{equation}
```
where ``\rr`` and ``\pp`` represent the electron's coordinates at any time subsequent to the termination of the laser pulse.
This scheme applies to electrons with positive energy, which are capable of ultimately escaping the parent ion and being detected.
Electrons with negative energy are presumed to be captured into Rydberg orbitals.

Electrons exhibiting comparable final momenta, i.e., those residing within the same bin of the final momentum grid, are aggregated by summing their respective probabilities: ``W_{\pp} = \sum_i W_{i}``. This aggregation results in the final momentum spectrum denoted by ``W_{\pp}``.

[^ShvetsovShilovski_2012]: N. I. Shvetsov-Shilovski, D. Dimitrovski, and L. B. Madsen, Ionization in elliptically polarized pulses: Multielectron polarization effects and asymmetry of photoelectron momentum distributions, *Phys. Rev. A* **85**, 023428 (2012). DOI: [10.1103/PhysRevA.85.023428](https://doi.org/10.1103/PhysRevA.85.023428).

---------------------------

## [Quantum Trajectory Monte Carlo (QTMC)](@id QTMC)

In contrast to the CTMC, the QTMC incorporates a quantum phase into each electron trajectory by employing the Feynman path-integral formalism [^Li_2014]. This quantum phase is equivalent to the action ``S_{\rm{traj}}`` defined in Eq. (7) in [SFA](@ref SFA), with due consideration of the Coulomb potential.
The phase gets accumulated during the electron's excursion and is expressed as
```math
\begin{equation}
    \Phi^{\rm{QTMC}} = S_{\rm{traj}} = - \int_{\tr}^{\infty} \left[ \frac{k^2}{2} + V(\rr) + \Ip \right] \dd t,
\end{equation}
```
where ``\tr`` signifies the instant when the electron exits the tunnel, and ``\kk=\dot{\rr}`` represents the momentum.
Ultimately, the PMD is determined by coherently summing the probability amplitudes corresponding to identical final momenta and taking the squared modulus of the resultant sum:
```math
\begin{equation}
    W_{\pp} = \abs{\sum_i \sqrt{W_{i}} \ee^{\ii \tilde{S}_i}}^2,
\end{equation}
```
where ``\tilde{S}_i`` encompasses the initial phase of the prefactor, the phase accumulated throughout the tunneling process, as well as the trajectory motion:
```math
\begin{equation}
    \tilde{S} = \arg \mathcal{P}_{\pp} + \Re S_{\rm{tun}} + S_{\rm{traj}}.
\end{equation}
```
The phase ``\Re S_{\rm{tun}}`` can be evaluated numerically but may be simplified if we adhere to the non-adiabatic-expansion scheme outlined in [SFA-SPANE](@ref SFAAE):
```math
\begin{equation}
\begin{aligned}
    \Re S_{\rm{tun}}
    &\approx \Re \!\! \int_{\ts}^{\tr} \!\!\!\! \dd\tau \left\{ \frac12 \left[ \pp + \AA(\tr) - \ii\ti\FF(\tr) + \frac12 \ti^2 \FF'(\tr)\right]^2 \!\!\! + \Ip \right\} \\
    &= - [\kk(\tr)\cdot\FF(\tr)] \frac{\ti^2}{2} + o(\ti^2) \\
    &\approx -\kk_0\cdot\rr_0.
\end{aligned}
\end{equation}
```
The last line of Eq. (7), i.e., ``-\kk_0\cdot\rr_0``, vanishes for the SFA-SPANE and ADK initial condition methods due to vanishing longitudinal initial momentum (``k_\parallel = 0``).

It is also noteworthy that in practical implementation, the upper limit of the integral in Eq. (4) need not extend to infinity.
Since electrons arriving at the same final momentum share the same energy after the laser vanishes (at ``\tf``), the integral
```math
\begin{equation}
    \int_{\tf}^{\infty} \left[ \frac{k^2}{2} + V(\rr) + \Ip \right] \dd t
\end{equation}
```
is uniform for electrons with the same final momentum.
Consequently, in numerical implementations, the upper limit of the phase integral in Eq. (4) can be set as the end of the laser pulse, i.e., ``\tf``, leading to
```math
\begin{equation}
    S^{\rm{QTMC}}_{\rm{traj}} = - \int_{\tr}^{\tf} \left[ \frac{k^2}{2} + V(\rr) + \Ip \right] \dd t.
\end{equation}
```

[^Li_2014]: M. Li, J.-W. Geng, H. Liu, Y. Deng, C. Wu, L.-Y. Peng, Q. Gong, and Y. Liu, Classical-quantum correspondence for above-threshold ionization, *Phys. Rev. Lett.* **112**, 113002 (2014). DOI: [10.1103/PhysRevLett.112.113002](https://doi.org/10.1103/PhysRevLett.112.113002)

---------------------------

## [Semiclassical Two-Step (SCTS) Model](@id SCTS)

The SCTS model [^ShvetsovShilovski_2016] improves the quantum phase in the QTMC scheme, giving
```math
\begin{equation}
    \Phi^{\rm{SCTS}} = \Re S_{\rm{tun}} + S_{\rm{traj}}
    = \underbrace{-\kk_0\cdot\rr_0}_{\Re S_{\rm{tun}}} \underbrace{- \int_{t_0}^\infty \left[ \frac{k^2}{2} + V(\rr) - \rr\cdot\bm{\nabla}V(\rr) + \Ip \right] \dd t}_{S_{\rm{traj}}}.
\end{equation}
```
Compared to the QTMC phase given by Eq. (4), the SCTS phase from Eq. (10) differs in two key aspects:
Firstly, there is the initial phase ``-\kk_0\cdot\rr_0``, which results from the tunneling process and is non-zero for non-adiabatic tunneling where ``k_{\parallel} \neq 0``.
Secondly, the integrand includes the ``\rr\cdot\grad V(\rr)`` term, which is absent in the QTMC formulation.
These differences stem from the distinct theoretical foundations of the two methods; the QTMC phase is derived within the framework of first-order perturbation theory, whereas the SCTS formulation extends beyond this approximation.
It should be noted that only the trajectory phase component of the SCTS model, represented by ``S_{\rm{traj}}`` in Eq. (10), is adopted here. The tunneling phase ``\Re S_{\rm{tun}}`` is intended to be incorporated during the preparation of initial conditions.

For the SCTS model, the phase integral in Eq. (10) over the interval ``[\tf, \infty)`` cannot be simply neglected due to the presence of the ``\rr\cdot\bm{\nabla}V(\rr)`` term.
However, the integral of this term can be reduced to an analytical expression in case of Coulomb potential, called the post-pulse Coulomb phase:
```math
\begin{equation}
\begin{aligned}
    S_{\rm{traj,f}}^{\rm{C}}(\tf)
    = \int_{\tf}^{\infty} \rr\cdot\bm{\nabla}V(\rr) \dd t
    = Z \int_{\tf}^{\infty} \frac{\dd t}{r}
    = - n^* \left[ \ln{g} + \sinh^{-1}\left( \frac{\kappa}{g}\rr_{\rm{f}}\cdot\pp_{\rm{f}} \right) \right],
\end{aligned}
\end{equation}
```
where ``\rr_{\rm{f}}=\rr(\tf)``, ``\pp_{\rm{f}}=\pp(\tf)`` and ``g=\sqrt{1+2\kappa^2 L^2}=\sqrt{1+2\kappa^2 (\rr_{\rm{f}}\times\pp_{\rm{f}})^2}``.

In this manner, we derive the expression for the SCTS trajectory phase suitable for numerical implementation:
```math
\begin{equation}
\begin{aligned}
    S^{\rm{SCTS}}_{\rm{traj}}
    = \Ip \tr - \int_{\tr}^{\tf} \left[ \frac{k^2}{2} + V(\rr) - \rr\cdot\bm{\nabla}V(\rr) \right] \dd t + S_{\rm{traj,f}}^{\rm{C}}(\tf).
\end{aligned}
\end{equation}
```

[^ShvetsovShilovski_2016]: N. I. Shvetsov-Shilovski, M. Lein, L. B. Madsen, E. Räsänen, C. Lemell, J. Burgdörfer, D. G. Arbó, and K. Tőkési, Semiclassical two-step model for strong-field ionization, *Phys. Rev. A* **94**, 013415 (2016). DOI: [10.1103/PhysRevA.94.013415](https://dx.doi.org/10.1103/PhysRevA.94.013415)
