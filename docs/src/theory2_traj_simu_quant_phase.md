# [Theory: Trajectory Simulation and Quantum Phase](@id theory_traj_phase)

Given the initial conditions, the tunneled electrons evolve classically in the combination of Coulomb and laser fields, following a classical trajectory, the scheme is named after the *Classical Trajectory Monte Carlo (CTMC)*.
Apart from the position and momentum, phase methods like the *Quantum Trajectory Monte Carlo (QTMC)* and *Semiclassical Two-Step (SCTS) Model* endow an additional quantum phase property to the classical trajectories, which are capable of reproducing more details in the final momentum spectrum than the full-classical CTMC.

In this section we review the scheme of trajectory simulation and introduce the quantum phase methods available in `eTraj`.

```@contents
Pages = ["theory2_traj_simu_quant_phase.md"]
Depth = 2
```

---------------------------

## [Classical Trajectory Monte Carlo (CTMC)](@id CTMC)

In the CTMC, each electron carries a probability ``W``, following a classical trajectory, and finally ends up with a final momentum ``\pp_\infty = \pp\rvert_{t=\infty}``, which is our interested physical quantity.

The tunneled electrons, each having different tunneling times, initial positions and momenta, evolve under the Hamiltonian equation:
```math
\begin{equation}
    \dot{\rr} = \bm{\nabla}_{\pp} H, \quad \dot{\pp} = - \bm{\nabla}_{\rr} H,
\end{equation}
```
and we use the Hamiltonian under the LG.

After the laser ends, the electron interacts only with the residual parent ion.
At a distance from the parent ion, the electron interacts with the potential's Coulomb tail, and its Runge-Lenz vector
```math
\begin{equation}
    \aa = \pp \times \LL - Z\rr/r
\end{equation}
```
can be viewed as asymptotically conserved.
Taking advantage of the conservation of ``\aa`` as well as the angular momentum and energy, we obtain the expression of the final momentum [^ShvetsovShilovski_2012]:
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
where ``\rr`` and ``\pp`` are quantities of the electron at any time after the laser ends.
This scheme applies for electrons with positive energy, which are able to finally escape the parent ion and reach the detector.
For electrons with negative energy, we assume that they finally fall on Rydberg states.

Finally, electrons with similar final momenta (i.e., in the same box of the final momentum grid) would be collected by summing up the probability they carry: ``W_{\pp} = \sum_i W_{i}``, and the final momentum spectrum is given by ``W_{\pp}``.

[^ShvetsovShilovski_2012]: N. I. Shvetsov-Shilovski, D. Dimitrovski, and L. B. Madsen, Ionization in elliptically polarized pulses: Multielectron polarization effects and asymmetry of photoelectron momentum distributions, *Phys. Rev. A* **85**, 023428 (2012). DOI: [10.1103/PhysRevA.85.023428](https://doi.org/10.1103/PhysRevA.85.023428).

---------------------------

## [Quantum Trajectory Monte Carlo (QTMC)](@id QTMC)

Compared with the CTMC, the QTMC scheme endows each electron trajectory with a quantum phase based on the Feynman path-integral approach [^Li_2014], which is actually the ``S_{\rm{traj}}`` in Eq. (7) in [SFA](@ref SFA) with additional account of the Coulomb potential.
The phase gets accumulated during the electron's excursion and is expressed as
```math
\begin{equation}
    \Phi^{\rm{QTMC}} = S_{\rm{traj}} = - \int_{\tr}^{\infty} \left[ \frac{k^2}{2} + V(\rr) + \Ip \right] \dd t,
\end{equation}
```
where ``\tr`` is the time when the electron tunnels, and ``\kk=\dot{\rr}`` denotes the momentum.
Finally, the momentum spectrum is given by coherently summing up the probability amplitude, and taking the square modulus of the summation result:
```math
\begin{equation}
    W_{\pp} = \abs{\sum_i \sqrt{W_{i}} \ee^{\ii \tilde{S}_i}}^2,
\end{equation}
```
where ``\tilde{S}_i`` contains the initial phase of the prefactor and the phase accumulated during the tunneling process, as well as the trajectory motion:
```math
\begin{equation}
    \tilde{S} = \arg \mathcal{P}_{\pp} + \Re S_{\rm{tun}} + S_{\rm{traj}}.
\end{equation}
```
The ``\Re S_{\rm{tun}}`` can be evaluated numerically, but can be simplified if we follow the non-adiabatic-expansion scheme in [SFA-SPANE](@ref SFAAE):
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
The last line of Eq. (7), i.e., ``-\kk_0\cdot\rr_0``, vanishes for the SFA-SPANE and ADK initial condition methods due to zero longitudinal initial momentum (``k_\parallel = 0``).

It's also worthwhile noting that in practical implementation, the upper limit of the integral in Eq. (4) doesn't have to be infinity.
Since electrons arrive at the same final momentum share the same energy after the laser ends (at ``\tf``), the integral
```math
\begin{equation}
    \int_{\tf}^{\infty} \left[ \frac{k^2}{2} + V(\rr) + \Ip \right] \dd t
\end{equation}
```
is same for electrons with the same final momentum.
Therefore, in numerical implementation, the upper limit of the phase integral in Eq. (4) can be simply set as the end of the laser, i.e., the ``\tf``, giving
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
The difference between the SCTS phase [``\Phi^{\rm{SCTS}}`` in Eq. (10)] and the QTMC phase [``\Phi^{\rm{QTMC}}`` in Eq. (4)] lies in two aspects:
The first is the initial phase ``-\kk_0\cdot\rr_0`` obtained in the tunneling process, which is non-zero for ``k_{\parallel}\neq 0`` w.r.t. the non-adiabatic tunneling process.
The second is the ``\rr\cdot\bm{\nabla}V(\rr)`` term in the integrand, which is omitted in the QTMC scheme.
The difference arises from formulation schemes of the two methods, as the QTMC obtained the phase under the first-order perturbation theory, while the SCTS's formulation went beyond the perturbation theory.
We note that we would only adopt the trajectory phase of the SCTS model, i.e., the ``S_{\rm{traj}}`` in Eq. (10), because the tunneling phase ``\Re S_{\rm{tun}}`` is supposed to be determined by the initial condition methods.

For the SCTS model, the phase integral in Eq. (10) in the interval ``[\tf, \infty)`` cannot be simply neglected due to the presence of the ``\rr\cdot\bm{\nabla}V(\rr)`` term.
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

In this way we obtain the expression of the SCTS trajectory phase that is suitable for numerical implementation:
```math
\begin{equation}
\begin{aligned}
    S^{\rm{SCTS}}_{\rm{traj}}
    = \Ip \tr - \int_{\tr}^{\tf} \left[ \frac{k^2}{2} + V(\rr) - \rr\cdot\bm{\nabla}V(\rr) \right] \dd t + S_{\rm{traj,f}}^{\rm{C}}(\tf).
\end{aligned}
\end{equation}
```

[^ShvetsovShilovski_2016]: N. I. Shvetsov-Shilovski, M. Lein, L. B. Madsen, E. Räsänen, C. Lemell, J. Burgdörfer, D. G. Arbó, and K. Tőkési, Semiclassical two-step model for strong-field ionization, *Phys. Rev. A* **94**, 013415 (2016). DOI: [10.1103/PhysRevA.94.013415](https://dx.doi.org/10.1103/PhysRevA.94.013415)
