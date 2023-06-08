# Initial Conditions

*This section reviews commonly-used theories used to provide initial conditions in the trajectory simulations.*

A number of theories can be adapted to provide initial conditions of the classical electrons in the trajectory simulation scheme.
The initial condition usually consists of three properties:

- Initial position ``\bm{r}_0`` (i.e., the tunneling exit position);
- Initial momentum ``\bm{p}_0``, we note that in the trajectory simulation schemes, initial momentum are usually denoted using symbol ``k``;
- The corresponding ionization probability ``W`` carried by each electron sample, depending on the time-dependent laser field and the properties of the target atoms/molecules.

In the following we will give a brief review on the available theories we implemented in *SemiclassicalSFI.jl*.
Atomic units (a.u.) are used throughout unless stated otherwise.

```@contents
Pages = ["theory_initial_conditions.md"]
```

## Ammosov-Delone-Krainov (ADK)

``\alpha = 1+Z/\sqrt{2I_{\mathrm{p}}}``

## Strong-Field Approximation (SFA)

The Strong-Field Approximation (SFA) is originated from the Keldysh theory of strong-field ionization.
Compared with the pertubative methods and adiabatic tunneling theories, the SFA is able to predict both the multi-photon and the tunneling process during the laser-atom interaction, as well as high-order non-pertubative phenomenona such as the above-threshold ionization (ATI).
The broad scope of SFA has contributed to its widespread application in theoretical investigations of strong-field ionization.

Considering the electron evolving under a combined field of the Coulomb field ``V(\bm{r})`` and the laser field ``\bm{F}(t)=-\partial_t \bm{A}(t)``, under the length gauge (LG), its Hamiltonian reads
```math
H = \frac12 \bm{p}^2 + V(\bm{r}) + \bm{F}(t)\cdot\bm{r}.
```
Denoting ``\mathinner{|\Psi_0\rangle} = \mathinner{|\psi_0\rangle} \mathrm{e}^{\mathrm{i}I_{\mathrm{p}}t}`` as the unperturbed initial state with ionization potential of ``I_{\mathrm{p}}``, ``\mathinner{|\Psi_{\bm{p}}\rangle}`` as the continuum state of momentum ``\bm{p}``, and
```math
U(t_{\mathrm{f}},t_{\mathrm{i}}) = \exp \left[ -\mathrm{i} \int_{t_{\mathrm{i}}}^{t_{\mathrm{f}}} H(\tau) \mathrm{d}\tau \right]
```
the time-evolution operator, the transition amplitude between the initial state and the final state of momentum ``\bm{p}`` is written as
```math
M_{\bm{p}} = \mathinner{\langle \Psi_{\bm{p}} | U(t_{\mathrm{f}},t_{\mathrm{i}}) | \Psi_0 \rangle}.
```

Here lies the key idea of SFA: when the influence of the Coulomb field to the ionized electrons is weak compared with that of the external laser field, we may neglect the influence of the Coulomb field in the expression of ``M_{\bm{p}}`` by replacing the time-evolution operator with a Coulomb-free one ``U_{\mathrm{f}}``, and meanwhile replacing the continuum state with the Volkov state ``\mathinner{| \Psi^{\mathrm{V}}_{\bm{p}} \rangle}`` which represents a free electron evolving under the same laser field:
```math
M_{\bm{p}} = \mathinner{\langle \Psi^{\mathrm{V}}_{\bm{p}} | U_{\mathrm{f}}(t_{\mathrm{f}},t_{\mathrm{i}}) | \Psi_0 \rangle},
```
where the Volkov state under the LG is the product of a plane wave and a phase factor:
```math
\mathinner{| \Psi^{\mathrm{V}}_{\bm{p}} \rangle} = \mathinner{| \bm{p}+\bm{A}(t) \rangle} \mathrm{e}^{-\mathrm{i}S_{\bm{p}}(t)},
```
and the phase has the expression:
```math
S_{\bm{p}}(t) = \int^{t} \frac12 [\bm{p}+\bm{A}(\tau)]^2 \mathrm{d}\tau.
```
In this way the ``M_{\bm{p}}`` is expressed as
```math
M_{\bm{p}} = -\mathrm{i} \int_{t_{\mathrm{i}}}^{t_{\mathrm{f}}} \mathinner{\langle \bm{p}+\bm{A}(\tau) | \bm{F}(\tau)\cdot\bm{r} | \psi_0 \rangle} \mathrm{e}^{\mathrm{i}\tilde{S}_{\bm{p}}(\tau)} \mathrm{d}\tau,
```
and we note that here we have extracted the phase factor of ``\mathinner{|\Psi_0\rangle}`` and combined it with the former ``\mathrm{e}^{\mathrm{i}S_{\bm{p}}(t)}``, giving
```math
\tilde{S}_{\bm{p}}(t) = \int^{t} \left[ \frac12 [\bm{p}+\bm{A}(\tau)]^2 + I_{\mathrm{p}} \right] \mathrm{d}\tau.
```

Utilizing the saddle-point approximation (SPA) would give a more consise expression of ``M_{\bm{p}}``.
The variation of ``\tilde{S}_{\bm{p}}(t)`` is much more sensitive than that of ``\mathinner{\langle \bm{p}+\bm{A}(t) | \bm{F}(t)\cdot\bm{r} | \psi_0 \rangle}`` as ``t`` varies, which leads to a fact that the whole integrand in our latest expression of ``M_{\bm{p}}`` oscillates in its complex phase and its values cancel out each other in most cases, except when the variation of the phase ``\tilde{S}_{\bm{p}}(t)`` becomes stable, i.e., at the saddle points of ``\tilde{S}_{\bm{p}}(t)``. The saddle points ``t_{\mathrm{s}}=t_{\mathrm{r}}+\mathrm{i}t_{\mathrm{i}}`` are the zeroes of the derivative of the complex function ``\tilde{S}_{\bm{p}}(t)``, which satisfy
```math
\partial_t \tilde{S}_{\bm{p}}(t) |_{t=t_{\mathrm{s}}} = \frac12 [\bm{p}+\bm{A}(t_{\mathrm{s}})]^2 + I_{\mathrm{p}} = 0.
```
The integral can be approximated by a summation over the saddle points:
```math
M_{\bm{p}} \approx \sum_{t_{\mathrm{s}}} P_{\bm{p}}(t_{\mathrm{s}}) \mathrm{e}^{\mathrm{i}\tilde{S}_{\bm{p}}(t_{\mathrm{s}})},
```
where ``P_{\bm{p}}(t_{\mathrm{s}})`` denotes the prefactor.
Here we use a modified version of SFA which takes account of the Coulomb potential, which gives the prefactor
```math
P_{\bm{p}}(t_{\mathrm{s}}) = \{ [\bm{p}+\bm{A}(t_{\mathrm{s}})] \cdot \bm{F}(t_{\mathrm{s}}) \}^{-\alpha/2},
```
and we recall that ``\alpha = 1+Z/\sqrt{2I_{\mathrm{p}}}``.
The phase ``\tilde{S}_{\bm{p}}(t_{\mathrm{s}})``, is obtained by solving the integral
```math
\begin{aligned}
    \tilde{S}_{\bm{p}}(t_{\mathrm{s}})
    &= \int_{-\infty}^{t_{\mathrm{s}}} \left[ \frac12 [\bm{p}+\bm{A}(\tau)]^2 + I_{\mathrm{p}} \right] \mathrm{d}\tau \\
    &= \left( -\int_{t_{\mathrm{s}}}^{t_{\mathrm{r}}} -\int_{t_{\mathrm{r}}}^{\infty} \right) \left[ \frac12 [\bm{p}+\bm{A}(\tau)]^2 + I_{\mathrm{p}} \right] \mathrm{d}\tau \\
    &= \Phi_{\mathrm{tun}} + \Phi_{\mathrm{traj}},
\end{aligned}
```
where ``\Phi_{\mathrm{tun}}`` represents the complex phase accumulation during the tunneling process, whose real part denotes the quantum phase, while its imaginary part, is related to the ionization probability; the ``\Phi_{\mathrm{traj}}``, is the phase accumulation during the electron trajectory motion in the continuum.

The SFA provides the final momentum distribution, while the trajectory simulation requires initial conditions of the eletrons.
To utilize the SFA to give initial conditions, we suppose that the classical electron is ejected at time ``t_{\mathrm{r}}`` at tunneling exit ``\bm{r}_0`` with momentum ``\bm{p}_0``.
The initial momentum ``\bm{p}_0``, neglecting the Coulomb potential, is related to the final momentum ``\bm{p}`` through
```math
\bm{p} = \bm{p}_0 + \int_{t_{\mathrm{r}}}^{\infty} \bm{F}(\tau) \mathrm{d}\tau = \bm{p}_0 - \bm{A}(t_{\mathrm{r}}).
```
The initial position ``\bm{r}_0``, i.e., the tunneling exit, is found by constructing a quantum tunneling trajectory.
The beginning of the trajectory, i.e., the tunneling entrance, has a real part of zero; the electron tunnels through the barrier during the time interval ``t_{\mathrm{s}}`` to ``t_{\mathrm{r}}`` and emerges as a classical electron at the tunneling exit with real position and momentum.
In this way we obtain the expression of the initial position:
```math
\bm{r}_0 = \mathrm{Re} \int_{t_{\mathrm{s}}}^{t_{\mathrm{r}}} \bm{A}(\tau) \mathrm{d}\tau.
```
The probablity density (in the final momentum space) carried by the electron sample is
```math
\mathrm{d}W/\mathrm{d}\bm{p} = \lvert P_{\bm{p}}(t_{\mathrm{s}}) \rvert^2 \exp(-2\ \mathrm{Im}\ \Phi_{\mathrm{tun}}).
```

## SFA Adiabatic Expansion (SFA-AE)

## Molecular ADK (MOADK)

## Weak-Field Asymptotic Theory (WFAT)
