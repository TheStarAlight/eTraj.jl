# Theory: Trajectory Simulation and Phase Methods

*This section reviews the trajectory simulation procedure and the phase methods within.*

Given the initial conditions, the tunneled electrons evolve classically in the combined field of Coulomb and laser, following a classical trajectory, the scheme is named after the *Classical Trajectory Monte Carlo (CTMC)*.
Apart from the position and momentum, phase methods like the *Quantum Trajectory Monte Carlo (QTMC)* and *Semiclassical Two-Step (SCTS) Model* give an additional quantum phase property to the classical trajectories, which are capable of reproducing more details in the final momentum spectrum than the full-classical CTMC.

In the following we review the scheme of trajectory simulation and introduce the quantum phase methods available at present.
Note that atomic units (a.u.) is used throughout.

```@contents
Pages = ["theory2_trajectory_simulation_phase_methods.md"]
```


## Classical Trajectory Monte Carlo (CTMC)

In the CTMC, each sample electron carries a probability ``W``, following a classical trajectory, and finally ends up with a final momentum ``\bm{p}_\infty = \bm{p}|_{t=\infty}``, which is our interested physical quantity.

The tunneled electrons, each having different tunneling time, initial positions and momenta, evolve under the Newtonian equation of motion:
```math
\ddot{\bm{r}} = - \bm{F}(t) - \bm{\nabla}V(\bm{r}),
```
where ``V(\bm{r})`` denotes the potential of the parent ion.

After the laser ends, the electron interacts only with the residual parent ion.
At a distance from the parent ion, the electron interacts with the potential's Coulomb tail, and its Runge-Lenz vector ``\bm{a} = \bm{p}\times\bm{L} - Z\bm{r}/r`` can be viewed as approximately conserved. Taking advantage of this conserved quantity, combining with the conservation of angular momentum and energy, we obtain the expression of the final momentum:
```math
\begin{aligned}
    \bm{p}_\infty &= p_\infty \frac{p_\infty(\bm{L}\times\bm{a})-\bm{a}}{1+p_\infty^2 L^2}, \\
    p_\infty^2/2 &= p^2/2 - Z/r, \\
    \bm{L} &= \bm{r}\times\bm{p}, \\
    \bm{a} &= \bm{p}\times\bm{L} - Z\bm{r}/r, \\
\end{aligned}
```
where ``\bm{r},\bm{p}`` are quantities of the electron at any time after the laser ends.
This scheme applies for electrons with positive energy, which are able to finally escape the parent ion and reach the detector.
For electrons with negative energy, we assume that they finally become rydberg states.

Finally, electrons with similar final momenta (i.e., in the same small box of the final momentum grid) would be collected by summing up the probabilities they carry: ``W_{\bm{p}} = \sum_i{W_i}``, and the final momentum spectrum is given by ``W_{\bm{p}}``.


## Quantum Trajectory Monte Carlo (QTMC)

Compared with the CTMC, the QTMC scheme endows each electron trajectory with a quantum phase ``\Phi`` based on the Feynman path-integral approach.
The phase gets acculmulated during the electron's excursion and is expressed as
```math
\Phi = - \int_{t_0}^\infty \left[ \frac{p^2}{2} + V(\bm{r}) + I_{\mathrm{p}} \right] \mathrm{d}t.
```
where ``t_0`` is the time when the electron tunneled.
Finally the momentum spectrum is given by coherently summing up the probability amplitude, and taking the square modulus of the summation result:
```math
W_{\bm{p}} = \left| \sum_i \sqrt{W_i}\ \mathrm{e}^{\mathrm{i}\Phi_i} \right|^2.
```

It's also worthwhile noting that in practical implementation, the upper limit of the integral of the quantum phase ``\Phi`` doesn't have to be infinity.
Since electrons which arrived at the same final momentum share the same energy after the laser ends (at ``t_{\mathrm{f}}``), the integral
```math
\int_{t_{\mathrm{f}}}^\infty \left[ \frac{p^2}{2} + V(\bm{r}) + I_{\mathrm{p}} \right] \mathrm{d}t
```
is same for electrons with the same final momentum.
Therefore, in numerical implementation, the upper limit of the phase integral can be simply set as the end of the laser, i.e., the ``t_{\mathrm{f}}``.

## Semiclassical Two-Step (SCTS) Model

The SCTS model improves the quantum phase in the QTMC scheme, giving
```math
\Phi = - \bm{k}_0\cdot\bm{r}_0 - \int_{t_0}^\infty \left[ \frac{p^2}{2} + V(\bm{r}) - \bm{r}\cdot\bm{\nabla}V(\bm{r}) + I_{\mathrm{p}} \right] \mathrm{d}t.
```
The difference between the SCTS phase and the QTMC lies in two aspects:
The first is the initial phase ``\bm{k}_0\cdot\bm{r}_0``, which is non-zero for non-zero initial longitutinal momentum ``k_\parallel`` w.r.t. the non-adiabatic tunneling process.
The second is the ``\bm{r}\cdot\bm{\nabla}V(\bm{r})`` term in the integrand which is omitted in the QTMC scheme.

For the SCTS model, the phase integral in the interval ``[t_{\mathrm{f}},\infty)`` cannot be simply neglected due to the presence of the ``\bm{r}\cdot\bm{\nabla}V(\bm{r})`` term in the integrand.
However, the integral of this term can be reduced to an analytical expression in case of Coulomb potential (``V(r)=Z/r``):
```math
\begin{aligned}
    \Phi_{\mathrm{f}}^{\mathrm{C}}(t_{\mathrm{f}})
    &= \int_{t_{\mathrm{f}}}^\infty \bm{r}\cdot\bm{\nabla}V(\bm{r}) \mathrm{d}t \\
    &= Z \int_{t_{\mathrm{f}}}^\infty \frac{\mathrm{d}t}{r} \\
    &= - \frac{Z}{\kappa} \left[ \ln{g} + \sinh^{-1} \left( \frac{\kappa}{g}\bm{r}_{\mathrm{f}}\cdot\bm{p}_{\mathrm{f}} \right) \right],
\end{aligned}
```
where ``\bm{r}_{\mathrm{f}}=\bm{r}(t_{\mathrm{f}})``, ``\bm{p}_{\mathrm{f}}=\bm{p}(t_{\mathrm{f}})`` and ``g = \sqrt{1+2\kappa^2 L^2} = \sqrt{1+2\kappa^2 (\bm{r}_{\mathrm{f}}\times\bm{p}_{\mathrm{f}})^2}``.
In this way we obtain the expression of the SCTS phase that is suitable for numerical implementation:
```math
\Phi = - \bm{k}_0\cdot\bm{r}_0 + I_{\mathrm{p}}t_0 - \int_{t_0}^{t_{\mathrm{f}}} \left[ \frac{p^2}{2} + V(\bm{r}) - \bm{r}\cdot\bm{\nabla}V(\bm{r}) \right] \mathrm{d}t + \Phi_{\mathrm{f}}^{\mathrm{C}}(t_{\mathrm{f}}).
```
