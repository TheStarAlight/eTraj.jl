# Theory: Trajectory Simulation and Phase Methods

*This section reviews the trajectory simulation procedure and the phase methods within.*

Given the initial conditions, the tunneled electrons evolve classically in the combined field of Coulomb and laser, following a classical trajectory, the scheme is named after the *Classical Trajectory Monte Carlo (CTMC)*.
Apart from the position and momentum, phase methods like the *Quantum Trajectory Monte Carlo (QTMC)* and *Semiclassical Two-Step Model (SCTS)* give an additional quantum phase property to the classical trajectories, which are capable of reproducing more details in the final momentum spectrum than the full-classical CTMC.

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

Finally, electrons with similar final momenta would be collected by summing up the probabilities they carry: ``W_{\bm{p}} = \sum_i{W_i}``, and the final momentum spectrum is given by ``W_{\bm{p}}``.



## Quantum Trajectory Monte Carlo (QTMC)

## Semiclassical Two-Step Model (SCTS)
