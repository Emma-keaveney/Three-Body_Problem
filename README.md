# The Three-Body Problem

This repository looks at finding solutions to the three-body problem using numerical methods in Python and C++. The project is based off a module undertaken during my MPhys at Durham University, in which I modelled the three-body problem in Python.

While no analytical, or closed-form, solution to the three-body problem exists, we can use approximations and numerical methods to model the orbits of the three bodies. The key approximation made by in this case is that one of the bodies is sufficiently light to have a neglible effect on the orbits of the two other bodies. Therefore, the two heavier bodies orbit their centre of mass in circular motion, and we can use numerical methods to plot the orbit of the light body. A three-body problem with this set-up is called the *restricted three-body problem* (https://orbital-mechanics.space/the-n-body-problem/circular-restricted-three-body-problem.html).

## Numerical Methods
The force felt by the light body orbiting the two heavier bodies in the reference frame of the barycentre of the two heavier bodies is,
```math
F = -G m M_1 \frac{r - R_1}{|r- R_1|^{3/2}} - -G m M_2 \frac{r - R_2}{|r- R_2|^{3/2}}.
```
Where $`m`$ is the mass of the light body, $`M_1`$ and $`M_2`$ the masses of the two heavier bodies. $`G`$ is the gravitational constant, and $`r, R_1, R_2`$ the distance of the light body and two heavier bodies respectively to the barycentre.

From this equation for force we yield a second-order ordinary differential equation which describes the acceleration of the light body. As no analytical solution exists for this differential equation we must use numerical methods to solve for the motion of the satellite at each point. Here we use two methods: the *Taylor expansion method* (http://www.met.reading.ac.uk/~sws02hs/teaching/TaylorSeries/TaylorSeriesNotes.pdf), and the *4th order Runge Kutta method*  (https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods). Numerical methods use information about an object at time $`t`$ to calculate its characteristics at time $`t+a`$ where $`a`$ is some small time-step. 

In the Taylor expansion, the position $`x`$ and velocity $`\dot x`$ of the satellite at step $`n+1`$ are found using the position, velocity and acceleration at the previous step,
```math
x_{n+1} = x_n + a \dot x_n + \frac{a^2}{2} \ddot x_n + O(a^3)
```
```math
\dot x_{n+1} = \dot x_n + a \ddot x_n + O(a^2)
```

The Runge Kutta method is more complex as the velocity and acceleration are evaluated at multiple points across each time step and a weighted average used. The fourth-order Runge-Kutta method (RK4) uses four points, with a resulting error term of the order $`O(a^5)`$. RK4 greatly improves the accuracy and stability of the solutions as it uses more points than a Taylor expansion method with the same $`a`$.

## Example 1 - a satellite orbiting the Earth-Moon L2 point

The first example used to simulate the restricted three-body is that of a satellite at the Lunar Lagrangian 2 (L2) point (https://en.wikipedia.org/wiki/Lagrange_point#). Lagrangian points are where the gravitational forces from the two bodies and the centrifugal force balance, producing equilibrium points. The Lunar L2 point lies on the line connecting the Earth and the Moon, positioned behind the Moon, and a satellite at this point will orbit the Earth-Moon barycentre with the same period as the Earth and Moon. An advantage of having a satellite in this location is that the Moon acts as a shield against radio radiation from the Earth, such as man-made radio emissions, or radiation from the Earth’s magnetosphere. This makes the L2 position ideal for satellites detecting low-frequency radio waves from the universe.

The Python and C++ files for simulating this example are 'Three_body.ipynb' and 'three_body.cpp'. The accompanying 'Three_body_visualisation.ipynb' file is a notebook which produces visualisations of the orbits produced in the C++ file.

## Example 2 - galactic collisions and fly-bys

Following on from this simple simulation, we use the code developed to model two galaxies colliding or undergoing close-encounters. Using the method described in Toomre and Toomre's seminal 1972 paper 'Galactic Bridges and Tails' (https://ui.adsabs.harvard.edu/abs/1972ApJ...178..623T/abstract), galactic collisions can be modelled by assuming no self-gravity or interstellar forces within each galaxy, with each star feeling a force from the centre of mass of both galaxies. At each time-step, the new position and velocity of each star is calculated, and then the centre of masses of both galaxies adjusted. We show that, as in Toomre and Toomre, by varying the speed and mass ratios of the two galaxies, the shape and size of the 'tails' and 'bridges' formed during collisions and fly-bys will change. This is of course an very simplified simulation. We only consider infinitely thin galaxies moving on a 2D plane, with no self-gravity, and for only one type of matter. 

We break this example into two stages, in the first we have a galaxy colliding with a point perturber. The Python notebook simple_perturber.ipynb contains the simulation and visualisations in Python. simple_perturber.cpp contains the code for simulating in C++, with the accompanying file simple_perturber_visualisation.ipynb to visualise the results.

Once we have shown how the restricted three-body problem can be used for the collision between a galaxy and a point perturber, we scale this up to have two galaxies colliding. The Python notebook galactic_collisions.ipynb contains the simulation and visualisations in Python. Galactic_collisions.cpp contains the code for simulating in C++, with the accompanying file galactic_collisions_visualisation.ipynb to visualise the results
