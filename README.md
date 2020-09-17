# Simulating Lennard-Jones (LJ) particles in a periodic boundary condition (PBC), using Kinetic Monte-Carlo (KMC)

The KMC code here is written for simulating
LJ particles in a cubic PBC box. The code uses a normal
distribution to generate new configurations in the PBC box by moving
the LJ particles in x, y, and z directions, and then uses Metropolis
criterion to either accept or reject the new configuration or move.
The code employs the gradient of the energy (force)
to bias the Monte-Carlo algorithm toward a more favorable direction.
