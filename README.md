# Calculation of spin dynamics in itinerant, dipolar ferromagnets
### Some notes and insights on parameters
- standard setting for calculation
--- g = 5.0
--- xi = 1.0
--- qs = np.logspace(-1, 3, 161) (101 steps probably fine as well)
--- j = 10.0
--- rho_upper_limit = 3.0 (exponent --> 10^rho_upper_limit)
--- nrho = 101
--- neta = 21
--- gamma0 = 5.1326
--- iterations = 50

- Physical parameters:
--- g: relative strength of the dipolar interaction to the exchange interaction
----- results for small g, which should correspond to Resibois-Piette function do not look like it
----- results for high g, seem like they could match calculations from Frey and Schwabl
----- for g => 100 (not for g = 60), there is a wiggle feature close to x = 0. Artefact or actual prediction?
--- xi: magnetic correlation length
--- j: exchange interaction strength (does not do anything in the scaled version of the dynamical scaling function)
--- qs: momentum space points for which the dynamical scaling function is calculated

- Calculation parameters:
--- 'gamma0': start value of the dynamical scaling function
----- unless super of, it seems to converge to the same result
--- rho_upper_limit: cut off of the rho integration variable in momentum space (physically from 0 to +ininity) 
----- changes the value for low 'x' quite drastically and not in a ragular monotonous manner
----- keeping number of rho points per exponent intervall constant does not help
--- nrho: number of points for numerical integration over rho
----- extending number of integration points in a does not lead to a convergence to a certain curve, which is convcerning
----- the large x behaviour is almost independent of nrho
----- the small x behaviour is sort of eratic, no monotonous trend
--- neta: number of points for numerical integration over eta
--- iterations: number of iterations in the 'fix point iteration calculation approach'
----- after 50 iterations the change to the curves is marginally

### Credit and Acknowledgement
Thanks goes to Steffen Säubert from the TUM Physik-Department's Chair for the Topology of Correlated Systems who wrote the original code for the Holstein-Primakoff dispersion and the resolution convolution with a SANS resolution function.