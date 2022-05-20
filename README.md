# Main programs

## jacobi.f90

This program is primarily designed to compare numerical integration of a finite region using either mixed or single Jacobi coordinates. The region is defined in mixed Jacobi coordinates, and the integrands for the two cases are hardcoded to correspond to the same function. The most direct comparison can be seen through the integrate_triple_Simpson function, which uses Simpson's rule to evaluate the triple integral in both coordinate systems.

Additionally, a separated mixed Jacobi integral is calculated and expected to agree with the unseparated triple integral. Monte Carlo integration can also used to evaluate the single Jacobi integral, and a further, unrelated Cartesian integral can be performed to test the triple integration with variable limits.

Unless otherwise stated:
* (R_a, r_a and gamma_a) are the single Jacobi coordinates for channel A
* (R_a, R_b and gamma_ab) are mixed Jacobi coordinates
* mass_a, mass_b and mass_c are the masses of atoms A, B and C
* mass_total is the sum of the masses of A, B and C
* mu_a and mu_b are reduced channel masses
* m_a and m_b are internal reduced masses


### module precisn

This defines the precision of real numbers throughout the program, currently such that real numbers are stored up to 12 decimal places.


### subroutine main

This sets up the main integration parameters used throughout the program and calls the integration functions, printing their results and the calculation time. When all code is uncommented, this leads to the evaluation of:

1. Separated mixed Jacobi integral using Simpson's rule
2. Unseparated mixed Jacobi integral using Simpson's rule
3. Unseparated single Jacobi integral using Simpson's rule
4. Unseparated single Jacobi integral using Monte Carlo integration
5. Unseparated Cartesian test integral using Simpson's rule

Note:
* (1), (2), (3), and (4) share the same integration limits
* (1) and (2) share the same integrand: mixed_jacobi_integrand_func
* (3) and (4) share the same integrand: single_jacobi_integrand_func
* single_jacobi_integrand_func and mixed_jacobi_integrand_func are hardcoded to represent the same integrand

Options:
* save_inputs: If true, saves the values of masses and Jacobi integration limits in an unformatted file (outputs/values)
* mixed_jacobi_n: Number of Simpson cells used for mixed Jacobi integration using Simpson's rule
* single_jacobi_n: Number of Simpson cells used for single Jacobi integration using Simpson's rule
* test_coord_n: Number of Simpson cells used for the test integration using Simpson's rule
* mc_n: Number of points to generate for Monte Carlo integration
* mixed_jacobi_lims: Integration limits used for mixed and single Jacobi integration, in terms of mixed Jacobi coordinates
* test_coord_lims: Integration limits used for the test integration (may be used as placeholders for variable limits)


### subroutine calculate_masses

This returns hardcoded definitions of three masses, as well as their sum, internal reduced masses, and channel reduced masses.

Inputs:
* mass_a (undefined)
* mass_b (undefined)
* mass_c (undefined)
* mass_total (undefined)
* mu_a (undefined)
* mu_b (undefined)
* m_a (undefined)
* m_b (undefined)

Outputs:
* mass_a
* mass_b
* mass_c
* mass_total
* mu_a
* mu_b
* m_a
* m_b


### function integrate_single_Simpson

This uses Simpson's rule to integrate the function defined in "mixed_jacobi_integrand_func". Currently, this allows the integration of x^4 or sin(x)cos^2(x).

Inputs:
* n: Number of Simpson intervals, which must be positive and even
* a: Lower limit of integration
* b: Upper limit of integration
* is_x_squared: Flag to determine whether to integrate x^4 (true) or sin(x)cos^2(x) (false)

Output:
* integral: Value of calculated integral

Errors:
* If n is not positive and even, the program will stop


### function mixed_jacobi_integrand_func

This defines and evaluates an integrand for a single variable, x. Currently, this evaluates either x^4 or sin(x)cos^2(x), based on a flag that is passed. Mixed Jacobi coordinates have a volume element R_a^2 R_b^2 sin(gamma_ab), so this currently allows for integration of the function R_a^2 R_b^2 cos^2(gamma_ab), one variable at a time.

Inputs:
* x: The current value of the variable in the integrand to be evaluated, i.e. R_a, R_b or gamma_ab
* is_x_squared: Flag to determine whether to integrate x^4 (true) or sin(x)cos^2(x) (false)

Output:
* integrand: The evaluated integrand


### function single_jacobi_integrand_func

This defines and evaluates an integrand for three variables. Currently, this evaluates the same function as in "mixed_jacobi_integrand_func", but for a different volume element: R_a^2 r_a^2 sin(gamma_a).

Inputs:
* x: The current value of the first variable in the integrand to be evaluated, R_a
* y: The current value of the second variable in the integrand to be evaluated, r_a
* z: The current value of the third variable in the integrand to be evaluated, gamma_a
* mu_a
* mass_c
* m_b

Output:
* integrand: The evaluated integrand


### function calc_r_a

This calculates the value of r_a, given known mixed Jacobi coordinates (R_a, R_b and gamma_ab).

Inputs:
* R_a
* R_b
* gamma_ab
* mu_a
* m_b
* mass_c

Output:
* coord: r_a


### function calc_gamma_a

This calculate gamma_a from known values of R_a, r_a and gamma_ab. Not all values of R_a, r_a and gamma_ab are consistent with the allowed range of R_b, so R_b is calculated first, and a placeholder value (-100._wp) is returned if it is not valid.

Inputs:
* small_r_a: r_a
* mu_a
* R_a
* mass_c
* gamma_ab
* m_b
* lims: The mixed Jacobi integration limits (R_a_min, R_a_max, R_b_min, R_b_max, gamma_ab_min, gamma_ab_max)
* rand_root: Flag to determine whether to use a random number generator to select either the positive or negative root when calculating R_b
* pos_root: Flag to determine whether to use the positive root when calculating R_b, if rand_root is false

Options:
* gamma_tol: A tolerance for the value of cos(gamma_a), which may be slightly outside the valid range (-1 - 1) due to numerical errors
* R_b_tol: A tolerance for the value of R_b, which may be slightly outside its defined range (R_b_min - R_b_max) due to numerical errors

Output:
* gamma_a

Errors:
* If cos(gamma_a) is outside its valid range, including the defined tolerance, the program will stop. This may be because the tolerance is too small, or there may be an error with the inputs


### function integrate_triple_Simpson

This uses Simpson's rule to integrate three variables. This is designed to be run in three different ways:

* If the mixed_jacobi flag is set to true, an integrand defined in "mixed_jacobi_integrand_func" for each coordinate is evaluated between the limits defined by mixed_lims

* If the single_jacobi flag is set to true, an integrand defined in "single_jacobi_integrand_func" is evaluated for all three coordinates together. This currently evaluates the same function, over the same volume, as when mixed_jacobi is true, through variable limits and a scaling factor

* If test_coords is set to true, an integrand defined by "test_integrand_func" for all three variables is evaluated for all three coordinates together. This integrand is unrelated to the mixed_jacobi and single_jacobi integrands, but as for "single_jacobi", the limits are variable

MPI is used to divide the outermost integration equally between processes.

Inputs:
* n: Simpson cells to split the three variables into (x_n, y_n, z_n)
* mixed_lims: The mixed Jacobi integration limits (R_a_min, R_a_max, R_b_min, R_b_max, gamma_ab_min, gamma_ab_max)
* mu_a
* mu_b
* m_a
* m_b
* mass_c
* mixed_jacobi: Flag to determine whether to integrate using mixed Jacobi coordinates if true
* single_jacobi: Flag to determine whether to integrate using single Jacobi coordinates if true
* test_coords: Flag to determine whether to integrate using Cartesian coordinates if true
* sub_comm: MPI communicator

Options:
* verbose: If true, progress of the integration and the number of limits not found is printed
* save_lims: If single_jacobi, saves the limits used for r_a and gamma_a in unformatted files (outputs/r_lims and outputs/gamma_lims)

Output:
* total_integral: Integral calculated

Errors:
* If more than one flag out of mixed_jacobi, single_jacobi and test_coords is set to true, the program will stop


### function integrate_MC

This aims to carry out the same integral as "integrate_triple_Simpson" for single_jacobi, using Monte Carlo integration, involving evaluating the integrand at randomly selected coordinates within the volume.

Inputs:
* n: Number of Monte Carlo samples
* mixed_lims: The mixed Jacobi integration limits (R_a_min, R_a_max, R_b_min, R_b_max, gamma_ab_min, gamma_ab_max)
* mu_a
* mu_b
* m_a
* m_b
* mass_c
* sub_comm: MPI communicator

Options:
* verbose: If true, progress of the integration and the number of limits not found is printed
* save_lims: If single_jacobi, saves the limits used for r_a and gamma_a in unformatted files (outputs/r_lims and outputs/gamma_lims)

Output:
* total_integral: Integral calculated


### function estimate_jacobi_lims

This estimates the minimum and maximum values of r_a, based on known R_a, or gamma_a, based on known R_a and r_a. This can be done through generation of random, valid coordinates, or uniform evaluation across the range of variables. The smallest and largest coordinates found will be returned. If no valid coordinates are found in this search, both limits are set to 0.

Inputs:
* R_a
* mixed_lims: The mixed Jacobi integration limits (R_a_min, R_a_max, R_b_min, R_b_max, gamma_ab_min, gamma_ab_max)
* mu_a
* m_b
* mass_c
* estimating_r_a: Estimate r_a limits if true
* estimating_gamma_a: Estimate gamma_a limits if true
* small_r_a: r_a (unused when estimating r_a limits)

Options:
* The number of test coordinates sampled, hardcoded via the size of test_coord in its declaration
* random_range: Flag to determine whether to sample randomly or uniformly across the variable range

Output:
* single_lims: An array of two values containing the estimated minimum and maximum r_a or gamma_a


### function get_limits

This gets the limits of r_a for a known R_a, or the limits of gamma_a for a known R_a and r_a. The limits can be estimated or calcualted analytically in some cases. The validity of the analytic expressions must be considered carefully before use. Both limits are set to 0._wp if valid coordinates cannot be found function.

Inputs:
* get_r_lims: Flag to determine whether the r_a limits are returned
* get_gamma_lims: Flag to determine whether the gamma_a limits are returned
* estimate_lims: Flag to determine whether estimate_jacobi_lims is used to estimate the limits, or whether an analytic expression is used
* R_a
* small_r_a: r_a (unused when estimating r_a limits)
* mixed_lims: The mixed Jacobi integration limits (R_a_min, R_a_max, R_b_min, R_b_max, gamma_ab_min, gamma_ab_max)
* mu_a
* m_b
* mass_c

Options:
* gamma_tol: A tolerance for the value of cos(gamma_a), which may be slightly outside the valid range (-1 - 1) due to numerical errors

Outputs:
* lims: An array of two values containing the estimated or analytic minimum and maximum r_a or gamma_a

Errors:
* If both get_r_lims and get_gamma_lims are set to true, the program will stop


### function get_test_limits

This calculates the limits needed for the test integrand. The y limits are currently defined as 0 - 1-x, and the z limits are currently defined as 0 - 1-x-y.

Inputs:
* get_y_lims: Flag to determine whether to get the limits for y
* get_z_lims: Flag to determine whether to get the limits for z
* x
* y

Outputs:
* lims: An array of two values containing the minimum and maximum y or z


Errors:
* If both get_y_lims and get_z_lims are set to true, the program will stop


### function test_integrand_func

This defines and evaluates an integrand for a three variables: x, y and z. Currently, this returns (x+y+z)^-0.5.

Inputs:
* x
* y
* z

Output:
* integrand: The evaluated integrand