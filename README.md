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



## transform.f90

This program is designed to transform surface amplitudes at the boudary between the inner region, defined by mixed Jacobi coordinates (R_a, R_b, gamma_ab), and the outer regions, defined by single Jacobi coordinates ((R_a, r_a, gamma_a) or (R_b, r_b, gamma_b)), where the boundary is defined by constant R_a or R_b respectively.

This requires the calculation of the inner products of wavefunctions within the inner region, over the two variables not held fixed at the relevant boundary, which is performed using Simpson's rule to carry out the double integration. The integration is carried out using both mixed Jacobi coordinates and single Jacobi coordinates for comparison, although the mixed Jacobi coordinates are expected to be preferable, as the integration limits are fixed in this basis.

Note:

* Prior to transforming the amplitudes, the program currently also calculates the "old" amplitudes in the mixed Jacobi basis. This should not be required in the future, as these amplitudes may be calculated much more efficiently, but this functionality allows for testing of the more efficient method, and provides placeholder amplitudes to be transformed

* The wavefunctions used to calculate the "old" and transformed amplitudes are arbitrary placeholders that do not meet all requirements. However, the channel functions have been defined such that they are zero at the other channel's boundary

* Many of the functions within this program are adapted from jacobi.f90

Unless otherwise stated:
* config_a defines whether currently considering channel A (true) or channel B (false)
* mixed_int is the flag to determine whether to carry out integrations using mixed Jacobi (true) or single Jacobi (false) coordinates

* (R_1, small_r_1 and gamma_1) are the single Jacobi coordinates for channel A (R_a, r_a, gamma_a) or B (R_b, r_b, gamma_b)
* (R_1, R_2 and gamma_ab) are mixed Jacobi coordinates, with R_2 representing R_b or R_a for channel A or B respectively

* mass_1 is the mass of atom A for channel A, or atom B for channel B
* mass_2 is the mass of atom B for channel A, or atom A for channel B
* mass_c is the mass of atom C
* mass_total is the sum of the masses of A, B and C
* mu_1 is the reduced channel mass for channel A (mu_a)
* mu_2 is the reduced channel mass for channel B (mu_b)
* m_1 is the internal reduced mass for channel A (m_b)
* m_2 is the internal reduced mass for channel B (m_a)

* simpson_n: Number of Simpson cells used for mixed and single Jacobi integration using Simpson's rule
* mixed_lims: Integration limits used for mixed and single Jacobi integration, in terms of mixed Jacobi coordinates
* boundary_val_1 is the boundary value of R_a or R_b for channel A or B respectively
* boundary_val_2 is the boundary value of R_b or R_a for channel A or B respectively

* coords: Grid of coordinates defining the Simpson cells to be used during integration
* trans_coords: Grid of coordinates in coordinates not being used during integration

* n_1 is the number of channel functions for configuration A
* nc_1 is the number of radial continuum basis orbitals for configuration A
* n_2 is the number of channel functions for configuration B
* nc_2 is the number of radial continuum basis orbitals for configuration B
* b is the number of quadratically integrable functions (not used currently)
* nt is the number of linearly independent basis functions

* channel_func_1 is a flag determining whether the bra is a channel function (true) or basis function (false)
* single_func_1 is a flag determining whether the bra is in single Jacobi (true) or mixed Jacobi (false) coordinates
* channel_func_2 is a flag determining whether the ket is a channel function (true) or basis function (false)
* single_func_2 is a flag determining whether the ket is in single Jacobi (true) or mixed Jacobi (false) coordinates

* sub_comm is a (subset of) the MPI communicator
* sub_comm_size is the size of sub_comm
* rank is the MPI process rank


### module precisn

This defines the precision of real numbers throughout the program, currently such that real numbers are stored up to 12 decimal places.


### subroutine main

This sets up the main integration parameters used throughout the program and calls the function that calculates and transforms the surface amplitudes. When all code is uncommented, these calculations are performed in both the mixed and single Jacobi coordinates, and the total time for this in each set of coordinates is output.

Options:
* save_inputs: If true, saves the values of masses, Jacobi integration limits and coordinate system for the integration in an unformatted file (outputs/values)
* simpson_n
* mixed_lims
* mixed_int: When all code is uncommented, this simply defines the order these are performed, as the transformation is performed using both coordinates
* n_1
* nc_1
* n_2
* nc_2
* b


### subroutine calculate_masses

See jacobi.f90


### function get_single_phi

This evaluates the ith channel function in single Jacobi coordinates for configuration A or B. Currently arbitrarily defined, other than to be zero at the opposite channel boundary.

Inputs:
* i: Index of channel function
* boundary_val_2
* R_1
* small_r_1
* gamma_1
* mu_1
* m_1
* mass_c
* config_a
* conj: Flag to specify complex conjugate (currently unused)

Outputs:
* func: Evaluated function

### function get_mixed_phi

This evaluates the ith channel function in mixed Jacobi coordinates for configuration A or B. Currently arbitrarily defined, other than to be zero at the opposite channel boundary.

Inputs:
* i: Index of channel function
* boundary_val_2
* R_1
* R_2
* gamma_ab
* config_a
* conj: Flag to specify complex conjugate (currently unused)

Output:
* func: Evaluated function


### function get_single_psi

This evaluates the kth basis function in single Jacobi coordinates for configuration A or B, using the relevant channel functions, radial continuum orbital functions and coefficients defined in their own functions.

Inputs:
* k: Index of basis function
* n_1
* nc_1
* boundary_val_2
* R_1
* small_r_1
* gamma_1
* mu_1
* m_1
* mass_c
* config_a
* conj: Flag to specify complex conjugate (currently unused)

Output:
* func: Evaluated function


### function get_mixed_psi

This evaluates the kth basis function in mixed Jacobi coordinates for configuration A or B, using the relevant channel functions, radial continuum orbital functions and coefficients defined in their own functions.

Inputs:
* k: Index of basis function
* n_1
* nc_1
* boundary_val_2
* R_1
* R_2
* gamma_ab
* config_a
* conj: Flag to specify complex conjugate (currently unused)

Output:
* func: Evaluated function


### function get_single_radial_func

This evaluates the radial continuum basis orbital (i, j) in single Jacobi coordinates. Currently arbitrarily defined.

Inputs:
* R_1
* i: Index of radial function
* j: Index of radial function
* config_a

Output:
* radial_func: Evaluated function


### function get_mixed_radial_func

This evaluates the radial continuum basis orbital (i, j) in mixed Jacobi coordinates. Currently arbitrarily defined.

Inputs:
* R_1
* i: Index of radial function
* j: Index of radial function
* config_a

Output:
* radial_func: Evaluated function


### function get_single_coeff

This evaluates the coefficient (i, j, k) in single Jacobi coordinates. Currently arbitrarily defined.

Inputs:
* i: Index of coefficient
* j: Index of coefficient
* k: Index of coefficient
* config_a

Output:
* coeff: Coefficient


### function get_mixed_coeff

This evaluates the coefficient (i, j, k) in mixed Jacobi coordinates. Currently arbitrarily defined.

Inputs:
* i: Index of coefficient
* j: Index of coefficient
* k: Index of coefficient
* config_a

Output:
* coeff: Coefficient


### function integrand_func

This evaluates the integrand at (x, y, z), including the volume element R_2^2 sin(gamma_ab) or small_r_2^2 sin(gamma_1), using a number of flags and indices to determine the relevant bra and ket functions being that are integrated. See also: jacobi.f90: single_jacobi_integrand_func

Inputs:
* x: Jacobi coordinate, representing R_1
* y: Jacobi coordinate, representing R_2 (mixed_int) or small_r_1 (not mixed_int)
* z: Jacobi coordinate, representing gamma_ab (mixed_int) or gamma_1 (not mixed_int)
* boundary_val_2
* mu_1
* m_1
* mass_c
* trans_coords
* idx_1: Bra function index
* idx_2: Ket function index
* n_1
* nc_1
* config_a
* mixed_int
* channel_func_1
* single_func_1
* channel_func_2
* single_func_2

Output:
* total_integrand: Evaluated integrand


### function integrate_double_Simpson

This uses Simpson's rule to integrate two variables. See also jacobi.f90: integrate_triple_Simpson

Inputs:
* simpson_n: Number of Simpson cells used for mixed and single Jacobi integration using Simpson's rule
* x: R_1
* boundary_val_2
* n_1
* nc_1
* idx_1: Bra function index
* idx_2: Ket function index
* mu_1
* mu_2
* m_1
* m_2
* mass_c
* coords
* trans_coords
* config_a
* mixed_int
* channel_func_1
* single_func_1
* channel_func_2
* single_func_2

Output:
* total_integral: Integral calculated


### function estimate_jacobi_lims

Estimates the range of small_r_1 or gamma_1. See jacobi.f90 (note: parameters passed in different order).


### function get_jacobi_lims

Gets the limits of small_r_1 or gamma_1, through estimation or analytically. See jacobi.f90:get_limits (note: parameters passed in different order, generalised for channels A and B).


### function calc_small_r_from_mixed

Calculates small_r_1 from mixed Jacobi coordinates.

Inputs:
* R_1
* R_2
* gamma_ab
* mu_1
* m_1
* mass_c

Output:
* r

Errors:
* If small_r_1 is calculated to be less than 0, a warning will be printed, and small_r_1 will be set to 0


### function calc_R_from_single

Calculates R_2 from single Jacobi coordinates.

Inputs:
* R_1
* small_r_1
* gamma_1
* mu_1
* m_1
* mass_c
* config_a

Output:
* R_2


### function calc_gamma_from_mixed

Calculates gamma_1 from mixed Jacobi coordinates. See jacobi.f90:calc_gamma_a (note: generalised for channels A and B).


### function calc_gamma_from_single

Calculates gamma_ab from single Jacobi coordinates.

Inputs:
* R_1
* small_r_1
* gamma_1
* mu_1
* m_1
* mass_c
* config_a

Options:
* gamma_tol: A tolerance for the value of cos(gamma_a), which may be slightly outside the valid range (-1 - 1) due to numerical errors

Output
* gamma_ab

Errors:
* If cos(gamma_a) is outside its valid range, including the defined tolerance, the program will stop. This may be because the tolerance is too small, or there may be an error with the inputs


### function calc_gamma_from_integral

Calculates gamma_1 from R_1, r_1 and gamma_ab.  See jacobi.f90:calc_gamma_a (note: generalised for channels A and B).


### function transform_mixed_to_single

Calculates a single Jacobi coordinate from a mixed Jacobi coordinate.

Inputs:
* R_1
* R_2
* gamma_ab
* mu_1
* m_1
* mass_c
* config_a

Outputs:
single_coords: (small_r_1, gamma_1)


### function transform_single_to_mixed

Calculates a mixed Jacobi coordinate from single Jacobi coordinate.

Inputs:
* R_1
* small_r_1
* gamma_1
* mu_1
* m_1
* mass_c
* config_a

Output:
* mixed_coords: (R_2, gamma_ab)


### function transform_grid

Transforms a grid of mixed Jacobi coordinates to single Jacobi coordinates or vice versa, with the direction of the transformation determined the coordinates being used in the integration (via the mixed_int flag).

Inputs:
* x: R_1
* coords
* n: Number of coordinates for each of the two variables to transform
* mu_1
* m_1
* mass_c
* config_a
* mixed_int
* sub_comm
* sub_comm_size
* rank

Options:
* verbose: If true, progress is printed

Output:
* trans_coords


### function get_grid

Determines the grid of coordinates defining the Simpson cells to be used during integration, with the choice of single or mixed Jacobi coordinates determined the mixed_int flag.

Inputs:
* x: R_1
* boundary_val_2
* n: Number of coordinates for each of the two variables defining the grid
* mixed_lims
* mu_1
* m_1
* mass_c
* config_a
* mixed_int
* sub_comm
* sub_comm_size
* rank

Options:
* verbose: If true, progress is printed

Output:
* coords


### function transform_amps

Transform surface amplitudes w_ik from mixed ("old") to single Jacobi coordinates.

Inputs:
* old_amps: Surface amplitudes in mixed Jacobi coordinates
* n_1
* nc_1
* nt
* simpson_n: Number of Simpson cells used for mixed and single Jacobi integration using Simpson's rule
* boundary_val_1
* boundary_val_2
* mu_1
* mu_2
* m_1
* m_2
* mass_c
* coords
* trans_coords
* config_a
* mixed_int
* sub_comm
* sub_comm_size
* rank

Options:
* verbose: If true, progress is printed

Output:
* amps: Transformed amplitudes (in the single Jacobi basis)


### function calc_old_amps

Calculates the surface amplitudes in the mixed Jacobi basis. For appropriately chosen channel functions, this should not be required in the future.

Inputs:
* n: n_1 or n_2
* nc: nc_1 or nc_2
* nt
* simpson_n
* coords
* trans_coords
* boundary_val_1
* boundary_val_2
* mu_1
* mu_2
* m_1
* m_2 mass_c
* config_a
* mixed_int
* sub_comm
* sub_comm_size
* rank

Options:
* verbose: If true, progress is printed

Output
* amps: Surface amplitudes in mixed Jacobi basis ("old")


### subroutine transform_all_amps

This calculates surface amplitudes in mixed Jacobi coordinates ("old") and transforms these into single Jacobi coordinates, printing both the "old" and "new" amplitudes, and the time taken during various calculations. This includes calling functions to determine the grid to be used during integration in both mixed and single Jacobi coordinates, as well as functions to calculate the surface amplitudes in mixed Jacobi coordinates and transform these into single Jacobi coordinates. Calculations for channel A and B are performed seperately through the config_a flag.

Inputs:
* n_1
* nc_1
* n_2
* nc_2
* nt
* simpson_n
* mixed_lims
* mu_a
* mu_b
* m_a
* m_b
* mass_c
* mixed_int
* sub_comm

* Options:
* verbose: If true, progress is printed