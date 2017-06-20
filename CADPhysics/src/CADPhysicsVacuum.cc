#include "CADPhysicsVacuum.hh"

// Create IMFP table.
// BETWEEN 1 AND 2 ELECTRONVOLTS (see definition of energy_axis),
// THIS ALWAYS GIVES ZERO IMFP. OUTSIDE THAT RANGE, IT GIVES NAN
// BECAUSE THE INTERPOLATION ROUTINE DOES NOT LIKE EXTRAPOLATION
// WHEN THE TWO VALUES ARE INFINITE.
typename material::imfp_table_t CADPhysicsVacuum::get_vacuum_imfp()
{
	using real = material::fast_real;
	ax_logspace<real> energy_axis(1, 2, 2);
	std::vector<real> data_vector(2, -std::numeric_limits<real>::infinity());

	return material::imfp_table_t::base_type
	(
		energy_axis, data_vector
	);
}

// Create ICDF table that always returns zero (scattering angle, energy loss, ...)
// Does work for all possible energies.
typename material::icdf_table_t CADPhysicsVacuum::get_vacuum_icdf()
{
	using real = material::fast_real;
	ax_logspace<real> energy_axis(1, 2, 2);
	ax_linspace<real> prob_axis(0, 1, 2);
	std::vector<real> data_vector(4, 0);

	return material::icdf_table_t::base_type
	(
		energy_axis, prob_axis, data_vector
	);
}

// Create ionization table that always returns -1 binding energy
// Also works for all possible energies.
typename material::ionization_table_t CADPhysicsVacuum::get_vacuum_ionization()
{
	using real = material::fast_real;
	ax_logspace<real> energy_axis(1, 2, 2);
	ax_linspace<real> prob_axis(0, 1, 2);
	std::vector<real> data_vector(4, -1);

	return material::ionization_table_t::base_type
	(
		energy_axis, prob_axis, data_vector
	);
}
