#ifndef CADPhysicsVacuum_h
#define CADPhysicsVacuum_h

#include <limits>
#include <csread/material.h>

class CADPhysicsVacuum
{
public:
// Create IMFP table.
// BETWEEN 1 AND 2 ELECTRONVOLTS (see definition of energy_axis),
// THIS ALWAYS GIVES ZERO IMFP. OUTSIDE THAT RANGE, IT GIVES NAN
// BECAUSE THE INTERPOLATION ROUTINE DOES NOT LIKE EXTRAPOLATION
// WHEN THE TWO VALUES ARE INFINITE.
typename material::imfp_table_t get_vacuum_imfp();

// Create ICDF table that always returns zero (scattering angle, energy loss, ...)
// Does work for all possible energies.
typename material::icdf_table_t get_vacuum_icdf();
// Create ionization table that always returns -1 binding energy
// Also works for all possible energies.
typename material::ionization_table_t get_vacuum_ionization();

typename material::imfp_table_t get_vacuum_range();
};

#endif
