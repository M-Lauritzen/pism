/* Copyright (C) 2026 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef PISM_SIA_SSA_WEIGHT_H
#define PISM_SIA_SSA_WEIGHT_H

#include <cmath>

namespace pism {
namespace stressbalance {

/*!
 * Bueler-Brown velocity weight for the optional weighted SIA+SSA hybrid.
 *
 *   f(|u_SSA|) = 1 - (2/pi) * atan( (|u_SSA| / u_ref)^2 )
 *
 * f -> 1 in the slow interior (keep the full SIA shear) and f -> 0 at fast
 * sliding cells (remove the SIA), with the f = 0.5 crossover at |u_SSA| = u_ref.
 * The SIA contribution is scaled by f and the SSA contribution by (1 - f), so
 * the two are BLENDED rather than added -- this is the pre-2011 PISM `-super`
 * behavior, and it avoids counting the observed velocity twice at outlets where
 * tau_c was inverted under SSA alone.
 *
 * Defined here, in one place, because the same weight has to be applied in three
 * places that would otherwise drift apart: the mass-continuity flux
 * (GeometryEvolution::compute_interface_fluxes), the reported 3D velocity
 * (SIAFD::compute_3d_horizontal_velocity) and the `sia_ssa_weight` diagnostic.
 *
 * @param[in] ssa_speed          |u_SSA| in the same units as reference_velocity
 * @param[in] reference_velocity u_ref > 0, the f = 0.5 crossover speed
 */
inline double sia_ssa_velocity_weight(double ssa_speed, double reference_velocity) {
  double r = ssa_speed / reference_velocity;
  return 1.0 - (2.0 / M_PI) * std::atan(r * r);
}

} // end of namespace stressbalance
} // end of namespace pism

#endif /* PISM_SIA_SSA_WEIGHT_H */
