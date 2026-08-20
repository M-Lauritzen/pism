/* Copyright (C) 2018 PISM Authors
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

#ifndef LAPSE_RATES_H
#define LAPSE_RATES_H

namespace pism {

namespace array { class Scalar; }

void lapse_rate_correction(const array::Scalar &surface,
                           const array::Scalar &reference_surface,
                           double lapse_rate,
                           array::Scalar &result);

/*!
 * Same, but with a SPATIALLY VARYING (and possibly time-dependent) lapse rate.
 *
 * Needed for the ISMIP7 SMB-elevation feedback, which prescribes a per-cell,
 * time-varying vertical gradient (`dmrrodz`, runoff change with surface elevation)
 * rather than a single number.
 */
void lapse_rate_correction(const array::Scalar &surface,
                           const array::Scalar &reference_surface,
                           const array::Scalar &lapse_rate,
                           array::Scalar &result);

} // end of namespace pism

#endif /* LAPSE_RATES_H */
