/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace algebra::smatrix::index {

enum track_indices : unsigned int {
  // bound indices
  bound_loc0 = 0,
  bound_loc1 = 1,
  bound_phi = 2,
  bound_theta = 3,
  bound_qoverp = 4,
  bound_time = 5,
  bound_size = 6,

  // free indices
  free_pos0 = 0,
  free_pos1 = 1,
  free_pos2 = 2,
  free_time = 3,
  free_dir0 = 4,
  free_dir1 = 5,
  free_dir2 = 6,
  free_qoverp = 7,
  free_size = 8,
};

}  // namespace algebra::smatrix::index