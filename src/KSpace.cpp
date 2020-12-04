
#include "KSpace.h"
#include <complex>
#include <iostream>

void KSpace::compute(const std::vector<Eigen::Vector3d> &xyz,
                     const std::vector<double> &charges) {
  double energy_sum = 0.0;
  double se = 0.0;
  const std::complex<double> imaginary(0.0, 1.0);

  for (Index ix = -k_max; ix <= k_max; ++ix) {
    for (Index iy = -k_max; iy <= k_max; ++iy) {
      for (Index iz = -k_max; iz <= k_max; ++iz) {
        if (ix == 0 && iy == 0 && iz == 0) {
          continue;
        }
        if (ix * ix + iy * iy + iz * iz > k_max * k_max) {
          continue; // sum is conditionally converged we should sum spherically
        }
        std::complex<double> charge_sum = 0.0;
        Eigen::Vector3d kvector = unit_cell.getKVector(ix, iy, iz);
        for (Index j = 0; j < xyz.size(); j++) {
          charge_sum += charges[j] * std::exp(imaginary * kvector.dot(xyz[j]));
        }
        energy_sum += getAk(kvector) * std::norm(charge_sum);
      }
    }
  }
  energy = (2 * pi / unit_cell.getVolume()) * energy_sum;
}

double KSpace::getAk(const Eigen::Vector3d &k) const {
  /* Compute the A_k factor, i.e. k^(-2) exp(-k^2/(4\alpha^2)) */
  double k_squared = k.squaredNorm();
  return std::exp(-k_squared / (4 * alpha * alpha)) / k_squared;
}