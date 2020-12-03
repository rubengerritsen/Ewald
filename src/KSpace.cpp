
#include "KSpace.h"
#include <complex>
#include <iostream>

void KSpace::compute(const std::vector<Eigen::Vector3d> &xyz,
                     const std::vector<double> &charges) {
  double energy_sum = 0.0;
  double se = 0.0;
  const std::complex<double> imaginary(0.0, 1.0);

  Index k_max;



  for (Index ix = -k_max; ix <= k_max; ++ix) {
    for (Index iy = -k_max; iy <= k_max; ++iy) {
      for (Index iz = -k_max; iz <= k_max; ++iz) {
        if (ix == 0 && iy == 0 && iz == 0) {
          continue;
        }
        double charge_sum = 0.0;
        Eigen::Vector3d kvector = unit_cell.getKVector(ix, iy, iz);
        for (Index j = 0; j < xyz.size(); j++) {
          for (Index i = 0; i < xyz.size(); i++) {
            charge_sum += charges[i] * charges[j] * std::cos(kvector.dot(xyz[i] - xyz[j]));
          }
        }

        std::cout << "ChargeSum: " << charge_sum << std::endl;
        energy_sum += getAk(kvector) * charge_sum;
        std::cout << energy_sum << "   " << getAk(kvector) << std::endl;
      }
    }
  }

  for (const auto &q : charges) {
    se += q * q;
  }

  energy = (2 * pi / unit_cell.getVolume()) * energy_sum;
}