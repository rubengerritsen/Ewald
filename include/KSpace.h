

// Standard includes
#include <VotcaDefs.h>

#include "UnitCell.h"

class KSpace {
public:
  KSpace(const double alpha, const double k_max, const UnitCell &unit_cell)
      : alpha(alpha), k_max(k_max), unit_cell(unit_cell){};

  void compute(const std::vector<Eigen::Vector3d> &xyz,
               const std::vector<double> &charges);

  double getAk(const Eigen::Vector3d &k) {
    /* Compute the A_k factor, i.e. k^(-2) exp(-k^2/(4\alpha^2)) */
    double k_squared = k.squaredNorm();
    return std::exp(-k_squared / (4 * alpha * alpha)) / k_squared;
  }

  const double getEnergy() const { return energy; }

private:
  double alpha;
  double k_max;
  UnitCell unit_cell;
  const double pi = boost::math::constants::pi<double>();

  // RESULTS
  double energy;
};
