
// Standard includes
#include <VotcaDefs.h>
#include <iostream>

#include "UnitCell.h"

class RSpace {
public:
  RSpace(const double alpha, const double r_max, const UnitCell &unit_cell)
      : alpha(alpha), r_max(r_max), unit_cell(unit_cell) {
    if (r_max >= unit_cell.maxRCutOff()) {
      throw std::runtime_error("Real space cutoff is greater than half the "
                               "perpendicular width of the simulation box.");
    }
  }

  void compute(const std::vector<Eigen::Vector3d> &xyz,
               const std::vector<double> &charges);

  double getEnergy() const { return energy; }

  double getSelfEnergy() const { return selfEnergy; }

  void updateAlpha(double a) { alpha = a; }
  void updateR_max(double r) { r_max = r; }

private:
  double alpha;
  double r_max;
  UnitCell unit_cell;
  const double pi = boost::math::constants::pi<double>();
  const double sqrt_pi = std::sqrt(pi);

  double energy;
  double selfEnergy;
};