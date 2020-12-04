#ifndef KSPACE_H
#define KSPACE_H

// Standard includes
#include <VotcaDefs.h>

#include "UnitCell.h"

class KSpace {
public:
  KSpace(const double alpha, const double k_max, const UnitCell &unit_cell)
      : alpha(alpha), k_max(k_max), unit_cell(unit_cell){};

  void compute(const std::vector<Eigen::Vector3d> &xyz,
               const std::vector<double> &charges);

  double getAk(const Eigen::Vector3d &k) const;

  const double getEnergy() const { return energy; }

  void updateAlpha(double a) { alpha = a; }
  void updateK_max(double k) { k_max = k; }

private:
  double alpha;
  double k_max;
  UnitCell unit_cell;
  const double pi = boost::math::constants::pi<double>();

  // RESULTS
  double energy;
};


#endif