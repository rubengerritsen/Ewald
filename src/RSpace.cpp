
#include "RSpace.h"
#include <complex>
#include <iostream>
#include <cmath>

void RSpace::compute(const std::vector<Eigen::Vector3d> &xyz,
                     const std::vector<double> &charges) {

  energy = 0.0;

  for (Index i = 0; i < xyz.size(); ++i) {
    for (Index j = i + 1; j < xyz.size(); ++j) {
      Eigen::Vector3d dr(xyz[i] - xyz[j]);
      double dist = dr.norm();
      if(dist < r_max){
        energy += charges[i] * std::erfc(alpha * dist)/dist * charges[j];
      }
    }
  }

  // Self energy correction
  double U_self = 0.0;
  for(Index i = 0; i < xyz.size(); ++i){
    energy += -charges[i]*charges[i] * alpha / sqrt_pi;
  }
}