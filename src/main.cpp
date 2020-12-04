#include "KSpace.h"
#include "RSpace.h"
#include <VotcaDefs.h>
#include <complex>
#include <fstream>
#include <iostream>

int main(int argc, char **argv) {

  std::cout << "************* TEST 1: NaCl *************" << std::endl;
  Index cry_l = 16;
  Index N = cry_l * cry_l * cry_l;
  Index l = (double)cry_l / 2.0;
  std::vector<Eigen::Vector3d> xyz(N);
  std::vector<double> charges(N);
  for (Index i = 0; i < xyz.size(); ++i) {
    Index ix = i % cry_l;
    Index iy = (i % (cry_l * cry_l)) / cry_l;
    Index iz = i / (cry_l * cry_l);

    xyz[i][0] = (double)(ix)*0.5;
    xyz[i][1] = (double)(iy)*0.5;
    xyz[i][2] = (double)(iz)*0.5;

    charges[i] = (((ix + iy + iz) % 2) ? 1.0 : -1.0);
  }
  double alpha = 4;
  double r_max = 1.9;
  double k_max = 30;
  std::array<Eigen::Vector3d, 3> unitCell;
  unitCell[0] << l, 0, 0;
  unitCell[1] << 0, l, 0;
  unitCell[2] << 0, 0, l;
  UnitCell cell(unitCell[0], unitCell[1], unitCell[2]);
  RSpace rspace(alpha, r_max, cell);
  rspace.compute(xyz, charges);

  KSpace kspace(alpha, k_max, cell);
  kspace.compute(xyz, charges);

  std::cout << "RSpace Energy: "
            << (rspace.getEnergy() + rspace.getSelfEnergy()) << std::endl;
  std::cout << "KSpace energy: " << kspace.getEnergy() << std::endl;
  std::cout << "Total Energy: "
            << (rspace.getEnergy() + rspace.getSelfEnergy()) +
                   kspace.getEnergy()
            << std::endl;
  std::cout << "Madelung: "
            << -(rspace.getEnergy() + rspace.getSelfEnergy() +
                 kspace.getEnergy()) /
                   N
            << std::endl;
  std::cout << "Madelung Expected: " << 1.747565 << std::endl;
}
