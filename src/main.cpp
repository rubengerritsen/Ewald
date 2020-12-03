#include "KSpace.h"
#include "RSpace.h"
#include <VotcaDefs.h>
#include <complex>
#include <iostream>

int main(int argc, char **argv) {

  Index cry_l = 16;
  Index N = 16 * 16 * 16;
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

  double alpha = 1.2;
  double r_max = 0.3;
  double k_max = 4;

  std::array<Eigen::Vector3d, 3> unitCell;
  unitCell[0] << 8.0, 0, 0;
  unitCell[1] << 0, 8.0, 0;
  unitCell[2] << 0, 0, 8.0;
  UnitCell cell(unitCell[0], unitCell[1], unitCell[2]);

  KSpace kspace(alpha, k_max, cell);

	kspace.compute(xyz, charges);

  RSpace rspace(alpha, r_max, cell);

  rspace.compute(xyz, charges);




	std::cout << "KSpace energy: " << kspace.getEnergy() << std::endl;
  std::cout << "RSpace energy: " << rspace.getEnergy() << std::endl;

}
