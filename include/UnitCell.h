#ifndef UNITCELL_H
#define UNITCELL_H

// Standard includes
#include <VotcaDefs.h>
#include <iostream>

class UnitCell {
public:
  UnitCell(Eigen::Vector3d L1, Eigen::Vector3d L2, Eigen::Vector3d L3) {
    /**
     * Based on GROMACS's triclinic boxes, that need to satisfy
     * a_y = a_z = b_z = 0
     * a_x > 0, b_y > 0, c_z > 0
     * b_x < 0.5 a_x, c_x < 0.5 a_x, c_y < 0.5 b_y
     */
    bool is_gromacs_triclinic_box = L1[1] == 0 && L1[2] == 0 && L2[2] == 0 &&
                                    L1[0] > 0 && L2[1] > 0 && L3[2] > 0 &&
                                    L2[0] < 0.5 * L1[0] &&
                                    L3[0] < 0.5 * L1[0] && L3[1] < 0.5 * L2[1];
    assert(is_gromacs_triclinic_box);

    cell_matrix << L1, L2, L3;
    cell_matrix_inv = cell_matrix.inverse();
    cell_volume = L1.dot(L2.cross(L3));
  }

  Eigen::Vector3d getKVector(Index nx, Index ny, Index nz) const {
    Eigen::Vector3d n(nx, ny, nz);
    return getKVector(n);
  }

  Eigen::Vector3d getKVector(Eigen::Vector3d n) const {
    return 2 * boost::math::constants::pi<double>() * cell_matrix_inv * n;
  }

  Eigen::Vector3d getLVector(Index nx, Index ny, Index nz) const {
    Eigen::Vector3d n(nx, ny, nz);
    return getLVector(n);
  }

  Eigen::Vector3d getLVector(Eigen::Vector3d n) const {
    return cell_matrix * n;
  }

  Eigen::Vector3d minImage(Eigen::Vector3d v1, Eigen::Vector3d v2) const {
    Eigen::Vector3d r_tp = v1 - v2;
    Eigen::Vector3d r_dp =
        r_tp - cell_matrix.col(2) * std::round(r_tp.z() / cell_matrix(2, 2));
    Eigen::Vector3d r_sp =
        r_dp - cell_matrix.col(1) * std::round(r_dp.y() / cell_matrix(1, 1));
    return r_sp - cell_matrix.col(0) * std::round(r_sp.x() / cell_matrix(0, 0));
  }

  double maxRCutOff() const {
    // Calculates maximum cutoff for which the minimal image convention still
    // makes sense.
    Eigen::Vector3d AxB = cell_matrix.col(0).cross(cell_matrix.col(1));
    Eigen::Vector3d BxC = cell_matrix.col(1).cross(cell_matrix.col(2));
    Eigen::Vector3d CxA = cell_matrix.col(2).cross(cell_matrix.col(0));

    double Wa = std::abs(cell_matrix.col(0).dot(BxC)) / BxC.norm();
    double Wb = std::abs(cell_matrix.col(1).dot(CxA)) / CxA.norm();
    double Wc = std::abs(cell_matrix.col(2).dot(AxB)) / AxB.norm();

    return 0.5 * std::min(std::min(Wa, Wb), Wc);
  }

  const double getVolume() const { return cell_volume; }

private:
  double cell_volume;
  Eigen::Matrix3d cell_matrix;
  Eigen::Matrix3d cell_matrix_inv;
};

#endif