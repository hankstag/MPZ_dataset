#ifndef ANGLES_H
#define ANGLES_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include "tangent_frame.h"

typedef CVec3T<double> CVec3D;

double cos_angle_rescaled_tri(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& TT,
  int f, int k
);

double signed_angle(
  const Eigen::Vector3d& v1,
  const Eigen::Vector3d& v2,
  const Eigen::Vector3d& n
);

template <typename VectorType>
double signed_angle(const VectorType& v1, const VectorType& v2, const VectorType &normal);

void scaled_triangle(
  const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXi& TT, int f, int k,
  Eigen::Vector3d& p0, Eigen::Vector3d& p1, Eigen::Vector3d& p2
);

void scaled_triangle(
  const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXi& TT, int f, int k,
  CVec3D &p0, CVec3D &p1, CVec3D &p2
);

#endif