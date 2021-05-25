#ifndef TURNING_NUMBER_H
#define TURNING_NUMBER_H

#include "field.h"

int turning_number_dual_loop(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& TT,
  const std::vector<int>& loop,
  const std::vector<Property>& prop
);

#endif