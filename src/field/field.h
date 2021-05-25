#ifndef FIELD_H
#define FIELD_H

#include <Eigen/Core>

// collection of properties defined per face
class Property{

public:
  
  Property();
  Property(Eigen::Vector3i _pj, Eigen::Vector3d _kappa, Eigen::MatrixXd _TP, double r): pj(_pj), kappa(_kappa), TP(_TP), angle(r) {};
  ~Property(){};

  // period jumps
  Eigen::Vector3i pj;
  
  // diff between neighbor frames
  Eigen::Vector3d kappa;
  
  // reference frame (2 principle directions)
  Eigen::MatrixXd TP;

  // angle of the field
  double angle;
  
};

bool load_ffield(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const char *fname,
  std::vector<Property>& props,
  bool load_into_frames = false
);

void save_ffield_hdf5(
  std::string fname,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  std::vector<Property>& props
);

void load_ffield_hdf5(
  std::string fname,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  std::vector<Property>& props
);


#endif