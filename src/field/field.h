#ifndef FIELD_H
#define FIELD_H

#include <Eigen/Core>
#include <vector>

// collection of properties defined per face
class Property{

public:
  
  Property(){};
  Property(Eigen::Vector3i _pj, Eigen::Vector3d _kappa, \
            Eigen::MatrixXd _TP, double r): pj(_pj), kappa(_kappa), TP(_TP), angle(r) {};
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

double vector_to_theta(const Eigen::MatrixXd& frame, const Eigen::Vector3d& v);

template <typename Scalar>
void convert_vector_to_matrix(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::Matrix<Scalar,-1,1>& pj_v, Eigen::Matrix<Scalar,-1,-1>& pj);

template <typename Scalar>
void convert_matrix_to_vector(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::Matrix<Scalar,-1,-1>& pj, Eigen::Matrix<Scalar,-1,1>& pj_v);

void reference_frames(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXi& first_edges, std::vector<Eigen::MatrixXd>& frames);

void compute_angle_defect(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::VectorXd& A);

void compute_kappa(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXi& first_edges, const std::vector<Eigen::MatrixXd>& TPs, Eigen::MatrixXd& kappa);

bool load_ffield(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const char *fname, std::vector<Property>& props, bool load_into_frames = false);
void save_ffield_hdf5(std::string fname, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, std::vector<Property>& props);
void load_ffield_hdf5(std::string fname, Eigen::MatrixXd& V, Eigen::MatrixXi& F, std::vector<Property>& props);

void find_cones(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const std::vector<Property>& props, Eigen::VectorXd& S);

void theta_to_vector(const Eigen::VectorXd& theta, const std::vector<Eigen::MatrixXd>& frames, Eigen::MatrixXd& R);
void induce_theta_hat(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, std::vector<Property>& props, Eigen::VectorXd& theta_hat);
void write_theta_data(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, std::vector<Property>& props, Eigen::VectorXd& theta_hat, std::string out_dir, std::string model_name);

#endif