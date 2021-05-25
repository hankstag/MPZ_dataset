#include "angles.h"
#include <cfloat>

template <typename VectorType>
VectorType cross_prod(const VectorType& v1, const VectorType& v2)
{
  return VectorType(
    v1[1] * v2[2] - v1[2] * v2[1],
    v1[2] * v2[0] - v1[0] * v2[2],
    v1[0] * v2[1] - v1[1] * v2[0]);
}

template <typename VectorType>
double dot_prod(const VectorType& v1, const VectorType& v2)
{
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

double signed_angle(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& n)
{
  double s = n.dot(v1.cross(v2)); // compute the sin
  double c = v1.dot(v2);          // compute the cos
  
  if(s == 0 && c == 0){
    return 0.0;
  }else
    return atan2(s,c);
}

template <typename VectorType>
double signed_angle(const VectorType& v1, const VectorType& v2, const VectorType &normal)
{
    double s = dot_prod(normal, cross_prod(v1, v2));
    double c = dot_prod(v1, v2);
    const double angle = (s == 0 && c == 0) ? 0.0 : atan2(s, c);
    return angle;
}

void scaled_triangle(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& TT,
  int f, int k,
  Eigen::Vector3d& p0, 
  Eigen::Vector3d& p1, 
  Eigen::Vector3d& p2
){
  int u = F(f,k);
  int v = F(f,(k+1)%3);
  int t = F(f,(k+2)%3);
  p0 = Eigen::Vector3d(0, 0, 0);
  p1 = Eigen::Vector3d((V.row(u)-V.row(v)).norm(), 0, 0);
  // double cos_t = cos_angle_rescaled_tri(h->prev());
  double cos_t = cos_angle_rescaled_tri(V,F,TT,f,(k+2)%3);
  double sin_t = sqrt(1.0 - cos_t*cos_t);
  p2 = (V.row(t)-V.row(u)).norm() * Eigen::Vector3d(cos_t, sin_t, 0);
}

void scaled_triangle(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& TT,
  int f, int k,
  CVec3D &p0, CVec3D &p1, CVec3D &p2
){
  int u = F(f,k);
  int v = F(f,(k+1)%3);
  int t = F(f,(k+2)%3);
  p0 = CVec3D(0, 0, 0);
  p1 = CVec3D((V.row(u)-V.row(v)).norm(), 0, 0);
  // double cos_t = cos_angle_rescaled_tri(h->prev());
  double cos_t = cos_angle_rescaled_tri(V,F,TT,f,(k+2)%3);
  double sin_t = sqrt(1.0 - cos_t*cos_t);
  p2 = (V.row(t)-V.row(u)).norm() * CVec3D(cos_t, sin_t, 0);
}

// computes the cosine of the angle at h->vertex()
double cos_angle_rescaled_tri(
  // EP_base::Halfedge_const_handle h
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& TT,
  int f, int k
){
  // if (TT(f,k) == -1) return 1.0; // cos(0)
	// if ( h->facet_degree() != 3) return 1.0;

    //
    //     /\.
    //  l0/  \l2
    //   /____\.
    // 01  l1
    // Law of cosines:
    // l2^2 = l0^2 + l1^2 - 2*l0*l1*cos(theta01)
    // => cos(theta01) = (l0^2 + l1^2 - l2^2) / (2*l0*l1);
    int u = F(f,k);
    int v = F(f,(k+1)%3);
    int t = F(f,(k+2)%3);
    double l0 = (V.row(u)-V.row(v)).norm();// h->scaledLength();
    double l1 = (V.row(v)-V.row(t)).norm();
    double l2 = (V.row(t)-V.row(u)).norm();
    double cos01;
    if (l0 == 0 && l1 == 0 && l2 == 0) {
        cos01 = 0.5;
    }
    else {
        double denom01 = 2*l0*l1;
        cos01 = denom01 == 0 ? 1.0 : ( l0*l0 + l1*l1 - l2*l2 ) / denom01;
    }

    // Will (may?) get triangle inequality violations if the following
    // lines are uncommented.
    const double EPS = FLT_EPSILON;
    if (1.0 <= cos01 && cos01 < 1 + EPS) { cos01 = 1.0; }
    else if (-1.0 >= cos01 && cos01 > -1 - EPS) { cos01 = -1.0; }

    if (-1 > cos01 || cos01 > 1) {
      std::cerr << "Cosine " << cos01 << " is outside [-1, 1]\n";
    }
    return cos01;
}

template double signed_angle<CVec3T<double> >(CVec3T<double> const&, CVec3T<double> const&, CVec3T<double> const&);