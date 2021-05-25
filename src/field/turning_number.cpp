#include "turning_number.h"
#include "angles.h"
#include <igl/PI.h>
#include <vector>

// discrete geodesic curvature of a cyclic triangle strip
double dual_cycle_geodesic_curvature(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& TT,
  const std::vector<int>& strip,
  const std::vector<Property>& prop
){
  double curvature_g = 0;
  int N = strip.size();
  assert(N >= 3 && "dual strip has to be at least length 3");
  for(int i = 0;i < N;i++){
    int f_prev = strip[(i-1+N)%N];
    int f_curr = strip[i];
    int f_next = strip[(i+1)%N];
    
    // edge_i, i \in {0,1,2}
    // means edge of endpoints {F(f,i), F(f,(i+1)%3)}
    
    int e_in = -1;  // from which edge get into f_curr
    int e_out = -1; // from which edge get out of f_curr
    for(int e = 0; e < 3; e++){
      if(TT(f_curr,e) == f_next) e_out = e;
      if(TT(f_curr,e) == f_prev) e_in = e;
    }
    bool left_turn = (e_in+1)%3 == e_out ? false : true;
    assert(e_in != -1 && e_out != -1 && e_in != e_out);
    
    Eigen::Vector3d p0, p1, p2;
    scaled_triangle(V,F,TT,f_curr,e_in,p0,p1,p2); // p0 is F(f,e_in)
    
    Eigen::Vector3d normal(0,0,1);
    if(left_turn){
      curvature_g += signed_angle(p1-p0, p2-p0, normal);
    }else{
      curvature_g += signed_angle(p0-p1, p2-p1, normal);
    }
  }
  return curvature_g;
}

int turning_number_dual_loop(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& TT,
  const std::vector<int>& loop,
  const std::vector<Property>& prop
){
  auto angle_between_edges = [&](int f, int k0, int k1){
    Eigen::RowVector3d vec1, vec2;
    int multiplier = -1;
    if((k0+1)%3 != k1){
      std::swap(k0,k1);
      multiplier = 1;
    }
    int k2 = (k1+1)%3;
    vec1 = V.row(F(f,k0)) - V.row(F(f,k1));
    vec2 = V.row(F(f,k2)) - V.row(F(f,k1));
    Eigen::Vector3d n = vec1.cross(vec2).normalized();
    return multiplier * signed_angle(vec1, vec2, n);
  };
  
  double curvature_g = 0;
  int pj_sum = 0;
  double kappa_sum = 0;
  int n = loop.size();
  for(int i=0;i<n;i++){
    int fk = loop[(i-1+n)%n];
    int fi = loop[ i       ];
    int fj = loop[(i + 1)%n];
    int e0 = 0, e1 = 0;
    for(int k=0;k<3;k++){
      if(TT(fi,k) == fk) e0 = k;
      if(TT(fi,k) == fj){
        pj_sum -= prop[fi].pj(k);
        kappa_sum -= prop[fi].kappa(k);
        e1 = k;
      }
    }
    assert(e0 != e1);
    // std::cout<<"signed angle = "<<angle_between_edges(fi,e0,e1)<<std::endl;
    curvature_g += angle_between_edges(fi,e0,e1);
    
  }
  double c_g = dual_cycle_geodesic_curvature(V,F,TT,loop,prop);
  // spdlog::info("c_g diff = {0:f} - {1:f} = {2:f}", curvature_g, c_g, curvature_g-c_g);
  // std::cout<<std::setprecision(17)<<curvature_g-c_g<<std::endl;
  // std::cout<<kappa_sum<<std::endl;
  double tn_val = pj_sum + 2*(kappa_sum - c_g)/igl::PI;
  int tn = std::round(pj_sum + 2*(kappa_sum - c_g)/igl::PI);
  // std::cout<<"float(tn) = "<<tn_val<<std::endl;
  assert(std::abs(tn - tn_val) < 1e-9);
  return tn;
}