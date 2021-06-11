#include "edge_map.h"

void init_EMap(const Eigen::MatrixXi& F, EMap& em){
  em.clear();
  for(int i=0;i<F.rows();i++){
    for(int k=0;k<3;k++){
      int k_1 = (k+1)%3;
      int u = F(i,k);
      int v = F(i,k_1);
      em[std::make_pair(u,v)] = std::make_pair(i,k);
    }
  }
}
