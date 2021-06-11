#ifndef EDGE_MAP_H
#define EDGE_MAP_H

#include <unordered_map>
#include <Eigen/Core>

struct HASH{
  size_t operator()(const std::pair<int,int>&x)const{
    return std::hash<long long>()(((long long)x.first)^(((long long)x.second)<<32));
  }
};
using EMap = std::unordered_map<std::pair<int,int>, std::pair<int,int>, HASH>;

void init_EMap(const Eigen::MatrixXi& F, EMap& em);

#endif
