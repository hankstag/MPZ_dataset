#include "field/field.h"
#include "argh.h"
#include <igl/readOBJ.h>

int main(int argc, char *argv[])
{

  auto cmdl = argh::parser(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);

  std::string data_dir, out_dir;
  std::string model, name;
  cmdl("-m") >> model;
  cmdl("-d") >> data_dir;
  cmdl("-o") >> out_dir;
  
  if(out_dir == "") out_dir = data_dir; // use same dir as input data by default

  name = model.substr(0, model.find_last_of('.'));
  name = name.substr(name.find_last_of('/')+1);
  std::string field_name = name+".ffield";
  
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::readOBJ(data_dir + "/" + model, V, F);

  Eigen::VectorXd S;
  std::vector<Property> props;
  Eigen::VectorXd theta;

  std::string field_file = data_dir + "/" + name + ".ffield";
  if(V.cols() != 3) V.conservativeResize(V.rows(), 3); // some models have more than 3 cols
  load_ffield(V, F, field_file.c_str(), props);
  
  std::cout<<"load field done. props loaded: "<<props.size()<<"/"<<F.rows()<<std::endl;

  std::string out_file = out_dir + "/" + name + ".hdf5";
  std::cerr << "file saved to " << out_file << std::endl;
  save_ffield_hdf5(out_file, V, F, props);
  
  std::cout<<"save field done.\n";

}
