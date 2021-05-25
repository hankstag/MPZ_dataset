#include "field.h"
#include "util.h"
#include "turning_number.h"
#include "angles.h"
#include <highfive/H5Easy.hpp>
#include <vector>
#include <set>
#include <fstream>
#include <cfloat>
#include <igl/is_border_vertex.h>
#include <igl/HalfEdgeIterator.h>
#include <igl/edge_topology.h>
#include <igl/per_face_normals.h>
#include <igl/triangle_triangle_adjacency.h>

typedef CVec3T<double> CVec3D;

double pos_fmod(double x, double y)
{
  return (0 == y) ? x : x - y * floor(x/y);
}

// convert from an angle in local frame
// to a vector in world space
void theta_to_vector(const Eigen::VectorXd& theta, const std::vector<Eigen::MatrixXd>& frames, Eigen::MatrixXd& R){
  R.resize(theta.rows(),3);
  for(int i=0;i<theta.rows();i++){
    double a  = theta(i);
    Eigen::Vector2d vp(cos(a),sin(a));
    R.row(i) << vp.transpose() * frames[i];
  }
}

void find_cones(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const std::vector<Property>& props,
  Eigen::VectorXd& S
){
  S.setZero(V.rows());
  Eigen::MatrixXi TT, TTi;
  igl::triangle_triangle_adjacency(F,TT,TTi);
  std::vector<std::vector<int>> VF, VFi;
  igl::vertex_triangle_adjacency(V,F,VF,VFi);

  std::vector<bool> is_border = igl::is_border_vertex(F);
  for(int i=0;i<V.rows();i++){
    assert(!VF[i].empty());
    int f0 = VF[i][0];
    int e0 = VFi[i][0];
    int f1 = TT(f0,e0);
    if(is_border[i]) continue;
    auto ei = igl::HalfEdgeIterator<Eigen::MatrixXi,Eigen::MatrixXi,Eigen::MatrixXi>(F,TT,TTi,f0,e0);
    std::vector<int> one_ring;
    ei.flipF();
    do{
      one_ring.push_back(ei.Fi());
      ei.NextFE();
    }while(ei.Fi() != f1);
    assert(one_ring.size() == VF[i].size());
    S(i) = 1 + turning_number_dual_loop(V,F,TT,one_ring,props) / 4.0;
  }
}

// reference angles when computed on the scaled triangle (that actually
// lives in R^6)
void reference_angles_rescaled_tri(
  // Halfedge_const_handle he, double& d1, double& d2
  const Eigen::MatrixXd& V, 
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& TT,
  const Eigen::MatrixXi& TTi,
  const Eigen::VectorXi& first_edges,
  const std::vector<Eigen::MatrixXd>& TPs,
  int f, int k, double& d1, double& d2
){ 
    // assume two incident faces
    // MeshType::Facet_const_handle f1 = he->facet(); 
    // MeshType::Facet_const_handle f2 = he->opposite()->facet();
    int f0 = f;
    int f1 = TT(f,k);
    int k0 = k;
    int k1 = TTi(f,k);

    // compute a representation of each scaled triangle on the plane
    CVec3D p0, p1, p2, q0, q1, q2;
    // scaled_triangle(he, p0, p1, p2);
    scaled_triangle(V,F,TT,f0,k0,p0,p1,p2);
    scaled_triangle(V,F,TT,f1,k1,q0,q1,q2);

    // MeshType::Halfedge_const_handle ref1 = f1->halfedge();
    // MeshType::Halfedge_const_handle ref2 = f2->halfedge();
    // CVec3D refv1(TPs[f0](0,0), TPs[f0](0,1), TPs[f0](0,2));
    CVec3D refv1 = first_edges(f0) == k0 ? p0 - p2 :
                   first_edges(f0) == (k0+1)%3 ? p1 - p0 : p2 - p1;
        // ref1 == he         ? p1 - p0 :
        // ref1 == he->next() ? p2 - p1 : p0 - p2;
    // CVec3D refv2(TPs[f1](0,0), TPs[f1](0,1), TPs[f1](0,2));
    CVec3D refv2 = first_edges(f1) == k1 ? q0 - q2 :
                   first_edges(f1) == (k1+1)%3 ? q1 - q0 : q2 - q1;
        // ref2 == he->opposite()         ? q1 - q0 :
        // ref2 == he->opposite()->next() ? q2 - q1 : q0 - q2;
    // normals
    const CVec3D norm(0,0,1);
    // angles from the common edge to ref edges
    d1 = -signed_angle(refv1, p1 - p0, norm);
    d2 = -signed_angle(refv2, q1 - q0, norm);
    // spdlog::info("between faces {0:d} and {1:d}, d1 = {2:f}, d2 = {3:f}", f0, f1, d1, d2);
}

double diff_reference_angles(
  const Eigen::MatrixXd& V, 
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& TT,
  const Eigen::MatrixXi& TTi,
  const Eigen::VectorXi& first_edges,
  const Eigen::MatrixXi& has_priority,
  const std::vector<Eigen::MatrixXd>& TPs,
  int f, int k
){
  
  if (!has_priority(f,k))
    return -diff_reference_angles(V, F, TT, TTi, first_edges, has_priority, TPs, TT(f,k), TTi(f,k));

  // angle between reference directions in faces i and j
  double di, dj;
  reference_angles_rescaled_tri(V,F,TT,TTi,first_edges,TPs,f,k,di,dj);
  double kij = di-dj+M_PI;
  // to be consistent with mixed integer, ensure that kij \in (-pi, pi]
  kij = pos_fmod(kij, 2*M_PI);
  if (kij > M_PI) kij -= 2*M_PI;
	if ( !(-M_PI < kij && kij <= M_PI) ){
    // Log.error() << "assert(-M_PI < kij && kij <= M_PI) failed" << endl;
    std::cerr << "assert(-M_PI < kij && kij <= M_PI) failed!\n";
  }
  return kij;
}

void compute_kappa(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::VectorXi& first_edges,
  const std::vector<Eigen::MatrixXd>& TPs,
  Eigen::MatrixXd& kappa
){
  Eigen::MatrixXi TT, TTi;
  igl::triangle_triangle_adjacency(F,TT,TTi);
  kappa.setZero(F.rows(),3);
  Eigen::MatrixXi has_priority(F.rows(),3);
  has_priority.setZero();
  for(int f=0;f<F.rows();f++){
    for(int k=0;k<3;k++){
      int u = F(f,k);
      int v = F(f,(k+1)%3);
      if(u > v) has_priority(f,k) = 1;
    }
  }
  for(int f=0;f<F.rows();f++){
    for(int k=0;k<3;k++){
      if(TT(f,k) != -1)
        kappa(f,k) = diff_reference_angles(V,F,TT,TTi,first_edges,has_priority,TPs,f,k);
    }
  }
}

void reference_frames(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::VectorXi& first_edges,
  std::vector<Eigen::MatrixXd>& frames
){
  Eigen::MatrixXd N;
  
  // Generate normals per face
  igl::per_face_normals(V, F, N);
  frames.clear();
  
  // Generate reference frames
  for(unsigned f=0; f<F.rows(); ++f){
    // First edge
    int index = first_edges(f);
    int u = F(f, (index+2)%3);
    int v = F(f, index);
    // spdlog::info("first edge of face {0:d} is #{1:d} - ({2:d}, {3:d})", f, index, u, v);
    Eigen::Vector3d e1 = V.row(v) - V.row(u);
    e1.normalize();
    Eigen::Vector3d e2 = N.row(f);
    e2 = e2.cross(e1);
    e2.normalize();

    Eigen::MatrixXd TP(2,3);
    TP << e1.transpose(), e2.transpose();
    frames.push_back(TP);
  }
}

bool load_ffield(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const char *fname,
  std::vector<Property>& props,
  bool load_into_frames
){

  Eigen::VectorXd S;
  Eigen::VectorXd angles;
  
  Eigen::MatrixXi TT, TTi;
  igl::triangle_triangle_adjacency(F, TT, TTi);
  
  Eigen::MatrixXi sharp(F.rows(),3);
  sharp.setZero();
  angles.setZero(F.rows());
  Eigen::VectorXi fixed_shapeop_f(F.rows());
  Eigen::VectorXi fixed_shapeop_v(V.rows());
  fixed_shapeop_f.setZero();
  fixed_shapeop_v.setZero();
  
  std::ifstream fin( fname );
  if ( !fin ) {
    std::cerr << fname << ": Unable to open .ffield file\n";
    return false;
  }

  int version = 0;
  
  // read the line that indicates whether the data is curvature or target frame
  std::string line;
  Util::readline( fin, line ); // line 1 (optional version string)
  if ( line.substr( 0, 10 ) == "# version " ) {
      version = atoi( line.substr( 10 ).c_str() );
      line.clear();
  }
  
  std::cerr << "ffield ver."<<version<<": from " << fname << std::endl;

    // line 1 (in old format, this is the first line)
    if ( line.empty() || line[0] == '#' ) {
        Util::readline_nocomments( fin, line );
    }
    if ( line != "target_frame" && line != "shape_operator" ) {
        // Log.error( "'%s': Expected 'target_frame' or 'shape_operator'. Got '%s'.", fname, line.c_str() );
        // spdlog::error("'{0:s}': Expected 'target_frame' or 'shape_operator'. Got '{1:s}'.", fname, line.c_str());
        return false;
    }
    bool target_frame = ( line == "target_frame" );
    bool crossfield_angles_loaded = false;

    int nframes = 0;
    Util::readline_nocomments( fin, line ); // line 2
    std::istringstream is( line );
    is >> nframes;
    if ( nframes != F.rows() ) {
        // Log.error( "'%s': Number of facets should be %d. Got %d.", fname, m->size_of_facets(), nframes );
        // spdlog::error("'{0:s}': Number of facets should be {1:d}. Got {2:d}.", fname, F.rows(), nframes);
        return false;
    }

    // ignore the final header line
    Util::readline_nocomments( fin, line ); // line 3

    // LOAD THE FRAMES
    int idx = 0;
    bool suppress_tangent_frame_check = false;
    for(int i=0;i<F.rows();i++){
      double k1, k2, vk1[3], vk2[3];
      Util::readline_nocomments( fin, line );
      std::istringstream is( line );
      is >> k1 >> k2
          >> vk1[0] >> vk1[1] >> vk1[2]
          >> vk2[0] >> vk2[1] >> vk2[2];
      if ( !fin || !is ) {
          // Log.error( "'%s' facet %d: unable to read tensor. File partially loaded", fname, idx );
          // spdlog::error("'{0:s}': facet {1:d}: unable to read tensor. File partially loaded", fname, idx);
          return false;
      }
      CVec3T<double> ev1( vk1[0], vk1[1], vk1[2] );
      CVec3T<double> ev2( vk2[0], vk2[1], vk2[2] );
      TangentFrame frame( k1, k2, ev1, ev2 );
      // load into target frames
      if ( load_into_frames ) {
        if ( !target_frame ) // loading shape operator => transform
          frame.inplace_perp();
      }
      // load into shape operator
      else {
        if ( target_frame ) // loading target frames => transform
          frame.inplace_inverse_perp();
        if ( !frame.is_orthogonal() )
          std::cerr << fname << ": loading non-orthogonal frame into shape operator!\n";
      }
    }

    // read the next token
    line.clear();
    if ( !Util::readline_nocomments( fin, line ) || ( line != "fixed" && line != "crossfield_angles" ) ) {
        // spdlog::error( "'{0:s}': Expected 'crossfield_angles' or 'fixed'. Got '{1:s}'. File partially loaded.", fname, line.c_str() );
        return false;
    }

    // LOAD THE CROSSFIELD ANGLES
    // output crossfield angles
    if ( line == "crossfield_angles" ) {
        line.clear();
        if ( !Util::readline_nocomments( fin, line ) ) {
            // spdlog::error( "'{0:s}': Missing list of crossfield angles. File partially loaded.", fname );
            return false;
        }
        std::istringstream is( line );
        for(int i=0;i<F.rows();i++)
          is >> angles(i);
        // for_each_facet( fi, m ) {
        //     is >> fi->crossfield_angle();
        // }
        if ( !is ) {
            // spdlog::error( "'{0:s}': Couldn't read crossfield angles. File partially loaded.", fname);
            return false;
        }

        // read in the next token
        line.clear();
        if ( !Util::readline_nocomments( fin, line ) || line != "fixed" ) {
            // spdlog::error( "'{0:s}': Expected 'fixed'. Got '{1:s}'. File partially loaded.", fname, line.c_str() );
            return false;
        }
        crossfield_angles_loaded = true;
    }
    
    // LOAD THE INDICES OF THE FIXED FACETS
    if ( line != "fixed" ) {
        // spdlog::error( "'{0:s}': Expected 'fixed'. Got '{1:s}'. File partially loaded.", fname, line.c_str() );
        return false;
    }
    bool failed = !Util::readline_nocomments( fin, line );
    if ( !failed ) {
        int nfixed = -1;
        std::istringstream is( line );
        is >> nfixed;
        // assert( 0 <= nfixed && nfixed <= (int)m->size_of_facets() ); // better error?
        std::set<int> fixed_facet_indices;
        for ( int i = 0; i < nfixed; ++i ) {
            int idx = -1;
            is >> idx;
            // assert( 0 <= idx && idx < (int)m->size_of_facets() );
            fixed_facet_indices.insert( idx );
        }
        if ( !is )
            failed = true;
        else {
            int idx = 0;
            for(int i=0;i<F.rows();i++){
              bool fixed = ( fixed_facet_indices.find( idx ) != fixed_facet_indices.end() );
              fixed_shapeop_f(idx) = fixed;
            }
            // for ( Facet_iterator fi = m->facets_begin(); fi != m->facets_end(); ++fi, ++idx ) {
            //     bool fixed = ( fixed_facet_indices.find( idx ) != fixed_facet_indices.end() );
            //     fi->fixed_shapeop() = fixed;
            // }
        }
    }
    if ( failed ) {
        // spdlog::error( "'{0:s}': Couldn't read indices of fixed facets. File partially loaded.", fname);
        return false;
    }

    // LOAD THE INDICES OF THE FIXED VERTICES
    line.clear();
    if ( !Util::readline_nocomments( fin, line ) ) {
        // spdlog::error( "'{0:s}': Couldn't read next line. File partially loaded.", fname);
        return false;
    }
    if ( line == "fixed_vertices" ) {
        // old versions don't have this
        bool failed = !Util::readline_nocomments( fin, line );
        if ( !failed ) {
            int nfixed = -1;
            std::istringstream is( line );
            is >> nfixed;
            // assert( 0 <= nfixed && nfixed <= (int)m->size_of_vertices() );
            std::set<int> fixed_vertex_indices;
            for ( int i = 0; i < nfixed; ++i ) {
                int idx = -1;
                is >> idx;
                // assert( 0 <= idx && idx < (int)m->size_of_vertices() );
                fixed_vertex_indices.insert( idx );
            }
            if ( !is )
                failed = true;
            else {
                int idx = 0;
                for(int i=0;i<V.rows();i++){
                  bool fixed = ( fixed_vertex_indices.find( idx ) != fixed_vertex_indices.end() );
                  fixed_shapeop_v(idx) = fixed;
                }
                // for ( Vertex_iterator vi = m->vertices_begin(); vi != m->vertices_end(); ++vi, ++idx ) {
                //     bool fixed = ( fixed_vertex_indices.find( idx ) != fixed_vertex_indices.end() );
                //     vi->fixed_shapeop() = fixed;
                // }
            }
        }
        if ( failed ) {
            // spdlog::error( "'{0:s}': Couldn't read indices of fixed vertices. File partially loaded.", fname );
            return false;
        }

        line.clear();
        if ( !Util::readline_nocomments( fin, line ) ) {
            // spdlog::error( "'{0:s}': Couldn't read next line. File partially loaded.", fname);
            return false;
        }
    }

    // LOAD THE MATCHINGS
    if ( ( line != "matchings" && line != "matchings_and_sharp" &&
           line != "MI_matchings" && line != "MI_matchings_and_sharp" ) ) {
        //spdlog::error( "'{0:s}': Expected one of 'matchings', 'matchings_and_sharp', 'MI_matchings',\n"
        //            "'MI_matchings_and_sharp'. Got '{1:s}'. File partially loaded.", fname, line.c_str() );
        return false;
    }
    bool load_sharp = ( line == "matchings_and_sharp" || line == "MI_matchings_and_sharp" );
    // old style matching was not equal to mixed integer matching
    bool old_style_matching = ( line.substr( 0, 3 ) != "MI_" );

    Eigen::MatrixXi period_jump;
    Eigen::MatrixXi period_jump_bit;
    period_jump.setZero(F.rows(),3);
    period_jump_bit.setZero(F.rows(),3); // whether pj of halfedge loaded

    bool bUninitializedMatchings = false;
    bool floating_point_matchings = false;
    idx = 0;
    Eigen::VectorXi first_edges(F.rows()); // stores the first edge of faces
    first_edges.setZero();
    for(int k=0;k<F.rows();k++){
      // read in the vertex indices and matchings
      int id[3], sh[3];
      double m[3];
      double ipart;
      Util::readline_nocomments( fin, line );
      std::istringstream is( line );
      is >> id[0] >> id[1] >> id[2] >> m[0] >> m[1] >> m[2];
      if ( load_sharp )
          is >> sh[0] >> sh[1] >> sh[2];
      if ( !fin || !is ) {
          // spdlog::error( "'{0:s}' facet %d: unable to read matchings. File partially loaded.", fname, k );
          return false;
      }
      // check if any floating-point matchings were loaded
      floating_point_matchings = floating_point_matchings ||
                                   ( ( modf( m[0], &ipart ) != 0.0 ) || ( modf( m[1], &ipart ) != 0.0 ) || ( modf( m[2], &ipart ) != 0.0 ) );
      // check for unintialized matchings
      const double BAD_MATCHING = 999999;
      for ( int i = 0; i < 3; ++i )
        bUninitializedMatchings |= m[i] == BAD_MATCHING;
      
      // locate first halfedge
      // [ // CGAL_For_all( hc, hc_end ) {            ]
      // [ //     if ( hc->vertex()->tag() == id[0] ) ]
      // [ //         break;                          ]
      // [ // }                                       ]
      int index = -1;
      for(int i=0;i<3;i++){
        if(id[0] == F(k,i)){
          index = i;
          break;
        }
      }
      first_edges(k) = index;
      // verify the facet is correct and load in the matchings
      int vnum = 0;
      for(int i=0;i<3;i++){
        int idx = (index + i)%3;
        if(F(k,idx) != id[i]){
          // spdlog::error( "'{0:s}' facet %d: Facet vertex indices ({1:d}, {2:d}, {3:d}) don't match mesh. File partially loaded.",
          //            fname, idx, id[0], id[1], id[2] );
          return false;
        }
        
        if ( old_style_matching ){
          period_jump(k,(idx+2)%3) = -m[i];
          // backward compatibility: the old matching was negative of the mixed integer matching
          ; // MIM::set_fmatching( hc, -m[vnum] );
        }else{
          period_jump(k,(idx+2)%3) = m[i];
          // mixed integer matching
          ; // MIM::set_fmatching( hc, m[vnum] );
        }
        
        if ( load_sharp )
          sharp(k,idx) = (sh[i] != 0); // hc->set_sharp_edge( sh[i] != 0 );
    
        int f1 = TT(k,(idx+2)%3);
        int e1 = TTi(k,(idx+2)%3);
        if(f1 != -1 && period_jump_bit(f1, e1) == 1 && fmod(period_jump(f1,e1) + period_jump(k,(idx+2)%3), 4.0) != 0)
          // spdlog::warn("{0:s} facet {1:d}-{2:d}: Inconsistent matchings {3:d} and {4:d} across edge.",
          //              fname,k,f1,period_jump(k,(idx+2)%3), period_jump(f1,e1));
        
        period_jump_bit(k,(idx+2)%3) = true; // mark matching as loaded
        
      }
      
    }
    // spdlog::info("matchings loaded!");
    
    // compute kappas
    std::vector<Eigen::MatrixXd> TPs;
    reference_frames(V,F,first_edges,TPs);
    Eigen::MatrixXd kappa;
    compute_kappa(V,F,first_edges,TPs,kappa);
    
    Eigen::MatrixXd R(F.rows(),3);
    R.setZero();
    for(int i=0;i<F.rows();i++){
      Eigen::Vector3i pj_col = period_jump.row(i).transpose();
      Eigen::Vector3d ka_col = kappa.row(i).transpose();
      Property pr(pj_col, ka_col, TPs[i], angles(i));
      props.push_back(pr);
    }
    find_cones(V,F,props,S);
    
// #define LOG
#ifdef LOG
    int count = 0;
    std::map<int,int> histogram;
    std::ofstream nf;
    std::string xx = fname;
    xx = xx + "_csv";
    nf.open(xx,std::ios_base::app);
    for(int i=0;i<S.rows();i++){
      if(S(i) != 0){
        // spdlog::info("{0:d}: {1:f}", i, S(i));
        nf<<"c "<<i+1<<" "<<2*(1-S(i))<<std::endl;
        count++;
      }
      if(S(i) < 0)
        histogram[4*(1-S(i))]++;
      if(S(i) > 0)
        histogram[4*(1-S(i))]++;
    }
    spdlog::info("#singularities loaded: {0:d}", count);
    nf.close();
    std::ofstream mf;
    std::pair<int,int> gb = statistics(V,F);
    mf.open("singularity_histogram.csv",std::ios_base::app);
    mf<<fname<<","<<gb.first<<","<<count<<",";
    for(auto key_val: histogram){
      int key = key_val.first;
      int val = key_val.second;
      mf<<val<<"("<<key<<")"<<",";
    }
    mf<<std::endl;
    
    mf.close();
    // exit(0);
#endif

  return true;
}

void save_ffield_hdf5(
  std::string fname,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  std::vector<Property>& props
){
  H5Easy::File hd_file(fname, H5Easy::File::Overwrite);

  // input list: V, F, pj, kappa, theta, TP, S
  Eigen::MatrixXd TP_stk(F.rows()*2, 3);
  Eigen::MatrixXi pj(F.rows(),3);
  Eigen::MatrixXd kappa(F.rows(),3);
  Eigen::VectorXd theta(F.rows());
  for(int i=0;i<F.rows();i++){
    TP_stk.row(i*2)    << props[i].TP.row(0);
    TP_stk.row(i*2+1)  << props[i].TP.row(1);
    pj.row(i)          << props[i].pj.transpose();
    kappa.row(i)       << props[i].kappa.transpose();
    theta.row(i)       << props[i].angle;
  }
  H5Easy::dump(hd_file, "TP", TP_stk);
  H5Easy::dump(hd_file, "kappa", kappa);
  H5Easy::dump(hd_file, "theta", theta);
  H5Easy::dump(hd_file, "pj", pj);
  H5Easy::dump(hd_file, "V", V);
  H5Easy::dump(hd_file, "F", F);
}

void load_ffield_hdf5(
  std::string fname,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  std::vector<Property>& props
){
  H5Easy::File hd_file(fname, H5Easy::File::ReadOnly);

  Eigen::MatrixXi pj    = H5Easy::load<Eigen::MatrixXi>(hd_file, "pj");
  Eigen::MatrixXd kappa = H5Easy::load<Eigen::MatrixXd>(hd_file, "kappa");
  Eigen::VectorXd theta = H5Easy::load<Eigen::VectorXd>(hd_file, "theta");
  Eigen::MatrixXd TP_stk= H5Easy::load<Eigen::MatrixXd>(hd_file, "TP");
  V                     = H5Easy::load<Eigen::MatrixXd>(hd_file, "V");
  F                     = H5Easy::load<Eigen::MatrixXi>(hd_file, "F");
  
  std::vector<Eigen::MatrixXd> TP;
  for(int i = 0; i < pj.rows(); i++){
    Eigen::MatrixXd tp(2, 3);
    tp << TP_stk.row(2*i), TP_stk.row(2*i+1);
    TP.push_back(tp);
  }

  Eigen::MatrixXd R;
  theta_to_vector(theta, TP, R);

  props.clear();
  for(int i = 0; i < pj.rows(); i++){
    Eigen::Vector3i _pj = pj.row(i);
    Eigen::Vector3d _kappa = kappa.row(i);
    Eigen::MatrixXd _TP = TP[i];
    double rd = theta(i);
    Property pr(_pj, _kappa, _TP, rd);
    props.push_back(pr);
  }
}