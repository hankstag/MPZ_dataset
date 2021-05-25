// -*- Mode: c++ -*-
// Denis Zorin, version 8/2005 


#ifndef	__CVEC3T_H__
#define	__CVEC3T_H__

#include <iostream>
#include <cassert>
#include <cmath>

const double cvec3t_EPSILON = 1.0e-8;

// lines in this file are up to 160 char long

// list of non-member operations:
// IO: operator>>, istream& operator<<
// F = scalar class, V = vector class 
// when used as argument const ref is implied
//
// univariate: void normalize(V&)
//             F    lenSq (V)
//             F    len   (V)
//             F    l2    (V)
//             F    l1    (V)
//             F    linfty(V)
//             V    op-   (V)
//             V    dir   (V) 
//             int  largestAbsComp(V)
//             int  smallestAbsComp(V)
//
// bivariate:  V&   op=    (V&,V)
//             V&   op+=   (V&,V) 
//             V&   op-=   (V&,V)
//             V&   op*=   (V&,F)
//             V&   op/=   (V&,F)
//             V    op+    (V, V)
//             V    op-    (V, V)
//             V    compMult(V, V)
//             V    compDiv(V, V)
//             V    op*    (V, F) 
//             V    op*    (F, V) 
//             V    op/    (V, F) 
//             V    max    (V, V)
//             V    min    (V, V)
//             F    dot    (V, V)
//             F    dist   (V, V) 
//             F    angle  (V, V) 
//             V    cross  (V, V)
//             V    project(V, V)
//          bool    isCollinear(V,V)   
//
// trivariate: V    lerp   (V, V, F) 
//             V    rotate (V, V, F)
//             F    triProd(V, V, V) 

// there are no epsilons used anywhere: we do test for precise 0 before
// dividing, and for the range for  acos; in general one may need 
// to test if vectors are close to being equal etc
// accuracy of such tests is problem depedent and does not belong 
// here 



template <class F>
class CVec3T {
public:
  enum{ X = 0, Y = 1, Z = 2 };

  CVec3T( void )                              { v[X] = 0; v[Y] = 0; v[Z] = 0;}
  CVec3T( const CVec3T& c )                   { v[X] = c.v[X]; v[Y] = c.v[Y]; v[Z] = c.v[Z];} // copy constructor
  CVec3T( const F& a, const F& b, const F& c) { v[X] = a;      v[Y] = b;      v[Z] = c;     } // vector from 3 numbers
  explicit CVec3T( const F* a )                        { v[X] = a[X];   v[Y] = a[Y];   v[Z] = a[Z];  } // vector from an array

  // This breaks sorting vectors, because they get cast to pointers and then compared.
  // XXX AM - And that's why it SHOULDN'T BE ADDED
  //operator const F*( void ) const { return &v[X]; }  // cast a vector to an array

  ~CVec3T( void ) {}

  template <class G>
  CVec3T& operator=(const CVec3T<G>& c2) { v[0] =  c2(0); v[1]  = c2(1); v[2]  = c2(2); return (*this); }

  template <class G>
  bool operator==(const CVec3T<G>& c2) const {
      return v[0] == c2.v[0] && v[1] == c2.v[1] && v[2] == c2.v[2];
  }

  template <class G>
  bool operator!=(const CVec3T<G>& c2) const {
      return !(*this == c2);
  }


  // access components  
        F& x( void )       { return v[X]; }  
  const F& x( void ) const { return v[X]; }
        F& y( void )       { return v[Y]; }  
  const F& y( void ) const { return v[Y]; }
        F& z( void )       { return v[Z]; }  
  const F& z( void ) const { return v[Z]; }
  
        F& operator() (const unsigned int i)       { assert(i < 3); return v[i]; }
  const F& operator() (const unsigned int i) const { assert(i < 3); return v[i]; }
        F& operator[] (const unsigned int i)       { assert(i < 3); return v[i]; }
  const F& operator[] (const unsigned int i) const { assert(i < 3); return v[i]; }
  
 protected:
  F v[3];
};


// ----------------------------------------------------------------------------------------------------

// input/output
template <class F> std::ostream& operator<<(std::ostream& os, const CVec3T<F>& c) {  return os << c(0) << " " << c(1) << " " << c(2); }
template <class F> std::istream& operator>>(std::istream& is,       CVec3T<F>& c) {  is >> c(0);  is >> c(1);  is >> c(2);  return is;}


// ----------------------------------------------------------------------------------------------------
// univariate operations; most could be implemented as member functions, but do
// everything as non-member functions for consistency to simplify formula translation 

template <class F>  void      normalize(CVec3T<F>& c) { F l = len(c); assert(l != F(0)); c = c/l; }             // make the vector unit length

// saves on a sqrt() operation when the vector is unit
template <class F>  void safe_normalize(CVec3T<F>& c)     // make the vector unit length
{
    F l = lenSq(c);
    if (l != F(0) && std::abs(l-F(1)) >= 1.0e-16)
        c *= (F(1)/sqrt(l));
}
template <class F>  CVec3T<F> normalized(const CVec3T<F>& c)    // return the unit vector
{
    F l = lenSq(c);
    l = (l == F(0) || std::abs(l-F(1)) < 1.0e-16) ? F(1) : F(1)/sqrt(l);
    return c*l;
}

template <class F>  F         lenSq    (const CVec3T<F>& c) { return c(0)*c(0) + c(1)*c(1) + c(2)*c(2); }          // squared length
template <class F>  F         len      (const CVec3T<F>& c) { return sqrt(lenSq(c)); }                             // length
template <class F>  F         l2       (const CVec3T<F>& c) { return sqrt(lenSq(c)); }                             // synonym
template <class F>  F         l1       (const CVec3T<F>& c) { return fabs(c(0)) + fabs(c(1)) + fabs(c(2)); }       // l1 norm |x| + |y| + |z|
template <class F>  F         linfty   (const CVec3T<F>& c) { return std::max(std::max(fabs(c(0)),fabs(c(1))),fabs(c(2))); } // l infinity norm max(|x|,|y|,|z)|

template <class F>  CVec3T<F> operator-(const CVec3T<F>& c) { return CVec3T<F>( -c(0), -c(1), -c(2)); }            // unary minus
template <class F>  CVec3T<F> dir      (const CVec3T<F>& c) { F l = len(c); assert(l != F(0)); return c/l;}                // unit vector of the same direction


// find the index of the largest and smallest components of the vector
template <class F>  int largestAbsComp (const CVec3T<F>& c1)  {
  F a= fabs(c1(0)), b = fabs(c1(1)), c = fabs(c1(2));
  if      (a >= b && a >= c) return 0;
  else if (b >= a && b >= c) return 1;
  else return 2;
}

template <class F>  int smallestAbsComp (const CVec3T<F>& c1)  {
  F a= fabs(c1(0)), b = fabs(c1(1)), c = fabs(c1(2));
  if      (a <= b && a <= c) return 0;
  else if (b <= a && b <= c) return 1;
  else return 2;
}

// ----------------------------------------------------------------------------------------------------

// using multiple template params below is a mechanism to facilitate generation of
// additional operations combining different types of vectors
// e.g. as needed to support automatic differentition in an efficient way
// (so that constant vectors do not have to be cast to variable vectors)
// however, these templates are never instantiated automatically, because
// overload on return type is not allowed 
// so if we have a template with three parameters F,G,H, for
// operation F op(G g, H h) 
// it is insufficient to resolve F1 vf,vf1, vf2; vf = op(vf1, vf2).
// for this to be resolved, we need an explicit template specialization;
// the specializations for all parameters of the same type is return is 
// added at the end
// If we want to allow operations F0 vf0; F1 vf1; F2 vf2; 
// vf0 = op(vf1,vf2) (assuming op(f1,f2) is defined and can be cast to F0, 
// we need either a partial specialization 
// template <class G, class H> F0 op(G g, H h) { return op<F0,G,H>(g,h); } 
// or (more safely) a complete specialization 
// F0 op(F1 g, F2 h){ return  op<F0,F1,F2>(g,h) }

// ----------------------------------------------------------------------------------------------------
// bivariate operations

//assignment shortcuts

template <class F, class G, class H> 
inline CVec3T<F>& operator+=(CVec3T<G>& c1, const CVec3T<H>& c2) { c1(0) += c2(0); c1(1) += c2(1); c1(2) += c2(2); return c1; }
template <class F, class G, class H> 
inline CVec3T<F>& operator-=(CVec3T<G>& c1, const CVec3T<H>& c2) { c1(0) -= c2(0); c1(1) -= c2(1); c1(2) -= c2(2); return c1; }
template <class F, class G, class H> 
inline CVec3T<F>& operator*=(CVec3T<G>& c1, const H& s         ) { c1(0) *= s    ; c1(1) *= s    ; c1(2) *= s    ; return c1; }
template <class F, class G, class H> 
inline CVec3T<F>& operator/=(CVec3T<G>& c1, const H& s         ) { assert(s!=F(0)); 
                                                                   c1(0) /= s    ; c1(1) /= s    ; c1(2) /= s    ; return c1; }


// bivariate vector arithmetic operations +,-, component division and multiplication,
// multiplication by constant,division by constant, 

// a gcc problem: there is aliasing with a template operator defined in basic_string
// can be either avoided by requiring not to have using namespace std
// before this file is included, or using a different name

//template <class F, class G, class H> 
//inline CVec3T<F> operator+(const CVec3T<G>& c1, const CVec3T<H>& c2) { return CVec3T<F>( c1(0) + c2(0), c1(1) + c2(1), c1(2) + c2(2)); }
template <class F, class G, class H> 
inline CVec3T<F>    opplus(const CVec3T<G>& c1, const CVec3T<H>& c2) { return CVec3T<F>( c1(0) + c2(0), c1(1) + c2(1), c1(2) + c2(2)); }
template <class F, class G, class H> 
inline CVec3T<F> operator-(const CVec3T<G>& c1, const CVec3T<H>& c2) { return CVec3T<F>( c1(0) - c2(0), c1(1) - c2(1), c1(2) - c2(2)); }
template <class F, class G, class H> 
inline CVec3T<F> compMult (const CVec3T<G>& c1, const CVec3T<H>& c2) { return CVec3T<F>( c1(0) * c2(0), c1(1) * c2(1), c1(2) * c2(2)); }
template <class F, class G, class H> 
inline CVec3T<F> compDiv  (const CVec3T<G>& c1, const CVec3T<H>& c2) { assert(c2(0) != F(0) && c2(1) != F(0) && c2(2) != F(0));
                                                                       return CVec3T<F>( c1(0) / c2(0), c1(1) / c2(1), c1(2) / c2(2)); }
template <class F, class G, class H> 
inline CVec3T<F> operator*(const CVec3T<G>& c1, const        H&  s ) { return CVec3T<F>( c1(0) * s,     c1(1) * s    , c1(2) * s    ); }
template <class F, class G, class H> 
inline CVec3T<F> operator*(const        G&  s , const CVec3T<H>& c1) { return CVec3T<F>( s     * c1(0), s     * c1(1), s     * c1(2)); }
template <class F, class G, class H> 
inline CVec3T<F> operator/(const CVec3T<G>& c1, const        H&  s ) { assert(s != F(0) );
                                                                       return CVec3T<F>( c1(0) / s,     c1(1) / s    , c1(2) / s    ); }

// componentwise max/min
template <class F, class G, class H> 
inline CVec3T<F> max      (const CVec3T<G>& c1, const CVec3T<H>& c2) { return CVec3T<F>( max(c1(0),c2(0)), max(c1(1),c2(1)),max(c1(2),c2(2))); }
template <class F, class G, class H> 
inline CVec3T<F> min      (const CVec3T<G>& c1, const CVec3T<H>& c2) { return CVec3T<F>( min(c1(0),c2(0)), min(c1(1),c2(1)),min(c1(2),c2(2))); }

// geometric operations: dot product, distance len(u-v), angle between vectors,
// collinearity query
// cross product, projection (u dot v) v/|v|^2, 

template <class F, class G, class H> 
inline F dot ( const CVec3T<G>& c1, const CVec3T<H>& c2 ) { return ( c1.x()* c2.x() + c1.y() * c2.y() +  c1.z() * c2.z()); } 
template <class F, class G, class H> 
inline F dist( const CVec3T<G>& c1, const CVec3T<H>& c2 ) { return len(c1-c2); } // eliminate?


// angle between two vectors 0.. Pi; conservative version: asserts on
// argument to acos > 1, which may happen because of numerical innacuracy
// however calculations are rearranged to minimize probability of this
// not using truncation to [-1,1] because of issues with automatic differentiation
// probably the right approach to make this robust is to define a special
// version of acos for doubles which does the truncation

template <class F, class G, class H> 
inline F angle( const CVec3T<G>& c1, const CVec3T<H>& c2 ) { 
//F s   = lenSq(c1)*lenSq(c2); assert(s != F(0));
//F dtp = dot(c1,c2); 
//F dps = dtp*dtp/s; assert( dps <= F(1) && dps >= F(-1) );
//return acos(dtp > F(0)? sqrt( dps ):-sqrt(dps) );
  F s = len(cross(c1, c2));
  F c = dot(c1, c2);
  return atan2(s, c);
}


template <class F, class G, class H> 
inline CVec3T<F> cross( const CVec3T<G>& c1, const CVec3T<H>& c2 )  
{ return CVec3T<F>( c1.y() * c2.z() - c1.z() * c2.y(),  
		    c1.z() * c2.x() - c1.x() * c2.z(), 
		    c1.x() * c2.y() - c1.y() * c2.x() ); }


template <class F, class G, class H> 
inline CVec3T<F> project( const CVec3T<G>& c1, const CVec3T<H>& c2 ) { return c2*dot(c1,c2)/lenSq(c2); }

//template <class F, class G, class H> 
//inline bool isCollinear( const CVec3T<G>& c1, const CVec3T<H>& c2 ) { return (lenSq(cross(c1,c2)) == F(0) ); } // eliminate?

// ----------------------------------------------------------------------------------------------------

// trivariate geometric operations: linear interpolation (1-t)u + tv, rotation of u
// around v by angle a,  triple product u dot (v cross w)

template <class F, class G, class H, class K>
inline CVec3T<F> lerp      ( const CVec3T<G>& c1, const CVec3T<H>& c2, const K& s )        { return c1 + s*(c2-c1); } //eliminate?
template <class F, class G, class H, class K>
inline CVec3T<F> rotate    ( const CVec3T<G>& c1, const CVec3T<H>& c2, const K& s )        
{ 
  F c2l = len(c2);
  assert(c2l != F(0));
  CVec3T<F> unitc2 = c2/c2l;
  CVec3T<F> c1par = dot(c1,unitc2)*unitc2;
  CVec3T<F> c1perp = c1-c1par;
  return c1par + cos(s)*c1perp + sin(s)*cross(unitc2,c1perp);
}

 
template <class F, class G, class H, class K> 
inline         F tripleProd( const CVec3T<G>& c1, const CVec3T<H>& c2, const CVec3T<K>& c3){ return dot(c1,cross(c2,c3)); }

// Compute the counterclockwise angle from v1 to v2 when looking to the
// supplied unit normal orthogonal to both the vectors (i.e. when the normal
// is pointing at you).
template <class F>
inline F angle_between(const CVec3T<F> &v1, const CVec3T<F> &v2, const CVec3T<F> &normal)
{
    assert(std::abs(lenSq(normal) - 1.0) < cvec3t_EPSILON);
    assert(std::abs(dot(normal, v1)) < cvec3t_EPSILON);
    assert(std::abs(dot(normal, v2)) < cvec3t_EPSILON);
    F l1sq = lenSq(v1); if (l1sq == 0) return F(0);
    F l2sq = lenSq(v2); if (l2sq == 0) return F(0);
    F l1l2_inv = 1.0 / sqrt(l1sq*l2sq);
    F c = dot(v1, v2) * l1l2_inv;
    F s = dot(normal, cross(v1, v2)) * l1l2_inv;
    return atan2(s, c);
}

// ----------------------------------------------------------------------------------------------------
// specializations for all parameter types equal to the return type


template <class F>  inline CVec3T<F>& operator+= (      CVec3T<F>& c1, const CVec3T<F>& c2) { return  operator+=<F,F,F>( c1,c2); }
template <class F>  inline CVec3T<F>& operator-= (      CVec3T<F>& c1, const CVec3T<F>& c2) { return  operator-=<F,F,F>( c1,c2); }
template <class F>  inline CVec3T<F>& operator*= (      CVec3T<F>& c1, const F&         s ) { return  operator*=<F,F,F>( c1,s ); }
template <class F>  inline CVec3T<F>& operator/= (      CVec3T<F>& c1, const F&         s ) { return  operator/=<F,F,F>( c1,s ); }

template <class F>  inline CVec3T<F>  operator+  (const CVec3T<F>& c1, const CVec3T<F>& c2) { return      opplus<F,F,F>( c1,c2); }
template <class F>  inline CVec3T<F>  operator-  (const CVec3T<F>& c1, const CVec3T<F>& c2) { return   operator-<F,F,F>( c1,c2); }
template <class F>  inline CVec3T<F>  compMult   (const CVec3T<F>& c1, const CVec3T<F>& c2) { return   compMult <F,F,F>( c1,c2); }
template <class F>  inline CVec3T<F>  compDiv    (const CVec3T<F>& c1, const CVec3T<F>& c2) { return   compDiv  <F,F,F>( c1,c2); }
template <class F>  inline CVec3T<F>  operator*  (const CVec3T<F>& c1, const         F&  s) { return   operator*<F,F,F>( c1,s ); }
template <class F>  inline CVec3T<F>  operator*  (const F& s,          const CVec3T<F>& c1) { return   operator*<F,F,F>(  s,c1); }
template <class F>  inline CVec3T<F>  operator/  (const CVec3T<F>& c1, const         F&  s) { return   operator/<F,F,F>( c1,s ); }

template <class F>  inline CVec3T<F>  max        (const CVec3T<F>& c1, const CVec3T<F>& c2) { return         max<F,F,F>( c1,c2); }
template <class F>  inline CVec3T<F>  min        (const CVec3T<F>& c1, const CVec3T<F>& c2) { return         min<F,F,F>( c1,c2); }

template <class F>  inline F          dot        (const CVec3T<F>& c1, const CVec3T<F>& c2) { return         dot<F,F,F>( c1,c2); }
template <class F>  inline F          dist       (const CVec3T<F>& c1, const CVec3T<F>& c2) { return        dist<F,F,F>( c1,c2); }
template <class F>  inline F          angle      (const CVec3T<F>& c1, const CVec3T<F>& c2) { return       angle<F,F,F>( c1,c2); }
template <class F>  inline CVec3T<F>  cross      (const CVec3T<F>& c1, const CVec3T<F>& c2) { return       cross<F,F,F>( c1,c2); }
template <class F>  inline CVec3T<F>  project    (const CVec3T<F>& c1, const CVec3T<F>& c2) { return     project<F,F,F>( c1,c2); }
template <class F>  inline bool  isCollinear(const CVec3T<F>& c1, const CVec3T<F>& c2) { return isCollinear<F,F,F>( c1,c2); }

template <class F>  inline CVec3T<F>  lerp       (const CVec3T<F>& c1, const CVec3T<F>& c2, const F& s         ) { return lerp      <F,F,F,F>( c1,c2,s ); }
template <class F>  inline CVec3T<F>  rotate     (const CVec3T<F>& c1, const CVec3T<F>& c2, const F& s         ) { return rotate    <F,F,F,F>( c1,c2,s ); }
template <class F>  inline F          tripleProd (const CVec3T<F>& c1, const CVec3T<F>& c2, const CVec3T<F>& c3) { return tripleProd<F,F,F,F>( c1,c2,c3); }

#endif	/* __CVEC3T_H__ */

