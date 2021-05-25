#ifndef AUXMATH_H
#define AUXMATH_H

#include <cassert>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <limits> // for std::numeric_limits
#include "Claussen.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif

const double SQRT2_2 = 0.7071067811865; // sqrt(2)/2 = 1/sqrt(2)

// TODO remove
const double EPSILON = FLT_EPSILON;

template <class T> T sqr( const T& a ) { return a*a; }
inline bool logical_xor(bool a, bool b) { return (a != b); }
inline bool logical_eq(bool a, bool b) { return (a == b); }

// computes the positive number k in 0..abs(b)-1 such that (k-a) mod abs(b) = 0
inline int pos_mod(int a, int b)
{
    b = abs(b);
    int k = a < 0 ? b - (-a % b) : a % b;
    if (k == b) k = 0;
    return k;
}

inline double pos_fmod(double x, double y)
{
    return (0 == y) ? x : x - y * floor(x/y);
}

inline bool is_between(double u, double u1, double u2)
{
    double mn = std::min(u1, u2);
    double mx = std::min(u1, u2);
    return (mn <= u && u <= mx);
}

inline double bary_from_u(double u, double u1, double u2)
{
    return (u1 == u2) ? std::numeric_limits<double>::quiet_NaN()
        : (u == u1) ? 0.0 // for numerical exactness
        : (u == u2) ? 1.0
        : (u - u1) / (u2 - u1);
}

// round to fixed number of fractional bits
inline double smallest_positive_fixed(int BINARY_FRAC_DIGS)
{
    // Don't exceed size of integer (otherwise the shift will be unpredictable)
    // Don't exceed the size of the mantissa in double.

    assert(0 <= BINARY_FRAC_DIGS && BINARY_FRAC_DIGS < (int)std::min(sizeof(unsigned long long)*8, (size_t)(52ul)));

    const double SCALEf = (double)((long long unsigned int)(1) << BINARY_FRAC_DIGS);

    return 1.0 / SCALEf;
}
inline double to_fixed(double v, int BINARY_FRAC_DIGS)
{
    const double SCALEf = smallest_positive_fixed(BINARY_FRAC_DIGS);
    return floor(v / SCALEf + 0.5) * SCALEf;
}

inline double l2s( const double angle ){
  const double s = 2*fabs(sin(angle));
  return s < 1e-40 ? -92.103403719761827360719658 : log(s);
}

inline double Lob( const double angle ){
  return 0 <= angle && angle <= M_PI ? claussen(2*angle)/2 : -1e10;
}

inline double PropLength( const double l, const double an, const double ad ){
  return l*(sin(an)/sin(ad));
}

inline double cot( const double x ){
  const double a = fabs(fmod(x+M_PI,2*M_PI)-M_PI)<1e-40 ? 1e-40 : x;
  return 1./tan(a);
}

template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

typedef std::complex<double> Cmplx;

#endif // AUXMATH_H

