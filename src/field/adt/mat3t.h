#ifndef _MATRIX3_H_
#define _MATRIX3_H_

#include <cstring>
#include "cvec3t.h"

template <class F> class Mat3T
{
public:
	// make this an identity matrix
	Mat3T&  setIdentity() {
		memset(_m, 0, 9 * sizeof(F));
		(*this)(0, 0) = (*this)(1, 1) = (*this)(2, 2) = 1.0f;
		return *this;
	}
	// returns an identity matrix  
	static Mat3T Identity() {
		return Mat3T().setIdentity();
	}
	Mat3T() {
		setIdentity();
	}
	// create from elements
	Mat3T(F m00, F m01, F m02,
		F m10, F m11, F m12,
		F m20, F m21, F m22) {
			Mat3T& M = *this;
			M(0, 0) = m00; M(0, 1) = m01; M(0, 2) = m02;
			M(1, 0) = m10; M(1, 1) = m11; M(1, 2) = m12;
			M(2, 0) = m20; M(2, 1) = m21; M(2, 2) = m22;
	}
	// create diagonal matrix from elements
	Mat3T(F m00, F m11, F m22) {
		Mat3T& M = *this;
		M(0, 0) = m00; M(0, 1) = F(0); M(0, 2) = F(0);
		M(1, 0) = F(0); M(1, 1) = m11; M(1, 2) = F(0);
		M(2, 0) = F(0); M(2, 1) = F(0); M(2, 2) = m22;
	}
	Mat3T(const Mat3T& mat) {
		memcpy(_m, mat._m, 9 * sizeof(F));
	}

	// cast to array
	// entries are in column-major order, i.e. first 3 elements are the first column
	operator F*()
	{
		return _m;
	}
	// cast to const array
	operator const F*() const {
		return _m;
	}
	// conversion from array to matrix: require to be explicit to avoid
	// unexpected implicit casts 
	explicit Mat3T(F* m) {
		memcpy(_m, m, 9 * sizeof(F));
	}

#if 0
	// set to matrix for translation transform
	Mat3T& setTranslation( const CVec3T<F>& trans ) { 
		setIdentity(); 
		(*this)(0,3) = trans.x();
		(*this)(1,3) = trans.y();
		(*this)(2,3) = trans.z();
		return *this;
	}

	static Mat3T Translation( const CVec3T<F>& trans ) {
		return Mat3T().setTranslation(trans);
	}
#endif

	// set to matrix for nonuniform scale transform
	Mat3T& setScale(const CVec3T<F>& scale) {
		setIdentity();
		(*this)(0, 0) = scale.x();
		(*this)(1, 1) = scale.y();
		(*this)(2, 2) = scale.z();
		return *this;
	}

	static Mat3T Scale(const CVec3T<F>& scale) {
		return Mat3T().setScale(scale);
	}
	// set to a rotation matrix for axis v and angle a; unlike OpenGL the angle is in radians, not degrees!
	Mat3T& setRotation(F a, const CVec3T<F>& v) {
		//Vec3T<F> u = v.dir();
		CVec3T<F> u = v / l2(v);
		F u1 = u.x();
		F u2 = u.y();
		F u3 = u.z();

		Mat3T U = Mat3T(u1*u1, u1*u2, u1*u3,
			u2*u1, u2*u2, u2*u3,
			u3*u1, u3*u2, u3*u3);
		Mat3T S = Mat3T(0, -u3, u2,
			u3, 0, -u1,
			-u2, u1, 0);

		(*this) = U + (Identity() - U) * cos(a) + S * sin(a);
#if 0
		(*this)(3,3) = 1.0;
#endif
		return *this;
	}

	// set to a rotation matrix that rotates vector 'from' to 'to' (both should already be normalized!)
	// The core idea is from: http://www.cs.brown.edu/people/jfh/papers/Moller-EBA-1999/paper.pdf
	// Moller & Hughes: Efficiently building a matrix to rotate one vector to
	// another (1999?)
	// The code is available at http://jgt.akpeters.com/papers/MollerHughes99/code.html
	Mat3T& setRotation(const CVec3T<F> &from, const CVec3T<F> &to) {
		assert(std::abs(lenSq(from) - 1.0) < 1.0e-8);
		assert(std::abs(lenSq(to) - 1.0) < 1.0e-8);
		Mat3T &R = *this;

		CVec3T<double> v = cross(from, to);
		double c = dot(from, to);

		if (std::abs(c) >= 0.99) // "from" and "to"-vector almost parallel
		{
			// vector most nearly orthogonal to "from"
			CVec3T<double> x = CVec3T<double>(
				(from(0) > 0.0) ? from(0) : -from(0),
				(from(1) > 0.0) ? from(1) : -from(1),
				(from(2) > 0.0) ? from(2) : -from(2));

			if (x(0) < x(1)) {
				if (x(0) < x(2)) {
					x(0) = 1.0; x(1) = x(2) = 0.0;
				}
				else {
					x(2) = 1.0; x(0) = x(1) = 0.0;
				}
			}
			else {
				if (x(1) < x(2)) {
					x(1) = 1.0; x(0) = x(2) = 0.0;
				}
				else {
					x(2) = 1.0; x(0) = x(1) = 0.0;
				}
			}

			CVec3T<double> u = x - from;
			CVec3T<double> v = x - to;
			double c1 = 2.0 / dot(u, u);
			double c2 = 2.0 / dot(v, v);
			double c3 = c1 * c2 * dot(u, v);

			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					R(i, j) = -c1 * u(i) * u(j)
						- c2 * v(i) * v(j)
						+ c3 * v(i) * u(j);
				}
				R(i, i) += 1.0;
			}
		}
		else  // the most common case, unless "from"="to", or "from"=-"to"
		{
			double h = 1.0 / (1.0 + c);
			R(0, 0) = c + h * v(0) * v(0);
			R(0, 1) = h * v(0) * v(1) - v(2);
			R(0, 2) = h * v(0) * v(2) + v(1);

			R(1, 0) = h * v(0) * v(1) + v(2);
			R(1, 1) = c + h * v(1) * v(1);
			R(1, 2) = h * v(1) * v(2) - v(0);

			R(2, 0) = h * v(0) * v(2) - v(1);
			R(2, 1) = h * v(1) * v(2) + v(0);
			R(2, 2) = c + h * v(2) * v(2);
		}
		return R;
	}

	static Mat3T Rotation(F a, const CVec3T<F>& v) {
		return Mat3T().setRotation(a, v);
	}
	static Mat3T Rotation(const CVec3T<F> &from, const CVec3T<F> &to) {
		return Mat3T().setRotation(from, to);
	}

	//// set to a rotation matrix that rotates u into v
	//Mat3T& setRotation(CVec3T<F> u, CVec3T<F> v) {
	//  normalize(v);
	//  CVec3T<F> w   = cross(u, v); 
	//  if (l2(w)==0)
	//  {
	//    setIdentity();
	//  } else {
	//    normalize(w);
	//    CVec3T<F> uP  = cross(w, u); normalize(uP);
	//    F phi = atan2( dot(uP, v), dot(u, v) );
	//    setRotation(phi, w);
	//  }
	//  return *this;
	//}

	//static Mat3T Rotation( const CVec3T<F>& u, const CVec3T<F>& v)
	//{
	//  return Mat3T().setRotation(u, v);
	//}

	Mat3T tangentPlanePrincipalCurvatures(CVec3T<F> normal, CVec3T<F> edge)
	{
		CVec3T<F> t1 = edge - normal*dot(normal, edge); normalize(t1);
		CVec3T<F> t2 = cross(normal, t1);   normalize(t2);

		F a = dot(t1, (*this)*t1);
		F b = dot(t1, (*this)*t2);
		F c = dot(t2, (*this)*t2);

		F discr = sqrt((a - c)*(a - c) + 4 * b*b);
		F lambda1 = 0.5 * ((a + c) + discr);
		F lambda2 = 0.5 * ((a + c) - discr);

		CVec3T<F> v1, v2;

		if (b == 0) // the above formulas are wrong for b=0, but that 
		{
			v1 = t1;
			v2 = t2;
		}
		else {
			v1 = (-b * t1) + (a - lambda1)*t2; normalize(v1);
			v2 = (-b * t1) + (a - lambda2)*t2; normalize(v2);
		}

		v1 = lambda1 * v1;
		v2 = lambda2 * v2;

		return Mat3T(
			v1.x(), v1.y(), v1.z(),
			v2.x(), v2.y(), v2.z(),
			normal.x(), normal.y(), normal.z()
			);
	}

	// (i,j) element access
	// this is the only function for which the internal storage order matters
	F operator()(int i, int j) const {
		return _m[3 * j + i];
	}

	F& operator()(int i, int j)  {
		return _m[3 * j + i];
	}

	Mat3T& operator=(const Mat3T& mat) {
		memcpy(_m, mat._m, 9 * sizeof(F));
		return *this;
	}

	// extract column
	CVec3T<F> col(int i) const {
		return CVec3T<F>((*this)(0, i), (*this)(1, i), (*this)(2, i));
	}
	// extract row
	CVec3T<F> row(int i) const {
		return CVec3T<F>((*this)(i, 0), (*this)(i, 1), (*this)(i, 2));
	}

	// set row / col
	void set_row(const CVec3T<F> &v, int r) {
		for (int i = 0; i < 3; i++)
			(*this)(r, i) = v(i);
	}

	void set_rows(const CVec3T<F> &v0, const CVec3T<F> &v1, const CVec3T<F> &v2) {
		set_row(v0, 0);
		set_row(v1, 1);
		set_row(v2, 2);
	}

	void set_col(const CVec3T<F> &v, int c) {
		for (int i = 0; i < 3; i++)
			(*this)(i, c) = v(i);
	}

	void set_cols(const CVec3T<F> &v0, const CVec3T<F> &v1, const CVec3T<F> &v2) {
		set_col(v0, 0);
		set_col(v1, 1);
		set_col(v2, 2);
	}

	// set the matrix to the inverse of M; return true if inverse exists  
	bool setInverse(const Mat3T& M) {
		Mat3T A;
		int i, j, k;
		F V;

		A = M;
		setIdentity();


		for (i = 0; i < 3; i++) {
			V = A(i, i);              /* Find the new pivot. */
			k = i;
			for (j = i + 1; j < 3; j++)
			if (fabs(A(j, i)) > fabs(V)) {
				/* Find maximum on col i, row i+1..n */
				V = A(j, i);
				k = j;
			}
			j = k;


			F tmp;
			if (i != j)
			for (k = 0; k < 3; k++) {
				tmp = A(i, k); A(i, k) = A(j, k); A(j, k) = tmp;
				tmp = (*this)(i, k); (*this)(i, k) = (*this)(j, k); (*this)(j, k) = tmp;
			}


			for (j = i + 1; j < 3; j++) {   /* Eliminate col i from row i+1..n. */
				if (A(j, i) != 0) {
					V = A(j, i) / A(i, i);

					for (k = 0; k < 3; k++) {
						A(j, k) -= V * A(i, k);
						(*this)(j, k) -= V * (*this)(i, k);
					}
				}

			}
		}

		for (i = 2; i >= 0; i--) {             /* Back Substitution. */
			if (A(i, i) == 0)
				return false;             /* Error. */

			for (j = 0; j < i; j++) {   /* Eliminate col i from row 1..i-1. */
				V = A(j, i) / A(i, i);

				for (k = 0; k < 3; k++) {
					/* A[j][k] -= V * A[i][k]; */
					(*this)(j, k) -= V * (*this)(i, k);
				}
			}
		}

		for (i = 0; i < 3; i++)        /* Normalize the inverse Matrix. */
		for (j = 0; j < 3; j++)
			(*this)(i, j) /= A(i, i);

		return true;
	}

	// return inverse; WARNING: returns garbage if the matrix is not invertible
	Mat3T inverse() const {
		Mat3T M; M.setInverse(*this);
		return M;
	}
	// set the matrix to transpose of M
	Mat3T& setTranspose(const Mat3T& M)  {
		for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			(*this)(i, j) = M(j, i);
		return (*this);
	}

	Mat3T transpose() const {
		Mat3T M; M.setTranspose(*this);
		return M;
	}

	// matrix vector product Mv, vector regarded as a column  
	template <typename Vector>
	Vector mult(const Vector& v) const
	{
		const Mat3T& M = *this;

		return Vector(
			v.x() * M(0, 0) + v.y() * M(0, 1) + v.z() * M(0, 2),
			v.x() * M(1, 0) + v.y() * M(1, 1) + v.z() * M(1, 2),
			v.x() * M(2, 0) + v.y() * M(2, 1) + v.z() * M(2, 2)
			);
	}
	CVec3T<F> operator* (const CVec3T<F>& v) const {
		return this->mult(v);
	}

	// scale a matrix
	Mat3T operator* (F s) const {
		const Mat3T& M = *this;
		return Mat3T(s*M(0, 0), s*M(0, 1), s*M(0, 2),
			s*M(1, 0), s*M(1, 1), s*M(1, 2),
			s*M(2, 0), s*M(2, 1), s*M(2, 2)
			);
	}
	// scale a matrix
	Mat3T operator*=(F s) {
		Mat3T& M = *this;
		M(0, 0) *= s; M(0, 1) *= s; M(0, 2) *= s;
		M(1, 0) *= s; M(1, 1) *= s; M(1, 2) *= s;
		M(2, 0) *= s; M(2, 1) *= s; M(2, 2) *= s;
		return M;
	}
	// multiply two matrices  
	Mat3T operator* (const Mat3T& M2) const {
		const Mat3T& M1 = *this;
		return Mat3T(
			M1(0, 0)*M2(0, 0) + M1(0, 1)*M2(1, 0) + M1(0, 2)*M2(2, 0),
			M1(0, 0)*M2(0, 1) + M1(0, 1)*M2(1, 1) + M1(0, 2)*M2(2, 1),
			M1(0, 0)*M2(0, 2) + M1(0, 1)*M2(1, 2) + M1(0, 2)*M2(2, 2),


			M1(1, 0)*M2(0, 0) + M1(1, 1)*M2(1, 0) + M1(1, 2)*M2(2, 0),
			M1(1, 0)*M2(0, 1) + M1(1, 1)*M2(1, 1) + M1(1, 2)*M2(2, 1),
			M1(1, 0)*M2(0, 2) + M1(1, 1)*M2(1, 2) + M1(1, 2)*M2(2, 2),

			M1(2, 0)*M2(0, 0) + M1(2, 1)*M2(1, 0) + M1(2, 2)*M2(2, 0),
			M1(2, 0)*M2(0, 1) + M1(2, 1)*M2(1, 1) + M1(2, 2)*M2(2, 1),
			M1(2, 0)*M2(0, 2) + M1(2, 1)*M2(1, 2) + M1(2, 2)*M2(2, 2)
			);
	}

	Mat3T &operator+=(const Mat3T& M2) {
		Mat3T& M1 = *this;
		M1(0, 0) += M2(0, 0); M1(0, 1) += M2(0, 1); M1(0, 2) += M2(0, 2);
		M1(1, 0) += M2(1, 0); M1(1, 1) += M2(1, 1); M1(1, 2) += M2(1, 2);
		M1(2, 0) += M2(2, 0); M1(2, 1) += M2(2, 1); M1(2, 2) += M2(2, 2);
		return M1;
	}

	Mat3T &operator-=(const Mat3T& M2) {
		Mat3T& M1 = *this;
		M1(0, 0) -= M2(0, 0); M1(0, 1) -= M2(0, 1); M1(0, 2) -= M2(0, 2);
		M1(1, 0) -= M2(1, 0); M1(1, 1) -= M2(1, 1); M1(1, 2) -= M2(1, 2);
		M1(2, 0) -= M2(2, 0); M1(2, 1) -= M2(2, 1); M1(2, 2) -= M2(2, 2);
		return M1;
	}

	Mat3T operator+ (const Mat3T& M2) const {
		const Mat3T& M1 = *this;
		return Mat3T(
			M1(0, 0) + M2(0, 0), M1(0, 1) + M2(0, 1), M1(0, 2) + M2(0, 2),
			M1(1, 0) + M2(1, 0), M1(1, 1) + M2(1, 1), M1(1, 2) + M2(1, 2),
			M1(2, 0) + M2(2, 0), M1(2, 1) + M2(2, 1), M1(2, 2) + M2(2, 2));
	}

	Mat3T operator- (const Mat3T& M2) const {
		const Mat3T& M1 = *this;
		return Mat3T(
			M1(0, 0) - M2(0, 0), M1(0, 1) - M2(0, 1), M1(0, 2) - M2(0, 2),
			M1(1, 0) - M2(1, 0), M1(1, 1) - M2(1, 1), M1(1, 2) - M2(1, 2),
			M1(2, 0) - M2(2, 0), M1(2, 1) - M2(2, 1), M1(2, 2) - M2(2, 2)
			);

	}

	// unary version
	Mat3T operator- () const {
		const Mat3T& M1 = *this;
		return Mat3T(
			-M1(0, 0), -M1(0, 1), -M1(0, 2),
			-M1(1, 0), -M1(1, 1), -M1(1, 2),
			-M1(2, 0), -M1(2, 1), -M1(2, 2)
			);
	}

	// Frobenius norm, i.e. sqrt of the sum of the squares of the entries also defined as trace(M^2)
	F frobnorm() const {
		F fn = 0;
		for (int i = 0; i < 9; i++) fn += _m[i] * _m[i];
		return sqrt(fn);
	}

	// Determinant
	F det() const {
		const Mat3T& M = *this;
		return   M(0, 0) * (M(1, 1) * M(2, 2) - M(1, 2) * M(2, 1))
			- M(0, 1) * (M(1, 0) * M(2, 2) - M(1, 2) * M(2, 0))
			+ M(0, 2) * (M(1, 0) * M(2, 1) - M(1, 1) * M(2, 0));
	}

private:
	// data stored in column major order for OpenGL compatibility, i.e. first for elements are the first column of the matrix
	F _m[9];
};
template <class F> inline Mat3T<F> operator*(F s, const Mat3T<F>& M) {
	return M*s;
}
template <class F> inline std::ostream& operator<< (std::ostream& os, const Mat3T<F>& M) {
	os << "[ " << M(0, 0) << " " << M(0, 1) << " " << M(0, 2) << " " << "; ";
	os << M(1, 0) << " " << M(1, 1) << " " << M(1, 2) << " " << "; ";
	os << M(2, 0) << " " << M(2, 1) << " " << M(2, 2) << " " << "; ";
	os << "] ";
	return os;
}

template <class F>
inline Mat3T<F> outerProd(const CVec3T<F>& c1, const CVec3T<F>& c2)
{
	return Mat3T<F>(c1.x() * c2.x(), c1.x() * c2.y(), c1.x() * c2.z(),
		c1.y() * c2.x(), c1.y() * c2.y(), c1.y() * c2.z(),
		c1.z() * c2.x(), c1.z() * c2.y(), c1.z() * c2.z());
}

#endif
