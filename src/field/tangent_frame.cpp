
// #include "stdafx.h"

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cfloat>
#include <cmath>
#include <cstdio>

#include "tangent_frame.h"
#include "auxmath.h"

TangentFrame::TangentFrame()
{
    m_k[0] = 0;
    m_k[1] = 0;
}

TangentFrame::TangentFrame(double _k1, double _k2,
        CVec3T<double> _ev1, CVec3T<double> _ev2)
{
    // length must be either exactly 0 (uninitialized) or close to 1; and if I want ansio?
    //assert(fabs(lenSq(_ev1)) == 0 || fabs(lenSq(_ev1) - 1.0) < EPSILON);
    //assert(fabs(lenSq(_ev2)) == 0 || fabs(lenSq(_ev2) - 1.0) < EPSILON);
    initialize(_k1, _k2, _ev1, _ev2);
}

TangentFrame::TangentFrame(const TangentFrame &frame)
{
    initialize(frame.k1(), frame.k2(), frame.ev1(), frame.ev2());
}

inline void
TangentFrame::initialize(double _k1, double _k2,
        CVec3T<double> _ev1, CVec3T<double> _ev2)
{
    m_k[0] = _k1;
    m_k[1] = _k2;
    m_ev[0] = _ev1;
    m_ev[1] = _ev2;
}

void TangentFrame::setIdentity()
{
	initialize(1, 1, CVec3T<double>(1,0,0), CVec3T<double>(0,1,0));
}

bool
TangentFrame::initialized() const
{
    return ((m_ev[0].x() != 0.0 || m_ev[0].y() != 0.0 || m_ev[0].z() != 0.0) &&
            (m_ev[1].x() != 0.0 || m_ev[1].y() != 0.0 || m_ev[1].z() != 0.0));
}

// compute the normal
CVec3T<double>
TangentFrame::normal() const
{
    CVec3T<double> norm = cross(m_ev[0], m_ev[1]);
    safe_normalize(norm);
    return norm;
}

// returns if the frame represents a (nearly) symmetric frame.
bool
TangentFrame::is_orthogonal() const
{
    return fabs(dot(m_ev[0], m_ev[1])) < 1.0e-5;
}

bool TangentFrame::is_orthonormal() const
{
	return is_orthogonal() && abs(len(m_ev[0]) - 1) < 1e-5 && abs(len(m_ev[1]) - 1) < 1e-5;
}

// rotate the frame by some multiple of 90 degrees around its normal
TangentFrame
TangentFrame::rotate(int offset) const
{
    double rk1 = (offset % 2 == 0) ? m_k[0] : m_k[1];
    double rk2 = (offset % 2 == 0) ? m_k[1] : m_k[0]; 
    // This is to avoid issues with mod returning a negative number
    // when off is negative (apparently this is implementation
    // specific). Assume that off is negative, then, in modulo 4 math:
    //     (mod 4)  off = -(-off) = -(4 off - off) = -3 off.
    // If off is negative, (-3 off) is positive, which gives a positive
    // result according to the standard.
    int off1 = offset < 0 ? (-3*offset) % 4 : offset % 4;
    int off2 = (off1 + 1) % 4;
    return TangentFrame(rk1, rk2, dir(off1), dir(off2));
}

// returns the rotated frame field that is aligned against the supplied *unit* normal
TangentFrame
TangentFrame::aligned_frame(const CVec3T<double> &fromNormal, const CVec3T<double> &toNormal) const
{
    assert(fabs(lenSq(fromNormal) - 1.0) < EPSILON);
    assert(fabs(lenSq(toNormal) - 1.0) < EPSILON);
    assert(fabs(dot(fromNormal, ev1())) < EPSILON);
    assert(fabs(dot(fromNormal, ev2())) < EPSILON);

    TangentFrame aligned;
    Mat3T<double> R = Mat3T<double>::Rotation(fromNormal, toNormal);
    aligned = TangentFrame(m_k[0], m_k[1], R*m_ev[0], R*m_ev[1]);

    assert(fabs(dot(toNormal, normalized(aligned.ev1()))) < 1.0e-5);
    assert(fabs(dot(toNormal, normalized(aligned.ev2()))) < 1.0e-5);
    return aligned;
}

TangentFrame
TangentFrame::aligned_frame(const CVec3T<double> &toNormal) const
{
    return aligned_frame(this->normal(), toNormal);
}

// gives one of the four directions, based on index mod 4 (0 through 3)
CVec3T<double>
TangentFrame::dir(int index) const
{
    // This is to avoid issues with mod returning a negative number
    // when off is negative (apparently this is implementation
    // specific). When off is negative, then, in modulo 4 math:
    //     (mod 4)  off = -(-off) = -(4 off - off) = -3 off.
    // (-3 off) is positive, which gives a positive result with mod
    // according to the standard.
    int off = index < 0 ? (-3*index) % 4 : index % 4;
    if (0 <= off && off < 2)
        return m_ev[off];
    else if (2 <= off && off < 4)
        return -m_ev[off-2];
    else { // shouldn't happen
        assert(0);
        return CVec3T<double>(0.0, 0.0, 0.0);
    }
}

// how much one needs to rotate to get from one frame to the next
// ASSUMES that the frames are co-planar
// AM - Actually, the integrity of the computation is still good even if the
// frames are not co-planar. Still, I'm leaving the assertions in there.
int
TangentFrame::offset_to_coplanar_frame(const TangentFrame &t, CVec3T<double> common_edge, CVec3T<double> normal) const
{
    const TangentFrame &self = *this;

    // assert that the frames are co-planar
    assert(fabs(lenSq(normal) - 1.0) < EPSILON);
    assert(fabs(dot(normal, t.normal()) - 1.0) < EPSILON);
    assert(fabs(dot(normal, self.normal()) - 1.0) < EPSILON);

    // Note that K = [ev1,ev2]^t is the gradient of the parameterization.
    // Edge vector e maps to the domain edge vector K e.
    // Keep self constant and compare against the four possibilities in
    // CVec3T.
    Cmplx es( dot(self.ev1(), common_edge), dot(self.ev2(), common_edge) );
    Cmplx et( dot(   t.ev1(), common_edge), dot(   t.ev2(), common_edge) );

    // If no magnitude available, use directions to match
    if (self.k1() != 0 || self.k2() != 0) {
        es = Cmplx(self.k1()*real(es), self.k2()*imag(es));
    }
    //else {
    //    std::cout << "*** WARNING: tangent frame with 0 magnitude\n";
    //}
    if (t.k1() != 0 || t.k2() != 0) {
        et = Cmplx(t.k1()*real(et), t.k2()*imag(et));
    }
    //else {
    //    std::cout << "*** WARNING: tangent frame with 0 magnitude\n";
    //}
    const Cmplx I = Cmplx(0,1);
    int off = 0;
    double min = norm(es - et);
    for (int i = 0; i < 3; ++i)
    {
        et *= I; // rotate 90 degrees
        double d = norm(es - et);
        if (d < min) {
            off = i+1;
            min = d;
        }
    }

// XXX assertion not valid
//  double dt = dot(self.dir(off), t.ev1());
//  if (dt < 0) {
//      const CVec3T<double> z(0,0,1);
//      // align self and t to the 2D plane for easier debugging
//      TangentFrame a = self.aligned_frame(normal, z);
//      TangentFrame b = t.aligned_frame(normal, z);
//      double ab = dot(a.dir(off), b.ev1());
//      assert(dt >= 0);
//  }

    return off;
}


TangentFrame
TangentFrame::perp() const
{
    return TangentFrame(*this).inplace_perp();
}

TangentFrame
TangentFrame::inverse_perp() const
{
    return TangentFrame(*this).inplace_inverse_perp();
}

TangentFrame &
TangentFrame::inplace_perp()
{
    std::swap(m_k[0], m_k[1]);
    if (is_orthogonal()) { // avoid round-off error by swapping
        std::swap(m_ev[0], m_ev[1]);
        m_ev[1] = -m_ev[1];
    } else {
        CVec3T<double> norm = cross(m_ev[0], m_ev[1]);
        safe_normalize(norm);
        m_ev[0] = cross(norm, m_ev[0]);
        m_ev[1] = cross(norm, m_ev[1]);
    }
    return *this;
}

TangentFrame &
TangentFrame::inplace_inverse_perp()
{
    std::swap(m_k[0], m_k[1]);
    if (is_orthogonal()) { // avoid round-off error by swapping
        std::swap(m_ev[0], m_ev[1]);
        m_ev[0] = -m_ev[0];
    }
    else {
        CVec3T<double> norm = cross(m_ev[0], m_ev[1]);
        safe_normalize(norm);
        m_ev[0] = cross(m_ev[0], norm);
        m_ev[1] = cross(m_ev[1], norm);
    }
    return *this;
}


// Compute the frame orthogonal to norm with the prescribed angle from
// vector v0.
// (Not passing v0 by reference since it is normalized in the code.)
void
TangentFrame::from_angle(double angle, CVec3T<double> v0, CVec3T<double> norm)
{
    m_k[0] = 1;
    m_k[1] = 1;
    safe_normalize(norm);
    safe_normalize(v0);
    CVec3T<double> v1 = cross(norm, v0);
    assert(fabs(lenSq(v1)-1) < EPSILON);
    double c = cos(angle);
    double s = sin(angle);
    m_ev[0] =  c * v0 + s * v1;
    m_ev[1] = -s * v0 + c * v1;
}

// Convert the frame to an angle with respect to v0.
// (Not passing v0 by reference since it is normalized in the code.)
// Returns an angle in -Pi..Pi.
double
TangentFrame::to_angle(CVec3T<double> v0, CVec3T<double> norm) const
{
    return TangentFrame::to_angle(m_ev[0], v0, norm);
}

// (static function)
// Convert the vector cvec to an angle with respect to v0.
// (Not passing v0 by reference since it is normalized in the code.)
// Returns an angle in -Pi..Pi.
double
TangentFrame::to_angle(CVec3T<double> some_vec, CVec3T<double> v0, CVec3T<double> norm)
{
    safe_normalize(some_vec);
    safe_normalize(norm);
    safe_normalize(v0);
    assert(fabs(lenSq(v0)-1) < EPSILON);
    assert(fabs(lenSq(norm)-1) < EPSILON);
	bool bBad = false;
    if (fabs(dot(norm, v0)) >= EPSILON) {
		bBad = true;
        printf("WARNING: |dot(norm, v0)=%lg|  >=  EPSILON=%lg\n", dot(norm, v0), EPSILON);
    }
    //assert(fabs(dot(norm, v0)) < EPSILON);
    CVec3T<double> v1 = cross(norm, v0);
	if ( !bBad )
		assert(fabs(lenSq(v1)-1) < EPSILON);
    double x = dot(v0, some_vec);
    double y = dot(v1, some_vec);
    if (x == 0 && y == 0)
        return 0;
    else
        return atan2(y, x);
}


// Returns the least-square similarity transformation to an input Jacobian.
// (Ignoring the scale gives the least-square rotation.)
//
// INPUT
// * gradu and gradv are the u and v gradient vectors (the rows of the
//   Jacobian).
// * normal is the normal to both of them.
//
// OUTPUT
// * vectors rotu and rotv, which are the unit magnitude of the rows of the
//   least-square rotation (the rows of the target Jacobian if it were a
//   pure rotation).
// * scale_uv, which is what rotu and rotv should be scaled by to get the
//   closest rotation+scale transformation.
//
void
LSqRotation_from_Jacobian(const CVec3T<double> gradu, const CVec3T<double> gradv, const CVec3T<double> normal_in,
        CVec3T<double> &rotu, CVec3T<double> &rotv, double &scale_uv)
{
    // rotation that maps the triangle to the x-y plane
    // Mapping:
    //            Rotate       Param
    //    3D       to 2D   (shape usually
    //                       different)
    //    /\   ->   /\   ---->   /\.
    //   /__\      /__\         /__\.
    //    A         B            C
    //
    // Jacobian J: A -> C
    // Rotation R: A -> B
    //
    // So the Jacobian from B to C is  J R'. (' is the transpose)
    // Since the parameterization gradients are the rows of the
    // Jacobian, it's easier to work with J'.
    //
    // Say, LSQRot computes the least-square rotation of a 2x2 Jacobian.
    //
    // So, the least-square rotation LSQR to J is
    //   LSQR(J) = LSQRot(J R') R
    // The computation below is done in terms of the transpose.
    //   LSQR(J)' = R' LSQRot( (R J')' )'
    //
    //
	CVec3T<double> normal = normal_in;

	if ( 0 ) {
		// test normal
		CVec3T<double> fnormal2 = cross(gradu, gradv);
		safe_normalize(fnormal2);
		assert( abs(len(fnormal2) - 1) < 1e-7 );
		//double diff = fmin(linfty(normal - fnormal2), linfty(normal + fnormal2));
    double la = linfty(normal - fnormal2);
    double lb = linfty(normal + fnormal2);
    double diff = std::fmin(la, lb);
		if ( diff > EPSILON ) {
			; // Log.warn("normals diff %f > 1.0e-7", diff);
			//normal = fnormal2;
		}
	}

    // ZL no reason the normal isn't unit
    //assert( abs(len(normal) - 1) < 1e-7 );

    Mat3T<double> R = Mat3T<double>::Rotation(normal, CVec3T<double>(0,0,1));
    // ZL no reason the rotation matrix isn't a rotation matrix
    //assert( abs( ( R.transpose()*R - Mat3T<double>::Identity() ).frobnorm() ) < 1e-7 );

    CVec3T<double> gu2d = R * gradu;
    CVec3T<double> gv2d = R * gradv;
    // try test_uvs()
	static bool bWarn = false; // zl: warn only once (for the whole run..)
    if ( !bWarn && fabs(gu2d.z()) > 1.0e-5 ) {
        ; // Log.warn("[show once] LSqRotation_from_Jacobian(): fabs(gu2d.z()) = %f > 1.0e-5", fabs(gu2d.z()));
		bWarn = true;
	}
    if ( !bWarn && fabs(gv2d.z()) > 1.0e-5 ) {
        // Log.warn("fabs(gv2d.z()) = %f > 1.0e-5", fabs(gv2d.z()));
        ;
		bWarn = true;
	}
    // Convert it to the least-square rotation+scale
    //      [gu2d.x() gu2d.y()]  =  [a b]   ->  _1_ [a+d b-c]
    //      [gv2d.x() gv2d.y()]     [c d]        2  [c-b a+d]
    // For least-square rotation, need to normalize by the sqrt of the determinant.
    double c = 0.5*(gu2d.x()+gv2d.y());
    double s = 0.5*(gv2d.x()-gu2d.y());

    double sqrt_det = sqrt(c*c + s*s);
    if ( fabs(sqrt_det) < FLT_EPSILON ) {
        c = 1;
        s = 0;
    }  else {
        c /= sqrt_det;
        s /= sqrt_det;
    }
    CVec3T<double> rotu2d(c, -s, 0);
    CVec3T<double> rotv2d(s,  c, 0);

    Mat3T<double> Rt; Rt.setTranspose(R);
    rotu = Rt * rotu2d;
    rotv = Rt * rotv2d;
    scale_uv = sqrt_det;
}

