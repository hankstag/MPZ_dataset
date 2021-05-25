#ifndef TANGENT_FRAME_H
#define TANGENT_FRAME_H

#include <complex>

#include "adt/mat3t.h"


// more general -- allows for a non-orthogonal frame
class TangentFrame {
protected:
    double m_k[2];          // magnitudes
    CVec3T<double> m_ev[2]; // unit vectors

public:
    TangentFrame();
    TangentFrame(double _k1, double _k2,
            CVec3T<double> _ev1, CVec3T<double> _ev2);
    TangentFrame(const TangentFrame &frame);

    double &k1() { return m_k[0]; }
    double &k2() { return m_k[1]; }
    double k1() const { return m_k[0]; }
    double k2() const { return m_k[1]; }
    CVec3T<double> &ev1() { return m_ev[0]; }
    CVec3T<double> &ev2() { return m_ev[1]; }
    const CVec3T<double> &ev1() const { return m_ev[0]; }
    const CVec3T<double> &ev2() const { return m_ev[1]; }
    CVec3T<double> vec1() const { return m_ev[0] * m_k[0]; }
    CVec3T<double> vec2() const { return m_ev[1] * m_k[1]; }

    // initialize the data
    void initialize(double _k1, double _k2,
        CVec3T<double> _ev1, CVec3T<double> _ev2);

    // considered initialized if both ev1 and ev2 have non-zero magnitude
    bool initialized() const;

	void setIdentity();

    // compute the normal
    CVec3T<double> normal() const;

    // returns whether or not the frame is orthogonal
    bool is_orthogonal() const;
	bool is_orthonormal() const;

    // rotate the frame by some multiple of 90 degrees around its normal
    TangentFrame rotate(int offset) const;

    // returns the rotated frame field that is aligned against the supplied *unit* normal
    TangentFrame aligned_frame(const CVec3T<double> &fromNormal, const CVec3T<double> &toNormal) const;
    TangentFrame aligned_frame(const CVec3T<double> &toNormal) const;

    // gives one of the four directions, based on index mod 4 (0 through 3)
    CVec3T<double> dir(int index) const;

    // how much one needs to rotate to get from one frame to the next
    // ASSUMES that the frames are co-planar
    int offset_to_coplanar_frame(const TangentFrame &t, CVec3T<double> common_edge, CVec3T<double> normal) const;

    // how much one needs to rotate to get from one frame to the next
    // - computation is done after alignment w.r.t. the supplied unit normal
    int offset_to_frame(const TangentFrame &t, CVec3T<double> common_edge, CVec3T<double> normal) const;

    // the perpendicular operation allows one to go from the tensor field
    // to the shape operator (perp) and vice versa (inverse_perp). In both
    // cases, the frame's orientation is preserved so that the cross product
    // of the two vectors is the same.
    TangentFrame perp() const;
    TangentFrame inverse_perp() const;
    TangentFrame &inplace_perp();
    TangentFrame &inplace_inverse_perp();

    // Compute the frame orthogonal to norm with the prescribed angle from
    // vector v0.
    void from_angle(double angle, CVec3T<double> v0, CVec3T<double> norm);
    // Convert the frame to an angle with respect to v0.
    double to_angle(const CVec3T<double> v0, CVec3T<double> norm) const;

    // Convert the vector some_vec to an angle with respect to v0.
    static double to_angle(const CVec3T<double> some_vec, const CVec3T<double> v0, CVec3T<double> norm);

    // disallow some functions (for now) to be safe
};



// Returns the least-square similarity transformation to an input Jacobian
// (a tangent frame can be interpreted as such).
// (Tgnoring the scale gives the least-square rotation.)
//
// INPUT
// * gradu and gradv are the u and v gradient vectors (the rows of the
//   Jacobian).
// * norm is the normal to both of them.
//
// OUTPUT
// * vectors rotu and rotv, which are the unit magnitude of the rows of the
//   least-square rotation (the rows of the target Jacobian if it were a
//   pure rotation).
// * scale_uv, which is what rotu and rotv should be scaled by to get the
//   closest rotation+scale transformation.
//
void LSqRotation_from_Jacobian(const CVec3T<double> gradu, const CVec3T<double> gradv, const CVec3T<double> norm,
            CVec3T<double> &rotu, CVec3T<double> &rotv, double &scale_uv);


#endif // TANGENT_FRAME_H

