#ifndef _QUATERNION_HPP_
#define _QUATERNION_HPP_

#include "Vector3.hpp"
#include "Matrix3x3.hpp"

template<typename T>
class Quaternion {
public:
    Quaternion() : x(0), y(0), z(0), w(1) {}
    Quaternion(T x0, T y0, T z0, T w0) : x(x0), y(y0), z(z0), w(w0) {}

    Quaternion& operator+=(const Quaternion &r) { x += r.x; y += r.y; z += r.z; w += r.w; return *this; }
    Quaternion& operator-=(const Quaternion &r) { x -= r.x; y -= r.y; z -= r.z; w -= r.w; return *this; }
    Quaternion& operator*=(const Quaternion &r) { 
        T wn = w * r.w - x * r.x - y * r.y - z * r.z;
        T xn = w * r.x + r.w * x + y * r.z - z * r.y;
        T yn = w * r.y + r.w * y + z * r.x - x * r.z;
        T zn = w * r.z + r.w * z + x * r.y - y * r.x;

        w = wn;
        x = xn;
        y = yn;
        z = zn;

        return *this;
    }

    Quaternion operator+(const Quaternion &r) const { return Quaternion(*this) += r; }
    Quaternion operator-(const Quaternion &r) const { return Quaternion(*this) -= r; }
    Quaternion operator*(const Quaternion &r) const { return Quaternion(*this) *= r; }

    Quaternion& operator*=(const T r) { x *= r; y *= r; z *= r; w *= r; return *this; }
    Quaternion& operator/=(const T r) { x /= r; y /= r; z /= r; w /= r; return *this; }
    Quaternion operator*(const T r) const { return Quaternion(*this) *= r; }
    Quaternion operator/(const T r) const { return Quaternion(*this) /= r; }

    // bool operator==(const glm::quat &q) const {
    //     return (std::abs(x - q.x) < 1e-5 &&
    //             std::abs(y - q.y) < 1e-5 &&
    //             std::abs(z - q.z) < 1e-5 &&
    //             std::abs(w - q.w) < 1e-5);
    // }

    friend std::ostream& operator<<(std::ostream &out, const Quaternion q) {
        return (out << q.w << " " << q.x << " " << q.y << " " << q.z);
    }

    void normalize() {
        *this /= norm();
    }

    Real norm() {
        return std::sqrt(w * w + x * x + y * y + z * z);
    }

    Matrix3x3<T> to_rotation_mat() {
        return Matrix3x3<T>(
            1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w),
            2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w),
            2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)
        );
    }

    T x, y, z, w;

};

typedef Quaternion<Real> Quatf;

inline const Quatf operator*(const Real r, const Quatf &q) { return q * r; }

#endif /* _QUATERNION_HPP_ */
