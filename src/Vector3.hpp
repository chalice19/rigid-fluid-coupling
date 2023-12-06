// ----------------------------------------------------------------------------
// Vector3.hpp
//
//  Created on: 13 Feb 2020
//      Author: Kiwon Um
//        Mail: kiwon.um@telecom-paris.fr
//
// Description: 3D Vector (DO NOT DISTRIBUTE!)
//
// Copyright 2020-2023 Kiwon Um
//
// The copyright to the computer program(s) herein is the property of Kiwon Um,
// Telecom Paris, France. The program(s) may be used and/or copied only with
// the written permission of Kiwon Um or in accordance with the terms and
// conditions stipulated in the agreement/contract under which the program(s)
// have been supplied.
// ----------------------------------------------------------------------------

#ifndef _VECTOR3_HPP_
#define _VECTOR3_HPP_

#include <cassert>

#include "Vector.hpp"
#include "Matrix3x3.hpp"


template<typename T> class Matrix3x3;

template<typename T>
class Vector3 {
public:
  enum { D = 3 };

  typedef T ValueT;

  union {
    struct { T x; T y; T z; };
    struct { T i; T j; T k; };
    T v[D];
  };

  explicit Vector3(const T &value=0) : x(value), y(value), z(value) {}
  Vector3(const T &a, const T &b, const T &c=0): x(a), y(b), z(c) {}

  Vector2<T> get2d() {
    return Vector2<T>(x, y);
  }

  // assignment operators
  Vector3& operator+=(const Vector3 &r) { x+=r.x; y+=r.y; z+=r.z; return *this; }
  Vector3& operator-=(const Vector3 &r) { x-=r.x; y-=r.y; z-=r.z; return *this; }
  Vector3& operator*=(const Vector3 &r) { x*=r.x; y*=r.y; z*=r.z; return *this; }
  Vector3& operator/=(const Vector3 &r) { x/=r.x; y/=r.y; z/=r.z; return *this; }

  Vector3& operator+=(const T *r) { x+=r[0]; y+=r[1]; z+=r[2]; return *this; }
  Vector3& operator-=(const T *r) { x-=r[0]; y-=r[1]; z-=r[2]; return *this; }
  Vector3& operator*=(const T *r) { x*=r[0]; y*=r[1]; z*=r[2]; return *this; }
  Vector3& operator/=(const T *r) { x/=r[0]; y/=r[1]; z/=r[2]; return *this; }

  Vector3& operator+=(const T s) { x+=s; y+=s; z+=s; return *this; }
  Vector3& operator-=(const T s) { x-=s; y-=s; z-=s; return *this; }
  Vector3& operator*=(const T s) { x*=s; y*=s; z*=s; return *this; }
  Vector3& operator/=(const T s) {
    const T d=static_cast<T>(1)/s; return operator*=(d);
  }

  // unary operators
  Vector3 operator+() const { return *this; }
  Vector3 operator-() const { return Vector3(-x, -y, -z); }

  // binary operators
  Vector3 operator+(const Vector3 &r) const { return Vector3(*this)+=r; }
  Vector3 operator-(const Vector3 &r) const { return Vector3(*this)-=r; }
  Vector3 operator*(const Vector3 &r) const { return Vector3(*this)*=r; }
  Vector3 operator/(const Vector3 &r) const { return Vector3(*this)/=r; }

  Vector3 operator+(const T *r) const { return Vector3(*this)+=r; }
  Vector3 operator-(const T *r) const { return Vector3(*this)-=r; }
  Vector3 operator*(const T *r) const { return Vector3(*this)*=r; }
  Vector3 operator/(const T *r) const { return Vector3(*this)/=r; }

  Vector3 operator+(const T s) const { return Vector3(*this)+=s; }
  Vector3 operator-(const T s) const { return Vector3(*this)-=s; }
  Vector3 operator*(const T s) const { return Vector3(*this)*=s; }
  Vector3 operator/(const T s) const { return Vector3(*this)/=s; }

  // comparison operators
  bool operator==(const Vector3 &r) const {
    return ((x==r.x) && (y==r.y) && (z==r.z));
  }
  bool operator!=(const Vector3 &r) const { return !(*this==r); }
  bool operator<(const Vector3 &r) const {
    return (x!=r.x) ? x<r.x : (y!=r.y) ? y<r.y : z<r.z;
  }
  bool operator<=(const Vector3 &r) const {
    return (x!=r.x) ? x<=r.x : (y!=r.y) ? y<=r.y : z<=r.z;
  }

  // cast operator
  template<typename T2>
  operator Vector3<T2>() const {
    return
      Vector3<T2>(static_cast<T2>(x),static_cast<T2>(y),static_cast<T2>(z));
  }

  const T& operator[](const tIndex i) const { assert(i<D); return v[i]; }
  T& operator[](const tIndex i) {
    return const_cast<T &>(static_cast<const Vector3 &>(*this)[i]);
  }

  // special calculative functions

  Vector3& normalize() { return (x==0&&y==0&&z==0) ? (*this):(*this)/=length(); }
  Vector3 normalized() const { return Vector3(*this).normalize(); }

  T dotProduct(const Vector3 &r) const { return x*r.x + y*r.y + z*r.z; }
  Vector3 crossProduct(const Vector3 &r) const {
    return Vector3(y*r.z - z*r.y, z*r.x - x*r.z, x*r.y - y*r.x);
  }
  Matrix3x3<T> crossProductMatrix() const {
    return Matrix3x3<T>(0, -z, y,  z, 0, -x,  -y, x, 0);
  }
  Matrix3x3<T> outerProduct(const Vector3 &r) const {
    return Matrix3x3<T>(
      x*r.x, x*r.y, x*r.z, y*r.x, y*r.y, y*r.z, z*r.x, z*r.y, z*r.z);
  }

  T length() const { return std::sqrt(lengthSquare()); }
  T lengthSquare() const { return (x*x + y*y + z*z); }
  T distanceTo(const Vector3 &t) const { return (*this-t).length(); }
  T distanceSquareTo(const Vector3 &t) const { return (*this-t).lengthSquare(); }

  friend std::istream& operator>>(std::istream &in, Vector3 &vec) {
    return (in >> vec.x >> vec.y >> vec.z);
  }
  friend std::ostream& operator<<(std::ostream &out, const Vector3 &vec) {
    return (out << vec.x << " " << vec.y << " " << vec.z);
  }
};

typedef Vector3<Real> Vec3f;

inline const Vec3f operator*(const Real s, const Vec3f &r) { return r*s; }

#endif  /* _VECTOR3_HPP_ */
