// ----------------------------------------------------------------------------
// RigidSolver.hpp
//
//  Created on: 18 Dec 2020
//      Author: Kiwon Um
//        Mail: kiwon.um@telecom-paris.fr
//
// Description: Simple Rigid Body Solver (DO NOT DISTRIBUTE!)
//
// Copyright 2020-2023 Kiwon Um
//
// The copyright to the computer program(s) herein is the property of Kiwon Um,
// Telecom Paris, France. The program(s) may be used and/or copied only with
// the written permission of Kiwon Um or in accordance with the terms and
// conditions stipulated in the agreement/contract under which the program(s)
// have been supplied.
// ----------------------------------------------------------------------------

#ifndef _RIGIDSOLVER_HPP_
#define _RIGIDSOLVER_HPP_

#include "Vector3.hpp"
#include "Matrix3x3.hpp"
#include "Quaternion.hpp"

struct BodyAttributes {
  BodyAttributes() :
    X(0, 0, 0), R(Mat3f::I()), P(0, 0, 0), L(0, 0, 0),
    V(0, 0, 0), omega(0, 0, 0), F(0, 0, 0), tau(0, 0, 0), q(Quatf()) {}

  tIndex vertexCount() { return vpos.size(); }
  const Vec3f &position(const tIndex i) const { return vpos[i]; }
  const Vec2f &position2d(const tIndex i) const { return vpos2d[i]; }
  const float &get_color(const tIndex i) const { return color[i]; }

  // glm::mat4 worldMat() const
  // {
  //   return glm::mat4(           // column-major
  //     R(0,0), R(1,0), R(2,0), 0,
  //     R(0,1), R(1,1), R(2,1), 0,
  //     R(0,2), R(1,2), R(2,2), 0,
  //     X[0],   X[1],   X[2],   1);
  // }

  Real M;                       // mass
  Mat3f I0, I0inv;              // inertia tensor and its inverse in body space
  Mat3f Iinv;                   // inverse of inertia tensor

  // rigid body state
  Vec3f X;                      // position
  Mat3f R;                      // rotation
  Vec3f P;                      // linear momentum
  Vec3f L;                      // angular momentum
  Quatf q;                      // rotation represented by quaternion

  // auxiliary quantities
  Vec3f V;                      // linear velocity
  Vec3f omega;                  // angular velocity

  // force and torque
  Vec3f F;                      // force
  Vec3f tau;                    // torque

  // mesh's vertices in body space
  std::vector<float> color;
  std::vector<Vec3f> vpos;
  std::vector<Vec2f> vpos2d;

  // Vec3f center_of_mass;
};

class Box : public BodyAttributes {
public:
  explicit Box(
    const Real w = 1.0, const Real h = 1.0, const Real d = 1.0, const Real m = 10.0,
    const Vec3f v0 = Vec3f(0), const Vec3f omega0 = Vec3f(0)) :
    width(w), height(h), depth(d)
  {
    V = v0;                     // initial velocity
    omega = omega0;             // initial angular velocity
    M = m;

    Real M_div_12 = M / 12.0f;

    I0 = Mat3f(M_div_12 * (height * height + depth * depth), 0.0f, 0.0f,
              0.0f, M_div_12 * (width * width + depth * depth), 0.0f,
              0.0f, 0.0f, M_div_12 * (width * width + height * height));
    I0inv = I0.invert();
    Iinv = R * I0inv * R.transposed();

    // vertices data
    vpos.push_back(Vec3f(-0.5*w, -0.5*h, 0.0));
    vpos.push_back(Vec3f( 0.5*w, -0.5*h, 0.0));
    vpos.push_back(Vec3f( 0.5*w,  0.5*h, 0.0));
    vpos.push_back(Vec3f(-0.5*w,  0.5*h, 0.0));

    // center_of_mass = Vec3f(0.0);
    // for (int i = 0; i < vertexCount(); ++i) {
    //   center_of_mass += vpos[i];
    // }
    // center_of_mass /= vertexCount();

    color = std::vector<float>(vertexCount() * 4, 1.0); // RGBA

    for (int i = 0; i < vertexCount(); ++i) {
      color[i * 4 + 0] = 1.0;
      color[i * 4 + 1] = 0.6;
      color[i * 4 + 2] = 0.6;
    }
  }

  // rigid body property
  Real width, height, depth;
};

class RigidSolver {
public:
  explicit RigidSolver(
    BodyAttributes *body0=nullptr, const Vec3f g=Vec3f(0, 0, 0)) :
    body(body0), _g(g), _step(0), _sim_t(0) {}

  void init(BodyAttributes *body0, Vec3f pos)
  {
    body = body0;
    body->vpos2d = std::vector<Vec2f>(body->vertexCount());
    update_pos2d();

    body->X = pos;

    // body->F = _g * body->M;
    body->F = Vec3f(0.0);
    body->tau = Vec3f(0.0f);

    _step = 0;
    _sim_t = 0;
  }

  Vec2f point_velocity(Vec2f p) const {
    Vec3f p3d = Vec3f(p.x, p.y, 0.0);
    Vec3f v3d = body->V + (body->omega).crossProduct(p3d - body->X);
    return v3d.get2d();
  }

  Vec2f translate(Vec2f p) {
    Vec3f p3d = Vec3f(p.x, p.y, 0.0);
    p3d = body->R * p3d + body->X;
    return p3d.get2d();
  }

  void step(const Real dt, std::map<Vec2f, Vec2f> forces)
  {
    // std::cout << "body pos: " << body->X.x << " " << body->X.y << " " << body->X.z << "; t=" << _sim_t << " (dt=" << dt << ")" << std::endl;

    computeForceAndTorque(forces);

    body->P += dt * body->F;
    body->V = body->P / body->M;;
    body->X += dt * body->V;
 
    body->Iinv = body->R * body->I0inv * body->R.transposed();
    body->L += dt * body->tau;
    body->omega = body->Iinv * body->L;

    Quatf my_w_hat(body->omega.x, body->omega.y, body->omega.z, 0);
    body->q += 0.5f * dt * (my_w_hat * body->q);
    body->q.normalize();
    
    body->R = body->q.to_rotation_mat();
    
    update_pos2d();

    ++_step;
    _sim_t += dt;
  }

  BodyAttributes *body;

private:
  void computeForceAndTorque(std::map<Vec2f, Vec2f> forces)
  {
    body->F = _g * body->M;
    body->tau = Vec3f(0.0f);

    for (auto const& force : forces) {
      Vec3f F = Vec3f(force.second.x, force.second.y, 0.0);
      Vec3f p = Vec3f(force.first.x, force.first.y, 0.0);
      body->F += F;
      body->tau += (p - body->X).crossProduct(F);
    }

    // std::cout << "body force : " << body->F << "; and torque: " << body->tau << std::endl;
  }

  void update_pos2d() {
    for (int i = 0; i < body->vertexCount(); ++i) {
      body->vpos2d[i] = (body->R * body->vpos[i] + body->X).get2d();
    }
  }

  // simulation parameters
  Vec3f _g;                     // gravity
  tIndex _step;                 // step count
  Real _sim_t;                 // simulation time
};

#endif  /* _RIGIDSOLVER_HPP_ */
