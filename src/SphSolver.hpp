#ifndef _SPHSOLVER_HPP_
#define _SPHSOLVER_HPP_

#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <set>

#ifndef M_PI
    #define M_PI 3.141592
#endif

#include "Vector.hpp"
#include "RigidSolver.hpp"

// SPH Kernel function: cubic spline
class CubicSpline
{
public:
  explicit CubicSpline(const Real h = 1) : _dim(2)
  {
    setSmoothingLen(h);
  }
  void setSmoothingLen(const Real h)
  {
    const Real h2 = square(h), h3 = h2 * h;
    _h = h;
    _sr = 2e0 * h;
    _c[0] = 2e0 / (3e0 * h);
    _c[1] = 10e0 / (7e0 * M_PI * h2);
    _c[2] = 1e0 / (M_PI * h3);
    _gc[0] = _c[0] / h;
    _gc[1] = _c[1] / h;
    _gc[2] = _c[2] / h;
  }
  Real smoothingLen() const { return _h; }
  Real supportRadius() const { return _sr; }

  Real f(const Real l) const
  {
    const Real q = l / _h;
    if (q < 1e0)
      return _c[_dim - 1] * (1e0 - 1.5 * square(q) + 0.75 * cube(q));
    else if (q < 2e0)
      return _c[_dim - 1] * (0.25 * cube(2e0 - q));
    return 0;
  }
  Real derivative_f(const Real l) const
  {
    const Real q = l / _h;
    if (q <= 1e0)
      return _gc[_dim - 1] * (-3e0 * q + 2.25 * square(q));
    else if (q < 2e0)
      return -_gc[_dim - 1] * 0.75 * square(2e0 - q);
    return 0;
  }

  Real w(const Vec2f &rij) const { return f(rij.length()); }
  Vec2f grad_w(const Vec2f &rij) const { return grad_w(rij, rij.length()); }
  Vec2f grad_w(const Vec2f &rij, const Real len) const
  {
    if (len < 1e-6) {
      return {0.0, 0.0};
    }

    return derivative_f(len) * rij / len;
  }

private:
  unsigned int _dim;
  Real _h, _sr, _c[3], _gc[3];
};

class SphSolver
{
public:
  explicit SphSolver(
      const Real nu = 0.1, const Real h = 0.5, const Real density = 1e3,
      const Vec2f g = Vec2f(0, -9.8), const Real eta = 0.01, const Real gamma = 7.0)
      : _kernel(h), _nu(nu), _h(h), _d0(density), _g(g), _eta(eta), _gamma(gamma)
  {
    _dt = 0.001;
    _m0 = _d0 * _h * _h;
    _c = std::fabs(_g.y) / _eta;
    _k = _d0 * _c * _c / _gamma;
  }

  // assume an arbitrary grid with the size of res_x*res_y; a fluid mass fill up
  // the size of f_width, f_height; each cell is sampled with 2x2 particles.
  void initScene(
      const int res_x, const int res_y, const int f_width, const int f_height, RigidSolver *rigid = nullptr)
  {
    _resX = res_x;
    _resY = res_y;

    // set wall for boundary
    _l = 0.5 * _h;
    _r = static_cast<Real>(res_x) - 0.5 * _h;
    _b = 0.5 * _h;
    _t = static_cast<Real>(res_y) - 0.5 * _h;

    _fixed_solid_particles = 0;

    // // bottom wall
    // int j = 0;
    // for (int i = 0; i < resX(); ++i)
    // {
    //   _pos_fluid.push_back(Vec2f(i + 0.25, j + 0.25));
    //   _pos_fluid.push_back(Vec2f(i + 0.75, j + 0.25));
    //   _pos_fluid.push_back(Vec2f(i + 0.25, j + 0.75));
    //   _pos_fluid.push_back(Vec2f(i + 0.75, j + 0.75));

    //   _fixed_solid_particles += 4;
    // }

    // // left wall
    // int i = 0;
    // for (int j = 1; j < 15; ++j)
    // {
    //   _pos_fluid.push_back(Vec2f(i + 0.25, j + 0.25));
    //   _pos_fluid.push_back(Vec2f(i + 0.75, j + 0.25));
    //   _pos_fluid.push_back(Vec2f(i + 0.25, j + 0.75));
    //   _pos_fluid.push_back(Vec2f(i + 0.75, j + 0.75));

    //   _fixed_solid_particles += 4;
    // }

    // // right wall
    // i = resX() - 1;
    // for (int j = 1; j < 15; ++j)
    // {
    //   _pos_fluid.push_back(Vec2f(i + 0.25, j + 0.25));
    //   _pos_fluid.push_back(Vec2f(i + 0.75, j + 0.25));
    //   _pos_fluid.push_back(Vec2f(i + 0.25, j + 0.75));
    //   _pos_fluid.push_back(Vec2f(i + 0.75, j + 0.75));

    //   _fixed_solid_particles += 4;
    // }

    // sample a fluid mass
    for (int j = 3; j < f_height + 3; ++j)
    {
      for (int i = 0; i < f_width; ++i)
      {
        _pos_fluid.push_back(Vec2f(i + 0.25, j + 0.25));
        _pos_fluid.push_back(Vec2f(i + 0.75, j + 0.25));
        _pos_fluid.push_back(Vec2f(i + 0.25, j + 0.75));
        _pos_fluid.push_back(Vec2f(i + 0.75, j + 0.75));
      }
    }

    for (int j = 0; j < 3; ++j)
    {
      for (int i = 0; i < resX(); ++i)
      {
        _pos_fluid.push_back(Vec2f(i + 0.25, j + 0.25));
        _pos_fluid.push_back(Vec2f(i + 0.75, j + 0.25));
        _pos_fluid.push_back(Vec2f(i + 0.25, j + 0.75));
        _pos_fluid.push_back(Vec2f(i + 0.75, j + 0.75));
      }
    }

    // make sure for the other particle quantities
    _vel = std::vector<Vec2f>(_pos_fluid.size(), Vec2f(0, 0));
    _acc = std::vector<Vec2f>(_pos_fluid.size(), Vec2f(0, 0));
    _p = std::vector<Real>(_pos_fluid.size(), 0);
    _d = std::vector<Real>(_pos_fluid.size(), 0);

    _col = std::vector<float>(_pos_fluid.size() * 4, 1.0); // RGBA
    _vln = std::vector<float>(_pos_fluid.size() * 4, 0.0); // GL_LINES

    // don't need to update _pidxInGrid each step if keep updated on the positions changes
    for (int i = 0; i < particleCountFluid(); ++i) {
      tIndex cell_x = static_cast<int>(_pos_fluid[i].x);
      tIndex cell_y = static_cast<int>(_pos_fluid[i].y);

      _pidxInGridFluid[idx1d(cell_x, cell_y)].insert(i);

      if (i < _fixed_solid_particles) {
        _col[i * 4 + 0] = 1.0;
        _col[i * 4 + 1] = 0.6;
        _col[i * 4 + 2] = 0.6;
      } else {
        updateColor(i);
      }
    }

    _body = rigid;
    constructBoundaryParticles();
  }

  std::map<Vec2f, Vec2f> update(bool show_velocities)
  {
    updateBoundaryParticles();

    _acc = std::vector<Vec2f>(particleCountFluid(), Vec2f(0, 0));
    _neighbors_fluid = std::vector<std::vector<tIndex>>(particleCountFluid(), std::vector<tIndex>(0, -1));
    _neighbors_body = std::vector<std::vector<tIndex>>(particleCountFluid(), std::vector<tIndex>(0, -1));

    for (int i = 0; i < particleCountFluid(); ++i) {
      _neighbors_fluid[i] = buildNeighbor(i, _pidxInGridFluid, _pos_fluid, _pos_fluid);
      _neighbors_body[i] = buildNeighbor(i, _pidxInGridBody, _pos_fluid, _pos_body);

      computeDensity(i);
    }

    computePressure();

    for (int i = 0; i < particleCountFluid(); ++i) {
      applyBodyForce(i);
      applyPressureForce(i);
      applyViscousForce(i);
    }

    for (int i = _fixed_solid_particles; i < particleCountFluid(); ++i) {
      updateVelocity(i);
      updatePosition(i);
      
      resolveCollision(i);
      updateColor(i);
      if (show_velocities)
        updateVelLine(i);
    }

    std::map<Vec2f, Vec2f> forces;
    for (int i = 0; i < particleCountBody(); ++i) {
      if (_body_forces[i] != Vec2f(0,0) && forces[_pos_body[i]] == Vec2f(0,0)) {
        forces[_pos_body[i]] += _body_forces[i];
      }
    }

    return forces;
  }

  tIndex particleCountFluid() const { return _pos_fluid.size(); }
  tIndex particleCountBody() const { return _pos_body_initial.size(); }

  const Vec2f &positionFluid(const tIndex i) const { return _pos_fluid[i]; }
  const Vec2f &positionBody(const tIndex i) const { return _pos_body[i]; }
  
  const float &colorFluid(const tIndex i) const { return _col[i]; }
  const float &colorBody(const tIndex i) const { return _col2[i]; }
  const float &vline(const tIndex i) const { return _vln[i]; }

  int resX() const { return _resX; }
  int resY() const { return _resY; }

  Real get_dt() const { return _dt; }

private:
  void constructBoundaryParticles() {
    _pos_body_initial.clear();
    auto body_pos = _body->body->vpos2d;

    for (int i = 0; i < body_pos.size(); ++i) {
      Vec2f v1 = body_pos[i];
      Vec2f v2 = body_pos[0];
      if (i < body_pos.size() - 1) {
        v2 = body_pos[i + 1];
      }

      float l = (v2 - v1).length();
      float dl = 0.25;
      int n = l / dl;
      dl = l / n;

      for (int i = 0; i < n; ++i) {
        float t = dl * i / l;
        _pos_body_initial.push_back((1 - t) * v1 + t * v2);
      }
    }

    _pos_body = std::vector<Vec2f>(particleCountBody(), Vec2f(0.0));
    updateBoundaryParticles();
    computeBoundaryContribution();
  }

  void updateBoundaryParticles() {
    _body_forces = std::vector<Vec2f>(particleCountBody(), Vec2f(0.0));
    _pidxInGridBody.clear();

    // fill in _pidxInGridBody
    for (int i = 0; i < particleCountBody(); ++i) {

      _pos_body[i] = _body->translate(_pos_body_initial[i]);

      tIndex cell_x = static_cast<int>(_pos_body[i].x);
      tIndex cell_y = static_cast<int>(_pos_body[i].y);

      _pidxInGridBody[idx1d(cell_x, cell_y)].insert(i);
    }
  }

  void computeBoundaryContribution() {
    _boundary_contrib = std::vector<float>(particleCountBody(), 0.0);
    
    for (int i = 0; i < particleCountBody(); ++i) {
      std::vector<tIndex> neighs = buildNeighbor(i, _pidxInGridBody, _pos_body, _pos_body);
      Real sum_w = 0.0;
      for (int j : neighs) {
        sum_w += _kernel.w(_pos_body[i] - _pos_body[j]);
      }

      _boundary_contrib[i] = _d0 / sum_w;
    }
  }

  // search of neighbors, can be used both for body_parcticles and fluid particles 
  std::vector<tIndex> buildNeighbor(tIndex i,
                     std::map<tIndex, std::set<tIndex>> &idxInGrid,
                     std::vector<Vec2f> &pos_list_i,
                     std::vector<Vec2f> &pos_list_j)
  {

    // for (int i = 0; i < particleCountFluid(); ++i) {
      
      int r = ceil(_kernel.supportRadius());
      tIndex cell_x = static_cast<int>(pos_list_i[i].x);
      tIndex cell_y = static_cast<int>(pos_list_i[i].y);
      std::vector<tIndex> neighs;

      for (int dx = -r; dx < r + 1; ++dx) {
        for (int dy = -r; dy < r + 1; ++dy) {

          tIndex cell_neigh = idx1d(cell_x + dx, cell_y + dy);

          // check that cell is not outside of the grid
          if (cell_neigh == -1) {
            continue;
          }

          for (tIndex i_neigh : idxInGrid[cell_neigh]) {
            if ((pos_list_i[i] - pos_list_j[i_neigh]).length() <= _kernel.supportRadius()) {
              neighs.push_back(i_neigh);
            }
          }
        }
      }

      return neighs;
    // }
  }

  void computeDensity(tIndex i)
  {
    // for (int i = 0; i < particleCountFluid(); ++i) {
      _d[i] = 0.;
      for (int j : _neighbors_fluid[i]) {
        _d[i] += _m0 * _kernel.w(_pos_fluid[i] - _pos_fluid[j]);
      }
      for (int j : _neighbors_body[i]) {
        _d[i] += _boundary_contrib[j] * _kernel.w(_pos_fluid[i] - _pos_body[j]);
      }
      
    // }
  }

  void computePressure()
  {
    for (int i = 0; i < particleCountFluid(); ++i) {
      _p[i] = fmax(_k * (pow((_d[i] / _d0), 7.0) - 1.0), 0.0);
    }
  }

  void applyBodyForce(tIndex i)
  {
    // for (int i = 0; i < particleCountFluid(); ++i) {
      _acc[i] += _g;
    // }
  }

  void applyPressureForce(tIndex i)
  {
    // da1 ~ 2 * da2
    Vec2f da1(0.0);  // video 15
    Vec2f da2(0.0);  // video 14

    // for (int i = 0; i < particleCountFluid(); ++i) {
      Real dpi = _p[i] / (_d[i] * _d[i]);
      for (int j : _neighbors_fluid[i]) {
        da1 += _m0 * (dpi + _p[j] / (_d[j] * _d[j]))
                   * _kernel.grad_w(_pos_fluid[i] - _pos_fluid[j]);

        da2 += _m0 * ((_p[i] + _p[j]) / (_d[i] * _d[j])) / 2.
                   * _kernel.grad_w(_pos_fluid[i] - _pos_fluid[j]);
      }

      _acc[i] -= da1;

      // on the solid's boundary
      for (int j : _neighbors_body[i]) {
        Vec2f da = _boundary_contrib[j] * _p[i] / (_d[i] * _d[i])
                 * _kernel.grad_w(_pos_fluid[i] - _pos_body[j]);

        _acc[i] -= da;

        Vec2f F = _m0 * da;
        _body_forces[j] += F;
      }
    // }
  }

  void applyViscousForce(tIndex i)
  {
    // for (int i = 0; i < particleCountFluid(); ++i) {
      for (int j : _neighbors_fluid[i]) {
        Vec2f xij = _pos_fluid[i] - _pos_fluid[j];
        _acc[i] += 2.0 * _nu * _m0 / (_d[i] + _d[j])
                   * xij.dotProduct(_kernel.grad_w(xij))
                   / (xij.dotProduct(xij) + 0.01 * _h * _h)
                   * (_vel[i] - _vel[j]);
      }

      // on the solid's boundary (needs velocities)
      for (int j : _neighbors_body[i]) {
        Vec2f xij = _pos_fluid[i] - _pos_body[j];
        Vec2f da = 2.0 * _nu / _d[i] * _boundary_contrib[j]
                   * xij.dotProduct(_kernel.grad_w(xij))
                   / (xij.dotProduct(xij) + 0.01 * _h * _h)
                   * (_vel[i] - _body->point_velocity(_pos_body[j]));

        _acc[i] += da;

        Vec2f F = _m0 * da;
        _body_forces[j] -= F;
      }
    // }
  }

  void updateVelocity(tIndex i)
  {
    // for (int i = 0; i < particleCountFluid(); ++i) {
      _vel[i] += _dt * _acc[i];
    // }
  }

  void updatePosition(tIndex i)
  {
    // for (int i = 0; i < particleCountFluid(); ++i) {
      // _pos_fluid[i] += _dt * _vel[i];
      updatePositionValue(i, _pos_fluid[i] + _dt * _vel[i]);
    // }
  }

  void updatePositionValue(tIndex i, Vec2f value) {
    if (i < _fixed_solid_particles) {
      return;
    }

    tIndex cell_x = static_cast<int>(_pos_fluid[i].x);
    tIndex cell_y = static_cast<int>(_pos_fluid[i].y);

    _pos_fluid[i] = value;

    tIndex cell_x_new = static_cast<int>(_pos_fluid[i].x);
    tIndex cell_y_new = static_cast<int>(_pos_fluid[i].y);

    // if the particle has moved to another cell of the grid
    if (cell_x != cell_x_new || cell_y != cell_y_new) {
      _pidxInGridFluid[idx1d(cell_x, cell_y)].erase(i);
      _pidxInGridFluid[idx1d(cell_x_new, cell_y_new)].insert(i);
    }
  }

  void updatePositionX(tIndex i, float new_x) {
      updatePositionValue(i, {new_x, _pos_fluid[i].y});
  }

  void updatePositionY(tIndex i, float new_y) {
      updatePositionValue(i, {_pos_fluid[i].x, new_y});
  }

  // simple collision detection/resolution for each particle
  void resolveCollision(tIndex i)
  {
    // for (tIndex i = 0; i < particleCountFluid(); ++i)
    // {
      
      if (_pos_fluid[i].x < _l || _pos_fluid[i].x > _r) {
        // _pos_fluid[i].x = clamp(_pos_fluid[i].x, _l, _r);
        updatePositionX(i, clamp(_pos_fluid[i].x, _l, _r));
        _vel[i].x *= -0.8;
        _vel[i].y *= 0.8;
      }

      if (_pos_fluid[i].y < _b || _pos_fluid[i].y > _t) {
        // _pos_fluid[i].y = clamp(_pos_fluid[i].y, _b, _t);
        updatePositionY(i, clamp(_pos_fluid[i].y, _b, _t));
        _vel[i].y *= -0.8;
        _vel[i].x *= 0.8;
      }
    // }
  }

  void updateColor(tIndex i)
  {
    if (i < _fixed_solid_particles) {
      return;
    }

    // for (tIndex i = 0; i < particleCountFluid(); ++i)
    // {
      _col[i * 4 + 0] = 0.6;
      _col[i * 4 + 1] = 0.6;
      _col[i * 4 + 2] = _d[i] / _d0;
    // }
  }

  void updateVelLine(tIndex i)
  {
    // for (tIndex i = 0; i < particleCountFluid(); ++i)
    // {
      _vln[i * 4 + 0] = _pos_fluid[i].x;
      _vln[i * 4 + 1] = _pos_fluid[i].y;
      _vln[i * 4 + 2] = _pos_fluid[i].x + _vel[i].x;
      _vln[i * 4 + 3] = _pos_fluid[i].y + _vel[i].y;
    // }
  }

  inline tIndex idx1d(const int i, const int j) {
    if (i < 0 || j < 0) {
      return -1;
    }
    if (i >= resX() || j >= resY()) {
      return -1;
    }

    return i + j * resX();
  }

  const CubicSpline _kernel;
  RigidSolver *_body;

  // particle data
  std::vector<Vec2f> _pos_fluid;    // position of fluid particles
  std::vector<Vec2f> _pos_body;     // position of body boundary particles in world space
  std::vector<Vec2f> _pos_body_initial;     // position of body boundary particles in bod space

  std::vector<Vec2f> _vel; // velocity
  std::vector<Vec2f> _acc; // acceleration
  std::vector<Real> _p;    // pressure
  std::vector<Real> _d;    // density

  std::vector<Vec2f> _body_forces;        // forces acting on the body's particles

  std::vector<Real> _boundary_contrib;

  std::map<tIndex, std::set<tIndex>> _pidxInGridFluid; // to find fluid neighbor particles faster
  std::map<tIndex, std::set<tIndex>> _pidxInGridBody; // to find boundary neighbor particles faster
  std::vector<std::vector<tIndex>> _neighbors_fluid;
  std::vector<std::vector<tIndex>> _neighbors_body;

  std::vector<float> _col; // fluid particle color; just for visualization
  std::vector<float> _col2; // fluid particle color; just for visualization
  std::vector<float> _vln; // fluid particle velocity lines; just for visualization

  // simulation
  Real _dt; // time step
  tIndex _fixed_solid_particles; // number of fixed particles stored in the beginning of _pos_fluid

  int _resX, _resY; // background grid resolution

  // wall
  Real _l, _r, _b, _t; // wall (boundary)

  // SPH coefficients
  Real _nu; // viscosity coefficient
  Real _d0; // rest density
  Real _h;  // particle spacing (i.e., diameter)
  Vec2f _g; // gravity

  Real _m0; // rest mass
  Real _k;  // EOS coefficient

  Real _eta;
  Real _c;     // speed of sound
  Real _gamma; // EOS power factor
};


#endif  /* _SPHSOLVER_HPP_ */