#pragma once

// #include "GaussianIntegrationPoints.h"
// #include "IntegrationPointSet.h"
// #include "NumericalIntegration.h"
#include "base/SimulationClock.h"

#include "weakforms/WeakFormManager.h"
#include "weakforms/WeakFormTags.h"

#include "utils/MatrixTypes.h"

#include <cassert>
#include <memory>

class LaSystemBase
{
public:
  virtual ~LaSystemBase() = default;

  void discretize(SparseMatrixR &A, DynamicVectorR &b) const
  {
    lhs(A);
    rhs(b);
  }

  virtual void lhs(SparseMatrixR &A) const
  {
    std::vector<Triplet> lhs;
    sparse_lhs(lhs);
    A.setFromTriplets(lhs.begin(), lhs.end());
  }

  virtual void rhs(DynamicVectorR &b) const
  {
    static auto set_vector =
      [&](Triplet const &t) { b[t.row()] += t.value(); };

    std::vector<Triplet> rhs;
    sparse_rhs(rhs);
    for_each(rhs.begin(), rhs.end(), set_vector);
  }

protected:
  LaSystemBase() = default;
  virtual void sparse_lhs(std::vector<Triplet> &) const = 0;
  virtual void sparse_rhs(std::vector<Triplet> &) const = 0;
};


class LaSystem : public LaSystemBase
{
public:

  using BilinearEngine_t = typename Operator<DynamicMatrixR>::OperatorEngine;
  using LinearEngine_t = typename Operator<DynamicVectorR>::OperatorEngine;

  LaSystem(Mesh const &m, BoundaryMesh const &bm)
    : LaSystemBase(), _mesh(m)
  {
    initialize();
  }

  ~LaSystem() = default;

  void sparse_lhs(std::vector<Triplet> &data) const override
  {
    using point_t = GaussianIntegrationPoints<LocalPoint2d , 2>;
    using element_t = typename Mesh::it_t;
    using gauss_pts = IntegrationPointSet<point_t, element_t>;

    size_t dof = 0;

    for (auto &element : _mesh)
    {
      dof = element.dof.size();
      DynamicMatrixR x(DynamicMatrixR::Zero(dof, dof));
      gauss_pts gpts(element);
      for (auto &form : _bilinear_forms)
      {
        GaussianIntegration::integrate(gpts, *form.get(), x);
      }
      SparseUtils::matrixToTriplets(x, element, data);
    }
  }

  void sparse_stiffness(std::vector<Triplet> &data) const
  {
    using point_t = GaussianIntegrationPoints<LocalPoint2d , 2>;
    using element_t = typename Mesh::it_t;
    using gauss_pts = IntegrationPointSet<point_t, element_t>;

    size_t dof = 0;

    for (auto &element : _mesh)
    {
      dof = element.dof.size();
      DynamicMatrixR x(DynamicMatrixR::Zero(dof, dof));
      gauss_pts gpts(element);
      for (auto &form : _stiffness_forms)
      {
        if (dynamic_cast<ContributesToInteratia*>(form.get())) continue;
        GaussianIntegration::integrate(gpts, *form.get(), x);
      }
      SparseUtils::matrixToTriplets(x, element, data);
    }
  }

  void sparse_mass(std::vector<Triplet> &data) const
  {
    using point_t = GaussianIntegrationPoints<LocalPoint2d , 2>;
    using element_t = typename Mesh::it_t;
    using gauss_pts = IntegrationPointSet<point_t, element_t>;

    size_t dof = 0;

    for (auto &element : _mesh)
    {
      dof = element.dof.size();
      DynamicMatrixR x(DynamicMatrixR::Zero(dof, dof));
      gauss_pts gpts(element);
      for (auto &form : _mass_forms)
      {
        GaussianIntegration::integrate(gpts, *form.get(), x);
      }
      SparseUtils::matrixToTriplets(x, element, data);
    }
  }

  void sparse_rhs(std::vector<Triplet> &data) const override
  {
    using point_t = GaussianIntegrationPoints<LocalPoint2d , 2>;
    using element_t = typename Mesh::it_t;
    using gauss_pts = IntegrationPointSet<point_t, element_t>;

    size_t dof = 0;

    for (auto &element : _mesh)
    {
      dof = element.dof.size();
      DynamicVectorR x(DynamicVectorR::Zero(dof));
      gauss_pts gpts(element);
      for (auto &form : _linear_forms)
      {
        GaussianIntegration::integrate(gpts, *form.get(), x);
      }
      SparseUtils::vectorToTriplets(x, element, data);
    }
  }

protected:
  void initialize()
  {
    auto &wfm = WeakFormManager::instance();
    wfm.initialize();

    for (auto &form : wfm.bilinearForms())
    {
      _bilinear_forms.push_back(std::unique_ptr<BilinearEngine_t>(form->createBilinearEngine()));
    }
    for (auto &form : wfm.stiffnessForms())
    {
      _stiffness_forms.push_back(std::unique_ptr<BilinearEngine_t>(form->createBilinearEngine()));
    }
    for (auto &form : wfm.massForms())
    {
      _mass_forms.push_back(std::unique_ptr<BilinearEngine_t>(form->createBilinearEngine()));
    }
    for (auto &form : wfm.linearForms())
    {
      _linear_forms.push_back(std::unique_ptr<LinearEngine_t>(form->createLinearEngine()));
    }
  }

  Mesh const &_mesh;
  BoundaryMesh const &_bmesh;
  std::vector<std::unique_ptr<BilinearEngine_t>> _bilinear_forms;
  std::vector<std::unique_ptr<BilinearEngine_t>> _stiffness_forms;
  std::vector<std::unique_ptr<BilinearEngine_t>> _mass_forms;
  std::vector<std::unique_ptr<LinearEngine_t>> _linear_forms;
};

// template<class Mesh,
//          class BoundaryMesh,
//          typename T>
// class ConstrainedLaSystem : public LaSystemBase
// {
// public:
//   using Pair_t = std::pair<int, T const *>;
//   using BilinearEngine_t = typename Operator<DynamicMatrixR>::OperatorEngine;
//   using LinearEngine_t = typename Operator<DynamicVectorR>::OperatorEngine;

//   template<typename Eng>
//   using EnginePair_t = std::pair<int, std::unique_ptr<Eng>>;

//   ConstrainedLaSystem(Mesh const &m, BoundaryMesh const &bm,
//                      std::vector<Pair_t> const &point,
//                      std::vector<Pair_t> const &boundary = {})
//     : LaSystemBase(), _mesh(m), _bmesh(bm), _point(point), _boundary(boundary)
//   {
//     initialize();
//   }

//   ~ConstrainedLaSystem() {} 

//   virtual void lhs(SparseMatrixR &A) const
//   {
//     auto offset = 0;
//     static auto renumber = [&](Triplet &t) { t.setRow(offset++); };

//     std::vector<Triplet> lhs;
//     sparse_lhs(lhs);
//     for_each(lhs.begin(), lhs.end(), renumber);

//     auto cdof = lhs.size();
//     auto ndof = _mesh.vertices();
//     A.resize(cdof, ndof);
//     A.setFromTriplets(lhs.begin(), lhs.end());
//   }

//   virtual void rhs(DynamicVectorR &b) const
//   {
//     auto offset = 0;
//     static auto renumber = [&](Triplet &t) { t.setRow(offset++); };

//     static auto set_vector =
//       [&](Triplet const &t) { b[t.row()] += t.value(); };

//     std::vector<Triplet> rhs;
//     sparse_rhs(rhs);
//     for_each(rhs.begin(), rhs.end(), renumber);

//     auto cdof = rhs.size();
//     b.resize(cdof);
//     std::fill(b.data(), b.data()+cdof, 0);
//     for_each(rhs.begin(), rhs.end(), set_vector);
//   }

//   void sparse_lhs(std::vector<Triplet> &data) const override
//   {
//     auto &vertexes = _bmesh.vertexes();
//     LocalPoint2d lp;
//     lp.time = SimulationClock::instance().time();
//     for (auto const &pair : _point_lhs)
//     {
//       auto element = vertexes[pair.first];
//       auto eng = pair.second.get();
//       //  
//       lp.x = element.x;
//       lp.dof = element.dof;
//       //
//       auto x = (*eng)(lp);
//       SparseUtils::matrixToTriplets(x, element, data);
//     }
//   }

//   void sparse_rhs(std::vector<Triplet> &data) const override
//   {
//     auto &vertexes = _bmesh.vertexes();
//     LocalPoint2d lp;
//     lp.time = SimulationClock::instance().time();
//     for (auto const &pair : _point_rhs)
//     {
//       auto element = vertexes[pair.first];
//       auto eng = pair.second.get();
//       //  
//       lp.x = element.x;
//       lp.dof = element.dof;
//       //
//       auto x = (*eng)(lp);
//       SparseUtils::vectorToTriplets(x, element, data);
//     }
//   }

// protected:
//   void initialize()
//   {
//     for (auto &pair : _point)
//     {
//       _point_lhs.push_back({pair.first, std::unique_ptr<BilinearEngine_t>(pair.second->getJacobian()->createBilinearEngine())});
//       _point_rhs.push_back({pair.first, std::unique_ptr<LinearEngine_t>(pair.second->getResidual()->createLinearEngine())});
//     }

//     for (auto &pair : _boundary)
//     {
//       _boundary_lhs.push_back({pair.first, std::unique_ptr<BilinearEngine_t>(pair.second->getJacobian()->createBilinearEngine())});
//       _boundary_rhs.push_back({pair.first, std::unique_ptr<LinearEngine_t>(pair.second->getResidual()->createLinearEngine())});
//     }
//   }

//   Mesh const &_mesh;
//   BoundaryMesh const &_bmesh;
//   std::vector<Pair_t> const  &_point;
//   std::vector<Pair_t> const  &_boundary;
//   std::vector< EnginePair_t<BilinearEngine_t> > _point_lhs;
//   std::vector< EnginePair_t<LinearEngine_t> >   _point_rhs;
//   std::vector< EnginePair_t<BilinearEngine_t> > _boundary_lhs;
//   std::vector< EnginePair_t<LinearEngine_t> >   _boundary_rhs;

// };


// template<class Mesh,
//          class BoundaryMesh,
//          typename T>
// class FeLaSystem
// {
// public:
//   using Pair_t = std::pair<int, T const *>;
//   using Unconstrained_t = LaSystem<Mesh, BoundaryMesh>;
//   using Constrained_t = ConstrainedLaSystem<Mesh, BoundaryMesh, T>;

//   FeLaSystem(Mesh const &m, BoundaryMesh const &bm,
//                      std::vector<Pair_t> const &point,
//                      std::vector<Pair_t> const &boundary = {})
//     : _mesh(m), _bmesh(bm), _lasystem(m,bm), _constrained(m, bm, point, boundary)
//   {
//     initialize();
//   }

//   virtual ~FeLaSystem() {} 

//   void discretize(SparseMatrixR &A, DynamicVectorR &b) const
//   {
//     lhs(A);
//     rhs(b);
//   }

//   virtual void lhs(SparseMatrixR &A) const
//   {
//     size_t ndof = 0;
//     std::vector<Triplet> lhs;
//     sparse_lhs(lhs, ndof);

//     A.resize(ndof, ndof);
//     A.setFromTriplets(lhs.begin(), lhs.end());
//   }

//   void composed_lhs(SparseMatrixR &A, double a, double b) const
//   {
//     static auto scale_a = [&](Triplet &t) { t.setValue(a*t.value()); };
//     static auto scale_b = [&](Triplet &t) { t.setValue(b*t.value()); };

//     std::vector<Triplet> mass, stiffness, cdata;
//     _lasystem.sparse_mass(mass);
//     _lasystem.sparse_stiffness(stiffness);
//     for_each(mass.begin(), mass.end(), scale_a);
//     for_each(stiffness.begin(), stiffness.end(), scale_b);

//     size_t ndof = 0;
//     sparse_clhs(cdata, ndof);
//     ndof += _mesh.vertices();

//     mass.insert(mass.end(), stiffness.begin(), stiffness.end());
//     mass.insert(mass.end(), cdata.begin(), cdata.end());

//     A.resize(ndof, ndof);
//     A.setFromTriplets(mass.begin(), mass.end());
//   }

//   virtual void rhs(DynamicVectorR &b) const
//   {
//     static auto set_vector =
//       [&](Triplet const &t) { b[t.row()] += t.value(); };

//     size_t ndof = 0;
//     std::vector<Triplet> rhs;
//     sparse_rhs(rhs, ndof);
//     b.resize(ndof);
//     std::fill(b.data(), b.data()+ndof, 0);
//     for_each(rhs.begin(), rhs.end(), set_vector);
//   }


//   virtual void stiffness(SparseMatrixR &A) const
//   {
//     size_t ndof = 0;
//     std::vector<Triplet> lhs;
//     sparse_stiffness(lhs, ndof);

//     A.resize(ndof, ndof);
//     A.setFromTriplets(lhs.begin(), lhs.end());
//   }

//   virtual void mass(SparseMatrixR &A) const
//   {
//     size_t ndof = 0;
//     std::vector<Triplet> lhs;
//     sparse_mass(lhs, ndof);

//     A.resize(ndof, ndof);
//     A.setFromTriplets(lhs.begin(), lhs.end());
//   }

//   virtual void unconstrained_mass(SparseMatrixR &A) const
//   {
//     std::vector<Triplet> lhs;
//     _lasystem.sparse_mass(lhs);

//     auto ndof = _mesh.vertices();
//     A.resize(ndof, ndof);
//     A.setFromTriplets(lhs.begin(), lhs.end());
//   }

//   virtual void unconstrained_stiffness(SparseMatrixR &A) const
//   {
//     std::vector<Triplet> lhs;
//     _lasystem.sparse_stiffness(lhs);

//     auto ndof = _mesh.vertices();
//     A.resize(ndof, ndof);
//     A.setFromTriplets(lhs.begin(), lhs.end());
//   }

//   void sparse_lhs(std::vector<Triplet> &data, size_t &ndof) const
//   {
//     std::vector<Triplet> cdata;

//     _lasystem.sparse_lhs(data);

//     sparse_clhs(cdata, ndof);
//     ndof += _mesh.vertices();

//     data.insert(data.end(), cdata.begin(), cdata.end());
//   }

//   void sparse_stiffness(std::vector<Triplet> &data, size_t &ndof) const
//   {
//     std::vector<Triplet> cdata;

//     _lasystem.sparse_stiffness(data);

//     sparse_clhs(cdata, ndof);
//     ndof += _mesh.vertices();

//     data.insert(data.end(), cdata.begin(), cdata.end());
//   }

//   void sparse_mass(std::vector<Triplet> &data, size_t &ndof) const
//   {
//     std::vector<Triplet> cdata;

//     _lasystem.sparse_mass(data);

//     sparse_clhs(cdata, ndof);
//     ndof += _mesh.vertices();

//     data.insert(data.end(), cdata.begin(), cdata.end());
//   }

//   void sparse_rhs(std::vector<Triplet> &data, size_t &ndof) const
//   {
//     std::vector<Triplet> cdata;
//     _lasystem.sparse_rhs(data);

//     sparse_crhs(cdata, ndof);
//     ndof += _mesh.vertices();

//     data.insert(data.end(), cdata.begin(), cdata.end());
//   }

//   Unconstrained_t const &unconstrained_system() const { return _lasystem; }
//   Constrained_t const &constrained_system() const { return _constrained; }

// private:

//   void sparse_crhs(std::vector<Triplet> &cdata, size_t &cdof) const
//   {
//     _constrained.sparse_rhs(cdata);
//     cdof = cdata.size();

//     auto offset = _mesh.vertices();
//     static auto renumber = [&](Triplet &t) { t.setRow(offset++); };

//     for_each(cdata.begin(), cdata.end(), renumber);
//   }

//   void sparse_clhs(std::vector<Triplet> &cdata, size_t &cdof) const
//   {
//     _constrained.sparse_lhs(cdata);
//     cdof = cdata.size();

//     auto offset = _mesh.vertices();
//     static auto renumber = [&](Triplet &t) { t.setRow(offset++); };
//     static auto transpose = [&](Triplet &t) { t.transpose(); };

//     for_each(cdata.begin(), cdata.end(), renumber);

//     auto cdataT = cdata; 
//     for_each(cdataT.begin(), cdataT.end(), transpose);

//     cdata.insert(cdata.end(), cdataT.begin(), cdataT.end());
//   }

//   void initialize()
//   {
//   }

//   Mesh const &_mesh;
//   BoundaryMesh const &_bmesh;
//   Unconstrained_t const _lasystem;
//   Constrained_t const _constrained;
// };