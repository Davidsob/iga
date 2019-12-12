#pragma once

#include "base/SimulationClock.h"

#include "splines/GeometricObject.h"

#include "iga/GeometricDofManager.h"
#include "iga/Quadrature.h"
#include "iga/QuadratureMesh.h"

#include "iga/constraints/ConstraintManager.h"
#include "iga/constraints/ConstraintBase.h"

#include "iga/weakforms/WeakFormManager.h"
#include "iga/weakforms/WeakFormTags.h"

#include "utils/MatrixTypes.h"

#include <cassert>
#include <memory>
#include <algorithm>


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

template<typename PrimaryVariable>
class LaSystem : public LaSystemBase
{
public:
  using BilinearEngine_t = typename Operator<IntegrationPoint,DynamicMatrixR>::OperatorEngine;
  using LinearEngine_t   = typename Operator<IntegrationPoint,DynamicVectorR>::OperatorEngine;

  LaSystem(GeometricObject const *obj)
    : LaSystemBase()
    , _obj(obj)
  {
    initialize();
  }

  ~LaSystem() = default;

  void sparse_lhs(std::vector<Triplet> &data) const override
  {
    auto &mgr = GeometricDofManager::instance();
    auto const dof = mgr.dofForShape(_obj,PrimaryVariable::ndof);
    auto const size = dof.size();
    DynamicMatrixR x(DynamicMatrixR::Zero(size,size));

    for (auto const ip : QuadratureMesh(_obj))
    {
      for (auto &form : _bilinear_forms)
      {
        quadrature::gauss(ip,x,*form.get());
      }
    }
    convert::to<void>(x,dof,dof,data);
  }

  void sparse_stiffness(std::vector<Triplet> &data) const
  {
    auto &mgr = GeometricDofManager::instance();
    auto const dof = mgr.dofForShape(_obj,PrimaryVariable::ndof);
    auto const size = dof.size();
    DynamicMatrixR x(DynamicMatrixR::Zero(size,size));
    
    for (auto const ip : QuadratureMesh(_obj))
    {
      for (auto &form : _stiffness_forms)
      {
        quadrature::gauss(ip,x,*form.get());
      }
    }
    convert::to<void>(x,dof,dof,data);
  }

  void sparse_mass(std::vector<Triplet> &data) const
  {
    auto &mgr = GeometricDofManager::instance();
    auto const dof = mgr.dofForShape(_obj,PrimaryVariable::ndof);
    auto const size = dof.size();
    DynamicMatrixR x(DynamicMatrixR::Zero(size,size));
    
    for (auto const ip : QuadratureMesh(_obj))
    {
      for (auto &form : _mass_forms)
      {
        quadrature::gauss(ip,x,*form.get());
      }
    }
    convert::to<void>(x,dof,dof,data);
  }

  void sparse_rhs(std::vector<Triplet> &data) const override
  {
    auto &mgr = GeometricDofManager::instance();
    auto const dof = mgr.dofForShape(_obj,PrimaryVariable::ndof);
    auto const size = dof.size();
    DynamicVectorR x(DynamicVectorR::Zero(size));

    for (auto const &ip : QuadratureMesh(_obj))
    {
      for (auto &form : _linear_forms)
      {
        quadrature::gauss(ip,x,*form.get());
      }
    }
    convert::to<void>(x,dof,data); // converts to triplets
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

    {
      auto &mgr = GeometricDofManager::instance(); 
      mgr.addShape(_obj);
    }
  }

  GeometricObject const * _obj;
  std::vector<std::unique_ptr<BilinearEngine_t>> _bilinear_forms;
  std::vector<std::unique_ptr<BilinearEngine_t>> _stiffness_forms;
  std::vector<std::unique_ptr<BilinearEngine_t>> _mass_forms;
  std::vector<std::unique_ptr<LinearEngine_t>> _linear_forms;
};

template<typename PrimaryVariable>
class ConstrainedLaSystem : public LaSystemBase
{
public:
  using BilinearEngine_t = typename Operator<IntegrationPoint, DynamicMatrixR>::OperatorEngine;
  using LinearEngine_t   = typename Operator<IntegrationPoint, DynamicVectorR>::OperatorEngine;

  template<typename Eng>
  using EnginePair_t = std::pair<GeometricObject const *, std::unique_ptr<Eng>>;

  ConstrainedLaSystem()
    : LaSystemBase()
  {
    initialize();
  }

  ~ConstrainedLaSystem() {} 

  virtual void lhs(SparseMatrixR &A) const
  {
    assert(0);
    // auto offset = 0;
    // static auto renumber = [&](Triplet &t) { t.setRow(offset++); };

    // std::vector<Triplet> lhs;
    // sparse_lhs(lhs);
    // for_each(lhs.begin(), lhs.end(), renumber);

    // auto cdof = lhs.size();
    // auto ndof = _dofManager.ids.size();
    // A.resize(cdof, ndof);
    // A.setFromTriplets(lhs.begin(), lhs.end());
  }

  virtual void rhs(DynamicVectorR &b) const
  {
    auto offset = 0;
    static auto renumber = [&](Triplet &t) { t.setRow(offset++); };

    static auto set_vector =
      [&](Triplet const &t) { b[t.row()] += t.value(); };

    std::vector<Triplet> rhs;
    sparse_rhs(rhs);
    for_each(rhs.begin(), rhs.end(), renumber);

    auto cdof = rhs.size();
    b.resize(cdof);
    std::fill(b.data(), b.data()+cdof, 0);
    for_each(rhs.begin(), rhs.end(), set_vector);
  }

  void sparse_lhs(std::vector<Triplet> &data) const override
  {
    auto &mgr = GeometricDofManager::instance();

    for (auto const &pair : _bilinear_forms) // loop of shape engine pairs
    {
      auto const shape = pair.first;
      auto const form = pair.second.get();
      auto const dof = mgr.dofForShape(shape,PrimaryVariable::ndof);
      auto const size = dof.size();
      DynamicMatrixR x(DynamicMatrixR::Zero(size,size));
      for (auto const &ip : QuadratureMesh(shape))
      {
        quadrature::gauss(ip,x,*form);
      }
      convert::to<void>(x,dof,dof,data);
    }

  }

  void sparse_rhs(std::vector<Triplet> &data) const override
  {
    auto &mgr = GeometricDofManager::instance();

    for (auto const &pair : _linear_forms) // loop of shape engine pairs
    {
      auto const shape = pair.first;
      auto const form = pair.second.get();
      auto const dof = mgr.dofForShape(shape,PrimaryVariable::ndof);
      auto const size = dof.size();
      DynamicVectorR x(DynamicVectorR::Zero(size));
      for (auto const &ip : QuadratureMesh(shape))
      {
        quadrature::gauss(ip,x,*form);
      }
      convert::to<void>(x,dof,data);
    }
  }

protected:
  void initialize()
  {
    auto &mgr = GeometricDofManager::instance();
    for (auto const &pair : ConstraintManager::instance())
    {
      _linear_forms  .push_back({pair.first, std::unique_ptr<LinearEngine_t  >(pair.second->getResidual()->createLinearEngine())});
      _bilinear_forms.push_back({pair.first, std::unique_ptr<BilinearEngine_t>(pair.second->getJacobian()->createBilinearEngine())});
      mgr.addShape(pair.first); // handled by the element mappers!
    }
  }

  std::vector< EnginePair_t<BilinearEngine_t> > _bilinear_forms;
  std::vector< EnginePair_t<LinearEngine_t> >   _linear_forms;
};


template<typename PrimaryVariable>
class FeLaSystem
{
public:
  using Unconstrained_t = LaSystem<PrimaryVariable>;
  using Constrained_t = ConstrainedLaSystem<PrimaryVariable>;

  FeLaSystem(GeometricObject const *obj)
    : _obj(obj), _lasystem(obj), _constrained()
  {
  }

  virtual ~FeLaSystem() {} 

  void discretize(SparseMatrixR &A, DynamicVectorR &b) const
  {
    lhs(A);
    rhs(b);
  }

  virtual void lhs(SparseMatrixR &A) const
  {
    size_t ndof = 0;
    std::vector<Triplet> lhs;
    sparse_lhs(lhs, ndof);
    auto cmp = [](auto const &a, auto const &b) { return a.row() < b.row(); };
    std::sort(lhs.begin(),lhs.end(),cmp);
    A.resize(ndof, ndof);
    A.setFromTriplets(lhs.begin(), lhs.end());
  }

  void composed_lhs(SparseMatrixR &A, double a, double b) const
  {
    static auto scale_a = [&](Triplet &t) { t.setValue(a*t.value()); };
    static auto scale_b = [&](Triplet &t) { t.setValue(b*t.value()); };

    std::vector<Triplet> mass, stiffness, cdata;
    _lasystem.sparse_mass(mass);
    _lasystem.sparse_stiffness(stiffness);
    for_each(mass.begin(), mass.end(), scale_a);
    for_each(stiffness.begin(), stiffness.end(), scale_b);

    size_t ndof = 0;
    sparse_clhs(cdata, ndof);
    ndof += unconstrainedDof();

    mass.insert(mass.end(), stiffness.begin(), stiffness.end());
    mass.insert(mass.end(), cdata.begin(), cdata.end());

    A.resize(ndof, ndof);
    A.setFromTriplets(mass.begin(), mass.end());
  }

  virtual void rhs(DynamicVectorR &b) const
  {
    static auto set_vector =
      [&](Triplet const &t) { b[t.row()] += t.value(); };
    size_t ndof = 0;
    std::vector<Triplet> rhs;
    sparse_rhs(rhs, ndof);
    b.resize(ndof);
    std::fill(b.data(), b.data()+ndof, 0);
    for_each(rhs.begin(), rhs.end(), set_vector);
  }

  virtual void stiffness(SparseMatrixR &A) const
  {
    size_t ndof = 0;
    std::vector<Triplet> lhs;
    sparse_stiffness(lhs, ndof);

    A.resize(ndof, ndof);
    A.setFromTriplets(lhs.begin(), lhs.end());
  }

  virtual void mass(SparseMatrixR &A) const
  {
    size_t ndof = 0;
    std::vector<Triplet> lhs;
    sparse_mass(lhs, ndof);

    A.resize(ndof, ndof);
    A.setFromTriplets(lhs.begin(), lhs.end());
  }

  virtual void unconstrained_mass(SparseMatrixR &A) const
  {
    std::vector<Triplet> lhs;
    _lasystem.sparse_mass(lhs);

    auto ndof = unconstrainedDof();
    A.resize(ndof, ndof);
    A.setFromTriplets(lhs.begin(), lhs.end());
  }

  virtual void unconstrained_stiffness(SparseMatrixR &A) const
  {
    std::vector<Triplet> lhs;
    _lasystem.sparse_stiffness(lhs);

    auto ndof = unconstrainedDof();
    A.resize(ndof, ndof);
    A.setFromTriplets(lhs.begin(), lhs.end());
  }

  void sparse_lhs(std::vector<Triplet> &data, size_t &ndof) const
  {
    std::vector<Triplet> cdata;

    _lasystem.sparse_lhs(data);

    size_t cdof;
    sparse_clhs(cdata, cdof);
    ndof += unconstrainedDof();
    ndof += cdof;
    data.insert(data.end(), cdata.begin(), cdata.end());
  }

  void sparse_stiffness(std::vector<Triplet> &data, size_t &ndof) const
  {
    std::vector<Triplet> cdata;

    _lasystem.sparse_stiffness(data);

    sparse_clhs(cdata, ndof);
    ndof += unconstrainedDof();

    data.insert(data.end(), cdata.begin(), cdata.end());
  }

  void sparse_mass(std::vector<Triplet> &data, size_t &ndof) const
  {
    std::vector<Triplet> cdata;

    _lasystem.sparse_mass(data);

    sparse_clhs(cdata, ndof);
    ndof += unconstrainedDof();

    data.insert(data.end(), cdata.begin(), cdata.end());
  }

  void sparse_rhs(std::vector<Triplet> &data, size_t &ndof) const
  {
    std::vector<Triplet> cdata;
    _lasystem.sparse_rhs(data);

    sparse_crhs(cdata, ndof);
    ndof += unconstrainedDof();

    data.insert(data.end(), cdata.begin(), cdata.end());
  }

  Unconstrained_t const &unconstrained_system() const { return _lasystem; }
  Constrained_t   const &constrained_system()   const { return _constrained; }

private:

  size_t unconstrainedDof() const
  { return GeometricDofManager::instance().ids.size()*PrimaryVariable::ndof; }

  size_t unique_dof(std::vector<Triplet> const &data) const
  {
    static auto same = [](auto const &a, auto const &b) { return a.row() == b.row(); };
    // assumes sorted data
    auto cpy = data;
    auto it = std::unique(cpy.begin(), cpy.end(), same);
    return std::distance(cpy.begin(),it);
  }

  void sparse_crhs(std::vector<Triplet> &cdata, size_t &cdof) const
  {
    static auto cmp = [](auto const &a, auto const &b)
    {
      return a.row() < b.row();
    };

    _constrained.sparse_rhs(cdata);
    std::sort(cdata.begin(),cdata.end(),cmp);
    cdof = unique_dof(cdata);
    // next renumber all starting with zero
    auto const offset = unconstrainedDof();
    int index = offset-1;
    int last  = -1;

    static auto renumber = [&index,&last](auto &t)
    {
      if ( t.row()!= last) 
      {
        last = t.row();
        index++;
      }
      t.setRow(index);
    };
    // // static auto renumber = [&](Triplet &t) { t.setRow(offset++); };
    // static auto renumber = [&](Triplet &t) { t.setRow(t.row()+offset); };

    for_each(cdata.begin(), cdata.end(), renumber);
  }

  void sparse_clhs(std::vector<Triplet> &cdata, size_t &cdof) const
  {
    static auto cmp = [](auto const &a, auto const &b)
    {
      return a.row() < b.row();
    };

    _constrained.sparse_lhs(cdata);
    std::sort(cdata.begin(),cdata.end(),cmp);
    cdof = unique_dof(cdata);

    // next renumber all starting with zero
    auto const offset = unconstrainedDof();
    int index = offset-1;
    int last  = -1;

    static auto renumber = [&index,&last](auto &t)
    {
      if ( t.row()!= last) 
      {
        last = t.row();
        index++;
      }
      t.setRow(index);
    };

    static auto transpose = [&](Triplet &t) { t.transpose(); };

    for_each(cdata.begin(), cdata.end(), renumber);

    auto cdataT = cdata; 
    for_each(cdataT.begin(), cdataT.end(), transpose);

    cdata.insert(cdata.end(), cdataT.begin(), cdataT.end());
  }

  GeometricObject const * _obj;
  Unconstrained_t  _lasystem;
  Constrained_t    _constrained;
};