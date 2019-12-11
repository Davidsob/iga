#ifndef WeakForms_h
#define WeakForms_h


#include "base/NamedObject.h"
#include "base/Operator.h"

#include "iga/IntegrationPoints.h"

#include "utils/MatrixTypeTraits.h"

#include <type_traits>

class LinearFormTag {};
class BilinearFormTag {};

class WeakFormBase
  : public NamedObject
{
public:

  ~WeakFormBase() = default; 

  using index_t           = IntegrationPoint;
  using linear_engine_t   = typename Operator<IntegrationPoint,DynamicVectorR>::OperatorEngine;
  using bilinear_engine_t = typename Operator<IntegrationPoint,DynamicMatrixR>::OperatorEngine;

  virtual linear_engine_t   * createLinearEngine() const = 0;
  virtual bilinear_engine_t * createBilinearEngine() const = 0;

protected:
  explicit WeakFormBase(std::string const &name = "WeakFormBase")
    : NamedObject(name)
  {} 
};


/**
 * @brief      Class for weak form.
 *
 * @tparam     Return_t  { description }
 */
template<typename Wf, typename Return_t, class Enable = void>
class WeakForm : public WeakFormBase, public BilinearFormTag
{
public:

  ~WeakForm() = default;

  template<typename ...Args>
  Return_t operator()(Args ...args) const
  {
    return static_cast<Wf const *>(this)->eval(args...);
  }

  class FormEngine : public Operator<IntegrationPoint, DynamicMatrixR>::OperatorEngine
  {
  public:
    using base_t = typename Operator<IntegrationPoint, DynamicMatrixR>::OperatorEngine;
    using var_t = WeakForm<Wf, Return_t>;
    using typename base_t::value_t;
    using typename base_t::index_t;

    FormEngine(var_t const &var)
      : _var(var)
    {}

    ~FormEngine() = default;

    value_t operator()(index_t const &idx) const override
    {
      return _var(idx);
    }

  private:
    var_t const &_var;
  };

  using typename WeakFormBase::linear_engine_t;
  using typename WeakFormBase::bilinear_engine_t;

  linear_engine_t * createLinearEngine() const override
  { 
    return nullptr;
  }

  bilinear_engine_t * createBilinearEngine() const override
  { 
    return new FormEngine(*this); 
  }

protected:
  explicit WeakForm(std::string const &className)
    : WeakFormBase(className)
  {}
};

template<typename Wf, typename Return_t>
class WeakForm<Wf, Return_t, typename std::enable_if<is_vector_type<Return_t>::value()>::type >
  : public WeakFormBase, public LinearFormTag
{
public:

  ~WeakForm() = default;

  template<typename ...Args>
  Return_t operator()(Args ...args) const
  {
    return static_cast<Wf const *>(this)->eval(args...);
  }

  class FormEngine : public Operator<IntegrationPoint, DynamicVectorR>::OperatorEngine
  {
  public:
    using base_t = typename Operator<IntegrationPoint, DynamicVectorR>::OperatorEngine;
    using var_t = WeakForm<Wf, Return_t>;
    using typename base_t::value_t;
    using typename base_t::index_t;

    FormEngine(var_t const &var)
      : _var(var)
    {}

    ~FormEngine() = default;

    value_t operator()(index_t const &idx) const override
    {
      return _var(idx);
    }

  private:
    var_t const &_var;
  };

  using typename WeakFormBase::linear_engine_t;
  using typename WeakFormBase::bilinear_engine_t;

  linear_engine_t * createLinearEngine() const override
  { 
    return new FormEngine(*this); 
  }

  bilinear_engine_t * createBilinearEngine() const override
  {
    return nullptr; 
  }

protected:
  explicit WeakForm(std::string const &className)
    : WeakFormBase(className)
  {}
};

/**
 * @brief      Class for linear weak form.
 *
 * @tparam     A         { description }
 * @tparam     F         { description }
 * @tparam     Return_t  { description }
 */
template<typename A, typename F, typename Return_t>
class LinearWeakForm
  : public WeakForm<LinearWeakForm<A, F, Return_t>, Return_t>
{
public:

  using return_t = Return_t;

  LinearWeakForm(A const &a, F const &f)
    : WeakForm<LinearWeakForm<A, F, Return_t>, Return_t>("LinearWeakForm")
    , _a(a), _f(f)
  {}

  template<typename ...Args>
  Return_t eval(Args ...args) const
  { 
    return _a(args...).transpose()*_f(args...);
  }

protected: 
  LinearWeakForm(A const &a, F const &f, std::string const &className)
    : WeakForm<LinearWeakForm<A, F, Return_t>, Return_t>(className)
    , _a(a), _f(f)
  {}

  A const &_a;
  F const &_f;
};

/**
 * @brief      Class for bilinear weak form.
 *
 * @tparam     A         { description }
 * @tparam     B         { description }
 * @tparam     Return_t  { description }
 */
template<typename A, typename B, typename Return_t>
class BilinearWeakForm
  : public LinearWeakForm<A, B, Return_t>
{
public:

  using LinearWeakForm<A, B, Return_t>::return_t;

  BilinearWeakForm(A const &a, B const &b)
    : LinearWeakForm<A, B, Return_t>(a, b, "BilinearWeakForm")
  {}

protected: 
  BilinearWeakForm(A const &a, B const &b, std::string const &className)
    : LinearWeakForm<A, B, Return_t>(a, b, className)
  {}
};

/**
 * @brief      Class for ad b weak form.
 *
 * @tparam     A         { description }
 * @tparam     D         { description }
 * @tparam     B         { description }
 * @tparam     Return_t  { description }
 */
template<typename A, typename D, typename B, typename Return_t>
class AdBWeakForm
  : public WeakForm<AdBWeakForm<A, D, B, Return_t>, Return_t>
{
public:

  using return_t = Return_t;

  AdBWeakForm(A const &a, D const &d, B const &b) 
    : WeakForm<AdBWeakForm<A, D, B, Return_t>, Return_t>("AdBWeakForm")
    , _a(a), _d(d), _b(b)
  {}

  template<typename ...Args>
  Return_t eval(Args ...args) const
  { 
    return _a(args...).transpose()*(_d(args...)*_b(args...));
  }

protected:
  AdBWeakForm(A const &a, D const &d, B const &b, std::string const &className)
    : WeakForm<AdBWeakForm<A, D, B, Return_t>, Return_t>(className)
    , _a(a), _d(d), _b(b) {}

  A const &_a;
  D const &_d;
  B const &_b;
};

/**
 * @brief      Class for bb weak form.
 *
 * @tparam     B         { description }
 * @tparam     Return_t  { description }
 */
template<typename B, typename Return_t>
class BBWeakForm : public BilinearWeakForm<B, B, Return_t>
{
public:
  using  BilinearWeakForm<B, B, Return_t>::return_t;

  BBWeakForm(B const &b)
    : BilinearWeakForm<B, B, Return_t>(b, b, "BBWeakForm")
  {}
};

/**
 * @brief      Class for bd b weak form.
 *
 * @tparam     B         { description }
 * @tparam     D         { description }
 * @tparam     Return_t  { description }
 */
template<typename B, typename D, typename Return_t>
class BdBWeakForm : public AdBWeakForm<B, D, B, Return_t>
{
public:

  using AdBWeakForm<B, D, B, Return_t>::return_t;

  BdBWeakForm(B const &b, D const &d, std::string const &className = "BdBWeakForm")
    : AdBWeakForm<B, D, B, Return_t>(b, d, b, className)
  {}
};

template<typename Wf, typename Tag>
class TaggedWeakForm : public Wf, public Tag
{
public:
  using Wf::return_t;

  template<typename ...Args>
  TaggedWeakForm(Args const&...args)
    : Wf(args...)
  {}
};

/**
 * @brief      Class for operator weak form.
 *
 * @tparam     A         { Operator class }
 * @tparam     Return_t  { value_t of operator class }
 */
template<class A>
class OperatorWeakForm
  : public WeakForm<OperatorWeakForm<A>, typename A::value_t>
{
public:

  using return_t = typename A::value_t;

  OperatorWeakForm(A const &a, std::string const &className = "OperatorWeakForm")
    : WeakForm<OperatorWeakForm<A>, return_t>(className)
    , _a(a)
  {}

  template<typename ...Args>
  typename A::value_t eval(Args ...args) const
  { 
    return _a(args...);
  }

private:
  A const &_a;
};

#endif //WeakForms_h