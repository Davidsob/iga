#pragma once

class GeometricStiffnessOperator_gradu 
{
public:
  using value_t = DynamicMatrixR;
  GeometricStiffnessOperator_gradu(){}
  ~GeometricStiffnessOperator_gradu() = default;

  template<typename Point>
  value_t operator()(Point const &p) const
  {
    auto const grad = p.mapper.grad(p.para);
    auto const dof = grad.cols();
    DynamicMatrixR b(DynamicMatrixR::Zero(6,6*dof));

    for (int i = 0; i < dof; i++)
    {
      b(0,i*6) = grad(0,i);
      b(1,i*6) = grad(1,i);

      b(2,i*6+1) = grad(0,i);
      b(3,i*6+1) = grad(1,i);

      b(4,i*6+2) = grad(0,i);
      b(5,i*6+2) = grad(1,i);
    }

    return b;
  }
};
