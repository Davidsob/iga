#pragma once

class GeometricStiffnessOperator_grada 
{
public:
  using value_t = DynamicMatrixR;
  GeometricStiffnessOperator_grada(){}
  ~GeometricStiffnessOperator_grada() = default;

  template<typename Point>
  value_t operator()(Point const &p) const
  {
    auto const grad = p.mapper.grad(p.para);
    auto const dof = grad.cols();
    DynamicMatrixR b(DynamicMatrixR::Zero(6,6*dof));

    for (int i = 0; i < dof; i++)
    {
      b(0,i*6+3) = grad(0,i);
      b(1,i*6+3) = grad(1,i);

      b(2,i*6+4) = grad(0,i);
      b(3,i*6+4) = grad(1,i);

      b(4,i*6+5) = grad(0,i);
      b(5,i*6+5) = grad(1,i);
    }

    return b;
  }
};
