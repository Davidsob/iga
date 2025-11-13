#pragma once

class GeometricStiffnessOperator_ida 
{
public:
  using value_t = DynamicMatrixR;
  GeometricStiffnessOperator_ida(){}
  ~GeometricStiffnessOperator_ida() = default;

  template<typename Point>
  value_t operator()(Point const &p) const
  {
    auto const shape = p.mapper.shape(p.para);
    auto const dof = shape.size();
    DynamicMatrixR b(DynamicMatrixR::Zero(6,6*dof));

    for (int i = 0; i < dof; i++)
    {
      b(0,i*6+3) = shape(i);
      b(1,i*6+3) = shape(i);

      b(2,i*6+4) = shape(i);
      b(3,i*6+4) = shape(i);

      b(4,i*6+5) = shape(i);
      b(5,i*6+5) = shape(i);
    }

    return b;
  }
};
