#include "TensorUtils.h"

StaticMatrixR<3,3> const TensorUtils::eye = StaticMatrixR<3,3>::Identity();

StaticMatrixR<6,6> const TensorUtils::eye6 = StaticMatrixR<6,6>::Identity();

StaticMatrixR<9,9> const TensorUtils::eye9 = StaticMatrixR<9,9>::Identity();

const size_t TensorUtils::ik[] = {0, 1, 2, 1, 0, 0};

const size_t TensorUtils::jl[] = {0, 1, 2, 2, 2, 1};

