#pragma once

#include "base/Values.h"

#include "iga/operators/Multiply.h"
#include "iga/operators/VectorizedVariable.h"

#include "LinearMembraneStrainOperator.h"
#include "LinearBendingStrainOperator.h"
#include "LinearTransverseStrainOperator.h"


using LinearMembraneStrain
    = BinaryMappedValue<
        LinearMembraneStrainOperator,
        VectorizedVariable<SixDofDisplacementVariable>,
        Multiply,
        StaticVectorR<3>
        >;

using LinearBendingStrain
    = BinaryMappedValue<
        LinearBendingStrainOperator,
        VectorizedVariable<SixDofDisplacementVariable>,
        Multiply,
        StaticVectorR<3>
        >;

using LinearTransverseShearStrain
    = BinaryMappedValue<
        LinearTransverseStrainOperator,
        VectorizedVariable<SixDofDisplacementVariable>,
        Multiply,
        StaticVectorR<2>
        >;
