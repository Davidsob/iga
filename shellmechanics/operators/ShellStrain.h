#pragma once

// #include "iga/operators/InterpolatedVariable.h"
#include "iga/Values.h"
#include "iga/operators/VectorizedVariable.h"

#include "LinearMembraneStrainOperator.h"
#include "LinearBendingStrainOperator.h"
#include "LinearTransverseStrainOperator.h"


using LinearMembraneStrain
    = BinaryMappedValue<
        LinearMembraneStrainOperator,
        VectorizedVariable<SixDofDisplacementVariable>,
        Multiply,
        StaticVectorR<6>
        >;

using LinearBendingStrain
    = BinaryMappedValue<
        LinearBendingStrainOperator,
        VectorizedVariable<SixDofDisplacementVariable>,
        Multiply,
        StaticVectorR<6>
        >;

using LinearTransverseShearStrain
    = BinaryMappedValue<
        LinearTransverseStrainOperator,
        VectorizedVariable<SixDofDisplacementVariable>,
        Multiply,
        StaticVectorR<6>
        >;
