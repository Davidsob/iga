#pragma once

#include "base/Values.h"

#include "iga/operators/Multiply.h"
#include "iga/operators/VectorizedVariable.h"

#include "LinearMembraneStrainOperator.h"
#include "LinearBendingStrainOperator.h"
#include "LinearTransverseStrainOperator.h"

#include "MembraneGreensStrainOperator.h"
#include "BendingGreensStrainOperator.h"
#include "TransverseGreensStrainOperator.h"


using LinearMembraneStrain
    = BinaryMappedValue<
        LinearMembraneStrainOperator,
        VectorizedVariable<SixDofDisplacementVariable>,
        Multiply,
        StaticVectorR<3>
        >;

using MembraneGreensStrain
    = BinaryMappedValue<
        MembraneGreensStrainOperator,
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

using BendingGreensStrain
    = BinaryMappedValue<
        BendingGreensStrainOperator,
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

using TransverseShearGreensStrain
    = BinaryMappedValue<
        TransverseGreensStrainOperator,
        VectorizedVariable<SixDofDisplacementVariable>,
        Multiply,
        StaticVectorR<2>
        >;
