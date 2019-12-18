#pragma once

#include "base/Values.h"

#include "iga/operators/Multiply.h"
#include "iga/operators/VectorizedVariable.h"

#include "materials/linear_elastic/MembraneMaterialTangent.h"
#include "materials/linear_elastic/TransverseShearMaterialTangent.h"

#include "ShellStrain.h"


using MembraneStress
    = BinaryMappedValue<
        MembraneMaterialTangent,
        LinearMembraneStrain,
        Multiply,
        StaticVectorR<3>
        >;

using BendingStress
    = BinaryMappedValue<
        MembraneMaterialTangent,
        LinearBendingStrain,
        Multiply,
        StaticVectorR<3>
        >;

using TransverseShearStress
    = BinaryMappedValue<
        TransverseShearMaterialTangent,
        LinearTransverseShearStrain,
        Multiply,
        StaticVectorR<2>
        >;
