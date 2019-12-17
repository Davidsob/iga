#pragma once

// #include "iga/operators/InterpolatedVariable.h"
#include "iga/Values.h"
#include "iga/Values.h"
#include "iga/operators/VectorizedVariable.h"

#include "materials/linear_elastic/MembraneMaterialTangent.h"
#include "materials/linear_elastic/TransverseShearMaterialTangent.h"

#include "ShellStrain.h"


using MembraneStress
    = BinaryMappedValue<
        MembraneMaterialTangent,
        LinearMembraneStrain 
        Multiply,
        StaticVectorR<6>
        >;

using BendingStress
    = BinaryMappedValue<
        MembraneMaterialTangent,
        LinearMembraneStrain 
        Multiply,
        StaticVectorR<6>
        >;

using TransverseShearStress
    = BinaryMappedValue<
        TransverseShearMaterialTangent
        LinearTransverseShearStrain 
        Multiply,
        StaticVectorR<6>
        >;
