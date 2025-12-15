#pragma once

#if USE_DOUBLE
using FloatL = double;
#else
using FloatL = float;
#endif

#if SOA
#define INDEX(arrayID,sturctureID,arraySize,structureSize) (arrayID+sturctureID*(arraySize))
#else
#define INDEX(arrayID,sturctureID,arraySize,structureSize) (arrayID*(structureSize)+sturctureID)
#endif

