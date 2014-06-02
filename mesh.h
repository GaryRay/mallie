#ifndef __MESH_H__
#define __MESH_H__

#include <cstdio>
#include "common.h"

struct FarPatchParam {
///  Field      | Bits | Content                                              
///  -----------|:----:|------------------------------------------------------
///  level      | 4    | the subdivision level of the patch                   
///  nonquad    | 1    | whether the patch is the child of a non-quad face    
///  rotation   | 2    | patch rotations necessary to match CCW face-winding  
///  v          | 10   | log2 value of u parameter at first patch corner      
///  u          | 10   | log2 value of v parameter at first patch corner      
///  reserved1  | 5    | padding                                              
    unsigned int faceIndex:32; // Ptex face index
    struct BitField {
        unsigned int field:32;
    } bitfield;
};

typedef struct {
  size_t numVertices;
#ifndef ENABLE_OSD_PATCH
  size_t numFaces;
#endif
  real *vertices;              /// [xyz] * numVertices
  real *facevarying_normals;   /// [xyz] * 3(triangle) * numFaces
  real *facevarying_tangents;  /// [xyz] * 3(triangle) * numFaces
  real *facevarying_binormals; /// [xyz] * 3(triangle) * numFaces
  real *facevarying_uvs;       /// [xyz] * 3(triangle) * numFaces
#ifndef ENABLE_OSD_PATCH
  unsigned int *faces;         /// triangle x numFaces
#endif
  unsigned int *materialIDs;   /// index x numFaces

#ifdef ENABLE_OSD_PATCH
  size_t numBezierPatches;
  real *bezierVertices;        /// [xyz] * 16 * numBezierPatches
  FarPatchParam *patchParams;  /// [FarPatchParam] * numBezierPatches
#endif

} Mesh;

#endif // __MESH_H__
