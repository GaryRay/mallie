#ifndef __MESH_H__
#define __MESH_H__

#include <cstdio>
#include "common.h"

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
  size_t numRegularPatches;
  unsigned int *regularPatchIndices; /// 16 x numRegularPatches
  real *bezierVertices;        /// [xyz] * 16 * numRegularPatches
#endif

} Mesh;

#endif // __MESH_H__
