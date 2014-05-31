#include <string>
#include <iostream>
#include <cassert>
#include <cstdio>

#include "common.h"
#include "patch_loader.h"
#include "eson.h"

namespace mallie{

bool
PatchLoader::LoadESON(
    std::vector< bezier_patch<vector3> >& patch_array, const char* filename)
{
  // @todo { Use mmap() }
  std::vector<uint8_t> buf;

  FILE* fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Failed to load file: %s\n", filename);
  }

  fseek(fp, 0, SEEK_END);
  size_t len = ftell(fp);
  rewind(fp);
  buf.resize(len);
  len = fread(&buf[0], 1, len, fp);
  fclose(fp);

  eson::Value v;

  std::string err = eson::Parse(v, &buf[0]);
  if (!err.empty()) {
    std::cout << "Err: " << err << std::endl;
    exit(1);
  }

  //std::cout << "[LoadESON] # of shapes in .obj : " << shapes.size() << std::endl;

  if (!v.Has("num_vertices")) {
    fprintf(stderr, "Mallie:error\tmsg:\"num_vertices\" field not found.\n");
    return false;
  }
  int64_t num_vertices = v.Get("num_vertices").Get<int64_t>();
  printf("# of vertices: %lld\n", num_vertices);

#ifndef ENABLE_OSD_PATCH
  if (!v.Has("num_faces")) {
    fprintf(stderr, "Mallie:error\tmsg:\"num_faces\" field not found.\n");
    return false;
  }

  int64_t num_faces = v.Get("num_faces").Get<int64_t>();
  printf("# of faces   : %lld\n", num_faces);
#endif

  if (!v.Has("vertices")) {
    fprintf(stderr, "Mallie:error\tmsg:\"vertices\" field not found.\n");
    return false;
  }

  eson::Binary vertices_data = v.Get("vertices").Get<eson::Binary>();
  const float* vertices = reinterpret_cast<float*>(const_cast<uint8_t*>(vertices_data.ptr));

#ifndef ENABLE_OSD_PATCH
  if (!v.Has("faces")) {
    fprintf(stderr, "Mallie:error\tmsg:\"faces\" field not found.\n");
    return false;
  }

  eson::Binary faces_data = v.Get("faces").Get<eson::Binary>();
  const int* faces = reinterpret_cast<int*>(const_cast<uint8_t*>(faces_data.ptr));
#endif

  const float* facevarying_normals = NULL;
  if (v.Has("facevarying_normals")) {
    eson::Binary facevarying_normals_data = v.Get("facevarying_uvs").Get<eson::Binary>();
    facevarying_normals = reinterpret_cast<float*>(const_cast<uint8_t*>(facevarying_normals_data.ptr));
  }

  const float* facevarying_uvs = NULL;
  if (v.Has("facevarying_uvs")) {
    eson::Binary facevarying_uvs_data = v.Get("facevarying_uvs").Get<eson::Binary>();
    facevarying_uvs = reinterpret_cast<float*>(const_cast<uint8_t*>(facevarying_uvs_data.ptr));
  }

  const unsigned short* material_ids = NULL;
  if (v.Has("material_ids")) {
    eson::Binary material_ids_data = v.Get("material_ids").Get<eson::Binary>();
    material_ids = reinterpret_cast<unsigned short*>(const_cast<uint8_t*>(material_ids_data.ptr));
  }

#ifdef ENABLE_OSD_PATCH
  int64_t num_regular_patches = v.Get("num_regular_patches").Get<int64_t>();
  printf("# of regular patches   : %lld\n", num_regular_patches);

  eson::Binary regular_indices_data = v.Get("regular_indices").Get<eson::Binary>();
  const int* regular_indices = reinterpret_cast<int*>(const_cast<uint8_t*>(regular_indices_data.ptr));
#endif
    
#if 0

  // ESON -> Mesh
#ifdef ENABLE_OSD_PATCH
  mesh.numRegularPatches     = num_regular_patches;
  mesh.numVertices  = num_vertices;
  mesh.vertices = new real[num_vertices * 3];
  mesh.regularPatchIndices    = new unsigned int[num_regular_patches * 16];
  mesh.materialIDs = new unsigned int[num_regular_patches];
  mesh.bezierVertices = new real[num_regular_patches * 16 * 3];

  for (size_t i = 0; i < 3*num_vertices; i++) {
    mesh.vertices[i] = vertices[i];
  }

  for (size_t i = 0; i < 16*num_regular_patches; i++) {
    mesh.regularPatchIndices[i] = regular_indices[i];
  }

  // convert bspline to bezier
  const float Q[4][4] = {
      { 1.f/6.f, 4.f/6.f, 1.f/6.f, 0.f },
      { 0.f,     4.f/6.f, 2.f/6.f, 0.f },
      { 0.f,     2.f/6.f, 4.f/6.f, 0.f },
      { 0.f,     1.f/6.f, 4.f/6.f, 1.f/6.f } };

  real *pbv = mesh.bezierVertices;
  for (size_t i = 0; i < num_regular_patches; i++) {
      for (int j = 0; j < 4; j++) {
          for (int k = 0; k < 4; k++) {
              real3 H[4];
              for (int l = 0; l < 4; l++) {
                  H[l][0] = H[l][1] = H[l][2] = 0;
                  for (int m = 0; m < 4; m++) {
                      int vert = mesh.regularPatchIndices[i*16 + l*4 + m];
                      H[l][0] += Q[j][m] * mesh.vertices[vert*3 + 0];
                      H[l][1] += Q[j][m] * mesh.vertices[vert*3 + 1];
                      H[l][2] += Q[j][m] * mesh.vertices[vert*3 + 2];
                  }
              }
              real3 cp;
              cp[0] = cp[1] = cp[2] = 0;
              for (int m = 0; m < 4; m++) {
                  cp[0] += Q[k][m] * H[m][0];
                  cp[1] += Q[k][m] * H[m][1];
                  cp[2] += Q[k][m] * H[m][2];
              }
              *pbv++ = cp[0];
              *pbv++ = cp[1];
              *pbv++ = cp[2];
          }
      }
  }

  if (material_ids) {
    assert(0); // @todo
    //for (size_t i = 0; i < num_faces; i++) {
    //  mesh.materialIDs[i] = material_ids[i];
    //}
  } else {
    //for (size_t i = 0; i < num_faces; i++) {
    //  mesh.materialIDs[i] = 0; // 0 = default material.
    //}
  }
#else
  mesh.numFaces     = num_faces;
  mesh.numVertices  = num_vertices;
  mesh.vertices = new real[num_vertices * 3];
  mesh.faces    = new unsigned int[num_faces * 3];
  mesh.materialIDs = new unsigned int[num_faces];

  for (size_t i = 0; i < 3*num_vertices; i++) {
    mesh.vertices[i] = vertices[i];
  }

  for (size_t i = 0; i < 3*num_faces; i++) {
    mesh.faces[i] = faces[i];
  }

  if (material_ids) {
    for (size_t i = 0; i < num_faces; i++) {
      mesh.materialIDs[i] = material_ids[i];
    }
  } else {
    for (size_t i = 0; i < num_faces; i++) {
      mesh.materialIDs[i] = 0; // 0 = default material.
    }
  }
#endif

  mesh.facevarying_normals = NULL;
  mesh.facevarying_uvs = NULL;
  mesh.facevarying_tangents = NULL;
  mesh.facevarying_binormals = NULL;
#endif

  return true;
}
    
}
