#ifndef __PATCH_LOADER_H__
#define __PATCH_LOADER_H__

#include "vector3.h"
#include "bezier/bezier_patch.hpp"
#include <vector>

namespace mallie{

class
PatchLoader
{
 public:
  PatchLoader();

  ///< Load wavefront obj data from a file.
  ///< Allocated memory for the mesh must be free'ed by the application.
  //static bool LoadObj(Mesh& mesh, const char* filename);

  ///< Load ESON mesh from a file.
  ///< Allocated memory for the mesh must be free'ed by the application.
  ///< Implicitly loads material infos whose name is basename(filename) + '.material.json'
    static bool LoadESON(std::vector< bezier_patch<vector3> >& path_array, const char* filename);

};
    
}

#endif  // __PATCH_LOADER_H__
// vim:set sw=2 ts=2 expandtab:
