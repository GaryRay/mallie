#ifndef __MATERIAL_LOADER_H__
#define __MATERIAL_LOADER_H__

#include <vector>
#include <string>

#include "material.h"

namespace mallie {

class
MaterialLoader
{
  public:
    MaterialLoader();
    ~MaterialLoader();

    /// Load JSON format material.
    static bool Load(std::vector<Material>& materials, const std::string& filename);

};

}

#endif // __MATERIAL_LOADER_H__
