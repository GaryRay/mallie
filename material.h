#ifndef __MALLIE_MATERIAL_H__
#define __MALLIE_MATERIAL_H__

#include "common.h"

struct Material {
  real3 diffuse;
  real3 reflection;
  real  reflection_glossiness;
  int   fresnel; // 1 = Fresnel reflection on, 0 = off.
  float ior; // index of reflection.
  int   id;
};

#endif // __MALLIE_MATERIAL_H__
