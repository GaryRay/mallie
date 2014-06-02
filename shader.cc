#include <cassert>

#include "shader.h"
#include "scene.h"
#include "vector3.h"

using namespace mallie;

//
// Some predefined shaders.
//
void
ShowNormal(float rgba[4], const Scene& scene, const Intersection& isect, const Ray& ray)
{
  rgba[0] = 0.5 * isect.normal[0] + 0.5;
  rgba[1] = 0.5 * isect.normal[1] + 0.5;
  rgba[2] = 0.5 * isect.normal[2] + 0.5;
  rgba[3] = 1.0;
}

void
EyeDotN(float rgba[4], const Scene& scene, const Intersection& isect, const Ray& ray)
{
  real3 I;
#if 0
  I[0] = scene.GetCamera().eye_[0] - isect.position[0];
  I[1] = scene.GetCamera().eye_[1] - isect.position[1];
  I[2] = scene.GetCamera().eye_[2] - isect.position[2];
#else
  I = ray.dir.neg();
#endif
  I.normalize();

  real IdotN = vdot(I, isect.normal); 
  IdotN = std::max(0.2, IdotN);

  rgba[0] = IdotN;
  rgba[1] = IdotN;
  rgba[2] = IdotN;
  rgba[3] = 1.0;
}

void
YourShader(float rgba[4], const Scene& scene, const Intersection& isect)
{

}

Shader::Shader()
{
  currentShaderIndex_ = 0;
  maxRegisteredShaders_ = 0;

  shaderList_.resize(kMaxShaders);
  for (size_t i = 0; i < kMaxShaders; i++) {
    shaderList_[i] = NULL;
  }

  RegisterShader(0, ShowNormal);
  RegisterShader(1, EyeDotN);

  // @note { Add your shader here! }
}

void Shader::RegisterShader(
int index, ShaderFunction func)
{
  assert(index < kMaxShaders);
  shaderList_[index] = func;

  if ((index+1) > maxRegisteredShaders_) {
    maxRegisteredShaders_ = (index+1);
  }
}

void Shader::SetShader(int index)
{
  assert(index < kMaxShaders);

  currentShaderIndex_ = index;
}

void Shader::Shade(float rgba[4], const Scene& scene, const Intersection& isect, const Ray& ray)
{
  if (shaderList_[currentShaderIndex_]) {
    shaderList_[currentShaderIndex_](rgba, scene, isect, ray);
  }
}
