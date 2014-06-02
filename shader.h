#ifndef __MALLIE_SHADER_H__
#define __MALLIE_SHADER_H__

#include <vector>
#include "common.h"
#include "intersection.h"

namespace mallie {

class Scene;

const int kMaxShaders = 16;

// Simple shader class
class Shader {
public:
  typedef void (*ShaderFunction)(float rgba[4], const Scene &scene,
                                 const Intersection &isect, const Ray &ray);

  Shader();
  ~Shader() {};

  /// Register custom shader function.
  void RegisterShader(int index, ShaderFunction func);

  /// Specify shader.
  void SetShader(int index);

  int NumShaders() const { return maxRegisteredShaders_; }

  /// Call shader specified by SetShader()
  void Shade(float rgba[4], const Scene &scene, const Intersection &isect,
             const Ray &ray);

private:
  std::vector<ShaderFunction> shaderList_;
  int currentShaderIndex_;
  int maxRegisteredShaders_;
};
}

#endif // __MALLIE_SHADER_H__
