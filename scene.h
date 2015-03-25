#ifndef __MALLIE_SCENE_H__
#define __MALLIE_SCENE_H__

#include <cassert>
#include <string>

#include "bvh_accel.h"
#include "material.h"
#include "camera.h"
#include "shader.h"
#include "texture.h"

#include "bezier/patch_accel.h"

namespace mallie {

typedef enum {
  PRIMITIVE_POLYGON,
} PrimitiveType;

typedef struct {
  PrimitiveType type;
  BVHAccel *accel;
} Primitive;

//< Scene graph node class
class Node {
public:
  Node() {}
  ~Node() {}

  //< Recompute transform matrix(includes child nodes).
  void UpdateTransform();

private:
  double xform_[16];               //< Transform of this node.
  double xformInverse_;            //< Inverse transform
  double xformInverse33_;          //< Upper 3x3 of inverse transform
  double xformInverseTranspose33_; //< Upper 3x3 of inverse transpose transform
  std::vector<Node *> children_;   //< Children
};

//< Scene class. Usually only one instance of Scene class in rendering.
class Scene {
public:
  Scene();
  ~Scene();

  // Fixme.
  bool Init(const std::string &objFilename, const std::string &esonFilename,
            const std::string &materialFilename,
            const std::string &ptexFilename,
            const std::string &envmapFilename, const std::string &envmapCoord,
            double sceneScale = 1.0,
            int numLeafPrimitives = 16);

  bool Trace(Intersection &isect, Ray &ray) const;

  void BoundingBox(real3 &bmin, real3 &bmax);

  real3 GetBackgroundRadiance(real3 &dir);

  const Material &GetMaterial(int matID) const {
    assert(matID < materials_.size());
    return materials_[matID];
  }

  size_t AddMaterial(const Material& mat) {
    size_t idx = materials_.size();
    materials_.push_back(mat);
    return idx;
  }

  void SetCamera(Camera &cam) { camera_ = cam; }

  const Camera &GetCamera() const { return camera_; }

  Shader &GetShader() { return shader_; }

  bool LoadEnvMap(const std::string &filename, const std::string& coord);

  // JSON format material loader.
  bool LoadMaterial(const std::string& filename);

  const Texture &GetEnvMap() const { return envMap_; }

  const PTexture *GetPTexture() const { return ptex_; }

protected:
  BVHAccel accel_;
  Mesh mesh_;
  std::vector<Material> materials_;
  Camera camera_;
  Shader shader_;

  Texture envMap_; // Environment map(Assume angular map coordinte)

  PTexture *ptex_;

  PatchAccel *patch_accel_;
};
}

#endif // __MALLIE_SCENE_H__
