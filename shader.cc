#include <cassert>

#include "shader.h"
#include "scene.h"
#include "vector3.h"

using namespace mallie;

namespace {

vector3 reflect(const vector3 &in, const vector3 &n) {
  float d = dot(in, n);
  return in - n * (2.0 * d);
}

vector3 refract(bool &tir, const vector3 &in, const vector3 &n, float eta) {
  vector3 ret;
  vector3 N;
  double e = eta;
  double cos1 = dot(in, n);
  if (cos1 < 0.0) { // entering
    N = n;
  } else { // outgoing
    cos1 = -cos1;
    e = 1.0f / eta;
    N = -n;
  }

  double k = 1.0 - (e * e) * (1.0 - cos1 * cos1);
  if (k <= 0.0) {
    // Toral internal reflection.
    ret = reflect(in, n);
    tir = true;
    ret.normalize();
    return ret;
  }

  k = -e * cos1 - sqrt(k);

  tir = false;
  ret = k * N + e * in;
  ret.normalize();

  return ret;
}

void fresnel(vector3 &refl, vector3 &refr, float &kr, float &kt,
             const vector3 &in, const vector3 &n, float eta) {
  float d = dot(in, n);

  refl = reflect(in, n);

  bool tir;
  refr = refract(tir, in, n, eta);

  if (tir) {
    kr = 1.0;
    kt = 0.0;
    return;
  }

  float cos_r = dot(refl, n);
  float cos_t = -dot(refr, n);

  float rp = (cos_t - eta * cos_r) / (cos_t + eta * cos_r);
  float rs = (cos_r - eta * cos_t) / (cos_r + eta * cos_t);
  kr = (rp * rp + rs * rs) * 0.5f;

  if (kr < 0.0f)
    kr = 0.0f;
  if (kr > 1.0f)
    kr = 1.0f;

  kt = 1.0f - kr;
}
}

//
// Some predefined shaders.
//
void ShowNormal(float rgba[4], const Scene &scene, const Intersection &isect,
                const Ray &ray) {
  rgba[0] = 0.5 * isect.normal[0] + 0.5;
  rgba[1] = 0.5 * isect.normal[1] + 0.5;
  rgba[2] = 0.5 * isect.normal[2] + 0.5;
  rgba[3] = 1.0;
}

void EyeDotN(float rgba[4], const Scene &scene, const Intersection &isect,
             const Ray &ray) {
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

void EnvShader(float rgba[4], const Scene &scene, const Intersection &isect,
               const Ray &ray) {
  vector3 in, n;

  in[0] = ray.dir[0];
  in[1] = ray.dir[1];
  in[2] = ray.dir[2];
  in.normalize();

  n[0] = isect.normal[0];
  n[1] = isect.normal[1];
  n[2] = isect.normal[2];
  n.normalize();

  vector3 r = reflect(in, n);

  float dir[3];
  dir[0] = r[0];
  dir[1] = r[1];
  dir[2] = r[2];

  if (scene.GetEnvMap().IsValid()) {
    if (scene.GetEnvMap().coordinate() == Texture::COORDINATE_LONGLAT) {
      LongLatMapSampler::Sample(rgba, dir, &(scene.GetEnvMap()));
    } else if (scene.GetEnvMap().coordinate() == Texture::COORDINATE_ANGULAR) {
      AngularMapSampler::Sample(rgba, dir, &(scene.GetEnvMap()));
    } else {
      rgba[0] = 1.0;
      rgba[1] = 1.0;
      rgba[2] = 1.0;
      rgba[3] = 1.0;
    }


    rgba[0] *= 3.14;
    rgba[1] *= 3.14;
    rgba[2] *= 3.14;

  } else {

    rgba[0] = 0.5 * r[0] + 0.5;
    rgba[1] = 0.5 * r[1] + 0.5;
    rgba[2] = 0.5 * r[2] + 0.5;
    rgba[3] = 1.0;
  }
}

void PtexShader(float rgba[4], const Scene &scene, const Intersection &isect,
                const Ray &ray) {

  float illum[4];
  EyeDotN(illum, scene, isect, ray);

  scene.GetPTexture()->Eval(rgba, 0, 3, isect.faceID, isect.u, isect.v, 0, 0, 0, 0);

  float scale = 4.0f;
  rgba[0] *= scale*illum[0];
  rgba[1] *= scale*illum[1];
  rgba[2] *= scale*illum[2];
  rgba[3] = 1;
}

void YourShader(float rgba[4], const Scene &scene, const Intersection &isect,
                const Ray &ray) {}

Shader::Shader() {
  currentShaderIndex_ = 0;
  maxRegisteredShaders_ = 0;

  shaderList_.resize(kMaxShaders);
  for (size_t i = 0; i < kMaxShaders; i++) {
    shaderList_[i] = NULL;
  }

  RegisterShader(0, ShowNormal);
  RegisterShader(1, EyeDotN);
  RegisterShader(2, EnvShader);
  RegisterShader(3, PtexShader);

  // @note { Add your shader here! }
}

void Shader::RegisterShader(int index, ShaderFunction func) {
  assert(index < kMaxShaders);
  shaderList_[index] = func;

  if ((index + 1) > maxRegisteredShaders_) {
    maxRegisteredShaders_ = (index + 1);
  }
}

void Shader::SetShader(int index) {
  assert(index < kMaxShaders);

  currentShaderIndex_ = index;
}

void Shader::Shade(float rgba[4], const Scene &scene, const Intersection &isect,
                   const Ray &ray) {
  if (shaderList_[currentShaderIndex_]) {
    shaderList_[currentShaderIndex_](rgba, scene, isect, ray);
  }
}
