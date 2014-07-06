#include <cassert>

#include "shader.h"
#include "scene.h"
#include "vector3.h"

using namespace mallie;

namespace mallie {

real3 PathTrace(const Scene &scene, const Camera &camera,
                const Intersection& s, const Ray& inRay); // render.cc

double randomreal(); // render.cc

}

namespace {

float vavg(real3 x)
{
  return (x[0] + x[1] + x[2]) / 3;
}


real3 vclamp01(real3 x)
{
  real3 ret;
  ret[0] = x[0];
  ret[1] = x[1];
  ret[2] = x[2];
}

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

void GenerateBasis(real3 &tangent, real3 &binormal,
                          const real3 &normal) {
  // Find the minor axis of the vector
  int i;
  int index = -1;
  double minval = 1.0e+6;
  double val = 0;

  for (int i = 0; i < 3; i++) {
    val = fabsf(normal[i]);
    if (val < minval) {
      minval = val;
      index = i;
    }
  }

  if (index == 0) {

    tangent.x = 0.0;
    tangent.y = -normal.z;
    tangent.z = normal.y;
    tangent.normalize();

    binormal = vcross(tangent, normal);
    binormal.normalize();

  } else if (index == 1) {

    tangent.x = -normal.z;
    tangent.y = 0.0;
    tangent.z = normal.x;
    tangent.normalize();

    binormal = vcross(tangent, normal);
    binormal.normalize();

  } else {

    tangent.x = -normal.y;
    tangent.y = normal.x;
    tangent.z = 0.0;
    tangent.normalize();

    binormal = vcross(tangent, normal);
    binormal.normalize();
  }
}

// (Modified) Ward glossy BRDF
// http://www.graphics.cornell.edu/~bjw/wardnotes.pdf
// Some from OpenShadingLanguage
void WardBRDF(
    real3* omega_in, // output
    float* pdf,     // output
    float* weight,  // output
    float ax,
    float ay,
    real3  omega_out,
    real3  normal,               // shading normal
    real3  geometric_normal)     // geometric normal
{
    float cosNO = vdot(normal, omega_out);

    (*pdf) = 0.0f;
    (*weight) = 0.0f;
    (*omega_in)[0] = 0.0f;
    (*omega_in)[1] = 0.0f;
    (*omega_in)[2] = 0.0f;

    if (cosNO > 0.0f) {
        // @todo { Supply tangent vector for true aniso-brdf. }
        real3 tangent, binormal;

        GenerateBasis(tangent, binormal, normal);

        float randu = randomreal();
        float randv = randomreal();

        float alphaRatio = ay / ax;
        float cosPhi, sinPhi;

        if (randu < 0.25f) {
            float val = 4 * randu;
            float tanPhi = alphaRatio * tanf((float) M_PI_2 * val);
            cosPhi = 1 / sqrtf(1 + tanPhi * tanPhi);
            sinPhi = tanPhi * cosPhi;
        } else if (randu < 0.5) {
            float val = 1 - 4 * (0.5f - randu);
            float tanPhi = alphaRatio * tanf((float) M_PI_2 * val);
            // phi = (float) M_PI - phi;
            cosPhi = -1 / sqrtf(1 + tanPhi * tanPhi);
            sinPhi = -tanPhi * cosPhi;
        } else if (randu < 0.75f) {
            float val = 4 * (randu - 0.5f);
            float tanPhi = alphaRatio * tanf((float) M_PI_2 * val);
            //phi = (float) M_PI + phi;
            cosPhi = -1 / sqrtf(1 + tanPhi * tanPhi);
            sinPhi = tanPhi * cosPhi;
        } else {
            float val = 1 - 4 * (1 - randu);
            float tanPhi = alphaRatio * tanf((float) M_PI_2 * val);
            // phi = 2 * (float) M_PI - phi;
            cosPhi = 1 / sqrtf(1 + tanPhi * tanPhi);
            sinPhi = -tanPhi * cosPhi;
        }

        // eq. 6
        // we take advantage of cos(atan(x)) == 1/sqrt(1+x^2)
        //                  and sin(atan(x)) == x/sqrt(1+x^2)
        float thetaDenom = (cosPhi * cosPhi) / (ax * ax) + (sinPhi * sinPhi) / (ay * ay);
        float tanTheta2 = -logf(1 - randv) / thetaDenom;
        float cosTheta  = 1 / sqrtf(1 + tanTheta2);
        float sinTheta  = cosTheta * sqrtf(tanTheta2);

        real3 h; // already normalized becaused expressed from spherical coordinates
        h.x = sinTheta * cosPhi;
        h.y = sinTheta * sinPhi;
        h.z = cosTheta;
        // compute terms that are easier in local space
        float dotx = h.x / ax;
        float doty = h.y / ay;
        float dotn = h.z;
        // transform to world space
        h = h.x * tangent + h.y * binormal + h.z * normal;
        // generate the final sample
        float oh = vdot(h, omega_out);
        //omega_in->x = 2 * oh * h.x - omega_out.x;
        //omega_in->y = 2 * oh * h.y - omega_out.y;
        //omega_in->z = 2 * oh * h.z - omega_out.z;
        omega_in->x = omega_out.x - 2 * oh * h.x;
        omega_in->y = omega_out.y - 2 * oh * h.y;
        omega_in->z = omega_out.z - 2 * oh * h.z;

        float ng_dot_wi = vdot(geometric_normal, (*omega_in));
        if (ng_dot_wi > 0.0f) {
            float cosNI = vdot(normal, (*omega_in));
            if (cosNI > 0.0f) {
                // eq. 9
                float exp_arg = (dotx * dotx + doty * doty) / (dotn * dotn);
                float denom = 4 * (float) M_PI * ax * ay * oh * dotn * dotn * dotn;
                (*pdf) = expf(-exp_arg) / denom;
                // compiler will reuse expressions already computed
                denom = (4 * (float) M_PI * ax * ay * sqrtf(cosNO * cosNI));
                float power = cosNI * expf(-exp_arg) / denom;
                (*weight) = power;
                //domega_in_dx = (2 * m_N.dot(domega_out_dx)) * m_N - domega_out_dx;
                //domega_in_dy = (2 * m_N.dot(domega_out_dy)) * m_N - domega_out_dy;
                // Since there is some blur to this reflection, make the
                // derivatives a bit bigger. In theory this varies with the
                // roughness but the exact relationship is complex and
                // requires more ops than are practical.
                //domega_in_dx *= 10;
                //domega_in_dy *= 10;
            }
        }
    }
}


}

//
// Some predefined shaders.
//

// Physically-based shader
void PBS(float rgba[4], const Scene &scene, const Intersection &isect,
                const Ray &ray) {

  // Currently available: diffuse + reflection(+glossy reflection)

  int matID = isect.matID;
  const Material& mat = scene.GetMaterial(matID);

  // Preserve Energy conservation for each channel.
  real3 diffuse = mat.diffuse;
  real3 reflection = mat.reflection;
  float reflectionGlossiness = mat.reflection_glossiness;

  real3 one(1.0, 1.0, 1.0);
  real3 ksRGB = reflection;
  real3 kdRGB = vclamp01((one - ksRGB) * diffuse);

  float ks = vavg(ksRGB); ks = std::min(1.0f, std::max(0.0f, ks));
  float kd = vavg(kdRGB); kd = std::min(1.0f, std::max(0.0f, kd));

  real3 kdRet(0.0, 0.0, 0.0);
  real3 ksRet(0.0, 0.0, 0.0);
  if (kd > 0.0) {
    kdRet[0] = kdRGB[0];
    kdRet[1] = kdRGB[1];
    kdRet[2] = kdRGB[2];
  }

  real3 in, n;

  in[0] = ray.dir[0];
  in[1] = ray.dir[1];
  in[2] = ray.dir[2];
  in.normalize();

  n[0] = isect.normal[0];
  n[1] = isect.normal[1];
  n[2] = isect.normal[2];
  n.normalize();


  // reflection
  if (ks > 0.0) {
    real3 r;

    float weight = 1.0;

    if (reflectionGlossiness < 1.0) {
      // glossy reflection. WardBRDF
      
      float pdf;

      // larget = sharper.
      float ax = 1000.0f * reflectionGlossiness; // isotropic
      float ay = 1000.0f * reflectionGlossiness;

      real3 wi = in.neg();

      WardBRDF(&r, &pdf, &weight, ax, ay, wi, n, n);
      //printf("w = %f, r = %f, %f, %f, n = %f, %f, %f\n",
      //  weight, r[0], r[1], r[2], n[0], n[1], n[2]);

    } else {
      // perfect specular.
      vector3 in_0(in[0], in[1], in[2]);
      vector3 n_0(n[0], n[1], n[2]);
      vector3 rr = reflect(in_0, n_0);
      r[0] = rr[0];
      r[1] = rr[1];
      r[2] = rr[2];
    }

    float dir[3];
    dir[0] = r[0];
    dir[1] = r[1];
    dir[2] = r[2];

    if ((weight > 0.0) && scene.GetEnvMap().IsValid()) {
      float envcol[4];
      if (scene.GetEnvMap().coordinate() == Texture::COORDINATE_LONGLAT) {
        LongLatMapSampler::Sample(envcol, dir, &(scene.GetEnvMap()));
      } else if (scene.GetEnvMap().coordinate() == Texture::COORDINATE_ANGULAR) {
        AngularMapSampler::Sample(envcol, dir, &(scene.GetEnvMap()));
      } else {
        // @todo
        envcol[0] = 1.0;
        envcol[1] = 1.0;
        envcol[2] = 1.0;
      }

      ksRet[0] = envcol[0] * ksRGB[0];
      ksRet[1] = envcol[1] * ksRGB[1];
      ksRet[2] = envcol[2] * ksRGB[2];

    } else {

      // @todo: Sunsky?
      ksRet[0] = 0.0;
      ksRet[1] = 0.0;
      ksRet[2] = 0.0;
    }
  }

  // @fixme.
  rgba[0] = kdRet[0] + ksRet[0];
  rgba[1] = kdRet[1] + ksRet[1];
  rgba[2] = kdRet[2] + ksRet[2];
  rgba[3] = 1.0;
}

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

void PathTraceShader(float rgba[4], const Scene &scene, const Intersection &isect,
                const Ray &ray) {

  real3 col = PathTrace(scene, scene.GetCamera(), isect, ray);

  // ???: negate makes us happy
  rgba[0] = 1.0f - col[0];
  rgba[1] = 1.0f - col[1];
  rgba[2] = 1.0f - col[2];
  rgba[3] = 1.0;
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
  //RegisterShader(0, PBS);
  RegisterShader(1, EyeDotN);
  RegisterShader(2, EnvShader);
  RegisterShader(3, PtexShader);
  RegisterShader(4, PathTraceShader);
  RegisterShader(5, PBS);

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
