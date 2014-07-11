#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include <cstdio>
#include <iostream>
#include <ctime>
#include <algorithm>
#include <limits>

#if defined(_WIN32) && !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES
#endif
#include <cmath>

#include "hashgrid.h"
#include "render.h"
#include "camera.h"
#include "timerutil.h"
#include "scene.h"
#include "script_engine.h"
#include "shader.h"
#include "prim-plane.h"
#include "feature-lines.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef ENABLE_PTEX
#include <Ptexture.h>
#endif

#ifdef _WIN32
#define THREAD_TLS __declspec(thread)
#else // Assume gcc-like compiler
#define THREAD_TLS __thread
#endif

extern bool CheckRenderCancel();

namespace mallie {

const double kFar = 1.0e+30;
const double kEPS = 1.0e-3;
const int kMaxPathLength = 16;
const int kMinPathLength = 2;

const int kTileSize = 8;

const int kPtexMaxMem = 1024*1024; // @fixme.

struct PathVertex {
  real3 P;          // Position
  real3 N;          // Normal
  real3 wi;         // Incident vector
  real3 throughput; // Path throughput(RGB)
  int matID;        // Material ID
};

typedef std::vector<PathVertex> Path;

// HACK: SGA14 TechBrief
Plane gPlane;
bool gPlaneInitialied = false;

void init_plane(const real3& sceneBMin, const real3& sceneBMax, int matID)
{
  float ymin  = sceneBMin[1];
  float ysize = sceneBMax[1] - sceneBMin[1];
  //plane.set(0, 1, 0, -(zmin - zsize * plane_distscale));  
  gPlane.Set(0, 1, 0, -(ymin - ysize * 0.0001), matID);  

  gPlaneInitialied = true;
}

unsigned int gSeed[1024][4];

void init_randomreal(void) {
#if _OPENMP
  assert(omp_get_max_threads() < 1024);

  for (int i = 0; i < omp_get_max_threads(); i++) {
    gSeed[i][0] = 123456789 + i;
    gSeed[i][1] = 362436069;
    gSeed[i][2] = 521288629;
    gSeed[i][3] = 88675123;
  }
#else
#endif
}

double randomreal(void) {
// xorshift RNG
#ifdef _OPENMP
  int tid = omp_get_thread_num();
  unsigned int x = gSeed[tid][0];
  unsigned int y = gSeed[tid][1];
  unsigned int z = gSeed[tid][2];
  unsigned int w = gSeed[tid][3];
  unsigned t = x ^ (x << 11);
  x = y;
  y = z;
  z = w;
  w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));

  gSeed[tid][0] = x;
  gSeed[tid][1] = y;
  gSeed[tid][2] = z;
  gSeed[tid][3] = w;
  return w * (1.0 / 4294967296.0);
#else
  // @fixme { don't use __thread keyword? }
  static unsigned int THREAD_TLS x = 123456789, y = 362436069, z = 521288629,
                                 w = 88675123;
  unsigned t = x ^ (x << 11);
  x = y;
  y = z;
  z = w;
  w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
  return w * (1.0 / 4294967296.0);
#endif
}


namespace {

#ifdef ENABLE_PTEX
PtexCache*
InitPtex()
{
  PtexCache* c = PtexCache::create(0, kPtexMaxMem);

  return c;
}

PtexTexture*
LoadPtex(
  PtexCache* cache,
  const char* filename)
{
  Ptex::String err;
  PtexTexture* r = PtexTexture::open(filename, err, /* premult */ 0);

  printf("Mallie:info\tmsg:PtexTexture: %p\n", r);

  if (!r) {
    std::cerr << "Mallie:error\tmsg:" << err.c_str() << std::endl;
    return NULL;
  }

  return r;
}

void PtexTest(PtexTexture* r)
{
  // @todo
  PtexFilter::Options opts(PtexFilter::f_bicubic, 0, 1.0);
  PtexPtr<PtexFilter> f(PtexFilter::getFilter(r, opts));

  float result[4];
  int faceid = 0;
  float u=0, v=0, uw=.125, vw=.125;

  for (v = 0; v <= 1; v += .125) {
    for (u = 0; u <= 1; u += .125) {
      f->eval(result, 0, 1, faceid, u, v, uw, 0, 0, vw);
      printf("%8f %8f -> %8f\n", u, v, result[0]);
    }
  }
}
#endif // ENABLE_PTEX


static void GenerateBasis(real3 &tangent, real3 &binormal,
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

// Mis power (1 for balance heuristic)
double Mis(double aPdf) { return aPdf; }

// Mis weight for 2 pdfs
double Mis2(double aSamplePdf, double aOtherPdf) {
  return Mis(aSamplePdf) / (Mis(aSamplePdf) + Mis(aOtherPdf));
}

void GenEyePath(const Scene &scene, int x, int y) {

  double u0 = randomreal();
  double u1 = randomreal();
}

void GenLightPath(Scene &scene, int numPhotons) {
  std::vector<Path> paths;

  real3 lightPos = real3(0.0, 20.0, 0.0);
  real3 lightDir = real3(0.0, -1.0, 0.0);

  for (int i = 0; i < numPhotons; i++) {
    Path path;

    Ray ray;

    real3 dir = lightDir;
    dir.normalize();

    ray.dir = dir;
    ray.org = lightPos;

    Intersection isect;
    bool hit = scene.Trace(isect, ray);

    paths.push_back(path);
  }
}

void EnvCol(float rgba[4], const Scene &scene, const real3& dir)
{
  float d[3];
  d[0] = dir[0];
  d[1] = dir[1];
  d[2] = dir[2];

  if (scene.GetEnvMap().IsValid()) {
    if (scene.GetEnvMap().coordinate() == Texture::COORDINATE_LONGLAT) {
      LongLatMapSampler::Sample(rgba, d, &(scene.GetEnvMap()));
    } else if (scene.GetEnvMap().coordinate() == Texture::COORDINATE_ANGULAR) {
      AngularMapSampler::Sample(rgba, d, &(scene.GetEnvMap()));
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

    rgba[0] = 0.0;
    rgba[1] = 0.0;
    rgba[2] = 0.0;
    rgba[3] = 1.0;
  }
}

} // namespace

// Importance sample diffuse BRDF.
double SampleDiffuseIS(real3 &dir, const real3 &normal) {
  real3 tangent, binormal;

  GenerateBasis(tangent, binormal, normal);

  double theta = acos(sqrt(1.0 - randomreal()));
  double phi = 2.0 * M_PI * randomreal();

  double cosTheta = cos(theta);

  /* D = T*cos(phi)*sin(theta) + B*sin(phi)*sin(theta) + N*cos(theta) */
  double cos_theta = cos(theta);
  real3 T = tangent * cos(phi) * sin(theta);
  real3 B = binormal * sin(phi) * sin(theta);
  real3 N = normal * (cos_theta);

  dir = T + B + N;

  return cos_theta; // PDF = weight
}


bool TraceRay(Intersection& isect, const Scene &scene, Ray &ray) {

  isect.t = std::numeric_limits<real>::max(); // far

  int hit = 0;
  hit = (int)scene.Trace(isect, ray);

  // plane hit test
  if (0) {
    Intersection planeIsect; planeIsect.t = std::numeric_limits<real>::max();
    bool planeHit = gPlane.Intersect(&planeIsect, ray);
    if (planeHit && (planeIsect.t < isect.t)) {
      isect = planeIsect;
    } 

    hit |= (int)planeHit;
  }

  return (hit ? true : false);

}


#if 0
real3 PathTrace(const Scene &scene, const Camera &camera,
                int px, int py) {
  //
  // 1. Sample eye(E0)
  //
  float u = randomreal() - 0.5;
  float v = randomreal() - 0.5;

  // Ray ray = camera.GenerateRay(px + u + step / 2.0f, py + v + step / 2.0f);
  Ray ray = camera.GenerateRay(px + u, py + v);

  Intersection isect;
  isect.t = kFar;

  real3 throughput;
  real3 radiance = real3(0.0, 0.0, 0.0);
  unsigned int pathLength = 1;
  bool lastSpecular = true;
  double lastPdfW = 1.0;

  real3 kLight = real3(0,-1,0);

  for (;; ++pathLength) {
    bool hit = scene.Trace(isect, ray);
    if (!hit) {

      if (pathLength < kMinPathLength) {
        // eye -> background hit.
        break;
      }

      real3 nf = isect.normal;
      if(vdot(nf,ray.dir)>0)nf=nf.neg();
      real kl = vdot(kLight, nf)+0.5;
      real3 nf2 = (nf+real3(1,1,1))*0.5;
      radiance = real3(kl, kl, kl)*nf2;
      return radiance;

      // Hit background.
#if 0
      real3 kd = real3(0.5, 0.5, 0.5);
#else
      real is = std::max(0.05, ray.dir.y);
      real3 kd = real3(is, is, is);
#endif
      radiance += kd / real3(pathLength, pathLength, pathLength);
    }

    // ptex
    //return real3(isect.u, isect.v, isect.faceID*0.001);

    if (pathLength >= kMaxPathLength) {
      break;
    }

    real3 hitP = ray.org + isect.t * ray.dir;

    // 2. Next event estimation
    {}

    // 3. Continue path tracing.
    {
      double r = randomreal();
      real3 sampledDir;

      // faceforward.
      real3 n = isect.normal;
      double ndoti = vdot(isect.normal, ray.dir.neg());
      if (ndoti < 0.0) {
        n = n.neg();
      }
      double pdf = SampleDiffuseIS(sampledDir, n);

      // throughput *= factor * (cosThetaOut / pdf);

      ray.org = hitP + kEPS * sampledDir;
      ray.dir = sampledDir;

      isect.t = kFar;
    }
  }

  return radiance;
}
#endif
real3 PathTrace(const Scene &scene, const Camera &camera,
                const Intersection& s, const Ray& inRay) {

  real3 throughput;
  real3 radiance = real3(0.0, 0.0, 0.0);
  unsigned int pathLength = 1;
  bool lastSpecular = true;
  double lastPdfW = 1.0;

  Intersection isect;
  Ray ray;

  {
    isect.t = kFar;

    double r = randomreal();
    real3 sampledDir;

    // faceforward.
    real3 n = s.normal;
    double ndoti = vdot(s.normal, inRay.dir.neg());
    if (ndoti < 0.0) {
      n = n.neg();
    }
    double pdf = SampleDiffuseIS(sampledDir, n);

    ray.org = s.position + kEPS * sampledDir;
    ray.dir = sampledDir;
  }

  real3 kLight = real3(0,1,0);

  for (;; ++pathLength) {
    bool hit = scene.Trace(isect, ray);
    if (!hit) {

      if (pathLength < kMinPathLength) {
        // eye -> background hit.
        break;
      }

      //real3 nf = isect.normal;
      //if(vdot(nf,ray.dir)>0)nf=nf.neg();
      //real kl = vdot(kLight, nf)+0.5;
      //real3 nf2 = (nf+real3(1,1,1))*0.5;
      //radiance = real3(kl, kl, kl)*nf2;
      //return radiance;

      // Hit background.
#if 1
      real3 kd = real3(0.5, 0.5, 0.5);
#else
      real is = std::max(0.05, ray.dir.y);
      real3 kd = real3(is, is, is);
#endif
      radiance += kd / real3(pathLength, pathLength, pathLength);
    }

    // ptex
    //return real3(isect.u, isect.v, isect.faceID*0.001);

    if (pathLength >= kMaxPathLength) {
      break;
    }

    real3 hitP = ray.org + isect.t * ray.dir;

    // 2. Next event estimation
    {}

    // 3. Continue path tracing.
    {
      double r = randomreal();
      real3 sampledDir;

      // faceforward.
      real3 n = isect.normal;
      double ndoti = vdot(isect.normal, ray.dir.neg());
      if (ndoti < 0.0) {
        n = n.neg();
      }
      double pdf = SampleDiffuseIS(sampledDir, n);

      // throughput *= factor * (cosThetaOut / pdf);

      ray.org = hitP + kEPS * sampledDir;
      ray.dir = sampledDir;

      isect.t = kFar;
    }
  }

  return radiance;
}

void Render(Scene &scene, const RenderConfig &config,
            std::vector<float> &image, // RGB
            std::vector<int> &count, const double eye[3],
            const double lookat[3], const double up[3], const double quat[4],
            int step) {
  int width = config.width;
  int height = config.height;
  double fov = config.fov;

  // For Posprocess.
  std::vector<int>   geomIDBuffer(width*height);
  memset(&geomIDBuffer.at(0), -1, width*height*sizeof(int)); // -1 = background
  std::vector<float> normalBuffer(width*height*3);
  std::vector<float> depthBuffer(width*height);

  double origin[3], corner[3], du[3], dv[3];
  Camera camera(eye, lookat, up);
  camera.BuildCameraFrame(origin, corner, du, dv, fov, quat, width, height);
  // printf("[Mallie] origin = %f, %f, %f\n", gOrigin[0], gOrigin[1],
  // gOrigin[2]);
  // printf("[Mallie] corner = %f, %f, %f\n", gCorner[0], gCorner[1],
  // gCorner[2]);
  // printf("[Mallie] du     = %f, %f, %f\n", gDu[0], gDu[1], gDu[2]);
  // printf("[Mallie] dv     = %f, %f, %f\n", gDv[0], gDv[1], gDv[2]);

  assert(image.size() >= 3 * width * height);
  // memset(&image.at(0), 0, sizeof(float) * width * height * 3);

  init_randomreal();

  if (!gPlaneInitialied) {
    real3 sceneBMin, sceneBMax;
    scene.BoundingBox(sceneBMin, sceneBMax);

    Material planeMat;
    planeMat.diffuse[0] = 0.0;
    planeMat.diffuse[1] = 0.0;
    planeMat.diffuse[2] = 0.0;
    planeMat.reflection[0] = 0.8;
    planeMat.reflection[1] = 0.8;
    planeMat.reflection[2] = 0.8;
    planeMat.fresnel = true;
    planeMat.ior = 1.33;
    planeMat.reflection_glossiness = 0.99;
    size_t planeMatID = scene.AddMaterial(planeMat);

    init_plane(sceneBMin, sceneBMax, planeMatID);

  }

  mallie::timerutil t;
  mallie::timerutil tEventTimer;

  ClearRenderCancel();

  t.start();
  tEventTimer.start();

  //
  // Clear background with gradation.
  //
  memset(&image[0], 0, sizeof(float) * width * height * 3);

  int y_count = 0;
  int y_count_iter = 16;
  bool canceled = false;

  for (int y = 0; y < height; y += step) {

    {
      y_count++;
      if (y_count > y_count_iter) {
        canceled = CheckRenderCancel();
        y_count = 0;
      }
    } 

    if (canceled) {
      //printf("cancel\n");
      break;
    }

    //if ((y % 100) == 0) {
      //printf("\rMallie:info\tRender %d of %d", y, height);
      //fflush(stdout);
    //}

#if 0
    for (int x = 0; x < width; x += step) {

      // random sample pixel position in [step x step] sized tile.
      int px = x + (int)(randomreal() * step);
      int py = y + (int)(randomreal() * step);
      px = std::min(px, (width - 1));
      py = std::min(py, (height - 1));

      real3 radiance =
          PathTrace(scene, camera, config, image, count, px, py, 1);

      image[3 * (py * width + px) + 0] = radiance[0];
      image[3 * (py * width + px) + 1] = radiance[1];
      image[3 * (py * width + px) + 2] = radiance[2];
      count[py * width + px]++;
    }

#else

    #pragma omp parallel for schedule(dynamic, 1)
    for (int x = 0; x < width; x += step) {

      float u = randomreal() - 0.5;
      float v = randomreal() - 0.5;

      Ray ray = camera.GenerateRay(x + u + step / 2.0f, y + v + step / 2.0f);
      ray.depth = 0;

      Intersection isect;
      //bool hit = scene.Trace(isect, ray);
      bool hit = TraceRay(isect, scene, ray);

      if (hit) {

#if 0
        double dotNI = fabs(vdot(isect.normal, ray.dir.neg()));

        image[3 * (y * width + x) + 0] = isect.normal[0];
        image[3 * (y * width + x) + 1] = isect.normal[1];
        image[3 * (y * width + x) + 2] = isect.normal[2];
#else
        float rgba[4];
        scene.GetShader().Shade(rgba, scene, isect, ray);
        image[3 * (y * width + x) + 0] = rgba[0];
        image[3 * (y * width + x) + 1] = rgba[1];
        image[3 * (y * width + x) + 2] = rgba[2];

        geomIDBuffer[(y*width+x)] = isect.faceID;
        normalBuffer[3*(y*width+x)+0] = isect.normal[0];
        normalBuffer[3*(y*width+x)+1] = isect.normal[1];
        normalBuffer[3*(y*width+x)+2] = isect.normal[2];
        depthBuffer[(y*width+x)] = isect.t;

#endif

        // block fill
        for (int v = 0; v < step; v++) {
          if (y+v >= height) continue;
          for (int u = 0; u < step; u++) {
            if (x+u >= width) continue;
            for (int k = 0; k < 3; k++) {
              image[((y + v) * width * 3 + (x + u) * 3) + k] =
                  image[3 * (y * width + x) + k];
              count[(y+v) * width + (x+u)]++;
            }

            geomIDBuffer[((y+v)*width+(x+u))] = isect.faceID;
            normalBuffer[3*((y+v)*width+(x+u))+0] = isect.normal[0];
            normalBuffer[3*((y+v)*width+(x+u))+1] = isect.normal[1];
            normalBuffer[3*((y+v)*width+(x+u))+2] = isect.normal[2];
            depthBuffer[((y+v)*width+(x+u))] = isect.t;
          }
        }
      } else {
        // fetch envcolor.
        float rgba[4];
        EnvCol(rgba, scene, ray.dir);
        image[3 * (y * width + x) + 0] = rgba[0];
        image[3 * (y * width + x) + 1] = rgba[1];
        image[3 * (y * width + x) + 2] = rgba[2];
        
        // block fill
        for (int v = 0; v < step; v++) {
          if (y+v >= height) continue;
          for (int u = 0; u < step; u++) {
            if (x+u >= width) continue;
            for (int k = 0; k < 3; k++) {
              image[((y + v) * width * 3 + (x + u) * 3) + k] =
                  image[3 * (y * width + x) + k];
              count[(y+v) * width + (x+u)]++;
            }
          }
        }
      }
    }

#endif
  }


  // Postprocess
  if (config.wireframe) {
    FeatureLineImageSpace featureLineDrawer;

    // filter param
    float featureLineWidth = 1.0;
    int N = 2;
    int h = featureLineWidth;

    std::vector<float> img(width * height);
    featureLineDrawer.Filter(&img.at(0), &geomIDBuffer.at(0), &normalBuffer.at(0), &depthBuffer.at(0), width, height, N, h);

    // overlay
    for (int i = 0; i < width * height; i++) {
      image[3 * i + 0] += img[i];
      image[3 * i + 1] += img[i];
      image[3 * i + 2] += img[i];
    }
  }

  t.end();

  double fps = 1000.0 / (double)t.msec();
  printf("\r[Mallie] Render time: %f sec(s) | %f fps",
         (double)t.msec() / 1000.0, fps);
  fflush(stdout);
}

} // namespace
