#include "feature-lines.h"

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <cstring>
#include <cstdio>

#include <vector>

#include "common.h"

using namespace mallie;

typedef struct {
  float x;
  float y;
} Point;

static inline int clamp(float x, int w) {
  int i = (int)x;
  if (i < 0)
    i = 0.;
  if (i > (w - 1))
    i = w - 1;
  return i;
}

FeatureLineImageSpace::FeatureLineImageSpace() {}

FeatureLineImageSpace::~FeatureLineImageSpace() {}

void FeatureLineImageSpace::Filter(float *outputs, // output image.
                                   int *geomIDs,   // int
                                   float *normals, // [xyz] x float
                                   float *depths,  // float
                                   int width,      // image width,
                                   int height,     // image height,
                                   int N, float h) {
  memset(outputs, 0, sizeof(float) * width * height);

  // Stencil point to compute gradient
  float stencilDir[4][2];
  stencilDir[0][0] = h;
  stencilDir[0][1] = 0.0f;
  stencilDir[1][0] = 0.0f;
  stencilDir[1][1] = h;
  stencilDir[2][0] = -h;
  stencilDir[2][1] = 0.0f;
  stencilDir[3][0] = 0.0f;
  stencilDir[3][1] = -h;

  // Precompute sampling point.
  std::vector<Point> samplePoints;
  for (int i = 1; i <= N; i++) {
    // Concentric circles.
    float radius = (i / (float)N);
    for (int k = 0; k < 8; k++) {
      float theta = 2.0 * M_PI * (k / 8.0);
      Point p;
      p.x = radius * cos(theta);
      p.y = radius * sin(theta);
      samplePoints.push_back(p);
    }
  }
  float mHalf = samplePoints.size() * 0.5f;

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {

      // s = center of sampling point.
      Point s;
      s.x = x + 0.5f;
      s.y = y + 0.5f;

      // 1) compute m.
      //    m is the number of samples(rays) r in filter region for which g_r !=
      //    g_s
      //    (i.e., the number of samples hitting different geom ID)
      int m = 0;

      int g_s = geomIDs[y * width + x];
      float normalDiffMax = 1.0f;
      float distDiffMax = 0.0f;
      float distDiffThres = 1.5f; // @fixme. scene scale dependent
      float t_s = depths[y * width + x];
      real3 normal;
      normal[0] = normals[3 * (y * width + x) + 0];
      normal[1] = normals[3 * (y * width + x) + 1];
      normal[2] = normals[3 * (y * width + x) + 2];
      for (int k = 0; k < samplePoints.size(); k++) {
        int px = clamp(s.x + samplePoints[k].x, width);
        int py = clamp(s.y + samplePoints[k].y, height);
        int g_r = geomIDs[py * width + px];
        float t_r = depths[py * width + px];
        if (g_s != g_r) {
          m++;
        }
        real3 n;
        n[0] = normals[3 * (py * width + px) + 0];
        n[1] = normals[3 * (py * width + px) + 1];
        n[2] = normals[3 * (py * width + px) + 2];
        // same dir = near to 1.0, diff dir = near to 0.0, so we take min() to
        // compute diffMax.
        float normalDiff = fabs(vdot(n, normal));
        normalDiffMax = std::min(normalDiff, normalDiffMax);

        if (fabsf(t_s - t_r) > 0.1f && (t_s > 0.1f) && (t_r > 0.1f)) {
          distDiffMax = std::max(distDiffMax, fabsf(t_s - t_r));
        }
      }

      outputs[y * width + x] = 0.0f;

      if (m == 0) {

        // Geom ID are all same.

        // 2.2) Check crease-edge.
        //      If any normal is significantly different, it is crease-edge.
        // if (normalDiffMax < 0.25) {
        //  outputs[y*width+x] = 4.0f * (0.25 - normalDiffMax);
        //}

        // 2.2) Check self-occlusion.
        // Disabled for now because visual quality is not good.
        // if (distDiffMax > distDiffThres) {
        //  outputs[y*width+x] = 1.0f;
        //}

      } else {
        float e = 1.0 - fabs(m - mHalf) / mHalf;
        outputs[y * width + x] = e;
      }
    }
  }
}
