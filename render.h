#ifndef __MALLIE_RENDER_H__
#define __MALLIE_RENDER_H__

#include <vector>
#include <string>

#include "scene.h"
#include "tinythread.h"

namespace mallie {

extern void ClearRenderCancel();
extern void PostRenderCancel();
extern bool CheckRenderCancel();

struct RenderConfig {
  double fov;
  int width;
  int height;
  double eye[3];
  double lookat[3];
  double up[3];
  double quat[4];

  double scene_scale;

  bool plane; // Infinite plane for debugging purpose

  int num_passes;
  int num_photons; // # of photon to shoot per pass.

  bool wireframe;

  std::string obj_filename;
  std::string eson_filename;
  std::string material_filename;
  std::string envmap_filename;
  std::string ptex_filename;
  std::string envmap_coord; // "longlat" or "angularmap"

  float display_gamma;

  RenderConfig()
      : fov(60.0), width(512), height(512), scene_scale(1.0), plane(false),
        num_passes(10), num_photons(10000), display_gamma(2.2f), wireframe(false) {

    eye[0] = 0.0;
    eye[1] = 0.0;
    eye[2] = -5.0;
    lookat[0] = 0.0;
    lookat[1] = 0.0;
    lookat[2] = 0.0;
    up[0] = 0.0;
    up[1] = 1.0;
    up[2] = 0.0;
    quat[0] = quat[1] = quat[2] = quat[3] = 0.0;

    envmap_coord = "longlat";
  }
};

extern void Render(Scene &scene, const RenderConfig &config,
                   std::vector<float> &image, // out image
                   std::vector<int> &count,   // per-pixel counter
                   const double eye[3], const double lookat[3],
                   const double up[3], const double quat[4], int step);
}

#endif // __MALLIE_RENDER_H__
