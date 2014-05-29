#include <cmath>
#include <limits>
#include <iostream>
#include "texture.h"
#include "vector3.h"

#ifdef ENABLE_PTEX
#include <Ptexture.h>
#endif

#ifdef ENABLE_PTEX
PtexCache *InitPtex() {
  int maxMem = 1024 * 1024;
  PtexCache *c = PtexCache::create(0, maxMem);

  return c;
}

PtexTexture *LoadPtex(PtexCache *cache, const char *filename) {
  Ptex::String err;
  PtexTexture *r = PtexTexture::open(filename, err, /* premult */ 0);

  printf("[Mallie] PtexTexture: %p\n", r);

  if (!r) {
    std::cerr << err.c_str() << std::endl;
    return NULL;
  }

  return r;
}
#endif

namespace {

#ifndef M_PI
#define M_PI 3.141592
#endif

int inline fasterfloor( const float x ) {
  if (x >= 0) {
    return (int)x;
  }

  int y = (int)x;
  if (std::abs(x - y) <= std::numeric_limits<float>::epsilon()) {
    // Do nothing.
  } else {
    y = y - 1;
  }

  return y;
}


bool myisnan(float a)
{
    volatile float d = a;
    return d != d;
}

inline void FilterByte(
  float* rgba,
  const unsigned char* image,
  int i00, int i10, int i01, int i11,
  float w[4], // weight
  int stride)
{
  unsigned char texel[4][4];

  const float inv = 1.0f / 255.0f;
  if (stride == 4) {

    for (int i = 0; i < 4; i++) {
      texel[0][i] = image[i00+i];
      texel[1][i] = image[i10+i];
      texel[2][i] = image[i01+i];
      texel[3][i] = image[i11+i];
    }

    for (int i = 0; i < 4; i++) {
      rgba[i] = w[0] * texel[0][i] +
                w[1] * texel[1][i] +
                w[2] * texel[2][i] +
                w[3] * texel[3][i];
      // normalize.
      rgba[i] *= inv;
    }

  } else {

    for (int i = 0; i < stride; i++) {
        texel[0][i] = image[i00+i];
        texel[1][i] = image[i10+i];
        texel[2][i] = image[i01+i];
        texel[3][i] = image[i11+i];
    }

    for (int i = 0; i < stride; i++) {
      rgba[i] = w[0] * texel[0][i] +
                w[1] * texel[1][i] +
                w[2] * texel[2][i] +
                w[3] * texel[3][i];
      // normalize.
      rgba[i] *= inv;
    }
  }

  if (stride < 4) {
    rgba[3] = 0.0;
  }

}

inline void FilterFloat(
  float* rgba,
  const float* image,
  int i00, int i10, int i01, int i11,
  float w[4], // weight
  int stride)
{
  float texel[4][4];

  if (stride == 4) {

    for (int i = 0; i < 4; i++) {
      texel[0][i] = image[i00+i];
      texel[1][i] = image[i10+i];
      texel[2][i] = image[i01+i];
      texel[3][i] = image[i11+i];
    }

    for (int i = 0; i < 4; i++) {
      rgba[i] = w[0] * texel[0][i] +
                w[1] * texel[1][i] +
                w[2] * texel[2][i] +
                w[3] * texel[3][i];
    }

  } else {

    for (int i = 0; i < stride; i++) {
        texel[0][i] = image[i00+i];
        texel[1][i] = image[i10+i];
        texel[2][i] = image[i01+i];
        texel[3][i] = image[i11+i];
    }

    for (int i = 0; i < stride; i++) {
      rgba[i] = w[0] * texel[0][i] +
                w[1] * texel[1][i] +
                w[2] * texel[2][i] +
                w[3] * texel[3][i];
    }
  }

  if (stride < 4) {
    rgba[3] = 0.0;
  }

}

} // namespace

using namespace mallie;

void
Texture::fetch(
    float *rgba,
    float u, float v) const
{
    if (!IsValid()) {
      if (rgba) {
        rgba[0] = 0.0f;
        rgba[1] = 0.0f;
        rgba[2] = 0.0f;
        rgba[3] = 0.0f;
      }
      return;
    }

    float sx = fasterfloor(u);
    float sy = fasterfloor(v);

    float uu = u - sx;
    float vv = v - sy;

    // clamp
    uu = std::max(uu, 0.0f); uu = std::min(uu, 1.0f);
    vv = std::max(vv, 0.0f); vv = std::min(vv, 1.0f);

    float px = (m_width  - 1) * uu;
    float py = (m_height - 1) * vv;

    int x0 = (int)px;
    int y0 = (int)py;
    int x1 = ((x0 + 1) >= m_width ) ? (m_width  - 1) : (x0 + 1);
    int y1 = ((y0 + 1) >= m_height) ? (m_height - 1) : (y0 + 1);

    float dx = px - (float)x0;
    float dy = py - (float)y0;

    float w[4];

    w[0] = (1.0f - dx) * (1.0 - dy);
    w[1] = (1.0f - dx) * (      dy);
    w[2] = (       dx) * (1.0 - dy);
    w[3] = (       dx) * (      dy);

    
    //unsigned char texel[4][4];

    int stride = m_components;

    int i00 = stride * (y0 * m_width + x0);
    int i01 = stride * (y0 * m_width + x1);
    int i10 = stride * (y1 * m_width + x0);
    int i11 = stride * (y1 * m_width + x1);

#if 0
    const float inv = 1.0f / 255.0f;

    if (stride == 4) {

      for (int i = 0; i < 4; i++) {
          texel[0][i] = m_image[i00+i];
          texel[1][i] = m_image[i10+i];
          texel[2][i] = m_image[i01+i];
          texel[3][i] = m_image[i11+i];
      }

      for (int i = 0; i < 4; i++) {
          rgba[i] = w[0] * texel[0][i] +
                    w[1] * texel[1][i] +
                    w[2] * texel[2][i] +
                    w[3] * texel[3][i];
          // normalize.
          rgba[i] *= inv;
      }

    } else {

      for (int i = 0; i < stride; i++) {
          texel[0][i] = m_image[i00+i];
          texel[1][i] = m_image[i10+i];
          texel[2][i] = m_image[i01+i];
          texel[3][i] = m_image[i11+i];
      }

      for (int i = 0; i < stride; i++) {
          rgba[i] = w[0] * texel[0][i] +
                    w[1] * texel[1][i] +
                    w[2] * texel[2][i] +
                    w[3] * texel[3][i];
          // normalize.
          rgba[i] *= inv;
      }
    }

    if (stride < 4) {
        rgba[3] = 0.0;
    }
#else
    if (m_format == FORMAT_BYTE) {
      FilterByte(rgba, m_image, i00, i10, i01, i11, w, stride);
    } else if (m_format == FORMAT_FLOAT) {
      FilterFloat(rgba, reinterpret_cast<float*>(m_image), i00, i10, i01, i11, w, stride);
    } else { // unknown

    }
#endif
}

void
Texture::fetchD(
    float *rgba0,
    float *rgba1,
    float *rgba2,
    float u, float v) const
{
    // @todo { optimize! }

    // fetch (i, j)
    fetch(rgba0, u, v);
    
    // fetch (i+1, j)
    float u1 = u + m_invWidth;
    fetch(rgba1, u1, v);

    // fetch (i, j+1)
    float v1 = v + m_invHeight;
    fetch(rgba2, u, v1);
}

void
AngularMapSampler::Sample(
    float* rgba,
    float  dir[3],
    const  Texture* texture)
{
    //float u, v;

    // (x,y,z) -> (u, v)

    vector3 d(dir[0], dir[1], dir[2]);
    d.normalize();

    float r;
    if (d[2] >= -1.0 && d[2] < 1.0) { // for safety
        r = (1.0 / M_PI) * acosf(d[2]);
    } else {
        r = 0.0;
    }

    float norm2 = d[0] * d[0] + d[1] * d[1];
    if (norm2 > 1.0e-6) {
        r /= sqrtf(norm2);
    }

    float uu = d[0] * r;
    float vv = d[1] * r;

    uu = 0.5 * uu + 0.5;
    vv = 0.5 * vv + 0.5;

    texture->fetch(rgba, uu, vv);
}

