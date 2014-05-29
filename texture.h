#ifndef __MALLIE_TEXTURE_H__
#define __MALLIE_TEXTURE_H__

namespace mallie {

class
Texture
{

public:

  typedef enum {
    FORMAT_BYTE,
    FORMAT_FLOAT,
  } Format;

  Texture() {
    // Make invalid texture
    m_width = -1;
    m_height = -1;
    m_image = NULL;
    m_components = -1;
  }   

  Texture(unsigned char *image, int width, int height, int components, Format format, float gamma = 1.0f) {
    Set(image, width, height, components, format, gamma);
  }

  ~Texture() { }

  void Set(unsigned char *image, int width, int height, int components, Format format, float gamma = 1.0f) {
    m_width         = width;
    m_height        = height;
    m_image         = image;
    m_invWidth      = 1.0f / width;
    m_invHeight     = 1.0f / height;
    m_components    = components;
    m_format        = format;
    m_invGamma      = 1.0f / gamma; // Take a inv for faster computation.
  }

  int width() const {
    return m_width;
  }

  int height() const {
    return m_height;
  }

  int components() const {
    return m_components;
  }

  unsigned char* image() const {
    return m_image;
  } 
  
  // Bilinear textel fetch.
  void fetch(float *rgba, float u, float v) const;

  // Fetch filtered (i, j), (i+1, j), (i, j+1) texel.
  void fetchD(float *rgba0, float *rgba1, float *rgba2, float u, float v) const;

  bool IsValid() const {
    return (m_image != NULL) && (m_width > 0) && (m_height > 0);
  }

private:

  int             m_width;
  int             m_height;
  float           m_invWidth;
  float           m_invHeight;
  int             m_components;
  unsigned char*  m_image;
  float           m_invGamma;
  int             m_format;    
};

//
// A sampler class for angular map texture.
//
class
AngularMapSampler
{
 public:
  AngularMapSampler() {};
  ~AngularMapSampler() {};

  static void Sample(float* rgba, float dir[3], const Texture* texture);
};

// PTex
class PTexture
{
public:
  PTexture();
  ~PTexture();

private:
#ifdef ENABLE_PTEX

#endif
};

} // namespace

#endif // __MALLIE_TEXTURE_H__
