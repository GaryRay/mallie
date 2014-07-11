#ifndef __MALLIE__FEATURE_LINES_H__
#define __MALLIE__FEATURE_LINES_H__

// Image space feature line drawing.

namespace mallie {

class FeatureLineImageSpace {
public:
  FeatureLineImageSpace();
  ~FeatureLineImageSpace();

  /// Draw feature lines in outputs,
  /// Input: geom ID image, normal image, depth image
  ///        filter params: N, h
  /// geomID = -1 -> background
  void Filter(float *outputs, // output image.
              int *geomIDs,   // int
              float *normals, // [xyz] x float
              float *depths,  // float
              int width,      // image width,
              int height,     // image height,
              int N, float h);
};
}

#endif // __MALLIE__FEATURE_LINES_H__
