#include "main_sdl.h"

#include <cassert>
#include <vector>
#include <cassert>
#include <iostream>


#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef ENABLE_SDL

#include <SDL.h> // SDL2

#include "trackball.h"
#include "camera.h"
#include "timerutil.h"
#include "render.h"
#include "tinythread.h"
#include "script_engine.h"
#include "jpge.h"

#if defined(_WIN32) && !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES
#endif
#include <cmath>

#define PIXELSTEP_COARSE (8)

namespace mallie {

//
// -- GUI variables --
//
static bool gViewChanged = true; // initial is true.
static bool gRenderInteractive = false;
static int gMouseX = -1, gMouseY = -1;
static bool gMouseMoving = false;
static int gMouseButton = 0;
static bool gNeedRedraw = true;
static bool gShiftPressed = false;
static bool gCtrlPressed = false;
static bool gLightEditing = false;

static double gEye[3] = {0.0, 0.0, -5.0};
static double gUp[3] = {0.0, 1.0, 0.0};
static double gLookat[3] = {0.0, 0.0, 0.0};
static double gFov = 45.0;
static double gScale = 0.1;
static double gPrevQuat[4] = {0.0, 0.0, 0.0, 0.0};
static double gCurrQuat[4] = {0.0, 0.0, 0.0, 0.0};
static double gInitQuat[4] = {0.0, 0.0, 0.0, 1.0};

static double gRotate[3] = {0.0f, 0.0f, 0.0f};
static double gOrigin[3], gCorner[3], gDu[3], gDv[3];

static double gIntensity = 1.0;

static int gShaderIndex = 0;

// Progressive render param
static int gRenderPixelStep = 1; //PIXELSTEP_COARSE;
static int gRenderPasses = 1;

static tthread::mutex gRenderThreadMutex;
static time_t gRenderClock = 0;
static bool gRenderQuit = false; // Only become true when we quit app.
static bool gRenderCancel = false;

int gWidth = 256;
int gHeight = 256;
float gGamma = 1.0f;

// SDL_Surface* gSurface = NULL;
SDL_Window *gWindow = NULL;
SDL_Surface *gSurface = NULL;
SDL_Renderer *gSDLRenderer = NULL;
SDL_mutex *gMutex = NULL;

std::vector<float> gImage;
std::vector<int> gCount;
std::vector<float> gFramebuffer; // HDR framebuffer
RenderConfig gRenderConfig;

typedef struct {
  Scene *scene;
  const RenderConfig *config;
} RenderContext;

inline unsigned char fclamp(float x) {
  int i = x * 255.5;
  if (i < 0)
    return 0;
  if (i > 255)
    return 255;
  return (unsigned char)i;
}

inline float gamma_correct(float x, float inv_gamma) {
  return pow(x, inv_gamma);
}


void HDRToLDR(std::vector<unsigned char> &out, const std::vector<float> &in,
              int width, int height) {
  out.resize(width * height * 3);
  assert(in.size() == (width * height * 3));

  float inv_gamma = 1.0 / gGamma;

  // Simple [0, 1] -> [0, 255]
  for (int i = 0; i < width * height * 3; i++) {
    out[i] = fclamp(gamma_correct(in[i], inv_gamma));
  }
}

void SaveAsJPEG(const char *filename)
{

  tthread::lock_guard<tthread::mutex> guard(gRenderThreadMutex);

  int ret = SDL_LockMutex(gMutex);
  assert(ret == 0);

  std::vector<unsigned char> ldr;

  ldr.resize(gWidth * gHeight * 3);

  SDL_LockSurface(gSurface);

  unsigned char *data = (unsigned char *)gSurface->pixels;

  for (int i = 0; i < gWidth * gHeight; i++) {
    ldr[3*i+0] = data[4*i+2];
    ldr[3*i+1] = data[4*i+1];
    ldr[3*i+2] = data[4*i+0];
  }

  SDL_UnlockSurface(gSurface);
  SDL_UnlockMutex(gMutex);

  jpge::params comp_params;
  comp_params.m_quality = 100;
  ret = jpge::compress_image_to_jpeg_file(filename, gWidth, gHeight, 3,
                                               &ldr.at(0), comp_params);
  assert(ret);

  printf("Save %s\n", filename);
}

void RequestRedraw() {
  tthread::lock_guard<tthread::mutex> guard(gRenderThreadMutex);
  gNeedRedraw = true;
  return;
}

static void EulerToQuatRad(double quat[4], double x, double y,
                           double z) // in radian. yaw, pitch, roll
{
  double rx = x;
  double ry = y;
  double rz = z;

  double hx = 0.5 * rx;
  double hy = 0.5 * ry;
  double hz = 0.5 * rz;

  double cosHx = cos(hx);
  double cosHy = cos(hy);
  double cosHz = cos(hz);

  double sinHx = sin(hx);
  double sinHy = sin(hy);
  double sinHz = sin(hz);

  quat[0] = cosHx * cosHy * cosHz + sinHx * sinHy * sinHz;
  quat[1] = sinHx * cosHy * cosHz - cosHx * sinHy * sinHz;
  quat[2] = cosHx * sinHy * cosHz + sinHx * cosHy * sinHz;
  quat[3] = cosHx * cosHy * sinHz - sinHx * sinHy * cosHz;
}

static void EulerToQuatZYX(double quat[4], double x, double y,
                           double z) // in radian. yaw, pitch, roll
{
  double rx = x;
  double ry = y;
  double rz = z;

  double hx = 0.5 * rx;
  double hy = 0.5 * ry;
  double hz = 0.5 * rz;

  double cosHx = cos(hx);
  double cosHy = cos(hy);
  double cosHz = cos(hz);

  double sinHx = sin(hx);
  double sinHy = sin(hy);
  double sinHz = sin(hz);

  quat[0] = cosHx * cosHy * cosHz + sinHx * sinHy * sinHz;
  quat[1] = sinHx * cosHy * cosHz - cosHx * sinHy * sinHz;
  quat[2] = cosHx * sinHy * cosHz + sinHx * cosHy * sinHz;
  quat[3] = cosHx * cosHy * sinHz - sinHx * sinHy * cosHz;
}

static void AccumImage(std::vector<float> &dst, const std::vector<float> &src) {
  assert(dst.size() == src.size());
  for (size_t i = 0; i < src.size(); i++) {
    dst[i] += src[i];
  }
}

static void ClearImage(std::vector<float> &img) {
  for (size_t i = 0; i < img.size(); i++) {
    img[i] = 0.0f;
  }
}

static void ClearCount(std::vector<int> &img) {
  for (size_t i = 0; i < img.size(); i++) {
    img[i] = 0;
  }
}

void SaveCamera(const std::string &filename) {
  FILE *fp = fopen(filename.c_str(), "w");

  fprintf(fp, "%f %f %f\n", gEye[0], gEye[1], gEye[2]);
  fprintf(fp, "%f %f %f\n", gLookat[0], gLookat[1], gLookat[2]);
  fprintf(fp, "%f %f %f %f\n", gCurrQuat[0], gCurrQuat[1], gCurrQuat[2],
          gCurrQuat[3]);

  fclose(fp);

  std::cout << "Mallie:info\tSave camera data to: " << filename << std::endl;

  char buf[1024];
  sprintf(buf, "eye = [%f, %f, %f];\n", gEye[0], gEye[1], gEye[2]);
  printf("\n");
  ScriptEngine::Eval(buf);
  printf("\n");
}

void LoadCamera(const std::string &filename) {
  FILE *fp = fopen(filename.c_str(), "r");
  if (!fp) {
    std::cerr << "Mallie:error\tmsg:camera file not found " << filename << std::endl;
    return;
  }

  fscanf(fp, "%lf %lf %lf\n", &gEye[0], &gEye[1], &gEye[2]);
  fscanf(fp, "%lf %lf %lf\n", &gLookat[0], &gLookat[1], &gLookat[2]);
  fscanf(fp, "%lf %lf %lf %lf\n", &gCurrQuat[0], &gCurrQuat[1], &gCurrQuat[2],
         &gCurrQuat[3]);

  fclose(fp);

  std::cout << "Mallie:info\tLoad camera data" << std::endl;
  RequestRedraw();
  gViewChanged = true;
}

void HandleMouseButton(SDL_Event e) {

  if (e.type == SDL_MOUSEBUTTONUP) {
    gMouseMoving = false;
    gViewChanged = true;
    gMouseButton = 0;
    gRenderInteractive = false;
    gRenderPasses = 1;
    gRenderPixelStep = 1;

    RequestRedraw();
  
  } else if (e.type == SDL_MOUSEBUTTONDOWN) {

    PostRenderCancel();

    gMouseX = e.motion.x;
    gMouseY = e.motion.y;
    gMouseMoving = true;
    gRenderInteractive = true;
    gRenderPasses = 1;
    gRenderPixelStep = PIXELSTEP_COARSE;

    if (SDL_GetMouseState(NULL, NULL) & SDL_BUTTON(1)) {
      gMouseButton = 1; // left
    } else if (SDL_GetMouseState(NULL, NULL) & SDL_BUTTON(2)) {
      gMouseButton = 2; // middle
    } else if (SDL_GetMouseState(NULL, NULL) & SDL_BUTTON(3)) {
      gMouseButton = 3; // right
    }
    trackball(gPrevQuat, 0.0, 0.0, 0.0, 0.0);
  }
}

void HandleMouseMotion(SDL_Event e) {
  float rotScale = 1.0;

  if (gMouseMoving) {

    int x = e.motion.x;
    int y = e.motion.y;

    gViewChanged = true;
    gRenderInteractive = true;
    gRenderPasses = 1;
    RequestRedraw();
    gRenderPixelStep = PIXELSTEP_COARSE;

    if (gCtrlPressed || (gMouseButton == 3)) {

      gEye[2] -= gScale * (gMouseY - y);
      gLookat[2] -= gScale * (gMouseY - y);

    } else if (gShiftPressed || (gMouseButton == 2)) {

      //printf("shift = %f\n", gScale * (gMouseX - x));
      gEye[0] += gScale * (gMouseX - x);
      gLookat[0] += gScale * (gMouseX - x);

      gEye[1] -= gScale * (gMouseY - y);
      gLookat[1] -= gScale * (gMouseY - y);

    } else {
      // trackball(gPrevQuat,
      //    0.0 * (2.0f * gMouseX - gWidth) / (float)gWidth,
      //    rotScale * (gHeight - 2.0f * gMouseY) / (float)gHeight,
      //    0.0 * (2.0f * x - gWidth) / (float)gWidth,
      //    rotScale * (gHeight - 2.0f * y) / (float)gHeight);

      // trackball(gPrevQuat,
      //    rotScale * (2.0f * gMouseX - gWidth) / (float)gWidth,
      //    rotScale * (gHeight - 2.0f * gMouseY) / (float)gHeight,
      //    rotScale * (2.0f * x - gWidth) / (float)gWidth,
      //    rotScale * (gHeight - 2.0f * y) / (float)gHeight);
      // trackball(gPrevQuat,
      //    0.0f,
      //    rotScale * (gHeight - 2.0f * gMouseY) / (float)gHeight,
      //    0.0f,
      //    rotScale * (gHeight - 2.0f * y) / (float)gHeight);

      double xx = (x - gMouseX) / (double)gWidth;
      double yy = (y - gMouseY) / (double)gHeight;
      double zz = 0.0;
      // EulerToQuatRad(gPrevQuat, xx, yy, zz);
      // printf("quat = %f, %f, %f, %f\n", gPrevQuat[0], gPrevQuat[1],
      // gPrevQuat[2], gPrevQuat[3]);

      double scale = M_PI * 2.0; // Heuristic value
      gRotate[0] += scale * xx;
      gRotate[1] += scale * yy;
      // clamp
      double eps = 1.0e-3;
      if (gRotate[1] <= -(0.5 * M_PI - eps))
        gRotate[1] = -0.5 * M_PI + eps;
      if (gRotate[1] >= (0.5 * M_PI - eps))
        gRotate[1] = 0.5 * M_PI - eps;

      add_quats(gPrevQuat, gCurrQuat, gCurrQuat);
    }
  }

  gMouseX = e.motion.x;
  gMouseY = e.motion.y;
}

bool HandleKey(Scene& scene, SDL_Event e) {
  if (e.type == SDL_KEYUP) {
    gShiftPressed = false;
    gCtrlPressed = false;
    gRenderInteractive = false;
    gRenderPasses = 1;
  } else if (e.type == SDL_KEYDOWN) {
    //gRenderInteractive = true;
    //gRenderPasses = 1;
    //gRenderPixelStep = PIXELSTEP_COARSE;
    switch (e.key.keysym.sym) {
    case SDLK_ESCAPE:
    case 'q':
      // exit(-1);
      return true;
      break;
    case 'd':
      // exit(-1);
      ScriptEngine::Eval("dumpCamera(eye, lookat, up);\n");
      break;
    case SDLK_SPACE:
      // reset rotation
      gEye[0] = gRenderConfig.eye[0];
      gEye[1] = gRenderConfig.eye[1];
      gEye[2] = gRenderConfig.eye[2];
      trackball(gCurrQuat, 0.0f, 0.0f, 0.0f, 0.0f);
      trackball(gPrevQuat, 0.0f, 0.0f, 0.0f, 0.0f);
      gRotate[0] = gRotate[1] = gRotate[2] = 0.0f;
      RequestRedraw();
      break;
    case 'i':
      gIntensity += 0.1f;
      RequestRedraw();
      break;
    case 'o':
      gIntensity -= 0.1f;
      if (gIntensity < 0.1f) {
        gIntensity = 0.1f;
      }
      RequestRedraw();
      break;
    case 'j':
      SaveAsJPEG("output.jpg");
      break;

    case SDLK_LSHIFT:
      gShiftPressed = true;
      break;
    case SDLK_TAB:
    case SDLK_LCTRL:
      gCtrlPressed = true;
      break;
    case 'c':
      SaveCamera("camera.dat");
      break;
    case 'v':
      LoadCamera("camera.dat");
      break;
    case 'm':
      {
        int numShaders = scene.GetShader().NumShaders();
        gShaderIndex++;
      
        scene.GetShader().SetShader(gShaderIndex % numShaders);
        printf("shader = %d\n", gShaderIndex % numShaders);
        RequestRedraw();
      }
      break;
    default:
      break;
    }
  }

  return false;
}

void Display(SDL_Surface *surface, const std::vector<float> &image,
             const std::vector<int> &counts, int width,
             int height, float gamma) {

  // Write to backbuffer.
  tthread::lock_guard<tthread::mutex> guard(gRenderThreadMutex);

  int ret = SDL_LockMutex(gMutex);
  assert(ret == 0);

  SDL_SetRenderDrawColor(gSDLRenderer, 0, 0, 0, 255);
  SDL_RenderClear(gSDLRenderer);

  SDL_LockSurface(surface);

  // ARGB
  unsigned char *data = (unsigned char *)surface->pixels;

  float inv_gamma = 1.0f / gamma;

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {

      // per-pixel count
      int c = counts[y*width+x];
      if (c < 1) c = 1;
      float scale = 1.0f / (float)c;

      //unsigned char col[3];
      //col[0] = x % 255;
      //col[1] = y % 255;
      //col[2] = 127;

#ifdef __APPLE__
      // BGRA?
      data[4 * (y * width + x) + 0] =
          fclamp(gamma_correct(scale * image[3 * (y * width + x) + 2], inv_gamma));
      data[4 * (y * width + x) + 1] =
          fclamp(gamma_correct(scale * image[3 * (y * width + x) + 1], inv_gamma));
      data[4 * (y * width + x) + 2] =
          fclamp(gamma_correct(scale * image[3 * (y * width + x) + 0], inv_gamma));
      data[4 * (y * width + x) + 3] = 255;
#else
      // BGRA?
      data[4 * (y * width + x) + 2] =
          fclamp(gamma_correct(scale * image[3 * (y * width + x) + 0], inv_gamma));
      data[4 * (y * width + x) + 1] =
          fclamp(gamma_correct(scale * image[3 * (y * width + x) + 1], inv_gamma));
      data[4 * (y * width + x) + 0] =
          fclamp(gamma_correct(scale * image[3 * (y * width + x) + 2], inv_gamma));
      data[4 * (y * width + x) + 3] = 255;
#endif
    }
  }

  SDL_UnlockSurface(surface);
  // SDL_RenderPresent(gSDLRenderer);

  SDL_UnlockMutex(gMutex);
}

Uint32 TimeLeft(int interval) {
  static Uint32 next_time = 0;
  Uint32 now;

  now = SDL_GetTicks();
  if (next_time <= now) {
    next_time = now + interval;
    return (0);
  }

  return (next_time - now);
}

static bool CheckSDLEvent() {
  SDL_Event event;
  SDL_PumpEvents();
  // if (SDL_PeepEvents(&event, 1, SDL_PEEKEVENT, SDL_EVENTMASK
  //(SDL_MOUSEBUTTONDOWN) | SDL_EVENTMASK(SDL_KEYDOWN)) > 0) {
  if (SDL_PeepEvents(&event, 1, SDL_PEEKEVENT,
                     SDL_MOUSEBUTTONDOWN | SDL_KEYDOWN,
                     SDL_MOUSEBUTTONDOWN | SDL_KEYDOWN) > 0) {
    return true;
  }
  return false;
}

static void Init(const RenderConfig &config) {
  // Save
  gRenderConfig = config;

  gFov = config.fov;

  gEye[0] = config.eye[0];
  gEye[1] = config.eye[1];
  gEye[2] = config.eye[2];

  gLookat[0] = config.lookat[0];
  gLookat[1] = config.lookat[1];
  gLookat[2] = config.lookat[2];

  gUp[0] = config.up[0];
  gUp[1] = config.up[1];
  gUp[2] = config.up[2];

  trackball(gCurrQuat, 0.0f, 0.0f, 0.0f, 0.0f);
  add_quats(gInitQuat, gCurrQuat, gCurrQuat);

  Camera camera(gEye, gLookat, gUp);
  camera.BuildCameraFrame(gOrigin, gCorner, gDu, gDv, gFov, gCurrQuat, gWidth,
                          gHeight);
  // printf("[Mallie] eye    = %f, %f, %f\n", gEye[0], gEye[1], gEye[2]);
  // printf("[Mallie] lookat = %f, %f, %f\n", gLookat[0], gLookat[1],
  // gLookat[2]);
  // printf("[Mallie] up     = %f, %f, %f\n", gUp[0], gUp[1], gUp[2]);

  gGamma = config.display_gamma;
}

time_t GetCurrentRenderClock() {
  time_t clk;
  tthread::lock_guard<tthread::mutex> guard(gRenderThreadMutex);
  clk = gRenderClock;

  return clk;
}

void NotifyRenderClock() {
  tthread::lock_guard<tthread::mutex> guard(gRenderThreadMutex);
  timerutil t;
  gRenderClock = t.current();
  return;
}

bool GetRenderQuitRequest() {
  bool ret;
  tthread::lock_guard<tthread::mutex> guard(gRenderThreadMutex);
  ret = gRenderQuit;

  return ret;
}

void NotifyRenderQuit() {
  tthread::lock_guard<tthread::mutex> guard(gRenderThreadMutex);
  gRenderQuit = true;
  return;
}

void RenderThread(void *arg) {
  RenderContext ctx = *(reinterpret_cast<RenderContext *>(arg));

  time_t prevRenderClock = 0;

  while (!GetRenderQuitRequest()) {

    time_t currentRenderClock = GetCurrentRenderClock();
    if ((gRenderPasses >= ctx.config->num_passes) &&
        (currentRenderClock <= prevRenderClock)) {
#ifdef _WIN32
      Sleep(33);
#else
      usleep(1000 * 33);
#endif
      continue;
    }

    if ((gRenderPasses >= ctx.config->num_passes)) {
      // printf("Render finished\n");
      // render finished
      // Display(gSurface, gFramebuffer, gRenderPasses, config.width,
      // config.height);
      continue;
    }

    // redraw request check.
    {
      tthread::lock_guard<tthread::mutex> guard(gRenderThreadMutex);
      if (gNeedRedraw) {
        ClearImage(gFramebuffer);
        ClearCount(gCount);

        gNeedRedraw = false;
      }
    }

    // Use Euler rotation.
    // printf("rot = %f, %f, %f\n", 180*gRotate[0]/M_PI, 180*gRotate[1]/M_PI,
    // 180*gRotate[2]/M_PI);
    EulerToQuatRad(gCurrQuat, gRotate[2], gRotate[0], gRotate[1] + M_PI);
    // printf("quat = %f, %f, %f, %f\n", gCurrQuat[0], gCurrQuat[1],
    // gCurrQuat[2], gCurrQuat[3]);

    Render(*(ctx.scene), *(ctx.config), gImage, gCount, gEye, gLookat, gUp,
           gCurrQuat, gRenderPixelStep);

    // redraw request check.
    {
      tthread::lock_guard<tthread::mutex> guard(gRenderThreadMutex);
      if (gNeedRedraw) {
        // Discard render result and do refresh rendering.
        ClearImage(gFramebuffer);
        ClearCount(gCount);

        gNeedRedraw = false;
        continue;
      }
    }


    AccumImage(gFramebuffer, gImage);

    Display(gSurface, gFramebuffer, gCount, ctx.config->width,
            ctx.config->height, ctx.config->display_gamma);

    //if (gMouseMoving) {
    //  ClearImage(gFramebuffer);
    //  ClearCount(gCount);
    //}

// printf("step = %d, interactive = %d\n", gRenderPixelStep,
// gRenderInteractive);

#if 0
    // Increment render pass.
    if (!gRenderInteractive && gRenderPixelStep == 1) {
      gRenderPasses++;
    }

    if (!gRenderInteractive) {
      gRenderPixelStep >>= 1;

      if (gRenderPixelStep == 1) {
        ClearImage(gFramebuffer);
      }

      if (gRenderPixelStep < 1) {
        gRenderPixelStep = 1;
      }
    } else {
      gRenderPasses = 1;
    }
#else
    gRenderPasses++;
#endif

    prevRenderClock = currentRenderClock;
  }
}

void DoMainSDL(Scene &scene, const RenderConfig &config) {
  printf("Mallie:info\tSDL window mode.\n");

  gWidth = config.width;
  gHeight = config.height;

  gWindow = SDL_CreateWindow("Mallie", SDL_WINDOWPOS_UNDEFINED,
                             SDL_WINDOWPOS_UNDEFINED, gWidth, gHeight, 0);
  if (!gWindow) {
    printf("Mallie:error\tSDL err: %s\n", SDL_GetError());
    exit(1);
  }

  gSDLRenderer = SDL_CreateRenderer(gWindow, -1, SDL_RENDERER_SOFTWARE);
  if (!gSDLRenderer) {
    printf("Mallie:error\tSDL err: %s\n", SDL_GetError());
    exit(1);
  }

  gSurface = SDL_GetWindowSurface(gWindow);
  if (!gSurface) {
    printf("Mallie:error\tSDL err: %s\n", SDL_GetError());
    exit(1);
  }

  gMutex = SDL_CreateMutex();

  gFramebuffer.resize(gWidth * gHeight * 3); // RGB
  ClearImage(gFramebuffer);

  gImage.resize(gWidth * gHeight * 3); // RGB

  gCount.resize(gWidth * gHeight);
  ClearCount(gCount);

  Init(config);

  RenderContext renderCtx;
  renderCtx.scene = &scene;
  renderCtx.config = &config;

  tthread::thread renderThread(RenderThread, (void *)&renderCtx);

  SDL_Event event;

  bool done = false;
  while (!done) {

    bool hasEvent = false;

    while (SDL_PollEvent(&event)) {
      switch (event.type) {
      case SDL_QUIT:
        done = true;
        break;
      case SDL_KEYUP:
        hasEvent = true;
      case SDL_KEYDOWN:
        done = HandleKey(scene, event);
        break;
      case SDL_MOUSEBUTTONUP:
        hasEvent = true;
      case SDL_MOUSEBUTTONDOWN:
        HandleMouseButton(event);
        break;
      case SDL_MOUSEMOTION:
        HandleMouseMotion(event);
        break;
      }

      if (done) {
        break;
      }
    }

    if (done) {
      break;
    }

#if 0
    if ((gRenderPasses >= config.num_passes)) {
      //Display(gSurface, gFramebuffer, gRenderPasses, config.width,
      //config.height);
      continue;
    }

    // Use Euler rotation.
    //printf("rot = %f, %f, %f\n", 180*gRotate[0]/M_PI, 180*gRotate[1]/M_PI,
    //180*gRotate[2]/M_PI);
    EulerToQuatRad(gCurrQuat, gRotate[2], gRotate[0], gRotate[1] + M_PI);
    //printf("quat = %f, %f, %f, %f\n", gCurrQuat[0], gCurrQuat[1],
    //gCurrQuat[2], gCurrQuat[3]);

    Render(scene, config, gImage, gEye, gLookat, gUp, gCurrQuat,
           gRenderPixelStep);

    // Always clear framebuffer for intermediate result
    if (gRenderPixelStep > 1) {
      ClearImage(gFramebuffer);
    }

    AccumImage(gFramebuffer, gImage);

    Display(gSurface, gFramebuffer, gRenderPasses, config.width, config.height);

    //printf("step = %d, interactive = %d\n", gRenderPixelStep,
    //gRenderInteractive);

    // Increment render pass.
    if (!gRenderInteractive && gRenderPixelStep == 1) {
      gRenderPasses++;
    }

    if (!gRenderInteractive) {
      gRenderPixelStep >>= 1;

      if (gRenderPixelStep == 1) {
        ClearImage(gFramebuffer);
      }

      if (gRenderPixelStep < 1) {
        gRenderPixelStep = 1;
      }
    } else {
      gRenderPasses = 1;
    }
#else

    if (hasEvent) {
      NotifyRenderClock();
    }

    SDL_Delay(33);
    {
      // Ensure render thread deson't write to a framebuffer.
      int ret = SDL_LockMutex(gMutex);
      assert(ret == 0);                // 0 = success
      SDL_RenderPresent(gSDLRenderer); // bitblit
      SDL_UnlockMutex(gMutex);
    }

#endif
  }

  NotifyRenderQuit();
  // renderThread.detach();
  renderThread.join();

  printf("\n");
  fflush(stdout); // for safety
}
}
#else  // ENABLE_SDL
namespace mallie {

void DoMainSDL(Scene &scene, const RenderConfig &config) {

  return;
}
}
#endif // ENABLE_SDL

namespace mallie {

#ifdef ENABLE_SDL
void PostRenderCancel()
{
  tthread::lock_guard<tthread::mutex> guard(gRenderThreadMutex);
  gRenderCancel = true;
}

void ClearRenderCancel()
{
  tthread::lock_guard<tthread::mutex> guard(gRenderThreadMutex);
  gRenderCancel = false;
}

bool CheckRenderCancel() {
  tthread::lock_guard<tthread::mutex> guard(gRenderThreadMutex);
  
  bool ret = gRenderCancel;

  return ret;
}
#else
void PostRenderCancel() {
}

void ClearRenderCancel() {
}

bool CheckRenderCancel() {
  return false;
}
#endif
}
