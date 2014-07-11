-- Unit test
newoption {
   trigger     = "with-unittest",
   description = "Build unit test."
}

-- OSD patch eson
newoption {
   trigger     = "with-osd-patch",
   description = "Use OSD patch format."
}


-- ARM target
newoption {
   trigger     = "arm",
   description = "Compile for ARM."
}

-- OpenMP
newoption {
   trigger     = "with-openmp",
   description = "Use OpenMP."
}

-- MPI
newoption {
   trigger     = "with-mpi",
   description = "Use MPI."
}

-- Ptex
newoption {
   trigger     = "with-ptex",
   description = "Use Ptex library."
}

-- FITS
newoption {
   trigger     = "with-cfitsio",
   description = "Use C-binding of FITS IO."
}

-- SDL
newoption {
   trigger     = "with-sdl",
   description = "Use SDL2.0."
}

sources = {
   "main.cc",
   "main_console.cc",
   "main_sdl.cc",
   "material.h",
   "mmm_io.cc",
   "mmm_io.h",
   "render.cc",
   "render.h",
   "camera.cc",
   "camera.h",
   "light.cc",
   "matrix.cc",
   "trackball.cc",
   "importers/tiny_obj_loader.cc",
   "importers/tiny_obj_loader.h",
   "importers/eson.cc",
   "importers/eson.h",
   "importers/mesh_loader.cc",
   "importers/mesh_loader.h",
   "importers/material_loader.h",
   "importers/material_loader.cc",
   "bvh_accel.cc",
   "bvh_accel.h",
   "scene.cc",
   "scene.h",
   "spectrum.cc",
   "jpge.cc",
   "deps/parson/parson.c",
   "tasksys.cc",
   "texture.cc",
   "vcm.cc",
   "script_engine.cc",
   "script_engine.h",
   "stb_image.c",
   "shader.cc",
   "shader.h",
   "deps/TinyThread++-1.1/source/tinythread.cpp",
   "prim-plane.h",
   "prim-plane.cc",
   "feature-lines.h",
   "feature-lines.cc",
}

tinyjs_sources = {
   "deps/tinyjs/TinyJS.cpp",
   "deps/tinyjs/TinyJS_Functions.cpp",
   "deps/tinyjs/TinyJS_MathFunctions.cpp",
}

test_sources = {
   "test/cctest/test-atomic.cc"
}

gtest_sources = {
   "deps/gtest-1.7.0/src/gtest-all.cc",
   "deps/gtest-1.7.0/src/gtest_main.cc"
}

bezier_sources = {
   "bezier/bezier.cc",
   "bezier/bilinear_patch_intersection.cc",
   "bezier/bezier_patch_intersection.cc",
   "bezier/bvh_composite_intersection.cc",
   "bezier/patch_accel.cc",

   "importers/patch_loader.cc",
}


newaction {
   trigger     = "install",
   description = "Install the software",
   execute = function ()
      -- copy files, etc. here
   end
}

-- premake4.lua
solution "MallieSolution"
   configurations { "Release", "Debug" }

   if (os.is("windows")) then
      platforms { "native", "x32", "x64" }
   else
      platforms { "native", "x32", "x64" }
   end

   -- A project defines one build target
   project "Mallie"
      kind "ConsoleApp"
      language "C++"

      files { sources, tinyjs_sources, bezier_sources }

      includedirs {
         "./",
         "deps/parson/",
         "deps/TinyThread++-1.1/source/",
         "deps/tinyjs/",
         "bezier/",
      }

      if _OPTIONS['with-osd-patch'] then
         defines { 'ENABLE_OSD_PATCH' }
      end

      -- MacOSX. Guess we use gcc.
      configuration { "macosx", "gmake" }

         defines { '_LARGEFILE_SOURCE', '_FILE_OFFSET_BITS=64' }

         -- SDL
         if _OPTIONS['with-sdl'] then
            defines { 'ENABLE_SDL' }
            buildoptions { "`sdl2-config --cflags`" }
            buildoptions { "-msse2" }
            linkoptions { "`sdl2-config --libs`" }
         end

         -- gcc openmp
         if _OPTIONS['with-openmp'] then
            buildoptions { "-fopenmp" }
            linkoptions { "-fopenmp" }
         end

         -- gcc mpi
         if _OPTIONS['with-mpi'] then
            defines { 'WITH_MPI' }
         end

         -- Ptex
         if _OPTIONS['with-ptex'] then
            defines { 'ENABLE_PTEX' }
            includedirs { "./deps/ptex-master/install/include" }
            libdirs { "./deps/ptex-master/install/lib" }
            links   { "Ptex" }
         end

      -- Windows general
      configuration { "windows" }

         includedirs { "./compat" } -- stdint

         if _OPTIONS['with-sdl'] then
            defines { 'ENABLE_SDL' }
            includedirs { "extlibs/windows/SDL2-2.0.3/include/" }
            links { "SDL2", "SDL2main" }
            libdirs { "extlibs/windows/SDL2-2.0.3/lib/x64/" }
         end

         defines { '_USE_MATH_DEFINES' }
         defines { 'NOMINMAX', '_LARGEFILE_SOURCE', '_FILE_OFFSET_BITS=64' }

      -- Windows + gmake specific
      configuration { "windows", "gmake" }

         defines { '__STDC_CONSTANT_MACROS', '__STDC_LIMIT_MACROS' } -- c99

         links { "stdc++", "msvcrt", "ws2_32", "winmm" }

      -- Linux specific
      configuration {"linux", "gmake"}
         defines { '__STDC_CONSTANT_MACROS', '__STDC_LIMIT_MACROS' } -- c99

         -- gcc openmp
         if _OPTIONS['with-openmp'] then
            buildoptions { "-fopenmp" }
            linkoptions { "-fopenmp" }
         end

         -- gcc mpi
         if _OPTIONS['with-mpi'] then
            defines { 'WITH_MPI' }
         end

         -- Ptex
         if _OPTIONS['with-ptex'] then
            defines { 'ENABLE_PTEX' }
            includedirs { "./deps/ptex-master/install/include" }
            libdirs { "./deps/ptex-master/install/lib" }
            links   { "Ptex" }
         end

         if _OPTIONS['with-sdl'] then
            defines { 'ENABLE_SDL' }
            buildoptions { "`sdl2-config --cflags`" }
            linkoptions { "`sdl2-config --libs`" }
         end

      configuration "Debug"
         defines { "DEBUG" } -- -DDEBUG
         flags { "Symbols" }
         targetdir "bin/"
         targetname "mallie_debug"

      configuration "Release"
         defines { "NDEBUG" }
         flags { "Symbols", "Optimize" }
         if _OPTIONS['arm'] then
            -- buildoptions { "NEON" }
         else 
            flags { "EnableSSE2" }
         end
         targetdir "bin/"
         targetname "mallie"

   if _OPTIONS['with-unittest'] then
      dofile "premake4-test.lua"
   end

