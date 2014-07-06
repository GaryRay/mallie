#include <iostream>
#include <fstream>
#include <cassert>

#include "material_loader.h"

#include "parson.h" // deps/parson

using namespace mallie;

bool
MaterialLoader::Load(
  std::vector<Material>& materials,
  const std::string& filename)
{

  { // file check
    std::ifstream is(filename.c_str());

    if (!is) {
      std::cerr << "File not found: " << filename << std::endl;
      return false;
    }
  }

  JSON_Value *root = json_parse_file(filename.c_str());
  if (json_value_get_type(root) != JSONArray) {
    std::cerr << "Root element must be array: " << filename << std::endl;
    return false;
  }

  JSON_Array* arr = json_value_get_array(root);
  int numMaterials = json_array_get_count(arr);

  materials.resize(numMaterials);

  for (int i = 0; i < numMaterials; i++) {
      Material mat;

      JSON_Object* m = json_array_get_object(arr, i);
      assert(m);

      int id = -1;
      if (json_value_get_type(json_object_dotget_value(m, "index")) ==
          JSONNumber) {
        id = json_object_dotget_number(m, "index");
      }
      assert(id > -1);
      assert(id < numMaterials);

      mat.id = id;
      printf("Mallie:debug\tid:%d\n", mat.id);

      if (json_value_get_type(json_object_dotget_value(m, "diffuse")) ==
          JSONArray) {
        JSON_Array *array = json_object_dotget_array(m, "diffuse");
        if (json_array_get_count(array) == 3) {
          mat.diffuse[0] = json_array_get_number(array, 0);
          mat.diffuse[1] = json_array_get_number(array, 1);
          mat.diffuse[2] = json_array_get_number(array, 2);
          printf("Mallie:debug\tdiffuse:%f, %f, %f\n", mat.diffuse[0], mat.diffuse[1], mat.diffuse[2]);
        }
      }

      if (json_value_get_type(json_object_dotget_value(m, "reflection")) ==
          JSONArray) {
        JSON_Array *array = json_object_dotget_array(m, "reflection");
        if (json_array_get_count(array) == 3) {
          mat.reflection[0] = json_array_get_number(array, 0);
          mat.reflection[1] = json_array_get_number(array, 1);
          mat.reflection[2] = json_array_get_number(array, 2);
          printf("Mallie:debug\treflection:%f, %f, %f\n", mat.reflection[0], mat.reflection[1], mat.reflection[2]);
        }
      }

      if (json_value_get_type(json_object_dotget_value(m, "reflection_glossiness")) ==
          JSONNumber) {
        mat.reflection_glossiness = json_object_dotget_number(m, "reflection_glossiness");
        printf("Mallie:debug\treflection_glossiness:%f\n", mat.reflection_glossiness);
      }

      materials[id] = mat;
  }

#if 0
  if (json_value_get_type(json_object_dotget_value(object, "index")) ==
      JSONNumber) {
    config.obj_filename = json_object_dotget_string(object, "obj_filename");
  }

  if (json_value_get_type(json_object_dotget_value(object, "eson_filename")) ==
      JSONString) {
    config.eson_filename = json_object_dotget_string(object, "eson_filename");
  }

  if (json_value_get_type(json_object_dotget_value(
          object, "material_filename")) == JSONString) {
    config.material_filename =
        json_object_dotget_string(object, "material_filename");
  }

  if (json_value_get_type(json_object_dotget_value(object, "scene_scale")) ==
      JSONNumber) {
    config.scene_scale = json_object_dotget_number(object, "scene_scale");
  }

  if (json_value_get_type(json_object_dotget_value(object, "eye")) ==
      JSONArray) {
    JSON_Array *array = json_object_dotget_array(object, "eye");
    if (json_array_get_count(array) == 3) {
      config.eye[0] = json_array_get_number(array, 0);
      config.eye[1] = json_array_get_number(array, 1);
      config.eye[2] = json_array_get_number(array, 2);
    }
  }

  if (json_value_get_type(json_object_dotget_value(object, "up")) ==
      JSONArray) {
    JSON_Array *array = json_object_dotget_array(object, "up");
    if (json_array_get_count(array) == 3) {
      config.up[0] = json_array_get_number(array, 0);
      config.up[1] = json_array_get_number(array, 1);
      config.up[2] = json_array_get_number(array, 2);
    }
  }

  if (json_value_get_type(json_object_dotget_value(object, "lookat")) ==
      JSONArray) {
    JSON_Array *array = json_object_dotget_array(object, "lookat");
    if (json_array_get_count(array) == 3) {
      config.lookat[0] = json_array_get_number(array, 0);
      config.lookat[1] = json_array_get_number(array, 1);
      config.lookat[2] = json_array_get_number(array, 2);
    }
  }

  if (json_value_get_type(json_object_dotget_value(object, "resolution")) ==
      JSONArray) {
    JSON_Array *array = json_object_dotget_array(object, "resolution");
    if (json_array_get_count(array) == 2) {
      config.width = json_array_get_number(array, 0);
      config.height = json_array_get_number(array, 1);
    }
  }

  if (json_value_get_type(json_object_dotget_value(object, "num_passes")) ==
      JSONNumber) {
    config.num_passes = json_object_dotget_number(object, "num_passes");
  }

  if (json_value_get_type(json_object_dotget_value(object, "num_photons")) ==
      JSONNumber) {
    config.num_passes = json_object_dotget_number(object, "num_photons");
  }

  if (json_value_get_type(json_object_dotget_value(object, "plane")) ==
      JSONBoolean) {
    config.plane = json_object_dotget_boolean(object, "plane");
  }

  if (json_value_get_type(json_object_dotget_value(object, "envmap_filename")) ==
      JSONString) {
    config.envmap_filename = json_object_dotget_string(object, "envmap_filename");
  }

  if (json_value_get_type(json_object_dotget_value(object, "envmap_coord")) ==
      JSONString) {
    config.envmap_coord = json_object_dotget_string(object, "envmap_coord");
  }

  if (json_value_get_type(json_object_dotget_value(object, "display_gamma")) ==
      JSONNumber) {
    config.display_gamma = json_object_dotget_number(object, "display_gamma");
  }
#endif

  json_value_free(root);

  return true;
}

