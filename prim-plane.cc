#include <cmath>
#include <limits>

#include "prim-plane.h"

using namespace mallie;

bool
Plane::Intersect(
    Intersection *info,
    const Ray& ray)
{
    real3 n(m_a, m_b, m_c);
    real3 v = ray.dir;
    v.normalize();
    float vn = vdot(v, n);

    if (std::abs(vn) > std::numeric_limits<float>::epsilon() * 1024.0f) {
        float on_d = vdot(ray.org, n) + m_d;
        //float t = -on_d / vn;
        float inv_vn = 1.0f / vn;
        float t = -on_d * inv_vn;

        if ((t > 0) && (t < info->t)) {

            info->t = t;
            info->position = ray.org + t * v;

            n.normalize();
            info->geometricNormal = n;
            info->normal          = n;
            info->tangent[0] = 1.0; // @fixme
            info->tangent[1] = 0.0;
            info->tangent[2] = 0.0;
            info->binormal[0] = 0.0; // @fixme
            info->binormal[1] = 0.0;
            info->binormal[2] = -1.0;
            //info->st[0]         = 0.0; 
            //info->st[1]         = 0.0; 
            info->matID       = m_matID; 

            return true;
        }

        return false;

    } else {
        return false;
    }
}
