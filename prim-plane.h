#ifndef __MALLIE_PRIM_PLANE_HPP__
#define __MALLIE_PRIM_PLANE_HPP__

#include "intersection.h"
#include "common.h"

namespace mallie {

class
Plane
{
public:
    Plane() {
    }

    ~Plane() {
    }

    void Set(float a, float b, float c, float d, int matID = 0) {
        m_a = a; m_b = b; m_c = c; m_d = d;
        m_matID = matID;
    }
 
    bool Intersect(Intersection *info, const Ray &ray);

    float m_a, m_b, m_c, m_d;
    int   m_matID;

};

}

#endif // __MALLIE_PLANE_HPP__
