#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>
#include "geometry.h"

class Model {
private:
	std::vector<vec3> verts_;
	std::vector<vec2> uv_;        // array of tex coords
    std::vector<vec3> norms_;     // array of normal vectors
	std::vector<int> facet_vrt_; 
    std::vector<int> facet_tex_;  // indices in the above arrays per triangle
    std::vector<int> facet_nrm_;
public:
	Model(const char *filename);
	int nverts() const;
    int nfaces() const;
    vec3 normal(const int iface, const int nthvert) const;  // per triangle corner normal vertex
    vec3 normal(const vec2 &uv) const;                      // fetch the normal vector from the normal map texture
    vec3 vert(const int i) const;
    vec3 vert(const int iface, const int nthvert) const;
    vec2 uv(const int iface, const int nthvert) const;
};

#endif //__MODEL_H__