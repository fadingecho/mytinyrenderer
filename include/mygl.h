#include "../include/tgaimage.h"
#include <vector>
#include <cmath>
#include "../include/geometry.h"
#include "../include/model.h"
#include <string>
#include <limits>

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor blue = TGAColor(0, 0, 255, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
const TGAColor black = TGAColor(0, 0, 0, 255);

struct IShader {
    mat<4, 4> u_viewport;   
    mat<4, 4> u_model;
    mat<4, 4> u_view;
    mat<4, 4> u_proj;
    mat<4, 4> u_vpI;    //viewport_invert
    mat<4, 4> u_projI;  //projection_invert
    virtual vec3 vertex(int iface, int nthvert) = 0;
    virtual bool fragment(vec3 bar, TGAColor &color) = 0;
};

vec<4> to_vec4(const vec<3> v3);
vec3 to_vec3(const vec4 v4);
vec3 trans_vec3(mat<4, 4> trans, vec3 coord, double fill);

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color);
vec3 barycentric(const vec2 pts[3], const vec2 P);
void triangle(const vec3 pts[3],IShader &shader, double *zbuff, TGAImage &image);
void my_clear(TGAImage& image, const TGAColor color);

mat<4, 4> get_scale(const double ratio);
mat<4, 4> get_translation(const vec3 world_coords);

mat<4, 4> get_viewport(int x, int y, int w, int h, int d);
mat<4, 4> get_model_trans(const double ratio, const vec3 world_coords);
mat<4, 4> get_view(const vec3 eye_pos, const vec3 center, const vec3 up);
mat<4, 4> get_projection(double eye_fov, double aspect_ratio, double zNear, double zFar);