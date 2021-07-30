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

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color);
vec3 barycentric(const vec2 pts[3], const vec2 P);
void triangle(const vec3 pts[3], const vec2 uvs[3], double *zbuff, TGAImage &image, TGAImage &texture, vec3 ill_indensity);
void my_clear(TGAImage& image, const TGAColor color);
