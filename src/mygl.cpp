#include "../include/tgaimage.h"
#include <vector>
#include <cmath>
#include "../include/geometry.h"
#include "../include/model.h"
#include "../include/mygl.h"
#include <string>
#include <limits>

void my_clear(TGAImage& image, const TGAColor color){
    for(int i = 0;i < image.get_width();i ++ )
        for(int j = 0;j < image.get_height();j ++ ){
            image.set(i, j, color);
        }
}

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color){
    bool steep = false;
    if(std::abs(x0-x1) < std::abs(y0-y1))//if steep, transpose the image
    {
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    
    if(x0>x1)//for symmetry
    {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }

    int dx = x1-x0;
    int dy = y1-y0;
    int derror = std::abs(dy)*2;
    int error = 0;
    int y = y0;

    for(int x=x0;x<=x1;x++)
    {
        if(steep)
        {
            image.set(y, x, color);
        }
        else{
            image.set(x, y, color);
        }

        error += derror;
        if(error > dx){
            y += (y1>y0?1:-1);
            error -= dx*2;
        }
    }
}

double interpolation(vec3 bcoords, double v1, double v2, double v3){
    // std::cout << bcoords[0]+bcoords[1]+bcoords[2] << std::endl;
    double ret = bcoords[0]*v1;
    ret += bcoords[1]*v2;
    ret += bcoords[2]*v3;
    return ret;
}

template<int n> vec<n> interpolation(vec3 bcoords, vec<n> v1, vec<n> v2, vec<n> v3){
    assert(n > 0);
    vec<n> ret = bcoords[0]*v1;
    ret = ret + bcoords[1]*v2;
    ret = ret + bcoords[2]*v3;
    return ret;
}

vec3 barycentric(const vec2 pts[3], const vec2 P){
    // vec2 v0 = pts[1] - pts[0], v1 = pts[2] - pts[0], v2 = P - pts[0];
    // double d00 = v0*v0;
    // double d01 = v0*v1;
    // double d11 = v1*v1;
    // double d20 = v2*v0;
    // double d21 = v2*v1;
    // double denom = d00*d11 - d01*d01;
    // double v = (d11 * d20 - d01 * d21) / denom;
    // double w = (d00 * d21 - d01 * d20) / denom;
    // double u = 1.0f - v - w;
    // return vec3(u, v, w);

    vec2 A = pts[0], B = pts[1], C = pts[2];
    vec3 x(C[0] - A[0], B[0] - A[0], A[0] - P[0]);
    vec3 y(C[1] - A[1], B[1] - A[1], A[1] - P[1]);

    vec3 u = cross(x, y); // u 向量和 x y 向量的点积为 0，所以 x y 向量叉乘可以得到 u 向量

    // 由于 A, B, C, P 的坐标都是 int 类型，所以 u[2] 必定是 int 类型
    // 如果 u[2] 为 0，则表示三角形 ABC 退化了（退还为直线 or 一个点），需要对其舍弃
    if (std::abs(u[2]) > 1e-2) {
        return vec3(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
    }

    // in this case generate negative coordinates, it will be thrown away by the rasterizator
    return vec3(-1, 1, 1);
}

void triangle(const vec3 pts[3], const vec2 uvs[3], double *zbuff, TGAImage &image, TGAImage &texture, vec3 ill_indensities){
    vec2 pts2[3] = {vec2(pts[0].x, pts[0].y),vec2(pts[1].x, pts[1].y),vec2(pts[2].x, pts[2].y)};

    vec2 bboxmin(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    vec2 bboxmax(-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());
    vec2 clamp(image.get_width() - 1, image.get_height() - 1);
    for(int i = 0;i < 3;i ++ )
        for(int j = 0;j < 2;j ++ )
        {
            bboxmin[j] = std::min(bboxmin[j] , pts[i][j]);
            bboxmax[j] = std::max(bboxmax[j], pts[i][j]);
        }
    for(int i = 0;i < 2;i ++ ) {
        bboxmin[i] = std::max(0., bboxmin[i]);
        bboxmax[i] = std::min(clamp[i], bboxmax[i]);
    }
    
    vec2 uv;
    int tw = texture.get_width(), th = texture.get_height();
    for(P.x = bboxmin.x;P.x < bboxmax.x;P.x ++ )
        for(P.y = bboxmin.y;P.y < bboxmax.y;P.y ++ )
        {
            vec3 bc_screen= barycentric(pts2, vec2(P.x, P.y));
            if(bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;
            P.z = interpolation(bc_screen, pts[0].z, pts[1].z, pts[2].z);

            if(P.z > zbuff[(int)(P.x + P.y*image.get_width())]){
                // image.set(P.x, P.y, white);
                // image.set(P.x, P.y, TGAColor(255*ill_indensity, 255*ill_indensity, 255*ill_indensity, 255));
                double ill_indensity = interpolation(bc_screen, ill_indensities[0], ill_indensities[1], ill_indensities[2]);
                TGAColor color = black;
                vec3 color_v = vec3(0, 0, 0);
                if(ill_indensity > 0){
                    uv = interpolation(bc_screen, uvs[0], uvs[1], uvs[2]);
                    uv.x *= th;
                    uv.y *= tw;
                    color = texture.get(uv.x, uv.y)*ill_indensity;
                }
              
                image.set(P.x, P.y, color);
                zbuff[(int)(P.x + P.y*image.get_width())] = P.z;
            } 
        }
}