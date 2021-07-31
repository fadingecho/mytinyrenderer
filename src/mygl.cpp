#include "mygl.h"
#define MY_PI 3.1415926

void my_clear(TGAImage& image, const TGAColor color){
    for(int i = 0;i < image.get_width();i ++ )
        for(int j = 0;j < image.get_height();j ++ ){
            image.set(i, j, color);
        }
}

vec<4> to_vec4(const vec<3> v3){
    vec4 v4;
    v4[0] = v3[0];
    v4[1] = v3[1];
    v4[2] = v3[2];
    v4[3] = 1;
    return v4;
}

vec3 to_vec3(const vec4 v4){
    double w = 1/v4[3];
    return vec3(v4[0]*w, v4[1]*w, v4[2]*w);
}

vec3 trans_vec3(mat<4, 4> trans, vec3 coord){
    return to_vec3(trans*to_vec4(coord));
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

void triangle(const vec3 pts[3],IShader &shader, double *zbuff, TGAImage &image){
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
    vec3 P;
    for(P.x =bboxmin.x;P.x < bboxmax.x;P.x ++ )
        for(P.y = bboxmin.y;P.y < bboxmax.y;P.y ++ )
        {
            vec3 bc_screen= barycentric(pts2, vec2(P.x, P.y));
            if(bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;
            P.z = interpolation(bc_screen, pts[0].z, pts[1].z, pts[2].z);
            if(P.z > zbuff[(int)(P.x + P.y*image.get_width())]){
                TGAColor color;
                bool discard = shader.fragment(bc_screen, color);
                if(discard) continue;
                image.set(P.x, P.y, color);
                zbuff[(int)(P.x + P.y*image.get_width())] = P.z;
            } 
        }
}

mat<4, 4> get_viewport(int x, int y, int w, int h){
    mat<4, 4> m = mat<4, 4>::identity();
    m[0][3] = x+w/2.f;
    m[1][3] = y+h/2.f;
    m[2][3] = 255/2.f;

    m[0][0] = w/2.f;
    m[1][1] = h/2.f;
    m[2][2] = 255/2.f;
    return m;
}

mat<4, 4> get_model_trans(const mat<4, 4> model_pos){
    mat<4, 4> m = mat<4, 4>::identity();

    mat<4, 4> scale = mat<4, 4>::identity()*100;
    scale[3][3] = 1;

    double angle = 0;
    angle = angle * MY_PI / 180.f;
    mat<4, 4> rot = mat<4, 4>::identity();
    rot[0][0] = cos(angle);
    rot[0][2] = sin(angle);
    rot[2][0] = -sin(angle);
    rot[2][2] = cos(angle);
    return rot*scale*m;
}

mat<4, 4> get_view(const vec3 eye_pos, const vec3 center, const vec3 up){
    vec3 z = (center - eye_pos).normalize(),
         x = cross(up, z).normalize(),
         y = cross(z, x).normalize();

    mat<4, 4> to_center = mat<4, 4>::identity(),
              base = mat<4, 4>::identity();

    for(int i = 0;i < 3;i ++ ){
        to_center[i][3] = -center[i];
        base[0][i] = x[i];
        base[1][i] = y[i];
        base[2][i] = z[i];
    }
    
    return base*to_center;
}

mat<4, 4> get_projection(double eye_fov, double aspect_ratio, double zNear, double zFar){
    mat<4, 4> projection = mat<4, 4>::identity();
    //the way I create these translation is same as games101, not tinyrenderer
    double l, r, b, t;
    t = -std::tan(eye_fov/360.0*MY_PI) * std::abs(zNear);
    r = aspect_ratio * t;
    l = -r;
    b = -t;

    mat<4,4> orth_l = mat<4,4>::identity();
    orth_l[0][0] = 2/(r - l);
    orth_l[1][1] = 2/(b - t);
    orth_l[2][2] = 2/(zNear - zFar);
    
    mat<4,4> orth_r = mat<4,4>::identity();
    orth_r[0][3] = -(r + l)/2;
    orth_r[1][3] = -(t + b)/2;
    orth_r[2][3] = -(zNear + zFar)/2;

    mat<4,4> to_orth = mat<4,4>::zero();
    to_orth[0][0] = zNear;
    to_orth[1][1] = zNear;
    to_orth[2][2] = zNear + zFar;
    to_orth[2][3] = -zNear*zFar;
    to_orth[3][2] = 1;
    return orth_l * orth_r * to_orth * projection;
}