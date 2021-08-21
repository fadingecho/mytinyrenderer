#include "mygl.h"
#define MY_PI 3.1415926

void my_clear(TGAImage& image, const TGAColor color){
    for(int i = 0;i < image.get_width();i ++ )
        for(int j = 0;j < image.get_height();j ++ ){
            image.set(i, j, color);
        }
}

// vec<4> to_vec4(const vec<3> v3){
//     vec4 v4;
//     v4[0] = v3[0];
//     v4[1] = v3[1];
//     v4[2] = v3[2];
//     v4[3] = 1;
//     return v4;
// }

// vec3 to_vec3(const vec4 v4){
//     double w = 1/v4[3];
//     return vec3(v4[0]*w, v4[1]*w, v4[2]*w);
// }

vec3 trans_vec3(const mat<4, 4> trans, const vec3 coord, double fill){
    vec4 v4 = trans*embed<4>(coord, fill);
    if(fill > std::numeric_limits<double>::epsilon()){
        double w = 1/v4[3];
        return proj<3>(v4)*w;
    }
    return proj<3>(v4);
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

template<int n> vec<n> interpolation(vec3 bcoords, vec<n> v1, vec<n> v2, vec<n> v3){
    assert(n > 0);
    vec<n> ret = bcoords[0]*v1;
    ret = ret + bcoords[1]*v2;
    ret = ret + bcoords[2]*v3;
    return ret;
}

vec3 barycentric(const vec2 pts[3], const vec2 P){
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

void triangle(const vec3 screen_coords[3],IShader &shader, double *zbuff, TGAImage &image){
    vec3 world_coords[3];// world coords in the paras are coords after camera trans, actually
    for(int i = 0;i < 3;i ++ ){
        world_coords[i] = trans_vec3(shader.u_projI*shader.u_vpI, screen_coords[i], 1.0);
    }
    vec2 wc2d[3] = {vec2(world_coords[0].x, world_coords[0].y),vec2(world_coords[1].x, world_coords[1].y),vec2(world_coords[2].x, world_coords[2].y)};
    vec2 sc2d[3] = {vec2(screen_coords[0].x, screen_coords[0].y),vec2(screen_coords[1].x, screen_coords[1].y),vec2(screen_coords[2].x, screen_coords[2].y)};

    vec2 clamp(image.get_width() - 1, image.get_height() - 1);
    vec2 bboxmin(clamp);
    vec2 bboxmax(std::numeric_limits<double>::epsilon(), std::numeric_limits<double>::epsilon());//initting it by 0 causes unknow error when light pos is (0, y, 0) etc.
    
    for(int i = 0;i < 3;i ++ )
        for(int j = 0;j < 2;j ++ )
        {
            bboxmin[j] = std::min(bboxmin[j] , sc2d[i][j]);
            bboxmax[j] = std::max(bboxmax[j], sc2d[i][j]);
        }
    
    vec2 uv;
    for(int x =(int)bboxmin.x;x < (int)bboxmax.x+1;x ++ )
        for(int y = (int)bboxmin.y;y < (int)bboxmax.y+1;y ++ )
        {
            vec3 bc_screen= barycentric(sc2d, vec2(x, y));
            double z = vec3(screen_coords[0][2], screen_coords[1][2], screen_coords[2][2])*bc_screen;
            vec3 v = trans_vec3(shader.u_projI*shader.u_vpI, vec3(x, y, z), 1.0);
            vec3 bc_world = barycentric(wc2d, vec2(v.x, v.y));
            
            if(z < zbuff[x + y*image.get_width()] || bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;
            TGAColor color;
            bool discard = shader.fragment(bc_world, color);
            if(discard) continue;
            image.set(x, y, color);
            zbuff[x + y*image.get_width()] = z;
        }
}

mat<4, 4> get_viewport(int x, int y, int w, int h, int d){
    mat<4, 4> m = mat<4, 4>::identity();
    m[0][3] = x+w/2.f;
    m[1][3] = y+h/2.f;
    m[2][3] = d/2.f;

    m[0][0] = w/2.f;
    m[1][1] = h/2.f;
    m[2][2] = d/2.f;
    return m;
}

mat<4, 4> get_scale(const double ratio){
    mat<4, 4> scale = mat<4, 4>::identity()*ratio;
    scale[3][3] = 1;
    return scale;
}

mat<4, 4> get_translation(const vec3 world_coords){
    mat<4, 4> translation = mat<4, 4>::identity();
    translation.set_col(3, embed<4>(world_coords));
    return translation;
}

mat<4, 4> get_model_trans(const double ratio, const vec3 world_coords){
    double angle = 0;
    angle = angle * MY_PI / 180.f;
    mat<4, 4> rot = mat<4, 4>::identity();
    rot[0][0] = cos(angle);
    rot[0][2] = sin(angle);
    rot[2][0] = -sin(angle);
    rot[2][2] = cos(angle);
    return get_translation(world_coords)*get_scale(ratio);
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
    r = -aspect_ratio * t;
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