#include "tgaimage.h"
#include "geometry.h"
#include "model.h"
#include "global.h"
#include "mygl.h"

const int height = 1024;
const int width = 1024;
const int depth = 255;

vec3 light_dir(0, 0, -1);
vec3 eye_pos(0, 0, 0);
vec3 center(0, 0, 300);
// vec3 center(100, 100, 270);
vec3 up(0, 1, 0);

mat<4, 4> model_trans;
mat<4, 4> view_trans;
mat<4, 4> projection_trans;
mat<4, 4> mvp;
mat<4, 4> view_port;

Model *model = NULL;
const std::string output_dir = "../output/";
TGAImage tex0;
TGAImage tex1;

TGAColor show_nv3(vec3 v3){
    return TGAColor(v3.x*255, v3.y*255, v3.z*255);
}

struct GouraudShader : public IShader{
    public:
    mat<2, 3> v_uv;//written by vertex shader, read by fragment shader
    mat<4, 4> u_vp_mvp;
    mat<4, 4> u_vp_mvp_it;

    mat<4, 4> u_model;
    mat<4, 4> u_model_it;

    vec3 v_ver[3];  //for barycentric coords in world space

    vec2 v_t0;
    vec2 v_t1;

    virtual vec3 vertex(int iface, int nthvert){
        v_uv.set_col(nthvert, model->uv(iface, nthvert));
        vec3 gl_vertex = model->vert(iface, nthvert);

        v_t0 = vec2(tex0.get_width(), tex0.get_height());
        v_t1 = vec2(tex1.get_width(), tex1.get_height());

        v_ver[nthvert] = gl_vertex;
        return trans_vec3(u_vp_mvp, gl_vertex);
    }

    virtual bool fragment(vec3 bar, TGAColor &color){
        vec3 l, n;
        vec2 b_uv = v_uv*bar;
        TGAColor c = tex1.get(v_t1.x*(b_uv.x), v_t1.y*(1-b_uv.y));
        n = vec3(((double)c[2])/255*2-1, ((double)c[1])/255*2-1, ((double)c[0]/255*2-1));   //according to tinyredner/model.cpp, cant understand it
        c = tex0.get(v_t0.x*(b_uv.x), v_t0.y*(1-b_uv.y));

        l = trans_vec3(u_model, light_dir).normalize(); 
        n = trans_vec3(u_model_it, n).normalize();
        double intensity = std::max<double>(0.0, -(l*n));
        color = c*intensity;
        return false;   //discard(mask...)
    }

};

int main(const int argc, const char** argv)
{
    model = new Model("../obj/diablo3_pose/object.obj");
    double *zbuffer = new double[width*height];
    for(int i = 0;i < width*height;i ++ ) zbuffer[i] = -std::numeric_limits<double>::max();

    double  eye_fov = 45, 
            aspect_ratio = width/height,
            zNear = -0.1,
            zFar = -1000,
            f1 = (zNear - zFar)/2.0,
            f2 = (zFar + zNear)/2.0;
 
    model_trans         = get_model_trans(model_trans);
    view_trans          = get_view(eye_pos, center, up);
    projection_trans    = get_projection(eye_fov, aspect_ratio, zNear, zFar);
    mvp = projection_trans*view_trans*model_trans;
    view_port           = get_viewport(0, 0, width, height, depth);

    /*
    write model
    */
    TGAImage image(width, height, TGAImage::RGB), tex;
    my_clear(image, white);
    tex0.read_tga_file("../obj/diablo3_pose/_diffuse.tga");
    tex1.read_tga_file("../obj/diablo3_pose/_nm.tga");

    GouraudShader shader;
    shader.u_vp_mvp = view_port*mvp;
    shader.u_vp_mvp_it = (view_port*mvp).invert_transpose();

    shader.u_model =(view_trans*model_trans);
    shader.u_model_it = (view_trans*model_trans).invert_transpose();
    for(int i = 0;i <model->nfaces();i ++ ){
        //for everyface, get vertex and draw lines
        vec3 screen_coords[3];
        for(int j = 0;j < 3;j ++ ){
            screen_coords[j] = shader.vertex(i, j);
            screen_coords[j].x = (int)screen_coords[j].x;
            screen_coords[j].y = (int)screen_coords[j].y;
        }
        triangle(screen_coords, shader, zbuffer, image);
    }

    vec3 origin = trans_vec3(view_port*mvp, vec3(0, 0, 0));
    vec3 xaxis = trans_vec3(view_port*mvp, vec3(1, 0, 0));
    vec3 yaxis = trans_vec3(view_port*mvp, vec3(0, 1, 0));
    vec3 zaxis = trans_vec3(view_port*mvp, vec3(0, 0, 1));
    line(int(origin.x), int(origin.y), int(xaxis.x), int(xaxis.y), image, red);
    line(int(origin.x), int(origin.y), int(yaxis.x), int(yaxis.y), image, green);
    line(int(origin.x), int(origin.y), int(zaxis.x), int(zaxis.y), image, blue);

    if(image.write_tga_file(output_dir + "output.tga"))
        std::cout << "file output to " << output_dir << std::endl;
    else 
        std::cout << "can not output file" << std::endl;

    /*
    write zbuffer
    */
    TGAImage buffer(width, height, TGAImage::RGBA);
    for(int i = 0;i < width;i ++ )
        for(int j = 0;j < height;j ++ )
        {
            double ratio = zbuffer[i + j*width];
            buffer.set(i, j, TGAColor(255*ratio, 255*ratio, 255*ratio, 255));
        }
    if(buffer.write_tga_file(output_dir + "zbuffer.tga"))
        std::cout << "zbuffer output to " << output_dir << std::endl;
    else 
        std::cout << "can not output zbuffer" << std::endl;
    
    delete model;
    return 0;
}
