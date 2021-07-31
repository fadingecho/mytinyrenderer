#include "tgaimage.h"
#include "geometry.h"
#include "model.h"
#include "global.h"
#include "mygl.h"

const int height = 800;
const int width = 800;
const int depth = 255;

vec3 light_dir(0, 0, 1);
vec3 eye_pos(0, 0, 0);
vec3 center(150, 300, 300);
vec3 up(0, 1, 0);

mat<4, 4> model_trans;
mat<4, 4> view_trans;
mat<4, 4> projection_trans;
mat<4, 4> mvp;
mat<4, 4> view_port;

Model *model = NULL;
const std::string output_dir = "../output/";

struct GouraudShader : public IShader{
    vec3 varying_intensity;//written by vertex shader, read by fragment shader

    virtual vec3 vertex(int iface, int nthvert){
        varying_intensity[nthvert] = std::max(0., model->normal(iface, nthvert)*light_dir);
        vec3 gl_vertex = model->vert(iface, nthvert);
        return trans_vec3(mvp, gl_vertex);
    }

    virtual bool fragment(vec3 bar, TGAColor &color){
        double intensity = varying_intensity*bar;
        color = TGAColor(255, 255, 255)*intensity;
        return false;   //discard(mask...)
    }

};

void show_model_origin(mat<4, 4> mvp, TGAImage image){
    vec3 origin = trans_vec3(mvp, vec3(0, 0, 0));
    vec3 xaxis = trans_vec3(mvp, vec3(1, 0, 0));
    vec3 yaxis = trans_vec3(mvp, vec3(0, 1, 0));
    vec3 zaxis = trans_vec3(mvp, vec3(0, 0, 1));
    line(int((origin.x + 1.)*width/2), int((origin.y + 1.)*height/2), int((xaxis.x + 1.)*width/2), int((xaxis.y + 1.)*height/2), image, red);
    line(int((origin.x + 1.)*width/2), int((origin.y + 1.)*height/2), int((yaxis.x + 1.)*width/2), int((yaxis.y + 1.)*height/2), image, green);
    line(int((origin.x + 1.)*width/2), int((origin.y + 1.)*height/2), int((zaxis.x + 1.)*width/2), int((zaxis.y + 1.)*height/2), image, blue);
}

int main(const int argc, const char** argv)
{
    model = new Model("../obj/head.obj");
    double *zbuffer = new double[width*height];
    for(int i = 0;i < width*height;i ++ ) zbuffer[i] = -std::numeric_limits<double>::max();

    double  eye_fov = 45, 
            aspect_ratio = width/height,
            zNear = -0.1,
            zFar = -1000,
            f1 = (zNear - zFar)/2.0,
            f2 = (zFar + zNear)/2.0;
 
    model_trans = get_model_trans(model_trans);
    view_trans = get_view(eye_pos, center, up);
    projection_trans = get_projection(eye_fov, aspect_ratio, zNear, zFar);
    mvp = projection_trans*view_trans*model_trans;
    view_port = get_viewport(0, 0, width, height);
    mvp = view_port*mvp;

    /*
    write model
    */
    TGAImage image(width, height, TGAImage::RGB), tex;
    my_clear(image, white);
    tex.read_tga_file("../obj/_diffuse.tga");

    GouraudShader shader;
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

    show_model_origin(mvp, image);
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
