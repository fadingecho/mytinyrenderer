#include "tgaimage.h"
#include "geometry.h"
#include "model.h"
#include "global.h"
#include "mygl.h"

const double MY_PI = 3.1415926;

const int height = 800;
const int width = 800;

Model *model = NULL;
const std::string output_dir = "../output/";

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

vec3 mvp_trans(mat<4, 4> mvp, vec3 coord){
    return to_vec3(mvp*to_vec4(coord));
}

mat<4, 4> get_model_trans(const mat<4, 4> model_pos){
    mat<4, 4> m = mat<4, 4>::identity();

    mat<4, 4> scale = mat<4, 4>::identity()*100;
    scale[3][3] = 1;

    double angle = 40;
    angle = angle * MY_PI / 180.f;
    mat<4, 4> rot = mat<4, 4>::identity();
    rot[0][0] = cos(angle);
    rot[0][2] = sin(angle);
    rot[2][0] = -sin(angle);
    rot[2][2] = cos(angle);
    return rot*scale*m;
}

mat<4, 4> get_viewport(const vec3 eye_pos){
    mat<4, 4> m = mat<4, 4>::identity();
    m[0][3] = -eye_pos[0];
    m[1][3] = -eye_pos[1];
    m[2][3] = -eye_pos[2];
    return m;
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

void draw(){
    model = new Model("../obj/head.obj");
    double *zbuffer = new double[width*height];
    for(int i = 0;i < width*height;i ++ ) zbuffer[i] = -std::numeric_limits<double>::max();

    // vec3 light_dir(pow(3., -0.5), pow(3., -0.5), pow(3., -0.5));   //shoot from infinite far place
    vec3 light_dir(0, 0, 1);
    vec3 eye_pos(0, 0, 400);

    double  eye_fov = 45, 
            aspect_ratio = width/height,
            zNear = -0.1,
            zFar = -1000,
            f1 = (zNear - zFar)/2.0,
            f2 = (zFar + zNear)/2.0;
 
    mat<4, 4> model_trans = get_model_trans(model_trans);
    mat<4, 4> viewport_trans = get_viewport(eye_pos);
    mat<4, 4> projection_trans = get_projection(eye_fov, aspect_ratio, zNear, zFar);
    mat<4, 4> mvp = mat<4, 4>::identity();
    mvp = projection_trans*viewport_trans*model_trans*mvp;

    /*
    write model
    */
    TGAImage image(width, height, TGAImage::RGB);
    my_clear(image, white);
    TGAImage tex;
    tex.read_tga_file("../obj/_diffuse.tga");
    for(int i = 0;i <model->nfaces();i ++ ){
        //for everyface, get vertex and draw lines
        vec3 screen_coords[3];
        vec3 world_coords[3];
        vec2 tex_coords[3];
        for(int j = 0;j < 3;j ++ ){
            world_coords[j] = model->vert(i, j);
        }
        vec3 ill_indensities;
        for(int j = 0;j < 3;j ++ ){
                world_coords[j] = mvp_trans(mvp, world_coords[j]);
                tex_coords[j] = model->uv(i, j);

                vec3 v = world_coords[j];
                //pay attention to the convert
                //it is really important
                screen_coords[j] = vec3(int((v.x + 1.)*width/2), int((v.y + 1)*height/2), (v.z*f1 + f2)); 
                //std::cout << v.z*f1 + f2 << std::endl;
                //note that the sampled coords are upset-down
                //only when light indensity above zero
                tex_coords[j].x = 1 - tex_coords[j].x;
                tex_coords[j].y = 1 - tex_coords[j].y;

                //gouraud shading
                ill_indensities[j] = model->normal(i, j).normalize() * light_dir;
        }
        triangle(screen_coords, tex_coords, zbuffer, image, tex, ill_indensities);
    }

    vec3 origin = mvp_trans(mvp, vec3(0, 0, 0));
    vec3 xaxis = mvp_trans(mvp, vec3(1, 0, 0));
    vec3 yaxis = mvp_trans(mvp, vec3(0, 1, 0));
    vec3 zaxis = mvp_trans(mvp, vec3(0, 0, 1));
    line(int((origin.x + 1.)*width/2), int((origin.y + 1.)*height/2), int((xaxis.x + 1.)*width/2), int((xaxis.y + 1.)*height/2), image, red);
    line(int((origin.x + 1.)*width/2), int((origin.y + 1.)*height/2), int((yaxis.x + 1.)*width/2), int((yaxis.y + 1.)*height/2), image, green);
    line(int((origin.x + 1.)*width/2), int((origin.y + 1.)*height/2), int((zaxis.x + 1.)*width/2), int((zaxis.y + 1.)*height/2), image, blue);

    printf("model grid done\n");

    //image.flip_vertically();
    if(image.write_tga_file(output_dir + "output.tga"))
        std::cout << "file output to " << output_dir << std::endl;
    else 
        std::cout << "can not output file" << std::endl;

    /*
    write zbuffer
    */
    std::cout << "depth buffer" << std::endl;
    TGAImage buffer(width, height, TGAImage::RGBA);
    for(int i = 0;i < width;i ++ )
        for(int j = 0;j < height;j ++ )
        {
            double ratio = zbuffer[i + j*width]*10;
            buffer.set(i, j, TGAColor(255*ratio, 255*ratio, 255*ratio, 255));
        }
    if(buffer.write_tga_file(output_dir + "zbuffer.tga"))
        std::cout << "zbuffer output to " << output_dir << std::endl;
    else 
        std::cout << "can not output zbuffer" << std::endl;
    
    delete model;
}

int main(const int argc, const char** argv)
{
    printf("start render\n");
    draw();
    return 0;
}
