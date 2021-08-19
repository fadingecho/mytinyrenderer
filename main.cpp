#include "tgaimage.h"
#include "geometry.h"
#include "model.h"
#include "global.h"
#include "mygl.h"
#include "string"

const int height = 800;
const int width = 800;
const int depth = 255;

vec3 light_dir(-1, 0, 0);
vec3 light_pos(300, 300, 0);
vec3 eye_pos(0, 0, 0);
vec3 center(180, 150, 300);
vec3 up(0, 1, 0);

double model_scale[2] = {100, 100};
vec3 model_poi[2] = {
    vec3(0, 0, 0),
    vec3(0, 0, -100)
};

std::string texture_path[2][3] = {
    {
        "../obj/diablo3_pose/_diffuse.tga",
        "../obj/diablo3_pose/_nm_tangent.tga",
        "../obj/diablo3_pose/_spec.tga"
    },
    {
        "../obj/floor_diffuse.tga",
        "../obj/floor_nm_tangent.tga",
        ""
    }
}; 

std::string model_path[2] = {
    "../obj/diablo3_pose/object.obj",
    "../obj/floor.obj"
};

mat<4, 4> model_trans;
mat<4, 4> view_trans;
mat<4, 4> projection_trans;
mat<4, 4> mvp;
mat<4, 4> view_port;

Model *model = NULL;
const std::string output_dir = "../output/";
TGAImage tex[3];

TGAColor show_nv3(vec3 v3){
    double s = v3.x + v3.y + v3.z;
    return TGAColor(v3.x*255/s, v3.y*255/s, v3.z*255/s);
}

struct DepthShader : public IShader{
    public:
    mat<3, 3> v_tri;
    
    virtual vec3 vertex(int iface, int nthvert){
        vec3 gl_vertex = trans_vec3(u_view*u_model, model->vert(iface, nthvert), 1.0);
        v_tri.set_col(nthvert, trans_vec3(u_viewport*u_proj, gl_vertex, 1.0));
        return gl_vertex;
    }

    virtual bool fragment(vec3 bar, TGAColor &color){
        double ratio = (v_tri*bar).z;
        // std::cout << ratio << std::endl;
        color = TGAColor(255*ratio, 255*ratio, 255*ratio, 255);
        return false;   //discard
    } 
};

struct BlingPhongShader : public IShader{
    public:
    mat<3, 3> v_tri;
    mat<2, 3> v_uv;//written by vertex shader, read by fragment shader
    mat<3, 3> v_nrm;//fs will interpolate normal
    mat<3, 3> ndc_tri;

    virtual vec3 vertex(int iface, int nthvert){
        vec3 gl_vertex = trans_vec3(u_view*u_model, model->vert(iface, nthvert), 1.0);

        v_uv.set_col(nthvert, model->uv(iface, nthvert));
        v_nrm.set_col(nthvert, proj<3>((u_view*u_model).invert_transpose()*embed<4>(model->normal(iface, nthvert), 0.0)));
        ndc_tri.set_col(nthvert, gl_vertex);

        v_tri.set_col(nthvert, gl_vertex);
        return gl_vertex;
    }

    virtual bool fragment(vec3 bar, TGAColor &color){
        vec3 l, n, r; //light dir, normal, reflected light dir
        vec2 b_uv = v_uv*bar;

        TGAColor c = tex[1].get(tex[1].get_width()*(b_uv.x), tex[1].get_height()*(1-b_uv.y));
        n = vec3((double)c[2]/255*2-1, (double)c[1]/255*2-1, (double)c[0]/255*2-1).normalize();   //according to tinyredner/model.cpp, cant understand it
        {//TBN
            vec3 b_n = (v_nrm*bar).normalize();
        
            mat<3,3> A;
            A[0] = ndc_tri.col(1) - ndc_tri.col(0);
            A[1] = ndc_tri.col(2) - ndc_tri.col(0);
            A[2] = b_n;//bn

            mat<3,3> AI = A.invert();
            vec3 i = AI * vec3(v_uv[0][1] - v_uv[0][0], v_uv[0][2] - v_uv[0][0], 0);
            vec3 j = AI * vec3(v_uv[1][1] - v_uv[1][0], v_uv[1][2] - v_uv[1][0], 0);
            
            mat<3, 3> B;
            B.set_col(0, i.normalize());
            B.set_col(1, j.normalize());
            B.set_col(2, b_n);

            n = (B*n).normalize();
        }
        
        // l = trans_vec3(u_view*u_model, (light_dir).normalize(), 0.0).normalize(); 
        l = trans_vec3(u_view*u_model, (ndc_tri*bar - light_pos).normalize(), 0.0).normalize();
        double nl = -(l*n);
    
        r = (n*(nl*2.0) + l).normalize();

        double diff = std::max<double>(nl, 0.0);
        double spec = 1-pow(1-std::max(r.z, 0.0), (double)tex[2].get(tex[2].get_width()*(b_uv.x), tex[2].get_height()*(1-b_uv.y))[0]);
        double ambi = 5;

        c = tex[0].get(tex[0].get_width()*(b_uv.x), tex[0].get_height()*(1-b_uv.y));
        for (int i=0; i<3; i++) color[i] = std::min<double>(ambi + (double)c[i]*(diff + 0.76*spec), 255);
      
        return false;   //discard
    }
};

void load(uint8_t cnt){
    for(int i = 0;i < 3;i ++ ){
        if(texture_path[cnt][i].empty()){
            std::cout << "texture#" << i << " is empty" << std::endl;
            continue;
        }
        if(!tex[i].read_tga_file(texture_path[cnt][i])) std::cout << "cant read texture#" << i << "in" << texture_path[cnt][i] << std::endl;
    }
    model = new Model(model_path[cnt].c_str());
}

int main(const int argc, const char** argv)
{
    TGAImage image(width, height, TGAImage::RGB);
    my_clear(image, black);
    double *zbuffer = new double[width*height];
    for(int i = 0;i < width*height;i ++ ) zbuffer[i] = -std::numeric_limits<double>::max();
    double  eye_fov = 45, 
            aspect_ratio = width/height,
            zNear = -1,
            zFar = -800,
            f1 = (zNear - zFar)/2.0,
            f2 = (zFar + zNear)/2.0;
    
    view_trans          = get_view(eye_pos, center, up);
    projection_trans    = get_projection(eye_fov, aspect_ratio, zNear, zFar);
    view_port           = get_viewport(0, 0, width, height, depth);
    for(uint8_t cnt = 0;cnt < 2;cnt ++ ){
        std::cout << "process model #" << char(cnt+'0') << std::endl;
        load(cnt);
        model_trans         = get_model_trans(model_scale[cnt], model_poi[cnt]);
        mvp                 = projection_trans*view_trans*model_trans;
        
        if(true){
            DepthShader shader;
            shader.u_viewport   = view_port;
            shader.u_model      = model_trans;   //to_center leads to strange diff
            shader.u_view       = view_trans;
            shader.u_proj       = projection_trans;
            shader.u_projI      = projection_trans.invert();
            shader.u_vpI        = view_port.invert();
            
            for(int i = 0;i <model->nfaces();i ++ ){
                //for everyface, get vertex and draw lines
                vec3 screen_coords[3];
                for(int j = 0;j < 3;j ++ ){
                    screen_coords[j] = shader.vertex(i, j);
                }
                triangle(screen_coords, shader, zbuffer, image);
            }
        }else{
            for(int i = 0;i <model->nfaces();i ++ ){
                vec3 screen_coords[3];
                for(int j = 0;j < 3;j ++ ){
                    screen_coords[j] = trans_vec3(view_port*mvp, model->vert(i, j), 1.0);
                    screen_coords[j].x = (int)screen_coords[j].x;
                    screen_coords[j].y = (int)screen_coords[j].y;
                }
                for(int i = 0;i < 3;i ++ ){
                    line(screen_coords[i].x, screen_coords[i].y, screen_coords[(i+1)%3].x, screen_coords[(i+1)%3].y, image, black);
                }
            }
        }
        std::cout << std::endl;
    }

    double ratio = 50;
    vec3 o = trans_vec3(view_port*projection_trans*view_trans*get_model_trans(ratio, vec3(0, 0, 0)), vec3(0, 0, 0), 1.0);
    vec3 x = trans_vec3(view_port*projection_trans*view_trans*get_model_trans(ratio, vec3(0, 0, 0)), vec3(1, 0, 0), 1.0);
    vec3 y = trans_vec3(view_port*projection_trans*view_trans*get_model_trans(ratio, vec3(0, 0, 0)), vec3(0, 1, 0), 1.0);
    vec3 z = trans_vec3(view_port*projection_trans*view_trans*get_model_trans(ratio, vec3(0, 0, 0)), vec3(0, 0, 1), 1.0);
    line(int(o.x), int(o.y), int(x.x), int(x.y), image, red);
    line(int(o.x), int(o.y), int(y.x), int(y.y), image, green);
    line(int(o.x), int(o.y), int(z.x), int(z.y), image, blue);

    if(image.write_tga_file(output_dir + "output.tga"))
        std::cout << "file output to " << output_dir + "output.tga" << std::endl;
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
        std::cout << "zbuffer output to " << output_dir + "zbuffer.tga" << std::endl;
    else 
        std::cout << "can not output zbuffer" << std::endl;
    
    delete model;
    return 0;
}
