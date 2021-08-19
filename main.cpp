#include "tgaimage.h"
#include "geometry.h"
#include "model.h"
#include "global.h"
#include "mygl.h"

const std::string shaders[4] = {"PhongShader", "NormalShader", "BlingPhongShader", "NoShader"};

const int height = 800;
const int width = 800;
const int depth = 255;

// vec3 light_dir(1, 0, 0);
// vec3 light_dir(-1, -1, -1);
vec3 light_dir(-1, 0, 0);
vec3 light_pos(300, 300, 0);
vec3 eye_pos(0, 0, 0);
vec3 center(300, 200, 0);
// vec3 center(100, 100, 270);
// vec3 center(-300, 0, 0);
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
TGAImage tex2;

TGAColor show_nv3(vec3 v3){
    return TGAColor(v3.x*255, v3.y*255, v3.z*255);
}

struct PhongShader : public IShader{
    public:
    mat<2, 3> v_uv;//written by vertex shader, read by fragment shader
    mat<3, 3> v_nrm;//fs will interpolate normal

    mat<4, 4> u_vp_mvp;
    mat<4, 4> u_model;

    vec2 u_t0;

    virtual vec3 vertex(int iface, int nthvert){
        v_uv.set_col(nthvert, model->uv(iface, nthvert));
        v_nrm.set_col(nthvert, proj<3>((u_model).invert_transpose()*embed<4>(model->normal(iface, nthvert), 0.0)));
        vec3 gl_vertex = model->vert(iface, nthvert);

        return trans_vec3(u_vp_mvp, gl_vertex, 1.0);
    }

    virtual bool fragment(vec3 bar, TGAColor &color){
        vec3 l, n; //light dir, normal
        vec2 b_uv = v_uv*bar;
        
        l = trans_vec3(u_model, light_dir.normalize(), 1.0).normalize(); 
        n = (v_nrm*bar).normalize();
        double diff = std::max<double>(-(l*n), 0.0);
       
        TGAColor c = tex0.get(u_t0.x*(b_uv.x), u_t0.y*(1-b_uv.y));
        color = c*diff;
        
        return false;   //discard
    }

};

mat<4, 4> getv(const vec3 eye_pos, const vec3 center, const vec3 up){
    vec3 z = (center - eye_pos).normalize(),
         x = cross(up, z).normalize(),
         y = cross(z, x).normalize();

    mat<4, 4> to_center = mat<4, 4>::identity(),
              base = mat<4, 4>::identity();

    for(int i = 0;i < 3;i ++ ){
        // to_center[i][3] = -center[i];
        base[0][i] = x[i];
        base[1][i] = y[i];
        base[2][i] = z[i];
    }
    
    return base*to_center;
}

struct NormalShader : public IShader{
    public:
    vec3 v_tri[3];
    mat<2, 3> v_uv;//written by vertex shader, read by fragment shader
    mat<3, 3> v_nrm;//fs will interpolate normal

    mat<4, 4> u_viewport;
    mat<4, 4> u_model;
    mat<4, 4> u_view;
    mat<4, 4> u_proj;
    mat<4, 4> u_model_it;

    vec2 u_t0;
    vec2 u_t1;
    vec2 u_t2;

    mat<3, 3> ndc_tri;

    virtual vec3 vertex(int iface, int nthvert){
        vec3 gl_vertex = trans_vec3(u_view*u_model, model->vert(iface, nthvert), 1.0);

        v_uv.set_col(nthvert, model->uv(iface, nthvert));
        v_nrm.set_col(nthvert, proj<3>((u_view*u_model).invert_transpose()*embed<4>(model->normal(iface, nthvert), 0.0)));
        ndc_tri.set_col(nthvert, gl_vertex);

        v_tri[nthvert] = trans_vec3(u_viewport*u_proj, gl_vertex, 1.0);
        return gl_vertex;
    }

    virtual bool fragment(vec3 bar, TGAColor &color){
        vec3 l, n, r; //light dir, normal, reflected light dir
        vec2 b_uv = v_uv*bar;
        // color = TGAColor(b_uv.x*255, b_uv.y*255, 0);
        // return false;

        TGAColor c = tex1.get(u_t1.x*(b_uv.x), u_t1.y*(1-b_uv.y));
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
        double diff = std::max<double>(nl, 0.0);

        c = tex0.get(u_t0.x*(b_uv.x), u_t0.y*(1-b_uv.y));
        color = c*diff;
        return false;   //discard
    }

};

struct BlingPhongShader : public IShader{
    public:
    vec3 v_tri[3];
    mat<2, 3> v_uv;//written by vertex shader, read by fragment shader
    mat<3, 3> v_nrm;//fs will interpolate normal

    mat<4, 4> u_viewport;
    mat<4, 4> u_model;
    mat<4, 4> u_view;
    mat<4, 4> u_proj;
    mat<4, 4> u_model_it;

    vec2 u_t0;
    vec2 u_t1;
    vec2 u_t2;

    mat<3, 3> ndc_tri;

    virtual vec3 vertex(int iface, int nthvert){
        vec3 gl_vertex = trans_vec3(u_view*u_model, model->vert(iface, nthvert), 1.0);

        v_uv.set_col(nthvert, model->uv(iface, nthvert));
        v_nrm.set_col(nthvert, proj<3>((u_view*u_model).invert_transpose()*embed<4>(model->normal(iface, nthvert), 0.0)));
        ndc_tri.set_col(nthvert, gl_vertex);

        v_tri[nthvert] = trans_vec3(u_viewport*u_proj, gl_vertex, 1.0);
        return gl_vertex;
    }

    virtual bool fragment(vec3 bar, TGAColor &color){
        vec3 l, n, r; //light dir, normal, reflected light dir
        vec2 b_uv = v_uv*bar;

        TGAColor c = tex1.get(u_t1.x*(b_uv.x), u_t1.y*(1-b_uv.y));
        n = vec3((double)c[2]/255*2-1, (double)c[1]/255*2-1, (double)c[0]/255*2-1).normalize();   //according to tinyredner/model.cpp, cant understand it
        // color = show_nv3(n);
        // return false;
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
        double spec = 1-pow(1-std::max(r.z, 0.0), (double)tex2.get(u_t2.x*(b_uv.x), u_t2.y*(1-b_uv.y))[0]);
        double ambi = 5;

        c = tex0.get(u_t0.x*(b_uv.x), u_t0.y*(1-b_uv.y));
        for (int i=0; i<3; i++) color[i] = std::min<double>(ambi + (double)c[i]*(diff + 0.76*spec), 255);
        // color = show_nv3(n);
        // color = white*diff;
       
        return false;   //discard
    }

};

void load(){
    // tex0.read_tga_file("../obj/diablo3_pose/_diffuse.tga");
    // tex1.read_tga_file("../obj/diablo3_pose/_nm_tangent.tga");
    // tex2.read_tga_file("../obj/diablo3_pose/_spec.tga");
    // model = new Model("../obj/diablo3_pose/object.obj");

    // tex0.read_tga_file("../obj/african_head/_diffuse.tga");
    // tex1.read_tga_file("../obj/african_head/normal.tga");
    // tex2.read_tga_file("../obj/african_head/_spec.tga");
    // model = new Model("../obj/african_head/object.obj");
    // tex1.read_tga_file("../obj/grid.tga");
    // tex0.read_tga_file("../obj/african_head/_diffuse.tga");
    // tex1.read_tga_file("../obj/african_head/normal.tga");
    // tex2.read_tga_file("../obj/african_head/_spec.tga");
    // model = new Model("../obj/african_head/object.obj");
    tex0.read_tga_file("../obj/floor_diffuse.tga");
    tex1.read_tga_file("../obj/floor_nm_tangent.tga");
    model = new Model("../obj/floor.obj");
}

int main(const int argc, const char** argv)
{
    TGAImage image(width, height, TGAImage::RGB);
    my_clear(image, white);
    load();

    double *zbuffer = new double[width*height];
    for(int i = 0;i < width*height;i ++ ) zbuffer[i] = -std::numeric_limits<double>::max();

    double  eye_fov = 45, 
            aspect_ratio = width/height,
            zNear = -1,
            zFar = -500,
            f1 = (zNear - zFar)/2.0,
            f2 = (zFar + zNear)/2.0;
 
    model_trans         = get_model_trans(model_trans);
    view_trans          = get_view(eye_pos, center, up);
    projection_trans    = get_projection(eye_fov, aspect_ratio, zNear, zFar);
    mvp                 = projection_trans*view_trans*model_trans;
    view_port           = get_viewport(0, 0, width, height, depth);

    /*
    write model
    */
    std::string shader_name = shaders[1];
    std::cout << shader_name << std::endl;
    
    if(shader_name == "BlingPhongShader"){
        BlingPhongShader shader;
        shader.u_viewport   = view_port;
        shader.u_model      = model_trans;   //to_center leads to strange diff
        shader.u_view       = view_trans;
        shader.u_proj       = projection_trans;
        shader.u_model_it   = (model_trans).invert_transpose();
        
        shader.u_t0 = vec2(tex0.get_width(), tex0.get_height());
        shader.u_t1 = vec2(tex1.get_width(), tex1.get_height());
        shader.u_t2 = vec2(tex2.get_width(), tex2.get_height());

        for(int i = 0;i <model->nfaces();i ++ ){
            //for everyface, get vertex and draw lines
            for(int j = 0;j < 3;j ++ ){
                shader.vertex(i, j);
            }
            triangle(shader.v_tri, shader, zbuffer, image);
        }
    }else if(shader_name == "NormalShader"){
        NormalShader shader;
        shader.u_viewport   = view_port;
        shader.u_model      = model_trans;   //to_center leads to strange diff
        shader.u_view       = view_trans;
        shader.u_proj       = projection_trans;
        shader.u_model_it   = (model_trans).invert_transpose();
        
        shader.u_t0 = vec2(tex0.get_width(), tex0.get_height());
        shader.u_t1 = vec2(tex1.get_width(), tex1.get_height());
        shader.u_t2 = vec2(tex2.get_width(), tex2.get_height());

        for(int i = 0;i <model->nfaces();i ++ ){
            //for everyface, get vertex and draw lines
            vec3 screen_coords[3];
            for(int j = 0;j < 3;j ++ ){
                shader.vertex(i, j);
            }
            triangle(shader.v_tri, shader, zbuffer, image);
        }
    }else if(shader_name == "PhongShader"){
        PhongShader shader;
        shader.u_vp_mvp = view_port*mvp;
        shader.u_model =(model_trans);   //to_center leads to strange diff
        shader.u_t0 = vec2(tex0.get_width(), tex0.get_height());

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

    vec3 origin = trans_vec3(view_port*mvp, vec3(0, 0, 0), 1.0);
    vec3 xaxis = trans_vec3(view_port*mvp, vec3(1, 0, 0), 1.0);
    vec3 yaxis = trans_vec3(view_port*mvp, vec3(0, 1, 0), 1.0);
    vec3 zaxis = trans_vec3(view_port*mvp, vec3(0, 0, 1), 1.0);
    line(int(origin.x), int(origin.y), int(xaxis.x), int(xaxis.y), image, red);
    line(int(origin.x), int(origin.y), int(yaxis.x), int(yaxis.y), image, green);
    line(int(origin.x), int(origin.y), int(zaxis.x), int(zaxis.y), image, blue);

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
