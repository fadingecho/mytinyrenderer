#include "tgaimage.h"
#include "geometry.h"
#include "model.h"
#include "global.h"
#include "mygl.h"
#include "string"

const int height = 800;
const int width = 800;
const int depth = 255;

vec3 light_pos(200, 350, 300);
vec3 eye_pos(0, 0, 0);
vec3 center(180, 150, 300);
vec3 up(0, 1, 0);
double *shadow_buffer;
TGAImage shadow_map;

const int model_cnt = 2;
const int texture_cnt = 3;
double model_scale[model_cnt] = {100, 200};
vec3 model_poi[model_cnt] = {
    vec3(0, 0, 0),
    vec3(-100, 50, -200)};
std::string model_path[model_cnt] = {
    // "../obj/african_head/object.obj",
    "../obj/diablo3_pose/object.obj",
    "../obj/floor.obj"};
std::string texture_path[model_cnt][texture_cnt] = {
    // {"../obj/african_head/_diffuse.tga",
    //  "../obj/african_head/normal.tga",
    //  "../obj/african_head/_spec.tga"},
    {"../obj/diablo3_pose/_diffuse.tga",
     "../obj/diablo3_pose/_nm_tangent.tga",
     "../obj/diablo3_pose/_spec.tga"},
    {"../obj/floor_diffuse.tga",
     "../obj/floor_nm_tangent.tga",
     ""}};

mat<4, 4> model_trans;
mat<4, 4> view_trans;
mat<4, 4> projection_trans;
mat<4, 4> mvp;
mat<4, 4> view_port;

Model *model = NULL, *models[model_cnt];
const std::string output_dir = "../output/";
TGAImage *tex, textures[model_cnt][texture_cnt];

TGAColor show_nv3(vec3 v3)
{
    double s = v3.norm();
    return TGAColor(v3.x * 255 / s, v3.y * 255 / s, v3.z * 255 / s);
}

void init_buffer(double *buffer, const int w, const int h)
{
    for (int i = 0; i < w * h; i++)
    {
        buffer[i] = -std::numeric_limits<double>::max();
    }
}

struct DepthShader : public IShader
{
public:
    mat<3, 3> v_tri;

    virtual vec3 vertex(int iface, int nthvert)
    {
        vec3 gl_vertex = trans_vec3(u_viewport * u_mvp, model->vert(iface, nthvert), 1.0);
        v_tri.set_col(nthvert, gl_vertex);
        return gl_vertex;
    }

    virtual bool fragment(vec3 bar, TGAColor &color)
    {
        return false; //discard
    }
};

struct BlingPhongShader : public IShader
{
public:
    mat<4, 4> u_shadow;
    mat<2, 3> v_uv;  //written by vertex shader, read by fragment shader
    mat<3, 3> v_nrm; //fs will interpolate normal
    mat<3, 3> v_tri;

    virtual vec3 vertex(int iface, int nthvert)
    {
        vec3 gl_vertex = trans_vec3(u_view * u_model, model->vert(iface, nthvert), 1.0);

        v_uv.set_col(nthvert, model->uv(iface, nthvert));
        v_nrm.set_col(nthvert, proj<3>(u_camera_normal_trans * embed<4>(model->normal(iface, nthvert), 0.0)));
        v_tri.set_col(nthvert, gl_vertex);
        gl_vertex = trans_vec3(u_viewport * u_proj, gl_vertex, 1.0);
        return gl_vertex;
    }

    virtual bool fragment(vec3 bar, TGAColor &color)
    {
        double ambi, diff, spec;
        vec3 l, n, r; //light dir, normal, reflected light dir
        vec2 b_uv = v_uv * bar;

        vec3 bc_camera_pt = v_tri * bar;
        vec3 bc_model_pt = trans_vec3(u_to_model, bc_camera_pt, 1.0);
        vec3 b_shadow = trans_vec3(u_shadow, bc_model_pt, 1.0);

        double shadow_z = (double)shadow_map.get(b_shadow.x, b_shadow.y)[0] / 255;
        double shadow = 0.3 + 0.7 * (shadow_z < b_shadow.z + 0.01);

        TGAColor c = tex[1].get(tex[1].get_width() * (b_uv.x), tex[1].get_height() * (1 - b_uv.y));
        n = vec3((double)c[2] / 255 * 2 - 1, (double)c[1] / 255 * 2 - 1, (double)c[0] / 255 * 2 - 1).normalize(); //according to tinyredner/model.cpp, cant understand it
        {                                                                                                         //TBN
            vec3 b_n = (v_nrm * bar).normalize();
            mat<3, 3> A;
            A[0] = v_tri.col(1) - v_tri.col(0);
            A[1] = v_tri.col(2) - v_tri.col(0);
            A[2] = b_n; //bn

            mat<3, 3> AI = A.invert();
            vec3 i = AI * vec3(v_uv[0][1] - v_uv[0][0], v_uv[0][2] - v_uv[0][0], 0);
            vec3 j = AI * vec3(v_uv[1][1] - v_uv[1][0], v_uv[1][2] - v_uv[1][0], 0);

            mat<3, 3> B;
            B.set_col(0, i.normalize());
            B.set_col(1, j.normalize());
            B.set_col(2, b_n);

            n = (B * n).normalize();
        }

        l = trans_vec3(u_view * u_model, (v_tri * bar - light_pos).normalize(), 0.0).normalize();
        double nl = -(l * n);

        r = (n * (nl * 2.0) + l).normalize();

        ambi = 5;
        diff = std::max<double>(nl, 0.0);
        spec = 1 - pow(1 - std::max(r.z, 0.0), (double)tex[2].get(tex[2].get_width() * (b_uv.x), tex[2].get_height() * (1 - b_uv.y))[0]);

        c = tex[0].get(tex[0].get_width() * (b_uv.x), tex[0].get_height() * (1 - b_uv.y));
        for (int i = 0; i < 3; i++)
            color[i] = std::min<double>(ambi + (double)c[i] * (diff + spec) * shadow, 255);

        return false; //discard
    }
};

void load()
{
    for (int m = 0; m < model_cnt; m++)
    {
        for (int t = 0; t < texture_cnt; t++)
        {
            if (texture_path[m][t].empty())
            {
                std::cout << "texture#" << t << " is empty" << std::endl;
                textures[m][t] = TGAImage(width, height, 0);
                textures[m][t].clear();
                continue;
            }
            if (!textures[m][t].read_tga_file(texture_path[m][t]))
                std::cout << "cant read texture#" << t << "in" << texture_path[m][t] << std::endl;
        }
        models[m] = new Model(model_path[m].c_str());
        std::cout << std::endl;
    }
}

int main(const int argc, const char **argv)
{
    MyTimer timer;
    TGAImage frame(width, height, TGAImage::RGB);
    shadow_map = TGAImage(width / 2, height / 2, 1);
    my_clear(frame, black);
    my_clear(shadow_map, black);
    double *zbuffer = new double[frame.get_width() * frame.get_height()];
    shadow_buffer = new double[shadow_map.get_width() * shadow_map.get_height()];
    init_buffer(zbuffer, frame.get_width(), frame.get_height());
    init_buffer(shadow_buffer, shadow_map.get_width(), shadow_map.get_height());

    double eye_fov = 45,
           aspect_ratio = width / height,
           zNear = -1,
           zFar = -1100;

    mat<4, 4> shadow_trans;
    load();

    { //shadow mapping
        view_trans = get_view(eye_pos, light_pos, up);
        projection_trans = get_projection(eye_fov, aspect_ratio, zNear, zFar);
        view_port = get_viewport(0, 0, shadow_map.get_width(), shadow_map.get_height(), depth);

        for (uint8_t cnt = 0; cnt < model_cnt; cnt++)
        {
            model = models[cnt];
            tex = textures[cnt];
            
            model_trans = get_model_trans(model_scale[cnt], model_poi[cnt]);
            shadow_trans = view_port * projection_trans * view_trans;

            DepthShader shader;
            shader.u_viewport = view_port;
            shader.u_model = model_trans; //to_center leads to strange diff
            shader.u_view = view_trans;
            shader.u_proj = projection_trans;
            shader.u_to_camera = projection_trans.invert() * view_port.invert();
            shader.u_mvp = projection_trans * view_trans * model_trans;

            for (int i = 0; i < model->nfaces(); i++)
            {
                //for everyface, get vertex and draw lines
                vec3 screen_coords[3];
                for (int j = 0; j < 3; j++)
                {
                    screen_coords[j] = shader.vertex(i, j);
                }
                triangle(screen_coords, shader, shadow_buffer, shadow_map);
            }
        }

        int w = shadow_map.get_width(), h = shadow_map.get_height();
        for (int i = 0; i < w; i++)
            for (int j = 0; j < h; j++)
            {
                double ratio = shadow_buffer[i + j * w];
                shadow_map.set(i, j, TGAColor(255*ratio));
            }
    }

    timer.start();
    { //rendering frame
        view_trans = get_view(eye_pos, center, up);
        projection_trans = get_projection(eye_fov, aspect_ratio, zNear, zFar);
        view_port = get_viewport(0, 0, frame.get_width(), frame.get_height(), depth);

        for (uint8_t cnt = 0; cnt < 2; cnt++)
        {
            model = models[cnt];
            tex = textures[cnt];
            
            model_trans = get_model_trans(model_scale[cnt], model_poi[cnt]);
            mvp = projection_trans * view_trans * model_trans;

            BlingPhongShader shader;
            shader.u_viewport = view_port;
            shader.u_model = model_trans; //to_center leads to strange diff
            shader.u_view = view_trans;
            shader.u_proj = projection_trans;
            shader.u_to_camera = projection_trans.invert() * view_port.invert();
            shader.u_mvp = projection_trans * view_trans * model_trans;
            shader.u_shadow = shadow_trans * model_trans;
            shader.u_to_model = model_trans.invert() * view_trans.invert();
            shader.u_camera_normal_trans = (view_trans * model_trans).invert_transpose();
            for (int i = 0; i < model->nfaces(); i++)
            {
                //for everyface, get vertex and draw lines
                vec3 screen_coords[3];
                for (int j = 0; j < 3; j++)
                {
                    screen_coords[j] = shader.vertex(i, j);
                }
                triangle(screen_coords, shader, zbuffer, frame);
            }
        }
        timer.end();
        if (frame.write_tga_file(output_dir + "output.tga"))
            std::cout << "file output to " << output_dir + "output.tga" << std::endl;
        else
            std::cout << "can not output file" << std::endl;

        // int w = frame.get_width(), h = frame.get_height();
        // TGAImage buffer(w, h, TGAImage::RGBA);
        // for (int i = 0; i < w; i++)
        //     for (int j = 0; j < h; j++)
        //     {
        //         double ratio = zbuffer[i + j * w];
        //         buffer.set(i, j, TGAColor(255 * ratio, 255 * ratio, 255 * ratio, 255));
        //     }
        // if (buffer.write_tga_file(output_dir + "zbuffer.tga"))
        //     std::cout << "zbuffer output to " << output_dir + "zbuffer.tga" << std::endl;
        // else
        //     std::cout << "can not output zbuffer" << std::endl;
    }
    delete model;
    return 0;
}

// double ratio = 50;
//     vec3 o = trans_vec3(view_port*projection_trans*view_trans*get_model_trans(ratio, vec3(0, 0, 0)), vec3(0, 0, 0), 1.0);
//     vec3 x = trans_vec3(view_port*projection_trans*view_trans*get_model_trans(ratio, vec3(0, 0, 0)), vec3(1, 0, 0), 1.0);
//     vec3 y = trans_vec3(view_port*projection_trans*view_trans*get_model_trans(ratio, vec3(0, 0, 0)), vec3(0, 1, 0), 1.0);
//     vec3 z = trans_vec3(view_port*projection_trans*view_trans*get_model_trans(ratio, vec3(0, 0, 0)), vec3(0, 0, 1), 1.0);
//     line(int(o.x), int(o.y), int(x.x), int(x.y), frame, red);
//     line(int(o.x), int(o.y), int(y.x), int(y.y), frame, green);
//     line(int(o.x), int(o.y), int(z.x), int(z.y), frame, blue);
