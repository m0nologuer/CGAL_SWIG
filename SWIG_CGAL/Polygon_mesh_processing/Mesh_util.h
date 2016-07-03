#ifndef SWIG_CGAL_MESH_UTIL_3_H
#define SWIG_CGAL_MESH_UTIL_3_H

#include <SWIG_CGAL/Common/Wrapper_iterator_helper.h>
#include <SWIG_CGAL/Kernel/Point_3.h>
#include <SWIG_CGAL/Kernel/Vector_3.h>
#include <SWIG_CGAL/Polyhedron_3/Polyhedron_3.h>

#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/bounding_box.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/algorithm.h>
#include <CGAL/Side_of_triangle_mesh.h>

//Create a bounding box for the polyhedron, aligned around the center of mass
template<class Polyhedron_wrapper>
struct MeshBoundingBox
{
    Point_3 min_point;
    Vector_3 x_axis;
    Vector_3 y_axis;
    Vector_3 z_axis;

    Point_3 get_origin(){return min_point;}
    Vector_3 get_x_axis(){return x_axis;}
    Vector_3 get_y_axis(){return y_axis;}
    Vector_3 get_z_axis(){return z_axis;}

    MeshBoundingBox(Vector_3 min, Vector_3 x, Vector_3 y, Vector_3 z){
        min_point = Point_3(min.x(), min.y(), min.z());
        x_axis = x;
        y_axis = y;
        z_axis = z;
    }

    MeshBoundingBox(Polyhedron_wrapper& poly){
        //get bounds
        Polyhedron_3_ mesh = poly.get_data();
        CGAL::Iso_cuboid_3<CGAL::Epick> c = CGAL::bounding_box(mesh.points_begin(),mesh.points_end());
        
        //orient correctly...
        min_point = c[0];
        x_axis = c[1] - c[0];
        y_axis = c[3] - c[0];
        z_axis = c[5] - c[0];
        
        double x_axis_sq_length = x_axis.squared_length ();
        double y_axis_sq_length = y_axis.squared_length ();
        double z_axis_sq_length = z_axis.squared_length ();

        Point_3 center_of_mass= Point_3(0,0,0);
        int vertex_count = 0;
        for (Polyhedron_3_::Vertex_iterator vert = mesh.vertices_begin(); 
                vert != mesh.vertices_end(); ++vert)
        {
            center_of_mass = Point_3(vert->point().x()+center_of_mass.x(),
                    vert->point().y() +center_of_mass.y(),
                    vert->point().z() + center_of_mass.z());
            vertex_count++;
        };
        Vector_3 shifted_com = center_of_mass - min_point;
        Vector_3 normalized_com = Vector_3((shifted_com *x_axis)/x_axis_sq_length,
                    (shifted_com *y_axis)/y_axis_sq_length,
                    (shifted_com *z_axis)/z_axis_sq_length);
        bool flip_x = normalized_com.x() > 0.5;
        bool flip_y = normalized_com.y() > 0.5;
        bool flip_z = normalized_com.z() > 0.5;

        //Re-align the box
        if (flip_x){
            min_point = min_point + x_axis;
            x_axis = -x_axis;
        };
        if (flip_y){
            min_point = min_point + y_axis;
            y_axis = -y_axis;
        };
        if (flip_z){
            min_point = min_point + z_axis;
            z_axis = -z_axis;
        };
    };
};


template<class Point>
struct Point_transform {

    Point_transform(Vector_3 t, Vector_3 x, Vector_3 y,
        Vector_3 z): translation(t), x_axis(x), y_axis(y), z_axis(z){};

    Point operator()( Point& p) { 
        Vector_3 p_vec = Vector_3(p.x(), p.y(), p.z())+ translation;
        p = Point(p_vec*x_axis, p_vec*y_axis, p_vec*z_axis);
        return p;
    };
    private:
        Vector_3 translation;
        Vector_3 x_axis;
        Vector_3 y_axis;
        Vector_3 z_axis;
};


template <class Polyhedron_wrapper>
class Polygon_mesh_util_wrapper
{
    typedef Polygon_mesh_util_wrapper<Polyhedron_wrapper> Self;
    //disable deep copy
    Self deepcopy();
    void deepcopy(const Self&);

public:
    #ifndef SWIG
    Polyhedron_wrapper& poly;
    #endif

    Polygon_mesh_util_wrapper(Polyhedron_wrapper& p): poly(p) {}

    Polyhedron_wrapper stitch_borders(){
        Polyhedron_3_ mesh = poly.get_data();
        CGAL::Polygon_mesh_processing::stitch_borders(mesh);
        Polyhedron_wrapper p_wrap(mesh);
        return mesh;
    }

    typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron_3_> Primitive;
    typedef CGAL::AABB_traits<CGAL::Epick, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;
    typedef CGAL::Side_of_triangle_mesh<Polyhedron_3_, CGAL::Epick> Point_inside;

    std::vector<std::vector<std::vector<float> > > voxelize(int dimension){
        // Construct AABB tree with a KdTree
        Polyhedron_3_ polyhedron = poly.get_data();
        Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
        tree.accelerate_distance_queries();
        // Initialize the point-in-polyhedron tester
        Point_inside inside_tester(tree);

        //Initialize 3D grif of booleans
        std::vector<std::vector<std::vector<float> > > voxel_grid = std::vector<std::vector<std::vector<float> > >();
        for (int i = 0; i < dimension; ++i){
            voxel_grid.push_back(std::vector<std::vector<float> >());
            for (int j = 0; j < dimension; ++j){
                voxel_grid[i].push_back(std::vector<float>());
                for (int k = 0; k < dimension; ++k)
                    voxel_grid[i][j].push_back(0.0);
            }
        }      

        MeshBoundingBox<Polyhedron_wrapper> box = MeshBoundingBox<Polyhedron_wrapper>(poly);

        for (int i = 0; i < dimension; ++i)
            for (int j = 0; j < dimension; ++j)
                for (int k = 0; k < dimension; ++k)
                {
                    float x = (float)i/(float)dimension;
                    float y = (float)j/(float)dimension;
                    float z = (float)k/(float)dimension;

                    Vector_3 offset = Vector_3(box.x_axis.x()*x + box.y_axis.x()*y + box.z_axis.x()*z,
                        box.x_axis.y()*x + box.y_axis.y()*y + box.z_axis.y()*z,
                        box.x_axis.z()*x + box.y_axis.z()*y + box.z_axis.z()*z);

                    CGAL::Point_3<CGAL::Epick> point = CGAL::Point_3<CGAL::Epick>(box.min_point.x() + offset.x(),
                        box.min_point.y() + offset.y(),box.min_point.z() + offset.z());

                    // Determine the side and return true if inside!
                    bool inside = (inside_tester(point) == CGAL::ON_BOUNDED_SIDE);
                    voxel_grid[i][j][k] = inside ? 1.0 : 0.0;
                }
        
        return voxel_grid;
    };

    Polyhedron_wrapper concatenate_mesh(std::vector<Polyhedron_wrapper> polylist){
        
        std::vector<Polyhedron_3_::Facet_const_handle> total_faces = std::vector<Polyhedron_3_::Facet_const_handle>();

        //Go through all the meshes
        for (typename std::vector<Polyhedron_wrapper>::iterator mesh = polylist.begin() ; mesh != polylist.end(); ++mesh)
        {
            //Collect all the faces
            for(Polyhedron_3_::Facet_iterator facet_it = mesh->get_data().facets_begin();
                facet_it != mesh->get_data().facets_end(); ++facet_it) {

                total_faces.push_back(facet_it);
            }
        };

        //Build faces into mesh
        Poly_Builder mesh_builder = Poly_Builder(total_faces);
        
        //wrap
        Polyhedron_wrapper p_wrap;
        p_wrap.get_data().delegate(mesh_builder);

        return p_wrap;
    };

    Polyhedron_wrapper transform_mesh(Vector_3 translation, Vector_3 x_axis,
        Vector_3 y_axis, Vector_3 z_axis){

        Point_transform<CGAL::Point_3<CGAL::Epick> > p_transform(translation, x_axis, y_axis, z_axis);

        Polyhedron_3_ mesh = poly.get_data();
        std::transform( mesh.points_begin(), mesh.points_end(), mesh.points_begin(), 
            p_transform) ;

        Polyhedron_wrapper p_wrap(mesh);
        return p_wrap;
    };
    
};
#endif //SWIG_CGAL_MESH_UTIL_3_H
