// ------------------------------------------------------------------------------
// Copyright (c) 2013 GeometryFactory (FRANCE)
// Distributed under the Boost Software License, Version 1.0. (See accompany-
// ing file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
// ------------------------------------------------------------------------------

#ifndef SWIG_CGAL_MESH_SEGMENTATION_3_H
#define SWIG_CGAL_MESH_SEGMENTATION_3_H

#include <SWIG_CGAL/Common/Wrapper_iterator_helper.h>
#include <SWIG_CGAL/Kernel/Point_3.h>
#include <SWIG_CGAL/Kernel/Vector_3.h>
#include <SWIG_CGAL/Polyhedron_3/Polyhedron_3.h>

#include <CGAL/mesh_segmentation.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/bounding_box.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/algorithm.h>
#include <CGAL/Side_of_triangle_mesh.h>


// Property map associating a facet with an integer as id to an
// element in a vector stored internally
template<class ValueType>
struct Facet_with_id_pmap
    : public boost::put_get_helper<ValueType&,
             Facet_with_id_pmap<ValueType> >
{
    typedef Polyhedron_3_::Facet_const_handle key_type;
    typedef ValueType value_type;
    typedef value_type& reference;
    typedef boost::lvalue_property_map_tag category;
    Facet_with_id_pmap(
      std::vector<ValueType>& internal_vector
    ) : internal_vector(internal_vector) { }
    reference operator[](key_type key) const
    { return internal_vector[key->id()]; }
private:
    std::vector<ValueType>& internal_vector;
};

typedef CGAL::Polyhedron_incremental_builder_3<Polyhedron_3_::HalfedgeDS> Mesh_Builder;
template <class HDS>
class Polyhedron_from_facelist_builder : public CGAL::Modifier_base<HDS> {
    std::vector<Polyhedron_3_::Facet_const_handle> faces;
    public:
        Polyhedron_from_facelist_builder(std::vector<Polyhedron_3_::Facet_const_handle> f) :faces(f){}
        
        void operator()( HDS& hds) {
            Mesh_Builder mesh_builder = Mesh_Builder(hds);
            int facet_count = faces.size();
            int total_vertices = 0;
            for (int j = 0; j < facet_count; j++)
                total_vertices += faces[j]->facet_degree(); 

            int vertex_counter = 0;
            mesh_builder.begin_surface( total_vertices, facet_count);
            for (int j = 0; j < facet_count; j++){

                Polyhedron_3_::Facet_const_handle face = faces[j];
                Polyhedron_3_::Halfedge_around_facet_const_circulator halfedge = face->facet_begin();
                halfedge++;
                
                //add rest of vertices
                int degree = 0;
                do{
                    mesh_builder.add_vertex(halfedge->vertex ()->point());
                    degree++;
                }while(halfedge++ != face->facet_begin());
                
                //add vertices to face
                mesh_builder.begin_facet();
                for (int k = 0; k < degree; k++)
                    mesh_builder.add_vertex_to_facet( vertex_counter + k);
                mesh_builder.end_facet();
                vertex_counter = vertex_counter+degree;
            }
            mesh_builder.end_surface();

        }
    };
typedef Polyhedron_from_facelist_builder<Polyhedron_3_::HalfedgeDS> Poly_Builder;

//Create similarity signature for a polyhedron, used for computing distance
template<class Polyhedron_wrapper>
struct MeshSignature
    {
        double weights[15625];
        Point_3 points[15625];

        double weight_at(int i){
            return weights[i];
        };
        Point_3 point_at(int i){
            return points[i];
        };

        MeshSignature(Polyhedron_wrapper& poly){

            //PCA bounding box orientation
            Polyhedron_3_ mesh = poly.get_data();
            CGAL::Iso_cuboid_3<CGAL::Epick> c = CGAL::bounding_box(mesh.points_begin(),mesh.points_end());
            
            Point_3 min_point = c[0];
            Vector_3 x_axis = c[1] - c[0];
            Vector_3 y_axis = c[3] - c[0];
            Vector_3 z_axis = c[5] - c[0];
            double x_axis_sq_length = x_axis.squared_length ();
            double y_axis_sq_length = y_axis.squared_length ();
            double z_axis_sq_length = z_axis.squared_length ();

            //Orient around center of mass
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

            //initialize signature to zero
            for (int i = 0; i < 15625; ++i){
                weights[i] = 0;
                points[i] = Point_3(0,0,0);
            }

            //Iterate through vertices
            for (Polyhedron_3_::Vertex_iterator vert = mesh.vertices_begin(); 
                vert != mesh.vertices_end(); ++vert)
            {
                Point_3 p = vert->point();
                Vector_3 direction = p - min_point;

                Point_3 scaled_pos = Point_3((direction *x_axis)/x_axis_sq_length,
                    (direction *y_axis)/y_axis_sq_length,
                    (direction *z_axis)/z_axis_sq_length);

                if (flip_x) {scaled_pos = Point_3(1-scaled_pos.x(), scaled_pos.y(), scaled_pos.z()); }
                if (flip_y) {scaled_pos = Point_3(scaled_pos.x(), 1-scaled_pos.y(), scaled_pos.z()); }
                if (flip_z) {scaled_pos = Point_3(scaled_pos.x(), scaled_pos.y(), 1-scaled_pos.z()); }

                //Place in box
                int x_box = (int)(fmin(scaled_pos.x()*25,24));
                int y_box = (int)(fmin(scaled_pos.y()*25,24));
                int z_box = (int)(fmin(scaled_pos.z()*25,24));

                int box_index = x_box*25*25+y_box*25+z_box;

                //Calculate the area of every incident face
                double vertex_angle = 0;
                double vertex_area = 0;
                    
                //Circulate around vertex
                Polyhedron_3_::Halfedge_around_vertex_circulator halfedge_it = vert->vertex_begin(); 
                do
                {  
                    //calculate halfedge angle
                    Point_3 point1 = halfedge_it->opposite()->vertex()->point();
                    Point_3 point2 = halfedge_it->next_on_vertex()->opposite()->vertex()->point();
                    Vector_3 vec1 = p - point1;
                    Vector_3 vec2 = p - point2;
                    double dot_product = vec1*vec2/sqrt(vec1.squared_length()*vec2.squared_length());
                    vertex_angle += acos(dot_product);
                    
                    //calculate facet area
                    double facet_area = 0;
                    int degree = 3;

                    Vector_3 dir1 = vec1;
                    Vector_3 dir2 = vec2;
                    Vector_3 cross_product = Vector_3(dir1.y()*dir2.z()-dir2.y()*dir1.z(),
                        dir1.z()*dir2.x()-dir2.z()*dir1.x(),dir1.x()*dir2.y()-dir2.x()*dir1.y());
                    facet_area += sqrt(cross_product.squared_length())*0.5;
                    
                    Point_3 point3;
                    Polyhedron_3_::Halfedge_around_vertex_circulator halfedge_face_it = halfedge_it->vertex_begin();
                    //sum area via triangles
                    for ( ; ++halfedge_face_it != halfedge_it->vertex_begin();){
                        point1 = point2;
                        point2 = point3;
                        point3 = halfedge_face_it->vertex()->point();
                        dir1 = point1 -point2;
                        dir2 = point1 -point3;
                        cross_product = Vector_3(dir1.y()*dir2.z()-dir2.y()*dir1.z(),
                            dir1.z()*dir2.x()-dir2.z()*dir1.x(),dir1.x()*dir2.y()-dir2.x()*dir1.y());

                        facet_area += sqrt(cross_product.squared_length())*0.5;
                        degree++;
                    }
                    
                    vertex_area += facet_area/float(degree); 

                }while(++halfedge_it != vert->vertex_begin());
                
                //update signature by most curved point
                double curvature_measure = (6.28-vertex_angle)/vertex_area;
                double measure = 1/(1+curvature_measure);

                if (measure > weights[box_index]){
                    weights[box_index] = measure;
                    points[box_index] = scaled_pos;
                }

            }
        }
    };

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
class Polygon_mesh_segmentation_wrapper
{
    typedef Polygon_mesh_segmentation_wrapper<Polyhedron_wrapper> Self;
    //disable deep copy
    Self deepcopy();
    void deepcopy(const Self&);

public:
    #ifndef SWIG
    Polyhedron_wrapper& poly;
    #endif

    Polygon_mesh_segmentation_wrapper(Polyhedron_wrapper& p): poly(p) {}

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
    }

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
        }

        //Build faces into mesh
        Poly_Builder mesh_builder = Poly_Builder(total_faces);
        
        //wrap
        Polyhedron_wrapper p_wrap;
        p_wrap.get_data().delegate(mesh_builder);

        return p_wrap;
    }

    Polyhedron_wrapper transform_mesh(Vector_3 translation, Vector_3 x_axis,
        Vector_3 y_axis, Vector_3 z_axis){

        Point_transform<CGAL::Point_3<CGAL::Epick> > p_transform(translation, x_axis, y_axis, z_axis);

        Polyhedron_3_ mesh = poly.get_data();
        std::transform( mesh.points_begin(), mesh.points_end(), mesh.points_begin(), 
            p_transform) ;

        Polyhedron_wrapper p_wrap(mesh);
        return p_wrap;
    }
    
    std::vector< std::vector<Point_3> > get_border_loops(){

        Polyhedron_3_ mesh = poly.get_data();
        std::vector< std::vector<Point_3> > border_loops = std::vector< std::vector<Point_3> >();
        std::vector<Polyhedron_3_::Halfedge_handle> used_halfedges = std::vector<Polyhedron_3_::Halfedge_handle>();

        // loop through all border halfedges
        for(Polyhedron_3_::Halfedge_iterator he_it = mesh.border_halfedges_begin();
          he_it != mesh.halfedges_end(); ++he_it) {

            //if we haven't found this already
            if (std::find(used_halfedges.begin(), used_halfedges.end(), he_it) == used_halfedges.end()) {
                //for each halfedge, extract a halfedge loop
                std::vector<Point_3> point_loop = std::vector<Point_3>();
                for (Polyhedron_3_::Halfedge_handle loop = he_it->next(); loop != he_it; loop = loop->next()){
                    point_loop.push_back(loop->vertex()->point());
                    used_halfedges.push_back(loop);
                }

                border_loops.push_back(point_loop);
            }
        }
        return border_loops;
    }

    std::vector< std::vector<Polyhedron_3_::Halfedge_handle> > get_cuts(){

        Polyhedron_3_ mesh = poly.get_data();
        boost::unordered_map<Polyhedron_3_::Facet_handle, bool> marked_faces = boost::unordered_map<Polyhedron_3_::Facet_handle, bool>();
        boost::unordered_map<Polyhedron_3_::Vertex_handle, bool> marked_vertices = boost::unordered_map<Polyhedron_3_::Vertex_handle, bool>();
        boost::unordered_map<Polyhedron_3_::Halfedge_handle,bool> marked_edges = boost::unordered_map<Polyhedron_3_::Halfedge_handle,bool> ();

        marked_faces.insert(std::make_pair(mesh.facets_begin(), true));
        bool keep_cutting = true;

        //while there remains an edge e adjacent to only one triangle t Remove e and t
        while(keep_cutting){

            keep_cutting = false; //assume we're not going to find anything
            
            //search through halfedges (except borders)
            for(Polyhedron_3_::Halfedge_iterator he_it = mesh.halfedges_begin();
              he_it != mesh.border_halfedges_begin(); ++he_it) {

                //if this edge is not marked
                if(marked_edges.find(he_it) == marked_edges.end()){

                    bool incident_facet_marked = marked_faces.find(he_it->facet()) != marked_faces.end();
                    bool opposite_facet_marked = marked_faces.find(he_it->opposite()->facet()) != marked_faces.end();

                    //If this is the kind of edge we want
                    if (incident_facet_marked != opposite_facet_marked){

                        //Mark this edge pair
                        marked_edges.insert(std::make_pair(he_it,true));
                        marked_edges.insert(std::make_pair(he_it->opposite(),true));

                        //Mark the correct face as well
                        if (!incident_facet_marked)
                            marked_faces.insert(std::make_pair(he_it->facet(), true));
                        else
                            marked_faces.insert(std::make_pair(he_it->opposite()->facet(), true));
                    
                        keep_cutting = true;
                    }
                }
            }
        }

        //while there remains a vertex v adjacent to only one edge e Remove v and e
        while(keep_cutting){

            keep_cutting = false; //assume we're not going to find anything
            
            //search through vertices
            for(Polyhedron_3_::Vertex_iterator v_it = mesh.vertices_begin();
              v_it != mesh.vertices_end(); ++v_it) {

                //if this vertex is not marked
                if(marked_vertices.find(v_it) == marked_vertices.end()){

                    //check how many non-marked edges are connected
                    int non_marked_edge_count = 0;
                    Polyhedron_3_::Halfedge_handle lone_edge;

                    for(Polyhedron_3_::Halfedge_handle he_it = v_it->vertex_begin();
                        ++he_it != v_it->vertex_begin();){
                        //if this edge is not marked
                        if(marked_edges.find(he_it) != marked_edges.end()){
                            non_marked_edge_count++;
                            lone_edge= he_it;
                        }
                    } 

                    //If this is a lone vertex
                    if (non_marked_edge_count == 1){

                        //Mark this edge pair
                        marked_edges.insert(std::make_pair(lone_edge,true));
                        marked_edges.insert(std::make_pair(lone_edge->opposite(),true));

                        //Mark this vertex
                        marked_vertices.insert(std::make_pair(lone_edge,true));

                        keep_cutting = true;
                    }
                }
            }
        }

        //Pull out the cuts
        std::vector< std::vector<Polyhedron_3_::Halfedge_handle> > cuts = std::vector< std::vector<Polyhedron_3_::Halfedge_handle> >();
        for(Polyhedron_3_::Vertex_iterator v_it = mesh.vertices_begin();
              v_it != mesh.vertices_end(); ++v_it) {

            //if this vertex is not marked
            if(marked_vertices.find(v_it) == marked_vertices.end()){

                //start a new loop
                cuts.push_back(std::vector<Polyhedron_3_::Halfedge_handle>());
                bool finished_loop = false;
                Polyhedron_3_::Vertex_handle loop_vertex = v_it;

                while (!finished_loop){
                    //look for the direction to go in
                    for(Polyhedron_3_::Halfedge_handle he_it = v_it->vertex_begin();
                            ++he_it != v_it->vertex_begin();){
                       
                        //if this edge is not marked
                        if(marked_edges.find(he_it) == marked_edges.end()){
                           loop_vertex = he_it->opposite()->vertex();
                           //if the vertex is not marked
                            if (marked_vertices.find(loop_vertex)  == marked_vertices.end()){
                                cuts.back().push_back(he_it);

                                //mark the vertex
                                marked_vertices.insert(std::make_pair(loop_vertex,true));
                            }
                            else{
                                finished_loop = true;
                            }

                        }
                    }    
                }
            }
        }

        return cuts;
    }
    
    std::vector<Polyhedron_wrapper> segmentation()
    {
        Polyhedron_3_ mesh = poly.get_data();

        // assign id field for each facet
        std::size_t facet_id = 0;
        for(Polyhedron_3_::Facet_iterator facet_it = mesh.facets_begin();
          facet_it != mesh.facets_end(); ++facet_it, ++facet_id) {
            facet_it->id() = facet_id;
        }

        // create a property-map for SDF values
        std::vector<double> sdf_values(mesh.size_of_facets());
        Facet_with_id_pmap<double> sdf_property_map(sdf_values);
        CGAL::sdf_values(mesh, sdf_property_map);

        // create a property-map for segment-ids
        std::vector<std::size_t> segment_ids(mesh.size_of_facets());
        Facet_with_id_pmap<std::size_t> segment_property_map(segment_ids);
        size_t segments = CGAL::segmentation_from_sdf_values(mesh, sdf_property_map, segment_property_map);
      
        //Create list of lists to store faces for each segment
        std::vector<std::vector<Polyhedron_3_::Facet_const_handle> > polyhedron_facets = std::vector<std::vector<Polyhedron_3_::Facet_const_handle> > ();
        for (int i = 0; i < segments; ++i)
        {
           std::vector<Polyhedron_3_::Facet_const_handle> new_vector = std::vector<Polyhedron_3_::Facet_const_handle>();
            polyhedron_facets.push_back(new_vector);
        }
        
        //Split the mesh into seperate lists of faces per segment
        for(Polyhedron_3_::Facet_iterator facet_it = mesh.facets_begin();
          facet_it != mesh.facets_end(); ++facet_it) {

            size_t index = segment_property_map[facet_it];
            polyhedron_facets[index].push_back(facet_it);
        }

        //Combine face lists into Polyhedra
        std::vector<Polyhedron_wrapper> output_polyhedron_list = std::vector<Polyhedron_wrapper>();
        for (int i = 0; i < segments; ++i)
        {   

            Poly_Builder mesh_builder = Poly_Builder(polyhedron_facets[i]);
            Polyhedron_wrapper p_wrap;
            p_wrap.get_data().delegate(mesh_builder);

            output_polyhedron_list.push_back(p_wrap);

        }

        return output_polyhedron_list;
    }
};


#endif //SWIG_CGAL_POINT_SET_PROCESSING_3_H
