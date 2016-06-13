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
          facet_it != mesh.facets_end(); ++facet_it, ++facet_id) {

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
