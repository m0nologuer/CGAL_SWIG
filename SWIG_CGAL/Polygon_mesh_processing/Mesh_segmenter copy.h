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

#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/bounding_box.h>

//Create similarity signature for a polyhedron, used for computing distance
template<class Polyhedron_wrapper>
struct MeshSignature
    {
        double weights[1000];
        Point_3 points[1000];

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
            for (int i = 0; i < 1000; ++i){
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
                int x_box = (int)(fmin(scaled_pos.x()*10,9));
                int y_box = (int)(fmin(scaled_pos.y()*10,9));
                int z_box = (int)(fmin(scaled_pos.z()*10,9));

                int box_index = x_box*10*10+y_box*10+z_box;

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
    Point_3 get_x_axis(){return x_axis;}
    Point_3 get_y_axis(){return y_axis;}
    Point_3 get_z_axis(){return z_axis;}

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
            x_axis = -1 * x_axis;
        }
        if (flip_y){
            min_point = min_point + y_axis;
            y_axis = -1 * y_axis;
        }
        if (flip_z){
            min_point = min_point + z_axis;
            z_axis = -1 * z_axis;
        }
    }
}


template <class Polyhedron_wrapper>
class Poly_mesh_processor
{
    typedef Poly_mesh_processor<Polyhedron_wrapper> Self;
    //disable deep copy
    Self deepcopy();
    void deepcopy(const Self&);

public:
    #ifndef SWIG
    Polyhedron_wrapper& poly;
    #endif

    Poly_mesh_processor(Polyhedron_wrapper& p): poly(p) {}

    Polyhedron_wrapper stitch_borders(){
        Polyhedron_3_ mesh = poly.get_data();
        CGAL::Polygon_mesh_processing::stitch_borders(mesh);
        Polyhedron_wrapper p_wrap(mesh);
        return mesh;
    }
    
};


#endif //SWIG_CGAL_POINT_SET_PROCESSING_3_H
