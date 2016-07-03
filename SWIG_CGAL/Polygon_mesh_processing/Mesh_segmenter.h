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
