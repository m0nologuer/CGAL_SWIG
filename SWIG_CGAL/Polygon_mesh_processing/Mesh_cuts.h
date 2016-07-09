#ifndef SWIG_CGAL_MESH_CUTS_3_H
#define SWIG_CGAL_MESH_CUTS_3_H

#include <SWIG_CGAL/Common/Wrapper_iterator_helper.h>
#include <SWIG_CGAL/Kernel/Point_3.h>
#include <SWIG_CGAL/Kernel/Vector_3.h>
#include <SWIG_CGAL/Polyhedron_3/Polyhedron_3.h>

#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/parameterize.h>
#include <CGAL/Parameterization_mesh_patch_3.h>

#include <CGAL/Polygon_2_algorithms.h>



typedef CGAL::Parameterization_polyhedron_adaptor_3<Polyhedron_3_> Parameterization_polyhedron_adaptor;
typedef CGAL::Parameterization_mesh_patch_3<Parameterization_polyhedron_adaptor> Mesh_patch_polyhedron;

typedef CGAL::Parameterizer_traits_3<Mesh_patch_polyhedron> Parameterizer;
typedef boost::unordered_map<Polyhedron_3_::Vertex_handle,std::pair<float, float> > Parameterization;

Polyhedron_3_::Halfedge_handle prev_halfedge;
bool sorter(Polyhedron_3_::Halfedge_handle const & a, Polyhedron_3_::Halfedge_handle const & b) {
                Vector_3 dist_a = a->vertex()->point() - prev_halfedge->vertex()->point();
                Vector_3 dist_b = b->vertex()->point() - prev_halfedge->vertex()->point();
                return dist_a.squared_length() < dist_b.squared_length();
            }


template <class Polyhedron_wrapper>
class Polygon_mesh_cuts_wrapper
{
    typedef Polygon_mesh_cuts_wrapper<Polyhedron_wrapper> Self;
    //disable deep copy
    Self deepcopy();
    void deepcopy(const Self&);

public:
    #ifndef SWIG
    Polyhedron_wrapper& poly;
    Polyhedron_3_ mesh;
    std::vector<Polyhedron_3_::Vertex_handle> seam;
    std::vector< std::vector<Polyhedron_3_::Halfedge_handle> > sub_paths;
    Parameterization parameterization;
    #endif

    Polygon_mesh_cuts_wrapper(Polyhedron_wrapper& p): poly(p) { 
        mesh = p.get_data();
        parameterization = Parameterization();
    }
    
    std::vector< std::vector<Point_3> > get_border_loops(){
        
        std::vector< std::vector<Point_3> > border_loops = std::vector< std::vector<Point_3> >();
        std::vector<Polyhedron_3_::Halfedge_handle> used_halfedges = std::vector<Polyhedron_3_::Halfedge_handle>();
        
        // loop through all border halfedges
        for(Polyhedron_3_::Halfedge_iterator he_it = mesh.halfedges_begin();
          he_it != mesh.halfedges_end(); he_it++) {

            //if we haven't found this already
            if (he_it->is_border_edge() && std::find(used_halfedges.begin(), used_halfedges.end(), he_it) == used_halfedges.end()) {
                //for each halfedge, extract a halfedge loop
                std::vector<Point_3> point_loop = std::vector<Point_3>();
                for (Polyhedron_3_::Halfedge_handle loop = he_it->next(); loop != he_it; loop = loop->next()){
                    point_loop.push_back(loop->vertex()->point());
                    used_halfedges.push_back(loop);
                };

                border_loops.push_back(point_loop);
            }
        }; 
        return border_loops;
        
    };

    void add_path_between_two_vertices(Polyhedron_3_::Halfedge_handle start, Polyhedron_3_::Halfedge_handle end){

        //can't have halfedges used twice
        boost::unordered_map<Polyhedron_3_::Halfedge_handle,bool> marked_edges;
        std::vector <Polyhedron_3_::Halfedge_handle>* path = new std::vector <Polyhedron_3_::Halfedge_handle>();

        Polyhedron_3_::Vertex_handle current_vertex = start->opposite()->vertex();

        while(current_vertex != end->vertex()){

            float min_distance;
            bool found = false;

            Point_3 p = current_vertex->point();

            printf("%f %f %f\n", p.x(), p.y(), p.z());

            //go around the vertex, decide which direction to go in
            Polyhedron_3_::Halfedge_around_vertex_circulator he_it = current_vertex->vertex_begin();
            Polyhedron_3_::Halfedge_around_vertex_circulator best_halfedge = he_it;
            do{
                //we only consider if this edge or it's opposite isn't already used
                if(marked_edges.find(he_it) == marked_edges.end()
                    && marked_edges.find(he_it->opposite()) == marked_edges.end()){

                     Point_3 q = he_it->opposite()->vertex()->point();
                    printf("check %f %f %f -> %f %f %f \n", p.x(), p.y(), p.z(), q.x(), q.y(), q.z());

                    //greedy algorithm - pick the closest to the other point
                    float distance = (he_it->opposite()->vertex()->point() - end->vertex()->point()).squared_length();
                    if (!found || (distance < min_distance)){
                        printf("maxed \n");
                        min_distance = distance;
                        best_halfedge = he_it;
                        found = true;
                    };
                }

            }while(++he_it != current_vertex->vertex_begin());

            if (!found){
                //backtrack
                printf("not found\n");
                path->pop_back();
                current_vertex = path->back()->vertex();
            }
            else{
                printf("found\n");
                //add edge to list
                path->push_back(best_halfedge);
                printf("?\n");

                //set current vertex
                current_vertex = best_halfedge->opposite()->vertex();
                printf("found\n");
            }
        }

        sub_paths.push_back(*path);
    }


    void get_cuts(){

        boost::unordered_map<Polyhedron_3_::Facet_handle, bool> marked_faces;
        boost::unordered_map<Polyhedron_3_::Vertex_handle, bool> marked_vertices;
        boost::unordered_map<Polyhedron_3_::Halfedge_handle,bool> marked_edges;

        marked_faces.insert(std::make_pair(mesh.facets_begin(), true));
        bool keep_cutting = true;

        printf("Begin edge mark \n");

        std::vector<Polyhedron_3_::Halfedge_handle> halfedges = std::vector<Polyhedron_3_::Halfedge_handle>();
        for(Polyhedron_3_::Halfedge_iterator he_it = ++mesh.halfedges_begin();
              ++he_it != mesh.halfedges_end(); ++he_it)
            halfedges.push_back(he_it);
        prev_halfedge = halfedges[0];

        //while there remains an edge e adjacent to only one triangle t Remove e and t
        while(keep_cutting){

            keep_cutting = false; //assume we're not going to find anything

            std::sort(halfedges.begin(), halfedges.end(), sorter );

            printf("search through halfedges \n");

            //search through halfedges
             for(int i = 0; i < halfedges.size(); i++) {
                Polyhedron_3_::Halfedge_iterator he_it = halfedges[i];

                printf("checking edge: border %o, marked %o \n",he_it->is_border_edge(), marked_edges.find(he_it) == marked_edges.end());

                //if not a border and this edge is not marked
                if(!he_it->is_border_edge() && marked_edges.find(he_it) == marked_edges.end()){

                    bool incident_facet_marked = marked_faces.find(he_it->facet()) != marked_faces.end();
                    bool opposite_facet_marked = marked_faces.find(he_it->opposite()->facet()) != marked_faces.end();

                    printf("checking edge: indicident %o, oppo %o \n", incident_facet_marked, opposite_facet_marked);

                    //If this is the kind of edge we want
                    if (incident_facet_marked != opposite_facet_marked){

                        //Mark this edge pair
                        marked_edges.insert(std::make_pair(he_it,true));
                        marked_edges.insert(std::make_pair(he_it->opposite(),true));

                    Point_3 p = he_it->vertex()->point();
                    Point_3 q = he_it->opposite()->vertex()->point();

                    printf("%f %f %f -> %f %f %f \n", p.x(), p.y(), p.z(), q.x(), q.y(), q.z());

                        //Mark the correct face as well
                        if (!incident_facet_marked)
                            marked_faces.insert(std::make_pair(he_it->facet(), true));
                        else
                            marked_faces.insert(std::make_pair(he_it->opposite()->facet(), true));
                    
                        printf("marked edge pair\n");
                        prev_halfedge = he_it;

                        keep_cutting = true;
                        break;
                    };
                };
                printf("End halfedge looop\n");
            };
        };

        keep_cutting = true;
        printf("Begin vertex mark\n");
        
        
        //while there remains a vertex v adjacent to only one edge e Remove v and e
        while(keep_cutting){

            keep_cutting = false; //assume we're not going to find anything
           
            printf("search through vertices\n");
           
            //search through vertices
            for(Polyhedron_3_::Vertex_iterator v_it = mesh.vertices_begin();
              v_it != mesh.vertices_end(); ++v_it) {

                Point_3 p = v_it->point();
                printf("%f %f %f -> \n", p.x(), p.y(), p.z());


                printf("Vertex marked? %o\n",marked_vertices.find(v_it) == marked_vertices.end());
                
                //if this vertex is not marked
                if(marked_vertices.find(v_it) == marked_vertices.end()){

                    //check how many non-marked edges are connected
                    int non_marked_edge_count = 0;
                    Polyhedron_3_::Halfedge_handle lone_edge;

                    printf("counting non marked edge\n");

                    Polyhedron_3_::Halfedge_around_vertex_circulator he_it = v_it->vertex_begin();
                    do{
                        //if this edge is not marked

                        p = he_it->vertex()->point();
                        Point_3 q = he_it->opposite()->vertex()->point();

                        printf("%f %f %f -> %f %f %f \n", p.x(), p.y(), p.z(), q.x(), q.y(), q.z());

                        printf("check edge\n");
                        if(marked_edges.find(he_it) == marked_edges.end()){
                            printf("non marked edge\n");
                            non_marked_edge_count++;
                            lone_edge= he_it;
                        }

                    } while(++he_it != v_it->vertex_begin());

                    printf("non marked edge count: %d\n", non_marked_edge_count);

                    //If this is a lone vertex
                    if (non_marked_edge_count == 1){

                        printf("marking pair\n");
                        p = lone_edge->vertex()->point();
                        Point_3 q = lone_edge->opposite()->vertex()->point();

                        printf("%f %f %f -> %f %f %f \n", p.x(), p.y(), p.z(), q.x(), q.y(), q.z());

                        //Mark this edge pair
                        marked_edges.insert(std::make_pair(lone_edge,true));
                        marked_edges.insert(std::make_pair(lone_edge->opposite(),true));

                        //Mark this vertex
                        marked_vertices.insert(std::make_pair(v_it,true));

                        keep_cutting = true;
                    }
                }

                printf("end vertex iteration\n");
            };
        };


        printf("pulling cuts\n");

        //Pull out the cuts
        for(Polyhedron_3_::Vertex_iterator v_it = mesh.vertices_begin();
              v_it != mesh.vertices_end(); ++v_it) { //Attempt to start loop at each vertex
 

            printf("point: %f %f %f\n",v_it->point().x(), v_it->point().y(), v_it->point().z());

            printf("Vertex marked? %o\n",marked_vertices.find(v_it) == marked_vertices.end());

            //if this vertex is not marked
            if(marked_vertices.find(v_it) == marked_vertices.end()){
                //start a new loop
                std::vector<Polyhedron_3_::Halfedge_handle>* cut = new std::vector<Polyhedron_3_::Halfedge_handle>();
                bool finished_loop = false;
                Polyhedron_3_::Vertex_handle loop_vertex = v_it;
                marked_vertices.insert(std::make_pair(loop_vertex,true)); //mark vertex

                printf("start loop\n");

                while (!finished_loop){

                    finished_loop = true; 

                    printf("look in each direction\n");

                    //look for the direction to go in
                    Polyhedron_3_::Halfedge_around_vertex_circulator he_it = loop_vertex->vertex_begin();
                    
                    do{

                        printf("opp point: %f %f %f\n",he_it->opposite()->vertex()->point().x(), he_it->opposite()->vertex()->point().y(), he_it->opposite()->vertex()->point().z());

                        printf("edge marked? %o\n",marked_edges.find(he_it) == marked_edges.end());

                        //if not a border edge and this edge is not marked
                        if(marked_edges.find(he_it) == marked_edges.end()){

                            //mark edges (shouldn't be able to go backwards)
                            marked_edges.insert(std::make_pair(he_it,true));
                            marked_edges.insert(std::make_pair(he_it->opposite(),true));

                            printf("mark edges?\n");

                            //if the next vertex is not marked
                            if (marked_vertices.find(he_it->opposite()->vertex()) == marked_vertices.end()){
                                printf("add loop vertex\n");

                                loop_vertex = he_it->opposite()->vertex();
                                //mark the vertex
                                marked_vertices.insert(std::make_pair(loop_vertex,true)); 
                                //record the edge
                                cut->push_back(he_it);

                                finished_loop = false; //keep going
                                break; //avoid dupes
                            }
                            else{
                                printf("do not add loop vertex, finish loop\n");
                                finished_loop = true;  //otherwise we've reached the end
                            }
                            
                        } 
                    }while(++he_it != loop_vertex->vertex_begin());
                     
                }
                printf("add cut\n");

                for (int i = 0; i < cut->size(); ++i)
                {
                    Point_3 p = (*cut)[i]->vertex()->point();
                    Point_3 q = (*cut)[i]->opposite()->vertex()->point();
                    printf("hi\n");
                    printf("%f %f %f -> %f %f %f \n", p.x(), p.y(), p.z(), q.x(), q.y(), q.z());
                }

                sub_paths.push_back((*cut));
            }
        }

        //straighten cuts
       //for (int i = 0; i < cuts.size(); i++){
         //   cuts[i] = straighten_path(cuts[i], 10);
        //}

        
        //create a cut that passes through every other cut
        int cuts_count = sub_paths.size()-1;
        for (int i = 0; i < cuts_count; i++){
            Polyhedron_3_::Halfedge_handle start = sub_paths[i].back();
            Polyhedron_3_::Halfedge_handle end = sub_paths[i+1].front();

            add_path_between_two_vertices(start,end);
        }

        //concatenate cuts into one seam
        boost::unordered_map< std::pair<Polyhedron_3_::Vertex_handle, std::pair<int, int> >, bool> used_vertices = boost::unordered_map< std::pair<Polyhedron_3_::Vertex_handle, std::pair<int, int> >, bool>();
       
        int cut_id = 0;
        int vertex_id = 0;
 
        //while we haven't ended the seam
        while (vertex_id < sub_paths[cut_id].size() ){

            printf("%d %d", cut_id, vertex_id);

            //log current vertex
            Polyhedron_3_::Vertex_handle vertex = sub_paths[cut_id][vertex_id]->vertex();
            seam.push_back(vertex);
            used_vertices.insert(std::make_pair(std::make_pair(vertex, 
                std::make_pair(cut_id, vertex_id)), true));

            //search all nodes for dupilicate use of this vertex
            for (int i = 0; i < sub_paths.size(); ++i)
                for (int j = 0; j < sub_paths[i].size(); ++j)     
                    if (sub_paths[i][j]->vertex() == vertex){  //if we find a repeated vertex
                        if (used_vertices.find(std::make_pair(vertex,std::make_pair(i,j))) == used_vertices.end()){
                            //that hasn't been used before
                            cut_id = i;
                            vertex_id = j; //switch to that vertex
                            break;
                        }
                    }

            vertex_id++; //continue along this path
        }

        if (seam.size() == 0){
            //If empty, add in the first vertex
            Polyhedron_3_::Halfedge_handle he = mesh.halfedges_begin();
            seam.push_back(he->vertex());
            seam.push_back(he->opposite()->vertex());

            std::vector<Polyhedron_3_::Halfedge_handle>* cut = new std::vector<Polyhedron_3_::Halfedge_handle>();
            cut->push_back(he);
            sub_paths.push_back((*cut));
        }

    };

    float path_length(std::vector<Polyhedron_3_::Vertex_handle> path){
        float distance = 0;

        for (int i = 1; i < path.size(); ++i)
        {
            distance += sqrt((path[i-1]->point()-path[i]->point()).squared_length());
        }
        return distance;
    }

    float path_length(std::vector<Polyhedron_3_::Halfedge_handle> path){
        float distance = 0;

        for (int i = 0; i < path.size(); ++i)
        {
                printf("%d..",i);

            distance += sqrt((path[i]->vertex()->point()-path[i]->opposite()->vertex()->point()).squared_length());
        }
        return distance;
    }

    //brute force search between vertices on a polyhedron
    std::vector<Polyhedron_3_::Halfedge_handle> search_for_optimal_path(Polyhedron_3_::Halfedge_handle he, Polyhedron_3_::Halfedge_handle he_2, int max_size){

        boost::unordered_map<Polyhedron_3_::Halfedge_handle,bool> marked_edges;
        std::vector<Polyhedron_3_::Halfedge_handle> best_path;
        bool path_found = false;
        bool all_paths_searched = false;

        printf("Start optimal path search");

        while(!all_paths_searched){

            float min_distance = 0;
            std::vector<Polyhedron_3_::Halfedge_handle> halfedge_list;
            Polyhedron_3_::Vertex_handle current_vertex = he->vertex();
            
            printf("Considering new path");

            //steps limited by max size
            for (int i = 0; i < max_size; ++i)
            {
                Polyhedron_3_::Halfedge_handle next_edge;

                //decide which direction to travel in
                Polyhedron_3_::Halfedge_around_vertex_circulator he_it = current_vertex->vertex_begin();
                bool edge_found = false;
                do{
                    //if this edge is not marked
                    if(marked_edges.find(he_it) == marked_edges.end()){
                        next_edge = he_it;
                        edge_found = true;
                    }
                }while(++he_it != current_vertex->vertex_begin());

                //if no edge available, backtrack & mark previous edge
                if (!edge_found){
                    printf("backtracking");
                    marked_edges.insert(std::make_pair(halfedge_list.back(), true));
                    halfedge_list.pop_back();
                    i -= 2; 
                    if (i < 0){ //can't backtrack anymore
                        all_paths_searched = true;
                        break;
                    }
                }
                //otherwise, add new edge to list
                else{
                    current_vertex = next_edge->opposite()->vertex();
                    halfedge_list.push_back(next_edge);

                    Point_3 p = current_vertex->point();
                    printf("point: %f %f %f\n",p .x(), p .y(), p .z());

                    //If we reach end of a path
                    if (current_vertex == he_2->opposite()->vertex()){

                    printf("end of path");

                        float distance = path_length(best_path);

                        //calculate path length
                        if (!path_found || distance < min_distance){
                            min_distance = distance;
                            path_found = true;
                            best_path = halfedge_list;
                            printf("best path");
                        }

                        printf("backtracking");

                        //Backtrack
                        marked_edges.insert(std::make_pair(halfedge_list.back(), true));
                        halfedge_list.pop_back();
                        i -= 2; 
                        if (i < 0){ //can't backtrack anymore
                            all_paths_searched = true;
                            break;
                        }
                    }
                }
            }

            //if we reach limit without anything
            if (halfedge_list.size() == max_size){
                printf("backtracking");
                marked_edges.insert(std::make_pair(halfedge_list.back(), true));
            }
            
        }

        return best_path;
    }

    //applies brute force down a path
    std::vector<Polyhedron_3_::Halfedge_handle> straighten_path(std::vector<Polyhedron_3_::Halfedge_handle>& seam, int max_step){

        for (int i = 0; i < seam.size(); i+= max_step){
            std::vector<Polyhedron_3_::Halfedge_handle> path = search_for_optimal_path(seam[i], seam[i+max_step],max_step*2);
            for (int j=0; j< max_step; j++)
                seam[i+j] = path[i];
        }        

        for (int i = max_step/2; i < seam.size(); i+= max_step){
            std::vector<Polyhedron_3_::Halfedge_handle> path = search_for_optimal_path(seam[i], seam[i+max_step],max_step*2);
            for (int j=0; j< max_step; j++)
                seam[i+j] = path[i];
        }
        return seam;
    }

    Parameterization fit_border(int resolution){

       printf("cuts done \n");

        //allocate distances
        std::vector<float> distance_allocations  = std::vector<float>();
        float total_cut_distance = 0;

        for (int i = 0; i < seam.size(); i++){
            printf("%f..",seam[i]->point().x());

        }

        for (int i = 0; i < sub_paths.size(); i++){
            float dist = path_length(sub_paths[i]);
            total_cut_distance += dist;                
            distance_allocations.push_back(dist);
        }

       printf("rounding distance \n");

        //round distances
        std::vector<int> pixel_allocations;
        int points_left =  4*(resolution-1);
       printf("%d \n", distance_allocations.size());

        for (int i = 0; i < distance_allocations.size()-1; ++i){
            int allocation = distance_allocations[i]*4*(resolution-1)/total_cut_distance;
            pixel_allocations.push_back(allocation);
            points_left -= allocation;
        }
        pixel_allocations.push_back(points_left);

       printf("fit around border \n");

        //fit around border
        float x = 0; float y = 0;
        float x_dir = 1; float y_dir = 0; 
        for (int i = 0; i < sub_paths.size(); ++i)
        {
            printf("fit path %d \n", sub_paths[i].size());

            if (sub_paths[i].size() == 0)
                break;

            //place parameterization down every path
            for (int j = 0; j < sub_paths[i].size(); ++j)
            {
                printf("INSERT \n");
                printf("%f \n", sub_paths[i][j]->vertex()->point().x());

                std::pair<float, float> coord = std::make_pair(x, y);
                parameterization.insert(std::make_pair(sub_paths[i][j]->vertex(),coord));
                printf("legnth \n");

                std::vector<Polyhedron_3_::Halfedge_handle> vec; vec.push_back(sub_paths[i][j]);
                float distance_change = path_length(vec)*4.0/total_cut_distance;///(pixel_allocations[i]/resolution);
                x += x_dir * distance_change;
                y += y_dir * distance_change;

                //(0,0) -> (1,0) -> (1,1) -> (0,1) -> (0, 0)
                if (x > 1) { y = (x-1); y_dir = 1; x = 1; x_dir = 0;}
                if (y > 1) { x = (y-1); x_dir = -1; x = 1; y_dir = 0;}
                if (x < 0) { y = (1-x); y_dir = -1; x = 0; x_dir = 0;}
                if (y < 0) { x = (1-y); x_dir = 1; x = 0; y_dir = 0;}
            }

            std::pair<float, float> coord = std::make_pair(x, y);
            int j = sub_paths[i].size() -1;
            parameterization.insert(std::make_pair(sub_paths[i][j]->opposite()->vertex(),coord));
        }

        //split degenerate triangles
        //split corner edges
        //rotate

       printf("initial parameterization \n");

        //come up with a basic parameterization for the rest of the mesh
        std::queue<Polyhedron_3_::Vertex_handle> vertex_list;
        for (int i = 0; i < seam.size(); ++i)
            vertex_list.push(seam[i]);

        while (vertex_list.size() > 0){
            Polyhedron_3_::Vertex_handle current_vertex = vertex_list.front();
            vertex_list.pop();

            //Traverse neighbouring vertices
            Polyhedron_3_::Halfedge_around_vertex_circulator he_it = current_vertex->vertex_begin();
            float total_distance = 0;
            float weighted_x = 0;
            float weighted_y = 0;
             do{
                //if this vertex is already parametrized
                //coords are the average of neighbors, weighted by distance
                if(parameterization.find(he_it->opposite()->vertex()) != parameterization.end()){
                        float distance = sqrt((he_it->opposite()->vertex()->point()- current_vertex->point()).squared_length());
                        std::pair<float,float> coords = (parameterization.find(he_it->opposite()->vertex()))->second;
                       
                    weighted_x += distance* std::get<0>(coords);
                    weighted_y += distance* std::get<1>(coords);
                    total_distance += distance;
                }   
                else{
                    vertex_list.push(he_it->opposite()->vertex());
                }
            }while(++he_it != current_vertex->vertex_begin());

            //if this (original) vertex is not yet marked
            if(parameterization.find(current_vertex) == parameterization.end()){
                std::pair<float,float> coords = std::make_pair(weighted_x/total_distance, weighted_y/total_distance);
                parameterization.insert(std::make_pair(current_vertex, coords));
            }
        }
        
        return parameterization; 
    }

    float facet_area(Polyhedron_3_::Facet_iterator facet){
        
        //Circulate around facet
        Polyhedron_3_::Halfedge_around_facet_circulator halfedge_it = facet->facet_begin(); 
         
        double facet_area = 0;

        Polyhedron_3_::Halfedge_around_vertex_circulator halfedge_face_it = halfedge_it->vertex_begin();
        Point_3 point1 = halfedge_face_it->vertex()->point(); ++halfedge_face_it;
        Point_3 point2 = halfedge_face_it->vertex()->point(); ++halfedge_face_it;
        Point_3 point3;

        do
        {  
            //sum area via triangles
            point3 = halfedge_face_it->vertex()->point();
            Vector_3 dir1 = point1 -point2;
            Vector_3 dir2 = point1 -point3;
                
            Vector_3 cross_product = Vector_3(dir1.y()*dir2.z()-dir2.y()*dir1.z(),
                            dir1.z()*dir2.x()-dir2.z()*dir1.x(),dir1.x()*dir2.y()-dir2.x()*dir1.y());

            facet_area += sqrt(cross_product.squared_length())*0.5;

            point1 = point2;
            point2 = point3;

        }while(++halfedge_face_it != halfedge_it->vertex_begin());

        return facet_area;
    }

    void invert_coords(Polyhedron_3_::Vertex_handle v1, Polyhedron_3_::Vertex_handle v2,
        Polyhedron_3_::Vertex_handle v3, float& dir_1_edge_1,
        float& dir_1_edge_2, float& dir_2_edge_1, float& dir_2_edge_2){
           
            std::pair<float,float> coords1 = parameterization.find(v1)->second;
            std::pair<float,float> coords2 = parameterization.find(v2)->second;
            std::pair<float,float> coords3 = parameterization.find(v3)->second;

            //Solve for the linear combination of two edges, needed to parametrize either way
            float a_x = std::get<0>(coords2) - std::get<0>(coords1);
            float a_y = std::get<1>(coords2) - std::get<1>(coords1);
            float b_x = std::get<0>(coords3) - std::get<0>(coords1);
            float b_y = std::get<1>(coords3) - std::get<1>(coords1);

            //Simple 2x2 matrix inversion
            float det = a_x*b_y - a_y*b_x;
            dir_1_edge_1 = b_y/det;
            dir_1_edge_2 = -a_y/det;
            dir_2_edge_1 = -b_x/det;
            dir_2_edge_2 = a_y/det;
    }

    float geometric_stretch(Parameterization parameterization){
        
       printf("geometric_stretch \n");

        float stretch = 0;

        for(Polyhedron_3_::Facet_iterator f_it = mesh.facets_begin();
              f_it != mesh.facets_end(); ++f_it) {
                
            printf("Jacobian \n");

            //Compute Jacobian for each face
            Polyhedron_3_::Halfedge_around_facet_circulator he = f_it->facet_begin ();

            printf("vertices \n");
            printf("%f \n", he->vertex()->point().x());

            //Invert coords
            Polyhedron_3_::Vertex_handle v1 = he->vertex();
            Polyhedron_3_::Vertex_handle v2 = (++he)->vertex();
            Polyhedron_3_::Vertex_handle v3 = (++he)->vertex();

            float dir_1_edge_1;
            float dir_1_edge_2;
            float dir_2_edge_1;
            float dir_2_edge_2;


            invert_coords(v1, v2, v3,
                dir_1_edge_1, dir_1_edge_2, dir_2_edge_1, dir_2_edge_2);
            
            printf("Vecots \n");

            //rows of Jacobian
            Vector_3 df_dx = dir_1_edge_1 * (v2->point() - v1->point()) + dir_1_edge_2 * (v3->point() - v1->point());
            Vector_3 df_dy = dir_2_edge_1 * (v2->point() - v1->point()) + dir_2_edge_2 * (v3->point() - v1->point());

            float facet_stretch = df_dx.squared_length() + df_dy.squared_length();

            printf("Area \n");

            //Compute face area
            stretch += facet_stretch * facet_area(f_it);

        }
        return stretch;
    }

    void optimize_vertex(Polyhedron_3_::Vertex_handle vertex, Parameterization& parameterization){

        float best_error = 0;
        float new_error = 0;

        std::pair<float,float> coords = parameterization.find(vertex)->second;
        float x = coords.first;
        float y = coords.second;
       
        printf("optimize_vertex \n");
        
        do {

            printf("iteration \n");

            //random search direction
            float vector_dir_x = rand();
            float vector_dir_y = rand();
            
            //make sure it is in the "cone"?

            //binary line search
            float pertubation_length = 1;
            for (int i = 0; i < 20; ++i)
            {
                printf("%d...",i);

                pertubation_length *= 0.5;
                float new_x = x+vector_dir_x*pertubation_length;
                float new_y = y+vector_dir_y*pertubation_length;

                parameterization.erase(vertex);
                parameterization.insert(std::make_pair(vertex, std::make_pair(x,y)));

                //compute the error, update coords if improvement
                new_error = geometric_stretch(parameterization);
                if (new_error < best_error){
                    x = new_x;
                    y = new_y;
                    best_error = new_error;
                }
            }
        }while(best_error == new_error); //keep going, whilst we're reducing error
        

    }

    Parameterization parameterize(int resolution, int iterations){

        fit_border(resolution);
        printf("start \n");

        for (int i = 0; i < iterations; ++i)
        {
            printf("%d..", i);

            //Optimize each vertex in a line
            for(Polyhedron_3_::Vertex_iterator v_it = mesh.vertices_begin();
              v_it != mesh.vertices_end(); ++v_it) {
                optimize_vertex(v_it, parameterization);
            }
        }

        return parameterization;
    }

    Point_3 sample_geometry(float x, float y){

        CGAL::Epick::Point_2 sample_point(x,y);

        //Find facet that contains these coordinates
        for(Polyhedron_3_::Facet_iterator f_it = mesh.facets_begin();
              f_it != mesh.facets_end(); ++f_it) {

            //Build list of coordinates
            std::vector<CGAL::Epick::Point_2> coord_points;
            Polyhedron_3_::Halfedge_around_facet_circulator halfedge_it = f_it->facet_begin(); 
            do{

                std::pair<float,float> coords = parameterization.find(halfedge_it->vertex())->second;
                CGAL::Epick::Point_2 point(std::get<0>(coords), std::get<1>(coords));
                coord_points.push_back(point);
            }while(++halfedge_it != f_it->facet_begin());

            //If this facet contains the point
            if (CGAL::bounded_side_2(coord_points.begin(), coord_points.end(), 
                sample_point, CGAL::Epick())==CGAL::ON_BOUNDED_SIDE ){
            
                //Consider 3 vertices only    
                Polyhedron_3_::Vertex_handle v1 = halfedge_it->vertex();
                Polyhedron_3_::Vertex_handle v2 = (++halfedge_it)->vertex();
                Polyhedron_3_::Vertex_handle v3 = (++halfedge_it)->vertex();

                // Interpolate, start by inverting coord system
                std::pair<float,float> coords = parameterization.find(v1)->second;
                float x_dir = x - std::get<0>(coords);
                float y_dir = y - std::get<1>(coords);

                float dir_1_edge_1; float dir_1_edge_2;
                float dir_2_edge_1; float dir_2_edge_2;
                invert_coords(v1, v2, v3,
                    dir_1_edge_1, dir_1_edge_2, dir_2_edge_1, dir_2_edge_2);

                float x_edge = x_dir*dir_1_edge_1 + y_dir*dir_2_edge_1;
                float y_edge = x_dir*dir_2_edge_1 + y_dir*dir_2_edge_2;

                Point_3 geometry = v1->point() + (v2->point()-v1->point())*x_edge
                         + (v3->point()-v1->point())*y_edge;
                return geometry;
            }
        }

        //If facet not found
        return Point_3(0,0,0);
    }
};
#endif //SWIG_CGAL_MESH_CUTS_3_H
