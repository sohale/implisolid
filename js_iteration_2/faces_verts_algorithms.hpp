#pragma once


std::vector< std::vector<int>> make_neighbour_faces_of_vertex(const faces_t& faces, vertexindex_type  max_vert_index) {
    /*
    note: max_vert_index could be derived from maximum index in faces, but 1-It is safer, may prevent future capacity increase? (not really needed though) 2- It will not be much smaller enyway
    This will be replaced by the more efficient "sparse" matrices (or equivalent data structures) anyway.
    */
    vertexindex_type num_verts = max_vert_index;  // verts.shape()[0];
    std::vector< std::vector<int>> neighbour_faces_of_vertex;
    for (vertexindex_type vi=0; vi < num_verts; vi++) {
        neighbour_faces_of_vertex.push_back(std::vector<int>());
    }

    // todo: initialise using C++ RAII principle:
    // std::vector< std::vector<int>> neighbour_faces_of_vertex(max_vert_index);

    // static_assert(vertexindex_type == neighbour_faces_of_vertex::index_type);
    // static_assert(faces_t::value_type == vertexindex_type);

    int num_faces = faces.shape()[0];
    for (int fi=0; fi < num_faces; fi++) {
        for (int side_i=0; side_i < 3; side_i++) {
            vertexindex_type v1 = faces[fi][side_i];
            neighbour_faces_of_vertex[v1].push_back(fi);
        }
    }

    return neighbour_faces_of_vertex;
}

/*
    std::vector<REAL> result_verts;
    std::vector<int> result_faces;
*/

vectorized_vect convert_vectorverts_to_vectorized_vect(const std::vector<REAL> & result_verts) {

    boost::array<int, 2> verts_shape = { static_cast<int>( result_verts.size() / 3 ), 3 };
    vectorized_vect  verts(verts_shape);

    int output_verts=0;
    auto i = result_verts.begin();
    auto e = result_verts.end();
    for(; i != e; i++, output_verts++) {
        verts[output_verts][0] = *(i);
        verts[output_verts][1] = *(i + 1);
        verts[output_verts][2] = *(i + 2);
        i++;
        i++;
    }
    return verts;
}

vectorized_faces  convert_vectorfaces_to_vectorized_faces(const std::vector<int> & result_faces) {

    vectorized_faces_shape  faces_shape = { static_cast<int>( result_faces.size()/3 ) , 3 };
    vectorized_faces  faces(faces_shape);

    int output_faces = 0;
    auto i_f = result_faces.begin();
    auto e_f = result_faces.end();
    for(; i_f!=e_f; i_f++, output_faces++) {
        faces[output_faces][0] = *(i_f);
        i_f++;
        faces[output_faces][1] = *(i_f);
        i_f++;
        faces[output_faces][2] = *(i_f);
    }
    return faces;
}

void replace_vectorverts_from_vectorized_vect(std::vector<REAL> & result_verts, const vectorized_vect & new_verts) {
    #if ASSERT_USED
        REAL cumul_abs_displacement = 0;
    #endif
        auto n = new_verts.shape()[0];
    for (int i=0; i < n; i++) {
        // std::clog << "result_verts  = new_verts: " <<result_verts[i*3+0] << " = " << new_verts[i][0] << std::endl;
        #if ASSERT_USED
            for(int j=0; j < 3; j++) {
                cumul_abs_displacement += std::abs( result_verts[i*3 + j] - new_verts[i][j]);
            }
        #endif
        result_verts[i*3 + 0] = new_verts[i][0];
        result_verts[i*3 + 1] = new_verts[i][1];
        result_verts[i*3 + 2] = new_verts[i][2];
    }
    #if ASSERT_USED
        std::clog << "<cumul |displacement|> = " <<  cumul_abs_displacement/((REAL)(new_verts.shape()[0]))/3 << std::endl;
    #endif
}
