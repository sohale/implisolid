#pragma once

//#pragma message ( "*******************************" )

vectorized_faces f2f(const std::vector<std::vector<vertexindex_type>>& f) {
    vectorized_faces_shape
        shape =
            vectorized_faces_shape{static_cast<vectorized_faces::index>(f.size()), 3};

    vectorized_faces faces {shape};
    for(int i=0; i < f.size(); ++i) {
        for(int j=0; j < 3; ++j) {
            faces[i][j] = f[i][j];
        }
    }
    return faces;
}

vectorized_vect v2v(const std::vector<std::vector<REAL>>& v, REAL scale) {
    vectorized_vect_shape
        shape =
            vectorized_vect_shape{static_cast<vectorized_vect::index>(v.size()), 3};

    vectorized_vect vects {shape};
    for(int i=0; i < v.size(); ++i) {
        for(int j=0; j < 3; ++j) {
            vects[i][j] = v[i][j] * scale;
        }
    }
    return vects;
}

vectorized_vect  vects2vects(
    const std::vector<REAL>& result_verts
) {
    boost::array<vectorized_vect::index, 2> verts_shape = { (vectorized_vect::index)(result_verts.size()/3) , 3 };
    vectorized_vect  verts(verts_shape);

    int output_verts = 0;
    auto i = result_verts.begin();
    auto e = result_verts.end();
    for ( ; i != e; i++, output_verts++) {
        verts[output_verts][0] = (*i);
        i++;
        verts[output_verts][1] = (*i);
        i++;
        verts[output_verts][2] = (*i);
    }
    return verts;
}

void set_vectorverts_from_vectorised_verts(
    std::vector<REAL>& result_verts,
    const vectorized_vect & verts
) {
    auto n = verts.shape()[0];
    for (int i=0; i < n; i++) {
        result_verts[i*3+0] = verts[i][0];
        result_verts[i*3+1] = verts[i][1];
        result_verts[i*3+2] = verts[i][2];
    }
}


// converts std::vector<> into boost::multi_array<>
boost::multi_array<vertexindex_type, 2> copy_faces_from_vectorfaces(
    const std::vector<vertexindex_type> & mesh_faces
) {
    // boost::multi_array<int, 2> faces = ;

    vectorized_vect::index  num_faces = static_cast<vectorized_vect::index>(mesh_faces.size()/3);

    boost::array<vectorized_vect::index, 2> faces_shape = { num_faces , 3 };
    vectorized_faces faces(faces_shape);

    int output_faces = 0;
    auto i_f = mesh_faces.begin();
    auto e_f = mesh_faces.end();
    for ( ; i_f != e_f; i_f++, output_faces++) {
        faces[output_faces][0] = (*i_f);
        i_f++;
        faces[output_faces][1] = (*i_f);
        i_f++;
        faces[output_faces][2] = (*i_f);
    }
    assert(num_faces == output_faces);

    return faces;
}


/**
 * Note: This function does require the target to be of exact same size,
 * but the target ,ust have enough elements.
 * The rest are left untouched.
 */
void copy_vertices_into_stdvector_incomplete(
    const vectorized_vect & verts_source,
          std::vector<REAL> & target_verts_vector
) {
    assert( target_verts_vector.size() >= verts_source.shape()[0] * 3 );
    //auto num_verts = verts_source.shape()[0];
    auto target_begin = target_verts_vector.begin();
    auto source_begin = verts_source.begin();
    auto source_end   = verts_source.end();
    for ( ; source_begin != source_end; ++source_begin) {
        *target_begin = (*source_begin)[0];
        ++target_begin;
        *target_begin = (*source_begin)[1];
        ++target_begin;
        *target_begin = (*source_begin)[2];
        ++target_begin;
    }
    assert(target_begin <= target_verts_vector.end());
}
