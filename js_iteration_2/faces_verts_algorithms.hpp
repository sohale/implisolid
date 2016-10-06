#pragma once


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

vectorized_faces  convert_vectorfaces_to_vectorized_faces(const std::vector<vertexindex_type> & result_faces) {

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
