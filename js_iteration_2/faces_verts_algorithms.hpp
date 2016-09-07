#pragma once

std::vector< std::vector<int>> make_neighbour_faces_of_vertex(const verts_t& verts, const faces_t& faces) {
    int nt = faces.shape()[0];
    int vt = verts.shape()[0];
    std::vector< std::vector<int>> neighbour_faces_of_vertex;
    for (int fi=0; fi < vt; fi++) {
        neighbour_faces_of_vertex.push_back(std::vector<int>());
    }
    for (int fi=0; fi < nt; fi++) {
        for (int vi=0; vi < 3; vi++) {
            int v1 = faces[fi][vi];
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
