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
