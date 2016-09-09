/*
This file contains functions as building blocks ("code blocks") that implement algorithms and operations on meshes.

These functions are not organised in classes. An OOP design can call these, but should not include the bodies.

No function should use any "verts". These are operations that are done on a "faces_t" variable.

These are efficient graph theoretic algorithms. Most of them involve building data structures involved in graph adjancency structure.

@author sohale
@author marcFraysse
@author schauvier
*/

#pragma once
namespace mp5_implicit {

void make_edge_lookup(faces_t faces, faces_t& edges_of_faces, faces_t& faces_of_edges){
  int nfaces = faces.shape()[0];
  cout << "nfaces is : " << nfaces << endl;
  assert(nfaces % 2 == 0);
  int num_edges = nfaces*3./2.;

  long modulo = long(num_edges);
  long lookup_array_size = modulo*num_edges + num_edges;
  map<int, int> eulookup;
  map<int, int>::iterator iter;
  int edge_counter = 0;
  for (int fi=0; fi<nfaces; fi++){
    for (int vj=0; vj<3; vj++){
      assert(fi<nfaces);
      bool new_edge = false;
      int v2j = (vj+1)%3;
      int e1 = faces[fi][vj];
      int e2 = faces[fi][v2j];
      int eu_pair_int;
      if (e2 > e1){
        eu_pair_int = int(e1 + e2*modulo);
      }
      else{
        eu_pair_int = int(e2 + e1*modulo);
      }

      iter = eulookup.find(eu_pair_int);
      if (iter== eulookup.end()){
        new_edge = true;
      }
      if (new_edge){
        int e_id = edge_counter;
        edges_of_faces[fi][vj] = e_id;
        faces_of_edges[e_id][0] = fi;
        assert (vj!= v2j);

        eulookup.insert(pair<int,int>(eu_pair_int,e_id));
        edge_counter ++;
        assert(edge_counter <= num_edges);
      }
      else{
        int e_id = iter->second;
        assert (e_id >= 0);
        edges_of_faces[fi][vj] = e_id;
        faces_of_edges[e_id][1] = fi;
        assert (vj!= v2j);
        int other_fi = faces_of_edges[e_id][0];

      }
    }
  }

}

void build_faces_of_faces(const faces_t& edges_of_faces, const faces_t& faces_of_edges, faces_t& faces_of_faces){
  for(int face = 0; face < edges_of_faces.shape()[0]; face++){
    for(int edge = 0; edge < 3; edge++){
      if(faces_of_edges[edges_of_faces[face][edge]][0]!=face)
        faces_of_faces[face][edge]=faces_of_edges[edges_of_faces[face][edge]][0];
      else{
        assert(faces_of_edges[edges_of_faces[face][edge]][1] != face);
        faces_of_faces[face][edge]=faces_of_edges[edges_of_faces[face][edge]][1];}
    }
  }
  for(int face = 0; face < faces_of_faces.shape()[0]; face ++){
    for(int faces = 0; faces<3; faces++){
      assert(face==faces_of_faces[faces_of_faces[face][faces]][0] ||
          face==faces_of_faces[faces_of_faces[face][faces]][1] ||
          face==faces_of_faces[faces_of_faces[face][faces]][2]);
    }
  }
}

}   // namespace

