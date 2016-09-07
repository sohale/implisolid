#pragma once

ofstream debug_part_1(string filename, std::vector<REAL>& result_verts) {
    ofstream f_out(filename);

    f_out << "0ld vertex :" << endl;
    for (int i=0; i< result_verts.size()/3.; i++){
        f_out << result_verts[3*i];
        f_out << " ";
        f_out << result_verts[3*i+1];
        f_out << " ";
        f_out << result_verts[3*i+2];
        f_out <<  "\n";
    }
    f_out << endl;
    return f_out;
}


void debug_part_2(ofstream& f_out, std::vector<REAL>& result_verts, std::vector<int>& result_faces, vectorized_vect & centroids) {
      f_out << "n3w vertex :" << endl;
      for (int i=0; i< result_verts.size()/3.; i++){
          f_out << result_verts[3*i];
          f_out << " ";
          f_out << result_verts[3*i+1];
          f_out << " ";
          f_out << result_verts[3*i+2];
          f_out <<  "\n";
      }

      f_out << "faces :" << endl;
      for (int i=0; i< result_faces.size()/3.; i++){
          f_out << result_faces[3*i];
          f_out << " ";
          f_out << result_faces[3*i+1];
          f_out << " ";
          f_out << result_faces[3*i+2];
          f_out <<  "\n";
      }

      f_out << "centroids:" << endl;
      for (int i=0; i< centroids.shape()[0]; i++){
          f_out << centroids[i][0];
          f_out << " ";
          f_out << centroids[i][1];
          f_out << " ";
          f_out << centroids[i][2];
          f_out <<  "\n";
      }
      f_out.close();
}




void writing_test_file_(
        string  output_file_name,
        std::vector<REAL>&result_verts, std::vector<int>&result_faces,
        vectorized_vect & new_verts, vectorized_vect & centroids
    )
{

    ofstream f_out(output_file_name);

    f_out << "0ld vertex :" << endl;
    for (int i=0; i< result_verts.size()/3.; i++){
        f_out << result_verts[3*i];
        f_out << " ";
        f_out << result_verts[3*i+1];
        f_out << " ";
        f_out << result_verts[3*i+2];
        f_out <<  "\n";
    }
    f_out << endl;

    /*
    boost::array<int, 2> verts_shape = { (int)result_verts.size()/3 , 3 };
    vectorized_vect  verts(verts_shape);
    boost::array<int, 2> faces_shape = { (int)result_faces.size()/3 , 3 };
    boost::multi_array<int, 2> faces(faces_shape);
    vectorized_vect  new_verts (verts_shape);
    vectorized_vect  centroids (faces_shape);
    process2_vertex_resampling_relaxation(new_verts, faces, verts, centroids, object, c);
    */
    //new_verts, centroids

    for (int i=0; i<new_verts.shape()[0]; i++){
        result_verts[i*3+0] = new_verts[i][0];
        result_verts[i*3+1] = new_verts[i][1];
        result_verts[i*3+2] = new_verts[i][2];
    }

    f_out << "n3w vertex :" << endl;
    for (int i=0; i< result_verts.size()/3.; i++){
        f_out << result_verts[3*i];
        f_out << " ";
        f_out << result_verts[3*i+1];
        f_out << " ";
        f_out << result_verts[3*i+2];
        f_out <<  "\n";
    }

    f_out << "faces :" << endl;
    for (int i=0; i< result_faces.size()/3.; i++){
        f_out << result_faces[3*i];
        f_out << " ";
        f_out << result_faces[3*i+1];
        f_out << " ";
        f_out << result_faces[3*i+2];
        f_out <<  "\n";
    }

    f_out << "centroids:" << endl;
    for (int i=0; i< centroids.shape()[0]; i++){
        f_out << centroids[i][0];
        f_out << " ";
        f_out << centroids[i][1];
        f_out << " ";
        f_out << centroids[i][2];
        f_out <<  "\n";
    }
    f_out.close();
}
