/*
A demo combining Emscripen + ThreeJS + MC + "implicit surfaces" (a sphere).

mc2_pb.cpp: MC based on Paul Bourke's code, version 2.

Runs a quick-and-dirty (and buggy) Marching Cubes algorithm based on Paul Bourke's code.
Usage:
    1- Make sure the boost (1_61) is extracted in the right folder (check -I cmdline option in mc1_pb.bat )
    2- Run the mc1_pb.bat
    3- Be careful: Don't add the generated js file (mc1_pb.js) to git.
    4- Run the mc1_pb_3js.html from your browser.
    5- Enjoy the sphere!

The batch file generates the "mc1_pb.js" which is used inside the html file.

Issues:
    This demo does not work when optimization is applied.

Careful:
    em++/emcc may generate js and html files. It may replace the human-coded html ad js files with the nsame name.

Versions tested:
  C++:  C++14
  emc++: 1.35.0
  Boost: 1_61_0
  ThreeJS: r77

@author Sohail

MC algorithm based on http://paulbourke.net/geometry/polygonise/
*/


#include <iostream>
using namespace std;

#define ASSERTS 1
#define VERBOSE  1
typedef float REAL;
typedef struct {
   REAL x, y, z;
} XYZ;

//#define ABS(x) ((x<0)?(-x):(x))   //Why is this wrong??
REAL ABS(REAL x){
  if(x<0)
    return -x;
  return x;
}

// -s ASSERTIONS=1 ??


//#include <algorithm>

std::ostream& operator<<( std::ostream& sout, XYZ p)
{
    sout << p.x << ",";
    sout << p.y << ",";
    sout << p.z << "";
    return sout;
}

/*
   Linearly interpolate the position where an isosurface cuts
   an edge between two vertices, each with their own scalar value
*/
XYZ VertexInterp(
    REAL isolevel,
    XYZ p1, XYZ p2,
    REAL valp1,
    REAL valp2)
{
  bool verbose = false;
   if(verbose)
      cout <<"("<< p1 << "-"<< p2<<")";
   REAL mu;
   XYZ p;

   if (ABS(isolevel-valp1) < 0.00001)
      return(p1);
   if (ABS(isolevel-valp2) < 0.00001)
      return(p2);
   if (ABS(valp1-valp2) < 0.00001)
      return(p1);
   mu = (isolevel - valp1) / (valp2 - valp1);
   p.x = p1.x + mu * (p2.x - p1.x);
   p.y = p1.y + mu * (p2.y - p1.y);
   p.z = p1.z + mu * (p2.z - p1.z);
   if(verbose)
      cout <<"->"<< p<<"  ";

   return(p);
}


typedef struct {
   XYZ p[3];
} TRIANGLE;

typedef struct {
   XYZ p[8];
   REAL val[8];
} GRIDCELL;



//template<class CharT, class Traits>
//template<class CharT>
std::ostream& operator<<(
      std::ostream& sout,
      TRIANGLE tr)
//      TRIANGLE const& tr)
{
    for(int i=0; i<3; i++){
      sout << " |";
      sout << tr.p[i].x << ",";
      sout << tr.p[i].y << ",";
      sout << tr.p[i].z << "";
      sout << "| ";
    }
    /*
    std::string s = "tr";
    std::copy(s.begin(),
                        s.end(),
                        std::ostream_iterator<CharT>(sout) );
    */
    return sout;
}
/*
   Given a grid cell and an isolevel, calculate the triangular
   facets required to represent the isosurface through the cell.
   Return the number of triangular facets, the array "triangles"
   will be loaded up with the vertices at most 5 triangular facets.
    0 will be returned if the grid cell is either totally above
   of totally below the isolevel.
*/
int Polygonise(const GRIDCELL& grid, REAL isolevel, TRIANGLE *triangles, bool verbose)
{
   int i,ntriang;
   int cubeindex;
   XYZ vertlist[12];

const int edgeTable[256]={
0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };
const int triTable[256][16] =
{{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
{3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
{3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
{3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
{9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
{9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
{8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
{9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
{3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
{4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
{5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
{9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
{0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
{10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
{5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
{9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
{0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
{1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
{2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
{7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
{11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
{11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
{1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
{9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
{2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
{6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
{6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
{5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
{1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
{6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
{3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
{10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
{10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
{1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
{0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
{10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
{3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
{7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
{10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
{2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
{7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
{2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
{10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
{6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
{8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
{6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
{8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
{5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
{9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
{7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
{6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
{0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
{5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
{11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
{2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
{1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
{9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
{8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
{11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
{4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
{0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
{9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
{1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

   /*
      Determine the index into the edge table which
      tells us which vertices are inside of the surface
   */
   cubeindex = 0;
   if (grid.val[0] < isolevel) cubeindex |= 1;
   if (grid.val[1] < isolevel) cubeindex |= 2;
   if (grid.val[2] < isolevel) cubeindex |= 4;
   if (grid.val[3] < isolevel) cubeindex |= 8;
   if (grid.val[4] < isolevel) cubeindex |= 16;
   if (grid.val[5] < isolevel) cubeindex |= 32;
   if (grid.val[6] < isolevel) cubeindex |= 64;
   if (grid.val[7] < isolevel) cubeindex |= 128;

   /* Cube is entirely in/out of the surface */
   if (edgeTable[cubeindex] == 0)
      return(0);

   /* Find the vertices where the surface intersects the cube */
   if (edgeTable[cubeindex] & 1)
      vertlist[0] =
         VertexInterp(isolevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);
   if (edgeTable[cubeindex] & 2)
      vertlist[1] =
         VertexInterp(isolevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);
   if (edgeTable[cubeindex] & 4)
      vertlist[2] =
         VertexInterp(isolevel,grid.p[2],grid.p[3],grid.val[2],grid.val[3]);
   if (edgeTable[cubeindex] & 8)
      vertlist[3] =
         VertexInterp(isolevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
   if (edgeTable[cubeindex] & 16)
      vertlist[4] =
         VertexInterp(isolevel,grid.p[4],grid.p[5],grid.val[4],grid.val[5]);
   if (edgeTable[cubeindex] & 32)
      vertlist[5] =
         VertexInterp(isolevel,grid.p[5],grid.p[6],grid.val[5],grid.val[6]);
   if (edgeTable[cubeindex] & 64)
      vertlist[6] =
         VertexInterp(isolevel,grid.p[6],grid.p[7],grid.val[6],grid.val[7]);
   if (edgeTable[cubeindex] & 128)
      vertlist[7] =
         VertexInterp(isolevel,grid.p[7],grid.p[4],grid.val[7],grid.val[4]);
   if (edgeTable[cubeindex] & 256)
      vertlist[8] =
         VertexInterp(isolevel,grid.p[0],grid.p[4],grid.val[0],grid.val[4]);
   if (edgeTable[cubeindex] & 512)
      vertlist[9] =
         VertexInterp(isolevel,grid.p[1],grid.p[5],grid.val[1],grid.val[5]);
   if (edgeTable[cubeindex] & 1024)
      vertlist[10] =
         VertexInterp(isolevel,grid.p[2],grid.p[6],grid.val[2],grid.val[6]);
   if (edgeTable[cubeindex] & 2048)
      vertlist[11] =
         VertexInterp(isolevel,grid.p[3],grid.p[7],grid.val[3],grid.val[7]);

   if(verbose)
   {
   cout << endl;

   for(int i=0;i<8;i++)
   {
       cout << " ["<<i<<"]="<< grid.p[i]<< "  " ;
   }
   cout << endl;

   cout << "edgeTable[cubeindex] = " << edgeTable[cubeindex] << " ";
   cout << endl;
   cout << "triTable[cubeindex]=" << triTable[cubeindex] << " "; // address
   cout << endl;

   cout << "grid:";
   for (int i=0;i<8;i++)
   {
        cout << "("<<grid.p[i] << ") ";
   }
   cout << endl;

   cout << "vertlist[i]:";
   for (int i=0;i<12;i++)
   {
        cout << "(" << vertlist[i] << ") ";
   }
   cout << endl;
   //for (int i=0;triTable[cubeindex][i]!=-1;i++)
   }


   /* Create the triangle */
   ntriang = 0;
   for (int i=0;triTable[cubeindex][i]!=-1;i+=3) {
      triangles[ntriang].p[0] = vertlist[triTable[cubeindex][i  ]];
      triangles[ntriang].p[1] = vertlist[triTable[cubeindex][i+1]];
      triangles[ntriang].p[2] = vertlist[triTable[cubeindex][i+2]];
      if (verbose){
      cout << "triTable[cubeindex][i="<< i <<" .. i+2]= ("
          << triTable[cubeindex][i] << " "
          << triTable[cubeindex][i+1] << " "
          << triTable[cubeindex][i+2] << ") "
          ;
      }
      ntriang++;
   }
   if(verbose){
   cout << endl;
   cout << endl;
 }

   return(ntriang);
}

void test1()
{

    GRIDCELL grid;
    for (int i=0;i<8;i++)
    {
        grid.p[i] = XYZ({1,2,3});
        grid.val[i] = 1.0;
    }


    REAL isolevel = 0;
    TRIANGLE triangles[10000];
    int r = Polygonise(grid, isolevel, triangles, true);
    cout << r;
}


#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include <math.h>
//void make_grid( boost::multi_array<REAL, 4>& grid_out, boost::multi_array<REAL, 3>& values_out )


typedef boost::multi_array<REAL, 2> verts_t;
typedef boost::multi_array<int, 2> faces_t;

//typedef int verts_t;
//typedef int faces_t;
typedef pair<verts_t, faces_t> vf_t;
//typedef pair<boost::multi_array<REAL, 2>, boost::multi_array<int, 2>> vf_t;


/*
boost::multi_array<REAL, 3>& mc(const boost::multi_array<REAL, 3>& values) {}
*/
vf_t mc(const boost::multi_array<REAL, 3>& values)
{
    vf_t vf;
    return vf;
}

//void make_XYZ(REAL pa[3], XYZ&output){
//const boost::multi_array<REAL, 1>
//void make_XYZ(const boost::array<int, 3> & pa, XYZ & output)
void make_XYZ(const boost::multi_array<REAL, 1> & pa, XYZ & output)
{
    output.x=pa[0];
    output.y=pa[1];
    output.z=pa[2];
}

void print_triangle_array(int count, TRIANGLE* tra)
{
    if(count > 0){
        cout << "Triangles: count: " << count;
        for(int ti=0; ti<count; ti++)
            cout << "("<<tra[ti] <<")"<< " ";
        cout << endl;
    }
}



int nsize = 20;
boost::array<int, 4> values_shape = {{ nsize, nsize, nsize }};
boost::multi_array<REAL, 3> values (values_shape);

//todo: Warning: may make a copy. (may call the copy constructor)

vector<TRIANGLE> make_grid()
//boost::multi_array<REAL, 3>& make_grid(const boost::multi_array<REAL, 3>& values)
{
    int nx = values.shape()[0];
    int ny = values.shape()[1];
    int nz = values.shape()[2];

    boost::array<int, 4> grid_shape = {{ nx, ny, nz, 3 }};
    boost::multi_array<REAL, 4> grid (grid_shape);

    // 8000 -> 99 -> 1325
    int expected_number_of_faces = (int)(sqrt(nx*ny*nz) * 15.)+10;
    //vector<TRIANGLE> result_triangles (expected_number_of_faces);
    vector<TRIANGLE> result_triangles;
    result_triangles.reserve(expected_number_of_faces);
    cout << "expected_number_of_faces: " << expected_number_of_faces << " ";
    cout << "grid total size: " << nx*ny*nz << " ";   //not the number of cells, but the number of nodes


    cout << endl;
    cout << endl;
    std::clog << "size: "     <<  result_triangles.size()     << '\n';
    std::clog << "capacity: " <<  result_triangles.capacity() << '\n';
    std::clog << "max_size: " <<  result_triangles.max_size() << '\n';


    //todo: send proper message if empty.
    //REAL grid_min = -1. -0.2;  //works
    REAL grid_min = -1. -0.2 -0.2 -0.1; //doesnt work
    //REAL grid_max = +1.;
    REAL grid_step = +0.2;

    //boost::multi_array<float, 2> verts;  =  {{0.1, 5.5, 0}, {0, 1, 0.5}, {1.5, 1, 0}, {1, 0, 0} };

    //float verts1[4][3] = {{0.1, 5.5, 0}, {0, 1, 0.5}, {1.5, 1, 0}, {1, 0, 0} };

    for(int di = 0; di < 3; di++)
    {
        for(int xi=0; xi < nx; xi++)
        for(int yi=0; yi < ny; yi++)
        for(int zi=0; zi < nz; zi++)
        {
            /*
            int ii;
            if(di==0){
                ii = xi;
            }else if(di==1){
                ii = yi;
            }else if(di==2){
                ii = zi;
            }
            */
            int ii = (di==0)?(xi):((di==1)?(yi):(zi));

            //not efficient
            REAL val = grid_min + grid_step * (REAL)ii;
            //REAL val = 0 + 1. * (REAL)ii;

            grid[xi][yi][zi][di] = val;
        }
    }
    //return values;

    if(ASSERTS){
        ;
    }


    for(int xi=0; xi < nx; xi++)
    for(int yi=0; yi < ny; yi++)
    for(int zi=0; zi < nz; zi++)
    {
        //78 msec version
        boost::multi_array<float, 1> c = grid[xi][yi][zi];
        //float f = 2.0 - (c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
        REAL x = c[0];
        REAL y = c[1];
        REAL z = c[2];
        //REAL f = 2.0 - (x*x+y*y+z*z);
        int shape = 1;
        if (shape == 1){
          REAL orb = exp(-abs(pow(z,2))*10*2)*5.+1.;
          REAL f = 2.0 - (pow(x,2) + pow(y,2) + pow(z,2))*orb;
          values[xi][yi][zi] = f;
          }
        else if(shape == 2){
          REAL f = pow(x,2) + pow(y,2) - pow(log(z + 3.20),2) - 0.02;
          values[xi][yi][zi] = f;
          }


    }

    cout << "first voxel:";
    for(int xi=0; xi < 2; xi++)
    for(int yi=0; yi < 2; yi++)
    for(int zi=0; zi < 2; zi++)
    {
      cout << values[xi][yi][zi] << " ";
    }
    cout << endl;

    int sign_counters[3]= {0,0,0};
    for(int xi=0; xi < nx; xi++)
    for(int yi=0; yi < ny; yi++)
    for(int zi=0; zi < nz; zi++)
    {
        REAL v = values[xi][yi][zi];
        int idx = v>0. ? 2 : (v<0. ? 1 : 0);  // (+) -> 2, (-) -> 1,  (0) -> 0
        sign_counters[idx]++;
    }
    cout << "zeros:"<< sign_counters[0] << "  (+):" << sign_counters[1] << " (-):" << sign_counters[2] << endl;


    for(int xi=0; xi < nx-1; xi++)
    for(int yi=0; yi < ny-1; yi++)
    for(int zi=0; zi < nz-1; zi++)
    {
        //XYZ x[2] = {XYZ(values[xi][yi][zi]), XYZ(values[xi+1][yi][zi])};
        XYZ p1;
        /*REAL pa[3] = grid[xi][yi][zi];
        p1.x=pa[0];
        p1.y=pa[1];
        p1.z=pa[2];*/
        make_XYZ(grid[xi][yi][zi], p1);

        /*
        XYZ p8[8];
        make_XYZ(grid[xi  ][yi  ][zi  ], p8[0]);
        make_XYZ(grid[xi+1][yi  ][zi  ], p8[1]);
        make_XYZ(grid[xi  ][yi+1][zi  ], p8[2]);
        make_XYZ(grid[xi+1][yi+1][zi  ], p8[3]);
        make_XYZ(grid[xi  ][yi  ][zi+1], p8[4]);
        make_XYZ(grid[xi+1][yi  ][zi+1], p8[5]);
        make_XYZ(grid[xi  ][yi+1][zi+1], p8[6]);
        make_XYZ(grid[xi+1][yi+1][zi+1], p8[7]);


        GRIDCELL g;
        make_XYZ (grid[xi  ][yi  ][zi  ], g.p[0] );
        make_XYZ (grid[xi+1][yi  ][zi  ], g.p[1] );
        make_XYZ (grid[xi  ][yi+1][zi  ], g.p[2] );
        make_XYZ (grid[xi+1][yi+1][zi  ], g.p[3] );
        make_XYZ (grid[xi  ][yi  ][zi+1], g.p[4] );
        make_XYZ (grid[xi+1][yi  ][zi+1], g.p[5] );
        make_XYZ (grid[xi  ][yi+1][zi+1], g.p[6] );
        make_XYZ (grid[xi+1][yi+1][zi+1], g.p[7] );
        //WRONG!!!
        */

        GRIDCELL g;

        make_XYZ (grid[xi  ][yi  ][zi  ], g.p[0] );
        make_XYZ (grid[xi+1][yi  ][zi  ], g.p[3] );
        make_XYZ (grid[xi  ][yi+1][zi  ], g.p[1] );
        make_XYZ (grid[xi+1][yi+1][zi  ], g.p[2] );
        make_XYZ (grid[xi  ][yi  ][zi+1], g.p[4] );
        make_XYZ (grid[xi+1][yi  ][zi+1], g.p[7] );
        make_XYZ (grid[xi  ][yi+1][zi+1], g.p[5] );
        make_XYZ (grid[xi+1][yi+1][zi+1], g.p[6] );


/*
        make_XYZV(grid[xi  ][yi  ][zi  ], g.p[0], g.val[0]);
        make_XYZV(grid[xi+1][yi  ][zi  ], g.p[1], g.val[0] );
        make_XYZV(grid[xi  ][yi+1][zi  ], g.p[2], g.val[0] );
        make_XYZV(grid[xi+1][yi+1][zi  ], g.p[3], g.val[0] );
        make_XYZV(grid[xi  ][yi  ][zi+1], g.p[4], g.val[0] );
        make_XYZV(grid[xi+1][yi  ][zi+1], g.p[5], g.val[0] );
        make_XYZV(grid[xi  ][yi+1][zi+1], g.p[6], g.val[0] );
        make_XYZV(grid[xi+1][yi+1][zi+1], g.p[7], g.val[0] );
*/
        g.val[0] = values[xi  ][yi  ][zi  ];
        g.val[3] = values[xi+1][yi  ][zi  ];
        g.val[1] = values[xi  ][yi+1][zi  ];
        g.val[2] = values[xi+1][yi+1][zi  ];
        g.val[4] = values[xi  ][yi  ][zi+1];
        g.val[7] = values[xi+1][yi  ][zi+1];
        g.val[5] = values[xi  ][yi+1][zi+1];
        g.val[6] = values[xi+1][yi+1][zi+1];

        TRIANGLE triangles[16];
        TRIANGLE*tra = triangles;
        REAL isolevel = 0.;
        //int Polygonise(const GRIDCELL& grid, REAL isolevel, TRIANGLE *triangles)
        int count = Polygonise(g, isolevel, tra, false);
        if (VERBOSE>2){
            print_triangle_array(count, tra);
        }
        if(count > 0){
            for(int ti=0; ti<count; ti++)
              result_triangles.push_back(triangles[ti]);
        }

        //return;

    }
    return result_triangles;
}


#define crossProduct_012(a,b,c) \
  (a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
  (a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
  (a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];

#define set_crossProduct_xyz(a,b,c) \
  (a).x = (b).y * (c).z - (c).y * (b).z; \
  (a).y = (b).z * (c).x - (c).z * (b).x; \
  (a).z = (b).x * (c).y - (c).x * (b).y;


#define set_subtract_xyz(a,b,c) \
  (a).x = (b).x - (c).x; \
  (a).y = (b).y - (c).y; \
  (a).z = (b).z - (c).z;

#define innerProduct_xyz(a, b) \
  (a).x * (b).x + \
  (a).y * (b).y + \
  (a).z * (b).z \
  ;

XYZ get_triangle_normal(TRIANGLE const& t){
    //XYZ a = t.p[0];
    //XYZ b = t.p[1];
    XYZ a,b;
    set_subtract_xyz(a, t.p[1], t.p[0]);
    set_subtract_xyz(b, t.p[2], t.p[0]);
    XYZ normal;
    //REAL Ax = b.x - a.x;
    set_crossProduct_xyz(normal, a, b);
    return normal;
    //fixme: reference or copy?
}
#define ASSERT(x)  {if(!(x)){cout<<"assertion error";}}

#define rnd() (rand() / (REAL)(RAND_MAX))


XYZ example_xyz(REAL x, REAL y, REAL z)
{
    XYZ a;
    a.x = x; a.y = y; a.z = z;
    return a;
}

void test_ext_prod()
{
    /*XYZ a,b;
    a.x = 1; a.y = 0; a.z = 0;
    b.x = 0; b.y = 1; b.z = 0;
    */
    XYZ a=example_xyz(1,0,0);
    XYZ b=example_xyz(0,1,0);
    XYZ normal;
    set_crossProduct_xyz(normal, a, b);
    cout << normal.x << "," << normal.y << "," << normal.z << endl;
}


int CUBE_MAP[8] = {0,3,1,2,4,7,5,6};

GRIDCELL make_grid_101()
{
  GRIDCELL g;
  int c = 0;
  for(int v3=0; v3<2; v3++)
  for(int v2=0; v2<2; v2++)
  for(int v1=0; v1<2; v1++)
  {
      XYZ v=example_xyz(v1,v2,v3);
      //g.p[c++] = v;
      g.p[CUBE_MAP[c]] = v;
      c++;
  }
  return g;
}

void test_gridcell1()
{
    //Tests a single grid cube (Marching) with one positive point in the corner. This should generate one single triangle.
    GRIDCELL g = make_grid_101();
    boost::array<REAL, 8> values =
      //{{ +10,-1, +10,-1,   -1,-1, -1,-1 }};
      {{ +2,-1, -1,-1,   -1,-1, -1,-1 }};
    for(int i=0; i<8; i++)
        g.val[i] = values[i];

    REAL isolevel = 0;
    TRIANGLE triangles[10];
    int count = Polygonise(g, isolevel, triangles, true);
    cout << " Count: " << count << endl;


    print_triangle_array(count, triangles);

    //todo: turn into a test
}


void test_gridcell2()
{
    //Tests a single grid cube (Marching) with two positive points. This should generate two single triangles that comprise a rectangle.
    GRIDCELL g = make_grid_101();
    cout << " --------------- t2 ---------------";
    boost::array<REAL, 8> values =
      //{{ +10,-1, +10,-1,   -1,-1, -1,-1 }};
      {{ +2,+2, -1,-1,   -1,-1, -1,-1 }};
    for(int i=0; i<8; i++)
        g.val[i] = values[i];

    REAL isolevel = 0;
    TRIANGLE triangles[10];
    int count = Polygonise(g, isolevel, triangles, true);
    cout << " Count: " << count << endl;

    print_triangle_array(count, triangles);

    //todo: turn into a test
    //The output should be similar to:
    // |0.666667,1,0|  |0.666667,0,0|  |0,0,0.666667| ) ( |0,1,0.666667|  |0.666667,1,0|  |0,0,0.666667| )
}

vf_t vector_to_vertsfaces(vector<TRIANGLE> const& ta)
{
    //pair<verts_t, faces_t> vf_t
    //verts_t
    //faces_t
    int nt = ta.size();
    int nv = ta.size()*3;

    boost::array<int, 2> v_shape = {{ nv, 3 }};
    boost::array<int, 2> f_shape = {{ nt, 3 }};
    boost::multi_array<REAL, 2> verts (v_shape);
    boost::multi_array<int, 2> faces (f_shape);

    /*
    const TRIANGLE& tr = ta[ti];

    paulbourke_mc1.cpp:698:15: error: binding of reference to type 'XYZ' to a value of type 'const XYZ' drops qualifiers
            XYZ & p0 = tr.p[0];
    */

    /*
    for(int ti=0; ti<nt; ti++){
        //inefficient
        TRIANGLE tr = ta[ti];
        XYZ & p0 = tr.p[0];
        //boost::array<REAL, 3> xxx = {p0.x, p0.y, p0.z};  // works here
        boost::array<int, 1> shape = {{ 3 }};
        boost::multi_array<REAL, 1> xxx(shape);
        xxx[0] = p0.x;
        xxx[1] = p0.y;
        xxx[2] = p0.z;
        verts[ti*3+0] = xxx;
        //verts[ti*3+0][0] = tr.p[0].x;
        //verts[ti*3+0][0] = tr.p[0].x;
    }
    */

    int FLIP[3] = {0, 2, 1};
    for(int ti=0; ti<nt; ti++)
    {
        //inefficient
        TRIANGLE tr = ta[ti];
        XYZ n = get_triangle_normal(tr);
        REAL sgn = innerProduct_xyz(n, tr.p[0]); //todo: centroid
        bool flip_verts =
              //0;
              //rnd() > 0.5;
              (sgn<0);
        for(int side=0; side<3; side++){
            int side2 = side;
            if(flip_verts){
                 side2 = FLIP[side];
            }
            else
                 side2 = side;

            const REAL NOISE_LEVEL= 0.1*0;
            verts[ti*3+side][0] = tr.p[side2].x+NOISE_LEVEL*rnd();
            verts[ti*3+side][1] = tr.p[side2].y+NOISE_LEVEL*rnd();
            verts[ti*3+side][2] = tr.p[side2].z+NOISE_LEVEL*rnd();
        }
    }
    ASSERT(nt*3 == verts.shape()[0]);

    for(int ti=0; ti<nt; ti++){
      bool flip = 0; //!(ti % 2);
      faces[ti][0] = ti*3 + 0;
      if(flip){
        faces[ti][1] = ti*3 + 2;
        faces[ti][2] = ti*3 + 1;
      }
      else{
        faces[ti][1] = ti*3 + 1;
        faces[ti][2] = ti*3 + 2;
      }
    }

    //Twitch
    for (int j=0;j<nt; j+= 1){
        boost::array<REAL, 3> centroid = {0,0,0};

        for (int s=0;s<3; s+= 1){
            for (int d=0;d<3; d+= 1){
                centroid[d] += verts[j*3+s][d] / 3.;
            }
        }
        /*
        {
        REAL dc = rnd()*0.1;
        for (int d=0;d<3; d+= 1){
            centroid[d] = dc;
        }
        }
        */

        // separation ! (delay) ==> accumulation!

        const REAL alpha = 0.90;
        const REAL beta = 1. - alpha;
        for (int s=0;s<3; s+= 1){
            for (int d=0;d<3; d+= 1){
                verts[j*3+s][d]  = verts[j*3+s][d] * alpha + centroid[d] * beta;
            }
        }
    }

    const bool ADD_THE_SPIKES = false;
    if(ADD_THE_SPIKES){
        // Add the spikes !
        int spike_triangle=100;
        for (int j=0;j<nt*3; j += 100){
            spike_triangle = j;
            verts[spike_triangle][0] *= 2;
            verts[spike_triangle][1] *= 2;
            verts[spike_triangle][2] *= 2;
        }
    }

    vf_t p2 = make_pair(verts, faces);
    return p2;
}



extern "C" {
    void make_object(float* verts, int *nv, int* faces, int *nf);
    int main();
}

#include "timer.hpp"

void make_object(float* verts, int *nv, int* faces, int *nf){

    timer timr;

    vector<TRIANGLE> ta = make_grid();
    timr.stop("make grid"); // 1789.14 1565.37 1443.99 msec

    //very inefficient
    vf_t vf = vector_to_vertsfaces(ta);
    timr.stop("vector_to_vertsfaces()");

    *nv = vf.first.shape()[0];
    *nf = vf.second.shape()[0];

    for(int vi=0; vi<*nv; vi++){
        for(int di=0; di<3; di++){
            verts[vi*3+di] = vf.first[vi][di];
        }
      }

    for(int fi=0; fi<*nf; fi++){
        for(int si=0; si<3; si++){
            faces[fi*3+si] = vf.second[fi][si];
        }
      }
    timr.stop("finished");
}

int main0()
{
    //boost::multi_array<float, 4> grid;
    //make_grid(grid, values)
    //boost::multi_array<float, 3> values = make_grid();
    vector<TRIANGLE> ta = make_grid();
    cout << endl;
    cout << endl;

    cout << "number of triangles: " << ta.size();
    cout << endl;

    std::clog << "size: "     <<  ta.size()     << '\n';
    std::clog << "capacity: " <<  ta.capacity() << '\n';
    std::clog << "max_size: " <<  ta.max_size() << '\n';

    //very inefficient
    vf_t vf = vector_to_vertsfaces(ta);

    return 0;
}


int main()
{
     test_ext_prod();
     if(false)
        test_gridcell1();
     test_gridcell2();

     cout << "staying alive...";
     return 0;
}


//CODE THAT IS NOT USED YET:
//IN PROGRESS

#include <tuple>

vf_t triangles_to_vertsfaces(vector<TRIANGLE> const& ta)
{
    //converts the expanded into condenced format.
    //extracts unique points


    int nt = ta.size();
    int nv = ta.size()*3;

    boost::array<int, 2> v_shape = {{ nv, 3 }};
    boost::array<int, 2> f_shape = {{ nt, 3 }};
    boost::multi_array<REAL, 2> verts (v_shape);
    boost::multi_array<int, 2> faces (f_shape);

    /*
    const TRIANGLE& tr = ta[ti];

    paulbourke_mc1.cpp:698:15: error: binding of reference to type 'XYZ' to a value of type 'const XYZ' drops qualifiers
            XYZ & p0 = tr.p[0];
    */

    /*
    for(int ti=0; ti<nt; ti++){
        //inefficient
        TRIANGLE tr = ta[ti];
        XYZ & p0 = tr.p[0];
        //boost::array<REAL, 3> xxx = {p0.x, p0.y, p0.z};  // works here
        boost::array<int, 1> shape = {{ 3 }};
        boost::multi_array<REAL, 1> xxx(shape);
        xxx[0] = p0.x;
        xxx[1] = p0.y;
        xxx[2] = p0.z;
        verts[ti*3+0] = xxx;
        //verts[ti*3+0][0] = tr.p[0].x;
        //verts[ti*3+0][0] = tr.p[0].x;
    }
    */

    int FLIP[3] = {0, 2, 1};
    for(int ti=0; ti<nt; ti++)
    {
        //inefficient
        TRIANGLE tr = ta[ti];
        XYZ n = get_triangle_normal(tr);
        REAL sgn = innerProduct_xyz(n, tr.p[0]); //todo: centroid
        bool flip_verts =
              //0;
              //rnd() > 0.5;
              (sgn<0);
        for(int side=0; side<3; side++){
            int side2 = side;
            if(flip_verts){
                 side2 = FLIP[side];
            }
            else
                 side2 = side;

            const REAL NOISE_LEVEL= 0.1*0;
            verts[ti*3+side][0] = tr.p[side2].x+NOISE_LEVEL*rnd();
            verts[ti*3+side][1] = tr.p[side2].y+NOISE_LEVEL*rnd();
            verts[ti*3+side][2] = tr.p[side2].z+NOISE_LEVEL*rnd();
        }
    }
    ASSERT(nt*3 == verts.shape()[0]);

    for(int ti=0; ti<nt; ti++){
      bool flip = 0; //!(ti % 2);
      faces[ti][0] = ti*3 + 0;
      if(flip){
        faces[ti][1] = ti*3 + 2;
        faces[ti][2] = ti*3 + 1;
      }
      else{
        faces[ti][1] = ti*3 + 1;
        faces[ti][2] = ti*3 + 2;
      }
    }

    //Twitch
    for (int j=0;j<nt; j+= 1){
        boost::array<REAL, 3> centroid = {0,0,0};

        for (int s=0;s<3; s+= 1){
            for (int d=0;d<3; d+= 1){
                centroid[d] += verts[j*3+s][d] / 3.;
            }
        }
        /*
        {
        REAL dc = rnd()*0.1;
        for (int d=0;d<3; d+= 1){
            centroid[d] = dc;
        }
        }
        */

        // separation ! (delay) ==> accumulation!

        const REAL alpha = 0.90;
        const REAL beta = 1. - alpha;
        for (int s=0;s<3; s+= 1){
            for (int d=0;d<3; d+= 1){
                verts[j*3+s][d]  = verts[j*3+s][d] * alpha + centroid[d] * beta;
            }
        }
    }
}
