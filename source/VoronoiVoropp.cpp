#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <iostream>
#include "voro++.hh"
using namespace voro;
using namespace std;



extern "C" {
  void voropp_generate_voronoi_tessellation(int, int, int, int *, int *, 
					    int *, double *, double *, 
					    double *, double *);
}

extern "C" {
  int voropp_find_voronoi_cell(double *, double *, double *);
}


extern "C" {
  void voropp_delete_container(void);
}


container *con;



//=============================================================================
// voropp_delete_container
// ..
//=============================================================================
void voropp_delete_container(void)
{
  printf("Deleting container\n");
  delete con;
  printf("Finished deleting container\n");
}



//=============================================================================
// voropp_generate_voronoi_tessellation
// ..
//=============================================================================
void voropp_generate_voronoi_tessellation
(int ndimaux,                           // Dimensionality
 int Naux,                              // No. of points/vertices
 int Nlinesmax,                         // Max. allowed number of connections
 int *Nlines,                           // No. of lines/connections
 int *i1,                               // List of connection vertices (1)
 int *i2,                               // "" (2)
 double *r,                             // Positions of vertices
 double *vol,                           // ..
 double *boxmin,                        // ..
 double *boxmax)                        // ..
{
  int dim = ndimaux;                    // Dimensionality
  int i;                                // Counter
  int ii;                               // Aux. facet counter
  int j;                                // Counter
  int jj;                               // Aux. facet counter
  int k;                                // Dimension counter
  int NauxLoc;
  int numpoints = Naux;                 // No. of vertices
  int plist[3];                         // List of vertices on facet
  double rmin[3];                       // Minimum extent of bounding box
  double rmax[3];                       // Maximum extent of bounding box
  double boxsize[3];                    // ..
  double cvol = 1.0;                    // ..

  int Nlinesaux = 0;

  printf("[voropp_generate_voronoi_tessellation]\n");
  printf("ndim : %d\n",ndimaux);
  printf("Naux : %d\n",Naux);
  printf("Nlinesmax : %d\n",Nlinesmax);

  *Nlines = 0;

  for (i=0; i<numpoints; i++) {
    for (k=0; k<dim; k++) {
      rmin[k] = std::min((double) r[dim*i + k],rmin[k]);
      rmax[k] = std::max((double) r[dim*i + k],rmax[k]);
    }
  }
  for (k=0; k<dim; k++) {
    boxsize[k] = rmax[k] - rmin[k];
    //rmin[k] -= 0.25*boxsize[k];
    //rmax[k] += 0.25*boxsize[k];
    rmin[k] = boxmin[k];
    rmax[k] = boxmax[k];
    cvol *= (rmax[k] - rmin[k]);
    cout << "Box : " << k << "   " << rmin[k] << "    " << rmax[k] << endl;
  }

  // Create object container
  //container con(rmin[0],rmax[0],rmin[1],rmax[1],rmin[2],rmax[2],
  //	1,1,1,false,false,false,Naux);
  //con = new container(rmin[0],rmax[0],rmin[1],rmax[1],rmin[2],rmax[2],
  //	      8,8,8,false,false,false,Naux);

  NauxLoc = Naux/(32*32*32);

  //  con = new container(rmin[0],rmax[0],rmin[1],rmax[1],rmin[2],rmax[2],
  //		      32,32,32,false,false,false,Naux);

  

  con = new container(rmin[0],rmax[0],rmin[1],rmax[1],rmin[2],rmax[2],
  		      32,32,32,false,false,false,NauxLoc);
  con->clear();

  // Add particles to container
  for (i=0; i<numpoints; i++) {
    con->put(i,r[dim*i],r[dim*i + 1],r[dim*i + 2]);
  }

  // Loop over all cells/particles
  c_loop_all vl(*con);
  voronoicell_neighbor c;
  int ijk,q,Nneib; double *pp;


  // Now loop over all particles inside cell
  //---------------------------------------------------------------------------
  if (vl.start()) do if(con->compute_cell(c,vl)) {
    
    std::vector<int> vi;
    ijk = vl.ijk; 
    q = vl.q; 
    pp = con->p[ijk] + con->ps*q;
    int icurrent = con->id[ijk][q];

    if (icurrent >= numpoints) {
      printf("Exceeded maximum allowed number of points (somehow) : %d %d\n",
	     icurrent,numpoints);
      exit(0);
    }

    Nneib = c.number_of_faces();
    c.neighbors(vi);


    for (std::vector<int>::iterator it = vi.begin() ; it != vi.end(); ++it) {
      if (*it >= 0) {

	if (*it >= numpoints) {
	  printf("Exceeded maximum allowed number of points (somehow2) : %d %d\n",
		 *it,numpoints);
	  exit(0);
	}

	// Check we've not filled memory buffer with connections
	if (Nlinesaux >= Nlinesmax) {
	  printf("Exceeded maximum allowed number of lines\n");
	  exit(0);
	}
	i1[Nlinesaux] = icurrent;
	i2[Nlinesaux] = *it;
	Nlinesaux++;
      }
    }
    vol[icurrent] = c.volume();


  } while(vl.inc());
  //---------------------------------------------------------------------------


  // Output the particle positions in gnuplot format
  con->draw_particles("points.gnu");
  
  // Output the Voronoi cells in gnuplot format
  con->draw_cells_gnuplot("vcells.gnu");

  double vvol=con->sum_cell_volumes();
  printf("Container volume : %g\n"
	 "Voronoi volume   : %g\n"
	 "Difference       : %g\n",cvol,vvol,vvol-cvol);

  *Nlines = Nlinesaux;

  return;
}



//=============================================================================
// voropp_find_voronoi_cell
// ..
//=============================================================================
int voropp_find_voronoi_cell
(double *x,
 double *y,
 double *z)
{
  bool cell_found;
  double xv,yv,zv;
  int iVoronoi;

  cout << "Searching for cell with point : " << *x << "   " << *y << "   " 
       << *z << "   " << endl;

  cell_found = con->find_voronoi_cell(*x, *y, *z, xv, yv, zv, iVoronoi);
				      

  cout << "Found in cell : " << iVoronoi << endl;

  if (!cell_found) iVoronoi = -1;

  return iVoronoi;
}
