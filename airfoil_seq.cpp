/*
 * Open source copyright declaration based on BSD open source template:
 * http://www.opensource.org/licenses/bsd-license.php
 *
 * This file is part of the OP2 distribution.
 *
 * Copyright (c) 2011, Mike Giles and others. Please see the AUTHORS file in
 * the main source directory for a full list of copyright holders.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The name of Mike Giles may not be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY Mike Giles ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL Mike Giles BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//
//     Nonlinear airfoil lift calculation
//
//     Written by Mike Giles, 2010-2011, based on FORTRAN code
//     by Devendra Ghate and Mike Giles, 2005
//

//
// standard headers
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <vector>

typedef double real;
// global constants

real gam, gm1, cfl, eps, mach, alpha, qinf[4];


int main(int argc, char **argv)
{

  int    *becell, *ecell,  *bound, *bedge, *edge, *cell;
  real  *x, *q, *qold, *adt, *res;

  int    nnode,ncell,nedge,nbedge,niter;
  real  rms;

  //timer
  double wall_t1, wall_t2;

  // read in grid

  printf("reading in grid \n");

  FILE *fp;
  if ( (fp = fopen("./new_grid.dat","r")) == NULL) {
    printf("can't open file new_grid.dat\n"); exit(-1);
  }

  if (fscanf(fp,"%d %d %d %d \n",&nnode, &ncell, &nedge, &nbedge) != 4) {
    printf("error reading from new_grid.dat\n"); exit(-1);
  }

  cell   = (int *) malloc(4*ncell*sizeof(int));
  edge   = (int *) malloc(2*nedge*sizeof(int));
  ecell  = (int *) malloc(2*nedge*sizeof(int));
  bedge  = (int *) malloc(2*nbedge*sizeof(int));
  becell = (int *) malloc(  nbedge*sizeof(int));
  bound  = (int *) malloc(  nbedge*sizeof(int));

  x      = (real *) malloc(2*nnode*sizeof(real));
  q      = (real *) malloc(4*ncell*sizeof(real));
  qold   = (real *) malloc(4*ncell*sizeof(real));
  res    = (real *) malloc(4*ncell*sizeof(real));
  adt    = (real *) malloc(  ncell*sizeof(real));

  //set variables for graph coloring
  int* cell2edge = (int *) malloc(4*ncell*sizeof(int));
  for(int i = 0; i< 4*ncell; ++i){
	  cell2edge[i] = -1;
  }
  int* edge2color = (int*) malloc(nedge*sizeof(int));
  for(int i = 0; i< nedge; ++i){
  	  edge2color[i] = -1;
   }
  for (int n=0; n<nnode; n++) {
    if (fscanf(fp,"%lf %lf \n",&x[2*n], &x[2*n+1]) != 2) {
      printf("error reading from new_grid.dat\n"); exit(-1);
    }
  }

  for (int n=0; n<ncell; n++) {
    if (fscanf(fp,"%d %d %d %d \n",&cell[4*n  ], &cell[4*n+1],
                                   &cell[4*n+2], &cell[4*n+3]) != 4) {
      printf("error reading from new_grid.dat\n"); exit(-1);
    }
  }

  for (int n=0; n<nedge; n++) {
    if (fscanf(fp,"%d %d %d %d \n",&edge[2*n], &edge[2*n+1],
                                   &ecell[2*n],&ecell[2*n+1]) != 4) {
      printf("error reading from new_grid.dat\n"); exit(-1);
    }
    for(int i = 0; i < 4;++i){
    	if(cell2edge[ecell[2*n]+i] == -1){
    		cell2edge[ecell[2*n]+i] = n;
    		break;
    	}
    }
    for(int i = 0; i < 4;++i){
       	if(cell2edge[ecell[2*n+1]+i] == -1){
            cell2edge[ecell[2*n+1]+i] = n;
       		break;
    	}
    }
  }

  for (int n=0; n<nbedge; n++) {
    if (fscanf(fp,"%d %d %d %d \n",&bedge[2*n],&bedge[2*n+1],
                                   &becell[n], &bound[n]) != 4) {
      printf("error reading from new_grid.dat\n"); exit(-1);
    }
  }
  fclose(fp);
  int max = 0;
  std::vector<std::vector<int> > color2edge;
  for(int edge_ind = 0; edge_ind < nedge; ++edge_ind){
	  int color = 0;
	  while(1){
		  bool valid_color = true;
		  for(int i = 0; i < 4; ++i){
			  if(edge2color[cell2edge[ecell[edge_ind]+i]] == color ||
					  edge2color[cell2edge[ecell[edge_ind+1]+i]] == color){
				  valid_color = false;
			  }
		  }
		  if(valid_color){
			  edge2color[edge_ind] = color;
			  if(color2edge.size() == color){
				  color2edge.push_back(std::vector<int>(1,edge_ind));
			  } else if(color < color2edge.size()){
				  color2edge[color].push_back(edge_ind);
			  } else {
				  printf("ismet para van\n");
			  }
			  if(color > max){
				  max = color;
				  printf("%d\n",max);
			  }
			  break;
		  }
		  ++color;
	  }
  }
  printf("%d %d %d %d %d %d", color2edge[0].size(),color2edge[1].size(),color2edge[2].size(),
		  color2edge[3].size(),color2edge[4].size(),color2edge[5].size());
  // set constants and initialise flow field and residual

  printf("initialising flow field \n");

  gam = 1.4f;
  gm1 = gam - 1.0f;
  cfl = 0.9f;
  eps = 0.05f;

  real mach  = 0.4f;
  real alpha = 3.0f*atan(1.0f)/45.0f;
  real p     = 1.0f;
  real r     = 1.0f;
  real u     = sqrt(gam*p/r)*mach;
  real e     = p/(r*gm1) + 0.5f*u*u;

  qinf[0] = r;
  qinf[1] = r*u;
  qinf[2] = 0.0f;
  qinf[3] = r*e;

  for (int n=0; n<ncell; n++) {
    for (int m=0; m<4; m++) {
        q[4*n+m] = qinf[m];
      res[4*n+m] = 0.0f;
    }
  }

  //initialise timers for total execution wall time
  wall_t1 = omp_get_wtime();

  // main time-marching loop

  niter = 1000;
  double save = 0, area = 0, update = 0, flux_res = 0, perem = 0, wall_t_a = 0, wall_t_b = 0;
  #pragma omp parallel
  {
  for(int iter=1; iter<=niter; iter++) {
	#pragma omp single
	{
		wall_t_b = omp_get_wtime();
	}
    // save old flow solution
	#pragma omp for
    for (int i = 0; i < ncell; i++) {
      for (int n=0; n<4; n++) qold[4*i+n] = q[4*i+n];
    }
	#pragma omp single
    {
    	wall_t_a = omp_get_wtime();
    	save += wall_t_a - wall_t_b;
    }
    // predictor/corrector update loop

    for(int k=0; k<2; k++) {

		#pragma omp single
    	{
    		wall_t_b = omp_get_wtime();
    	}
      // calculate area/timstep
	  #pragma omp for
      for (int i = 0; i < ncell; i++) {
        double dx,dy, ri,u,v,c;

        ri =  1.0f/q[4*i+0];
        u  =   ri*q[4*i+1];
        v  =   ri*q[4*i+2];
        c  = sqrt(gam*gm1*(ri*q[4*i+3]-0.5f*(u*u+v*v)));

        dx = x[cell[4*i+1]*2+0] - x[cell[4*i+0]*2+0];
        dy = x[cell[4*i+1]*2+1] - x[cell[4*i+0]*2+1];
        adt[i]  = fabs(u*dy-v*dx) + c*sqrt(dx*dx+dy*dy);

        dx = x[cell[4*i+2]*2+0] - x[cell[4*i+1]*2+0];
        dy = x[cell[4*i+2]*2+1] - x[cell[4*i+1]*2+1];
        adt[i] += fabs(u*dy-v*dx) + c*sqrt(dx*dx+dy*dy);

        dx = x[cell[4*i+3]*2+0] - x[cell[4*i+2]*2+0];
        dy = x[cell[4*i+3]*2+1] - x[cell[4*i+2]*2+1];
        adt[i] += fabs(u*dy-v*dx) + c*sqrt(dx*dx+dy*dy);

        dx = x[cell[4*i+0]*2+0] - x[cell[4*i+3]*2+0];
        dy = x[cell[4*i+0]*2+1] - x[cell[4*i+3]*2+1];
        adt[i] += fabs(u*dy-v*dx) + c*sqrt(dx*dx+dy*dy);

        adt[i] = adt[i] / cfl;

      }
	  #pragma omp single
      {
        wall_t_a = omp_get_wtime();
        area += (wall_t_a - wall_t_b);
      // calculate flux residual
      wall_t_b = omp_get_wtime();
      }

      for(size_t color = 0; color < color2edge.size(); ++color){
    	  #pragma omp for
    	  for(size_t ind = 0; ind < color2edge[color].size(); ++ind){
    		int i = color2edge[color][ind];
    		double dx,dy,mu, ri, p1,vol1, p2,vol2, f;

    		dx = x[edge[2*i+0]*2+0] - x[edge[2*i+1]*2+0];
        	dy = x[edge[2*i+0]*2+1] - x[edge[2*i+1]*2+1];

        	ri   = 1.0f/q[ecell[2*i+0]*4+0];
        	p1   = gm1*(q[ecell[2*i+0]*4+3]-0.5f*ri*(q[ecell[2*i+0]*4+1]*q[ecell[2*i+0]*4+1]+q[ecell[2*i+0]*4+2]*q[ecell[2*i+0]*4+2]));
        	vol1 =  ri*(q[ecell[2*i+0]*4+1]*dy - q[ecell[2*i+0]*4+2]*dx);

        	ri   = 1.0f/q[ecell[2*i+1]*4+0];
        	p2   = gm1*(q[ecell[2*i+1]*4+3]-0.5f*ri*(q[ecell[2*i+1]*4+1]*q[ecell[2*i+1]*4+1]+q[ecell[2*i+1]*4+2]*q[ecell[2*i+1]*4+2]));
        	vol2 =  ri*(q[ecell[2*i+1]*4+1]*dy - q[ecell[2*i+1]*4+2]*dx);

        	mu = 0.5f*(adt[ecell[2*i+0]]+adt[ecell[2*i+1]])*eps;

        	f = 0.5f*(vol1* q[ecell[2*i+0]*4+0]         + vol2* q[ecell[2*i+1]*4+0]        ) + mu*(q[ecell[2*i+0]*4+0]-q[ecell[2*i+1]*4+0]);
        	res[ecell[2*i+0]*4+0] += f;
        	res[ecell[2*i+1]*4+0] -= f;
        	f = 0.5f*(vol1* q[ecell[2*i+0]*4+1] + p1*dy + vol2* q[ecell[2*i+1]*4+1] + p2*dy) + mu*(q[ecell[2*i+0]*4+1]-q[ecell[2*i+1]*4+1]);
        	res[ecell[2*i+0]*4+1] += f;
        	res[ecell[2*i+1]*4+1] -= f;
        	f = 0.5f*(vol1* q[ecell[2*i+0]*4+2] - p1*dx + vol2* q[ecell[2*i+1]*4+2] - p2*dx) + mu*(q[ecell[2*i+0]*4+2]-q[ecell[2*i+1]*4+2]);
        	res[ecell[2*i+0]*4+2] += f;
        	res[ecell[2*i+1]*4+2] -= f;
        	f = 0.5f*(vol1*(q[ecell[2*i+0]*4+3]+p1)     + vol2*(q[ecell[2*i+1]*4+3]+p2)    ) + mu*(q[ecell[2*i+0]*4+3]-q[ecell[2*i+1]*4+3]);
        	res[ecell[2*i+0]*4+3] += f;
        	res[ecell[2*i+1]*4+3] -= f;
    	  }
      }
      //for (int i = 0; i < nedge; i++) {

      //}

	  #pragma omp single
      {
      wall_t_a = omp_get_wtime();
      flux_res += (wall_t_a - wall_t_b);
      // Apply boundary conditions
      wall_t_b = omp_get_wtime();

      for (int i = 0; i < nbedge; i++) {
        double dx,dy,mu, ri, p1,vol1, p2,vol2, f;

        dx = x[bedge[2*i+0]*2+0] - x[bedge[2*i+1]*2+0];
        dy = x[bedge[2*i+0]*2+1] - x[bedge[2*i+1]*2+1];

        ri = 1.0f/q[becell[i]*4+0];
        p1 = gm1*(q[becell[i]*4+3]-0.5f*ri*(q[becell[i]*4+1]*q[becell[i]*4+1]+q[becell[i]*4+2]*q[becell[i]*4+2]));

        if (bound[i]==1) { //Far-field
          res[becell[i]*4+1] += + p1*dy;
          res[becell[i]*4+2] += - p1*dx;
        }
        else {
          vol1 =  ri*(q[becell[i]*4+1]*dy - q[becell[i]*4+2]*dx);

          ri   = 1.0f/qinf[0];
          p2   = gm1*(qinf[3]-0.5f*ri*(qinf[1]*qinf[1]+qinf[2]*qinf[2]));
          vol2 =  ri*(qinf[1]*dy - qinf[2]*dx);

          mu = adt[becell[i]]*eps;

          f = 0.5f*(vol1* q[becell[i]*4+0]         + vol2* qinf[0]        ) + mu*(q[becell[i]*4+0]-qinf[0]);
          res[becell[i]*4+0] += f;
          f = 0.5f*(vol1* q[becell[i]*4+1] + p1*dy + vol2* qinf[1] + p2*dy) + mu*(q[becell[i]*4+1]-qinf[1]);
          res[becell[i]*4+1] += f;
          f = 0.5f*(vol1* q[becell[i]*4+2] - p1*dx + vol2* qinf[2] - p2*dx) + mu*(q[becell[i]*4+2]-qinf[2]);
          res[becell[i]*4+2] += f;
          f = 0.5f*(vol1*(q[becell[i]*4+3]+p1)     + vol2*(qinf[3]+p2)    ) + mu*(q[becell[i]*4+3]-qinf[3]);
          res[becell[i]*4+3] += f;
        }
      }
      wall_t_a = omp_get_wtime();
      perem += (wall_t_a - wall_t_b);

      // update flow field
      wall_t_b = omp_get_wtime();
      rms = 0.0;
      }
      #pragma omp for reduction(+:rms)
      for (int i = 0; i < ncell; i++) {
        double del, adti;

        adti = 1.0f/adt[i];

        for (int n=0; n<4; n++) {
          del    = adti*res[4*i+n];
          q[4*i+n]   = qold[4*i+n] - del;
          res[4*i+n] = 0.0f;
          rms  += del*del;
        }
      }
	#pragma omp single
      {
      wall_t_a = omp_get_wtime();
      update += (wall_t_a - wall_t_b);
      }

    }

    // print iteration history
	#pragma omp single
    {
    rms = sqrt(rms/(double) ncell);
    if (iter%100 == 0){
    	printf(" %d  %10.5e \n",iter,rms);
    	printf("\tsave: %f\n",save);
    	printf("\tarea: %f\n",area);
    	printf("\tflux_res: %f\n",flux_res);
    	printf("\tperem: %f\n",perem);
    	printf("\tupdate: %f\n",update);
    }
    }
  }
  }
  wall_t2 = omp_get_wtime();

  printf("Max total runtime = \n%f\n",wall_t2-wall_t1);


  free(cell);
  free(edge);
  free(ecell);
  free(bedge);
  free(becell);
  free(bound);
  free(x);
  free(q);
  free(qold);
  free(res);
  free(adt);


  free(cell2edge);
  free(edge2color);
}

