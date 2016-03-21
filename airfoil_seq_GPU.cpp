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
#include "writeVTK.h"
typedef float real;

// global constants

real gam, gm1, cfl, eps, mach, alpha, qinf[4];

void init(int &nnode, int &ncell, int &nedge, int &nbedge, int* &becell,
		int* &ecell, int* &bound, int* &bedge, int* &edge, int* &cell, real* &x,
		real* &q, real* &qold, real* &adt, real* &res) {
	FILE *fp;
	if ((fp = fopen("./new_grid.dat", "r")) == NULL) {
		printf("can't open file new_grid.dat\n");
		exit(-1);
	}

	if (fscanf(fp, "%d %d %d %d \n", &nnode, &ncell, &nedge, &nbedge) != 4) {
		printf("error reading from new_grid.dat\n");
		exit(-1);
	}
	cell = (int *) malloc(4 * ncell * sizeof(int));
	edge = (int *) malloc(2 * nedge * sizeof(int));
	ecell = (int *) malloc(2 * nedge * sizeof(int));
	bedge = (int *) malloc(2 * nbedge * sizeof(int));
	becell = (int *) malloc(nbedge * sizeof(int));
	bound = (int *) malloc(nbedge * sizeof(int));

	x = (real *) malloc(2 * nnode * sizeof(real));
	q = (real *) malloc(4 * ncell * sizeof(real));
	qold = (real *) malloc(4 * ncell * sizeof(real));
	res = (real *) malloc(4 * ncell * sizeof(real));
	adt = (real *) malloc(ncell * sizeof(real));

	for (int n = 0; n < nnode; n++) {
		if (fscanf(fp, "%f %f \n", &x[2 * n], &x[2 * n + 1]) != 2) {
			printf("error reading from new_grid.dat\n");
			exit(-1);
		}
	}
	for (int n = 0; n < ncell; n++) {
		if (fscanf(fp, "%d %d %d %d \n", &cell[4 * n], &cell[4 * n + 1],
				&cell[4 * n + 2], &cell[4 * n + 3]) != 4) {
			printf("error reading from new_grid.dat\n");
			exit(-1);
		}
	}
	for (int n = 0; n < nedge; n++) {
		if (fscanf(fp, "%d %d %d %d \n", &edge[2 * n], &edge[2 * n + 1],
				&ecell[2 * n], &ecell[2 * n + 1]) != 4) {
			printf("error reading from new_grid.dat\n");
			exit(-1);
		}
	}
	for (int n = 0; n < nbedge; n++) {
		if (fscanf(fp, "%d %d %d %d \n", &bedge[2 * n], &bedge[2 * n + 1],
				&becell[n], &bound[n]) != 4) {
			printf("error reading from new_grid.dat\n");
			exit(-1);
		}
	}
	fclose(fp);
	// set constants and initialise flow field and residual
	printf("initialising flow field \n");

	gam = 1.4f;
	gm1 = gam - 1.0f;
	cfl = 0.9f;
	eps = 0.05f;

	real mach = 0.4f;
	real alpha = 3.0f * atan(1.0f) / 45.0f;
	real p = 1.0f;
	real r = 1.0f;
	real u = sqrt(gam * p / r) * mach;
	real e = p / (r * gm1) + 0.5f * u * u;

	qinf[0] = r;
	qinf[1] = r * u;
	qinf[2] = 0.0f;
	qinf[3] = r * e;

	for (int n = 0; n < ncell; n++) {
		for (int m = 0; m < 4; m++) {
			q[4 * n + m] = qinf[m];
			res[4 * n + m] = 0.0f;
		}
	}
}

void coloring(int ncell, int* edge2cell, int nedge, int cellnumbyedge,
		std::vector<std::vector<int> > &color2edge) {
	printf("coloring with ncell: %d nedge: %d, cellsbyedges: %d\n", ncell, nedge, cellnumbyedge);
	//set variables for graph coloring
	int* cell2edge = (int *) malloc(4 * ncell * sizeof(int));//ind: 4*cellind+i val:edge next to cell
	int* edge2color = (int*) malloc(nedge * sizeof(int));//ind: edgeind val: color
	for (int i = 0; i < 4 * ncell; ++i) {
		cell2edge[i] = -1;
	}
	for (int i = 0; i < nedge; ++i) {
		edge2color[i] = -1;
	}
	//parse cells to edges
	for (unsigned n = 0; n < nedge; ++n) {
		for(unsigned cellnum = 0; cellnum < cellnumbyedge; cellnum++ ){
			bool placed = false;
			for (int i = 0; i < 4 && !placed; ++i) {
				if (cell2edge[4*edge2cell[cellnumbyedge * n + cellnum] + i] == -1) {
					cell2edge[4*edge2cell[cellnumbyedge * n + cellnum] + i] = n;
					placed = true;
				}
			}
		}
	}

	int max = 0;
	for (int edge_ind = 0; edge_ind < nedge; ++edge_ind) {
		int color = 0;
		while (true) {//this loop will break if we find the right color
			bool valid_color = true;
			for(int cellnum = 0; cellnum < cellnumbyedge; ++cellnum){
				for (int i = 0; i < 4; ++i) {
					/*printf("edge: %d, cellnum: %d, edge2cell: %d, cell2edge: %d, edge2color: %d\n",edge_ind,cellnum,
							edge2cell[cellnumbyedge * edge_ind + cellnum],
							cell2edge[4*edge2cell[cellnumbyedge * edge_ind + cellnum] + i],
							edge2color[cell2edge[4*edge2cell[cellnumbyedge * edge_ind + cellnum] + i]]);*/
					if (cell2edge[4*edge2cell[cellnumbyedge * edge_ind + cellnum] + i] != -1
							&& edge2color[cell2edge[4*edge2cell[cellnumbyedge * edge_ind + cellnum] + i]] == color) {
						valid_color = false;
					}
				}
			}

			if (valid_color) {
				edge2color[edge_ind] = color;
				if (color2edge.size() == color) {
					color2edge.push_back(std::vector<int>(1, edge_ind));
				} else if (color < color2edge.size()) {
					(color2edge[color]).push_back(edge_ind);
				} else {
					printf("edgeind:%d, nedge: %d color:%d\n",edge_ind, nedge, color);
					exit(1);
					printf("ismet para van\n");
				}
				if (color > max) {
					max = color;
					//printf("%d\n", max);
				}
				break;
			}
			++color;
		}
	}

	for(int i = 0; i<color2edge.size(); ++i){
		printf("Number of edges with color %d: %d\n",i, color2edge[i].size());
	}
	free(cell2edge);
	free(edge2color);
}

void celledgestatistic(int ncell,int nedge,int* ecell){
	std::vector<int> cells(ncell,0);
	for(int edgeid = 0; edgeid < nedge; edgeid++){
		int cellid = ecell[2*edgeid];
		cells[cellid]++;
		cellid = ecell[2*edgeid+1];
		cells[cellid]++;
	}
	int num3 = 0, num1 = 0, num4 = 0, num0 = 0, num2 = 0, numother = 0;
	for(int cellid = 0; cellid < cells.size(); cellid++){
		if(cells[cellid] == 1) num1++;
		else if(cells[cellid] == 0) num0++;
		else if(cells[cellid] == 4) num4++;
		else if(cells[cellid] == 2) num2++;
		else if(cells[cellid] == 3) num3++;
		else numother++;
	}
	printf("0: %d\n1: %d\n2: %d\n3: %d\n4: %d\no: %d\n", num0,num1,num2,num3,num4,numother);
}

int main(int argc, char **argv) {

	int *becell, *ecell, *bound, *bedge, *edge, *cell;
	real *x, *q, *qold, *adt, *res;

	int nnode, ncell, nedge, nbedge, niter;
	real rms;


	// read in grid
	printf("Reading in grid \n");
	init(nnode, ncell, nedge, nbedge, becell, ecell, bound, bedge, edge, cell,
			x, q, qold, adt, res);


	printf("Read ready!\nBegin coloring!\n");
	//create coloring
	std::vector<std::vector<int> > color2edge, color2boundedge;
	coloring(ncell,ecell,nedge,2,color2edge);
	coloring(ncell,becell,nbedge,1,color2boundedge);


	// main time-marching loop
	niter = 1000;

	//create coloring vecs:
	int innerColorNum = color2edge.size();
	int * innerColoringLengths  = malloc(innerColorNum*sizeof(int));//freeeeee
	int * color2edgeArray = malloc(nedge*sizeof(int));//freee 
	int currind = 0;
	for(int i = 0; i < color2edge.size();++i){
		innerColoringLengths[i] = color2edge[i].size();
		for(int j = 0; j < color2edge[i].size();++j){	
			color2edgeArray[currind] = color2edge[i][j];
			currind++;
		}
	}
	
	
	int boundColorNum = color2boundedge.size();
	int * boundColoringLengths = malloc(boundColorNum*sizeof(int));//free
	int * color2bedgeArray = malloc(nbedge*sizeof(int));//free
	currind = 0;
	for(int i = 0; i < color2boundedge.size();++i){
		for(int j = 0; j < color2boundedge[i].size();++j){
			color2bedgeArray[currind] = color2boundedge[i][j];
			currind++;
		}
		boundColoringLengths[i] = color2edge[i].size();
	}
	


	//timer
	//initialise timers for total execution wall time
	#pragma omp target data map(to:ncell,nedge,nbedge,gm1,eps,gam)\
	map(to: cell[:4*ncell],becell[:nbedge],ecell[:2*nedge],edge[:2*nedge],x[:2*nnode],qinf[:4],q[:4*ncell])\
	map(to:innerColorNum,innerColoringLengths[:innerColorNum],color2edgeArray[:nedge])\
	map(to:boundColorNum,boundColoringLengths[:boundColorNum],color2bedgeArray[:nbedge])\
	map(alloc: adt[:ncell], qold[:4*ncell],res[4*ncell],rms)
	{
	for(int iter=1; iter<=niter; iter++) {

		// save old flow solution
		#pragma omp target
		#pragma omp parallel for
		for (int i = 0; i < ncell; i++) {
			for (int n=0; n<4; n++) qold[4*i+n] = q[4*i+n];
	    }

		// predictor/corrector update loop
	    for(int k=0; k<2; k++) {
	      // calculate area/timstep
	      #pragma omp target
		  #pragma omp parallel for
	      for (int i = 0; i < ncell; i++) {
	        real dx,dy, ri,u,v,c;

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

	      // calculate flux residual TODO
	      int offset = 0;
	      for(int color = 0; color < innerColorNum; ++color){
		      offset+= color>0?innerColoringLengths[color-1]:0;
		  #pragma omp target map(to:offset)
	    	  #pragma omp parallel for
		      for(int relind = 0; relind < innerColoringLengths[color];++relind){
	    		int i = color2edgeArray[offset+relind];
	    		real dx,dy,mu, ri, p1,vol1, p2,vol2, f;

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

	      // Apply boundary conditions TODO
	      offset = 0;
	      for(int color = 0; color < boundColorNum; ++color){
		      offset+= color>0?boundColoringLengths[color-1]:0;
		  #pragma omp target map(to:offset)
	    	  #pragma omp parallel for
		      for(int relind = 0; relind < boundColoringLengths[color];++relind){
	    		  int i = color2bedgeArray[offset+relind];
	    		  real dx,dy,mu, ri, p1,vol1, p2,vol2, f;

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
	      }


		#pragma omp target
	      {
	      // update flow field
	      rms = 0.0;
	      #pragma omp parallel for reduction(+:rms)
	      for (int i = 0; i < ncell; i++) {
	        real del, adti;

	        adti = 1.0f/adt[i];

	        for (int n=0; n<4; n++) {
	          del    = adti*res[4*i+n];
	          q[4*i+n]   = qold[4*i+n] - del;
	          res[4*i+n] = 0.0f;
	          rms  += del*del;
	        }
	      }
	      }

	    }

	    // print iteration history

	      if (iter%100 == 0){
#pragma omp target update from(rms)
	    	rms = sqrt(rms/(real) ncell);
	    	printf(" %d  %10.5e \n",iter,rms);
	        //char buf[50];
	        //sprintf(buf,"out%d.vtk",iter);
	        //WriteMeshToVTKAscii(buf, x, nnode, cell, ncell, q);
	      }

	}
	}

	//printf("Max total runtime = \n%f\n", wall_t2 - wall_t1);

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
}

