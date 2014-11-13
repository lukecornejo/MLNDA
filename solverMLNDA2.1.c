#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>  /* clock_t, clock, CLOCKS_PER_SEC */
#include "Eigen/Sparse"
#include "Eigen/Dense"
#include "Eigen/IterativeLinearSolvers"
#include <vector>
#include "ioMLNDA2.1.h"
#include "scMLNDA2.1.h"

using namespace Eigen;
using namespace std;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

extern int maxiter_lo;
extern double epsilon_lo;

extern int eta_star;

extern int nx, ny;
extern double *x, *y, *hx, *hy, *xe, *ye; // grid arrays

extern double k_eff, delta; // K-effective
extern double ****sigmaT, *****sigmaS, ****nuSigmaF, ****chi, ****s_ext; // material arrays
extern double ****phi , ****j_x , ****j_y; // NDA scalar flux and current solution
extern double ***phiB_L, ***phiB_R, ***phiB_B, ***phiB_T; // edge scalar flux from NDA
extern double ****D_xP, ****D_xN, ****D_yP, ****D_yN; // Positive and Negative
extern double ***FL, ***FR, ***FB, ***FT;

// residual data
extern double *res_mbal, *res_ml, *res_mr, *res_mb, *res_mt;
extern int *i_mbal, *j_mbal, *j_ml, *j_mr, *i_mb, *i_mt;
extern int *g_mbal, *g_ml, *g_mr, *g_mb, *g_mt;

// iteration data
extern int num_sol;
extern vector< vector<double> > err_lo, dt_lo, dt_pc;
extern vector< vector<int> > num_logm;

extern ofstream temfile;

clock_t t_pc;


//======================================================================================//
//++ function to solve GNDA problem ++++++++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void GNDAsolutionFS(int preconditioner) {
	int i, j, p, c, N_ukn, N_c, N_xf, N_yf;
	int etaL, g;
	double res=0;
	double **phiG, *phiB_LG, *phiB_RG, *phiB_BG, *phiB_TG, *FLG, *FRG, *FBG, *FTG;
	double **D_xPG, **D_xNG, **D_yPG, **D_yNG, **sigmaTG, **sigmaSG, **nuSigmaFG, **s_extG;
	
	etaL=eta_star-1; g=0;
	phiG=phi[etaL][g];
	phiB_LG=phiB_L[etaL][g]; phiB_RG=phiB_R[etaL][g];
	phiB_BG=phiB_B[etaL][g]; phiB_TG=phiB_T[etaL][g];
	
	FLG=FL[etaL][g]; FRG=FR[etaL][g];
	FBG=FL[etaL][g]; FTG=FR[etaL][g];
	
	D_xPG=D_xP[etaL][g]; D_xNG=D_xN[etaL][g];
	D_yPG=D_yP[etaL][g]; D_yNG=D_yN[etaL][g];
	
	sigmaTG=sigmaT[etaL][g];
	sigmaSG=sigmaS[etaL][g][g];
	nuSigmaFG=nuSigmaF[etaL][g];
	s_extG=s_ext[etaL][g];
	
	temfile<<"Grey NDA Fixed Source Solution\n";
	
	N_c=nx*ny;
	N_xf=2*nx;
	N_yf=2*ny;
	N_ukn=nx*ny+2*nx+2*ny;
	
	VectorXd x(N_ukn);
	VectorXd b(N_ukn);
    SpMat A(N_ukn,N_ukn);
	
	p=0; // initialize solution vector guess to transport solution
	for (j=0; j<ny; j++) {
		for (i=0; i<nx; i++) { x[p]=phiG[i][j]; p++; } // cell centre flux
	}
	for (j=0; j<ny; j++) { x[p]=phiB_LG[j]; p++; } // Left   BC Flux
	for (i=0; i<nx; i++) { x[p]=phiB_BG[i]; p++; } // Bottom BC Flux
	for (j=0; j<ny; j++) { x[p]=phiB_RG[j]; p++; } // Right BC Flux
	for (i=0; i<nx; i++) { x[p]=phiB_TG[i]; p++; } // Top   BC Flux
	
	// Assign matrix A and vector b
	p=0; // ++++++++ Central Cell ++++++++
	for ( j=0; j<ny; j++ ) {
		for ( i=0; i<nx; i++ ) {
			A.insert(p,p)=hy[j]*D_xNG[i+1][j]/xe[i+1] + hy[j]*D_xPG[i][j]/xe[i] + hx[i]*D_yNG[i][j+1]/ye[j+1] + hx[i]*D_yPG[i][j]/ye[j] +
			hx[i]*hy[j]*(sigmaTG[i][j]-sigmaSG[i][j]-nuSigmaFG[i][j]/k_eff);
			b[p]=hx[i]*hy[j]*s_extG[i][j]; p++;
		}
	}
	p=0; // ++++++++ Periphery Cells ++++++++
	if ( nx==1 and ny==1 ) { // Single Cell
		A.insert(0,1)=-hy[0]*D_xNG[0][0]/xe[0];   A.insert(0,2)=-hx[0]*D_yNG[0][0]/ye[0];
		A.insert(0,3)=-hy[0]*D_xPG[1][0]/xe[1];   A.insert(0,4)=-hx[0]*D_yPG[0][1]/ye[1];
	}
	else if ( nx==1 ) { // One Cell Wide
		A.insert(p,N_c)      =-hy[0]*D_xNG[0][0]/xe[0];   A.insert(p,N_c+ny)=-hx[0]*D_yNG[0][0]/ye[0];
		A.insert(p,N_c+ny+nx)=-hy[0]*D_xPG[1][0]/xe[1];   A.insert(p,1)     =-hx[0]*D_yPG[0][1]/ye[1];
		p++; // Bottom Cell
		for ( j=1; j<ny-1; j++ ) { // Middle Cells
			A.insert(p,N_c+j)      =-hy[j]*D_xNG[0][j]/xe[0];   A.insert(p,j-1)=-hx[0]*D_yNG[0][j]  /ye[j];
			A.insert(p,N_c+ny+nx+j)=-hy[j]*D_xPG[1][j]/xe[1];   A.insert(p,j+1)=-hx[0]*D_yPG[0][j+1]/ye[j+1];
			p++;
		}
		A.insert(p,N_c+ny-1)     =-hy[ny-1]*D_xNG[0][ny-1]/xe[0];   A.insert(p,ny-2)           =-hx[0]*D_yNG[0][ny-1]/ye[ny-1];
		A.insert(p,N_c+2*ny+nx-1)=-hy[ny-1]*D_xPG[1][ny-1]/xe[1];   A.insert(p,N_c+2*ny+2*nx-1)=-hx[0]*D_yPG[0][ny]  /ye[ny];
		p++; // Top Cell	
	}
	else if ( ny==1 ) { // One Cell Tall
		A.insert(p,N_c)=-hy[0]*D_xNG[0][0]/xe[0];   A.insert(p,N_c+ny)     =-hx[0]*D_yNG[0][0]/ye[0];
		A.insert(p,1)  =-hy[0]*D_xPG[1][0]/xe[1];   A.insert(p,N_c+2*ny+nx)=-hx[0]*D_yPG[0][1]/ye[1];
		p++; // Left Cell
		for ( i=1; i<nx-1; i++ ) { // Middle Cell
			A.insert(p,i-1)=-hy[0]*D_xNG[i][0]  /xe[i];    A.insert(p,N_c+ny+i)     =-hx[i]*D_yNG[i][0]/ye[0];
			A.insert(p,i+1)=-hy[0]*D_xPG[i+1][0]/xe[i+1];  A.insert(p,N_c+2*ny+nx+i)=-hx[i]*D_yPG[i][1]/ye[1];
			p++;
		}
		A.insert(p,nx-2)       =-hy[0]*D_xNG[nx-1][0]/xe[nx-1];  A.insert(p,N_c+ny+nx-1)    =-hx[nx-1]*D_yNG[nx-1][0]/ye[0];
		A.insert(p,N_c+2*ny+nx-1)=-hy[0]*D_xPG[nx][0]  /xe[nx];    A.insert(p,N_c+2*ny+2*nx-1)=-hx[nx-1]*D_yPG[nx-1][1]/ye[1];
		p++; // Right Cell
	}
	else {
		A.insert(p,N_c)=-hy[0]*D_xNG[0][0]/xe[0];   A.insert(p,N_c+ny)=-hx[0]*D_yNG[0][0]/ye[0];
		A.insert(p,1)  =-hy[0]*D_xPG[1][0]/xe[1];   A.insert(p,nx)    =-hx[0]*D_yPG[0][1]/ye[1];
		p++; // Bottom Left Corner
		for ( i=1; i<nx-1; i++ ) { // Bottom
			A.insert(p,p-1)=-hy[0]*D_xNG[i][0]  /xe[i];    A.insert(p,N_c+ny+i)=-hx[i]*D_yNG[i][0]/ye[0];
			A.insert(p,p+1)=-hy[0]*D_xPG[i+1][0]/xe[i+1];  A.insert(p,nx+i)    =-hx[i]*D_yPG[i][1]/ye[1];
			p++;
		}
		A.insert(p,nx-2)     =-hy[0]*D_xNG[nx-1][0]/xe[nx-1]; A.insert(p,N_c+ny+nx-1)=-hx[nx-1]*D_yNG[nx-1][0]/ye[0];
		A.insert(p,N_c+ny+nx)=-hy[0]*D_xPG[nx][0]  /xe[nx];   A.insert(p,2*nx-1)     =-hx[nx-1]*D_yPG[nx-1][1]/ye[1];
		p++; // Bottom Right Corner
		for ( j=1; j<ny-1; j++ ) { // Middle
			i=0; // Left Side
			A.insert(p,N_c+j)=-hy[j]*D_xNG[i][j]  /xe[i];    A.insert(p,p-nx)=-hx[i]*D_yNG[i][j]  /ye[j];
			A.insert(p,p+1)  =-hy[j]*D_xPG[i+1][j]/xe[i+1];  A.insert(p,p+nx)=-hx[i]*D_yPG[i][j+1]/ye[j+1];
			p++;
			for ( i=1; i<nx-1; i++ ) { // Centre Cells
				A.insert(p,p-1)=-hy[j]*D_xNG[i][j]  /xe[i];    A.insert(p,p-nx)=-hx[i]*D_yNG[i][j]  /ye[j];
				A.insert(p,p+1)=-hy[j]*D_xPG[i+1][j]/xe[i+1];  A.insert(p,p+nx)=-hx[i]*D_yPG[i][j+1]/ye[j+1];
				p++;
			}
			i=nx-1; // Right Side
			A.insert(p,p-1)        =-hy[j]*D_xNG[i][j]  /xe[i];    A.insert(p,p-nx)=-hx[i]*D_yNG[i][j]  /ye[j];
			A.insert(p,N_c+ny+nx+j)=-hy[j]*D_xPG[i+1][j]/xe[i+1];  A.insert(p,p+nx)=-hx[i]*D_yPG[i][j+1]/ye[j+1];
			p++;
		}
		j=ny-1; // Top Left Corner
		A.insert(p,N_c+j)=-hy[j]*D_xNG[0][j]/xe[0];   A.insert(p,p-nx)       =-hx[0]*D_yNG[0][j]  /ye[j];
		A.insert(p,p+1)  =-hy[j]*D_xPG[1][j]/xe[1];   A.insert(p,N_c+2*ny+nx)=-hx[0]*D_yPG[0][j+1]/ye[j+1];
		p++;
		for ( i=1; i<nx-1; i++ ) { // Top
			A.insert(p,p-1)=-hy[j]*D_xNG[i][j]  /xe[i];    A.insert(p,p-nx)         =-hx[i]*D_yNG[i][j]/ye[j];
			A.insert(p,p+1)=-hy[j]*D_xPG[i+1][j]/xe[i+1];  A.insert(p,N_c+2*ny+nx+i)=-hx[i]*D_yPG[i][j+1]/ye[j+1];
			p++;
		}
		i=nx-1; // Top Right Corner
		A.insert(p,p-1)          =-hy[j]*D_xNG[i][j]  /xe[i];    A.insert(p,p-nx)         =-hx[i]*D_yNG[i][j]/ye[j];
		A.insert(p,N_c+2*ny+nx-1)=-hy[j]*D_xPG[i+1][j]/xe[i+1];  A.insert(p,N_c+2*ny+nx+i)=-hx[i]*D_yPG[i][j+1]/ye[j+1];
		p++;
	}
	// ++++++++ Boundary Conditions ++++++++
	p=N_c;
	for ( j=0; j<ny; j++ ) {
		A.insert(p,j*nx)= D_xPG[0][j]/xe[0];   A.insert(p,N_c+j)=FLG[j]-D_xNG[0][j]/xe[0];   b[p]=0; p++; // Left   BC
	}
	for ( i=0; i<nx; i++ ) {
		A.insert(p,i)= D_yPG[i][0]/ye[0];   A.insert(p,N_c+ny+i)=FBG[i]-D_yNG[i][0]/ye[0];   b[p]=0; p++; // Bottom BC
	}
	for ( j=0; j<ny; j++ ) {
		A.insert(p,(j+1)*nx-1)=-D_xNG[nx][j]/xe[nx];   A.insert(p,N_c+ny+nx+j)=FRG[j]+D_xPG[nx][j]/xe[nx];   b[p]=0; p++; // Right  BC
	}
	for ( i=0; i<nx; i++ ) {
		A.insert(p,(ny-1)*nx+i)=-D_yNG[i][ny]/ye[ny];   A.insert(p,N_c+2*ny+nx+i)=FTG[i]+D_yPG[i][ny]/ye[ny];   b[p]=0; p++; // Top    BC
	}
	
	A.prune(1e-17, 10);
	
	// Solve Ax=b iteratively with BiCGSTAB
	if ( preconditioner==2 ) {
		t_pc=clock(); // start low-order timer
		BiCGSTAB<SpMat> solver;
		solver.setTolerance(epsilon_lo);     // set convergence criteria 
		solver.setMaxIterations(maxiter_lo); // set the max number of lo iterations
		solver.compute(A);
		
		t_pc=clock()-t_pc; // stop low-order timer
		dt_pc[etaL].push_back(((double)t_pc)/CLOCKS_PER_SEC); // add high-order solution time to vector
		
		x = solver.solveWithGuess(b,x);
		err_lo[etaL].push_back(solver.error());      // error in lo solution
		num_logm[etaL].push_back(solver.iterations()); // number of lo iterations
	}
	else if ( preconditioner==3 ) {
		t_pc=clock(); // start low-order timer
		BiCGSTAB<SpMat,IncompleteLUT<double> > solver;
		solver.preconditioner().setFillfactor(11);
		solver.setTolerance(epsilon_lo);     // set convergence criteria 
		solver.setMaxIterations(maxiter_lo); // set the max number of lo iterations
		solver.compute(A);
		
		t_pc=clock()-t_pc; // stop low-order timer
		dt_pc[etaL].push_back(((double)t_pc)/CLOCKS_PER_SEC); // add high-order solution time to vector
		
		x = solver.solveWithGuess(b,x);
		err_lo[etaL].push_back(solver.error());      // error in lo solution
		num_logm[etaL].push_back(solver.iterations()); // number of lo iterations
	}
	else cout<<">>Preconditioner Type Error \n";
	
	
	
	num_sol++;
	
	// set solution back to problem values
	res_mbal[etaL]=0.0;
	res_ml[etaL]=0.0;
	res_mr[etaL]=0.0;
	res_mb[etaL]=0.0;
	res_mt[etaL]=0.0;
	p=0;
	for (j=0; j<ny; j++) {
		for (i=0; i<nx; i++) { phiG[i][j]=x[p]; p++; } // Cell Centre Flux
	}
	for (j=0; j<ny; j++) { phiB_LG[j]=x[p]; p++; } // Left   BC Flux
	for (i=0; i<nx; i++) { phiB_BG[i]=x[p]; p++; } // Bottom BC Flux
	for (j=0; j<ny; j++) { phiB_RG[j]=x[p]; p++; } // Right BC Flux
	for (i=0; i<nx; i++) { phiB_TG[i]=x[p]; p++; } // Top   BC Flux
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Calculate current ///////////////////////////////////////////////////////////////////////////////////////
	// Boundary currents ///////////////////////////////////////////////////////////////////////////////////////
	for (j=0; j<ny; j++) j_x[etaL][g][0][j]=FLG[j]*phiB_LG[j]; // Left boundary current J_x
	for (i=0; i<nx; i++) j_y[etaL][g][i][0]=FBG[i]*phiB_BG[i]; // Bottom boundary current J_y
	// Inner currents J_x J_y
	for (j=0; j<ny; j++) for (i=1; i<nx; i++) j_x[etaL][g][i][j]=-(D_xPG[i][j]*phiG[i][j]-D_xNG[i][j]*phiG[i-1][j])/xe[i];
	for (j=1; j<ny; j++) for (i=0; i<nx; i++) j_y[etaL][g][i][j]=-(D_yPG[i][j]*phiG[i][j]-D_yNG[i][j]*phiG[i][j-1])/ye[j];
	// Boundary currents
	for (j=0; j<ny; j++) j_x[etaL][g][nx][j]=FRG[j]*phiB_RG[j]; // Right boundary current
	for (i=0; i<nx; i++) j_y[etaL][g][i][ny]=FTG[i]*phiB_TG[i]; // Top boundary current
	
	// Calculate residuals of matrix equations
	i=0; j=0;
	res=abs(-hy[j]*D_xPG[i+1][j]*phiG[i+1][j]/xe[i+1]-hx[i]*D_yPG[i][j+1]*phiG[i][j+1]/ye[j+1]
	+(hy[j]*D_xNG[i+1][j]/xe[i+1]+hy[j]*D_xPG[i][j]/xe[i]+hx[i]*D_yNG[i][j+1]/ye[j+1]+hx[i]*D_yPG[i][j]/ye[j]
	+hx[i]*hy[j]*(sigmaTG[i][j]-sigmaSG[i][j]-nuSigmaFG[i][j]/k_eff))*phiG[i][j]
	-hy[j]*D_xNG[i][j]*phiB_LG[j]/xe[i]-hx[i]*D_yNG[i][j]*phiB_BG[i]/ye[j]-hx[i]*hy[j]*s_extG[i][j]);
	if (res>res_mbal[etaL]) { res_mbal[etaL]=res; i_mbal[etaL]=i; j_mbal[etaL]=j; }
	for (j=1; j<ny-1; j++) {
		i=0;
		res=abs(-hy[j]*D_xPG[i+1][j]*phiG[i+1][j]/xe[i+1]-hx[i]*D_yPG[i][j+1]*phiG[i][j+1]/ye[j+1]
		+(hy[j]*D_xNG[i+1][j]/xe[i+1]+hy[j]*D_xPG[i][j]/xe[i]+hx[i]*D_yNG[i][j+1]/ye[j+1]+hx[i]*D_yPG[i][j]/ye[j]
		+hx[i]*hy[j]*(sigmaTG[i][j]-sigmaSG[i][j]-nuSigmaFG[i][j]/k_eff))*phiG[i][j]
		-hy[j]*D_xNG[i][j]*phiB_LG[j]/xe[i]-hx[i]*D_yNG[i][j]*phiG[i][j-1]/ye[j]-hx[i]*hy[j]*s_extG[i][j]);
		if (res>res_mbal[etaL]) { res_mbal[etaL]=res; i_mbal[etaL]=i; j_mbal[etaL]=j; }
		for (i=1; i<nx-1; i++) {
			res=abs(-hy[j]*D_xPG[i+1][j]*phiG[i+1][j]/xe[i+1]-hx[i]*D_yPG[i][j+1]*phiG[i][j+1]/ye[j+1]
			+(hy[j]*D_xNG[i+1][j]/xe[i+1]+hy[j]*D_xPG[i][j]/xe[i]+hx[i]*D_yNG[i][j+1]/ye[j+1]+hx[i]*D_yPG[i][j]/ye[j]
			+hx[i]*hy[j]*(sigmaTG[i][j]-sigmaSG[i][j]-nuSigmaFG[i][j]/k_eff))*phiG[i][j]
			-hy[j]*D_xNG[i][j]*phiG[i-1][j]/xe[i]-hx[i]*D_yNG[i][j]*phiG[i][j-1]/ye[j]-hx[i]*hy[j]*s_extG[i][j]);
			if (res>res_mbal[etaL]) { res_mbal[etaL]=res; i_mbal[etaL]=i; j_mbal[etaL]=j; }
		}
		i=nx-1;
		res=abs(-hy[j]*D_xPG[i+1][j]*phiB_RG[j]/xe[i+1]-hx[i]*D_yPG[i][j+1]*phiG[i][j+1]/ye[j+1]
		+(hy[j]*D_xNG[i+1][j]/xe[i+1]+hy[j]*D_xPG[i][j]/xe[i]+hx[i]*D_yNG[i][j+1]/ye[j+1]+hx[i]*D_yPG[i][j]/ye[j]
		+hx[i]*hy[j]*(sigmaTG[i][j]-sigmaSG[i][j]-nuSigmaFG[i][j]/k_eff))*phiG[i][j]
		-hy[j]*D_xNG[i][j]*phiG[i-1][j]/xe[i]-hx[i]*D_yNG[i][j]*phiG[i][j-1]/ye[j]-hx[i]*hy[j]*s_extG[i][j]);
		if (res>res_mbal[etaL]) { res_mbal[etaL]=res; i_mbal[etaL]=i; j_mbal[etaL]=j; }
	}
	i=nx-1; j=ny-1;
	res=abs(-hy[j]*D_xPG[i+1][j]*phiB_RG[j]/xe[i+1]-hx[i]*D_yPG[i][j+1]*phiB_TG[i]/ye[j+1]
	+(hy[j]*D_xNG[i+1][j]/xe[i+1]+hy[j]*D_xPG[i][j]/xe[i]+hx[i]*D_yNG[i][j+1]/ye[j+1]+hx[i]*D_yPG[i][j]/ye[j]
	+hx[i]*hy[j]*(sigmaTG[i][j]-sigmaSG[i][j]-nuSigmaFG[i][j]/k_eff))*phiG[i][j]
	-hy[j]*D_xNG[i][j]*phiG[i-1][j]/xe[i]-hx[i]*D_yNG[i][j]*phiG[i][j-1]/ye[j]-hx[i]*hy[j]*s_extG[i][j]);
	if (res>res_mbal[etaL]) { res_mbal[etaL]=res; g_mbal[etaL]=g; i_mbal[etaL]=i; j_mbal[etaL]=j; }
	
	for (j=1; j<ny-1; j++) {
		res=abs( D_xPG[0][j] *phiG[0][j]   /xe[0] +(FLG[j]-D_xNG[0][j] /xe[0])*phiB_LG[j]);
		if (res>res_ml[etaL]) { res_ml[etaL]=res; g_ml[etaL]=g; j_ml[etaL]=j; }
		res=abs(-D_xNG[nx][j]*phiG[nx-1][j]/xe[nx]+(FRG[j]+D_xPG[nx][j]/xe[nx])*phiB_RG[j]);
		if (res>res_mr[etaL]) { res_mr[etaL]=res; g_ml[etaL]=g; j_mr[etaL]=j; }
	}
	for (i=1; i<nx-1; i++) {
		res=abs( D_yPG[i][0] *phiG[i][0]   /ye[0] +(FBG[i]-D_yNG[i][0]/ye[0] )*phiB_BG[i]);
		if (res>res_mb[etaL]) { res_mb[etaL]=res; g_ml[etaL]=g; i_mb[etaL]=j; }
		res=abs(-D_yNG[i][ny]*phiG[i][ny-1]/ye[ny]+(FTG[i]+D_yPG[i][ny]/ye[ny])*phiB_TG[i]);
		if (res>res_mt[etaL]) { res_mt[etaL]=res; g_ml[etaL]=g; i_mt[etaL]=j; }
	}
	
}
//======================================================================================//

//======================================================================================//
//++ function to solve GNDA problem ++++++++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void GNDAsolutionKE(int preconditioner, double k_effL, vector< vector<double> > phi_gL)
{
	int i, j, p, c, N_ukn, N_c, N_xf, N_yf;
	int etaL, g;
	double res=0;
	double **phiG, *phiB_LG, *phiB_RG, *phiB_BG, *phiB_TG, *FLG, *FRG, *FBG, *FTG;
	double **D_xPG, **D_xNG, **D_yPG, **D_yNG, **sigmaTG, **sigmaSG, **nuSigmaFG;
	
	etaL=eta_star-1; g=0;
	phiG=phi[etaL][g];
	phiB_LG=phiB_L[etaL][g]; phiB_RG=phiB_R[etaL][g];
	phiB_BG=phiB_B[etaL][g]; phiB_TG=phiB_T[etaL][g];
	
	FLG=FL[etaL][g]; FRG=FR[etaL][g];
	FBG=FL[etaL][g]; FTG=FR[etaL][g];
	
	D_xPG=D_xP[etaL][g]; D_xNG=D_xN[etaL][g];
	D_yPG=D_yP[etaL][g]; D_yNG=D_yN[etaL][g];
	
	sigmaTG=sigmaT[etaL][g];
	sigmaSG=sigmaS[etaL][g][g];
	nuSigmaFG=nuSigmaF[etaL][g];
	
	temfile<<"Grey NDA K Eigenvalue Solution\n";
	
	N_c=nx*ny;
	N_xf=2*nx;
	N_yf=2*ny;
	N_ukn=nx*ny+2*nx+2*ny;
	
	VectorXd x(N_ukn);
	VectorXd b(N_ukn);
    SpMat A(N_ukn,N_ukn);
	
	p=0; // initialize solution vector guess to transport solution
	for (j=0; j<ny; j++) {
		for (i=0; i<nx; i++) { x[p]=phiG[i][j]; p++; } // cell centre flux
	}
	for (j=0; j<ny; j++) { x[p]=phiB_LG[j]; p++; } // Left   BC Flux
	for (i=0; i<nx; i++) { x[p]=phiB_BG[i]; p++; } // Bottom BC Flux
	for (j=0; j<ny; j++) { x[p]=phiB_RG[j]; p++; } // Right BC Flux
	for (i=0; i<nx; i++) { x[p]=phiB_TG[i]; p++; } // Top   BC Flux
	
	// Assign matrix A and vector b
	p=0; // ++++++++ Central Cell ++++++++
	for ( j=0; j<ny; j++ ) {
		for ( i=0; i<nx; i++ ) {
			A.insert(p,p)=hy[j]*D_xNG[i+1][j]/xe[i+1] + hy[j]*D_xPG[i][j]/xe[i] + hx[i]*D_yNG[i][j+1]/ye[j+1] + hx[i]*D_yPG[i][j]/ye[j] +
			hx[i]*hy[j]*(sigmaTG[i][j]-sigmaSG[i][j]-(1-delta)*nuSigmaFG[i][j]/k_effL);
			b[p]=hx[i]*hy[j]*(delta*nuSigmaFG[i][j]*phi_gL[i][j]/k_effL); p++;
		}
	}
	p=0; // ++++++++ Periphery Cells ++++++++
	if ( nx==1 and ny==1 ) { // Single Cell
		A.insert(0,1)=-hy[0]*D_xNG[0][0]/xe[0];   A.insert(0,2)=-hx[0]*D_yNG[0][0]/ye[0];
		A.insert(0,3)=-hy[0]*D_xPG[1][0]/xe[1];   A.insert(0,4)=-hx[0]*D_yPG[0][1]/ye[1];
	}
	else if ( nx==1 ) { // One Cell Wide
		A.insert(p,N_c)      =-hy[0]*D_xNG[0][0]/xe[0];   A.insert(p,N_c+ny)=-hx[0]*D_yNG[0][0]/ye[0];
		A.insert(p,N_c+ny+nx)=-hy[0]*D_xPG[1][0]/xe[1];   A.insert(p,1)     =-hx[0]*D_yPG[0][1]/ye[1];
		p++; // Bottom Cell
		for ( j=1; j<ny-1; j++ ) { // Middle Cells
			A.insert(p,N_c+j)      =-hy[j]*D_xNG[0][j]/xe[0];   A.insert(p,j-1)=-hx[0]*D_yNG[0][j]  /ye[j];
			A.insert(p,N_c+ny+nx+j)=-hy[j]*D_xPG[1][j]/xe[1];   A.insert(p,j+1)=-hx[0]*D_yPG[0][j+1]/ye[j+1];
			p++;
		}
		A.insert(p,N_c+ny-1)     =-hy[ny-1]*D_xNG[0][ny-1]/xe[0];   A.insert(p,ny-2)           =-hx[0]*D_yNG[0][ny-1]/ye[ny-1];
		A.insert(p,N_c+2*ny+nx-1)=-hy[ny-1]*D_xPG[1][ny-1]/xe[1];   A.insert(p,N_c+2*ny+2*nx-1)=-hx[0]*D_yPG[0][ny]  /ye[ny];
		p++; // Top Cell	
	}
	else if ( ny==1 ) { // One Cell Tall
		A.insert(p,N_c)=-hy[0]*D_xNG[0][0]/xe[0];   A.insert(p,N_c+ny)     =-hx[0]*D_yNG[0][0]/ye[0];
		A.insert(p,1)  =-hy[0]*D_xPG[1][0]/xe[1];   A.insert(p,N_c+2*ny+nx)=-hx[0]*D_yPG[0][1]/ye[1];
		p++; // Left Cell
		for ( i=1; i<nx-1; i++ ) { // Middle Cell
			A.insert(p,i-1)=-hy[0]*D_xNG[i][0]  /xe[i];    A.insert(p,N_c+ny+i)     =-hx[i]*D_yNG[i][0]/ye[0];
			A.insert(p,i+1)=-hy[0]*D_xPG[i+1][0]/xe[i+1];  A.insert(p,N_c+2*ny+nx+i)=-hx[i]*D_yPG[i][1]/ye[1];
			p++;
		}
		A.insert(p,nx-2)       =-hy[0]*D_xNG[nx-1][0]/xe[nx-1];  A.insert(p,N_c+ny+nx-1)    =-hx[nx-1]*D_yNG[nx-1][0]/ye[0];
		A.insert(p,N_c+2*ny+nx-1)=-hy[0]*D_xPG[nx][0]  /xe[nx];    A.insert(p,N_c+2*ny+2*nx-1)=-hx[nx-1]*D_yPG[nx-1][1]/ye[1];
		p++; // Right Cell
	}
	else {
		A.insert(p,N_c)=-hy[0]*D_xNG[0][0]/xe[0];   A.insert(p,N_c+ny)=-hx[0]*D_yNG[0][0]/ye[0];
		A.insert(p,1)  =-hy[0]*D_xPG[1][0]/xe[1];   A.insert(p,nx)    =-hx[0]*D_yPG[0][1]/ye[1];
		p++; // Bottom Left Corner
		for ( i=1; i<nx-1; i++ ) { // Bottom
			A.insert(p,p-1)=-hy[0]*D_xNG[i][0]  /xe[i];    A.insert(p,N_c+ny+i)=-hx[i]*D_yNG[i][0]/ye[0];
			A.insert(p,p+1)=-hy[0]*D_xPG[i+1][0]/xe[i+1];  A.insert(p,nx+i)    =-hx[i]*D_yPG[i][1]/ye[1];
			p++;
		}
		A.insert(p,nx-2)     =-hy[0]*D_xNG[nx-1][0]/xe[nx-1]; A.insert(p,N_c+ny+nx-1)=-hx[nx-1]*D_yNG[nx-1][0]/ye[0];
		A.insert(p,N_c+ny+nx)=-hy[0]*D_xPG[nx][0]  /xe[nx];   A.insert(p,2*nx-1)     =-hx[nx-1]*D_yPG[nx-1][1]/ye[1];
		p++; // Bottom Right Corner
		for ( j=1; j<ny-1; j++ ) { // Middle
			i=0; // Left Side
			A.insert(p,N_c+j)=-hy[j]*D_xNG[i][j]  /xe[i];    A.insert(p,p-nx)=-hx[i]*D_yNG[i][j]  /ye[j];
			A.insert(p,p+1)  =-hy[j]*D_xPG[i+1][j]/xe[i+1];  A.insert(p,p+nx)=-hx[i]*D_yPG[i][j+1]/ye[j+1];
			p++;
			for ( i=1; i<nx-1; i++ ) { // Centre Cells
				A.insert(p,p-1)=-hy[j]*D_xNG[i][j]  /xe[i];    A.insert(p,p-nx)=-hx[i]*D_yNG[i][j]  /ye[j];
				A.insert(p,p+1)=-hy[j]*D_xPG[i+1][j]/xe[i+1];  A.insert(p,p+nx)=-hx[i]*D_yPG[i][j+1]/ye[j+1];
				p++;
			}
			i=nx-1; // Right Side
			A.insert(p,p-1)        =-hy[j]*D_xNG[i][j]  /xe[i];    A.insert(p,p-nx)=-hx[i]*D_yNG[i][j]  /ye[j];
			A.insert(p,N_c+ny+nx+j)=-hy[j]*D_xPG[i+1][j]/xe[i+1];  A.insert(p,p+nx)=-hx[i]*D_yPG[i][j+1]/ye[j+1];
			p++;
		}
		j=ny-1; // Top Left Corner
		A.insert(p,N_c+j)=-hy[j]*D_xNG[0][j]/xe[0];   A.insert(p,p-nx)       =-hx[0]*D_yNG[0][j]  /ye[j];
		A.insert(p,p+1)  =-hy[j]*D_xPG[1][j]/xe[1];   A.insert(p,N_c+2*ny+nx)=-hx[0]*D_yPG[0][j+1]/ye[j+1];
		p++;
		for ( i=1; i<nx-1; i++ ) { // Top
			A.insert(p,p-1)=-hy[j]*D_xNG[i][j]  /xe[i];    A.insert(p,p-nx)         =-hx[i]*D_yNG[i][j]/ye[j];
			A.insert(p,p+1)=-hy[j]*D_xPG[i+1][j]/xe[i+1];  A.insert(p,N_c+2*ny+nx+i)=-hx[i]*D_yPG[i][j+1]/ye[j+1];
			p++;
		}
		i=nx-1; // Top Right Corner
		A.insert(p,p-1)          =-hy[j]*D_xNG[i][j]  /xe[i];    A.insert(p,p-nx)         =-hx[i]*D_yNG[i][j]/ye[j];
		A.insert(p,N_c+2*ny+nx-1)=-hy[j]*D_xPG[i+1][j]/xe[i+1];  A.insert(p,N_c+2*ny+nx+i)=-hx[i]*D_yPG[i][j+1]/ye[j+1];
		p++;
	}
	// ++++++++ Boundary Conditions ++++++++
	p=N_c;
	for ( j=0; j<ny; j++ ) {
		A.insert(p,j*nx)= D_xPG[0][j]/xe[0];   A.insert(p,N_c+j)=FLG[j]-D_xNG[0][j]/xe[0];   b[p]=0; p++; // Left   BC
	}
	for ( i=0; i<nx; i++ ) {
		A.insert(p,i)= D_yPG[i][0]/ye[0];   A.insert(p,N_c+ny+i)=FBG[i]-D_yNG[i][0]/ye[0];   b[p]=0; p++; // Bottom BC
	}
	for ( j=0; j<ny; j++ ) {
		A.insert(p,(j+1)*nx-1)=-D_xNG[nx][j]/xe[nx];   A.insert(p,N_c+ny+nx+j)=FRG[j]+D_xPG[nx][j]/xe[nx];   b[p]=0; p++; // Right  BC
	}
	for ( i=0; i<nx; i++ ) {
		A.insert(p,(ny-1)*nx+i)=-D_yNG[i][ny]/ye[ny];   A.insert(p,N_c+2*ny+nx+i)=FTG[i]+D_yPG[i][ny]/ye[ny];   b[p]=0; p++; // Top    BC
	}
	
	A.prune(1e-17, 10);
	
	// Solve Ax=b iteratively with BiCGSTAB
	if ( preconditioner==2 ) {
		t_pc=clock(); // start low-order timer
		BiCGSTAB<SpMat> solver;
		solver.setTolerance(epsilon_lo);     // set convergence criteria 
		solver.setMaxIterations(maxiter_lo); // set the max number of lo iterations
		solver.compute(A);
		
		t_pc=clock()-t_pc; // stop low-order timer
		dt_pc[etaL].push_back(((double)t_pc)/CLOCKS_PER_SEC); // add high-order solution time to vector
		
		x = solver.solveWithGuess(b,x);
		err_lo[etaL].push_back(solver.error());      // error in lo solution
		num_logm[etaL].push_back(solver.iterations()); // number of lo iterations
	}
	else if ( preconditioner==3 ) {
		t_pc=clock(); // start low-order timer
		BiCGSTAB<SpMat,IncompleteLUT<double> > solver;
		solver.preconditioner().setFillfactor(11);
		solver.setTolerance(epsilon_lo);     // set convergence criteria 
		solver.setMaxIterations(maxiter_lo); // set the max number of lo iterations
		solver.compute(A);
		
		t_pc=clock()-t_pc; // stop low-order timer
		dt_pc[etaL].push_back(((double)t_pc)/CLOCKS_PER_SEC); // add high-order solution time to vector
		
		x = solver.solveWithGuess(b,x);
		err_lo[etaL].push_back(solver.error());      // error in lo solution
		num_logm[etaL].push_back(solver.iterations()); // number of lo iterations
	}
	else cout<<">>Preconditioner Type Error \n";
	num_sol++;
	
	//cout<<"x "<<x<<endl;
	
	// set solution back to problem values
	res_mbal[etaL]=0.0;
	res_ml[etaL]=0.0;
	res_mr[etaL]=0.0;
	res_mb[etaL]=0.0;
	res_mt[etaL]=0.0;
	p=0;
	for (j=0; j<ny; j++) {
		for (i=0; i<nx; i++) { phiG[i][j]=x[p]; p++; } // Cell Centre Flux
	}
	for (j=0; j<ny; j++) { phiB_LG[j]=x[p]; p++; } // Left   BC Flux
	for (i=0; i<nx; i++) { phiB_BG[i]=x[p]; p++; } // Bottom BC Flux
	for (j=0; j<ny; j++) { phiB_RG[j]=x[p]; p++; } // Right BC Flux
	for (i=0; i<nx; i++) { phiB_TG[i]=x[p]; p++; } // Top   BC Flux
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Calculate current ///////////////////////////////////////////////////////////////////////////////////////
	// Boundary currents ///////////////////////////////////////////////////////////////////////////////////////
	for (j=0; j<ny; j++) j_x[etaL][g][0][j]=FLG[j]*phiB_LG[j]; // Left boundary current J_x
	for (i=0; i<nx; i++) j_y[etaL][g][i][0]=FBG[i]*phiB_BG[i]; // Bottom boundary current J_y
	// Inner currents J_x J_y
	for (j=0; j<ny; j++) for (i=1; i<nx; i++) j_x[etaL][g][i][j]=-(D_xPG[i][j]*phiG[i][j]-D_xNG[i][j]*phiG[i-1][j])/xe[i];
	for (j=1; j<ny; j++) for (i=0; i<nx; i++) j_y[etaL][g][i][j]=-(D_yPG[i][j]*phiG[i][j]-D_yNG[i][j]*phiG[i][j-1])/ye[j];
	// Boundary currents
	for (j=0; j<ny; j++) j_x[etaL][g][nx][j]=FRG[j]*phiB_RG[j]; // Right boundary current
	for (i=0; i<nx; i++) j_y[etaL][g][i][ny]=FTG[i]*phiB_TG[i]; // Top boundary current
	
	// Calculate residuals of matrix equations
	i=0; j=0;
	res=abs(-hy[j]*D_xPG[i+1][j]*phiG[i+1][j]/xe[i+1]-hx[i]*D_yPG[i][j+1]*phiG[i][j+1]/ye[j+1]
	+(hy[j]*D_xNG[i+1][j]/xe[i+1]+hy[j]*D_xPG[i][j]/xe[i]+hx[i]*D_yNG[i][j+1]/ye[j+1]+hx[i]*D_yPG[i][j]/ye[j]
	+hx[i]*hy[j]*(sigmaTG[i][j]-sigmaSG[i][j]-(1-delta)*nuSigmaFG[i][j]/k_effL))*phiG[i][j]
	-hy[j]*D_xNG[i][j]*phiB_LG[j]/xe[i]-hx[i]*D_yNG[i][j]*phiB_BG[i]/ye[j]
	-hx[i]*hy[j]*(delta*nuSigmaFG[i][j]*phi_gL[i][j]/k_effL));
	if (res>res_mbal[etaL]) { res_mbal[etaL]=res; i_mbal[etaL]=i; j_mbal[etaL]=j; }
	for (j=1; j<ny-1; j++) {
		i=0;
		res=abs(-hy[j]*D_xPG[i+1][j]*phiG[i+1][j]/xe[i+1]-hx[i]*D_yPG[i][j+1]*phiG[i][j+1]/ye[j+1]
		+(hy[j]*D_xNG[i+1][j]/xe[i+1]+hy[j]*D_xPG[i][j]/xe[i]+hx[i]*D_yNG[i][j+1]/ye[j+1]+hx[i]*D_yPG[i][j]/ye[j]
		+hx[i]*hy[j]*(sigmaTG[i][j]-sigmaSG[i][j]-(1-delta)*nuSigmaFG[i][j]/k_effL))*phiG[i][j]
		-hy[j]*D_xNG[i][j]*phiB_LG[j]/xe[i]-hx[i]*D_yNG[i][j]*phiG[i][j-1]/ye[j]
		-hx[i]*hy[j]*(delta*nuSigmaFG[i][j]*phi_gL[i][j]/k_effL));
		if (res>res_mbal[etaL]) { res_mbal[etaL]=res; i_mbal[etaL]=i; j_mbal[etaL]=j; }
		for (i=1; i<nx-1; i++) {
			res=abs(-hy[j]*D_xPG[i+1][j]*phiG[i+1][j]/xe[i+1]-hx[i]*D_yPG[i][j+1]*phiG[i][j+1]/ye[j+1]
			+(hy[j]*D_xNG[i+1][j]/xe[i+1]+hy[j]*D_xPG[i][j]/xe[i]+hx[i]*D_yNG[i][j+1]/ye[j+1]+hx[i]*D_yPG[i][j]/ye[j]
			+hx[i]*hy[j]*(sigmaTG[i][j]-sigmaSG[i][j]-(1-delta)*nuSigmaFG[i][j]/k_effL))*phiG[i][j]
			-hy[j]*D_xNG[i][j]*phiG[i-1][j]/xe[i]-hx[i]*D_yNG[i][j]*phiG[i][j-1]/ye[j]
			-hx[i]*hy[j]*(delta*nuSigmaFG[i][j]*phi_gL[i][j]/k_effL));
			if (res>res_mbal[etaL]) { res_mbal[etaL]=res; i_mbal[etaL]=i; j_mbal[etaL]=j; }
		}
		i=nx-1;
		res=abs(-hy[j]*D_xPG[i+1][j]*phiB_RG[j]/xe[i+1]-hx[i]*D_yPG[i][j+1]*phiG[i][j+1]/ye[j+1]
		+(hy[j]*D_xNG[i+1][j]/xe[i+1]+hy[j]*D_xPG[i][j]/xe[i]+hx[i]*D_yNG[i][j+1]/ye[j+1]+hx[i]*D_yPG[i][j]/ye[j]
		+hx[i]*hy[j]*(sigmaTG[i][j]-sigmaSG[i][j]-(1-delta)*nuSigmaFG[i][j]/k_effL))*phiG[i][j]
		-hy[j]*D_xNG[i][j]*phiG[i-1][j]/xe[i]-hx[i]*D_yNG[i][j]*phiG[i][j-1]/ye[j]
		-hx[i]*hy[j]*delta*nuSigmaFG[i][j]*phi_gL[i][j]/k_effL);
		if (res>res_mbal[etaL]) { res_mbal[etaL]=res; i_mbal[etaL]=i; j_mbal[etaL]=j; }
	}
	i=nx-1; j=ny-1;
	res=abs(-hy[j]*D_xPG[i+1][j]*phiB_RG[j]/xe[i+1]-hx[i]*D_yPG[i][j+1]*phiB_TG[i]/ye[j+1]
	+(hy[j]*D_xNG[i+1][j]/xe[i+1]+hy[j]*D_xPG[i][j]/xe[i]+hx[i]*D_yNG[i][j+1]/ye[j+1]+hx[i]*D_yPG[i][j]/ye[j]
	+hx[i]*hy[j]*(sigmaTG[i][j]-sigmaSG[i][j]-(1-delta)*nuSigmaFG[i][j]/k_effL))*phiG[i][j]
	-hy[j]*D_xNG[i][j]*phiG[i-1][j]/xe[i]-hx[i]*D_yNG[i][j]*phiG[i][j-1]/ye[j]
	-hx[i]*hy[j]*delta*nuSigmaFG[i][j]*phi_gL[i][j]/k_effL);
	if (res>res_mbal[etaL]) { res_mbal[etaL]=res; g_mbal[etaL]=g; i_mbal[etaL]=i; j_mbal[etaL]=j; }
	
	for (j=1; j<ny-1; j++) {
		res=abs( D_xPG[0][j] *phiG[0][j]   /xe[0] +(FLG[j]-D_xNG[0][j] /xe[0])*phiB_LG[j]);
		if (res>res_ml[etaL]) { res_ml[etaL]=res; g_ml[etaL]=g; j_ml[etaL]=j; }
		res=abs(-D_xNG[nx][j]*phiG[nx-1][j]/xe[nx]+(FRG[j]+D_xPG[nx][j]/xe[nx])*phiB_RG[j]);
		if (res>res_mr[etaL]) { res_mr[etaL]=res; g_ml[etaL]=g; j_mr[etaL]=j; }
	}
	for (i=1; i<nx-1; i++) {
		res=abs( D_yPG[i][0] *phiG[i][0]   /ye[0] +(FBG[i]-D_yNG[i][0]/ye[0] )*phiB_BG[i]);
		if (res>res_mb[etaL]) { res_mb[etaL]=res; g_ml[etaL]=g; i_mb[etaL]=j; }
		res=abs(-D_yNG[i][ny]*phiG[i][ny-1]/ye[ny]+(FTG[i]+D_yPG[i][ny]/ye[ny])*phiB_TG[i]);
		if (res>res_mt[etaL]) { res_mt[etaL]=res; g_ml[etaL]=g; i_mt[etaL]=j; }
	}
	
}
//======================================================================================//

//======================================================================================//
//++ function to solve NDA problem ++++++++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void NDAsolution(int preconditioner, int etaL, int g, std::vector< std::vector<double> > B_up, std::vector< std::vector<double> > S_s, std::vector< std::vector<double> > S_f)
{
	int i, j, p, c, N_ukn, N_c, N_xf, N_yf;
	double res=0;
	double **phiG, *phiB_LG, *phiB_RG, *phiB_BG, *phiB_TG, *FLG, *FRG, *FBG, *FTG;
	double **D_xPG, **D_xNG, **D_yPG, **D_yNG, **sigmaTG, **sigmaSG, **nuSigmaFG, **chiG, **s_extG;
	
	phiG=phi[etaL][g];
	phiB_LG=phiB_L[etaL][g]; phiB_RG=phiB_R[etaL][g];
	phiB_BG=phiB_B[etaL][g]; phiB_TG=phiB_T[etaL][g];
	
	FLG=FL[etaL][g]; FRG=FR[etaL][g];
	FBG=FL[etaL][g]; FTG=FR[etaL][g];
	
	D_xPG=D_xP[etaL][g]; D_xNG=D_xN[etaL][g];
	D_yPG=D_yP[etaL][g]; D_yNG=D_yN[etaL][g];
	
	sigmaTG=sigmaT[etaL][g];
	sigmaSG=sigmaS[etaL][g][g];
	nuSigmaFG=nuSigmaF[etaL][g];
	chiG=chi[etaL][g];
	s_extG=s_ext[etaL][g];
	
	temfile<<"NDA Solution Grid "<<etaL<<" group "<<g<<endl;
	
	N_c=nx*ny;
	N_xf=2*nx;
	N_yf=2*ny;
	N_ukn=nx*ny+2*nx+2*ny;
	
	VectorXd x(N_ukn);
	VectorXd b(N_ukn);
    SpMat A(N_ukn,N_ukn);
	
	p=0; // initialize solution vector guess to transport solution
	for (j=0; j<ny; j++) {
		for (i=0; i<nx; i++) { x[p]=phiG[i][j]; p++; } // cell centre flux
	}
	for (j=0; j<ny; j++) { x[p]=phiB_LG[j]; p++; } // Left   BC Flux
	for (i=0; i<nx; i++) { x[p]=phiB_BG[i]; p++; } // Bottom BC Flux
	for (j=0; j<ny; j++) { x[p]=phiB_RG[j]; p++; } // Right BC Flux
	for (i=0; i<nx; i++) { x[p]=phiB_TG[i]; p++; } // Top   BC Flux
	
	// Assign matrix A and vector b
	p=0; // ++++++++ Central Cell ++++++++
	for ( j=0; j<ny; j++ ) {
		for ( i=0; i<nx; i++ ) {
			A.insert(p,p)=hy[j]*D_xNG[i+1][j]/xe[i+1] + hy[j]*D_xPG[i][j]/xe[i] + hx[i]*D_yNG[i][j+1]/ye[j+1] + hx[i]*D_yPG[i][j]/ye[j] +
			hx[i]*hy[j]*(sigmaTG[i][j]-sigmaSG[i][j]);
			
			//temfile<<D_xNG[i+1][j]<<" "<<D_xPG[i][j]<<" "<<D_yNG[i][j+1]<<" "<<D_yPG[i][j]<<endl;
			
			b[p]=hx[i]*hy[j]*(S_s[i][j]+B_up[i][j]+chiG[i][j]*S_f[i][j]/k_eff+s_extG[i][j]); p++;
		}
	}
	//temfile<<A<<endl;
	// ++++++++ Periphery Cells ++++++++
	p=0;
	if ( nx==1 and ny==1 ) { // Single Cell
		A.insert(0,1)=-hy[0]*D_xNG[0][0]/xe[0];   A.insert(0,2)=-hx[0]*D_yNG[0][0]/ye[0];
		A.insert(0,3)=-hy[0]*D_xPG[1][0]/xe[1];   A.insert(0,4)=-hx[0]*D_yPG[0][1]/ye[1];
	}
	else if ( nx==1 ) { // One Cell Wide
		A.insert(p,N_c)      =-hy[0]*D_xNG[0][0]/xe[0];   A.insert(p,N_c+ny)=-hx[0]*D_yNG[0][0]/ye[0];
		A.insert(p,N_c+ny+nx)=-hy[0]*D_xPG[1][0]/xe[1];   A.insert(p,1)     =-hx[0]*D_yPG[0][1]/ye[1];
		p++; // Bottom Cell 
		for ( j=1; j<ny-1; j++ ) { // Middle Cells
			A.insert(p,N_c+j)      =-hy[j]*D_xNG[0][j]/xe[0];   A.insert(p,j-1)=-hx[0]*D_yNG[0][j]  /ye[j];
			A.insert(p,N_c+ny+nx+j)=-hy[j]*D_xPG[1][j]/xe[1];   A.insert(p,j+1)=-hx[0]*D_yPG[0][j+1]/ye[j+1];
			p++;
		}
		A.insert(p,N_c+ny-1)     =-hy[ny-1]*D_xNG[0][ny-1]/xe[0];   A.insert(p,ny-2)           =-hx[0]*D_yNG[0][ny-1]/ye[ny-1];
		A.insert(p,N_c+2*ny+nx-1)=-hy[ny-1]*D_xPG[1][ny-1]/xe[1];   A.insert(p,N_c+2*ny+2*nx-1)=-hx[0]*D_yPG[0][ny]  /ye[ny];
		p++; // Top Cell	
	}
	else if ( ny==1 ) { // One Cell Tall
		A.insert(p,N_c)=-hy[0]*D_xNG[0][0]/xe[0];   A.insert(p,N_c+ny)     =-hx[0]*D_yNG[0][0]/ye[0];
		A.insert(p,1)  =-hy[0]*D_xPG[1][0]/xe[1];   A.insert(p,N_c+2*ny+nx)=-hx[0]*D_yPG[0][1]/ye[1];
		p++; // Left Cell
		for ( i=1; i<nx-1; i++ ) { // Middle Cell
			A.insert(p,i-1)=-hy[0]*D_xNG[i][0]  /xe[i];    A.insert(p,N_c+ny+i)     =-hx[i]*D_yNG[i][0]/ye[0];
			A.insert(p,i+1)=-hy[0]*D_xPG[i+1][0]/xe[i+1];  A.insert(p,N_c+2*ny+nx+i)=-hx[i]*D_yPG[i][1]/ye[1];
			p++;
		}
		A.insert(p,nx-2)       =-hy[0]*D_xNG[nx-1][0]/xe[nx-1];  A.insert(p,N_c+ny+nx-1)    =-hx[nx-1]*D_yNG[nx-1][0]/ye[0];
		A.insert(p,N_c+2*ny+nx-1)=-hy[0]*D_xPG[nx][0]  /xe[nx];    A.insert(p,N_c+2*ny+2*nx-1)=-hx[nx-1]*D_yPG[nx-1][1]/ye[1];
		p++; // Right Cell
	}
	else {
		A.insert(p,N_c)=-hy[0]*D_xNG[0][0]/xe[0];   A.insert(p,N_c+ny)=-hx[0]*D_yNG[0][0]/ye[0];
		A.insert(p,1)  =-hy[0]*D_xPG[1][0]/xe[1];   A.insert(p,nx)    =-hx[0]*D_yPG[0][1]/ye[1];
		p++; // Bottom Left Corner
		for ( i=1; i<nx-1; i++ ) { // Bottom
			A.insert(p,p-1)=-hy[0]*D_xNG[i][0]  /xe[i];    A.insert(p,N_c+ny+i)=-hx[i]*D_yNG[i][0]/ye[0];
			A.insert(p,p+1)=-hy[0]*D_xPG[i+1][0]/xe[i+1];  A.insert(p,nx+i)    =-hx[i]*D_yPG[i][1]/ye[1];
			p++;
		}
		A.insert(p,nx-2)     =-hy[0]*D_xNG[nx-1][0]/xe[nx-1]; A.insert(p,N_c+ny+nx-1)=-hx[nx-1]*D_yNG[nx-1][0]/ye[0];
		A.insert(p,N_c+ny+nx)=-hy[0]*D_xPG[nx][0]  /xe[nx];   A.insert(p,2*nx-1)     =-hx[nx-1]*D_yPG[nx-1][1]/ye[1];
		p++; // Bottom Right Corner
		for ( j=1; j<ny-1; j++ ) { // Middle
			i=0; // Left Side
			A.insert(p,N_c+j)=-hy[j]*D_xNG[i][j]  /xe[i];    A.insert(p,p-nx)=-hx[i]*D_yNG[i][j]  /ye[j];
			A.insert(p,p+1)  =-hy[j]*D_xPG[i+1][j]/xe[i+1];  A.insert(p,p+nx)=-hx[i]*D_yPG[i][j+1]/ye[j+1];
			p++;
			for ( i=1; i<nx-1; i++ ) { // Centre Cells
				A.insert(p,p-1)=-hy[j]*D_xNG[i][j]  /xe[i];    A.insert(p,p-nx)=-hx[i]*D_yNG[i][j]  /ye[j];
				A.insert(p,p+1)=-hy[j]*D_xPG[i+1][j]/xe[i+1];  A.insert(p,p+nx)=-hx[i]*D_yPG[i][j+1]/ye[j+1];
				p++;
			}
			i=nx-1; // Right Side
			A.insert(p,p-1)        =-hy[j]*D_xNG[i][j]  /xe[i];    A.insert(p,p-nx)=-hx[i]*D_yNG[i][j]  /ye[j];
			A.insert(p,N_c+ny+nx+j)=-hy[j]*D_xPG[i+1][j]/xe[i+1];  A.insert(p,p+nx)=-hx[i]*D_yPG[i][j+1]/ye[j+1];
			p++;
		}
		j=ny-1; // Top Left Corner
		A.insert(p,N_c+j)=-hy[j]*D_xNG[0][j]/xe[0];   A.insert(p,p-nx)       =-hx[0]*D_yNG[0][j]  /ye[j];
		A.insert(p,p+1)  =-hy[j]*D_xPG[1][j]/xe[1];   A.insert(p,N_c+2*ny+nx)=-hx[0]*D_yPG[0][j+1]/ye[j+1];
		p++;
		for ( i=1; i<nx-1; i++ ) { // Top
			A.insert(p,p-1)=-hy[j]*D_xNG[i][j]  /xe[i];    A.insert(p,p-nx)         =-hx[i]*D_yNG[i][j]/ye[j];
			A.insert(p,p+1)=-hy[j]*D_xPG[i+1][j]/xe[i+1];  A.insert(p,N_c+2*ny+nx+i)=-hx[i]*D_yPG[i][j+1]/ye[j+1];
			p++;
		}
		i=nx-1; // Top Right Corner
		A.insert(p,p-1)          =-hy[j]*D_xNG[i][j]  /xe[i];    A.insert(p,p-nx)         =-hx[i]*D_yNG[i][j]/ye[j];
		A.insert(p,N_c+2*ny+nx-1)=-hy[j]*D_xPG[i+1][j]/xe[i+1];  A.insert(p,N_c+2*ny+nx+i)=-hx[i]*D_yPG[i][j+1]/ye[j+1];
		p++;
	}
	
	// ++++++++ Boundary Conditions ++++++++
	p=N_c;
	for ( j=0; j<ny; j++ ) {
		A.insert(p,j*nx)= D_xPG[0][j]/xe[0];   A.insert(p,N_c+j)=FLG[j]-D_xNG[0][j]/xe[0];   b[p]=0; p++; // Left   BC
	}
	for ( i=0; i<nx; i++ ) {
		A.insert(p,i)= D_yPG[i][0]/ye[0];   A.insert(p,N_c+ny+i)=FBG[i]-D_yNG[i][0]/ye[0];   b[p]=0; p++; // Bottom BC
	}
	for ( j=0; j<ny; j++ ) {
		A.insert(p,(j+1)*nx-1)=-D_xNG[nx][j]/xe[nx];   A.insert(p,N_c+ny+nx+j)=FRG[j]+D_xPG[nx][j]/xe[nx];   b[p]=0; p++; // Right  BC
	}
	for ( i=0; i<nx; i++ ) {
		A.insert(p,(ny-1)*nx+i)=-D_yNG[i][ny]/ye[ny];  A.insert(p,N_c+2*ny+nx+i)=FTG[i]+D_yPG[i][ny]/ye[ny]; b[p]=0; p++; // Top    BC
	}
	
	A.prune(1e-17, 10);
	//temfile<<A<<endl;
	
	// Solve Ax=b iteratively with BiCGSTAB
	if ( preconditioner==2 ) {
		t_pc=clock(); // start low-order timer
		BiCGSTAB<SpMat> solver;
		solver.setTolerance(epsilon_lo);     // set convergence criteria 
		solver.setMaxIterations(maxiter_lo); // set the max number of lo iterations
		solver.compute(A);
		
		t_pc=clock()-t_pc; // stop low-order timer
		dt_pc[etaL].push_back(((double)t_pc)/CLOCKS_PER_SEC); // add high-order solution time to vector
		
		x = solver.solveWithGuess(b,x);
		err_lo[etaL].push_back(solver.error());      // error in lo solution
		num_logm[etaL].push_back(solver.iterations()); // number of lo iterations
	}
	else if ( preconditioner==3 ) {
		t_pc=clock(); // start low-order timer
		BiCGSTAB<SpMat,IncompleteLUT<double> > solver;
		solver.preconditioner().setFillfactor(11);
		solver.setTolerance(epsilon_lo);     // set convergence criteria 
		solver.setMaxIterations(maxiter_lo); // set the max number of lo iterations
		solver.compute(A);
		
		t_pc=clock()-t_pc; // stop low-order timer
		dt_pc[etaL].push_back(((double)t_pc)/CLOCKS_PER_SEC); // add high-order solution time to vector
		
		x = solver.solveWithGuess(b,x);
		err_lo[etaL].push_back(solver.error());      // error in lo solution
		num_logm[etaL].push_back(solver.iterations()); // number of lo iterations
	}
	else cout<<">>Preconditioner Type Error \n";
	num_sol++;
	
	//cout<<"x "<<x<<endl;
	
	// set solution back to problem values
	p=0;
	for (j=0; j<ny; j++) {
		for (i=0; i<nx; i++) { phiG[i][j]=x[p]; p++; } // Cell Centre Flux
	}
	for (j=0; j<ny; j++) { phiB_LG[j]=x[p]; p++; } // Left   BC Flux
	for (i=0; i<nx; i++) { phiB_BG[i]=x[p]; p++; } // Bottom BC Flux
	for (j=0; j<ny; j++) { phiB_RG[j]=x[p]; p++; } // Right BC Flux
	for (i=0; i<nx; i++) { phiB_TG[i]=x[p]; p++; } // Top   BC Flux
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Calculate current ///////////////////////////////////////////////////////////////////////////////////////
	// Boundary currents ///////////////////////////////////////////////////////////////////////////////////////
	for (j=0; j<ny; j++) j_x[etaL][g][0][j]=FLG[j]*phiB_LG[j]; // Left boundary current J_x
	for (i=0; i<nx; i++) j_y[etaL][g][i][0]=FBG[i]*phiB_BG[i]; // Bottom boundary current J_y
	// Inner currents J_x J_y
	for (j=0; j<ny; j++) for (i=1; i<nx; i++) j_x[etaL][g][i][j]=-(D_xPG[i][j]*phiG[i][j]-D_xNG[i][j]*phiG[i-1][j])/xe[i];
	for (j=1; j<ny; j++) for (i=0; i<nx; i++) j_y[etaL][g][i][j]=-(D_yPG[i][j]*phiG[i][j]-D_yNG[i][j]*phiG[i][j-1])/ye[j];
	// Boundary currents
	for (j=0; j<ny; j++) j_x[etaL][g][nx][j]=FRG[j]*phiB_RG[j]; // Right boundary current
	for (i=0; i<nx; i++) j_y[etaL][g][i][ny]=FTG[i]*phiB_TG[i]; // Top boundary current
	
	i=0; j=0; // Calculate residuals of matrix equations
	res=abs(-hy[j]*D_xPG[i+1][j]*phiG[i+1][j]/xe[i+1]-hx[i]*D_yPG[i][j+1]*phiG[i][j+1]/ye[j+1]
	+(hy[j]*D_xNG[i+1][j]/xe[i+1]+hy[j]*D_xPG[i][j]/xe[i]+hx[i]*D_yNG[i][j+1]/ye[j+1]+hx[i]*D_yPG[i][j]/ye[j]
	+hx[i]*hy[j]*(sigmaTG[i][j]-sigmaSG[i][j]))*phiG[i][j]-hy[j]*D_xNG[i][j]*phiB_LG[j]/xe[i]-hx[i]*D_yNG[i][j]*phiB_BG[i]/ye[j]
	-hx[i]*hy[j]*(S_s[i][j]+B_up[i][j]+chiG[i][j]*S_f[i][j]/k_eff+s_extG[i][j]));
	if (res>res_mbal[etaL]) { res_mbal[etaL]=res; i_mbal[etaL]=i; j_mbal[etaL]=j; }
	for (j=1; j<ny-1; j++) {
		i=0;
		res=abs(-hy[j]*D_xPG[i+1][j]*phiG[i+1][j]/xe[i+1]-hx[i]*D_yPG[i][j+1]*phiG[i][j+1]/ye[j+1]
		+(hy[j]*D_xNG[i+1][j]/xe[i+1]+hy[j]*D_xPG[i][j]/xe[i]+hx[i]*D_yNG[i][j+1]/ye[j+1]+hx[i]*D_yPG[i][j]/ye[j]
		+hx[i]*hy[j]*(sigmaTG[i][j]-sigmaSG[i][j]))*phiG[i][j]
		-hy[j]*D_xNG[i][j]*phiB_LG[j]/xe[i]-hx[i]*D_yNG[i][j]*phiG[i][j-1]/ye[j]
		-hx[i]*hy[j]*(S_s[i][j]+B_up[i][j]+chiG[i][j]*S_f[i][j]/k_eff+s_extG[i][j]));
		if (res>res_mbal[etaL]) { res_mbal[etaL]=res; i_mbal[etaL]=i; j_mbal[etaL]=j; }
		for (i=1; i<nx-1; i++) {
			res=abs(-hy[j]*D_xPG[i+1][j]*phiG[i+1][j]/xe[i+1]-hx[i]*D_yPG[i][j+1]*phiG[i][j+1]/ye[j+1]
			+(hy[j]*D_xNG[i+1][j]/xe[i+1]+hy[j]*D_xPG[i][j]/xe[i]+hx[i]*D_yNG[i][j+1]/ye[j+1]+hx[i]*D_yPG[i][j]/ye[j]
			+hx[i]*hy[j]*(sigmaTG[i][j]-sigmaSG[i][j]))*phiG[i][j]
			-hy[j]*D_xNG[i][j]*phiG[i-1][j]/xe[i]-hx[i]*D_yNG[i][j]*phiG[i][j-1]/ye[j]
			-hx[i]*hy[j]*(S_s[i][j]+B_up[i][j]+chiG[i][j]*S_f[i][j]/k_eff+s_extG[i][j]));
			if (res>res_mbal[etaL]) { res_mbal[etaL]=res; i_mbal[etaL]=i; j_mbal[etaL]=j; }
		}
		i=nx-1;
		res=abs(-hy[j]*D_xPG[i+1][j]*phiB_RG[j]/xe[i+1]-hx[i]*D_yPG[i][j+1]*phiG[i][j+1]/ye[j+1]
		+(hy[j]*D_xNG[i+1][j]/xe[i+1]+hy[j]*D_xPG[i][j]/xe[i]+hx[i]*D_yNG[i][j+1]/ye[j+1]+hx[i]*D_yPG[i][j]/ye[j]
		+hx[i]*hy[j]*(sigmaTG[i][j]-sigmaSG[i][j]))*phiG[i][j]
		-hy[j]*D_xNG[i][j]*phiG[i-1][j]/xe[i]-hx[i]*D_yNG[i][j]*phiG[i][j-1]/ye[j]
		-hx[i]*hy[j]*(S_s[i][j]+B_up[i][j]+chiG[i][j]*S_f[i][j]/k_eff+s_extG[i][j]));
		if (res>res_mbal[etaL]) { res_mbal[etaL]=res; i_mbal[etaL]=i; j_mbal[etaL]=j; }
	}
	i=nx-1; j=ny-1;
	res=abs(-hy[j]*D_xPG[i+1][j]*phiB_RG[j]/xe[i+1]-hx[i]*D_yPG[i][j+1]*phiB_TG[i]/ye[j+1]
	+(hy[j]*D_xNG[i+1][j]/xe[i+1]+hy[j]*D_xPG[i][j]/xe[i]+hx[i]*D_yNG[i][j+1]/ye[j+1]+hx[i]*D_yPG[i][j]/ye[j]
	+hx[i]*hy[j]*(sigmaTG[i][j]-sigmaSG[i][j]))*phiG[i][j]-hy[j]*D_xNG[i][j]*phiG[i-1][j]/xe[i]-hx[i]*D_yNG[i][j]*phiG[i][j-1]/ye[j]
	-hx[i]*hy[j]*(S_s[i][j]+B_up[i][j]+chiG[i][j]*S_f[i][j]/k_eff+s_extG[i][j]));
	if (res>res_mbal[etaL]) { res_mbal[etaL]=res; g_mbal[etaL]=g; i_mbal[etaL]=i; j_mbal[etaL]=j; }
	
	for (j=1; j<ny-1; j++) {
		res=abs( D_xPG[0][j] *phiG[0][j]   /xe[0] +(FLG[j]-D_xNG[0][j] /xe[0])*phiB_LG[j]);
		if (res>res_ml[etaL]) { res_ml[etaL]=res; g_ml[etaL]=g; j_ml[etaL]=j; }
		res=abs(-D_xNG[nx][j]*phiG[nx-1][j]/xe[nx]+(FRG[j]+D_xPG[nx][j]/xe[nx])*phiB_RG[j]);
		if (res>res_mr[etaL]) { res_mr[etaL]=res; g_mr[etaL]=g; j_mr[etaL]=j; }
	}
	for (i=1; i<nx-1; i++) {
		res=abs( D_yPG[i][0] *phiG[i][0]   /ye[0] +(FBG[i]-D_yNG[i][0]/ye[0] )*phiB_BG[i]);
		if (res>res_mb[etaL]) { res_mb[etaL]=res; g_mb[etaL]=g; i_mb[etaL]=j; }
		res=abs(-D_yNG[i][ny]*phiG[i][ny-1]/ye[ny]+(FTG[i]+D_yPG[i][ny]/ye[ny])*phiB_TG[i]);
		if (res>res_mt[etaL]) { res_mt[etaL]=res; g_mt[etaL]=g; i_mt[etaL]=j; }
	}
	
}
//======================================================================================//

