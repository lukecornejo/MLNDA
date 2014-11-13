#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <time.h>  /* clock_t, clock, CLOCKS_PER_SEC */
#include <vector>
#include "ioMLNDA2.1.h"
#include "scMLNDA2.1.h"
#include "solverMLNDA2.1.h"
using namespace std;

/*
  2D Transport Solution With Step Characteristics and Non-linear Diffusion Acceleration
  
  Programer : Luke Cornejo
  Version   : 1.0
  Date      : 2-19-14
  
                         Changes
  ******************************************************

  
  To complile
  Windows
c++ -I .\eigen SC_MLNDA2.1.c ioMLNDA2.1.c solverMLNDA2.1.c scMLNDA2.1.c -o SC_MLNDA2.1.exe -std=c++0x
  Linux
c++ -I ./eigen -g SC_MLNDA2.1.c ioMLNDA2.1.c solverMLNDA2.1.c scMLNDA2.1.c -o SC_MLNDA2.1.exe
  
  hpc
c++ -I .\eigen -O3 SC_MLNDA2.1.c ioMLNDA2.1.c solverMLNDA2.1.c scMLNDA2.1.c -o SC_MLNDA2.1b.exe
-O
-O2
-O3 better
-Os
-O3 -ffast-math


icpc -I .\eigen -O3 SC_MLNDA2.1.c ioMLNDA2.1.c solverMLNDA2.1.c scMLNDA2.1.c -o SC_MLNDA2.1a.exe
-O1
-O2
-O3 best
-xO
-fast

  Boundary type options
  1: incoming according to side
  2: incoming according to angle
  3: reflective on all sides
  4: reflective on Left and Bottom, incoming on Right and Top according to side
  5: reflective on Right and Top, incoming on Left and Bottom according to side
  
  BC input order
  according to side : Left, Bottom, Right, Top
  according to angle: quad 1, quad 2, quad 3, quad 4
  
  Quadratures
  S4, S6, S8, S12, S16
  Q20, Q20
  

*/

const double pi=3.141592654;
int maxiter_lo=1000;
bool o_angular=false, KE_problem=false, write_xs=false; // option variables
double epsilon_si=1e-5, epsilon_kho=1e-5, epsilon_lo=1e-10, *epsilon_phi, *epsilon_keff; // default tolerances
vector<int> stop_phi, preconditioner_type; // limit number of iterations
int    kbc; // Kind of BC
int nx, ny;
double *x, *y, *hx, *hy, *xe, *ye; // grid arrays
int sn, N=8;                     // quadrature
double *mu, *eta, *xi, *w;   // quadrature
// Cross-Section Data
int *Ng, eta_star, **omegaP;
double ****sigmaT, *****sigmaS, ****nuSigmaF, ****chi, ****s_ext; // material arrays
double **bcL, **bcR, **bcB, **bcT;
// Solution
double k_eff=1.0, delta=0.5; // K-effective
double ****psiL, ****psiR, ****psiB, ****psiT;
double ***phiT, ***phi_xT, ***phi_yT, ***j_xT, ***j_yT; // scalar flux and current from transport
double ****phi , ****j_x , ****j_y; // NDA scalar flux and current solution
double ***phiB_L, ***phiB_R, ***phiB_B, ***phiB_T; // edge scalar flux from NDA
double **phiH, ***phiL; // One-group flux and scalar flux from last iteration
// Problem Factors
double **f; // correction factor
double ***FL, ***FR, ***FB, ***FT;
double ***D, ***D_x, ***D_y; // Diffusion coefficients
double ***D_xT, ***D_yT; // Tilde, Positive and Negative
double ****D_xP, ****D_xN, ****D_yP, ****D_yN; // Positive and Negative
// residual data
double res_ho;
int i_ho, j_ho;
double *res_bal;
int *i_bal, *j_bal, *g_bal;
double *res_mbal, *res_ml, *res_mr, *res_mb, *res_mt;
int *i_mbal, *j_mbal, *j_ml, *j_mr, *i_mb, *i_mt;
int *g_mbal, *g_ml, *g_mr, *g_mb, *g_mt;
// iteration data
int n_iterations, num_sol;
vector<int> num_mtot;
vector<double> rho_ho, rho_kho, norm_ho, norm_kho, dt_ho;
vector< vector<double> > rho_phi, rho_keff, norm_phi, norm_keff, err_lo, dt_lo, dt_pc;
vector< vector<int> > num_losi, num_logm, num_grid;
ofstream outfile;
ofstream datfile;
ofstream temfile;
clock_t t;
clock_t t_lo;
clock_t t_ho;


//======================================================================================//
//++ function to initialize problem memory space +++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
int initialize()
{
	int i, j, m, g, k;
	// initialize memory space
	res_bal=new double[eta_star];
	i_bal=new int[eta_star]; j_bal=new int[eta_star]; g_bal=new int[eta_star];
	res_mbal=new double[eta_star]; i_mbal=new int[eta_star]; j_mbal=new int[eta_star]; g_mbal=new int[eta_star];
	res_ml=new double[eta_star]; j_ml=new int[eta_star]; g_ml=new int[eta_star];
	res_mr=new double[eta_star]; j_mr=new int[eta_star]; g_mr=new int[eta_star];
	res_mb=new double[eta_star]; i_mb=new int[eta_star]; g_mb=new int[eta_star];
	res_mt=new double[eta_star]; i_mt=new int[eta_star]; g_mt=new int[eta_star];
	
	for (k=0; k<eta_star; k++) {
		vector< double > rowd;
		rho_phi.push_back(rowd);
		rho_phi[k].push_back(0.5);
		rho_keff.push_back(rowd);
		rho_keff[k].push_back(0.5);
		norm_phi.push_back(rowd);
		norm_phi[k].push_back(1);
		norm_keff.push_back(rowd);
		norm_keff[k].push_back(1);
		err_lo.push_back(rowd);
		dt_lo.push_back(rowd);
		dt_pc.push_back(rowd);
		vector< int > rowi;
		num_losi.push_back(rowi);
		num_losi[k].push_back(0);
		num_logm.push_back(rowi);
		num_grid.push_back(rowi);
	}
	
	// Multi-group variables
	
	// cell average flux
	phiL=new double**[Ng[0]];
	D   =new double**[Ng[0]];
	for (g=0; g<Ng[0]; g++) {
		phiL[g]=new double*[nx];
		D[g]   =new double*[nx];
		for (i=0; i<nx; i++) {
			phiL[g][i]=new double[ny];
			D[g][i]   =new double[ny];
		}
	}
	
	// cell edge values on x grid
	D_x   =new double**[Ng[0]];
	D_xT  =new double**[Ng[0]];
	for (g=0; g<Ng[0]; g++) {
		D_x[g]   =new double*[nx+1];
		D_xT[g]  =new double*[nx+1];
		for (i=0; i<nx+1; i++) {
			D_x[g][i]   =new double[ny];
			D_xT[g][i]  =new double[ny];
		}
	}
	
	// cell edge flux on y grid
	D_y   =new double**[Ng[0]];
	D_yT  =new double**[Ng[0]];
	for (g=0; g<Ng[0]; g++) {
		D_y[g]   =new double*[nx];
		D_yT[g]  =new double*[nx];
		for (i=0; i<nx; i++) {
			D_y[g][i]   =new double[ny+1];
			D_yT[g][i]  =new double[ny+1];
		}
	}
	
	// Multi-grid variables 
	FL=new double**[eta_star];
	FR=new double**[eta_star];
	FB=new double**[eta_star];
	FT=new double**[eta_star];
	
	phiB_L=new double**[eta_star];
	phiB_R=new double**[eta_star];
	phiB_B=new double**[eta_star];
	phiB_T=new double**[eta_star];
	for (k=0; k<eta_star; k++) { // Energy Grids
		FL[k]=new double*[Ng[k]];
		FR[k]=new double*[Ng[k]];
		FB[k]=new double*[Ng[k]];
		FT[k]=new double*[Ng[k]];
	
		phiB_L[k]=new double*[Ng[k]];
		phiB_R[k]=new double*[Ng[k]];
		phiB_B[k]=new double*[Ng[k]];
		phiB_T[k]=new double*[Ng[k]];
		for (g=0; g<Ng[k]; g++) { // Energy Groups
			FL[k][g]=new double[ny];
			FR[k][g]=new double[ny];
			FB[k][g]=new double[nx];
			FT[k][g]=new double[nx];
			
			phiB_L[k][g]=new double[ny];
			phiB_R[k][g]=new double[ny];
			phiB_B[k][g]=new double[nx];
			phiB_T[k][g]=new double[nx];
		}
	}
	
	phi=new double***[eta_star];
	for (k=0; k<eta_star; k++) { // Energy Grids
		phi[k]=new double**[Ng[k]];
		for (g=0; g<Ng[k]; g++) { // Energy Groups
			phi[k][g]=new double*[nx];
			for (i=0; i<nx; i++) phi[k][g][i]=new double[ny];
		}
	}
	phiH=phi[eta_star-1][0]; // One-group flux
	
	f=new double*[eta_star];
	for (k=1; k<eta_star; k++) f[k]=new double[Ng[k]];
	
	j_x =new double***[eta_star];
	D_xP=new double***[eta_star];
	D_xN=new double***[eta_star];
	for (k=0; k<eta_star; k++) { // Energy Grids
		j_x[k] =new double**[Ng[k]];
		D_xP[k]=new double**[Ng[k]];
		D_xN[k]=new double**[Ng[k]];
		for (g=0; g<Ng[k]; g++) { // Energy Groups
			j_x[k][g] =new double*[nx+1];
			D_xP[k][g]=new double*[nx+1];
			D_xN[k][g]=new double*[nx+1];
			for (i=0; i<nx+1; i++) {
				j_x[k][g][i] =new double[ny];
				D_xP[k][g][i]=new double[ny];
				D_xN[k][g][i]=new double[ny];
			}
		}
	}
	
	j_y =new double***[eta_star];
	D_yP=new double***[eta_star];
	D_yN=new double***[eta_star];
	for (k=0; k<eta_star; k++) { // Energy Grids
		j_y[k] =new double**[Ng[k]];
		D_yP[k]=new double**[Ng[k]];
		D_yN[k]=new double**[Ng[k]];
		for (g=0; g<Ng[k]; g++) { // Energy Groups
			j_y[k][g] =new double*[nx];
			D_yP[k][g]=new double*[nx];
			D_yN[k][g]=new double*[nx];
			for (i=0; i<nx; i++) {
				j_y[k][g][i] =new double[ny+1];
				D_yP[k][g][i]=new double[ny+1];
				D_yN[k][g][i]=new double[ny+1];
			}
		}
	}
	
	// find values for D's
	for (g=0; g<Ng[0]; g++) {
		for ( i=0; i<nx; i++ ) for ( j=0; j<ny; j++ ) D[g][i][j]=1/(3*sigmaT[0][g][i][j]); // Diffusion Coefficient. Isotropic scattering so sigma_tr = sigma_t
		for ( i=1; i<nx; i++ ) for ( j=0; j<ny; j++ ) D_x[g][i][j]=2*xe[i]*D[g][i-1][j]*D[g][i][j]/(D[g][i-1][j]*hx[i]+D[g][i][j]*hx[i-1]); // X Grid D
		for ( j=0; j<ny; j++ ) { // Top and Bottom
			D_x[g][0][j] =2*xe[0] *D[g][0][j]   /hx[0];
			D_x[g][nx][j]=2*xe[nx]*D[g][nx-1][j]/hx[nx-1];
		}
		for ( i=0; i<nx; i++ ) for ( j=1; j<ny; j++ ) D_y[g][i][j]=2*ye[j]*D[g][i][j-1]*D[g][i][j]/(D[g][i][j-1]*hy[j]+D[g][i][j]*hy[j-1]); // Y Grid D
		for ( i=0; i<nx; i++ ) { // Left and Right
			D_y[g][i][0] =2*ye[0] *D[g][i][0]   /hy[0];
			D_y[g][i][ny]=2*ye[ny]*D[g][i][ny-1]/hy[ny-1];
		}
	}
	return 0;
}
//======================================================================================//


//======================================================================================//
//++ Recursive function to solve One Group +++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void oneGroupSolFS() {
	int i, j;
	double norm_p, res;
	int etaL, g;
	
	etaL=eta_star-1; g=0;
	res_mbal[etaL]=0.0;
	
	cout<<string((etaL+1)*5,' ')<<"Solve Low-order One-group Eq. \n";
	// Solve NDA on each energy group
	//write_grid_average_dat(phi, 16, temfile);
	GNDAsolutionFS(preconditioner_type[etaL]); //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//write_grid_average_dat(phi, 16, temfile);
	
	// Calculate Equation Residuals
	res_bal[etaL]=0.0;
	for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) { // Balance in Cell
			res=abs(hy[j]*(j_x[etaL][g][i+1][j]-j_x[etaL][g][i][j])+hx[i]*(j_y[etaL][g][i][j+1]-j_y[etaL][g][i][j])+
			hx[i]*hy[j]*(sigmaT[etaL][g][i][j]-sigmaS[etaL][g][g][i][j]-nuSigmaF[etaL][g][i][j]/k_eff)*phi[etaL][g][i][j]
			-hx[i]*hy[j]*s_ext[etaL][g][i][j]);
			if ( res>res_bal[etaL] ) { res_bal[etaL]=res; i_bal[etaL]=i; j_bal[etaL]=j;  g_bal[etaL]=g; }
		}
	}
	
	num_grid[etaL].back()++;
	num_losi[etaL].push_back(1);
	
	norm_phi[etaL].push_back(0.0);
	rho_phi[etaL].push_back(0.0);
	norm_keff[etaL].push_back(0.0);
	rho_keff[etaL].push_back(0.0);
	
}
//======================================================================================//


//======================================================================================//
//++ Normalize the Eigenfunction to 1 ++++++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void normalizeEigen(int etaL) {
	int i, j, g;
	double sum=0.0;
	double ***phiG, ***j_xG, ***j_yG, **phiB_BG, **phiB_TG, **phiB_LG, **phiB_RG;
	phiG=phi[etaL]; j_xG=j_x[etaL]; j_yG=j_y[etaL];
	phiB_BG=phiB_B[etaL]; phiB_TG=phiB_T[etaL]; phiB_LG=phiB_L[etaL]; phiB_RG=phiB_R[etaL];
	// Find coefficient
	for (g=0; g<Ng[etaL]; g++) for (i=0; i<nx; i++) for (j=0; j<ny; j++) sum+=hx[i]*hy[j]*phiG[g][i][j];
	// Normalize
	for (g=0; g<Ng[etaL]; g++) {
		for (i=0; i<nx; i++) for (j=0; j<ny; j++) phiG[g][i][j]/=sum;
		for (i=0; i<nx+1; i++) for (j=0; j<ny; j++) j_xG[g][i][j]/=sum;
		for (i=0; i<nx; i++) for (j=0; j<ny+1; j++) j_yG[g][i][j]/=sum;
		for (i=0; i<nx; i++) {
			phiB_BG[g][i]/=sum;
			phiB_TG[g][i]/=sum;
		}
		for (j=0; j<ny; j++) {
			phiB_LG[g][j]/=sum;
			phiB_RG[g][j]/=sum;
		}
	}
	
}
//======================================================================================//

//======================================================================================//
//++ Recursive function to solve One Group +++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void oneGroupSolKE() {
	int i, j, l;
	double rho_p, rho_k, norm_p, norm_k, norm_pL, norm_kL, k_effL, res;
	double **phiG, **j_xG, **j_yG;
	double int_fission, int_jx, int_jy, int_abs, L1_norm;
	double **sigmaTG, **sigmaSG, **nuSigmaFG;
	vector< vector<double> > phi_gL;
	int etaL, g;
	
	
	etaL=eta_star-1; g=0;
	sigmaTG=sigmaT[etaL][g];
	sigmaSG=sigmaS[etaL][g][g];
	nuSigmaFG=nuSigmaF[etaL][g];
	
	phiG=phi[etaL][0];
	phi_gL.resize(nx);
	for (i=0; i<nx; i++) phi_gL[i].resize(ny);
	j_xG=j_x[etaL][g]; j_yG=j_y[etaL][g];
	
	cout<<string((etaL+1)*5,' ')<<"Solve Low-order One-group Eq. \n";
	norm_p=norm_phi[etaL].back();
	norm_k=norm_keff[etaL].back();
	l=0;
	do { // Do While Loop
		
		k_effL=k_eff;
		for (i=0; i<nx; i++) for (j=0; j<ny; j++) phi_gL[i][j]=phiG[i][j];
		norm_kL=norm_k;
		norm_pL=norm_p;
		
		// Solve NDA on each energy group
		GNDAsolutionKE(preconditioner_type[etaL], k_effL, phi_gL); //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		
		int_fission=0; int_abs=0; int_jx=0; int_jy=0;
		for (i=0; i<nx; i++) for (j=0; j<ny; j++) int_fission+=nuSigmaFG[i][j]*phiG[i][j]*hx[i]*hy[j];
		for (i=0; i<nx; i++) for (j=0; j<ny; j++) int_abs+=(sigmaTG[i][j]-sigmaSG[i][j])*phiG[i][j]*hx[i]*hy[j];
		for (j=0; j<ny; j++) int_jx+=(j_xG[nx][j]-j_xG[0][j])*hy[j];
		for (i=0; i<nx; i++) int_jy+=(j_yG[i][ny]-j_yG[i][0])*hx[i];
		
		k_eff=int_fission/(int_jx+int_jy+int_abs);
		
		norm_k=abs(k_eff-k_effL);
		rho_k=norm_k/norm_kL;
		norm_p=0;
		for (i=0; i<nx; i++) for (j=0; j<ny; j++) if ( abs(phiG[i][j]-phi_gL[i][j])>norm_p ) norm_p=abs(phiG[i][j]-phi_gL[i][j]);
		
		//write_grid_average_dat(phi, 16, temfile);
		
		rho_p=norm_p/norm_pL;
		norm_phi[etaL].push_back(norm_p);
		rho_phi[etaL].push_back(rho_p);
		norm_keff[etaL].push_back(norm_k+1e-20);
		rho_keff[etaL].push_back(rho_k);
		
		// Calculate Equation Residuals
		res_bal[etaL]=0.0;
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) { // Balance in Cell
				res=abs(hy[j]*(j_x[etaL][g][i+1][j]-j_x[etaL][g][i][j])+hx[i]*(j_y[etaL][g][i][j+1]-j_y[etaL][g][i][j])+
				hx[i]*hy[j]*(sigmaT[etaL][g][i][j]-sigmaS[etaL][g][g][i][j]-(1-delta)*nuSigmaF[etaL][g][i][j]/k_effL)*phiG[i][j]
				-hx[i]*hy[j]*(delta*nuSigmaFG[i][j]*phi_gL[i][j]/k_effL));
				if ( res>res_bal[etaL] ) { res_bal[etaL]=res; i_bal[etaL]=i; j_bal[etaL]=j;  g_bal[etaL]=g; }
			}
		}
		
		// Normalize the solution to 1
		normalizeEigen(etaL);
		
		l++;
		num_grid[etaL].back()++;
	} while ( ( norm_p>epsilon_phi[etaL]*(1/rho_p-1) or norm_k>epsilon_keff[etaL]*(1/rho_k-1) ) and l<stop_phi[etaL] );
	
	num_losi[etaL].push_back(l);
}
//======================================================================================//

//======================================================================================//
//++ Calculate the Correction Factors and corrected flux +++++++++++++++++++++++++++++++//
//======================================================================================//
void fCorrection(int etaL, int p, int i, int j, double fprod, int eta_stop) { // Recursive Function to Correct flux 
	int g;
	
	for (g=omegaP[etaL][p]; g<omegaP[etaL][p+1]; g++) {
		phi[etaL-1][g][i][j]*=fprod*f[etaL][p];
		if ( etaL>eta_stop+1 ) fCorrection(etaL-1,g,i,j,fprod*f[etaL][p], eta_stop);
	}
}
//======================================================================================//
void calcCorrection(int etaL) {
	int g, i, j, k, p, n;
	double phi_sum, fprod;
	
	// Calculate correction factors
	for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) { 
			for (k=etaL+1; k<eta_star; k++) {
				for (p=0; p<Ng[k]; p++) {
					phi_sum=0.0;
					for (g=omegaP[k][p]; g<omegaP[k][p+1]; g++) phi_sum+=phi[k-1][g][i][j];
					f[k][p]=phi[k][p][i][j]/phi_sum;
				}
			}
			fCorrection(eta_star-1, 0, i, j, 1, etaL);
		}
	}
}
//======================================================================================//


//======================================================================================//
//++ Recursive function to solve Coarse Grid +++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void coarseGridSol(int etaL) {
	int g, gg, i, j, l, p, pp, n_in;
	double rho_p, rho_k, norm_p, norm_k, norm_pL, norm_kL, k_effL, res, phiS, phiP, phiN, Ssum;
	vector< vector<double> > phi_gL, B_up, S_s, S_f;
	double ***phiG;
	
	res_mbal[etaL]=0.0;
	res_ml[etaL]=0.0;
	res_mr[etaL]=0.0;
	res_mb[etaL]=0.0;
	res_mt[etaL]=0.0;
	phiG=phi[etaL];
	
	phi_gL.resize(nx);
	B_up.resize(nx);
	S_s.resize(nx);
	S_f.resize(nx);
	for (i=0; i<nx; i++) {
		phi_gL[i].resize(ny);
		B_up[i].resize(ny);
		S_s[i].resize(ny);
		S_f[i].resize(ny);
	}
	cout<<string((etaL+1)*5,' ')<<"Solve Low-order Eq. on Grid # "<<etaL<<endl;
	norm_p=norm_phi[etaL].back();
	norm_k=norm_keff[etaL].back();
	l=0;
	do { // Do While Loop
		
		if (l!=0) calcCorrection(etaL); // calculate correction factors and update fluxes
		
		temfile<<" -- LO Grid "<<etaL<<" Iteration # "<<l<<" -- "<<endl;
		//write_grid_average_dat(phi, 16, temfile);
		//temfile<<"phiT\n";
		//write_group_average_dat(phiT, 0, 16, temfile);
		//temfile<<"D_xT \n";
		//write_group_edge_x_dat(D_xT, 0, 16, temfile);
		
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) {
				S_f[i][j]=0.0; for (g=0; g<Ng[etaL]; g++) S_f[i][j]+=nuSigmaF[etaL][g][i][j]*phiG[g][i][j];
			}
		}
		
		norm_pL=norm_p;
		norm_kL=norm_k;
		norm_p=0;
		res_bal[etaL]=0.0;
		k_effL=k_eff;
		
		for (g=0; g<Ng[etaL]; g++) { 
			//temfile<<"D_xN 1\n";
			//write_cell_edge_x_dat(D_xN[etaL][g], 16, temfile);
			
			for (i=0; i<nx; i++) for (j=0; j<ny; j++) phi_gL[i][j]=phiG[g][i][j]; // Set phi_gL to previous iteration
			for (i=0; i<nx; i++) {
				for (j=0; j<ny; j++) {
					B_up[i][j]=0.0; // Find up scattering term from previous iteration
					for (gg=g+1; gg<Ng[etaL]; gg++) B_up[i][j]+=sigmaS[etaL][gg][g][i][j]*phiG[gg][i][j];
					S_s[i][j]=0.0; 	// Find down scattering and fission production terms from current iteration
					for (gg=0; gg<g; gg++) S_s[i][j]+=sigmaS[etaL][gg][g][i][j]*phiG[gg][i][j];
				}
			}
			
			//temfile<<"D_xN 1\n";
			//write_grid_edge_x_dat(D_xN, 16, temfile);
			
			// Solve NDA on each energy group
			NDAsolution(preconditioner_type[etaL], etaL, g, B_up, S_s, S_f); //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			
			//write_grid_average_dat(phi, 16, temfile);
			
			// Calculate Equation Residuals
			for (i=0; i<nx; i++) {
				for (j=0; j<ny; j++) { // Residual Balance in Cell
					res=abs(hy[j]*(j_x[etaL][g][i+1][j]-j_x[etaL][g][i][j])+hx[i]*(j_y[etaL][g][i][j+1]-j_y[etaL][g][i][j])+
					hx[i]*hy[j]*(sigmaT[etaL][g][i][j]-sigmaS[etaL][g][g][i][j])*phiG[g][i][j]
					-hx[i]*hy[j]*(S_s[i][j]+B_up[i][j]+chi[etaL][g][i][j]*S_f[i][j]/k_eff+s_ext[etaL][g][i][j]));
					if ( res>res_bal[etaL] ) { res_bal[etaL]=res; i_bal[etaL]=i; j_bal[etaL]=j; g_bal[etaL]=g; }
				}
			}
			
			// Calculate Flux K-Infinity Norm
			for (i=0; i<nx; i++) for (j=0; j<ny; j++) if ( abs(phiG[g][i][j]-phi_gL[i][j])>norm_p ) norm_p=abs(phiG[g][i][j]-phi_gL[i][j]);
		}
		
		
		// Calculate Group Averaged values ////////////////////////////////////////////////////////////
		// Cross sections
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) {
				for (p=0; p<Ng[etaL+1]; p++) {
					phiS=0.0;
					sigmaT[etaL+1][p][i][j]  =0.0;
					nuSigmaF[etaL+1][p][i][j]=0.0;
					chi[etaL+1][p][i][j]     =0.0;
					s_ext[etaL+1][p][i][j]   =0.0;
					for (g=omegaP[etaL+1][p]; g<omegaP[etaL+1][p+1]; g++) {
						phiS+=phi[etaL][g][i][j];
						sigmaT[etaL+1][p][i][j]  +=phiG[g][i][j]*sigmaT[etaL][g][i][j];
						nuSigmaF[etaL+1][p][i][j]+=phiG[g][i][j]*nuSigmaF[etaL][g][i][j];
						chi[etaL+1][p][i][j]     +=chi[etaL][g][i][j];
						s_ext[etaL+1][p][i][j]   +=s_ext[etaL][g][i][j];
					}
					sigmaT[etaL+1][p][i][j]  /=phiS;
					nuSigmaF[etaL+1][p][i][j]/=phiS;
				}
				for (pp=0; pp<Ng[etaL+1]; pp++) {
					for (p=0; p<Ng[etaL+1]; p++) {
						phiS=0.0;
						sigmaS[etaL+1][pp][p][i][j]=0.0;
						for (gg=omegaP[etaL+1][pp]; gg<omegaP[etaL+1][pp+1]; gg++) {
							Ssum=0.0;
							for (g=omegaP[etaL+1][p]; g<omegaP[etaL+1][p+1]; g++) Ssum+=sigmaS[etaL][gg][g][i][j];
							sigmaS[etaL+1][pp][p][i][j]+=Ssum*phiG[gg][i][j];
							phiS+=phiG[gg][i][j];
						}
						//cout<<"pp "<<pp<<" p "<<p<<" "<<sigmaS[etaL+1][pp][p][i][j]<<" "<<phiS<<endl;
						sigmaS[etaL+1][pp][p][i][j]/=phiS;
						//cout<<"pp "<<pp<<" p "<<p<<" "<<sigmaS[etaL+1][pp][p][i][j]<<endl;
					}
				}
			}
		}
		
		// X Grid Diffusion consistency terms
		for (i=1; i<nx; i++) {
			for (j=0; j<ny; j++) {
				for (p=0; p<Ng[etaL+1]; p++) {
					phiP=0.0; phiN=0.0;
					D_xP[etaL+1][p][i][j]=0.0;
					D_xN[etaL+1][p][i][j]=0.0;
					for (g=omegaP[etaL+1][p]; g<omegaP[etaL+1][p+1]; g++) {
						phiP+=phi[etaL][g][i][j];
						phiN+=phi[etaL][g][i-1][j];
						D_xP[etaL+1][p][i][j]+=D_xP[etaL][g][i][j]*phiG[g][i][j];
						D_xN[etaL+1][p][i][j]+=D_xN[etaL][g][i][j]*phiG[g][i-1][j];
					}
					D_xP[etaL+1][p][i][j]/=phiP;
					D_xN[etaL+1][p][i][j]/=phiN;
				}
			}
		}
		// Y Grid Diffusion consistency terms
		for (i=0; i<nx; i++) {
			for (j=1; j<ny; j++) {
				for (p=0; p<Ng[etaL+1]; p++) {
					phiP=0.0; phiN=0.0;
					D_yP[etaL+1][p][i][j]=0.0;
					D_yN[etaL+1][p][i][j]=0.0;
					for (g=omegaP[etaL+1][p]; g<omegaP[etaL+1][p+1]; g++) {
						phiP+=phi[etaL][g][i][j];
						phiN+=phi[etaL][g][i][j-1];
						D_yP[etaL+1][p][i][j]+=D_yP[etaL][g][i][j]*phiG[g][i][j];
						D_yN[etaL+1][p][i][j]+=D_yN[etaL][g][i][j]*phiG[g][i][j-1];
					}
					D_yP[etaL+1][p][i][j]/=phiP;
					D_yN[etaL+1][p][i][j]/=phiN;
				}
			}
		}
		// Left and Right Boundary Conditions
		for (j=0; j<ny+1; j++) {
			for (p=0; p<Ng[etaL+1]; p++) {
				phiP=0.0; phiN=0.0;
				FL[etaL+1][p][j]=0.0;
				FR[etaL+1][p][j]=0.0;
				for (g=omegaP[etaL+1][p]; g<omegaP[etaL+1][p+1]; g++) {
					phiN+=phiB_L[etaL][g][j];
					phiP+=phiB_R[etaL][g][j];
					FL[etaL+1][p][j]+=FL[etaL][g][j]*phiB_L[etaL][g][j];
					FR[etaL+1][p][j]+=FR[etaL][g][j]*phiB_R[etaL][g][j];
				}
				FL[etaL+1][p][j]/=phiN;
				FR[etaL+1][p][j]/=phiP;
				
				phiP=0.0; phiN=0.0;
				D_xP[etaL+1][p][0][j]=0.0;
				D_xN[etaL+1][p][0][j]=0.0;
				for (g=omegaP[etaL+1][p]; g<omegaP[etaL+1][p+1]; g++) {
					phiP+=phiG[g][0][j];
					phiN+=phiB_L[etaL][g][j];
					D_xP[etaL+1][p][0][j]+=D_xP[etaL][g][0][j]*phiG[g][0][j];
					D_xN[etaL+1][p][0][j]+=D_xN[etaL][g][0][j]*phiB_L[etaL][g][j];
				}
				D_xP[etaL+1][p][0][j]/=phiP;
				D_xN[etaL+1][p][0][j]/=phiN;
				
				phiP=0.0; phiN=0.0;
				D_xP[etaL+1][p][nx][j]=0.0;
				D_xN[etaL+1][p][nx][j]=0.0;
				for (g=omegaP[etaL+1][p]; g<omegaP[etaL+1][p+1]; g++) {
					phiP+=phiB_R[etaL][g][j];
					phiN+=phiG[g][nx-1][j];
					D_xP[etaL+1][p][nx][j]+=D_xP[etaL][g][nx][j]*phiB_R[etaL][g][j];
					D_xN[etaL+1][p][nx][j]+=D_xN[etaL][g][nx][j]*phiG[g][nx-1][j];
				}
				D_xP[etaL+1][p][nx][j]/=phiP;
				D_xN[etaL+1][p][nx][j]/=phiN;
			}
		}
		// Bottom and Top Boundary Conditions
		for (i=0; i<nx; i++) {
			for (p=0; p<Ng[etaL+1]; p++) {
				phiP=0.0; phiN=0.0;
				FB[etaL+1][p][i]=0.0;
				FT[etaL+1][p][i]=0.0;
				for (g=omegaP[etaL+1][p]; g<omegaP[etaL+1][p+1]; g++) {
					phiN+=phiB_B[etaL][g][i];
					phiP+=phiB_T[etaL][g][i];
					FB[etaL+1][p][i]+=FB[etaL][g][i]*phiB_B[etaL][g][i];
					FT[etaL+1][p][i]+=FT[etaL][g][i]*phiB_T[etaL][g][i];
				}
				FB[etaL+1][p][i]/=phiN;
				FT[etaL+1][p][i]/=phiP;
				
				phiP=0.0; phiN=0.0;
				D_yP[etaL+1][p][i][0]=0.0;
				D_yN[etaL+1][p][i][0]=0.0;
				for (g=omegaP[etaL+1][p]; g<omegaP[etaL+1][p+1]; g++) {
					phiP+=phiG[g][i][0];
					phiN+=phiB_B[etaL][g][i];
					D_yP[etaL+1][p][i][0]+=D_yP[etaL][g][i][0]*phiG[g][i][0];
					D_yN[etaL+1][p][i][0]+=D_yN[etaL][g][i][0]*phiB_B[etaL][g][i];
				}
				D_yP[etaL+1][p][i][0]/=phiP;
				D_yN[etaL+1][p][i][0]/=phiN;
				
				phiP=0.0; phiN=0.0;
				D_yP[etaL+1][p][i][ny]=0.0;
				D_yN[etaL+1][p][i][ny]=0.0;
				for (g=omegaP[etaL+1][p]; g<omegaP[etaL+1][p+1]; g++) {
					phiP+=phiB_T[etaL][g][i];
					phiN+=phiG[g][i][ny-1];
					D_yP[etaL+1][p][i][ny]+=D_yP[etaL][g][i][ny]*phiB_T[etaL][g][i];
					D_yN[etaL+1][p][i][ny]+=D_yN[etaL][g][i][ny]*phiG[g][i][ny-1];
				}
				D_yP[etaL+1][p][i][ny]/=phiP;
				D_yN[etaL+1][p][i][ny]/=phiN;
			}
		}
		
		// Go to the next coarser grid
		if (etaL==eta_star-2) {
			if ( KE_problem ) oneGroupSolKE();
			else oneGroupSolFS();
		}
		else coarseGridSol(etaL+1);
		
		//temfile<<"D_xN 3\n";
		//write_grid_edge_x_dat(D_xN, 16, temfile);
		
		// Calculate itteration data
		rho_p=norm_p/norm_pL;
		norm_k=abs(k_eff-k_effL);
		rho_k=norm_k/norm_kL;
		norm_phi[etaL].push_back(norm_p);
		rho_phi[etaL].push_back(rho_p);
		norm_keff[etaL].push_back(norm_k);
		rho_keff[etaL].push_back(rho_k);
		
		l++; // increment for next loop
		num_grid[etaL].back()++;
	} while ( ( norm_p>epsilon_phi[etaL]*(1/rho_p-1) or norm_k>epsilon_keff[etaL]*(1/rho_k-1) ) and l<stop_phi[etaL] );
	n_in=l;
	num_losi[etaL].push_back(n_in);
	
}
//======================================================================================//

//======================================================================================//
//++ Calculate the Diffusion Coefficients +++++++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void calcDiffusion() {
	int g, i, j;
	// Calculate D^tilde and F factors
	
	// D_x
	for (g=0; g<Ng[0]; g++) {
		for (i=1; i<nx; i++) for (j=0; j<ny; j++) 
			D_xT[g][i][j]=2*(j_xT[g][i][j]+D_x[g][i][j]*(phiT[g][i][j]-phiT[g][i-1][j])/xe[i])/(phiT[g][i][j]+phiT[g][i-1][j]); // X Grid
		for ( j=0; j<ny; j++ ) { // Left and Right
			D_xT[g][0][j] =2*(j_xT[g][0][j] +D_x[g][0][j] *(phiT[g][0][j]   -phi_xT[g][0][j]) /xe[0]) /(phiT[g][0][j]   +phi_xT[g][0][j]); // Left Side
			FL[0][g][j]=j_xT[g][0][j]/phi_xT[g][0][j]; // Left F
			D_xT[g][nx][j]=2*(j_xT[g][nx][j]+D_x[g][nx][j]*(phi_xT[g][nx][j]-phiT[g][nx-1][j])/xe[nx])/(phi_xT[g][nx][j]+phiT[g][nx-1][j]); // Right Side
			FR[0][g][j]=j_xT[g][nx][j]/phi_xT[g][nx][j]; // Right F
		}
		// D_y
		for (i=0; i<nx; i++) for (j=1; j<ny; j++) 
			D_yT[g][i][j]=2*(j_yT[g][i][j]+D_y[g][i][j]*(phiT[g][i][j]-phiT[g][i][j-1])/ye[j])/(phiT[g][i][j]+phiT[g][i][j-1]); // Y Grid
		for ( i=0; i<nx; i++ ) { // Bottom and Top
			D_yT[g][i][0] =2*(j_yT[g][i][0] +D_y[g][i][0] *(phiT[g][i][0]   -phi_yT[g][i][0]) /ye[0]) /(phiT[g][i][0]   +phi_yT[g][i][0]); // Bottom Side
			FB[0][g][i]=j_yT[g][i][0]/phi_yT[g][i][0]; // Bottom F
			D_yT[g][i][ny]=2*(j_yT[g][i][ny]+D_y[g][i][ny]*(phi_yT[g][i][ny]-phiT[g][i][ny-1])/ye[ny])/(phi_yT[g][i][ny]+phiT[g][i][ny-1]); // Top
			FT[0][g][i]=j_yT[g][i][ny]/phi_yT[g][i][ny]; // Top F
		}
		for (i=0; i<nx+1; i++) for (j=0; j<ny; j++) D_xP[0][g][i][j]=D_x[g][i][j]-0.5*D_xT[g][i][j]*xe[i]; // D_xP X Grid
		for (i=0; i<nx; i++) for (j=0; j<ny+1; j++) D_yP[0][g][i][j]=D_y[g][i][j]-0.5*D_yT[g][i][j]*ye[j]; // D_yP Y Grid
		for (i=0; i<nx+1; i++) for (j=0; j<ny; j++) D_xN[0][g][i][j]=D_x[g][i][j]+0.5*D_xT[g][i][j]*xe[i]; // D_xN X Grid
		for (i=0; i<nx; i++) for (j=0; j<ny+1; j++) D_yN[0][g][i][j]=D_y[g][i][j]+0.5*D_yT[g][i][j]*ye[j]; // D_yN Y Grid
		
		// In case of reflective BC set edge values to zero
		switch (kbc) {
			case 3: // BC type 3 Reflective on all sides 
				for (i=0; i<nx; i++) { FB[0][g][i]=0;	FT[0][g][i]=0; } // Bottom and Top
				for (j=0; j<ny; j++) { FL[0][g][j]=0;	FR[0][g][j]=0; } // Left and Right
				break;
			case 4: // BC type 4 Reflective on Bottom and Left
				for (i=0; i<nx; i++) FB[0][g][i]=0;
				for (j=0; j<ny; j++) FL[0][g][j]=0;
				break;
			case 5: // BC type 5 Reflective on Top and Right
				for (i=0; i<nx; i++) FT[0][g][i]=0;
				for (j=0; j<ny; j++) FR[0][g][j]=0;
				break;
			default:
				break;
		}
	}
}
//======================================================================================//

//======================================================================================//
//++ Set initial guess for LO solutions ++++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void initialLoGuess() {
	int i, j, g, gg, k;
	// Low-order fluxes
	// Fine group
	for (g=0; g<Ng[0]; g++) {
		for (i=0; i<nx; i++) for (j=0; j<ny; j++) phi[0][g][i][j]=phiT[g][i][j];
		for (i=0; i<nx; i++) {
			phiB_B[0][g][i]=phi_yT[g][i][0];
			phiB_T[0][g][i]=phi_yT[g][i][ny];
		}
		for (j=0; j<ny; j++) {
			phiB_L[0][g][j]=phi_xT[g][0][j];
			phiB_R[0][g][j]=phi_xT[g][nx][j];
		}
	}
	// Coarse Group flux
	for (k=1; k<eta_star; k++) {
		for (g=0; g<Ng[k]; g++) {
			for (i=0; i<nx; i++) {
				for (j=0; j<ny; j++) {
					phi[k][g][i][j]=0;
					for (gg=omegaP[k][g]; gg<omegaP[k][g+1]; gg++) phi[k][g][i][j]+=phi[k-1][gg][i][j];
				}
			}
			for (i=0; i<nx; i++) {
				phiB_B[k][g][i]=0;
				phiB_T[k][g][i]=0;
				for (gg=omegaP[k][g]; gg<omegaP[k][g+1]; gg++) {
					phiB_B[k][g][i]=phiB_B[k-1][gg][i];
					phiB_T[k][g][i]=phiB_T[k-1][gg][i];
				}
			}
			for (j=0; j<ny; j++) {
				phiB_L[k][g][j]=0;
				phiB_R[k][g][j]=0;
				for (gg=omegaP[k][g]; gg<omegaP[k][g+1]; gg++) {
					phiB_L[k][g][j]=phiB_L[k-1][g][j];
					phiB_R[k][g][j]=phiB_R[k-1][g][j];
				}
			}
		}
	}
}
//======================================================================================//

//======================================================================================//
//++ Iterate to converge on solution +++++++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void Iterations() {
	int g, gg, i, j, m, s, k;
	double norm_si, norm_siL, phi_sum, res, psi_const=1.0;
	double k_effL, norm_k, norm_kL, rho_k;
	double S_f, S_s;
	
	norm_si=4.0*pi;
	norm_siL=norm_si;
	norm_k=1;
	norm_kL=norm_k;
	s=0;
	rho_ho.push_back(0.5);
	rho_kho.push_back(0.5);
	
	// begin iterations
	while ( norm_si>epsilon_si*(1/rho_ho[s]-1) or norm_k>epsilon_kho*(1/rho_kho[s]-1) or s==0 ) { //========================================================================
		cout<<" -- Iteration # "<<s<<" -- "<<endl;
		temfile<<" -- Iteration # "<<s<<" -- "<<endl;
		
		// set previous iteration to phiLast
		for (g=0; g<Ng[0]; g++) for (i=0; i<nx; i++) for (j=0; j<ny; j++) phiL[g][i][j]=phiT[g][i][j];
		norm_siL=norm_si; // set previous norm to normLast
		norm_kL=norm_k;
		k_effL=k_eff;
		
		//temfile<<"phiT 1\n";
		//write_group_average_dat(phiT, 0, 16, temfile);
		t_ho=clock(); // start high-order timer
		if ( s==0 ) {
			initialAngleSweep(); // calculate initial transport solution from constant angular flux <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			initialLoGuess();    // creat initial guess for LO solutions to prevent devision by zero <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
		}
		else {
			calcCorrection(0);   // Correct Flux using correction factors <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			angleSweep();       // perform sweep through angles <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		}
		t_ho=clock()-t_ho; // stop high-order timer
		dt_ho.push_back(((double)t_ho)/CLOCKS_PER_SEC); // add high-order solution time to vector
		
		//temfile<<"phiT 2\n";
		//write_group_average_dat(phiT, 0, 16, temfile);
		
		//write_grid_average_dat(phi, 16, temfile);
		
		// Calculate D^tilde and F factors
		calcDiffusion(); // Calculate diffusion factors and Boundary conditions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		
		// Calculate Transport Residuals
		res_ho=0;
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) {
				S_f=0.0;
				for (g=0; g<Ng[0]; g++) S_f+=nuSigmaF[0][g][i][j]*phi[0][g][i][j];
				for (g=0; g<Ng[0]; g++) {
					S_s=0.0;
					for (gg=0; gg<Ng[0]; gg++) S_s+=sigmaS[0][gg][g][i][j]*phi[0][gg][i][j];
					res=abs((j_xT[g][i+1][j]-j_xT[g][i][j])*hy[j]+(j_yT[g][i][j+1]-j_yT[g][i][j])*hx[i]+
						(sigmaT[0][g][i][j]*phiT[g][i][j]-S_s-chi[0][g][i][j]*S_f/k_eff-s_ext[0][g][i][j])*hx[i]*hy[j]);
					if ( res>res_ho ) { res_ho=res; i_ho=i; j_ho=j; }
				}
			}
		}
		num_sol=0;
		for (k=0; k<eta_star; k++) num_grid[k].push_back(0);
		// Solve NDA on each grid
		t_lo=clock(); // start low-order timer
		coarseGridSol(0); // call MLNDA function for fine group <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		t_lo=clock()-t_lo; // stop low-order timer
		dt_lo[0].push_back(((double)t_lo)/CLOCKS_PER_SEC); // add high-order solution time to vector
		num_mtot.push_back(num_sol); // Total # of times LO matirx was solved
		
		
		cout<<"Iteration Completed in "<<dt_ho[s]+dt_lo[0][s]<<" sec\n";
		
		//write_grid_average_dat(phi, 16, temfile);
		
		// Find new norm
		norm_si=0;
		for (g=0; g<Ng[0]; g++) 
			for (i=0; i<nx; i++) 
				for (j=0; j<ny; j++) if ( abs(phiT[g][i][j]-phiL[g][i][j])>norm_si ) norm_si=abs(phiT[g][i][j]-phiL[g][i][j]);
		norm_ho.push_back(norm_si);
		rho_ho.push_back(norm_si/norm_siL);
		
		norm_k=abs(k_eff-k_effL);
		rho_k=norm_k/norm_kL;
		norm_kho.push_back(norm_k);
		rho_kho.push_back(rho_k);
		
		s++;
	}
	n_iterations=s;
	
	// for periodic case normalize solution
	if ( kbc==3 or KE_problem ) {
		for (k=0; k<eta_star; k++) normalizeEigen(k);
	}
}
//======================================================================================//

//======================================================================================//
//++ Main program ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
int main (int argc, char* argv[])
{
	char fn[50], outf[50]={""};
	int i, j, m;
	double nu1, c;
	int run_stat_input, run_stat_quad, run_stat_initial;
	strcpy (fn,argv[1]);
	if ( argc!=2 ) cout << "usage: " << argv[0] << " input file name";
	else cout << fn << endl;
	
	run_stat_input=input(fn); // get input data
	
	if ( run_stat_input==0 ) {  // If there are no input errors Continue to run
		run_stat_quad=quadSet(); // find quadrature 
		
		if ( run_stat_quad==0 ) { // If there are no quadrature errors continue to run
			run_stat_initial=initialize(); // initialize memory space for solution
			
			if ( run_stat_initial==0 ) { // if there are no initializtion errors continue to run
				strcat (outf,fn);
				strcat (outf,".temp.csv");
				cout<<outf<<endl;
				temfile.open(outf); // open temporary file
				cout<<"++++++++++++++++++++++\n";
				t = clock();    // start timer
				// ------------------------------------------
				Iterations(); // Call Iterations
				// ------------------------------------------
				t = clock() -t; // stop timer
				cout<<"++++++++++++++++++++++\n";
				temfile.close(); // close temporary file
				
				output(fn);     // write out solutions
			}
			else {
				outfile.close(); // Close output file
				cout<<"Initialization Error "<<run_stat_initial<<": Code failed from error in 'initialize' function\n";
			}
		}
		else {
			outfile.close(); // Close output file
			cout<<"Quadrature Error "<<run_stat_quad<<": Code failed from error in 'quadSet' function\n";
		}
	}
	else {
		outfile.close(); // Close output file
		cout<<"Input Error "<<run_stat_input<<": Code failed from error in 'input' function\n";
	}
	
	cout<<"Code Completed Successfully!!\n";
	
	return 0;
}
//======================================================================================//

