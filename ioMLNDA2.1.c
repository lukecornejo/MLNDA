#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <string>
#include <vector>

using namespace std;


double **sigma_gT, ***sigma_gS, **sigma_gF, **nu_gF, **chi_g, **s_gext; // material # arrays
int *matnum; // material #
vector<string> xsname; // material name
int **material; // material #

const double pi=3.141592654;

extern bool o_angular, KE_problem, write_xs; // option variables
extern int nx, ny, N;
extern double *x, *y;
extern ofstream outfile;
extern ofstream datfile;

extern double *mu, *eta, *xi, *w;
extern double *x, *y, *hx, *hy, *xe, *ye; // grid arrays
extern double **bcL, **bcR, **bcB, **bcT;
extern double epsilon_si, epsilon_kho, epsilon_lo, *epsilon_phi, *epsilon_keff; // default tolerances
extern double k_eff, delta;
extern vector<int> stop_phi, preconditioner_type; // limit number of iterations

extern int kbc;
extern int eta_star;
extern int *Ng, eta_star, **omegaP;
extern double ****phi , ****j_x , ****j_y; // NDA scalar flux and current solution
extern double ***phiT, ***phi_xT, ***phi_yT, ***j_xT, ***j_yT, ***phiL; // scalar flux and current from transport
extern double ***phiB_L, ***phiB_R, ***phiB_B, ***phiB_T; // edge scalar flux from NDA
extern double ****sigmaT, *****sigmaS, ****nuSigmaF, ****chi, ****s_ext; // material arrays
extern double ***FL, ***FR, ***FB, ***FT; // NDA Boundary conditions
extern double ***D, ***D_x, ***D_y; // Diffusion coefficients
extern double ***D_xT, ***D_yT; // Tilde
extern double ****D_xP, ****D_xN, ****D_yP, ****D_yN; // Tilde, Positive and Negative

// Variables for output
extern int n_iterations, sn;
// residual data
extern double res_ho;
extern int i_ho, j_ho;
extern double *res_bal;
extern int *i_bal, *j_bal, *g_bal;
extern double *res_mbal, *res_ml, *res_mr, *res_mb, *res_mt;
extern int *i_mbal, *j_mbal, *j_ml, *j_mr, *i_mb, *i_mt;
extern int *g_mbal, *g_ml, *g_mr, *g_mb, *g_mt;

// iteration data
extern vector<int> num_mtot;
extern vector<double> rho_ho, rho_kho, norm_ho, norm_kho, dt_ho;
extern vector< vector<double> > rho_phi, rho_keff, norm_phi, norm_keff, err_lo, dt_lo, dt_pc;
extern vector< vector<int> > num_losi, num_logm, num_grid; // # of and # of solver iterations 

extern clock_t t;
extern clock_t t_lo;
extern clock_t t_ho;
extern clock_t t_pc;

// Function to parse Integer input
int iparse(char str[], int& i, int nl) {
	int p, t, num=0;
	while ( isdigit(str[i+1]) ) {
		num+=str[i]-'0';
		num*=10;
		i++;
	}
	num+=str[i]-'0';
	return num;
}

// Function to parse Decimal input
double dparse(char str[], int& i, int nl) {
	double num=0.0;
	int p, t, power;
	while ( isdigit(str[i+1]) ) {
		t=str[i]-'0';
		num+=t;
		num*=10.0;
		i++;
	}
	t=str[i]-'0';
	num+=t;
	if ( str[i+1]=='.' ) {
		i++; p=1;
		while ( isdigit(str[i+1]) ) {
			i++;
			t=str[i]-'0';
			num+=t/pow(10.0,double(p));
			p++;
		}
	}
	if ( str[i+1]=='e' or str[i+1]=='E' ) {
		i++;
		if ( str[i+1]=='-' ) {
			i+=2;
			if ( isdigit(str[i]) ) power=iparse(str,i,nl);
			num/=pow(10.0,double(power));
		}
		else if ( str[i+1]=='+' ) {
			i+=2;
			if ( isdigit(str[i]) ) power=iparse(str,i,nl);
			num*=pow(10.0,double(power));
		}
		else if ( isdigit(str[i+1]) ) {
			i++;
			power=iparse(str,i,nl);
			num*=pow(10.0,double(power));
		}
	}
	return num;
}

string print_out(double value)
{
	string dest;
	char buffer[50];
	
    sprintf(buffer, "%.8e", value); // First print out using scientific notation with 0 mantissa digits
	
	dest=buffer;
	if ( value>=0.0 ) dest=" " + dest;
	dest=" " + dest;
	return dest;
}

string print_csv(double value)
{
	string dest;
	char buffer[50];
	
    sprintf(buffer, "%.8e", value); // First print out using scientific notation with 0 mantissa digits
	
	dest=buffer;
	if ( value>=0.0 ) dest=" " + dest;
	dest=dest+",";
	return dest;
}

// Function to find the material region
int matRegion(int *matnum, int *regnum, int i, int j, int nmat, int nreg, double xn[], double yn[], double xp[], double yp[]) {
	int k, m, reg, ret;
	double xc, yc;
	
	xc=(x[i+1]+x[i])/2.0;
	yc=(y[j+1]+y[j])/2.0;
	for (k=0; k<nreg; k++) {
		if ( xc>=xn[k] && xc<=xp[k] && yc>=yn[k] && yc<=yp[k] ) {
			reg=k;
			break;
		}
	}
	for (m=0; m<nmat; m++) {
		if ( regnum[reg]==matnum[m] ) {
			ret=m;
			break;
		}
	}
	return ret;
}
//======================================================================================//

//======================================================================================//
//++ Read Input ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
int XSinput(int k)
{
	string temp;
	char matf[52]={""};
	const int nl=250;
	char namel[nl], inl[nl];
	int g, p, i, G, m;
	
	temp=xsname[k];
	temp+=".xs"; // XS input file name
	
	for (i=0; i<temp.size(); i++) matf[i]=temp[i];
	// open file +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	ifstream xsfile;
	xsfile.open (matf); // open XS input file
	
	if ( !xsfile.is_open() ) { cout<<">> File "<<matf<<" could not be opened\n"; return 1; }
	else {
		// read in data
		xsfile.getline(namel, nl); // read in name of the input
		
		xsfile.getline(inl, nl); // get # of groups
		for (i=0; i<nl; i++) {
			if ( isdigit(inl[i]) ) {
				G=iparse(inl,i,nl); // find number of mesh regions
			}
			if ( inl[i]==';' ) break;
		}
		if ( G!=Ng[0] ) { cout<<">> Material "<<matf<<" has mismatched # of Groups\n"; return 2; }
		
		xsfile.getline(inl, nl); // space
		
		for (g=0; g<G; g++) {
			xsfile.getline(inl, nl); // Get line of XS
			p=0;
			for (i=0; i<nl; i++) {
				if ( isdigit(inl[i]) ) {
					if (p==0)      sigma_gT[k][g]=dparse(inl,i,nl);
					else if (p==1) sigma_gF[k][g]=dparse(inl,i,nl);
					else if (p==2)    nu_gF[k][g]=dparse(inl,i,nl);
					else if (p==3)    chi_g[k][g]=dparse(inl,i,nl);
					else break;
					p++;
				}
				if ( inl[i]==';' ) break;
			}
		}
		
		xsfile.getline(inl, nl); // space
		
		for (g=0; g<G; g++) {
			xsfile.getline(inl, nl); // Get line of XS
			p=0;
			for (i=0; i<nl; i++) {
				if ( isdigit(inl[i]) ) {
					sigma_gS[k][g][p]=dparse(inl,i,nl);
					p++;
				}
				if ( !(p<G) ) break;
				if ( inl[i]==';' ) break;
			}
		}
		xsfile.close(); // close XS input file ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	}
	return 0;
}
//======================================================================================//

//======================================================================================//
//++ Read Input ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
int input(char fn[]) //
{
	char inpf[25]={""}, outf[25]={""}, temp[25]={""};
	const int nl=250;
	char namel[nl], inl[nl];
	int *regnum, sum;
	double *xn, *yn, *xp, *yp, *st, *ss, *sf, *nf, *se, *xg, *yg;
	double *xbc_e, *ybc_e, *bc_l, *bc_r, *bc_b, *bc_t, center; // bc temp
	double dtemp;
	int i, j, g, gg, k, p, nmat, nreg, ngx=0, ngy=0,  *nxt, *nyt, nbc, xbc_n, ybc_n;
	
	// input file name
	strcat (inpf,fn); strcat (inpf,".inp"); // input file name
	cout << inpf << endl;
	
	// open file +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	ifstream infile;
	infile.open (inpf); // open input file
	
	// read in data
	infile.getline(namel, nl); // read in name of the input
	
	// read x grid data
	infile.getline(inl, nl); // x grid edges
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			if ( ngx==0 ) {
				ngx=iparse(inl,i,nl); // find number of mesh regions
				nxt=new int[ngx];
				xg=new double[ngx+1];
				p=0;
			}
			else {
				xg[p]=dparse(inl,i,nl); p++;
			}
		}
		if ( inl[i]==';' ) break;
	}
	if ( p!=ngx+1 )	{ cout<<">> x grid input error a\n"; return 1; }
	
	infile.getline(inl, nl); // # cells in x grid
	p=0;
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			nxt[p]=iparse(inl,i,nl); p++;
		}
		if ( inl[i]==';' )	break;
	}
	if ( p!=ngx ) { cout<<">> x grid input error b\n"; return 1; }
	
	// read y grid data
	infile.getline(inl, nl); // read in grid data
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			if ( ngy==0 ) {
				ngy=iparse(inl,i,nl);
				nyt=new int[ngy];
				yg=new double[ngy+1];
				p=0;
			}
			else {
				yg[p]=dparse(inl,i,nl); p++; // y grid zones
			}
		}
		if ( inl[i]==';' )	break;
	}
	if ( p!=ngy+1 )	{ cout<<">> y grid input error a\n"; return 2; }
	
	infile.getline(inl, nl); // read in grid data
	p=0;
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			nyt[p]=iparse(inl,i,nl); p++; // # cells in y grid zone
		}
		if ( inl[i]==';' ) break;
	}
	if ( p!=ngy ) { cout<<">> y grid input error b\n"; return 2; }
	
	
	// read in boundary conditions
	p=0; // Get Type of BC
	infile.getline(inl, nl);
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			if ( p==0 ) kbc=iparse(inl,i,nl); // kind of BC
			else       KE_problem=iparse(inl,i,nl); // Type of problem
			p++;
		}
		if ( inl[i]==';' ) break;
	}
	
	p=0; // Get Bottom and Top BC
	infile.getline(inl, nl);
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			if ( p==0 ) {
				nbc=iparse(inl,i,nl); // # of BC regions
				xbc_n=nbc; 
				xbc_e=new double[nbc+1];
				bc_b=new double[nbc];
				bc_t=new double[nbc];
			}
			else if ( p<nbc ) 	xbc_e[p]=dparse(inl,i,nl); // Left   BC or Quadrant 1
			else if ( p<2*nbc ) bc_b[p-nbc]  =dparse(inl,i,nl); // Bottom BC or Quadrant 2
			else                bc_t[p-2*nbc]=dparse(inl,i,nl); // Right  BC or Quadrant 3
			p++;
		}
		if ( inl[i]==';' ) break;
	}
	
	p=0; // Get Left and Right BC
	infile.getline(inl, nl);
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			if ( p==0 ) {
				nbc=iparse(inl,i,nl); // # of BC regions
				ybc_n=nbc;
				ybc_e=new double[nbc+1];
				bc_l=new double[nbc];
				bc_r=new double[nbc];
			}
			else if ( p<nbc )   ybc_e[p]=dparse(inl,i,nl); // Left   BC or Quadrant 1
			else if ( p<2*nbc ) bc_l[p-nbc]  =dparse(inl,i,nl); // Bottom BC or Quadrant 2
			else                bc_r[p-2*nbc]=dparse(inl,i,nl); // Right  BC or Quadrant 3
			p++;
		}
		if ( inl[i]==';' ) break;
	}
	
	infile.getline(inl, nl); // get number of energy grids
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			eta_star=iparse(inl,i,nl); break;
		}
		if ( inl[i]==';' ) break;
	}
	
	Ng=new int[eta_star]; 
	
	infile.getline(inl, nl); // get number of energy groups for fine mesh
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			Ng[0]=iparse(inl,i,nl); break;
		}
		if ( inl[i]==';' ) break;
	}
	
	omegaP=new int*[eta_star];
	omegaP[0]=new int[Ng[0]+1];
	
	omegaP[0][0]=0;
	for (g=0; g<Ng[0]; g++) omegaP[0][g+1]=omegaP[0][g]+1;
	
	for (k=1; k<eta_star-1; k++) {
		infile.getline(inl, nl); // get coarse grids
		p=0;
		for (i=0; i<nl; i++) {
			if ( isdigit(inl[i]) ) {
				if ( p==0 ) {
					Ng[k]=iparse(inl,i,nl); // 
					omegaP[k]=new int[Ng[k]+1];
					omegaP[k][p]==0;
				}
				else omegaP[k][p]=omegaP[k][p-1]+iparse(inl,i,nl);
				p++;
			}
			if ( inl[i]==';' ) break;
		}
	}
	
	
	Ng[eta_star-1]=1;
	omegaP[eta_star-1]=new int[2];
	omegaP[eta_star-1][0]=0;
	omegaP[eta_star-1][1]=Ng[eta_star-2];
	
	for (k=1; k<eta_star; k++) { if (omegaP[k][Ng[k]]!=Ng[k-1]) { cout<<">>Error in energy groups: Grid "<<k<<endl; return 3; } }
	
	infile.getline(inl, nl); // space
	
	infile.getline(inl, nl);  // number of materials
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			nmat=iparse(inl,i,nl); break;
		}
	}
	
	// cross sections
	matnum  =new int[nmat];
	sigma_gT=new double*[nmat];
	sigma_gS=new double**[nmat];
	sigma_gF=new double*[nmat];
	nu_gF   =new double*[nmat];
	chi_g   =new double*[nmat];
	s_gext  =new double*[nmat];
	for (k=0; k<nmat; k++) {
		sigma_gT[k]=new double[Ng[0]];
		sigma_gS[k]=new double*[Ng[0]];
		sigma_gF[k]=new double[Ng[0]];
		nu_gF[k]   =new double[Ng[0]];
		chi_g[k]   =new double[Ng[0]];
		s_gext[k]  =new double[Ng[0]];
		for (g=0; g<Ng[0]; g++) {
			s_gext[k][g]=0.0;
			sigma_gS[k][g]=new double[Ng[0]];
		}
	}
	
	for (k=0; k<nmat; k++) {
		p=0;
		infile.getline(inl, nl);
		for (i=0; i<nl; i++) {
			if ( isdigit(inl[i]) ) {
				if ( p==0 )  matnum[k]=iparse(inl,i,nl); // Material #
				else if ( (p%2) != 0 ) g=iparse(inl,i,nl); // group of external source
				else s_gext[k][g-1]=4.0*pi*dparse(inl,i,nl); // external source
				p++;
			}
			if ( inl[i]=='"' ) { // find XS file name
				i++;
				for (j=0; j<25; j++) {
					if (inl[i]=='"') break;
					temp[j]=inl[i];
					i++;
				}
				i++;
			}
			if ( inl[i]==';' ) break;
		}
		xsname.push_back(temp);
		for (i=0; i<25; i++) temp[i]='\0';
	}
	
	
	infile.getline(inl, nl);  // number of material regions
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			nreg=iparse(inl,i,nl); break;
		}
	}
	
	// material edges
	regnum=new int[nreg];
	xn=new double[nreg];  xp=new double[nreg];
	yn=new double[nreg];  yp=new double[nreg];
	
	for (k=0; k<nreg; k++) {
		p=0;
		infile.getline(inl, nl);
		for (i=0; i<nl; i++) {
			if ( isdigit(inl[i]) ) {
				if ( p==0 )  regnum[k]=iparse(inl,i,nl); // Material Number in Region (int)
				else if ( p==1 ) xn[k]=dparse(inl,i,nl); // Left   Boundary of Region (double)
				else if ( p==2 ) yn[k]=dparse(inl,i,nl); // Bottom Boundary of Region (double)
				else if ( p==3 ) xp[k]=dparse(inl,i,nl); // Right  Boundary of Region (double)
				else if ( p==4 ) yp[k]=dparse(inl,i,nl); // Top    Boundary of Region (double)
				else { cout<<">> Region Input Error \n"; return 4; }
				p++;
			}
			if ( inl[i]==';' ) break;
		}
	}
	
	p=0;
	infile.getline(inl, nl); // Read High order Data
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			if ( p==0 )      N=iparse(inl,i,nl); // Quadrature # (int)
			if ( p==1 )      epsilon_si=dparse(inl,i,nl); // get flux convergence criteria (double)
			if ( p==2 )      epsilon_kho=dparse(inl,i,nl); // High Order k_eff convergence criteria (double)
			p++;
		}
		if ( inl[i]==';' ) break;
	}
	// Input Error Handling
	if ( epsilon_si<1e-12 ) { cout<<">> Error in SI convergence criteria input!\n"; epsilon_si=1e-5; }
	if ( epsilon_lo<1e-14 or epsilon_lo>epsilon_si) { cout<<">> Error in low-order convergence criteria input!\n"; epsilon_lo=1e-10; }
	epsilon_kho=epsilon_si; // Default
	if ( N!=4 and N!=6 and N!=8 and N!=12 and N!=16 and N!=20 and N!=36 ) {
		cout<<">> Error in Quadrature input!\n"; N=36;
	}
	
	// Initialize Low-order convergence criteria
	epsilon_phi=new double[eta_star];
	epsilon_keff=new double[eta_star];
	stop_phi.resize(eta_star);
	preconditioner_type.resize(eta_star);
	for (k=0; k<eta_star; k++) stop_phi[k]=10000;
	epsilon_phi[0]=epsilon_si/10.0;
	for (k=1; k<eta_star; k++) epsilon_phi[k]=epsilon_phi[k-1]/10.0;
	for (k=0; k<eta_star; k++) epsilon_keff[k]=epsilon_phi[k]; // default
	for (k=0; k<eta_star; k++) preconditioner_type[k]=3;
	
	p=0; // Read LO Flux NDA stopping criteria
	infile.getline(inl, nl); // Read LO Flux NDA stopping criteria
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			if ( p==0 ) if( !iparse(inl,i,nl) ) break; // Read whether or not to get data
			if ( p>0 and p<eta_star+1 ) epsilon_phi[p-1]=dparse(inl,i,nl); // get low-order convergence criteria
			p++;
		}
		if ( inl[i]==';' ) break;
	}
	p=0;
	infile.getline(inl, nl); // Read LO Keff NDA stopping criteria
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			if ( p==0 ) if(iparse(inl,i,nl)==0) break; // Read whether or not to get data
			if ( p>0 and p<eta_star+1 ) epsilon_keff[p-1]=dparse(inl,i,nl); // get low-order convergence criteria 
			p++;
		}
		if ( inl[i]==';' ) break;
	}
	p=0;
	infile.getline(inl, nl); // Read Iteration # NDA stopping criteria
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			if ( p==0 ) if(iparse(inl,i,nl)==0) break; // Read whether or not to get data
			if ( p>0 and p<eta_star+1 ) stop_phi[p-1]=iparse(inl,i,nl); // get low-order convergence criteria
			p++;
		}
		if ( inl[i]==';' ) break;
	}
	
	p=0; // Read preconditioner type
	infile.getline(inl, nl); // 
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			if ( p==0 ) if(iparse(inl,i,nl)==0) break; // Read whether or not to get data
			if ( p>0 and p<eta_star+1 ) preconditioner_type[p-1]=iparse(inl,i,nl); // Preconditioner type
			p++;
		}
		if ( inl[i]==';' ) break;
	}
	
	for (k=0; k<eta_star; k++) {
		if ( preconditioner_type[k]>3 and preconditioner_type[k]<2 ) {
			cout<<">>Error in preconditioner Type input in grid "<<k<<" use default preconditioner 3"<<endl;
			preconditioner_type[k]=3;
			return 5;
		}
	}
	
	p=0; // Read k-eigenvalue iteration data
	infile.getline(inl, nl); // 
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			if (p==0)      delta=dparse(inl,i,nl); // Get Delta for Weilandt-Shift (double)
			if (p==1) epsilon_lo=dparse(inl,i,nl); // Get convergence for BiCGStab solution (double)
			p++;
		}
		if ( inl[i]==';' ) break;
	}
	
	p=0; // Read Output options
	infile.getline(inl, nl); // 
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			if (p==0)      o_angular=iparse(inl,i,nl); // Write Angular Flux (bool) 
			if (p==1) write_xs=iparse(inl,i,nl); // Write Cross Sections (bool)
			p++;
		}
		if ( inl[i]==';' ) break;
	}
	
	infile.close(); // close input file ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	if ( KE_problem ) {
		for (k=0; k<nmat; k++) for (g=0; g<Ng[0]; g++) s_gext[k][g]=0; // Set source to zero if it is a k-eigenvalue problem
	}
	else { k_eff=1.0; delta=0.0; }
	
	
	// total grid
	nx=0; for (i=0; i<ngx; i++) nx+=nxt[i]; // find total number of cells in x grid
	ny=0; for (i=0; i<ngy; i++) ny+=nyt[i]; // find total number of cells in y grid
	
	x=new (nothrow) double[nx+1]; hx=new (nothrow) double[nx]; xe=new double[nx+1];
	y=new (nothrow) double[ny+1]; hy=new (nothrow) double[ny]; ye=new double[ny+1];
	
	x[0]=xg[0]; p=0;
	for (i=0; i<ngx; i++) {
		for (j=p; j<p+nxt[i]; j++) {
			hx[j]=(xg[i+1]-xg[i])/nxt[i];
			x[j+1]=x[j]+hx[j];
		}
		p+=nxt[i];
	}
	y[0]=yg[0]; p=0;
	for (i=0; i<ngy; i++) {
		for (j=p; j<p+nyt[i]; j++) {
		hy[j]=(yg[i+1]-yg[i])/nyt[i];
		y[j+1]=y[j]+hy[j];
		}
		p+=nyt[i];
	}
	xe[0]=0.5*hx[0]; xe[nx]=0.5*hx[nx-1];
	ye[0]=0.5*hy[0]; ye[ny]=0.5*hy[ny-1];
	for (i=1; i<nx; i++) xe[i]=0.5*(hx[i-1]+hx[i]);
	for (j=1; j<ny; j++) ye[j]=0.5*(hy[j-1]+hy[j]);
	
	// BC
	bcB=new double*[Ng[0]]; // Initialize Bottom BC
	bcT=new double*[Ng[0]]; // Initialize Top BC
	bcL=new double*[Ng[0]]; // Initialize Left BC
	bcR=new double*[Ng[0]]; // Initialize Right BC
	for (g=0; g<Ng[0]; g++) {
		bcB[g]=new double[nx]; // Initialize Bottom BC
		bcT[g]=new double[nx]; // Initialize Top BC
		bcL[g]=new double[ny]; // Initialize Left BC
		bcR[g]=new double[ny]; // Initialize Right BC
	}
	xbc_e[0]=x[0]; xbc_e[xbc_n]=x[nx];
	ybc_e[0]=y[0]; ybc_e[ybc_n]=y[ny];
	
	for (i=0; i<nx; i++) {
		center=(x[i]+x[i+1])/2;
		for (p=0; p<xbc_n; p++) {
			if ( xbc_e[p]<center and xbc_e[p+1]>center ) {
				for (g=0; g<Ng[0]; g++) {
					bcB[g][i]=bc_b[p];
					bcT[g][i]=bc_t[p];
				}
			}
		}
	}
	for (j=0; j<ny; j++) {
		center=(y[j]+y[j+1])/2;
		for (p=0; p<ybc_n; p++) {
			if ( ybc_e[p]<center and ybc_e[p+1]>center ) {
				for (g=0; g<Ng[0]; g++) {
					bcL[g][j]=bc_l[p];
					bcR[g][j]=bc_r[p];
				}
			}
		}
	}
	
	for (k=0; k<nmat; k++) {
		int xs_stat=XSinput(k); // Get Cross-Sections for each Material
		if ( xs_stat!=0 ) return xs_stat;
	}
	// Ensure there are no errors in cross sections
	for (k=0; k<nmat; k++) {
		dtemp=0.0; for (g=0; g<Ng[0]; g++) dtemp+=chi_g[k][g]; // Make sure chi sums to 1
		if ( dtemp!=0.0 ) for (g=0; g<Ng[0]; g++) chi_g[k][g]/=dtemp;
	}
	
	material=new int*[nx];
	for (i=0; i<nx; i++) material[i] =new int[ny];
	
	// cross section data
	sigmaT=new double***[eta_star];
	sigmaS=new double****[eta_star];
	nuSigmaF=new double***[eta_star];
	chi    =new double***[eta_star];
	s_ext =new double***[eta_star];
	for (k=0; k<eta_star; k++) {
		sigmaT[k]=new double**[Ng[k]];
		sigmaS[k]=new double***[Ng[k]];
		nuSigmaF[k]=new double**[Ng[k]];
		chi[k]   =new double**[Ng[k]];
		s_ext[k] =new double**[Ng[k]];
		for (g=0; g<Ng[k]; g++) {
			sigmaT[k][g]=new double*[nx];
			sigmaS[k][g]=new double**[Ng[k]];
			nuSigmaF[k][g]=new double*[nx];
			chi[k][g]   =new double*[nx];
			s_ext[k][g] =new double*[nx];
			for (gg=0; gg<Ng[k]; gg++) {
				sigmaS[k][g][gg]=new double*[nx];
				for (i=0; i<nx; i++) sigmaS[k][g][gg][i]=new double[ny];
			}
			for (i=0; i<nx; i++) {
				sigmaT[k][g][i]=new double[ny];
				nuSigmaF[k][g][i]=new double[ny];
				chi[k][g][i]   =new double[ny];
				s_ext[k][g][i] =new double[ny];
			}
		}
	}
	
	for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) {
			k=matRegion(matnum, regnum, i, j, nmat, nreg, xn, yn, xp, yp);
			material[i][j]=matnum[k];
			for (g=0; g<Ng[0]; g++) {
				sigmaT[0][g][i][j]  =sigma_gT[k][g];
				nuSigmaF[0][g][i][j]=nu_gF[k][g]*sigma_gF[k][g];
				chi[0][g][i][j]     =chi_g[k][g];
				s_ext[0][g][i][j]   = s_gext[k][g];
				for (gg=0; gg<Ng[0]; gg++) sigmaS[0][g][gg][i][j]=sigma_gS[k][g][gg];
				if ( sigmaT[0][g][i][j]==0.0 ) sigmaT[0][g][i][j]=1e-23;
			}
		}
	}
	
	
	strcat (outf,fn);
	strcat (outf,".out");
	cout<<outf<<endl;
	outfile.open(outf); // open output file. closed in output function
	
	outfile<<"Output of File : "<<outf<<endl;
	outfile<<"2D Multi-level NDA with Step Characteristics Transport\n";
	outfile<<"Version: 2.0\n";
	outfile<<"Programer : Luke Cornejo\n";
	outfile<<"Case Name : "<<namel<<endl;
	// current date/time based on current system
	time_t now = time(0);
	// convert now to string form
	char* dt = ctime(&now);
	
	outfile<<"Program Ran On: "<<dt<<endl;
	if ( KE_problem ) outfile<<" K-eigenvalue Problem \n";
	else outfile<<" Fixed Source Problem \n";
	
	outfile<<"+-----------------------------------------------+\n";
	outfile<<"Iteration Convergence Criteria: "<<epsilon_si<<endl;
	outfile<<"HO K_eff  Convergence Criteria: "<<epsilon_kho<<endl;
	outfile<<"BiCGStab Solver Tolerance: "<<epsilon_lo<<endl;
	outfile<<  "          Grid #                : ";
	for (k=0; k<eta_star; k++ ) outfile<<setw(10)<<k<<setw(6)<<"  ";
	outfile<<"\nLow-order  Convergence  Criteria: ";
	for (k=0; k<eta_star; k++ ) outfile<<print_out(epsilon_phi[k]);
	outfile<<"\nLO K_eff   Convergence  Criteria: ";
	for (k=0; k<eta_star; k++ ) outfile<<print_out(epsilon_keff[k]);
	outfile<<"\nMax Num. of Low-order Iterations: ";
	for (k=0; k<eta_star; k++ ) outfile<<setw(10)<<stop_phi[k]<<setw(6)<<"  ";
	outfile<<"\nLow-order    Preconditioner Type: ";
	for (k=0; k<eta_star; k++ ) outfile<<setw(10)<<preconditioner_type[k]<<setw(6)<<"  ";
	outfile<<endl;
	outfile<<"Weilandt-Shift Iteration Perameter Delta "<<print_out(delta)<<endl;
	outfile<<"+-----------------------------------------------+\n";
	
	outfile<<"\n -- Code Options -- \n";
	outfile<<"+-----------------------------------------------+\n";
	outfile<<"Write Angular Flux : "<<o_angular<<endl;
	outfile<<"Write Cross Sections : "<<write_xs<<endl;
	outfile<<"+-----------------------------------------------+\n";
	
	outfile<<"\n -- Material Properties --\n";
	for (k=0; k<nmat; k++) {
		outfile<<"Material # "<<matnum[k]<<"  Cross-section name: "<<xsname[k]<<endl;
		outfile<<"Group# |    Total XS   |   Fission XS  |      nuF      |      chi      | Ext Source \n";
		outfile.precision(6);
		for (g=0; g<Ng[0]; g++) outfile<<setw(6)<<g+1<<print_out(sigma_gT[k][g])<<print_out(sigma_gF[k][g])
			<<print_out(nu_gF[k][g])<<print_out(chi_g[k][g])<<print_out(s_gext[k][g])<<endl;
		outfile<<" Scattering Matrix \n";
		outfile<<"  g \\ g'  ";
		for (gg=0; gg<Ng[0]; gg++) outfile<<setw(6)<<gg+1<<setw(10)<<" ";
		outfile<<endl;
		for (g=0; g<Ng[0]; g++) {
			outfile<<setw(6)<<g+1;
			for (gg=0; gg<Ng[0]; gg++) outfile<<print_out(sigma_gS[k][g][gg]);
			outfile<<endl;
		}
	}
	
	return 0;
}
//======================================================================================//

//======================================================================================//
//++ function to write cell average values in viewer friendly format +++++++++++++++++++//
//======================================================================================//
void write_cell_average_out(double **outp, int outw, ofstream& file) { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	int i, j, m, g;
	// cell average values
	for (m=0; m<int(nx/10); m++) {
		file<<setw(outw)<<"sol. grid "<<setw(6)<<" ";
		for (i=m*10; i<(m+1)*10; i++) file<<print_out((x[i]+x[i+1])/2);
		file<<endl;
		file<<setw(outw)<<" "<<setw(6)<<"index";
		for (i=m*10; i<(m+1)*10; i++) file<<setw(outw-6)<<i+1<<setw(6)<<" ";
		file<<endl;
		for (j=ny-1; j>=0; j--) {
			file<<print_out((y[j]+y[j+1])/2)<<setw(5)<<j+1<<" ";
			for (i=m*10; i<(m+1)*10; i++) file<<print_out(outp[i][j]);
			file<<endl;
		}
	}
	if ( m*10<nx ) {
		file<<setw(outw)<<"sol. grid "<<setw(6)<<" ";
		for (i=m*10; i<nx; i++) file<<print_out((x[i]+x[i+1])/2);
		file<<endl;
		file<<setw(outw)<<" "<<setw(6)<<"index";
		for (i=m*10; i<nx; i++) file<<setw(outw-6)<<i+1<<setw(6)<<" ";
		file<<endl;
		for (j=ny-1; j>=0; j--) {
			file<<print_out((y[j]+y[j+1])/2)<<setw(5)<<j+1<<" ";
			for (i=m*10; i<nx; i++) file<<print_out(outp[i][j]);
			file<<endl;
		}
	}
}
//======================================================================================//
//++ function to write cell edge x values in viewer friendly format ++++++++++++++++++++//
//======================================================================================//
void write_cell_edge_x_out(double **outp, int outw, ofstream& file) { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	int i, j, m, g;
	// cell edge values on x grid
	for (m=0; m<int(nx/10); m++) {
		file<<setw(outw)<<"sol. grid "<<setw(6)<<" ";
		for (i=m*10; i<(m+1)*10; i++) file<<print_out(x[i]); file<<endl;
		file<<setw(outw)<<" "<<setw(6)<<"index";
		for (i=m*10; i<(m+1)*10; i++) file<<setw(outw-6)<<i+1<<setw(6)<<" "; file<<endl;
		for (j=ny-1; j>=0; j--) {
			file<<print_out((y[j]+y[j+1])/2)<<setw(5)<<j+1<<" ";
			for (i=m*10; i<(m+1)*10; i++) file<<print_out(outp[i][j]);
			file<<endl;
		}
	}
	if ( m*10<nx+1 ) {
		file<<setw(outw)<<"sol. grid "<<setw(6)<<" ";
		for (i=m*10; i<nx+1; i++) file<<print_out(x[i]); file<<endl;
		file<<setw(outw)<<" "<<setw(6)<<"index";
		for (i=m*10; i<nx+1; i++) file<<setw(outw-6)<<i+1<<setw(6)<<" "; file<<endl;
		for (j=ny-1; j>=0; j--) {
			file<<print_out((y[j]+y[j+1])/2)<<setw(5)<<j+1<<" ";
			for (i=m*10; i<nx+1; i++) file<<print_out(outp[i][j]); file<<endl;
		}
	}
}
//======================================================================================//
//++ function to write cell edge y values in viewer friendly format ++++++++++++++++++++//
//======================================================================================//
void write_cell_edge_y_out(double **outp, int outw, ofstream& file) { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	int i, j, m, g;
	// cell edge values on y grid
	for (m=0; m<int(nx/10); m++) {
		file<<setw(outw)<<"sol. grid "<<setw(6)<<" ";
		for (i=m*10; i<(m+1)*10; i++) file<<print_out((x[i]+x[i+1])/2); file<<endl;
		file<<setw(outw)<<" "<<setw(6)<<"index";
		for (i=m*10; i<(m+1)*10; i++) file<<setw(outw-6)<<i+1<<setw(6)<<" "; file<<endl;
		for (j=ny; j>=0; j--) {
			file<<print_out(y[j])<<setw(5)<<j+1<<" ";
			for (i=m*10; i<(m+1)*10; i++) file<<print_out(outp[i][j]); file<<endl;
		}
	}
	if ( m*10<nx ) {
		file<<setw(outw)<<"sol. grid "<<setw(6)<<" ";
		for (i=m*10; i<nx; i++) file<<print_out((x[i]+x[i+1])/2); file<<endl;
		file<<setw(outw)<<" "<<setw(6)<<"index";
		for (i=m*10; i<nx; i++) file<<setw(outw-6)<<i+1<<setw(6)<<" "; file<<endl;
		for (j=ny; j>=0; j--) {
			file<<print_out(y[j])<<setw(5)<<j+1<<" ";
			for (i=m*10; i<nx; i++) file<<print_out(outp[i][j]); file<<endl;
		}
	}
}
//======================================================================================//
//++ function to write cell average values in data format ++++++++++++++++++++++++++++++//
//======================================================================================//
void write_cell_average_dat(double **outp, int outw, ofstream& file) { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	int i, j, g;
	// cell average values
	file<<setw(outw)<<" sol. grid ,"<<setw(6)<<" "<<",";
	for (i=0; i<nx; i++) file<<print_csv((x[i]+x[i+1])/2);
	file<<endl;
	file<<setw(outw)<<","<<setw(6)<<"index"<<",";
	for (i=0; i<nx; i++) file<<setw(outw-6)<<i+1<<setw(6)<<",";
	file<<endl;
	for (j=0; j<ny; j++) {
		file<<print_csv((y[j]+y[j+1])/2)<<setw(5)<<j+1<<" ,";
		for (i=0; i<nx; i++) file<<print_csv(outp[i][j]);
		file<<endl;
	}
}
//======================================================================================//
//++ function to write cell average values in data format ++++++++++++++++++++++++++++++//
//======================================================================================//
void write_cell_average_dat(vector< vector<double> > outp, int outw, ofstream& file) { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	int i, j, g;
	// cell average values
	file<<setw(outw)<<" sol. grid ,"<<setw(6)<<" "<<",";
	for (i=0; i<nx; i++) file<<print_csv((x[i]+x[i+1])/2);
	file<<endl;
	file<<setw(outw)<<","<<setw(6)<<"index"<<",";
	for (i=0; i<nx; i++) file<<setw(outw-6)<<i+1<<setw(6)<<",";
	file<<endl;
	for (j=0; j<ny; j++) {
		file<<print_csv((y[j]+y[j+1])/2)<<setw(5)<<j+1<<" ,";
		for (i=0; i<nx; i++) file<<print_csv(outp[i][j]);
		file<<endl;
	}
}
//======================================================================================//
//++ function to write cell edge x values in data format +++++++++++++++++++++++++++++++//
//======================================================================================//
void write_cell_edge_x_dat(double **outp, int outw, ofstream& file) { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	int i, j, g;
	// cell edge values on x grid
	file<<setw(outw)<<" sol. grid ,"<<setw(6)<<" "<<",";
	for (i=0; i<nx+1; i++) file<<print_csv(x[i]);
	file<<endl;
	file<<setw(outw)<<","<<setw(6)<<"index"<<",";
	for (i=0; i<nx+1; i++) file<<setw(outw-6)<<i+1<<setw(6)<<",";
	file<<endl;
	for (j=0; j<ny; j++) {
		file<<print_csv((y[j]+y[j+1])/2)<<setw(5)<<j+1<<" ,";
		for (i=0; i<nx+1; i++) file<<print_csv(outp[i][j]);
		file<<endl;
	}
}
//======================================================================================//
//++ function to write cell edge y values in data format +++++++++++++++++++++++++++++++//
//======================================================================================//
void write_cell_edge_y_dat(double **outp, int outw, ofstream& file) { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	int i, j, g;
	
	// cell edge values on y grid
	file<<setw(outw)<<" sol. grid ,"<<setw(6)<<" "<<",";
	for (i=0; i<nx; i++) file<<print_csv((x[i]+x[i+1])/2);
	file<<endl;
	file<<setw(outw)<<","<<setw(6)<<"index"<<",";
	for (i=0; i<nx; i++) file<<setw(outw-6)<<i+1<<setw(6)<<",";
	file<<endl;
	for (j=0; j<ny+1; j++) {
		file<<print_csv(y[j])<<setw(5)<<j+1<<" ,";
		for (i=0; i<nx; i++) file<<print_csv(outp[i][j]);
		file<<endl;
	}
}
//======================================================================================//

//======================================================================================//
//++ function to write multi group cell edge avg values in output format +++++++++++++++//
//======================================================================================//
void write_group_average_out(double ***outp, int etaL, int outw, ofstream& file) {
	int g;
	for (g=0; g<Ng[etaL]; g++) {
		file<<"Grid # "<<etaL<<" Energy Group # "<<g+1<<endl;
		write_cell_average_out(outp[g], outw, file);
	}
}
//======================================================================================//
//++ function to write multi group cell edge x values in viewer friendly format ++++++++//
//======================================================================================//
void write_group_edge_x_out(double ***outp, int etaL, int outw, ofstream& file) { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	int g;
	for (g=0; g<Ng[etaL]; g++) {
		file<<"Grid # "<<etaL<<" Energy Group # "<<g+1<<endl;
		write_cell_edge_x_out(outp[g], outw, file);
	}
}
//======================================================================================//
//++ function to write multi group cell edge y values in viewer friendly format ++++++++//
//======================================================================================//
void write_group_edge_y_out(double ***outp, int etaL, int outw, ofstream& file) { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	int g;
	for (g=0; g<Ng[etaL]; g++) {
		file<<"Grid # "<<etaL<<" Energy Group # "<<g+1<<endl;
		write_cell_edge_y_out(outp[g], outw, file);
	}
}
//======================================================================================//
//++ function to write multi group cell average values in data format ++++++++++++++++++//
//======================================================================================//
void write_group_average_dat(double ***outp, int etaL, int outw, ofstream& file) { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	int g;
	for (g=0; g<Ng[etaL]; g++) {
		file<<"Grid # "<<etaL<<" Energy Group # "<<g+1<<endl;
		write_cell_average_dat(outp[g], outw, file);
	}
}
//======================================================================================//
//++ function to write multi group cell edge x values in data format +++++++++++++++++++//
//======================================================================================//
void write_group_edge_x_dat(double ***outp, int etaL, int outw, ofstream& file) { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	int g;
	for (g=0; g<Ng[etaL]; g++) {
		file<<"Grid # "<<etaL<<" Energy Group # "<<g+1<<endl;
		write_cell_edge_x_dat(outp[g], outw, file);
	}
}
//======================================================================================//
//++ function to write multi group cell edge y values in data format +++++++++++++++++++//
//======================================================================================//
void write_group_edge_y_dat(double ***outp, int etaL, int outw, ofstream& file) { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	int g;
	for (g=0; g<Ng[etaL]; g++) {
		file<<"Grid # "<<etaL<<" Energy Group # "<<g+1<<endl;
		write_cell_edge_y_dat(outp[g], outw, file);
	}
}
//======================================================================================//


//======================================================================================//
//++ function to write multi grid cell edge avg values in output format +++++++++++++++//
//======================================================================================//
void write_grid_average_out(double ****outp, int outw, ofstream& file) {
	int k;
	for (k=0; k<eta_star; k++) {
		if (k==0 ) file<<" -- Fine Energy Grid --\n";
		else if (k==eta_star-1) file<<" -- Grey Energy Grid --\n";
		else file<<" -- Energy Grid # "<<k<<" --\n";
		write_group_average_out(outp[k], k, outw, file);
	}
}
//======================================================================================//
//++ function to write multi grid cell edge x values in viewer friendly format ++++++++++++++++++++//
//======================================================================================//
void write_grid_edge_x_out(double ****outp, int outw, ofstream& file) { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	int k;
	for (k=0; k<eta_star; k++) {
		if (k==0 ) file<<" -- Fine Energy Grid --\n";
		else if (k==eta_star-1) file<<" -- Grey Energy Grid --\n";
		else file<<" -- Energy Grid # "<<k<<" --\n";
		write_group_edge_x_out(outp[k], k, outw, file);
	}
}
//======================================================================================//
//++ function to write multi grid cell edge y values in viewer friendly format ++++++++++++++++++++//
//======================================================================================//
void write_grid_edge_y_out(double ****outp, int outw, ofstream& file) { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	int k;
	for (k=0; k<eta_star; k++) {
		if (k==0 ) file<<" -- Fine Energy Grid --\n";
		else if (k==eta_star-1) file<<" -- Grey Energy Grid --\n";
		else file<<" -- Energy Grid # "<<k<<" --\n";
		write_group_edge_y_out(outp[k], k, outw, file);
	}
}
//======================================================================================//
//++ function to write multi grid cell average values in data format ++++++++++++++++++++++++++++++//
//======================================================================================//
void write_grid_average_dat(double ****outp, int outw, ofstream& file) { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	int k;
	for (k=0; k<eta_star; k++) {
		if (k==0 ) file<<" -- Fine Energy Grid --\n";
		else if (k==eta_star-1) file<<" -- Grey Energy Grid --\n";
		else file<<" -- Energy Grid # "<<k<<" --\n";
		write_group_average_dat(outp[k], k, outw, file);
	}
}
//======================================================================================//
//++ function to write multi grid cell edge x values in data format +++++++++++++++++++++++++++++++//
//======================================================================================//
void write_grid_edge_x_dat(double ****outp, int outw, ofstream& file) { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	int k;
	for (k=0; k<eta_star; k++) {
		if (k==0 ) file<<" -- Fine Energy Grid --\n";
		else if (k==eta_star-1) file<<" -- Grey Energy Grid --\n";
		else file<<" -- Energy Grid # "<<k<<" --\n";
		write_group_edge_x_dat(outp[k], k, outw, file);
	}
}
//======================================================================================//
//++ function to write multi grid cell edge y values in data format +++++++++++++++++++++++++++++++//
//======================================================================================//
void write_grid_edge_y_dat(double ****outp, int outw, ofstream& file) { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	int k;
	for (k=0; k<eta_star; k++) {
		if (k==0 ) file<<" -- Fine Energy Grid --\n";
		else if (k==eta_star-1) file<<" -- Grey Energy Grid --\n";
		else file<<" -- Energy Grid # "<<k<<" --\n";
		write_group_edge_y_dat(outp[k], k, outw, file);
	}
}
//======================================================================================//

//======================================================================================//
//++ function to write multi grid cell average values in column format +++++++++++++++++//
//======================================================================================//
void write_column_average_data(double ***outp, int etaL, int outw, ofstream& file) {
	int i, j, g;
	
	file<<" x index , x grid , y index , y grid ,";
	for (g=0; g<Ng[etaL]; g++) file<<" group "<<g<<",";
	file<<endl;
	for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) {
			file<<i+1<<","<<print_csv((x[i]+x[i+1])/2)<<j+1<<","<<print_csv((y[j]+y[j+1])/2);
			for (g=0; g<Ng[etaL]; g++) file<<print_csv(outp[g][i][j]);
			file<<endl;
		}
	}
}
//======================================================================================//

//======================================================================================//
//++ Function to Recursively Write Iteration data ++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void rec_write_iteration_out(int it, int etaL) {
	int i, g, k, m;
	int num_total, num_start;
	
	k=Ng[etaL]*num_losi[etaL][it];
	outfile<<" # of Inner Iterations "<<num_losi[etaL][it+1]-num_losi[etaL][it]<<endl;
	
	num_start=num_losi[etaL][it];
	num_total=num_losi[etaL][it+1]-num_start;
	
	for (m=0; m<int(num_total/10); m++) {
		outfile<<string((etaL+1)*5,' ')<<"Inner   Iteration  #  :";
		for (i=num_start+m*10; i<num_start+(m+1)*10; i++) outfile<<setw(10)<<i-num_losi[etaL][it]+1<<setw(6)<<" ";
		outfile<<endl;
		outfile<<string((etaL+1)*5,' ')<<"---- Flux      rho    :";
		for (i=num_start+m*10; i<num_start+(m+1)*10; i++) outfile<<print_out(rho_phi[etaL][i]);
		outfile<<endl;
		outfile<<string((etaL+1)*5,' ')<<"Flux Diff. L-inf Norm :";
		for (i=num_start+m*10; i<num_start+(m+1)*10; i++) outfile<<print_out(norm_phi[etaL][i]);
		outfile<<endl;
		outfile<<string((etaL+1)*5,' ')<<"---- K_eff    rho     :";
		for (i=num_start+m*10; i<num_start+(m+1)*10; i++) outfile<<print_out(rho_keff[etaL][i]);
		outfile<<endl;
		outfile<<string((etaL+1)*5,' ')<<"K_eff Diff. L-inf Norm:";
		for (i=num_start+m*10; i<num_start+(m+1)*10; i++) outfile<<print_out(norm_keff[etaL][i]);
		outfile<<endl;
	}
	if ( m*10<num_total ) {
		outfile<<string((etaL+1)*5,' ')<<"Inner   Iteration  #  :";
		for (i=num_start+m*10; i<num_losi[etaL][it+1]; i++) outfile<<setw(10)<<i-num_losi[etaL][it]+1<<setw(6)<<" ";
		outfile<<endl;
		outfile<<string((etaL+1)*5,' ')<<"---- Flux      rho    :";
		for (i=num_start+m*10; i<num_losi[etaL][it+1]; i++) outfile<<print_out(rho_phi[etaL][i]);
		outfile<<endl;
		outfile<<string((etaL+1)*5,' ')<<"Flux Diff. L-inf Norm :";
		for (i=num_start+m*10; i<num_losi[etaL][it+1]; i++) outfile<<print_out(norm_phi[etaL][i]);
		outfile<<endl;
		outfile<<string((etaL+1)*5,' ')<<"---- K_eff    rho     :";
		for (i=num_start+m*10; i<num_losi[etaL][it+1]; i++) outfile<<print_out(rho_keff[etaL][i]);
		outfile<<endl;
		outfile<<string((etaL+1)*5,' ')<<"K_eff Diff. L-inf Norm:";
		for (i=num_start+m*10; i<num_losi[etaL][it+1]; i++) outfile<<print_out(norm_keff[etaL][i]);
		outfile<<endl;
	}
	if (etaL+1!=eta_star) {
		for (i=num_losi[etaL][it]; i<num_losi[etaL][it+1]; i++) {
			outfile<<string((etaL+2)*5,' ')<<"Iteration # "<<i-num_losi[etaL][it]+1;
			rec_write_iteration_out(i, etaL+1);
		}
	}
}
//======================================================================================//
void rec_write_iteration_dat(int it, int etaL) {
	int i, g, k;
	
	k=Ng[etaL]*num_losi[etaL][it];
	datfile<<" # of Inner Iterations ,"<<num_losi[etaL][it+1]-num_losi[etaL][it]<<","<<endl;
	datfile<<string((etaL+1)*5,' ')<<"Inner   Iteration  #  :,";
	for (i=num_losi[etaL][it]; i<num_losi[etaL][it+1]; i++) datfile<<setw(10)<<i-num_losi[etaL][it]+1<<setw(6)<<",";
	datfile<<endl;
	datfile<<string((etaL+1)*5,' ')<<"---- Flux      rho    :,";
	for (i=num_losi[etaL][it]; i<num_losi[etaL][it+1]; i++) datfile<<print_csv(rho_phi[etaL][i]);
	datfile<<endl;
	datfile<<string((etaL+1)*5,' ')<<"Flux Diff. L-inf Norm :,";
	for (i=num_losi[etaL][it]; i<num_losi[etaL][it+1]; i++) datfile<<print_csv(norm_phi[etaL][i]);
	datfile<<endl;
	datfile<<string((etaL+1)*5,' ')<<"---- K_eff    rho     :,";
	for (i=num_losi[etaL][it]; i<num_losi[etaL][it+1]; i++) datfile<<print_csv(rho_keff[etaL][i]);
	datfile<<endl;
	datfile<<string((etaL+1)*5,' ')<<"K_eff Diff. L-inf Norm:,";
	for (i=num_losi[etaL][it]; i<num_losi[etaL][it+1]; i++) datfile<<print_csv(norm_keff[etaL][i]);
	datfile<<endl;
	if (etaL+1!=eta_star) {
		for (i=num_losi[etaL][it]; i<num_losi[etaL][it+1]; i++) {
			datfile<<string((etaL+2)*5,' ')<<"Iteration # "<<i-num_losi[etaL][it]+1;
			rec_write_iteration_dat(i, etaL+1);
		}
	}
}
//======================================================================================//
void rec_write_iteration_long_dat(int it, int etaL) {
	int i, g, k;
	
	k=Ng[etaL]*num_losi[etaL][it];
	for (i=num_losi[etaL][it]; i<num_losi[etaL][it+1]; i++) {
		for (g=0; g<Ng[etaL]; g++) {
			datfile<<string((etaL+1)*5,' ')<<"Pre cond. time ,"<<print_csv(dt_pc[etaL][k])<<" # of BiCGStab , "
			<<num_logm[etaL][k]<<" , Err. in LO Sol. ,"<<print_csv(err_lo[etaL][k])<<endl;
			k++;
		}
		datfile<<string((etaL+1)*5,' ')<<"Convergence Rate, "<<print_csv(rho_phi[etaL][i])<<endl;
		if (etaL+1==eta_star) continue;
		else rec_write_iteration_long_dat(i, etaL+1);
	}
}
//======================================================================================//

//======================================================================================//
//++ function to output data +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void output(char fn[])
{
	char outf[50]={""}, temp[10]={""};
	int p, g, gg, i, j, m, k, outw=16;
	double phi_sum=0.0, psi_max=0.0, res;
	double conv, conv_p, conv_pl, conv_pr, conv_pb, conv_pt, conv_jx, conv_jy;
	double L2_norm, L2_p, L2_pl, L2_pr, L2_pb, L2_pt, L2_jx, L2_jy;
	double *res_x=0, *res_y=0, *res_lbc=0, *res_rbc=0, *res_bbc=0, *res_tbc=0;
	int *i_xx=0, *i_yy=0, *i_bbc=0, *i_tbc=0;
	int *j_xx=0, *j_yy=0, *j_lbc=0, *j_rbc=0;
	int *g_xx=0, *g_yy=0, *g_bbc=0, *g_tbc=0, *g_lbc=0, *g_rbc=0;
	double res_hoit, *res_loit, sFission, sScatter; // Iterative Residuals
	int i_hoit, j_hoit, g_hoit, *i_loit, *j_loit, *g_loit; // Iterative Residuals Location

	cout<<"Begin Output\n";
	
	outfile<<"\n -- Material Map -- \n";
	for (j=ny-1; j>=0; j--) {
		for (i=0; i<nx; i++) outfile<<setw(3)<<material[i][j];
		outfile<<endl;
	}
	
	// calculate residuals ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	res_x  =new double[eta_star]; res_y  =new double[eta_star];
	res_lbc=new double[eta_star]; res_rbc=new double[eta_star];
	res_bbc=new double[eta_star]; res_tbc=new double[eta_star];
	i_xx   =new int[eta_star]; j_xx   =new int[eta_star]; g_xx   =new int[eta_star];
	i_yy   =new int[eta_star]; j_yy   =new int[eta_star]; g_yy   =new int[eta_star];
	j_lbc=new int[eta_star]; g_lbc=new int[eta_star]; j_rbc=new int[eta_star]; g_rbc=new int[eta_star];
	i_bbc=new int[eta_star]; g_bbc=new int[eta_star]; i_tbc=new int[eta_star]; g_tbc=new int[eta_star];
	res_loit=new double[eta_star]; i_loit=new int[eta_star]; j_loit=new int[eta_star]; g_loit=new int[eta_star];
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Equation residuals ///////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	for (k=0; k<eta_star; k++) {
		res_x[k]=0.0; res_y[k]=0.0;	res_lbc[k]=0.0; res_rbc[k]=0.0; res_bbc[k]=0.0; res_tbc[k]=0.0;
		for (g=0; g<Ng[k]; g++) {
			for (i=1; i<nx-1; i++) { // First Moment X Grid
				for (j=0; j<ny; j++) {
					res=abs(j_x[k][g][i][j]+(D_xP[k][g][i][j]*phi[k][g][i][j]-D_xN[k][g][i][j]*phi[k][g][i-1][j])/xe[i]); 
					if ( res>res_x[k] ) { res_x[k]=res; i_xx[k]=i; j_xx[k]=j; g_xx[k]=g; }
				}
			}
			for (i=0; i<nx; i++) { // First Moment Y Grid
				for (j=1; j<ny-1; j++) {
					res=abs(j_y[k][g][i][j]+(D_yP[k][g][i][j]*phi[k][g][i][j]-D_yN[k][g][i][j]*phi[k][g][i][j-1])/ye[j]); 
					if ( res>res_y[k] ) { res_y[k]=res; i_yy[k]=i; j_yy[k]=j; g_yy[k]=g; }
				}
			}
			// Boundary Conditions
			for (j=0; j<ny; j++) { // Left and Right
				res=abs(FL[k][g][j]*phiB_L[k][g][j]-j_x[k][g][0][j]);  
				if ( res>res_lbc[k] ) { res_lbc[k]=res; j_lbc[k]=j; g_lbc[k]=g; } // Left BC residual
				res=abs(FR[k][g][j]*phiB_R[k][g][j]-j_x[k][g][nx][j]); 
				if ( res>res_rbc[k] ) { res_rbc[k]=res; j_rbc[k]=j; g_rbc[k]=g; } // Right BC residual
			}
			for (i=0; i<nx; i++) { // Bottom and Top
				res=abs(FB[k][g][i]*phiB_B[k][g][i]-j_y[k][g][i][0]);  
				if ( res>res_bbc[k] ) { res_bbc[k]=res; i_bbc[k]=i; g_bbc[k]=g; } // Bottom BC residual
				res=abs(FT[k][g][i]*phiB_T[k][g][i]-j_y[k][g][i][ny]); 
				if ( res>res_tbc[k] ) { res_tbc[k]=res; i_tbc[k]=i; g_tbc[k]=g; } // Top BC residual
			}
		}
		// Calculate Iterative Residuals // Residuals of the actual NDA Equaitons
		res_loit[k]=0.0;
		for (i=1; i<nx-1; i++) { 
			for (j=0; j<ny; j++) {
				sFission=0;
				for (g=0; g<Ng[k]; g++) sFission+=nuSigmaF[k][g][i][j]*phi[k][g][i][j];
				for (g=0; g<Ng[k]; g++) {
					sScatter=0;
					for (gg=0; gg<Ng[k]; gg++) sScatter+=sigmaS[k][gg][g][i][j]*phi[k][gg][i][j];
					res=abs((j_x[k][g][i+1][j]-j_x[k][g][i][j])*hy[j]+(j_y[k][g][i][j+1]-j_y[k][g][i][j])*hx[i]
					+hx[i]*hy[j]*(sigmaT[k][g][i][j]*phi[k][g][i][j]-sScatter-chi[k][g][i][j]*sFission/k_eff-s_ext[k][g][i][j]));
					if ( res>res_loit[k] ) { res_loit[k]=res; i_loit[k]=i; j_loit[k]=j; g_loit[k]=g; } // residual
				}
			}
		}
	}
	
	// Calculate Iterative Residuals for High-Order Problem // Residuals of the actual NDA Equaitons
	res_hoit=0.0; k=0;
	for (i=1; i<nx-1; i++) { 
		for (j=0; j<ny; j++) {
			sFission=0;
			for (g=0; g<Ng[k]; g++) sFission+=nuSigmaF[k][g][i][j]*phiT[g][i][j];
			for (g=0; g<Ng[k]; g++) {
				sScatter=0;
				for (gg=0; gg<Ng[k]; gg++) sScatter+=sigmaS[k][gg][g][i][j]*phiT[gg][i][j];
				res=abs((j_xT[g][i+1][j]-j_xT[g][i][j])*hy[j]+(j_yT[g][i][j+1]-j_yT[g][i][j])*hx[i]
				+hx[i]*hy[j]*(sigmaT[k][g][i][j]*phiT[g][i][j]-sScatter-chi[k][g][i][j]*sFission/k_eff-s_ext[k][g][i][j]));
				if ( res>res_hoit ) { res_hoit=res; i_hoit=i; j_hoit=j; g_hoit=g; } // residual
			}
		}
	}
	
	// file output ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Energy Solution Grid
	outfile<<"\n -- Energy Grid -- \n";
	outfile<<"# of Fine Energy Groups : "<<Ng[0]<<endl;
	for (k=1; k<eta_star-1; k++) {
		outfile<<"# of Coarse Energy Groups : "<<Ng[k]<<" with group combinations ";
		for (g=0; g<Ng[k]; g++) outfile<<" "<<omegaP[k][g+1]-omegaP[k][g];
		outfile<<endl;
	}
	outfile<<"One Group Grid "<<" with group combination "<<omegaP[eta_star-1][1]-omegaP[eta_star-1][0]<<endl;
	
	// boundary conditions
	outfile<<"\n -- Boundary Conditions -- \n";
	outfile<<"Type of BC: "<<kbc;   // Write type of BC
	
	// output grid
	outfile<<"\n -- Solution Grid -- \n";
	outfile<<" X grid \n"<<" Index     Cell Edge       Width Avg     Cell center      Cell Width   ";
	switch ( kbc ) {
	case 1:
		outfile<<"  Bottom BC In  "; // Write Bottom BC
		outfile<<"  Top    BC In  "; // Write Top BC
		break;
	case 2:
		outfile<<"  Quad 2 BC In  "; // Write Bottom BC
		outfile<<"  Quad 4 BC In  "; // Write Top BC
		break;
	case 3:
		outfile<<" Bottom BC REFL "; // Write Bottom BC
		outfile<<" Top    BC REFL "; // Write Top BC
		break;
	case 4:
		outfile<<" Bottom BC REFL "; // Write Bottom BC
		outfile<<"  Top    BC In  "; // Write Top BC
		break;
	case 5:
		outfile<<"  Bottom BC In  "; // Write Bottom BC
		outfile<<" Top    BC REFL "; // Write Top BC
		break;
	default:
		cout<<">>Invalid BC Type !!!!!\n";
		break;
	}
	outfile<<endl;
	for (i=0; i<nx; i++) outfile<<setw(6)<<i+1<<print_out(x[i])<<print_out(xe[i])<<print_out((x[i]+x[i+1])/2)<<print_out(hx[i])
	<<print_out(bcB[0][i])<<print_out(bcT[0][i])<<endl; // write x grid
	outfile<<setw(6)<<nx+1<<print_out(x[nx])<<print_out(xe[nx])<<endl;
	outfile<<" Y grid \n"<<" Index     Cell Edge       Width Avg     Cell center      Cell Width   ";
	switch ( kbc ) {
	case 1:
		outfile<<"  Left   BC In  "; // Write Left BC
		outfile<<"  Right  BC In  "; // Write Right BC
		break;
	case 2:
		outfile<<"  Quad 1 BC In  "; // Write Left BC
		outfile<<"  Quad 3 BC In  "; // Write Right BC
		break;
	case 3:
		outfile<<" Left   BC REFL "; // Write Left BC
		outfile<<" Right  BC REFL "; // Write Right BC
		break;
	case 4:
		outfile<<" Left   BC REFL "; // Write Left BC
		outfile<<"  Right  BC In  "; // Write Right BC
		break;
	case 5:
		outfile<<"  Left   BC In  "; // Write Left BC
		outfile<<" Right  BC REFL "; // Write Right BC
		break;
	default:
		cout<<">>Invalid BC Type !!!!!\n";
		break;
	}
	outfile<<endl;
	for (j=0; j<ny; j++) outfile<<setw(6)<<j+1<<print_out(y[j])<<print_out(ye[j])<<print_out((y[j]+y[j+1])/2)<<print_out(hy[j])
	<<print_out(bcL[0][j])<<print_out(bcR[0][j])<<endl; // write y grid
	outfile<<setw(6)<<ny+1<<print_out(y[ny])<<print_out(ye[ny])<<endl;
	
	// output quadrature
	outfile<<"\n -- Quadrature -- \n";
	if ( N==20 or N==36 ) outfile<<" Octant-Range Q";
	else outfile<<"Level Symmeteric S";
	outfile<<N<<" Normalized to Integrate to 4*pi \n";
	outfile<<"   m        mu              eta              xi            weight   \n";
	for (m=0; m<sn; m++) outfile<<setw(5)<<m+1<<print_out(mu[m])<<print_out(eta[m])<<print_out(xi[m])<<print_out(w[m])<<endl; // quadrature
	outfile<<endl;
	
	outfile<<endl;
	outfile<<" -------------- \n";
	outfile<<" -- Solution -- \n";
	outfile<<" -------------- \n";
	outfile<<"Program run time : "<<((float)t)/CLOCKS_PER_SEC<<" seconds"<<endl;
	outfile<<"\nNumber of iterations "<<n_iterations<<endl;
	
	for (k=0; k<eta_star; k++) {
		rho_phi[k].erase(rho_phi[k].begin());
		rho_keff[k].erase(rho_keff[k].begin());
		norm_phi[k].erase(norm_phi[k].begin());
		norm_keff[k].erase(norm_keff[k].begin());
	}
	
	outfile<<"\n -- Iteration Data -- \n";
	outfile<<"  Iter."<<"    Estimated   "<<"    High-order  "<<"     rho of     "<<"   k_effetive   "<<"  High-order    "<<"   Low-order   ";
	for (k=0; k<eta_star; k++) outfile<<" Grid # "<<setw(2)<<k<<" ";
	outfile<<"   Number of   \n";
	outfile<<"    #  "<<"  Spectral Rad. "<<" Diff L-inf Norm"<<"  k_effective   "<<" Diff L-inf Norm"<<" Sol. Time [s]  "<<" Sol. Time [s] ";
	for (k=0; k<eta_star; k++) outfile<<"iterations ";
	outfile<<"   Matrix Sol.  \n";
	for (i=0; i<n_iterations; i++) {
		outfile<<setw(4)<<i<<setw(2)<<" "<<print_out(rho_ho[i+1])<<print_out(norm_ho[i])<<print_out(rho_kho[i+1])<<print_out(norm_kho[i])<<
		print_out(dt_ho[i])<<print_out(dt_lo[0][i]);
		for (k=0; k<eta_star; k++) outfile<<setw(7)<<num_grid[k][i]<<setw(4)<<" ";
		outfile<<setw(10)<<num_mtot[i]<<setw(6)<<endl;
	}
	// Output Detailed Convergence data
	for (k=0; k<eta_star; k++) for (i=0; i<num_losi[k].size()-1; i++) num_losi[k][i+1]+=num_losi[k][i]; // Cumulative sum of all the iterations
	outfile<<" -- Nested Convergence Data -- \n";
	for (i=0; i<n_iterations; i++) {
		outfile<<"Transport Iteration # "<<i;
		rec_write_iteration_out(i, 0);
	}
	
	outfile<<"\n -- Residuals -- \n";
	outfile<<"High-order Residual \n";
	outfile<<"Balance Residual:"<<print_out(res_ho)<<" at "<<i_ho<<" , "<<j_ho<<endl;
	
	outfile<<"Low-order Residuals \n";
	for (k=0; k<eta_star; k++) {
		outfile<<"Energy Grid # "<<k<<endl;
		outfile<<"Residual of Equations Solved In Matrix \n";
		outfile<<"Cell Balance Residual:"<<print_out(res_mbal[k])<<" in group "<<g_mbal[k]<<" at "<<i_mbal[k]<<" , "<<j_mbal[k]<<endl; 
		outfile<<"Left    BC   Residual:"<<print_out(res_ml[k])<<" in group "<<g_ml[k]<<" at "<<j_ml[k]<<endl;
		outfile<<"Right   BC   Residual:"<<print_out(res_mr[k])<<" in group "<<g_mr[k]<<" at "<<j_mr[k]<<endl;
		outfile<<"Bottom  BC   Residual:"<<print_out(res_mb[k])<<" in group "<<g_mb[k]<<" at "<<i_mb[k]<<endl;
		outfile<<"Top     BC   Residual:"<<print_out(res_mt[k])<<" in group "<<g_mt[k]<<" at "<<i_mt[k]<<endl;
		
		outfile<<"Residual of General Equations \n";
		outfile<<"Cell Balance Residual:"<<print_out(res_bal[k])<<" in group "<<g_bal[k]<<" at "<<i_bal[k]<<" , "<<j_bal[k]<<endl; 
		outfile<<"X       Grid Residual:"<<print_out(res_x[k])<<" in group "<<g_xx[k]<<" at "<<i_xx[k]<<" , "<<j_xx[k]<<endl;
		outfile<<"Y       Grid Residual:"<<print_out(res_y[k])<<" in group "<<g_yy[k]<<" at "<<i_yy[k]<<" , "<<j_yy[k]<<endl;
		outfile<<"Left    BC   Residual:"<<print_out(res_lbc[k])<<" in group "<<g_lbc[k]<<" at "<<j_lbc[k]<<endl;
		outfile<<"Right   BC   Residual:"<<print_out(res_rbc[k])<<" in group "<<g_rbc[k]<<" at "<<j_rbc[k]<<endl;
		outfile<<"Bottom  BC   Residual:"<<print_out(res_bbc[k])<<" in group "<<g_bbc[k]<<" at "<<i_bbc[k]<<endl;
		outfile<<"Top     BC   Residual:"<<print_out(res_tbc[k])<<" in group "<<g_tbc[k]<<" at "<<i_tbc[k]<<endl;
	}
	
	outfile<<" -- Iterative Residuals -- \n";
	outfile<<"High-order Iterative Residual "<<print_out(res_hoit)<<" in group "<<g_hoit<<" at "<<i_hoit<<" , "<<j_hoit<<endl;
	outfile<<"Low-order Iterative Residuals \n";
	for (k=0; k<eta_star; k++) outfile<<"Grid "<<k<<" Residual "<<print_out(res_loit[k])<<" in group "<<g_loit[k]<<" at "<<i_loit[k]<<" , "<<j_loit[k]<<endl;
	
	// ++ Evaluate Consistency between solutions +++++++++++++++++++++++++++++
	// Find difference between Transport solution and NDA solution
	conv_p=0.0; conv_pl=0.0; conv_pr=0.0; conv_pb=0.0; conv_pt=0.0; conv_jx=0.0; conv_jy=0.0;
	for (g=0; g<Ng[0]; g++) {
		for ( i=0; i<nx; i++ ) {
			for ( j=0; j<ny; j++ ) { conv=abs((phi[0][g][i][j]-phiT[g][i][j])/phiT[g][i][j]); if ( conv>conv_p ) conv_p=conv; }
		}
		for ( j=0; j<ny; j++ ) {
			conv=abs((phiB_L[0][g][j]-phi_xT[g][0][j] )/phi_xT[g][0][j]);  if ( conv>conv_pl ) conv_pl=conv;
			conv=abs((phiB_R[0][g][j]-phi_xT[g][nx][j])/phi_xT[g][nx][j]); if ( conv>conv_pr ) conv_pr=conv;
		}
		for ( i=0; i<nx; i++ ) {
			conv=abs((phiB_B[0][g][i]-phi_yT[g][i][0] )/phi_yT[g][i][0]);  if ( conv>conv_pb ) conv_pb=conv;
			conv=abs((phiB_T[0][g][i]-phi_yT[g][i][ny])/phi_yT[g][i][ny]); if ( conv>conv_pt ) conv_pt=conv;
		}
		for ( i=0; i<nx+1; i++ ) {
			for ( j=0; j<ny; j++ ) { conv=abs((j_x[0][g][i][j]-j_xT[g][i][j])/j_xT[g][i][j]); if ( conv>conv_jx ) conv_jx=conv; }
		}

		for ( i=0; i<nx; i++ ) {
			for ( j=0; j<ny+1; j++ ) { conv=abs((j_y[0][g][i][j]-j_yT[g][i][j])/j_yT[g][i][j]); if ( conv>conv_jy ) conv_jy=conv; }
		}
	}
	outfile<<"\n -- L-Infinity Norm of Relative Difference Between Transport and NDA Solution -- \n";
	outfile<<"Cell      Averaged Flux : "<<print_out(conv_p)<<endl;
	outfile<<" L        Boundary Flux : "<<print_out(conv_pl)<<endl;
	outfile<<" R        Boundary Flux : "<<print_out(conv_pr)<<endl;
	outfile<<" B        Boundary Flux : "<<print_out(conv_pb)<<endl;
	outfile<<" T        Boundary Flux : "<<print_out(conv_pt)<<endl;
	outfile<<"X Grid Face Avg Current : "<<print_out(conv_jx)<<endl;
	outfile<<"Y Grid Face Avg Current : "<<print_out(conv_jy)<<endl;
	
	
	outfile<<"\n -- Consistency Between Low-order Solutions on Successive Grids -- \n";
	for (k=1; k<eta_star; k++) {
		conv_p=0.0; conv_pl=0.0; conv_pr=0.0; conv_pb=0.0; conv_pt=0.0; conv_jx=0.0; conv_jy=0.0;
		L2_p=0.0; L2_pl=0.0; L2_pr=0.0; L2_pb=0.0; L2_pt=0.0; L2_jx=0.0; L2_jy=0.0;
		for (p=0; p<Ng[k]; p++) {
			for ( i=0; i<nx; i++ ) {
				for ( j=0; j<ny; j++ ) {
					phi_sum=0.0; for (g=omegaP[k][p]; g<omegaP[k][p+1]; g++) phi_sum+=phi[k-1][g][i][j];
					conv=abs((phi[k][p][i][j]-phi_sum)/phi_sum); if (conv>conv_p) conv_p=conv;
					L2_p+=conv*conv;
				}
			}
			for ( j=0; j<ny; j++ ) {
				phi_sum=0.0; for (g=omegaP[k][p]; g<omegaP[k][p+1]; g++) phi_sum+=phiB_L[k-1][g][j];
				conv=abs((phiB_L[k][p][j]-phi_sum)/phi_sum); if (conv>conv_pl) conv_pl=conv;
				L2_pl+=conv*conv;
				phi_sum=0.0; for (g=omegaP[k][p]; g<omegaP[k][p+1]; g++) phi_sum+=phiB_R[k-1][g][j];
				conv=abs((phiB_R[k][p][j]-phi_sum)/phi_sum); if (conv>conv_pr) conv_pr=conv;
				L2_pr+=conv*conv;
			}
			for ( i=0; i<nx; i++ ) {
				phi_sum=0.0; for (g=omegaP[k][p]; g<omegaP[k][p+1]; g++) phi_sum+=phiB_B[k-1][g][i];
				conv=abs((phiB_B[k][p][i]-phi_sum)/phi_sum); if (conv>conv_pb) conv_pb=conv;
				L2_pb+=conv*conv;
				phi_sum=0.0; for (g=omegaP[k][p]; g<omegaP[k][p+1]; g++) phi_sum+=phiB_T[k-1][g][i];
				conv=abs((phiB_T[k][p][i]-phi_sum)/phi_sum); if (conv>conv_pt) conv_pt=conv;
				L2_pt+=conv*conv;
			}
			for ( i=0; i<nx+1; i++ ) {
				for ( j=0; j<ny; j++ ) {
					phi_sum=0.0; for (g=omegaP[k][p]; g<omegaP[k][p+1]; g++) phi_sum+=j_x[k-1][g][i][j]; phi_sum+=1e-15;
					conv=abs((j_x[k][p][i][j]-phi_sum)/phi_sum); if (conv>conv_jx) conv_jx=conv;
					L2_jx+=conv*conv;
				}
			}
			for ( i=0; i<nx; i++ ) {
				for ( j=0; j<ny+1; j++ ) {
					phi_sum=0.0; for (g=omegaP[k][p]; g<omegaP[k][p+1]; g++) phi_sum+=j_y[k-1][g][i][j]; phi_sum+=1e-15;
					conv=abs((j_y[k][p][i][j]-phi_sum)/phi_sum); if (conv>conv_jy) conv_jy=conv;
					L2_jy+=conv*conv;
				}
			}
		}
		
		outfile<<"Relative Difference Between Grid # "<<k-1<<" and Grid # "<<k<<endl;
		outfile<<"                            "<<"L-infinity Norm"<<"    L2 Norm     "<<endl;
		outfile<<"Cell      Averaged Flux : "<<print_out(conv_p)<<print_out(sqrt(L2_p))<<endl;
		outfile<<" L        Boundary Flux : "<<print_out(conv_pl)<<print_out(sqrt(L2_pl))<<endl;
		outfile<<" R        Boundary Flux : "<<print_out(conv_pr)<<print_out(sqrt(L2_pr))<<endl;
		outfile<<" B        Boundary Flux : "<<print_out(conv_pb)<<print_out(sqrt(L2_pb))<<endl;
		outfile<<" T        Boundary Flux : "<<print_out(conv_pt)<<print_out(sqrt(L2_pt))<<endl;
		outfile<<"X Grid Face Avg Current : "<<print_out(conv_jx)<<print_out(sqrt(L2_jx))<<endl;
		outfile<<"Y Grid Face Avg Current : "<<print_out(conv_jy)<<print_out(sqrt(L2_jy))<<endl;
	}
	if (KE_problem) {
		outfile<<"\n -- K-eigenvalue -- \n";
		outfile<<"K_eff = "<<print_out(k_eff)<<endl;
	}
	
	// output flux
	outfile<<"\n ------------------- ";
	outfile<<"\n -- NDA Solution -- "; // Write NDA solution
	outfile<<"\n ------------------- \n";
	outfile<<"\n -- Cell Averaged Scalar Flux -- \n";
	write_grid_average_out(phi, outw, outfile); // call function to write out cell average scalar flux
	
	outfile<<"\n -- X Face Average Normal Current J_x -- \n";
	write_grid_edge_x_out(j_x, outw, outfile); // call function to write out cell edge current on x grid
	
	outfile<<"\n -- Y Face Average Normal Current J_y -- \n";
	write_grid_edge_y_out(j_y, outw, outfile); // call function to write out cell edge scalar flux on y grid
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	outfile<<"\n ------------------------- ";
	outfile<<"\n -- High Order Solution -- "; // Write Step Characteristics solution
	outfile<<"\n ------------------------- \n";
	outfile<<"\n -- Cell Averaged Scalar Flux -- \n";
	write_group_average_out(phiT, 0, outw, outfile); // call function to write out cell average scalar flux
	
	//outfile<<"\n -- X Grid Face Averaged Scalar Flux -- \n";
	//write_group_edge_x_out(phi_xT, 0, outw, outfile); // call function to write out cell edge scalar flux on x grid
	
	//outfile<<"\n -- Y Grid Face Averaged Scalar Flux -- \n";
	//write_group_edge_y_out(phi_yT, 0, outw, outfile); // call function to write out cell edge scalar flux on y grid
	
	outfile<<"\n -- X Face Average Normal Current J_x -- \n";
	write_group_edge_x_out(j_xT, 0, outw, outfile); // call function to write out cell edge current on x grid
	
	outfile<<"\n -- Y Face Average Normal Current J_y -- \n";
	write_group_edge_y_out(j_yT, 0, outw, outfile); // call function to write out cell edge scalar flux on y grid
	
	outfile.close(); // close output file opened in input file +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	outf[0]='\0';
	strcat (outf,fn); strcat (outf,".csv"); // name output file
	cout<<outf<<endl;
	datfile.open(outf); // open output data file ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// output flux
	datfile<<"# of iterations ,"<<n_iterations<<endl;
	datfile<<" -- Iteration Data -- \n";
	datfile<<"  Iter.,   Convergence,  High-order   ,   Low-order   \n";
	datfile<<"    #  ,      Rate    , Sol. Time [s] ,  Sol. Time [s]\n";
	for (i=0; i<n_iterations; i++) datfile<<setw(4)<<i<<setw(2)<<" "<<","<<print_csv(rho_ho[i])<<
		print_csv(dt_ho[i])<<print_csv(dt_lo[0][i])<<endl;
	datfile<<"Number of High-order Iterations, "<<n_iterations<<endl;
	// Output Detailed Convergence data
	for (i=0; i<n_iterations; i++) {
		datfile<<"Transport Iteration # "<<i<<",";
		rec_write_iteration_dat(i, 0);
	}
	datfile<<"Number of High-order Iterations, "<<n_iterations<<endl;
	// Output Detailed Convergence data
	for (i=0; i<n_iterations; i++) {
		datfile<<"High-order Iteration #, "<<i<<endl;
		datfile<<"High-order Solution Time, "<<print_csv(dt_ho[i])<<endl;
		datfile<<"Low- order Solution Time, "<<print_csv(dt_lo[0][i])<<endl;
		rec_write_iteration_long_dat(i, 0);
	}
	
	
	datfile<<" # of x cells , # of y cells ,\n";
	datfile<<nx<<" , "<<ny<<", \n";
	
	// output grid
	datfile<<" -- Solution Grid -- \n";
	
	datfile<<" X grid \n"<<" Index   ,  Cell Edge   ,    Width Avg   ,  Cell Center   ,   Cell Width   "<<endl;
	for (i=0; i<nx; i++) datfile<<setw(6)<<i+1<<","<<print_csv(x[i])<<print_csv(xe[i])<<print_csv((x[i]+x[i+1])/2)<<print_csv(hx[i])<<endl; // write x grid
	datfile<<setw(6)<<nx+1<<","<<print_csv(x[nx])<<print_csv(xe[nx])<<endl;
	
	datfile<<" Y grid \n"<<" Index   ,  Cell Edge   ,    Width Avg   ,  Cell Center   ,   Cell Width   "<<endl;
	for (j=0; j<ny; j++) datfile<<setw(6)<<j+1<<","<<print_csv(y[j])<<print_csv(ye[j])<<print_csv((y[j]+y[j+1])/2)<<print_csv(hy[j])<<endl; // write y grid
	datfile<<setw(6)<<ny+1<<","<<print_csv(y[ny])<<print_csv(ye[ny])<<endl;
	
	datfile.close(); // close output data file ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	
	outf[0]='\0';
	strcat (outf,fn); strcat (outf,".lo.csv"); // name output file
	cout<<outf<<endl;
	datfile.open(outf); // open output data file ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	datfile<<" # of x cells , # of y cells ,\n";
	datfile<<nx<<" , "<<ny<<", \n";
	datfile<<"x edge grid , "; for (i=0; i<nx+1; i++) datfile<<print_csv(x[i]); datfile<<endl;
	datfile<<"x center grid , "; for (i=0; i<nx; i++) datfile<<print_csv((x[i]+x[i+1])/2); datfile<<endl;
	datfile<<"y edge grid , "; for (j=0; j<ny+1; j++) datfile<<print_csv(y[j]); datfile<<endl;
	datfile<<"y center grid , "; for (j=0; j<ny; j++) datfile<<print_csv((y[j]+y[j+1])/2); datfile<<endl;
	datfile<<"number of energy grids, "<<eta_star<<",\n";
	datfile<<"number of groups in grid, ";
	for (k=0; k<eta_star; k++) datfile<<Ng[k]<<" , ";
	datfile<<endl;
	
	datfile<<" -- Cell Averaged Scalar Flux -- \n";
	write_grid_average_dat(phi, outw, datfile);
	
	datfile<<" -- X Face Average Normal Current J_x -- \n";
	write_grid_edge_x_dat(j_x, outw, datfile); // call function to write out cell edge current on x grid
	
	datfile<<" -- Y Face Average Normal Current J_y -- \n";
	write_grid_edge_y_dat(j_y, outw, datfile); // call function to write out cell edge current on y grid
	
	datfile<<" ----------------------- \n";
	datfile<<" -- Consistency Terms -- \n"; // 
	datfile<<" ----------------------- \n";
	
	datfile<<" -- X Grid Positive Diffusion Coefficient D^+_x -- \n";
	write_grid_edge_x_dat(D_xP, outw, datfile);
	
	datfile<<" -- Y Grid Positive Diffusion Coefficient D^+_y -- \n";
	write_grid_edge_y_dat(D_yP, outw, datfile);
	
	datfile<<" -- X Grid Negative Diffusion Coefficient D^-_x -- \n";
	write_grid_edge_x_dat(D_xN, outw, datfile);
	
	datfile<<" -- Y Grid Negative Diffusion Coefficient D^-_y -- \n";
	write_grid_edge_y_dat(D_yN, outw, datfile);
	
	datfile<<"\n -- Left and Right Boundary Flux -- \n";
	for (k=0; k<eta_star; k++) {
		datfile<<" -- Energy Grid # "<<k<<endl;
		for (g=0; g<Ng[k]; g++) {
			datfile<<"Grid # "<<k<<" Energy Group # "<<g<<endl;
			datfile<<"  index    sol. grid     Left   Flux     Right  Flux  \n";
			for ( j=ny-1; j>=0; j-- ) datfile<<setw(6)<<j+1<<","<<print_csv((y[j]+y[j+1])/2)<<print_csv(phiB_L[k][g][j])<<print_csv(phiB_R[k][g][j])<<endl;
		}
	}
	
	datfile<<"\n -- Bottom and Top Boundary Flux -- \n";
	datfile<<"  index    sol. grid     Bottom  Flux     Top   Flux  \n";
	for (k=0; k<eta_star; k++) {
		datfile<<" -- Energy Grid # "<<k<<endl;
		for (g=0; g<Ng[k]; g++) {
			datfile<<"Grid # "<<k<<" Energy Group # "<<g<<endl;
			datfile<<"  index    sol. grid     Bottom  Flux     Top   Flux  \n";
			for ( i=nx-1; i>=0; i-- ) datfile<<setw(6)<<i+1<<","<<print_csv((x[i]+x[i+1])/2)<<print_csv(phiB_B[k][g][i])<<print_csv(phiB_T[k][g][i])<<endl;
		}
	}
	
	// Boundary conditions
	datfile<<"\n -- Bottom and Top Boundary Factors -- \n";
	for (k=0; k<eta_star; k++) {
		datfile<<" -- Energy Grid # "<<k<<endl;
		for (g=0; g<Ng[k]; g++) {
			datfile<<"Grid # "<<k<<" Energy Group # "<<g<<endl;
			datfile<<" index "<<"     x center    "<<"   F  Bottom    "<<"     F  Top     \n";
			for (i=0; i<nx; i++) datfile<<setw(6)<<i+1<<","<<print_csv((x[i]+x[i+1])/2)<<print_csv(FB[k][g][i])<<print_csv(FT[k][g][i])<<endl;
		}
	}
	
	datfile<<"\n -- Left and Right Boundary Factors -- \n";
	for (k=0; k<eta_star; k++) {
		datfile<<" -- Energy Grid # "<<k<<endl;
		for (g=0; g<Ng[k]; g++) {
			datfile<<"Grid # "<<k<<" Energy Group # "<<g<<endl;
			datfile<<" index "<<"     y center    "<<"    F  Left     "<<"    F  Right    \n";
			for (j=0; j<ny; j++) datfile<<setw(6)<<j+1<<","<<print_csv((y[j]+y[j+1])/2)<<print_csv(FL[k][g][j])<<print_csv(FR[k][g][j])<<endl;
		}
	}
	
	datfile.close(); // close output data file ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	outf[0]='\0';
	strcat (outf,fn); strcat (outf,".ho.csv"); // name output file
	cout<<outf<<endl;
	datfile.open(outf); // open output data file ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// output flux
	
	datfile<<" # of x cells , # of y cells ,\n";
	datfile<<nx<<" , "<<ny<<", \n";
	datfile<<"x edge grid , "; for (i=0; i<nx+1; i++) datfile<<print_csv(x[i]); datfile<<endl;
	datfile<<"x center grid , "; for (i=0; i<nx; i++) datfile<<print_csv((x[i]+x[i+1])/2); datfile<<endl;
	datfile<<"y edge grid , "; for (j=0; j<ny+1; j++) datfile<<print_csv(y[j]); datfile<<endl;
	datfile<<"y center grid , "; for (j=0; j<ny; j++) datfile<<print_csv((y[j]+y[j+1])/2); datfile<<endl;
	datfile<<"number of energy grids, "<<eta_star<<",\n";
	datfile<<"number of groups in grid, ";
	for (k=0; k<eta_star; k++) datfile<<Ng[k]<<" , ";
	datfile<<endl;
	
	datfile<<" -- Cell Averaged Scalar Flux -- \n";
	write_group_average_dat(phiT, 0, outw, datfile); // call function to write out cell average scalar flux
	
	datfile<<" -- X Vertical Cell Edge Scalar Flux -- \n";
	write_group_edge_x_dat(phi_xT, 0, outw, datfile); // call function to write out cell edge scalar flux on x grid
	
	datfile<<" -- Y Horizontal Cell Edge Scalar Flux -- \n";
	write_group_edge_y_dat(phi_yT, 0, outw, datfile); // call function to write out cell edge scalar flux on y grid
	
	datfile<<" -- X Face Average Normal Current J_x -- \n";
	write_group_edge_x_dat(j_xT, 0, outw, datfile); // call function to write out cell edge current on x grid
	
	datfile<<" -- Y Face Average Normal Current J_y -- \n";
	write_group_edge_y_dat(j_yT, 0, outw, datfile); // call function to write out cell edge current on y grid
	
	// Write D terms
	datfile<<" ---------------------------- \n";
	datfile<<" -- Diffusion Coefficients -- \n";
	datfile<<" ---------------------------- \n";
	
	datfile<<" -- Cell Average Diffusion Coefficient -- \n";
	write_group_average_dat(D, 0, outw, datfile);
	
	datfile<<" -- X Grid Cell-edge Diffusion Coefficients -- \n";
	write_group_edge_x_dat(D_x, 0, outw, datfile);
	
	datfile<<" -- Y Grid Cell-edge Diffusion Coefficients -- \n";
	write_group_edge_y_dat(D_y, 0, outw, datfile);
	
	datfile<<" ----------------------- \n";
	datfile<<" -- Consistency Terms -- \n"; // 
	datfile<<" ----------------------- \n";
	
	datfile<<" -- X Grid Consistency Term D_xTilde -- \n";
	write_group_edge_x_dat(D_xT, 0, outw, datfile);
	
	datfile<<" -- Y Grid Consistency Term D_yTilde -- \n";
	write_group_edge_y_dat(D_yT, 0, outw, datfile);
	
	datfile.close(); // close output data file ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	
	if ( write_xs ) {
		outf[0]='\0';
		strcat (outf,fn); strcat (outf,".xs.csv"); // name output file
		cout<<outf<<endl;
		datfile.open(outf); // open output data file ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		datfile<<" -------------------- \n";
		datfile<<" -- Cross-Sections -- \n"; // Write NDA solution
		datfile<<" -------------------- \n";
		datfile<<" -- Total XS --\n";
		write_grid_average_dat(sigmaT, outw, datfile);
		
		datfile<<" -- NuF x Fission XS --\n";
		write_grid_average_dat(nuSigmaF, outw, datfile);
		
		datfile<<" -- chi --\n";
		write_grid_average_dat(chi, outw, datfile);
		
		datfile<<" -- Scattering XS --\n";
		for (k=0; k<eta_star; k++) {
			datfile<<"Energy Grid # "<<k<<endl;
			for (gg=0; gg<Ng[k]; gg++) {
				for (g=0; g<Ng[k]; g++) {
					datfile<<" scattering from group "<<gg+1<<" to group "<<g+1<<endl;
					write_cell_average_dat(sigmaS[k][gg][g], outw, datfile);
				}
			}
		}
		
		datfile<<" -- External Source --\n";
		write_grid_average_dat(s_ext, outw, datfile);
		datfile.close(); // close output data file ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	}
	
}
//======================================================================================//

