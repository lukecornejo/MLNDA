#include <iostream>
#include <fstream>
#include <cmath>
#include "ioMLNDA2.1.h"
#include <stdlib.h>
#include <vector>

using namespace std;

// Transport values
extern bool o_angular; // option variables
extern int N, sn;
extern double *mu, *eta, *xi, *w;
extern double **bcL, **bcR, **bcB, **bcT;

extern double ****psiL, ****psiR, ****psiB, ****psiT;

double **psi , **psi_x, **psi_y; // angular flux

// Other variables
const double pi=3.141592654;
extern int kbc;
extern int nx, ny, *Ng;
extern double k_eff;
extern double *x, *y, *hx, *hy, *xe, *ye; // grid arrays
extern double ****phi; // NDA scalar flux and current solution
extern double ***phiT, ***phi_xT, ***phi_yT, ***j_xT, ***j_yT, ***phiL; // scalar flux and current from transport
extern double ***phiB_L, ***phiB_R, ***phiB_B, ***phiB_T; // edge scalar flux from NDA
extern double ****sigmaT, *****sigmaS, ****nuSigmaF, ****chi, ****s_ext; // material arrays
extern double ***FL, ***FR, ***FB, ***FT;
extern double ***D, ***D_x, ***D_y; // Diffusion coefficients
extern double ***D_xT, ***D_yT; // Tilde
extern double ****D_xP, ****D_xN, ****D_yP, ****D_yN; // Tilde, Positive and Negative
extern ofstream temfile;

//======================================================================================//
//++ function to find the quadrature +++++++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
int quadSet()
{
	double mup[N/2], c, wt[8], test=0.0;
	int i, j, m, g, nw, nt;
	int wi[N*(N+2)/8];
	double *wp[N*(N+2)/8];
	
	nt=N*(N+2)/8;
	
	if ( N<20 ) {
		switch ( N ) {
		case 2:
			nw=1;
			wi[0]=1;
			break;
		case 4: // S4 quadrature
			nw=1;
			wi[0]=1; wi[1]=1; 
				 wi[2]=1;
			mup[0]=0.3500212;
			wt[0]=0.333333333;
			break;
		case 6: // S6 quadrature
			nw=2;
			wi[0]=1; wi[1]=2; wi[2]=1; 
				 wi[3]=2; wi[4]=2; 
					wi[5]=1;
			mup[0]=0.2666355;
			wt[0]=0.1761263; wt[1]=0.1572071;
			break;
		case 8: // S8 quadrature
			nw=3;
			wi[0]=1; wi[1]=2; wi[2]=2; wi[3]=1; 
				  wi[4]=2; wi[5]=3; wi[6]=2; 
					  wi[7]=2; wi[8]=2; 
						  wi[9]=1;
			mup[0]=0.2182179;
			wt[0]=0.1209877; wt[1]=0.0907407; wt[2]=0.0925926;
			break;
		case 10:
			nw=4;
			wi[0]=1; wi[1]=2; wi[2]=3; wi[3]=2; wi[4]=1;
				  wi[5]=2; wi[6]=4; wi[7]=4; wi[8]=2; 
					  wi[9]=3; wi[10]=4; wi[11]=3; 
						  wi[12]=2; wi[13]=2; 
							  wi[14]=1;
			break;
		case 12: // S12 quadrature
			nw=5;
			wi[0]=1; wi[1]=2; wi[2]=3; wi[3]=3; wi[4]=2; wi[5]=1;
				 wi[6]=2; wi[7]=4; wi[8]=5; wi[9]=4; wi[10]=2;
					wi[11]=3; wi[12]=5; wi[13]=5; wi[14]=3;
						 wi[15]=3; wi[16]=4; wi[17]=3;
							 wi[18]=2; wi[19]=2;
								 wi[20]=1;
			mup[0]=0.1672126;
			wt[0]=0.0707626; wt[1]=0.0558811; wt[2]=0.0373377; wt[3]=0.0502819; wt[4]=0.0258513;
			break;
		case 14:
			nw=7;
			wi[0]=1; wi[1]=2; wi[2]=3; wi[3]=4; wi[4]=3; wi[5]=2; wi[6]=1;
			   wi[7]=2; wi[8]=5; wi[9]=6; wi[10]=6; wi[11]=5; wi[12]=2;
				   wi[13]=3; wi[14]=6; wi[15]=7; wi[16]=6; wi[17]=3;
					  wi[18]=4; wi[19]=6; wi[20]=6; wi[21]=4;
						   wi[22]=3; wi[23]=5; wi[24]=3;
							   wi[25]=2; wi[26]=2;
								   wi[27]=1;
			break;
		case 16: // S16 quadrature
			nw=8;
			wi[0]=1; wi[1]=2; wi[2]=3; wi[3]=4; wi[4]=4; wi[5]=3; wi[6]=2; wi[7]=1;
			  wi[8]=2; wi[9]=5; wi[10]=6; wi[11]=7; wi[12]=6; wi[13]=5; wi[14]=2;
				  wi[15]=3; wi[16]=6; wi[17]=8; wi[18]=8; wi[19]=6; wi[20]=3;
					 wi[21]=4; wi[22]=7; wi[23]=8; wi[24]=7; wi[25]=4;
						 wi[26]=4; wi[27]=6; wi[28]=6; wi[29]=4;
							 wi[30]=3; wi[31]=5; wi[32]=3;
								 wi[33]=2; wi[34]=2;
									 wi[35]=1;
			mup[0]=0.1389568;
			wt[0]=0.0489872; wt[1]=0.0413296; wt[2]=0.0212326; wt[3]=0.0256207; wt[4]=0.0360486; wt[5]=0.0144589; wt[6]=0.0344958; wt[7]=0.0085179;
			break;
		
		default:
			cout<<"invalid quadrature\n";
			break;
		}
		
		// find weights in quadrent
		for (i=0; i<nt; i++) wp[i]=wt+wi[i]-1; // pointers to angle weights
		
		c=2.0*(1-3*pow(mup[0],2))/(N-2);
		for (i=1; i<N/2; i++) mup[i]=sqrt(pow(mup[i-1],2)+c); // find cosines
		
		// quaderature
		sn=N*(N+2)/8; // total number of number combinations
		mu=new double[sn];
		eta=new double[sn];
		xi=new double[sn];
		w=new double[sn];
		
		// first quadrant
		m=0;
		for (i=N/2; i>0; i--) { // first quadrant
			for (j=0; j<i; j++) {
				mu[m] = mup[i-j-1]; // cosine on x  >
				eta[m]= mup[j];     // cosine on y  <
				xi[m] = mup[N/2-i]; // cosine on z  <
				w[m]=*wp[m]*pi;
				m++;
			}
		}
		
	}
	else {
		sn=N;
		mu=new double[sn];
		eta=new double[sn];
		xi=new double[sn];
		w=new double[sn];
		
		switch ( N ) {
		case 36: // Q36=q461214 quadrature
			
			w[ 0]=8.454511187252e-03; mu[ 0]=9.717784813336e-01; eta[ 0]=1.096881837272e-02; xi[ 0]=2.356401244281e-01;
			w[ 1]=1.913728513580e-02; mu[ 1]=9.701698603928e-01; eta[ 1]=5.695764868253e-02; xi[ 1]=2.356401244312e-01;
			w[ 2]=2.863542971348e-02; mu[ 2]=9.622473153642e-01; eta[ 2]=1.362124657777e-01; xi[ 2]=2.356401244295e-01;
			w[ 3]=3.648716160597e-02; mu[ 3]=9.410672772109e-01; eta[ 3]=2.426233944222e-01; xi[ 3]=2.356401244311e-01;
			w[ 4]=4.244873302980e-02; mu[ 4]=8.997294996538e-01; eta[ 4]=3.673697853806e-01; xi[ 4]=2.356401244316e-01;
			w[ 5]=4.642823955812e-02; mu[ 5]=8.335743322378e-01; eta[ 5]=4.996274255819e-01; xi[ 5]=2.356401244286e-01;
			w[ 6]=4.841339013884e-02; mu[ 6]=7.417637460141e-01; eta[ 6]=6.279014865859e-01; xi[ 6]=2.356401244320e-01;
			w[ 7]=4.841339013884e-02; mu[ 7]=6.279014865859e-01; eta[ 7]=7.417637460141e-01; xi[ 7]=2.356401244320e-01;
			w[ 8]=4.642823955812e-02; mu[ 8]=4.996274255819e-01; eta[ 8]=8.335743322378e-01; xi[ 8]=2.356401244286e-01;
			w[ 9]=4.244873302980e-02; mu[ 9]=3.673697853806e-01; eta[ 9]=8.997294996538e-01; xi[ 9]=2.356401244316e-01;
			w[10]=3.648716160597e-02; mu[10]=2.426233944222e-01; eta[10]=9.410672772109e-01; xi[10]=2.356401244311e-01;
			w[11]=2.863542971348e-02; mu[11]=1.362124657777e-01; eta[11]=9.622473153642e-01; xi[11]=2.356401244295e-01;
			w[12]=1.913728513580e-02; mu[12]=5.695764868253e-02; eta[12]=9.701698603928e-01; xi[12]=2.356401244312e-01;
			w[13]=8.454511187252e-03; mu[13]=1.096881837272e-02; eta[13]=9.717784813336e-01; xi[13]=2.356401244281e-01;
			w[14]=8.352354145856e-03; mu[14]=7.656319455497e-01; eta[14]=1.160393058611e-02; xi[14]=6.431742164832e-01;
			w[15]=1.873220073879e-02; mu[15]=7.633693960835e-01; eta[15]=5.995074957044e-02; xi[15]=6.431742164834e-01;
			w[16]=2.759429759588e-02; mu[16]=7.524467626583e-01; eta[16]=1.419535016004e-01; xi[16]=6.431742164829e-01;
			w[17]=3.442681426024e-02; mu[17]=7.241384940891e-01; eta[17]=2.488983098158e-01; xi[17]=6.431742164835e-01;
			w[18]=3.901232700510e-02; mu[18]=6.711819639118e-01; eta[18]=3.685670882907e-01; xi[18]=6.431742164829e-01;
			w[19]=4.130171453748e-02; mu[19]=5.909368760506e-01; eta[19]=4.869502395267e-01; xi[19]=6.431742164829e-01;
			w[20]=4.130171453748e-02; mu[20]=4.869502395267e-01; eta[20]=5.909368760506e-01; xi[20]=6.431742164829e-01;
			w[21]=3.901232700510e-02; mu[21]=3.685670882907e-01; eta[21]=6.711819639118e-01; xi[21]=6.431742164829e-01;
			w[22]=3.442681426024e-02; mu[22]=2.488983098158e-01; eta[22]=7.241384940891e-01; xi[22]=6.431742164835e-01;
			w[23]=2.759429759588e-02; mu[23]=1.419535016004e-01; eta[23]=7.524467626583e-01; xi[23]=6.431742164829e-01;
			w[24]=1.873220073879e-02; mu[24]=5.995074957044e-02; eta[24]=7.633693960835e-01; xi[24]=6.431742164834e-01;
			w[25]=8.352354145856e-03; mu[25]=1.160393058611e-02; eta[25]=7.656319455497e-01; xi[25]=6.431742164832e-01;
			w[26]=1.460888798152e-02; mu[26]=4.445439440056e-01; eta[26]=2.447911451942e-02; xi[26]=8.954225007226e-01;
			w[27]=2.995376809966e-02; mu[27]=4.288508824476e-01; eta[27]=1.196054590036e-01; xi[27]=8.954225007227e-01;
			w[28]=3.798783310581e-02; mu[28]=3.670788892962e-01; eta[28]=2.519357740235e-01; xi[28]=8.954225007226e-01;
			w[29]=3.798783310581e-02; mu[29]=2.519357740235e-01; eta[29]=3.670788892962e-01; xi[29]=8.954225007226e-01;
			w[30]=2.995376809966e-02; mu[30]=1.196054590036e-01; eta[30]=4.288508824476e-01; xi[30]=8.954225007227e-01;
			w[31]=1.460888798152e-02; mu[31]=2.447911451942e-02; eta[31]=4.445439440056e-01; xi[31]=8.954225007226e-01;
			w[32]=6.404244616724e-03; mu[32]=1.483114568272e-01; eta[32]=1.670387759191e-02; xi[32]=9.887996218887e-01;
			w[33]=1.162080754372e-02; mu[33]=1.293388490485e-01; eta[33]=7.447663982495e-02; xi[33]=9.887996218887e-01;
			w[34]=1.162080754372e-02; mu[34]=7.447663982495e-02; eta[34]=1.293388490485e-01; xi[34]=9.887996218887e-01;
			w[35]=6.404244616724e-03; mu[35]=1.670387759191e-02; eta[35]=1.483114568272e-01; xi[35]=9.887996218887e-01;
			
			break;
		case 20: // Q20=q2468 quadrature
			
			w[ 0]=2.419260514149E-02; mu[ 0]=9.713274064903E-01; eta[ 0]=3.157215799340E-02; xi[ 0]=2.356401244281E-01;
			w[ 1]=5.213067212540E-02; mu[ 1]=9.586898685237E-01; eta[ 1]=1.593344524838E-01; xi[ 1]=2.356401244307E-01;
			w[ 2]=7.185542471164E-02; mu[ 2]=9.028558915298E-01; eta[ 2]=3.596178122512E-01; xi[ 2]=2.356401244304E-01;
			w[ 3]=8.182604839076E-02; mu[ 3]=7.770210099715E-01; eta[ 3]=5.837054752370E-01; xi[ 3]=2.356401244296E-01;
			w[ 4]=8.182604839076E-02; mu[ 4]=5.837054752370E-01; eta[ 4]=7.770210099715E-01; xi[ 4]=2.356401244296E-01;
			w[ 5]=7.185542471164E-02; mu[ 5]=3.596178122512E-01; eta[ 5]=9.028558915298E-01; xi[ 5]=2.356401244304E-01;
			w[ 6]=5.213067212540E-02; mu[ 6]=1.593344524838E-01; eta[ 6]=9.586898685237E-01; xi[ 6]=2.356401244307E-01;
			w[ 7]=2.419260514149E-02; mu[ 7]=3.157215799340E-02; eta[ 7]=9.713274064903E-01; xi[ 7]=2.356401244281E-01;
			w[ 8]=2.998205782366E-02; mu[ 8]=7.645615896150E-01; eta[ 8]=4.210110375297E-02; xi[ 8]=6.431742164827E-01;
			w[ 9]=6.147460425028E-02; mu[ 9]=7.375714298063E-01; eta[ 9]=2.057068622698E-01; xi[ 9]=6.431742164831E-01;
			w[10]=7.796304620960E-02; mu[10]=6.313311043797E-01; eta[10]=4.332989313333E-01; xi[10]=6.431742164827E-01;
			w[11]=7.796304620960E-02; mu[11]=4.332989313333E-01; eta[11]=6.313311043797E-01; xi[11]=6.431742164827E-01;
			w[12]=6.147460425028E-02; mu[12]=2.057068622698E-01; eta[12]=7.375714298063E-01; xi[12]=6.431742164831E-01;
			w[13]=2.998205782366E-02; mu[13]=4.210110375297E-02; eta[13]=7.645615896150E-01; xi[13]=6.431742164827E-01;
			w[14]=2.932993043666E-02; mu[14]=4.424202396002E-01; eta[14]=4.982847370367E-02; xi[14]=8.954225007227E-01;
			w[15]=5.322055875020E-02; mu[15]=3.858240341629E-01; eta[15]=2.221674140412E-01; xi[15]=8.954225007227E-01;
			w[16]=5.322055875020E-02; mu[16]=2.221674140412E-01; eta[16]=3.858240341629E-01; xi[16]=8.954225007227E-01;
			w[17]=2.932993043666E-02; mu[17]=4.982847370367E-02; eta[17]=4.424202396002E-01; xi[17]=8.954225007227E-01;
			w[18]=1.802505216045E-02; mu[18]=1.409476441875E-01; eta[18]=4.908227124734E-02; xi[18]=9.887996218887E-01;
			w[19]=1.802505216045E-02; mu[19]=4.908227124734E-02; eta[19]=1.409476441875E-01; xi[19]=9.887996218887E-01;
			
			break;
		default:
			cout<<"invalid quadrature\n";
			break;
		}
		
		// first quadrant
		for (m=0; m<N; m++) { // first quadrant
			w[m]=w[m]*pi;
		}
		
	}
	
	for (m=0; m<sn; m++) test+=w[m];
	cout<<"weight test "<<test<<endl;
	// normalize weight
	for (m=0; m<sn; m++) w[m]*=pi/test;
	
	// define SC transport variable memory
	
	// angular flux
	psi=new double*[nx];
	for (i=0; i<nx; i++) psi[i]=new double[ny];
	psi_x=new double*[nx+1];
	for (i=0; i<nx+1; i++) psi_x[i]=new double[ny];
	psi_y=new double*[nx];
	for (i=0; i<nx; i++) psi_y[i]=new double[ny+1];
	
	// initialize boundary conditions
	psiB=new double***[Ng[0]]; 
	psiT=new double***[Ng[0]];
	psiL=new double***[Ng[0]];
	psiR=new double***[Ng[0]];
	for (g=0; g<Ng[0]; g++) {
		psiB[g]=new double**[nx]; 
		psiT[g]=new double**[nx];
		psiL[g]=new double**[ny];
		psiR[g]=new double**[ny];
		for (i=0; i<nx; i++) {
			psiB[g][i]=new double*[sn];
			psiT[g][i]=new double*[sn];
			for (m=0; m<sn; m++) {
				psiB[g][i][m]=new double[5];
				psiT[g][i][m]=new double[5];
				for (int p=0; p<5; p++) {
					psiB[g][i][m][p]=1.0;
					psiT[g][i][m][p]=1.0;
				}
			}
		}
		for (j=0; j<ny; j++) {
			psiL[g][j]=new double*[sn];
			psiR[g][j]=new double*[sn];
			for (m=0; m<sn; m++) {
				psiL[g][j][m]=new double[5];
				psiR[g][j][m]=new double[5];
				for (int p=0; p<5; p++) {
					psiL[g][j][m][p]=1.0;
					psiR[g][j][m][p]=1.0;
				}
			}
		}
	}
	
	// cell average flux
	phiT=new double**[Ng[0]];
	for (g=0; g<Ng[0]; g++) {
		phiT[g]=new double*[nx];
		for (i=0; i<nx; i++) phiT[g][i]=new double[ny];
	}
	
		// cell edge values on x grid
	phi_xT=new double**[Ng[0]];
	j_xT  =new double**[Ng[0]];
	for (g=0; g<Ng[0]; g++) {
		phi_xT[g]=new double*[nx+1];
		j_xT[g]  =new double*[nx+1];
		for (i=0; i<nx+1; i++) {
			phi_xT[g][i]=new double[ny];
			j_xT[g][i]  =new double[ny];
		}
	}
	
	// cell edge flux on y grid
	phi_yT=new double**[Ng[0]];
	j_yT  =new double**[Ng[0]];
	for (g=0; g<Ng[0]; g++) {
		phi_yT[g]=new double*[nx];
		j_yT[g]  =new double*[nx];
		for (i=0; i<nx; i++) {
			phi_yT[g][i]=new double[ny+1];
			j_yT[g][i]  =new double[ny+1];
		}
	}
	return 0;
}
//======================================================================================//

//======================================================================================//
//++ function to solve transport in a single general cell ++++++++++++++++++++++++++++++//
//======================================================================================//
void cellSolution(double psiInB, double psiInL, double SA, double sigma, double mut, double etat, double xit, 
double LT, double LR, double& psiOutT, double& psiOutR, double& psiA  )
{
	double epsilon, exp_epsilon, epsilon_2, mup, muc, mu, du;
	double psiOut1, psiOut2, psiOut3, psiA1, psiA2, psiA3;
	double A1, A2, Lout1, Lout2, Lout3;
	
	mup=sqrt(1-xit*xit);             // mu'
	muc=LT/sqrt(LT*LT+LR*LR);        // mu of cell
	mu=mut/sqrt(mut*mut+etat*etat); // projection onto x-y plane
	
	if ( abs(mu-muc)<1e-15 ) { // ray passes through both corners of cell
		du=sqrt(LT*LT+LR*LR);
		if ( sigma<1e-10 ) {
			// triangle A
			psiOutT=psiInL+SA*du/mup/2.0; // find out going angular flux
			psiA1  =psiInL+SA*du/mup/3; // find cell angular flux
			
			// triangle C
			psiOutR=psiInB+SA*du/mup/2.0; // find out going angular flux
			psiA3  =psiInB+SA*du/mup/3; // find cell angular flux
			
			psiA=0.5*(psiInL+psiInB)+SA*du/mup/3.0;
		}
		else {
			epsilon=sigma*du/mup; // optical thickness
			exp_epsilon=exp(-epsilon); // exponent of epsilon
			epsilon_2=epsilon*epsilon; // square of epsilon
			// triangle A
			psiOutT=(psiInL*(1-exp_epsilon)+SA*(epsilon+exp_epsilon-1)/sigma)/epsilon; // find out going angular flux
			psiA1=2*(psiInL*(epsilon+exp_epsilon-1)+SA*(1+0.5*epsilon_2-exp_epsilon-epsilon)/sigma)/epsilon_2; // find cell angular flux
			
			// triangle C
			psiOutR=(psiInB*(1-exp_epsilon)+SA*(epsilon+exp_epsilon-1)/sigma)/epsilon; // find out going angular flux
			psiA3=2*(psiInB*(epsilon+exp_epsilon-1)+SA*(1+0.5*epsilon_2-exp_epsilon-epsilon)/sigma)/epsilon_2; // find cell angular flux
			
			psiA=((psiInL+psiInB)*(epsilon+exp_epsilon-1.0)+2.0*SA*(1.0+0.5*epsilon_2-exp_epsilon-epsilon)/sigma)/epsilon_2;
		}
		
		// total
		//psiOutT=psiOut1;
		//psiOutR=psiOut3;
		//psiA=0.5*(psiA1+psiA3);
		
		//cout<<"balance 1"<<endl;
		//cout<<"T1 balance "<<SA*du/mup-2*(psiOutT-psiInL)-epsilon*psiA1<<endl; // triangle 1 balance
		//cout<<"T3 balance "<<SA*du/mup-2*(psiOutR-psiInB)-epsilon*psiA3<<endl; // triangle 3 balance
		//cout<<"cell angular balance "<<SA*LT*LR-mut*LR*(psiOutR-psiInL)-etat*LT*(psiOutT-psiInB)-LT*LR*sigma*psiA<<endl;
		//cout<<"cell ang bal "<<SA*LT*LR-LR*(psiOutR-psiInL)-LT*(psiOutT-psiInB)-LT*LR*sigma*psiA<<endl;
		//cout<<SA*du-mup*(psiOutR-psiInL)-mup*(psiOutT-psiInB)-sigma*du*psiA<<endl;
	}
	else if ( mu<muc ) { // ray splits the top and bottom of the cell
		Lout1=mu*LR/sqrt(1.0-mu*mu);
		du=Lout1/mu;
		A1=Lout1*LR/2.0; // Triangle 1 Area
		Lout2=LT-Lout1;
		A2=LT*LR-2*A1; // Parallelogram 2 Area
		if ( sigma<1e-10 ) {
			// triangle A
			psiOut1=psiInL+SA*du/mup/2; // find out going angular flux
			psiA1  =psiInL+SA*du/mup/3; // find cell angular flux
			
			// parallelogram B
			psiOut2=psiInB+SA*du/mup;   // find out going angular flux
			psiA2  =psiInB+SA*du/mup/2; // find cell angular flux
			
			// triangle C
			psiOutR=psiA2;              // find out going angular flux
			psiA3  =psiInB+SA*du/mup/3; // find cell angular flux
		}
		else {
			epsilon=sigma*du/mup; // optical thickness
			exp_epsilon=exp(-epsilon); // exponent of epsilon
			epsilon_2=epsilon*epsilon; // square of epsilon
			// triangle A
			psiOut1=(psiInL*(1-exp_epsilon)+SA*(epsilon+exp_epsilon-1)/sigma)/epsilon; // find out going angular flux
			psiA1=2*(psiInL*(epsilon+exp_epsilon-1)+SA*(1+0.5*epsilon_2-exp_epsilon-epsilon)/sigma)/epsilon_2; // find cell angular flux
			
			// parallelogram B
			psiOut2=psiInB*exp_epsilon+SA*(1-exp_epsilon)/sigma;                       // find out going angular flux
			psiA2  =(psiInB*(1-exp_epsilon)+SA*(epsilon+exp_epsilon-1)/sigma)/epsilon; // find cell angular flux
			
			// triangle C
			psiOutR=psiA2;                                                                                       // find out going angular flux
			psiA3  =2*(psiInB*(epsilon+exp_epsilon-1)+SA*(1+0.5*epsilon_2-exp_epsilon-epsilon)/sigma)/epsilon_2; // find cell angular flux
			
		}
		
		// total
		psiOutT=(Lout1*psiOut1+Lout2*psiOut2)/LT;
		psiA=(A1*psiA1+A2*psiA2+A1*psiA3)/(LT*LR);
		
		//cout<<"balance 2"<<endl;
		//cout<<"T1 balance "<<SA*du/mup-2*(psiOut1-psiInL)-epsilon*psiA1<<endl; // triangle 1 balance
		//cout<<"P2 balance "<<SA*du/mup-(psiOut2-psiInB)-epsilon*psiA2<<endl; // parallelogram 2 balance
		//cout<<"T3 balance "<<SA*du/mup-2*(psiOutR-psiInB)-epsilon*psiA3<<endl; // triangle 3 balance
		//cout<<"cell angular balance "<<SA*LT*LR-mut*LR*(psiOutR-psiInL)-etat*LT*(psiOutT-psiInB)-LT*LR*sigma*psiA<<endl;
	}
	else { // ray splits the right and left side of the cell
		du=LT/mu;
		A1=sqrt(du*du-LT*LT)*LT/2.0; // Triangle 1 Area
		Lout3=sqrt(du*du-LT*LT); // Triangle 3 Length
		Lout2=LR-Lout3;
		A2=LT*LR-2*A1; // Parallelogram 2 Area
		
		if ( sigma<1e-10 ) {
			// triangle A
			psiOutT=psiInL+SA*du/mup/2; // find out going angular flux
			psiA1  =psiInL+SA*du/mup/3; // find cell angular flux
			
			// parallelogram B
			psiOut2=psiInL+SA*du/mup; // find out going angular flux
			psiA2  =psiOutT;          // find cell angular flux
			
			// triangle C
			psiOut3=psiInB+SA*du/mup/2; // find out going angular flux
			psiA3  =psiInB+SA*du/mup/2; // find cell angular flux
			
		}
		else {
			epsilon=sigma*du/mup; // optical thickness
			exp_epsilon=exp(-epsilon); // exponent of epsilon
			epsilon_2=epsilon*epsilon; // square of epsilon
			// triangle A
			psiOutT=(psiInL*(1-exp_epsilon)+SA*(epsilon+exp_epsilon-1)/sigma)/epsilon; // find out going angular flux
			psiA1=2*(psiInL*(epsilon+exp_epsilon-1)+SA*(1+0.5*epsilon_2-exp_epsilon-epsilon)/sigma)/epsilon_2; // find cell angular flux
			
			// parallelogram B
			psiOut2=psiInL*exp_epsilon+SA*(1-exp_epsilon)/sigma; // find out going angular flux
			psiA2=psiOutT;                                       // find cell angular flux
			
			// triangle C
			psiOut3=(psiInB*(1-exp_epsilon)+SA*(epsilon+exp_epsilon-1)/sigma)/epsilon; // find out going angular flux
			psiA3=2*(psiInB*(epsilon+exp_epsilon-1)+SA*(1+0.5*epsilon_2-exp_epsilon-epsilon)/sigma)/epsilon_2; // find cell angular flux
		}
		
		// total
		psiOutR=(Lout3*psiOut3+Lout2*psiOut2)/LR;
		psiA=(A1*psiA1+A2*psiA2+A1*psiA3)/(LT*LR);
		
		//cout<<"balance 3"<<endl;
		//cout<<"T1 balance "<<SA*du/mup-2*(psiOutT-psiInL)-epsilon*psiA1<<endl; // triangle 1 balance
		//cout<<"P2 balance "<<SA*du/mup-(psiOut2-psiInL)-epsilon*psiA2<<endl; // parallelogram 2 balance
		//cout<<"T3 balance "<<SA*du/mup-2*(psiOut3-psiInB)-epsilon*psiA3<<endl; // triangle 3 balance
		//cout<<"cell angular balance "<<SA*LT*LR-mut*LR*(psiOutR-psiInL)-etat*LT*(psiOutT-psiInB)-LT*LR*sigma*psiA<<endl;
	}
	//cout<<"cell angular balance "<<SA*LT*LR-mut*LR*(psiOutR-psiInL)-etat*LT*(psiOutT-psiInB)-LT*LR*sigma*psiA<<endl;
}
//======================================================================================//

//======================================================================================//
//++ sweep through angles and cells in each angular quadrant +++++++++++++++++++++++++++//
//======================================================================================//
void quad1(int g, vector< vector<double> > SA) // solution in quadrant 1
{
	int i, j, m, outw=16;
	double psiA;
	double omega_x, omega_y;
	
	//temfile<<"g "<<g<<" phiT 1\n";
	//write_group_average_dat(phiT, 0, 16, temfile);
	
	for (m=0; m<sn; m++) { // first quadrant
		omega_x=mu[m]; omega_y=eta[m];
		
		for (i=0; i<nx; i++) psi_y[i][0]=psiB[g][i][m][1]; // Bottom In BC
		for (j=0; j<ny; j++) psi_x[0][j]=psiL[g][j][m][1]; // Left In BC
		
		for (j=0; j<ny; j++) { // bottom to top
			for (i=0; i<nx; i++) { // left to right
				//psiInL=psi_x[i][j]; // incoming angular flux on the left
				//psiInB=psi_y[i][j]; // incoming angular flux on the bottom
				//SA=((sigmaS[i][j] + nuF[i][j]*sigmaF[i][j])*phi[i][j]+s_ext[i][j])/(4*pi); // source in the cell
				cellSolution( psi_y[i][j], psi_x[i][j], SA[i][j], sigmaT[0][g][i][j], mu[m], eta[m], xi[m], hx[i], hy[j], psi_y[i][j+1], psi_x[i+1][j], psiA );
				//psi_x[i+1][j]=psiOutR; // outgoing angular flux on the right
				//psi_y[i][j+1]=psiOutT; // outgoing angular flux on the top
				psi[i][j]=psiA;        // cell average angular flux
				// Calculate cell centre values
				phiT[g][i][j]+=psiA*w[m];
			}
		}
		// Calculate Cell Edge values
		for (i=0; i<nx+1; i++) {
			for (j=0; j<ny; j++) {
				phi_xT[g][i][j]+=psi_x[i][j]*w[m];
				j_xT[g][i][j]  +=omega_x*psi_x[i][j]*w[m];
			}
		}
		for (i=0; i<nx; i++) {
			for (j=0; j<ny+1; j++) {
				phi_yT[g][i][j]+=psi_y[i][j]*w[m];
				j_yT[g][i][j]  +=omega_y*psi_y[i][j]*w[m];
			}
		}
		
		for (i=0; i<nx; i++) psiT[g][i][m][1]=psi_y[i][ny]; // Top Out BC
		for (j=0; j<ny; j++) psiR[g][j][m][1]=psi_x[nx][j]; // Right Out BC
		
		// option to print out angular flux
		if ( o_angular ) {
			temfile<<" Quadrant 1 Direction # ,"<<m<<",\n";
			temfile<<"  Qmega_x ,"<<print_csv(omega_x)<<"  Omega_y ,"<<print_csv(omega_y)<<"  Weight ,"<<print_csv(w[m])<<endl;
			temfile<<" -- Cell Averaged Angular Flux -- \n";
			write_cell_average_dat(psi, outw, temfile); // call function to write out cell average scalar flux
			temfile<<" -- X Vertical Cell Edge Angular Flux -- \n";
			write_cell_edge_x_dat(psi_x, outw, temfile); // call function to write out cell edge scalar flux on x grid
			temfile<<" -- Y Horizontal Cell Edge Angular Flux -- \n";
			write_cell_edge_y_dat(psi_y, outw, temfile); // call function to write out cell edge scalar flux on y grid
		}
	}
	//temfile<<"g "<<g<<" phiT 2\n";
	//write_group_average_dat(phiT, 0, 16, temfile);
}
//======================================================================================//
void quad2(int g, vector< vector<double> > SA) // solution in quadrant 2
{
	int i, j, m, outw=16;
	double psiA;
	double omega_x, omega_y;
	
	for (m=0; m<sn; m++) { // second quadrant
		omega_x=-mu[m]; omega_y=eta[m];
		
		for (i=0; i<nx; i++) psi_y[i][0]=psiB[g][i][m][2]; // Bottom In BC
		for (j=0; j<ny; j++) psi_x[nx][j]=psiR[g][j][m][2]; // Right In BC
		
		for (j=0; j<ny; j++) { // bottom to top
			for (i=nx-1; i>=0; i--) { // right to left
				//psiInL=psi_x[i+1][j]; // incoming angular flux on the left
				//psiInB=psi_y[i][j];   // incoming angular flux on the bottom
				cellSolution( psi_y[i][j], psi_x[i+1][j], SA[i][j], sigmaT[0][g][i][j], mu[m], eta[m], xi[m], hx[i], hy[j], psi_y[i][j+1], psi_x[i][j], psiA );
				//psi_x[i][j]=psiOutR;   // outgoing angular flux on the right
				//psi_y[i][j+1]=psiOutT; // outgoing angular flux on the top
				psi[i][j]=psiA;        // cell average angular flux
				// Calculate cell centre values
				phiT[g][i][j]+=psiA*w[m];
			}
		}
		// Calculate Cell Edge values
		for (i=0; i<nx+1; i++) {
			for (j=0; j<ny; j++) {
				phi_xT[g][i][j]+=psi_x[i][j]*w[m];
				j_xT[g][i][j]  +=omega_x*psi_x[i][j]*w[m];
			}
		}
		for (i=0; i<nx; i++) {
			for (j=0; j<ny+1; j++) {
				phi_yT[g][i][j]+=psi_y[i][j]*w[m];
				j_yT[g][i][j]  +=omega_y*psi_y[i][j]*w[m];
			}
		}
		
		// Out BC
		for (i=0; i<nx; i++) psiT[g][i][m][2]=psi_y[i][ny]; // Top Out BC
		for (j=0; j<ny; j++) psiL[g][j][m][2]=psi_x[0][j]; // Left Out BC
		
		// option to print out angular flux
		if ( o_angular ) {
			temfile<<" Quadrant 2 Direction # ,"<<m<<",\n";
			temfile<<"  Qmega_x ,"<<print_csv(omega_x)<<"  Omega_y ,"<<print_csv(omega_y)<<"  Weight ,"<<print_csv(w[m])<<endl;
			temfile<<" -- Cell Averaged Angular Flux -- \n";
			write_cell_average_dat(psi, outw, temfile); // call function to write out cell average scalar flux
			temfile<<" -- X Vertical Cell Edge Angular Flux -- \n";
			write_cell_edge_x_dat(psi_x, outw, temfile); // call function to write out cell edge scalar flux on x grid
			temfile<<" -- Y Horizontal Cell Edge Angular Flux -- \n";
			write_cell_edge_y_dat(psi_y, outw, temfile); // call function to write out cell edge scalar flux on y grid
		}
		
	}
}
//======================================================================================//
void quad3(int g, vector< vector<double> > SA) // solution in quadrant 3
{
	int i, j, m, outw=16;
	double psiA;
	double omega_x, omega_y;
	
	for (m=0; m<sn; m++) { // third quadrant
		omega_x=-mu[m]; omega_y=-eta[m];
		
		for (i=0; i<nx; i++) psi_y[i][ny]=psiT[g][i][m][3]; // Top In BC
		for (j=0; j<ny; j++) psi_x[nx][j]=psiR[g][j][m][3]; // Right In BC
		
		for (j=ny-1; j>=0; j--) { // top to bottom
			for (i=nx-1; i>=0; i--) { // right to left
				//psiInL=psi_x[i+1][j]; // incoming angular flux on the left
				//psiInB=psi_y[i][j+1]; // incoming angular flux on the bottom
				cellSolution( psi_y[i][j+1], psi_x[i+1][j], SA[i][j], sigmaT[0][g][i][j], mu[m], eta[m], xi[m], hx[i], hy[j], psi_y[i][j], psi_x[i][j], psiA );
				//psi_x[i][j]=psiOutR; // outgoing angular flux on the right
				//psi_y[i][j]=psiOutT; // outgoing angular flux on the top
				psi[i][j]=psiA;      // cell average angular flux
				// Calculate cell centre values
				phiT[g][i][j]+=psiA*w[m];
			}
		}
		// Calculate Cell Edge values
		for (i=0; i<nx+1; i++) {
			for (j=0; j<ny; j++) {
				phi_xT[g][i][j]+=psi_x[i][j]*w[m];
				j_xT[g][i][j]  +=omega_x*psi_x[i][j]*w[m];
			}
		}
		for (i=0; i<nx; i++) {
			for (j=0; j<ny+1; j++) {
				phi_yT[g][i][j]+=psi_y[i][j]*w[m];
				j_yT[g][i][j]  +=omega_y*psi_y[i][j]*w[m];
			}
		}
		
		for (i=0; i<nx; i++) psiB[g][i][m][3]=psi_y[i][0]; // Bottom Out BC
		for (j=0; j<ny; j++) psiL[g][j][m][3]=psi_x[0][j]; // Left Out BC
		
		// option to print out angular flux
		if ( o_angular ) {
			temfile<<" Quadrant 3 Direction # ,"<<m<<",\n";
			temfile<<"  Qmega_x ,"<<print_csv(omega_x)<<"  Omega_y ,"<<print_csv(omega_y)<<"  Weight ,"<<print_csv(w[m])<<endl;
			temfile<<" -- Cell Averaged Angular Flux -- \n";
			write_cell_average_dat(psi, outw, temfile); // call function to write out cell average scalar flux
			temfile<<" -- X Vertical Cell Edge Angular Flux -- \n";
			write_cell_edge_x_dat(psi_x, outw, temfile); // call function to write out cell edge scalar flux on x grid
			temfile<<" -- Y Horizontal Cell Edge Angular Flux -- \n";
			write_cell_edge_y_dat(psi_y, outw, temfile); // call function to write out cell edge scalar flux on y grid
		}
	}
}
//======================================================================================//
void quad4(int g, vector< vector<double> > SA) // solution in quadrant 4
{
	int i, j, m, outw=16;
	double psiA;
	double omega_x, omega_y;
	
	for (m=0; m<sn; m++) { // fourth quadrant
		omega_x=mu[m]; omega_y=-eta[m];
		
		for (i=0; i<nx; i++) psi_y[i][ny]=psiT[g][i][m][4]; // Top In BC
		for (j=0; j<ny; j++) psi_x[0][j]=psiL[g][j][m][4]; // Left In BC
		
		for (j=ny-1; j>=0; j--) { // top to bottom
			for (i=0; i<nx; i++) { // left to right
				//psiInL=psi_x[i][j];   // incoming angular flux on the left
				//psiInB=psi_y[i][j+1]; // incoming angular flux on the bottom
				cellSolution( psi_y[i][j+1], psi_x[i][j], SA[i][j], sigmaT[0][g][i][j], mu[m], eta[m], xi[m], hx[i], hy[j], psi_y[i][j], psi_x[i+1][j], psiA );
				//psi_x[i+1][j]=psiOutR; // outgoing angular flux on the right
				//psi_y[i][j]=psiOutT;   // outgoing angular flux on the top
				psi[i][j]=psiA;        // cell average angular flux
				// Calculate cell centre values
				phiT[g][i][j]+=psiA*w[m];
			}
		}
		// Calculate Cell Edge values
		for (i=0; i<nx+1; i++) {
			for (j=0; j<ny; j++) {
				phi_xT[g][i][j]+=psi_x[i][j]*w[m];
				j_xT[g][i][j]  +=omega_x*psi_x[i][j]*w[m];
			}
		}
		for (i=0; i<nx; i++) {
			for (j=0; j<ny+1; j++) {
				phi_yT[g][i][j]+=psi_y[i][j]*w[m];
				j_yT[g][i][j]  +=omega_y*psi_y[i][j]*w[m];
			}
		}
		
		for (i=0; i<nx; i++) psiB[g][i][m][4]=psi_y[i][0]; // Bottom Out BC
		for (j=0; j<ny; j++) psiR[g][j][m][4]=psi_x[nx][j]; // Right Out BC
		
		// option to print out angular flux
		if ( o_angular ) {
			temfile<<" Quadrant 4 Direction # ,"<<m<<",\n";
			temfile<<"  Qmega_x ,"<<print_csv(omega_x)<<"  Omega_y ,"<<print_csv(omega_y)<<"  Weight ,"<<print_csv(w[m])<<endl;
			temfile<<" -- Cell Averaged Angular Flux -- \n";
			write_cell_average_dat(psi, outw, temfile); // call function to write out cell average scalar flux
			temfile<<" -- X Vertical Cell Edge Angular Flux -- \n";
			write_cell_edge_x_dat(psi_x, outw, temfile); // call function to write out cell edge scalar flux on x grid
			temfile<<" -- Y Horizontal Cell Edge Angular Flux -- \n";
			write_cell_edge_y_dat(psi_y, outw, temfile); // call function to write out cell edge scalar flux on y grid
		}
	}
}
//======================================================================================//

//======================================================================================//
//++ Determine how to sweep thought cells ++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void angleSweep() {
	int g, gg, i, j, m;
	vector< vector<double> > SA, S_f;
	cout<<"Transport Sweep Started : ";
	temfile<<"Transport Solution \n";
	//temfile<<"phiT 1\n";
	//write_group_average_dat(phiT, 0, 16, temfile);
	
	S_f.resize(nx);
	SA.resize(nx);
	for (i=0; i<nx; i++) {
		S_f[i].resize(ny);
		SA[i].resize(ny);
		for (j=0; j<ny; j++) {
			S_f[i][j]=0.0;
			for (g=0; g<Ng[0]; g++) S_f[i][j]+=nuSigmaF[0][g][i][j]*phi[0][g][i][j];
		}
	}
	
	// Zero out Flux and Currents
	for (g=0; g<Ng[0]; g++) {
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) {
				phiT[g][i][j]=0.0;
				SA[i][j]=chi[0][g][i][j]*S_f[i][j]/k_eff+s_ext[0][g][i][j];
				for (gg=0; gg<Ng[0]; gg++) SA[i][j]+=sigmaS[0][gg][g][i][j]*phi[0][gg][i][j];
				SA[i][j]/=4*pi;
			}
		}
		for (i=0; i<nx+1; i++) {
			for (j=0; j<ny; j++) {
				phi_xT[g][i][j]=0.0; j_xT[g][i][j]=0.0;
			}
		}
		for (i=0; i<nx; i++) {
			for (j=0; j<ny+1; j++) {
				phi_yT[g][i][j]=0.0; j_yT[g][i][j]=0.0;
			}
		}
	
		switch ( kbc ) {
		case 1: // incoming boundary conditions on face only 111111111111111111111111111111111111111111111111111111111111111111
			for (m=0; m<sn; m++) {
				for (i=0; i<nx; i++) { // top and bottom boundary conditions
					psiB[g][i][m][1]=bcB[g][i]; // Quad 1 Bottom boundary condition
					psiB[g][i][m][2]=bcB[g][i]; // Quad 2 Bottom boundary condition
					psiT[g][i][m][3]=bcT[g][i]; // Quad 3 Top boundary condition
					psiT[g][i][m][4]=bcT[g][i]; // Quad 4 Top boundary condition
				}
				for (j=0; j<ny; j++) { // left and right boundary conditions
					psiL[g][j][m][1]=bcL[g][j]; // Quad 1 Left boundary condition
					psiL[g][j][m][4]=bcL[g][j]; // Quad 4 Left boundary condition
					psiR[g][j][m][2]=bcR[g][j]; // Quad 2 Right boundary condition
					psiR[g][j][m][3]=bcR[g][j]; // Quad 3 Right boundary condition
				}
			}
			// Start solution sweep///////////////////////////////////////
			quad1(g,SA); // solution sweep through quadrant 1
			quad2(g,SA); // solution sweep through quadrant 2
			quad3(g,SA); // solution sweep through quadrant 3
			quad4(g,SA); // solution sweep through quadrant 4
			break;
		case 2: // incoming boundary conditions in Quadrants 22222222222222222222222222222222222222222222222222222222222222222222
			for (m=0; m<sn; m++) {
				for (i=0; i<nx; i++) { // top and bottom boundary conditions
					psiB[g][i][m][1]=bcL[g][i]; // Quad 1 Bottom boundary condition
					psiB[g][i][m][2]=bcB[g][i]; // Quad 2 Bottom boundary condition
					psiT[g][i][m][3]=bcR[g][i]; // Quad 3 Top    boundary condition
					psiT[g][i][m][4]=bcT[g][i]; // Quad 4 Top    boundary condition
				}
				for (j=0; j<ny; j++) { // left and right boundary conditions
					psiL[g][j][m][1]=bcL[g][j]; // Quad 1 Left  boundary condition
					psiL[g][j][m][4]=bcT[g][j]; // Quad 4 Left  boundary condition
					psiR[g][j][m][2]=bcB[g][j]; // Quad 2 Right boundary condition
					psiR[g][j][m][3]=bcR[g][j]; // Quad 3 Right boundary condition
				}
			}
			
			// Start solution sweep///////////////////////////////////////
			quad1(g,SA); // solution sweep through quadrant 1
			quad2(g,SA); // solution sweep through quadrant 2
			quad3(g,SA); // solution sweep through quadrant 3
			quad4(g,SA); // solution sweep through quadrant 4
			break;
		case 3: // All Reflective BC 333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
			// Start solution sweep
			for (m=0; m<sn; m++) { // quad 1
				for (i=0; i<nx; i++) {
					psiB[g][i][m][1]=psiB[g][i][m][4]; // Quad 1 Bottom Reflective BC
					psiB[g][i][m][2]=psiB[g][i][m][3]; // Quad 2 Bottom Reflective BC
					psiT[g][i][m][3]=psiT[g][i][m][2]; // Quad 3 Top Reflective BC
					psiT[g][i][m][4]=psiT[g][i][m][1]; // Quad 4 Top Reflective BC
				}
				for (j=0; j<ny; j++) {
					psiL[g][j][m][1]=psiL[g][j][m][2]; // Quad 1 Left Reflective BC
					psiR[g][j][m][2]=psiR[g][j][m][1]; // Quad 2 Right Reflective BC
					psiR[g][j][m][3]=psiR[g][j][m][4]; // Quad 3 Right Reflective BC
					psiL[g][j][m][4]=psiL[g][j][m][3]; // Quad 4 Left Reflective BC
				}
			}
			//temfile<<"g "<<g<<" phiT 1\n";
			//write_group_average_dat(phiT, 0, 16, temfile);
			//temfile<<"SA\n";
			//write_cell_average_dat(SA, 16, temfile);
			quad1(g,SA); // solution sweep through quadrant 1
			//temfile<<"g "<<g<<" phiT 2\n";
			//write_group_average_dat(phiT, 0, 16, temfile);
			quad2(g,SA); // solution sweep through quadrant 2
			//temfile<<g<<"phiT 3\n";
			//write_group_average_dat(phiT, 0, 16, temfile);
			quad3(g,SA); // solution sweep through quadrant 3
			//temfile<<g<<"phiT 4\n";
			//write_group_average_dat(phiT, 0, 16, temfile);
			quad4(g,SA); // solution sweep through quadrant 4
			//temfile<<g<<"phiT 5\n";
			//write_group_average_dat(phiT, 0, 16, temfile);
				
			break;
		case 4: //Reflective BC on bottom and left, face BC 444444444444444444444444444444444444444444444444444444444444444444444444
			// Start solution sweep ///////////////////////////////////////////////////////////
			for (m=0; m<sn; m++) { // quad 3
				for (i=0; i<nx; i++) psiT[g][i][m][3]=bcT[g][i]; // Quad 3 Top BC
				for (j=0; j<ny; j++) psiR[g][j][m][3]=bcR[g][j]; // Quad 3 Right BC
			}
			quad3(g,SA); // solution sweep through quadrant 3
			
			for (m=0; m<sn; m++) { // quad 4
				for (i=0; i<nx; i++) psiT[g][i][m][4]=bcT[g][i];           // Quad 4 Top BC
				for (j=0; j<ny; j++) psiL[g][j][m][4]=psiL[g][j][m][3]; // Quad 4 Left Reflective BC
			}
			quad4(g,SA); // solution sweep through quadrant 4
			
			for (m=0; m<sn; m++) { // quad 2
				for (i=0; i<nx; i++) psiB[g][i][m][2]=psiB[g][i][m][3]; // Quad 2 Bottom Reflective BC
				for (j=0; j<ny; j++) psiR[g][j][m][2]=bcR[g][j];           // Quad 2 Right BC
			}
			quad2(g,SA); // solution sweep through quadrant 2
			
			for (m=0; m<sn; m++) { // quad 1
				for (i=0; i<nx; i++) psiB[g][i][m][1]=psiB[g][i][m][4]; // Quad 1 Bottom Reflective BC
				for (j=0; j<ny; j++) psiL[g][j][m][1]=psiL[g][j][m][2]; // Quad 1 Left Reflective BC
			}
			quad1(g,SA); // solution sweep through quadrant 1
			break;
		case 5: // Reflective BC on top and right, face BC 55555555555555555555555555555555555555555555555555555555555555555555555555
			// Start solution sweep ////////////////////////////////////////////////////////
			for (m=0; m<sn; m++) {
				for (i=0; i<nx; i++) psiB[g][i][m][1]=bcB[g][i]; // Quad 1 Bottom boundary condition
				for (j=0; j<ny; j++) psiL[g][j][m][1]=bcL[g][j]; // Quad 1 Left boundary condition
			}
			quad1(g,SA); // solution sweep through quadrant 1
			
			for (m=0; m<sn; m++) {
				for (i=0; i<nx; i++) psiB[g][i][m][2]=bcB[g][i];           // Quad 2 Bottom BC
				for (j=0; j<ny; j++) psiR[g][j][m][2]=psiR[g][j][m][1]; // Quad 2 Right reflective BC
			}
			quad2(g,SA); // solution sweep through quadrant 2
			
			for (m=0; m<sn; m++) { // quad 4
				for (i=0; i<nx; i++) psiT[g][i][m][4]=psiT[g][i][m][1]; // Quad 4 Top reflective BC
				for (j=0; j<ny; j++) psiL[g][j][m][4]=bcL[g][j];           // Quad 4 Left BC
			}
			quad4(g,SA); // solution sweep through quadrant 4
			
			for (m=0; m<sn; m++) { // quad 3
				for (i=0; i<nx; i++) psiT[g][i][m][3]=psiT[g][i][m][2]; // Quad 3 Top reflective BC
				for (j=0; j<ny; j++) psiR[g][j][m][3]=psiR[g][j][m][4]; // Quad 3 Right reflective BC
			}
			quad3(g,SA); // solution sweep through quadrant 3
			break;
		default:
			cout<<"bad boundary conditions: incorrect boundary type"<<endl;
			break;
		}
	}
	//temfile<<"phiT 7\n";
	//write_group_average_dat(phiT, 0, 16, temfile);
	cout<<"Completed \n";
}
//======================================================================================//

//======================================================================================//
//++ Initial Transport Solution Found ++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void initialAngleSweep()
{
	int g, gg, i, j, m, k;
	double psi_const=1.0;
	temfile<<"Initial Angular Flux Guess\n";
	// Find Initial Flux and Current from constant angular flux
	for (g=0; g<Ng[0]; g++) {
		for ( i=0; i<nx; i++ ) {
			for ( j=0; j<ny; j++ ) {
				phiT[g][i][j]=0.0;
				for ( m=0; m<sn; m++ ) phiT[g][i][j]+=psi_const*w[m];
				phiT[g][i][j]*=4.0;
			}
		}
		for ( i=0; i<nx+1; i++ ) { // X Grid
			for ( j=0; j<ny; j++ ) {
				phi_xT[g][i][j]=0.0; j_xT[g][i][j]=0.0;
				for ( m=0; m<sn; m++ ) {
					j_xT[g][i][j]+=(psi_const-psi_const)*w[m]*mu[m];
					phi_xT[g][i][j]+=psi_const*w[m];
				}
				phi_xT[g][i][j]*=4.0; j_xT[g][i][j]*=2.0;
			}
		}
		for ( i=0; i<nx; i++ ) { // Y Grid
			for ( j=0; j<ny+1; j++ ) {
				phi_yT[g][i][j]=0.0; j_yT[g][i][j]=0.0;
				for ( m=0; m<sn; m++ ) {
					j_yT[g][i][j]+=(psi_const-psi_const)*w[m]*mu[m];
					phi_yT[g][i][j]+=psi_const*w[m];
				}
				phi_yT[g][i][j]*=4.0; j_yT[g][i][j]*=2.0;
			}
		}
	}
	
	if (kbc==3) {
		for (g=0; g<Ng[0]; g++) for ( i=0; i<nx; i++ ) for ( j=0; j<ny; j++ ) phiT[g][i][j]+=rand()%5;
	}
}
//======================================================================================//

