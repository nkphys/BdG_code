#include <iostream>
#include <math.h>  //fabs(double x) =|x|
#include <algorithm>
#include <stdlib.h>  //for div(q,n).rem(quot),rand
#include <time.h>
#include <fstream>
#include <limits>
#include <iomanip>
#include <stdio.h>



#include "functions_real.h"
using namespace std;

#ifndef SkyrmionDef_
#define SkyrmionDef_
class SKYRMION{


public:
	SKYRMION(string &Skyrmion_Type_)
	    :Skyrmion_Type(Skyrmion_Type_)
	{
	}
	

	int Diameter;
	double Radius_x, Radius_y; //Radius of skyrmion
	double Beta; //decay coeffcient
		
	double OptimizedBeta;
	int No_of_Sk_x, No_of_Sk_y;
	int Lx, Ly;	
	 
	Mat_2_doub Theta_, Phi_;
	double Spin_Size;

	string BraviasLattice;
		
	double Skyrmion_number;
	double vorticity, helicity, polarity; 
	string Skyrmion_Type;
	
	void Initialize_Skyrmion();	
	void Create_Skyrmion();
	void Print_Skyrmion(string SkyrmionFilepath);
// Need to add skyrmion number routine
	void Skyrmion_Number_Calculate();
};



void SKYRMION::Initialize_Skyrmion(){


//Theta, Phi initialization-----------
Theta_.resize(Lx);Phi_.resize(Lx);
for(int i=0;i<Lx;i++){
Theta_[i].resize(Ly);
Phi_[i].resize(Ly);
}

for(int i=0;i<Lx;i++){
for(int j=0;j<Ly;j++){
Theta_[i][j]=0.0;
Phi_[i][j]=0.0;
}
}

No_of_Sk_x = int(Lx/Diameter);
No_of_Sk_y = int(Ly/Diameter);
//assert((Lx%(Diameter))==0);
//assert((Ly%(Diameter))==0);

Radius_x = (1.0*Diameter)/2.0;
Radius_y = (1.0*Diameter)/2.0;

if(Radius_x!=Radius_y){
cout<<" Radius_x -ne Radius_y  needs NEW EQUATION"<<endl;
assert(false);
}


Spin_Size=1.0;
OptimizedBeta=0.05;

}

void SKYRMION::Create_Skyrmion(){


double offset=0.0001;
double alpha = 1.0; //contribution decay constant used in "Fermi" function
double sk_cent_x, sk_cent_y;
double sk_cent_x_temp, sk_cent_y_temp;
double dis;
double dis_x, dis_y;

/*
vorticity = 1.0 (skyr) ; vorticity = -1.0 (antiskyr)
helicity = pi (inward Neel skyr) ; helicity = 0.0 (outward Neel skyr)
helicity = -pi/2 (anticlock Bloch skyr) ; helicity = pi/2 (Clock Bloch skyr)
helicity = pi/n (n = integer) for twisted skyrmion
polarity = 1.0 (up skyr center spin) ; vorticity = 1.0 (down skyr center spin)
*/

for(int ix=0;ix<Lx;ix++){
for(int iy=0;iy<Ly;iy++){

Theta_[ix][iy]=0;Phi_[ix][iy]=0;

for(int sk_ix=0;sk_ix<No_of_Sk_x;sk_ix++){
for(int sk_iy=0;sk_iy<No_of_Sk_y;sk_iy++){

sk_cent_x = ((  (2*Radius_x + 1.0)*sk_ix) + Radius_x );
sk_cent_y = ((  (2*Radius_y + 1.0)*sk_iy) + Radius_y );


//PBC for sites
//sk_cent_x = sk_cent_x_temp%Lx;
//sk_cent_y = sk_cent_y_temp%Ly;

//dis = Distance(ix*1.0, iy*1.0, sk_cent_x, sk_cent_y);


if(BraviasLattice=="SquareLattice"){
//dis = Distance(ix*1.0, iy*1.0, sk_cent_x, sk_cent_y);
dis_x = abs(ix*1.0 - sk_cent_x);
dis_y = abs(iy*1.0 - sk_cent_y);
}
if(BraviasLattice=="TriangularLattice"){
dis_x = abs((ix*1.0 - sk_cent_x) + ((iy*1.0 - sk_cent_y)*0.5));
dis_y = abs((sqrt(3.0)/2.0)*(iy*1.0 - sk_cent_y));
}
dis = Distance(dis_x, dis_y, 0, 0);



//Theta_[ix][iy] += ( 2.0*atan(Radius_x/dis)*exp(Beta*(-1.0*dis)) + acos(-polarity)) * Theta((Radius_x) -dis);


/*
Theta_[ix][iy] += ( 2.0*atan(Radius_x/dis)*exp(Beta*(-1.0*dis)) + acos(-polarity)) * Fermi(dis, Radius_x, alpha);

//Find reference of following equations
if ( iy<sk_cent_y ) {
Phi_[ix][iy] += (-acos(vorticity*(ix-sk_cent_x)/dis) - helicity) * Fermi(dis, Radius_x, alpha);
//Phi_[ix][iy] += ( -acos(vorticity*(ix-sk_cent_x)/dis) - helicity )*Theta((Radius_x) -dis);
}
else {
Phi_[ix][iy] += ( acos(vorticity*(ix-sk_cent_x)/dis) - helicity ) * Fermi(dis, Radius_x, alpha);
//Phi_[ix][iy] += ( acos(vorticity*(ix-sk_cent_x)/dis) - helicity )*Theta((Radius_x) -dis);
}

if ( ix==sk_cent_x && iy==sk_cent_y ) {
	Phi_[ix][iy] = helicity;
}
*/


if(Skyrmion_Type=="AntiSkyrmion"){
Theta_[ix][iy] += 2.0*atan(Radius_x/dis)*exp(Beta*(-1.0*dis))*Theta((Radius_x) -dis);

if( iy>=sk_cent_y && ix>=sk_cent_x ){
Phi_[ix][iy] += (atan(  (dis_x+offset)/(dis_y+offset)  ) ) * Theta((Radius_x) -dis);
}
else if (iy<sk_cent_y && ix>=sk_cent_x){
Phi_[ix][iy] += (PI_ - atan(  (abs(dis_x)+offset)/(abs(dis_y)+offset)  ) ) * Theta((Radius_x) -dis);
}
else if (iy<=sk_cent_y && ix<sk_cent_x){
Phi_[ix][iy] += (PI_ + atan(  (abs(dis_x)+offset)/(abs(dis_y)+offset)  ) ) * Theta((Radius_x) -dis);
}
else if (iy>sk_cent_y && ix<sk_cent_x) {
Phi_[ix][iy] += (2*PI_ - atan(  (abs(dis_x)+offset)/(abs(dis_y)+offset)  ) ) * Theta((Radius_x) -dis);
}
}


if(Skyrmion_Type=="BlochSkyrmion"){
// add here
}

if(Skyrmion_Type=="NeelSkyrmion"){

Theta_[ix][iy] += 2.0*atan(Radius_x/dis)*exp(Beta*(-1.0*dis))*Theta((Radius_x) -dis);

if( iy>=sk_cent_y && ix>=sk_cent_x ){
Phi_[ix][iy] += (atan(  (dis_y+offset)/(dis_x+offset)  ) ) * Theta((Radius_x) -dis);
}
else if (ix<sk_cent_x && iy>=sk_cent_y){
Phi_[ix][iy] += (PI_ -  atan(  (abs(dis_y)+offset)/(abs(dis_x)+offset)  ) ) * Theta((Radius_x) -dis);
}
else if (iy<sk_cent_y && ix<=sk_cent_x){
Phi_[ix][iy] += (PI_*1.0 + atan(  (abs(dis_y)+offset)/(abs(dis_x)+offset)  ) ) * Theta((Radius_x) -dis);
}
else if (iy<sk_cent_y && ix>sk_cent_x) {
Phi_[ix][iy] += (2*PI_ - atan(  (abs(dis_y)+offset)/(abs(dis_x)+offset)  ) ) * Theta((Radius_x) -dis);
}



//add Here
}

//similarly more Skyrmions if needed ex. bimeron etc.

/*
if(iy<sk_cent_y){
Phi_[ix][iy] -= acos((ix-sk_cent_x+offset)/(dis+offset))*Theta((Radius_x) -dis);
}
else{
Phi_[ix][iy] += acos((ix-sk_cent_x+offset)/(dis+offset))*Theta((Radius_x) -dis);
}
*/

// if(Theta((Radius_x) -dis)!=0.0){
// cout<<ix<<"   "<<iy<<"  "<<sk_ix<<"  "<<sk_iy<<"  "<<sk_cent_x<<"   "<<sk_cent_y<<endl;
// }

}
}

}
}

}



void SKYRMION::Print_Skyrmion(string SkyrmionFilepath){

double Sz, Sx, Sy;

ofstream outfile(SkyrmionFilepath.c_str());

outfile<<"#ix   iy   Sz   Sx   Sy   Theta[ix][iy]     Phi[ix][iy]"<<endl;
for(int ix=0;ix<Lx;ix++){
for(int iy=0;iy<Ly;iy++){

Sz=Spin_Size*cos(Theta_[ix][iy]);
Sx=Spin_Size*sin(Theta_[ix][iy])*cos(Phi_[ix][iy]);
Sy=Spin_Size*sin(Theta_[ix][iy])*sin(Phi_[ix][iy]);

outfile<<ix<<"  "<<iy<<"  "<< Sz<<"  "<<Sx<<"  "<<Sy<<"  "<< Theta_[ix][iy]<<"  "<<Phi_[ix][iy]<< "  "<<cos(Phi_[ix][iy])<<"  "<<sin(Phi_[ix][iy])<<endl;
}

outfile<<endl;

}

outfile<<"#Skyrmion number = "<<Skyrmion_number<<endl; 

}


void SKYRMION::Skyrmion_Number_Calculate(){

double chi1, chi2, chi;
int jx, jy, jx1, jy1, jx2, jy2;

Mat_2_doub sx; sx.resize(Lx);
for (int ix = 0; ix < Lx; ix++) {
		sx[ix].resize(Ly);
		for (int iy = 0; iy < Ly; iy++) {
				sx[ix][iy] = 0.0;
		}
}

double** sy = new double*[Lx];
for (int ix = 0; ix < Lx; ix++) {
		sy[ix] = new double[Ly];
		for (int iy = 0; iy < Ly; iy++) {
				sy[ix][iy] = 0.0;
		}
}

double** sz = new double*[Lx];
for (int ix = 0; ix < Lx; ix++) {
		sz[ix] = new double[Ly];
		for (int iy = 0; iy < Ly; iy++) {
				sz[ix][iy] = 0.0;
		}
}

for (int ix = 0 ; ix < Lx ; ix++){
	for (int iy = 0 ; iy < Ly ; iy++){
		sx[ix][iy] = Spin_Size * sin(Theta_[ix][iy]) * cos(Phi_[ix][iy]);
		sy[ix][iy] = Spin_Size * sin(Theta_[ix][iy]) * sin(Phi_[ix][iy]);
		sz[ix][iy] = Spin_Size * cos(Theta_[ix][iy]);
	}
}

chi = 0.0;
for (int ix = 0 ; ix < Lx ; ix++){
	for (int iy = 0 ; iy < Ly ; iy++){
		chi1 = 0.0 ; chi2 = 0.0;
		jy = iy      ; jx1 = ix + 1; if (jx1 >= Lx) {jx1 = 0;}
		jy = iy      ; jx2 = ix - 1; if (jx2 < 0)  {jx2 = Lx - 1;}
		jy1 = iy + 1 ; jx = ix     ; if (jy1 >= Ly) {jy1 = 0;}
		jy2 = iy - 1 ; jx = ix     ; if (jy2 < 0)  {jy2 = Ly - 1;}

		chi1 = (1/(8.0*PI_)) * (  (sx[ix][iy] * (sy[jx1][jy]*sz[jx][jy1] - sz[jx1][jy]*sy[jx][jy1])) + \
															(sy[ix][iy] * (sz[jx1][jy]*sx[jx][jy1] - sx[jx1][jy]*sz[jx][jy1])) + \
															(sz[ix][iy] * (sx[jx1][jy]*sy[jx][jy1] - sy[jx1][jy]*sx[jx][jy1]))  );

		chi2 = (1/(8.0*PI_)) * (  (sx[ix][iy] * (sy[jx2][jy]*sz[jx][jy2] - sz[jx2][jy]*sy[jx][jy2])) + \
															(sy[ix][iy] * (sz[jx2][jy]*sx[jx][jy2] - sx[jx2][jy]*sz[jx][jy2])) + \
															(sz[ix][iy] * (sx[jx2][jy]*sy[jx][jy2] - sy[jx2][jy]*sx[jx][jy2]))  );

		chi = chi + chi1 + chi2;

	}
}

Skyrmion_number=chi;

}

#endif
