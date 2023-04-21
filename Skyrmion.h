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

#ifndef Skyrmion_
#define Skyrmion_
class SKYRMION{


public:
/*	SKYRMION(int &Lx, int &Ly)
	    :Lx(Lx), Ly(Ly)
	{
	}
*/	

	int Diameter;
	double Radius_x, Radius_y; //Radius of skyrmion
	double Beta; //decay coeffcient
		
	int No_of_Sk_x, No_of_Sk_y;
	int Lx, Ly;	
	 
	Mat_2_doub Theta_, Phi_;
	double Spin_Size;

	
	void Initialize_Skyrmion();	
	void Create_Skyrmion();
	void Print_Skyrmion(string SkyrmionFilepath);
// Need to add skyrmion number routine
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

}

void SKYRMION::Create_Skyrmion(){


double offset=0.0001;
double sk_cent_x, sk_cent_y;
double sk_cent_x_temp, sk_cent_y_temp;
double dis;

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

dis = Distance(ix*1.0, iy*1.0, sk_cent_x, sk_cent_y);

Theta_[ix][iy] += 2.0*atan(Radius_x/dis)*exp(Beta*(-1.0*dis))*Theta((Radius_x) -dis);



if(ix>sk_cent_x ){
Phi_[ix][iy] += (0*PI_ + atan(  (iy-sk_cent_y+offset)/(ix-sk_cent_x+offset)  ) ) * Theta((Radius_x) -dis);
}
else{
Phi_[ix][iy] += (0*PI_+ atan(  (iy-sk_cent_y+offset)/(ix-sk_cent_x+offset)  ) ) * Theta((Radius_x) -dis);
}



/*
if(iy<sk_cent_y){
Phi_[ix][iy] -= acos((ix-sk_cent_x+offset)/(dis+offset))*Theta((Radius_x) -dis);
}
else{
Phi_[ix][iy] += acos((ix-sk_cent_x+offset)/(dis+offset))*Theta((Radius_x) -dis);
}
*/

if(Theta((Radius_x) -dis)!=0.0){
cout<<ix<<"   "<<iy<<"  "<<sk_ix<<"  "<<sk_iy<<"  "<<sk_cent_x<<"   "<<sk_cent_y<<endl;
}

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

}


#endif
