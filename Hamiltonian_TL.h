#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor_type.h"
#include "Parameters_TL.h"
#include "Coordinates_TL.h"
#include "Matrix.h"
#include "Skyrmion.h"
#define PI acos(-1.0)

#ifndef Hamiltonian_TL_class
#define Hamiltonian_TL_class


extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
                         std::complex<double> *,int *, double *, int *);

extern "C" void dgesdd_ (char *, int *, int *, double *, int *, double *, double *, int *, double *, int*,
                         double *, int *, int *, int *);


extern "C" void zgesdd_ (char *, int *, int *, std::complex<double> *, int *, double *, std::complex<double> *, int *, std::complex<double> *, int*,
                         std::complex<double> *, int *, double * , int *, int *);


class Hamiltonian_TL
{
public:
    Hamiltonian_TL(Parameters_TL  &Parameters__ , Coordinates_TL &Coordinates__, SKYRMION &Skyrmion__)
        :Parameters_(Parameters__),Coordinates_(Coordinates__),Skyrmion_(Skyrmion__)

    {
        Initialize();
        HTBCreate();
    }


    void Initialize();
    void HTBCreate();
    void Add_SpinFermionTerm();
    void Add_PairingTerm();
    void Add_ChemicalPotentialTerm();
    void HamilCreation();
    void Diagonalize(char option);



    bool SpinFermionTerm, PairingTerm, ChemicalPotentialTerm;

    Matrix<complex<double>> Pauli_z, Pauli_x, Pauli_y;

    Coordinates_TL &Coordinates_;
    Parameters_TL &Parameters_;
    SKYRMION &Skyrmion_;
    int lx_, ly_, ncells_, n_orbs_;

    Mat_1_doub eigs_;

    Matrix<complex<double>> Ham_, HTB_;



};


void Hamiltonian_TL::Diagonalize(char option){


    char jobz=option;
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3*n -2);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

//     Ham_.print();

}




void Hamiltonian_TL::Add_SpinFermionTerm(){



int col_,row_;
double Sz_i, Sx_i, Sy_i; 

int ix, iy;

for(int site=0;site<ncells_;site++){
ix=Coordinates_.indx_cellwise_[site];
iy=Coordinates_.indy_cellwise_[site];

Sz_i = Skyrmion_.Spin_Size*cos(Skyrmion_.Theta_[ix][iy]);
Sx_i = Skyrmion_.Spin_Size*sin(Skyrmion_.Theta_[ix][iy])*cos(Skyrmion_.Phi_[ix][iy]);
Sx_i = Skyrmion_.Spin_Size*sin(Skyrmion_.Theta_[ix][iy])*sin(Skyrmion_.Phi_[ix][iy]);


for(int spin_beta=0;spin_beta<2;spin_beta++){
for(int spin_alpha=0;spin_alpha<2;spin_alpha++){

col_ = site + ncells_*spin_beta;
row_ = site + ncells_*spin_alpha;


Ham_(row_,col_) += 0.5*0.5*Parameters_.J_Hund*(Pauli_z(spin_alpha,spin_beta)*Sz_i  + 
                   Pauli_x(spin_alpha,spin_beta)*Sx_i +
		   Pauli_y(spin_alpha,spin_beta)*Sy_i) ;


Ham_(row_ + 2*ncells_,col_+2*ncells_) += -0.5*0.5*Parameters_.J_Hund*(Pauli_z(spin_alpha,spin_beta)*Sz_i  +
                   Pauli_x(spin_alpha,spin_beta)*Sx_i +
                   Pauli_y(spin_alpha,spin_beta)*Sy_i) ;




}
}
}


}


void Hamiltonian_TL::Add_PairingTerm(){


int spin_up, spin_down, col_,row_;
spin_up=0;spin_down=1;
for(int site=0;site<ncells_;site++){


col_= site + ncells_*spin_down + 2*ncells_;
row_= site + ncells_*spin_up;
Ham_(row_,col_) += 0.5*Parameters_.Delta_s;


col_= site + ncells_*spin_up + 2*ncells_;
row_= site + ncells_*spin_down;
Ham_(row_,col_) += -0.5*Parameters_.Delta_s;

col_= site + ncells_*spin_up;
row_= site + ncells_*spin_down + 2*ncells_;
Ham_(row_,col_) += 0.5*Parameters_.Delta_s;

col_= site + ncells_*spin_down;
row_= site + ncells_*spin_up + 2*ncells_;
Ham_(row_,col_) += -0.5*Parameters_.Delta_s;

}




}


void Hamiltonian_TL::Add_ChemicalPotentialTerm(){

int a,b;
for(int l=0;l<ncells_;l++){
for (int spin=0;spin<2;spin++){

a = l + ncells_*spin;
b = l + ncells_*spin;

// assert(a!=b);
Ham_(b,a) += -1.0*Parameters_.mu*0.5;
Ham_(b+2*ncells_,a+2*ncells_) += 1.0*Parameters_.mu*0.5;


}
}

}


void Hamiltonian_TL::HamilCreation(){
Ham_=HTB_;


if(SpinFermionTerm){
Add_SpinFermionTerm();
}


if(PairingTerm){
Add_PairingTerm();
}

if(ChemicalPotentialTerm){
Add_ChemicalPotentialTerm();
}



ly_ = Parameters_.ly;
lx_ = Parameters_.lx;
ncells_ = lx_ * ly_;
n_orbs_ = 1;
int space = 2 * ncells_ *2 ;

// Ham_.print();
/*
for (int i=0; i<space; i++){
  for (int j=0; j<space; j++){
    cout<<i<<"    "<<j<<"    "<<Ham_(i,j)<<endl;
  }
}*/


}


void Hamiltonian_TL::Initialize()
{

    ly_ = Parameters_.ly;
    lx_ = Parameters_.lx;
    ncells_ = lx_ * ly_;
    n_orbs_ = 1;
    int space = 2 * ncells_ *2 ;

    HTB_.resize(space, space);
    Ham_.resize(space, space);


   Pauli_z.resize(2,2);
   Pauli_x.resize(2,2);
   Pauli_y.resize(2,2);


   Pauli_z.fill(0);
   Pauli_z(0,0)=1.0;Pauli_z(1,1)=-1.0;

   Pauli_x.fill(0);
   Pauli_x(0,1)=1.0;Pauli_x(1,0)=1.0;

   Pauli_y.fill(0);
   Pauli_y(0,1)=-1.0*iota_comp;Pauli_y(1,0)=1.0*iota_comp;


   SpinFermionTerm=Parameters_.SpinFermionTerm;
   PairingTerm=Parameters_.PairingTerm;
   ChemicalPotentialTerm=Parameters_.ChemicalPotentialTerm;

   
} // ----------

void Hamiltonian_TL::HTBCreate(){

Mat_1_int t_neighs;
Mat_1_doub t_hoppings;


t_neighs.push_back(0);t_neighs.push_back(2);t_neighs.push_back(7);
t_hoppings.push_back(Parameters_.t_nn); t_hoppings.push_back(Parameters_.t_nn);t_hoppings.push_back(Parameters_.t_nn);

int m,a,b;
for(int l=0;l<ncells_;l++){

for(int neigh=0;neigh<t_neighs.size();neigh++){
m=Coordinates_.getneigh(l,t_neighs[neigh]);

if( (!Coordinates_.HIT_X_BC) || Parameters_.PBC_X){
   if( (!Coordinates_.HIT_Y_BC) || Parameters_.PBC_Y){

for (int spin=0;spin<2;spin++){

a = l + ncells_*spin;
b = m + ncells_*spin;

assert(a!=b);
HTB_(b,a) = -1.0*t_hoppings[neigh]*0.5;
HTB_(a,b) = conj(HTB_(b,a)); 

HTB_(b+2*ncells_,a+2*ncells_) = 1.0*t_hoppings[neigh]*0.5;
HTB_(a+2*ncells_,b+2*ncells_) = conj(HTB_(b+2*ncells_,a+2*ncells_));


}
}}

}}



}

#endif
