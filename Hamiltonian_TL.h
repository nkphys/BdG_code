#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor_type.h"
#include "Parameters_TL.h"
#include "Coordinates_TL.h"
#include "Matrix.h"
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
    Hamiltonian_TL(Parameters_TL  &Parameters__ , Coordinates_TL &Coordinates__)
        :Coordinates_(Coordinates__),Parameters_(Parameters__)

    {
        Initialize();
        HTBCreate();
    }


    void Initialize();
    void HTBCreate();
    void HamilCreation();
    void Diagonalize(char option);
    void Print_DOS(string DOSfile);

    Coordinates_TL &Coordinates_;
    Parameters_TL &Parameters_;
    int lx_, ly_, ncells_, n_orbs_;

    Mat_1_doub eigs_;
    Matrix<complex<double>> HTB_;
    Matrix<complex<double>> Ham_; 



};


void Hamiltonian_TL::Print_DOS(string DOSfile){


ofstream DOSFILE_(DOSfile.c_str());
DOSFILE_<<"#w   DOS(w)"<<endl;

double w_min=eigs_[0] -1.0;
double w_max=eigs_[eigs_.size()-1] +1.0;
double dw=0.01;
double eta=0.05;
double dos_;

double w_val=w_min;
while(w_val<=w_max){
dos_=0.0;
for(int n=0;n<eigs_.size();n++){
dos_ += Lorentzian(eta, w_val - eigs_[n]); 
}

DOSFILE_<<w_val<< "   "<<dos_<<endl;

w_val +=dw;
}


}

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

    // Ham_.print();

}




void Hamiltonian_TL::HamilCreation(){

Ham_=HTB_;

}





void Hamiltonian_TL::Initialize()
{

    ly_ = Parameters_.ly;
    lx_ = Parameters_.lx;
    ncells_ = lx_ * ly_;
    n_orbs_ = 1;
    int space = 2 * ncells_;

    HTB_.resize(space, space);
    Ham_.resize(space, space);
   
} // ----------

void Hamiltonian_TL::HTBCreate(){

Mat_1_int t_neighs;
Mat_1_doub t_hoppings;


t_neighs.push_back(0);t_neighs.push_back(7);t_neighs.push_back(3);
t_hoppings.push_back(Parameters_.t_nn); t_hoppings.push_back(Parameters_.t_nn); t_hoppings.push_back(Parameters_.t_nn);

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

HTB_(b,a) = -1.0*t_hoppings[neigh];
HTB_(a,b) = conj(HTB_(b,a)); 
}
}}

}}



}

#endif
