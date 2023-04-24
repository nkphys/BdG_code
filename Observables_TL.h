#include "Parameters_TL.h"
#include "Coordinates_TL.h"
#include "Hamiltonian_TL.h"
#include "tensor_type.h"
#include "functions_real.h"

#ifndef OBSERVABLES_H
#define OBSERVABLES_H


//HERE
//dsytri (UPLO, N, A, LDA, IPIV, WORK, INFO )
extern "C" void  dsytri_(char *, int *, double *, int *, int *, double *, int *);


//dsytrf(UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO)
extern "C" void  dsytrf_(char *, int *, double *, int *, int *, double *, int *, int *);


class Observables_TL{
public:

    Observables_TL(Parameters_TL& Parameters__, Coordinates_TL& Coordinates__, Hamiltonian_TL& Hamiltonian__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), Hamiltonian_(Hamiltonian__)
    {
        Initialize();
    }




	void Initialize();
	double lx_, ly_, ns_;
        void Print_DOS(string DOSfile);
 	void Print_Eigvals(string EigValsfile);
        void Calculate_Conductance();

	Parameters_TL &Parameters_;
	Coordinates_TL &Coordinates_;
	Hamiltonian_TL &Hamiltonian_;


};

void Observables_TL::Initialize(){
lx_=Parameters_.lx;
ly_=Parameters_.ly;
ns_=Parameters_.ns;
}



void Observables_TL::Print_Eigvals(string EigValsfile){


ofstream EIGVALFILE_(EigValsfile.c_str());
EIGVALFILE_<<"#n   E(n)"<<endl;

for(int n=0;n<Hamiltonian_.eigs_.size();n++){
EIGVALFILE_<<n<<"  "<<Hamiltonian_.eigs_[n]<<endl;
}


}



void Observables_TL::Print_DOS(string DOSfile){


ofstream DOSFILE_(DOSfile.c_str());
DOSFILE_<<"#w   DOS(w)"<<endl;

double w_min=Hamiltonian_.eigs_[0] -1.0;
double w_max=Hamiltonian_.eigs_[Hamiltonian_.eigs_.size()-1] +1.0;
double dw=0.01;
double eta=0.005;
double dos_;

double w_val=w_min;
while(w_val<=w_max){
dos_=0.0;
for(int n=0;n<Hamiltonian_.eigs_.size();n++){
dos_ += Lorentzian(eta, w_val - Hamiltonian_.eigs_[n]);
}

DOSFILE_<<w_val<< "   "<<dos_<<endl;

w_val +=dw;
}


}




void Observables_TL::Calculate_Conductance(){

/*    string fileout_sigma = "Conductance.txt";
    ofstream file_sigma_out(fileout_sigma.c_str());
    file_sigma_out<<"#omega  "<<endl;
*/

    //--------------------------------------------------//
    double eta = 0.005;
    //---------------------------------------------------//



    Mat_1_string Directions_str;
    Directions_str.push_back("px");Directions_str.push_back("py");Directions_str.push_back("pxmy");


    Mat_1_int Directions_int; //values taken from Coordinates_TL class.
    Directions_int.push_back(1);Directions_int.push_back(2);Directions_int.push_back(7);


    int i_1_, i_2_;
    Mat_4_Complex_doub PSI_;
    Mat_1_Complex_doub Vals_;
   
    PSI_.resize(3);
    Vals_.resize(3);
    for(int i=0;i<3;i++){
    PSI_[i].resize(3);
    Vals_[i]=0.0;;
	for(int j=0;j<3;j++){
	PSI_[i][j].resize(2*ns_);
        for(int l=0;l<2*ns_;l++){
	PSI_[i][j][l].resize(2*ns_);
       for(int k=0;k<2*ns_;k++){
	PSI_[i][j][l][k]=0.0;
	}
    }
    }
    }

  


   //ofstream checkfilePSI("PSIcheck.txt");
   // checkfilePSI<<"#n  m  val"<<endl; 


	for(int dir1_=0;dir1_<3;dir1_++){
                for(int dir2_=0;dir2_<3;dir2_++){

    for(int n=0;n<2*ns_;n++){
        for(int m=0;m<2*ns_;m++){

            
	Vals_[dir1_]=0.0;
	Vals_[dir2_]=0.0;
                

            for(int i=0;i<ns_;i++){
		
                i_1_ = Coordinates_.neigh(i,Directions_int[dir1_]);
                i_2_ = Coordinates_.neigh(i,Directions_int[dir2_]);


                for(int spin=0;spin<2;spin++){
                    Vals_[dir1_] += ( conj(Hamiltonian_.Ham_(i_1_ + (ns_*spin),n))*Hamiltonian_.Ham_(i + (ns_*spin),m) )
                            -  ( conj(Hamiltonian_.Ham_(i + (ns_*spin),n))*Hamiltonian_.Ham_(i_1_ + (ns_*spin),m) );

                    Vals_[dir2_] += ( conj(Hamiltonian_.Ham_(i_2_ + (ns_*spin),n))*Hamiltonian_.Ham_(i + (ns_*spin),m) )
                            -  ( conj(Hamiltonian_.Ham_(i + (ns_*spin),n))*Hamiltonian_.Ham_(i_2_ + (ns_*spin),m) );

                }
            }


            PSI_[dir1_][dir2_][n][m] = Vals_[dir1_]*conj(Vals_[dir2_]);
           
//	if((dir1_==dir2_) && (dir1_==0)){
//	checkfilePSI<<n<<"  "<<m<<"  "<<PSI_[dir1_][dir2_][n][m].real()<<"  "<<PSI_[dir1_][dir2_][n][m].imag()<<endl;
//	}
 
        }

//	if((dir1_==dir2_) && (dir1_==0)){
//        checkfilePSI<<endl;
 //       }
    }

	cout<<"PSI_["<<Directions_str[dir1_]<<"]["<<Directions_str[dir2_]<<"] DONE"<<endl;

	}
	}



    Mat_2_Complex_doub sigma_;
    sigma_.resize(3);
    for(int i=0;i<3;i++){
    sigma_[i].resize(3);
    for(int j=0;j<3;j++){
	sigma_[i][j]=0.0;
    }
    }


	cout<<Parameters_.mus<<"   "<<Parameters_.beta_T<<endl;
	
	for(int dir1_=0;dir1_<3;dir1_++){
                for(int dir2_=0;dir2_<3;dir2_++){

        for(int n=0;n<2*ns_;n++){
            for(int m=0;m<2*ns_;m++){

                if(n!=m){
                    sigma_[dir1_][dir2_] += (PSI_[dir1_][dir2_][n][m])*
				(
                             ((1.0/( exp(1.0*(Hamiltonian_.eigs_[m]-Parameters_.mus)*Parameters_.beta_T ) + 1.0)))
                            -((1.0/( exp(1.0*(Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta_T ) + 1.0)))
				)
                            *Lorentzian(Hamiltonian_.eigs_[n] - Hamiltonian_.eigs_[m], eta);
                }
            }
        }

        sigma_[dir1_][dir2_] = sigma_[dir1_][dir2_]*PI_*(1.0/(ns_));




	
	cout<<"#sigma_["<<Directions_str[dir1_]<<"]["<<Directions_str[dir2_]<<"] = "<<sigma_[dir1_][dir2_]<<endl;

		}}


}







#endif
