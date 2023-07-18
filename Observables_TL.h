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
        void Calculate_Akxw_ribbon();
        void Calculate_LDOS_SiteResolved();

	Parameters_TL &Parameters_;
	Coordinates_TL &Coordinates_;
	Hamiltonian_TL &Hamiltonian_;


};

void Observables_TL::Initialize(){
lx_=Parameters_.lx;
ly_=Parameters_.ly;
ns_=Parameters_.ns;
}



void Observables_TL::Calculate_LDOS_SiteResolved(){



double w_min=Hamiltonian_.eigs_[0] -1.0;
double w_max=Hamiltonian_.eigs_[Hamiltonian_.eigs_.size()-1] +1.0;
double dw=0.01;
double eta=0.005;
double dos_;
int comp;

for(int SiteNo=0;SiteNo<Parameters_.TotalSites;SiteNo++){
string fileout = "Site_" +to_string(SiteNo) +"_SiteResolvedDOS.txt" ;
ofstream file_out(fileout.c_str());
comp = SiteNo + 0*(Parameters_.TotalSites);

double w_val=w_min;
while(w_val<=w_max){
dos_=0.0;
for(int n=0;n<Hamiltonian_.eigs_.size();n++){
dos_ += abs(Hamiltonian_.Ham_(comp, n))*abs(Hamiltonian_.Ham_(comp, n))*Lorentzian(eta, w_val - Hamiltonian_.eigs_[n]);
}

file_out<<w_val<< "   "<<dos_<<endl;

w_val +=dw;
}
}//Site

}


void Observables_TL::Calculate_Akxw_ribbon()
{

    Mat_1_doub eigs_ = Hamiltonian_.eigs_;

    //---------Read from input file-----------------------//
    string fileout = "Akx_w.txt" ;
    double omega_min, omega_max, d_omega;
    double eta;
    int Nby2_ = eigs_.size()/2;

    if(Hamiltonian_.BdG_bool){
//    omega_min = eigs_[Nby2_] - 0.5;
 //   omega_max = eigs_[Nby2_] + 0.5;
	omega_min = eigs_[0] - 0.5;
        omega_max = eigs_[eigs_.size()-1] + 0.5;
	}
	else{
	omega_min = eigs_[0] - 0.1;
	omega_max = eigs_[eigs_.size()-1] + 0.1;
	}


    d_omega = 0.005*(omega_max-omega_min);
    eta = 1.0*d_omega;
    //---------------------------------------------------//

    int omega_index_max = int((omega_max - omega_min) / (d_omega));

    ofstream file_Akw_out(fileout.c_str());

    int c1, c2;

    Mat_3_Complex_doub A_nk;
    A_nk.resize(2);
    for(int i=0;i<2;i++){
        A_nk[i].resize(eigs_.size());
        for(int n=0;n<eigs_.size();n++){
            A_nk[i][n].resize(ns_);

        }
    }


    complex<double> temp_doub;
    int local_dof;
    double kx_val, ky_val;
    //----- A_\alpha\sigma(n,kx,iy) = sum_{i} Psi_{i,alpha, sigma ;n}exp(iota*K_vec \cdot r_vec) --------
    for(int n=0;n<eigs_.size();n++){
        for(int kx=0;kx<lx_;kx++){
            kx_val = (2.0*PI_*kx)/(1.0*lx_) ;
            for(int iy=0;iy<ly_;iy++){

                    for(int spin=0;spin<2;spin++){
                        local_dof=spin;

                        temp_doub=0.0;

                        for(int ix=0;ix<lx_;ix++){

			   if(Parameters_.BdG_bool){
                            c1 = Coordinates_.Ncell(ix, iy) + spin*(lx_*ly_) + 2*ns_;
			    c2 = Coordinates_.Ncell(ix, iy) + spin*(lx_*ly_);
                            temp_doub += (Hamiltonian_.Ham_(c1, n) + Hamiltonian_.Ham_(c2, n))*exp(-1.0*iota_comp* (kx_val*ix));
				}
				else{
                            c2 = Coordinates_.Ncell(ix, iy) + spin*(lx_*ly_);
                            temp_doub += (Hamiltonian_.Ham_(c2, n))*exp(-1.0*iota_comp* (kx_val*ix));
				}
                        }
                        A_nk[local_dof][n][Coordinates_.Ncell(kx,iy)]=temp_doub;
                    }
                
            }
        }
        if(n%100==0){
            cout << "n = "<<n<<" done"<<endl;}

    }


    cout<< "A_{alpha,sigma}(n,kx,iy) is done"<<endl;


    Mat_3_Complex_doub A_kw;

    A_kw.resize(2);
    for(int i=0;i<2;i++){
        A_kw[i].resize(Coordinates_.lx_);
        for(int n=0;n<Coordinates_.lx_;n++){
            A_kw[i][n].resize(omega_index_max);
        }
    }



    complex<double> Nup_check(0, 0);
    complex<double> Ndn_check(0, 0);



    // A_alpha_sigma(kx,w) = 1/N^2 \sum_{n,iy} |A_alpha, sigma(n,kx,iy)|^2 delta(w-eps_{n})
    int k, kp, k2;
        for(int spin=0;spin<2;spin++){
            local_dof=spin;

            for(int kx=0;kx<lx_;kx++){

                for(int w_no=0;w_no<omega_index_max;w_no++){

                    A_kw[local_dof][kx][w_no]=0.0;

                    for(int iy=0;iy<ly_;iy++){
                        k=Coordinates_.Ncell(kx,iy);
                        for(int n=0;n<eigs_.size();n++){
                            A_kw[local_dof][kx][w_no] += A_nk[local_dof][n][k]*conj(A_nk[local_dof][n][k])*
                                                         Lorentzian(eta,omega_min + (w_no * d_omega) -eigs_[n]);

                        }
                    }
                    A_kw[local_dof][kx][w_no] = A_kw[local_dof][kx][w_no]*(1.0/(lx_*lx_));

                }
                if(kx%10==0){
                    cout << "kx = "<<kx<<" done"<<endl;
                }
            }

        }
    




    cout << "Nup_check = " << Nup_check << endl;
    cout << "Ndn_check = " << Ndn_check << endl;


    double kx, ky;
    int kx_i, ky_i;

    Mat_1_intpair k_path;
    k_path.clear();
    Mat_1_intpair k_path2;
    k_path2.clear();
    pair_int temp_pair;


    // ---k_path---------

    //-------- 1d path-----------------
    kx_i = 0;
    for (kx_i = 0; kx_i < lx_; kx_i++)
    {
        temp_pair.first = kx_i;
        temp_pair.second = 0;
        k_path.push_back(temp_pair);
    }
    //----------------------------------



    //                  because in gnuplot use "set pm3d corners2color c1"
    temp_pair.first = 0;
    temp_pair.second = 0;
    k_path.push_back(temp_pair);

    //----------------------------------

    //----k_path done-------





    file_Akw_out<<"#k_point     kx    omega_val    omega_ind        Akw[spin=0]     Akw[spin=1]"<<endl;
    for (int k_point = 0; k_point < k_path.size(); k_point++)
    {

        kx_i = k_path[k_point].first;
        ky_i = k_path[k_point].second;
        kx = (2.0 * PI * kx_i) / (1.0 * Parameters_.lx);
        //ky = (2.0 * PI * ky_i) / (1.0 * Parameters_.ly);

        for (int omega_ind = 0; omega_ind < omega_index_max; omega_ind++)
        {

            //Use 1:6:7----for gnuplot
            file_Akw_out << k_point << "   " << kx<< "    " << omega_min + (d_omega * omega_ind) << "   " << omega_ind << "    ";

                for(int spin=0;spin<2;spin++){
                    local_dof = spin;
                    file_Akw_out << A_kw[local_dof][kx_i][omega_ind].real()<<"     ";
                }
            
            file_Akw_out << endl;
        }
        file_Akw_out << endl;
    }

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
