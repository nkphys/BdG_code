#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include "functions_real.h"

#include "Skyrmion.h"
#include "Hamiltonian_TL.h"
#include "Coordinates_TL.h"
#include "Parameters_TL.h"
#include "Observables_TL.h"

int main(int argc, char *argv[]){


string ex_string_original = argv[0];
string model_inputfile = argv[1];



Parameters_TL Parameters_TL_;
Parameters_TL_.Initialize(model_inputfile);

string my_skyrmion_type=Parameters_TL_.Skyrmion_Type;
SKYRMION my_skyrmion(my_skyrmion_type);

my_skyrmion.Lx=Parameters_TL_.lx;
my_skyrmion.Ly=Parameters_TL_.ly;
my_skyrmion.Diameter=Parameters_TL_.Skyrmion_Diameter;

my_skyrmion.Spin_Size=1.0;

my_skyrmion.Initialize_Skyrmion();
my_skyrmion.Beta=Parameters_TL_.Skyrmion_Beta;
my_skyrmion.BraviasLattice=Parameters_TL_.Lattice_Type;

my_skyrmion.Create_Skyrmion();
my_skyrmion.Skyrmion_Number_Calculate();
my_skyrmion.Magnetization_Calculate();
my_skyrmion.Print_Skyrmion("MySkyrmion.txt");

 
//assert(false);

//-------------Check-------


Coordinates_TL Coordinates_TL_(Parameters_TL_.lx,Parameters_TL_.ly, 1);


Hamiltonian_TL Hamiltonian_TL_(Parameters_TL_, Coordinates_TL_, my_skyrmion);
Hamiltonian_TL_.HamilCreation();
char Dflag='V';
Hamiltonian_TL_.Diagonalize(Dflag);



Observables_TL Observables_TL_(Parameters_TL_, Coordinates_TL_, Hamiltonian_TL_);

Observables_TL_.Print_Eigvals("Eigenvals.txt");
Observables_TL_.Print_DOS("dos_check.txt");

Observables_TL_.Calculate_Akxw_ribbon();

//------------------------


return 0;
}
