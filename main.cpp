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
int main(int argc, char *argv[]){


string ex_string_original = argv[0];
string model_inputfile = argv[1];

string my_skyrmion_type="AntiSkyrmion";
SKYRMION my_skyrmion(my_skyrmion_type);

my_skyrmion.Lx=51;
my_skyrmion.Ly=51;
my_skyrmion.Diameter=50;

my_skyrmion.Spin_Size=1.0;

my_skyrmion.vorticity=-1.0;
my_skyrmion.helicity=0;
my_skyrmion.polarity=-1;

my_skyrmion.Initialize_Skyrmion();
my_skyrmion.Beta=0.08;
my_skyrmion.Create_Skyrmion();
my_skyrmion.Skyrmion_Number_Calculate();
my_skyrmion.Print_Skyrmion("MySkyrmion.txt");



//-------------Check-------
Parameters_TL Parameters_TL_;
Parameters_TL_.Initialize(model_inputfile);


Coordinates_TL Coordinates_TL_(Parameters_TL_.lx,Parameters_TL_.ly, 1);


Hamiltonian_TL Hamiltonian_TL_(Parameters_TL_, Coordinates_TL_);
Hamiltonian_TL_.HamilCreation();
char Dflag='V';
Hamiltonian_TL_.Diagonalize(Dflag);
Hamiltonian_TL_.Print_DOS("dos_check.txt");

//------------------------


return 0;
}
