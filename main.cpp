#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include "functions_real.h"

#include "Skyrmion.h"

int main(){



SKYRMION my_skyrmion;

my_skyrmion.Lx=21;
my_skyrmion.Ly=21;
my_skyrmion.Diameter=20;
my_skyrmion.Beta=0.2;
my_skyrmion.Spin_Size=1.0;

my_skyrmion.Initialize_Skyrmion();
my_skyrmion.Create_Skyrmion();
my_skyrmion.Print_Skyrmion("Skyrmion_1.txt");



return 0;
}
