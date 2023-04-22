#ifndef Parameters_TL_class
#define Parameters_TL_class
#include "tensor_type.h"

class Parameters_TL
{

public:
    int lx, ly, ns;
    
    double BoundaryConnection_X, BoundaryConnection_Y;
    bool PBC_X, PBC_Y;

    double t_nn;
    double V_attract, J_Hund;
    
    //For observse
    double beta_T, mus;
    double eta, domega;

    double Temperature;

    double Total_Particles;
    int IterMax;
    double Convergence_Error;
    int RandomSeed;
    double alpha_OP;
    bool Self_consistency;


    void Initialize(string inputfile_);
    double matchstring(string file, string match);
    string matchstring2(string file, string match);

};

void Parameters_TL::Initialize(string inputfile_)
{


    cout << "____________________________________" << endl;
    cout << "Reading the inputfile for specific model: " << inputfile_ << endl;
    cout << "____________________________________" << endl;



  
    t_nn = matchstring(inputfile_, "t_nn");
    cout<<"t_nn = "<<t_nn<<endl;

    V_attract = double(matchstring(inputfile_, "V_attraction"));
    J_Hund = double(matchstring(inputfile_, "J_Hund"));


    lx = int(matchstring(inputfile_, "Xsite"));
    ly = int(matchstring(inputfile_, "Ysite"));


    BoundaryConnection_X = double(matchstring(inputfile_, "PBC_X"));
    BoundaryConnection_Y = double(matchstring(inputfile_, "PBC_Y"));
    if(BoundaryConnection_X==1.0){
        PBC_X=true;
    }
    else{
        PBC_X=false;
    }
    if(BoundaryConnection_Y==1.0){
        PBC_Y=true;
    }
    else{
        PBC_Y=false;
    }


    ns = lx * ly;
    cout << "TotalNumberOf Unit cells = " << ns << endl;


    double Self_consistency_double;
    Self_consistency_double=double(matchstring(inputfile_,"Self_consistency"));
    if(Self_consistency_double==1){
        Self_consistency=true;}
    else{
        Self_consistency=false;
    }

    cout << "____________________________________" << endl;
}

double Parameters_TL::matchstring(string file, string match)
{
    string test;
    string line;
    ifstream readFile(file);
    double amount;
    bool pass = false;
    while (std::getline(readFile, line))
    {
        std::istringstream iss(line);
        if (std::getline(iss, test, '=') && pass == false)
        {
            // ---------------------------------
            if (iss >> amount && test == match)
            {
                // cout << amount << endl;
                pass = true;
            }
            else
            {
                pass = false;
            }
            // ---------------------------------
            if (pass)
                break;
        }
    }
    if (pass == false)
    {
        string errorout = match;
        errorout += "= argument is missing in the input file!";
        throw std::invalid_argument(errorout);
    }
    cout << match << " = " << amount << endl;
    return amount;
}

string Parameters_TL::matchstring2(string file, string match)
{

    string line;
    ifstream readFile(file);
    string amount;
    int offset;

    if (readFile.is_open())
    {
        while (!readFile.eof())
        {
            getline(readFile, line);

            if ((offset = line.find(match, 0)) != string::npos)
            {
                amount = line.substr(offset + match.length() + 1);
            }
        }
        readFile.close();
    }
    else
    {
        cout << "Unable to open input file while in the Parameters_TL class." << endl;
    }

    cout << match << " = " << amount << endl;
    return amount;
}

#endif