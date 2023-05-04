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
    double V_attract, J_Hund, mu;
    int JJ_width;
    
    //For observse
    double beta_T, mus;
    double eta, domega;

    double BdG_double;
    bool BdG_bool;
    double Temperature;
    double Boltzman_Const;

    double Total_Particles;
    int IterMax;
    double Convergence_Error;
    int RandomSeed;
    double alpha_OP;
    bool Self_consistency;


    double Delta_s;

    bool SpinFermionTerm, PairingTerm, JJ_Channel, ChemicalPotentialTerm;

    void Initialize(string inputfile_);
    double matchstring(string file, string match);
    string matchstring2(string file, string match);

};

void Parameters_TL::Initialize(string inputfile_)
{


    Boltzman_Const=1.0;
    mus=0.0;




    cout << "____________________________________" << endl;
    cout << "Reading the inputfile for specific model: " << inputfile_ << endl;
    cout << "____________________________________" << endl;



  
    t_nn = matchstring(inputfile_, "t_nn");
    cout<<"t_nn = "<<t_nn<<endl;

    V_attract = double(matchstring(inputfile_, "V_attraction"));
    J_Hund = double(matchstring(inputfile_, "J_Hund"));
    Delta_s = double(matchstring(inputfile_, "Delta_s_local"));
    mu = double(matchstring(inputfile_, "mu"));
    JJ_width = int(matchstring(inputfile_, "JJ_width"));


    Temperature = double(matchstring(inputfile_, "Temperature"));
    beta_T = Boltzman_Const/Temperature;
    

    lx = int(matchstring(inputfile_, "Xsite"));
    ly = int(matchstring(inputfile_, "Ysite"));



    BdG_double = double(matchstring(inputfile_,"BdG_bool"));
    if(BdG_double==1.0){
	BdG_bool=true;
	}
	else{
	BdG_bool=false;
	}

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


    double SpinFermionTerm_double;
    SpinFermionTerm_double=double(matchstring(inputfile_,"SpinFermionTerm"));
    if(SpinFermionTerm_double==1){
        SpinFermionTerm=true;}
    else{
        SpinFermionTerm=false;
    }

    double PairingTerm_double;
    PairingTerm_double=double(matchstring(inputfile_,"PairingTerm"));
    if(PairingTerm_double==1){
        PairingTerm=true;}
    else{
        PairingTerm=false;
    }
    
    double JJ_Channel_double;
    JJ_Channel_double=double(matchstring(inputfile_,"JJ_Channel"));
    if(JJ_Channel_double==1){
        JJ_Channel=true;}
    else{
        JJ_Channel=false;
    }

    double ChemicalPotentialTerm_double;
    ChemicalPotentialTerm_double=double(matchstring(inputfile_,"ChemicalPotentialTerm"));
    if(ChemicalPotentialTerm_double==1){
        ChemicalPotentialTerm=true;}
    else{
        ChemicalPotentialTerm=false;
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
