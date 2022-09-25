#ifndef Parameters_class
#define Parameters_class
#include "tensor_type.h"

class Parameters{

public:

    double V1_param, V2_param, V3_param, Psi_param, omega_param, a_moire;
    double TwistTheta, a_monolayer, MStar, eps_DE;
    int Grid_L1, Grid_L2, ns;
    int BZ_L1, BZ_L2;
    double Vz_;
    int q_slices, theta_slices, r_slices;


//	q_slices=200
//theta_slices=200
//r_sclices=181



    void Initialize(string inputfile_);
    double matchstring(string file,string match);
    string matchstring2(string file,string match);

};


void Parameters::Initialize(string inputfile_){


    cout << "____________________________________" << endl;
    cout << "Reading the inputfile name: " << inputfile_ << endl;
    cout << "____________________________________" << endl;

    V1_param = matchstring(inputfile_,"V1_param_in_meV");
    V2_param = matchstring(inputfile_,"V2_param_in_meV");
    V3_param = matchstring(inputfile_,"V3_param_in_meV");
    Psi_param = matchstring(inputfile_,"Psi_param_in_radians");
    omega_param = matchstring(inputfile_,"omega_param_in_meV");
    Grid_L1 = int(matchstring(inputfile_,"Grid_ReciprocalLattice_L1"));
    Grid_L2 = int(matchstring(inputfile_,"Grid_ReciprocalLattice_L2"));
    BZ_L1 = int(matchstring(inputfile_,"BZ_L1"));
    BZ_L2 = int(matchstring(inputfile_,"BZ_L2"));

    q_slices=int(matchstring(inputfile_,"q_slices"));
    theta_slices=int(matchstring(inputfile_,"theta_slices"));
    r_slices=int(matchstring(inputfile_,"r_slices"));


    a_monolayer = matchstring(inputfile_,"a_monolayer_in_angstorm");
    TwistTheta = matchstring(inputfile_,"Twist_Theta_in_radians");
    a_moire = a_monolayer/abs(TwistTheta);
    cout<<"a_moire (in Angstorm)= "<<a_moire<<endl;
    MStar = matchstring(inputfile_,"MStar_in_RestMass");
    eps_DE = matchstring(inputfile_, "eps_DE");
    Vz_=matchstring(inputfile_,"Layer_Potential_Diff_Vz_in_meV");

    ns=Grid_L1*Grid_L2;
}


double Parameters::matchstring(string file,string match) {
    string test;
    string line;
    ifstream readFile(file);
    double amount;
    bool pass=false;
    while (std::getline(readFile, line)) {
        std::istringstream iss(line);
        if (std::getline(iss, test, '=') && pass==false) {
            // ---------------------------------
            if (iss >> amount && test==match) {
                // cout << amount << endl;
                pass=true;
            }
            else {
                pass=false;
            }
            // ---------------------------------
            if(pass) break;
        }
    }
    if (pass==false) {
        string errorout=match;
        errorout+="= argument is missing in the input file!";
        throw std::invalid_argument(errorout);
    }
    cout << match << " = " << amount << endl;
    return amount;
}

string Parameters::matchstring2(string file,string match) {

    string line;
    ifstream readFile(file);
    string amount;
    int offset;

    if(readFile.is_open())
    {
        while(!readFile.eof())
        {
            getline(readFile,line);

            if ((offset = line.find(match, 0)) != string::npos) {
                amount = line.substr (offset+match.length()+1);				}

        }
        readFile.close();
    }
    else
    {cout<<"Unable to open input file while in the Parameters class."<<endl;}




    cout << match << " = " << amount << endl;
    return amount;
}

#endif



