#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <random>
#include <complex>
#include <cmath>
#include <cassert>
using namespace std;
#include "Matrix.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "Hamiltonian.h"


#include "random"


int main(int argc, char *argv[]) {

    string ex_string_original =argv[0];

    string ex_string;
    //ex_string.substr(ex_string_original.length()-5);
    ex_string=ex_string_original.substr (2);
    cout<<"'"<<ex_string<<"'"<<endl;



    if(ex_string=="MoireBands"){
        string model_inputfile = argv[1];

        if (argc<2) { throw std::invalid_argument("USE:: executable inputfile"); }

        Parameters Parameters_;
        Parameters_.Initialize(model_inputfile);

        Coordinates Coordinates_(Parameters_.Grid_L1, Parameters_.Grid_L2, 1);
        Hamiltonian Hamiltonian_(Parameters_, Coordinates_);


        string file_bands_out="Bands_energy.txt";
        ofstream FileBandsOut(file_bands_out.c_str());
        FileBandsOut<<"#index kx_value ky_value E0(k)  E1(k)   E2(k) ....."<<endl;


        int n1, n2;
        Mat_1_intpair k_path;
        k_path.clear();
        pair_int temp_pair;
        int L1_,L2_;
        L1_=Parameters_.BZ_L1;
        L2_=Parameters_.BZ_L2;


	//Gamma to  K=(n1=-L1/3, n2=L2/3)
        n1=0;
  	n2=0;
        while(n2<=int(L2_/3)){
	temp_pair.first = n1;
            temp_pair.second = n2;
            k_path.push_back(temp_pair);
	n1--;
	n2++;
        //cout<<n1<<"  "<<n2<<endl;
	}


	//K to M=(n1=0, n2 = L2/2)
	n2 = int(L2_/3)+1;
        n1 = -int(L1_/3) +2;
	while(n2<=int(L2_/2)){
	temp_pair.first = n1;
        temp_pair.second = n2;
        k_path.push_back(temp_pair);	
        n2 = n2+1; 
 	n1 = n1+2;
	//cout<<n1<<"  "<<n2<<endl;
	}
	
	//assert(n1==0);
	n1=0;
	n2=int(L2_/2);	
	//M to Gamma
	n2 = n2-1;
	while(n2>=0){
        temp_pair.first = n1;
        temp_pair.second = n2;
        k_path.push_back(temp_pair);    
        n2 = n2-1;
        }



	double E_gamma_max;
	//-------------------------------------
	n1=0;
	n2=0;
	Hamiltonian_.kx_=(2.0*PI/Parameters_.a_moire)*(n1*(1.0/(sqrt(3)*L1_))  +  n2*(1.0/(sqrt(3)*L2_)));
        Hamiltonian_.ky_=(2.0*PI/Parameters_.a_moire)*(n1*(-1.0/(L1_))  +  n2*(1.0/(L2_)));
        Hamiltonian_.HTBCreate();
        char Dflag='V';
        Hamiltonian_.Diagonalize(Dflag);
	E_gamma_max = Hamiltonian_.eigs_[Hamiltonian_.Ham_.n_col()-1]; 
	//-----------------------------

        for(int index=0;index<k_path.size();index++){
        n1=k_path[index].first;
        n2=k_path[index].second;
	cout<<n1<<"  "<<n2<<endl;

        Hamiltonian_.kx_=(2.0*PI/Parameters_.a_moire)*(n1*(1.0/(sqrt(3)*L1_))  +  n2*(1.0/(sqrt(3)*L2_)));
        Hamiltonian_.ky_=(2.0*PI/Parameters_.a_moire)*(n1*(-1.0/(L1_))  +  n2*(1.0/(L2_)));
        Hamiltonian_.HTBCreate();
        char Dflag='V';
        Hamiltonian_.Diagonalize(Dflag);


	double _kx, _ky , _fk;
	_kx = Hamiltonian_.kx_*Parameters_.a_moire;
	_ky = Hamiltonian_.ky_*Parameters_.a_moire;
        //cout <<Hamiltonian_.Ham_.n_col()<<endl;
        FileBandsOut<<index<<"  "<<Hamiltonian_.kx_<<"  "<<Hamiltonian_.ky_<<"   ";
        //FileBandsOut<<1.0*( abs(  exp(iota_complex*_ky)   +   (2.0*exp(-0.5*iota_complex*_ky)*cos(sqrt(3.0)*0.5*_kx))    )   )   <<"  "<<-1.0*( abs(  exp(iota_complex*_ky)   +   (2.0*exp(-0.5*iota_complex*_ky)*cos(sqrt(3.0)*0.5*_kx))    )   )   <<"  ";
	
	_fk  = (4.0*(cos(0.5*_ky)*cos(0.5*_ky))) + 4.0*(cos(0.5*_ky)*cos(sqrt(3.0)*0.5*_kx));	


	FileBandsOut<<1.0*(sqrt(1.0+_fk))<<"  "<<-1.0*(sqrt(1.0+_fk))<<"  ";

	for(int band=Hamiltonian_.Ham_.n_col()-1;band>=0;band--){
            FileBandsOut<<Hamiltonian_.eigs_[band]-E_gamma_max<<"  ";
        }
        FileBandsOut<<endl;
        }


	Hamiltonian_.Print_Moire_Potential();
	Mat_1_int bands_;bands_.clear();

	//Honeycomb-1orb
//	bands_.push_back(0);bands_.push_back(1);

	//Honeycomb- 2orb (px,py)
	bands_.push_back(2);bands_.push_back(3);bands_.push_back(4);bands_.push_back(5);
        Hamiltonian_.Get_Wannier_function(bands_);

    }






    cout << "--------THE END--------" << endl;
} // main
