#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#define PI acos(-1.0)

#ifndef Hamiltonian_class
#define Hamiltonian_class

extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
                         std::complex<double> *,int *, double *, int *);
//zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);

class Hamiltonian {
public:

    Hamiltonian(Parameters& Parameters__, Coordinates&  Coordinates__)
        :Parameters_(Parameters__),Coordinates_(Coordinates__)

    {
        Initialize();
    }


    void Initialize();    //::DONE
    void Hoppings();        //::DONE
    double GetCLEnergy();    //::DONE
    void InteractionsCreate();   //::DONE
    void Check_Hermiticity();  //::DONE
    void HTBCreate();   //::DONE
    double chemicalpotential(double muin,double Particles);    //::DONE
    void Get_Wannier_function(Mat_1_int bands_);
    void Update_Bloch_States_Using_Projection(Mat_2_Complex_doub & Psi_state_,int space_slices,double r1_min,double d_r1,double r2_min,double d_r2);
    double Gaussian(double rx_, double ry_, double Rx_center, double Ry_center, double std_dev);
    double P_orbital_WF(double rx_, double ry_, double Rx_center, double Ry_center, double alpha, int n);
    void Print_Moire_Potential();
    double TotalDensity();   //::DONE
    double E_QM();   //::DONE
    double NIEnergy(double kx_val, double ky_val);

    void Diagonalize(char option);   //::DONE
    void copy_eigs(int i);  //::DONE

    int convert_jm_to_int(string jm_val);

    Parameters &Parameters_;
    Coordinates &Coordinates_;
    int ns_, l1_, l2_;
    double kx_, ky_;
    double k_plusx, k_minusx, k_plusy, k_minusy;

    Mat_1_doub Wnr_center_1, Wnr_center_2;
    Mat_1_doub Wnr_center_x, Wnr_center_y;

    Matrix<complex<double>> HTB_;
    Matrix<complex<double>> Ham_;
    Matrix<double> Tx,Ty,Tpxpy,Tpxmy;
    vector<double> eigs_,eigs_saved_, eigs_bands, sx_,sy_,sz_;

    //real space  effective H params
    int L1_eff, L2_eff;
    Mat_4_Complex_doub Tij;
    Mat_2_Complex_doub Uij;
    Mat_2_Complex_doub U_onsite_inter;

};


void Hamiltonian::Initialize(){

    ns_=Parameters_.ns;
    l1_=Parameters_.Grid_L1;
    l2_=Parameters_.Grid_L2;

    int space=ns_;

    HTB_.resize(space,space);
    Ham_.resize(space,space);
    eigs_.resize(space);
    eigs_saved_.resize(space);

    k_plusx = (-1.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    k_plusy = (-1.0/3.0)*(2.0*PI/Parameters_.a_moire);
    k_minusx = (-1.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    k_minusy = (1.0/3.0)*(2.0*PI/Parameters_.a_moire);


    //real space  effective H params
    L1_eff=6;L2_eff=6;
    
    Tij.resize(4);
    for(int band1=0;band1<4;band1++){
    Tij[band1].resize(4);
    for(int band2=0;band2<4;band2++){
    Tij[band1][band2].resize(L1_eff*L2_eff);
    for(int i=0;i<L1_eff*L2_eff;i++){
        Tij[band1][band2][i].resize(L1_eff*L2_eff);
    }
    }
    }


    Uij.resize(L1_eff*L2_eff);
    for(int i=0;i<L1_eff*L2_eff;i++){
        Uij[i].resize(L1_eff*L2_eff);
    }



    Wnr_center_1.resize(4);Wnr_center_2.resize(4);
    Wnr_center_x.resize(4);Wnr_center_y.resize(4);

    Wnr_center_1[0]=(-1.0*Parameters_.a_moire)/((3.0));Wnr_center_2[0]=(-1.0*Parameters_.a_moire)/((3.0));
    Wnr_center_1[1]=(-1.0*Parameters_.a_moire)/((3.0));Wnr_center_2[1]=(-1.0*Parameters_.a_moire)/((3.0));
    Wnr_center_1[2]=(-2.0*Parameters_.a_moire)/((3.0));Wnr_center_2[2]=(1.0*Parameters_.a_moire)/((3.0));
    Wnr_center_1[3]=(-2.0*Parameters_.a_moire)/((3.0));Wnr_center_2[3]=(1.0*Parameters_.a_moire)/((3.0));
    

    for(int i=0;i<4;i++){	
	 Wnr_center_x[i] =  ((Wnr_center_1[i])*(sqrt(3.0)/2.0)) + ((Wnr_center_2[i])*(sqrt(3.0)/2.0));
	 Wnr_center_y[i] = ((Wnr_center_1[i])*((-1.0)/2.0)) + ((Wnr_center_2[i])*((1.0)/2.0));
	
	cout<<"Wnr_center_rxry["<<i<<"]=("<<Wnr_center_x[i]<<", "<<Wnr_center_y[i]<<")"<<endl;

	}


    eigs_bands.resize(4);	


} // ----------

double Hamiltonian::TotalDensity(){

    double n1=0.0;
    /*
    for(int j=0;j<eigs_.size();j++){
        n1 +=  1.0f/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
    }
    */
    return n1;

} // ----------



double Hamiltonian::E_QM(){

    return 0.0;

} // ----------

double Hamiltonian::NIEnergy(double kx_val, double ky_val){

    double energy_;
    //energy_ = -1.0*(Parameters_.RedPlanckConst*Parameters_.RedPlanckConst*(kx_val*kx_val  + ky_val*ky_val))*(0.5/Parameters_.MStar);
    energy_ = -1.0*(((3.809842*1000)/Parameters_.MStar)*(kx_val*kx_val  + ky_val*ky_val));

    return energy_;
}

double Hamiltonian::GetCLEnergy(){

    return 0.0;

} // ----------


int Hamiltonian::convert_jm_to_int(string jm_val){

    int val;
    if(jm_val=="3by2_m3by2"){val=0;}
    if(jm_val=="3by2_3by2"){val=1;}
    if(jm_val=="3by2_m1by2"){val=2;}
    if(jm_val=="3by2_1by2"){val=3;}
    if(jm_val=="1by2_m1by2"){val=4;}
    if(jm_val=="1by2_1by2"){val=5;}
    return val;
}

void Hamiltonian::Check_Hermiticity()

{
    complex<double> temp(0,0);
    complex<double>temp2;

    for(int i=0;i<Ham_.n_row();i++) {
        for(int j=0;j<Ham_.n_row();j++) {
            if(
                    abs(Ham_(i,j) - conj(Ham_(j,i)))>0.00001
                    ) {
                cout<<Ham_(i,j)<<endl;
                cout<<conj(Ham_(j,i))<<endl;

            }
            assert(
                        abs(Ham_(i,j) - conj(Ham_(j,i)))<0.00001
                        ); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            //temp +=temp2*conj(temp2);
        }
    }

    // cout<<"Hermiticity: "<<temp<<endl;
}





void Hamiltonian::Diagonalize(char option){

    //extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
    //                       std::complex<double> *,int *, double *, int *);


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

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}


}

void Hamiltonian::Print_Moire_Potential(){

string file_out =  "Moire_Potential.txt";
ofstream fl_out(file_out.c_str());


Mat_1_doub VParams;
    VParams.resize(3);
    VParams[0]=Parameters_.V1_param;
    VParams[1]=Parameters_.V2_param;
    VParams[2]=Parameters_.V3_param;



   double b1x_, b1y_, b2x_, b2y_;
    b1x_=(2.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    b1y_=(0.0)*(2.0*PI/Parameters_.a_moire);
    b2x_=(1.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    b2y_=(1.0)*(2.0*PI/Parameters_.a_moire);


    double rx_min, ry_min, rx_max, ry_max, d_rx, d_ry;
    rx_min=-3.0*Parameters_.a_moire;
    ry_min=-3.0*Parameters_.a_moire;
    rx_max=3.0*Parameters_.a_moire;
    ry_max=3.0*Parameters_.a_moire;
    int r_ind;
    int space_slices=200;
    d_rx=(rx_max-rx_min)/(space_slices);
    d_ry=(ry_max-ry_min)/(space_slices);



    Mat_2_int neigh_G_shell_1, neigh_G_shell_2;
    neigh_G_shell_1.resize(3);neigh_G_shell_2.resize(3);
    for(int s=0;s<3;s++){
        neigh_G_shell_1[s].resize(6);
        neigh_G_shell_2[s].resize(6);
        }

        //shell=1
        neigh_G_shell_1[0][0]=1;neigh_G_shell_2[0][0]=0;
        neigh_G_shell_1[0][1]=0;neigh_G_shell_2[0][1]=1;
        neigh_G_shell_1[0][2]=-1;neigh_G_shell_2[0][2]=1;
        neigh_G_shell_1[0][3]=-1;neigh_G_shell_2[0][3]=0;
        neigh_G_shell_1[0][4]=0;neigh_G_shell_2[0][4]=-1;
        neigh_G_shell_1[0][5]=1;neigh_G_shell_2[0][5]=-1;

        //shell=2
        neigh_G_shell_1[1][0]=1;neigh_G_shell_2[1][0]=1;
        neigh_G_shell_1[1][1]=-1;neigh_G_shell_2[1][1]=2;
        neigh_G_shell_1[1][2]=-2;neigh_G_shell_2[1][2]=1;
        neigh_G_shell_1[1][3]=-1;neigh_G_shell_2[1][3]=-1;
        neigh_G_shell_1[1][4]=1;neigh_G_shell_2[1][4]=-2;
        neigh_G_shell_1[1][5]=2;neigh_G_shell_2[1][5]=-1;

        //shell=3
        neigh_G_shell_1[2][0]=2;neigh_G_shell_2[2][0]=0;
        neigh_G_shell_1[2][1]=0;neigh_G_shell_2[2][1]=2;
        neigh_G_shell_1[2][2]=-2;neigh_G_shell_2[2][2]=2;
        neigh_G_shell_1[2][3]=-2;neigh_G_shell_2[2][3]=0;
        neigh_G_shell_1[2][4]=0;neigh_G_shell_2[2][4]=-2;
        neigh_G_shell_1[2][5]=2;neigh_G_shell_2[2][5]=-2;




double kx_, ky_;
double dx, dy;
complex<double> Del_r;

for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
                for(int ry_ind=0;ry_ind<space_slices;ry_ind++){
  dx= rx_min + rx_ind*d_rx;
  dy= ry_min + ry_ind*d_ry;

  Del_r=0.0;
for (int s=0;s<3;s++){

                        for(int neigh_ind=0;neigh_ind<6;neigh_ind++){
                        kx_ = neigh_G_shell_1[s][neigh_ind]*b1x_ + neigh_G_shell_2[s][neigh_ind]*b2x_ ;
			ky_ = neigh_G_shell_1[s][neigh_ind]*b1y_ + neigh_G_shell_2[s][neigh_ind]*b2y_ ;

			Del_r += VParams[s]*exp(iota_complex*( (kx_*dx) + (ky_*dy) + (Parameters_.Psi_param)   ));		

                    
                        }
                     }


fl_out<<dx<<"   "<<dy<<"   "<<Del_r.real()<<"   "<<Del_r.imag()<<endl;

                     }

fl_out<<endl;
         }





}

void Hamiltonian::HTBCreate(){


    Ham_.resize(ns_,ns_);
    double b1x_, b1y_, b2x_, b2y_;
    b1x_=(2.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    b1y_=(0.0)*(2.0*PI/Parameters_.a_moire);
    b2x_=(1.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    b2y_=(1.0)*(2.0*PI/Parameters_.a_moire);

    int Bottom_, Top_;
    Bottom_=0;Top_=1;

    //l1_/2,l2_/2 is the k-point

    Mat_2_int neigh_G_shell_1, neigh_G_shell_2;
    neigh_G_shell_1.resize(3);neigh_G_shell_2.resize(3);
    for(int s=0;s<3;s++){
	neigh_G_shell_1[s].resize(6);
	neigh_G_shell_2[s].resize(6);
	}

	//shell=1
	neigh_G_shell_1[0][0]=1;neigh_G_shell_2[0][0]=0;
	neigh_G_shell_1[0][1]=0;neigh_G_shell_2[0][1]=1;
	neigh_G_shell_1[0][2]=-1;neigh_G_shell_2[0][2]=1;
	neigh_G_shell_1[0][3]=-1;neigh_G_shell_2[0][3]=0;
	neigh_G_shell_1[0][4]=0;neigh_G_shell_2[0][4]=-1;
	neigh_G_shell_1[0][5]=1;neigh_G_shell_2[0][5]=-1;
	
	//shell=2
	neigh_G_shell_1[1][0]=1;neigh_G_shell_2[1][0]=1;
        neigh_G_shell_1[1][1]=-1;neigh_G_shell_2[1][1]=2;
        neigh_G_shell_1[1][2]=-2;neigh_G_shell_2[1][2]=1;
        neigh_G_shell_1[1][3]=-1;neigh_G_shell_2[1][3]=-1;
        neigh_G_shell_1[1][4]=1;neigh_G_shell_2[1][4]=-2;
        neigh_G_shell_1[1][5]=2;neigh_G_shell_2[1][5]=-1;	

	//shell=3
        neigh_G_shell_1[2][0]=2;neigh_G_shell_2[2][0]=0;
        neigh_G_shell_1[2][1]=0;neigh_G_shell_2[2][1]=2;
        neigh_G_shell_1[2][2]=-2;neigh_G_shell_2[2][2]=2;
        neigh_G_shell_1[2][3]=-2;neigh_G_shell_2[2][3]=0;
        neigh_G_shell_1[2][4]=0;neigh_G_shell_2[2][4]=-2;
        neigh_G_shell_1[2][5]=2;neigh_G_shell_2[2][5]=-2;




    Mat_1_doub VParams;
    VParams.resize(3);
    VParams[0]=Parameters_.V1_param;
    VParams[1]=Parameters_.V2_param;
    VParams[2]=Parameters_.V3_param;
   

    double kx_local, ky_local;

    int row, col;
    int i1_neigh, i2_neigh;
    for(int i1=0;i1<l1_;i1++){
        for(int i2=0;i2<l2_;i2++){
            kx_local = kx_ + (-(l1_/2)+i1)*(b1x_) + (-(l2_/2)+i2)*(b2x_);
            ky_local = ky_ + (-(l1_/2)+i1)*(b1y_) + (-(l2_/2)+i2)*(b2y_);
            
                row=Coordinates_.Nbasis(i1, i2, 0);
                

                    //1
                    col = row;
                    Ham_(row,col) += NIEnergy(kx_local, ky_local) + (0.0*Parameters_.Vz_);


		     for (int s=0;s<3;s++){
		   
			for(int neigh_ind=0;neigh_ind<6;neigh_ind++){
			i1_neigh = i1 + neigh_G_shell_1[s][neigh_ind];
                        i2_neigh = i2 + neigh_G_shell_2[s][neigh_ind];
			
			if( (i1_neigh<l1_)  && (i2_neigh<l2_)  && (i1_neigh>=0)  && (i2_neigh>=0)  ){
			col = Coordinates_.Nbasis(i1_neigh, i2_neigh, 0);
                        Ham_(row,col) += VParams[s]*exp(iota_complex*Parameters_.Psi_param);
			}
	
			}


		     }
            
            
        }
    }



} // ----------



double Hamiltonian::Gaussian(double rx_, double ry_, double Rx_center, double Ry_center, double std_dev){
double val;
val = sqrt(1.0/(PI*std_dev*std_dev))*exp(-1.0*(  ((rx_-Rx_center)*(rx_-Rx_center) + (ry_-Ry_center)*(ry_-Ry_center))/(2.0*std_dev*std_dev)) );
return val;
}


double Hamiltonian::P_orbital_WF(double rx_, double ry_, double Rx_center, double Ry_center, double alpha, int n){
double val;
double r_val;
double sin_theta, cos_theta; // theta is in b/w x and r

r_val = sqrt((rx_ - Rx_center)*(rx_ - Rx_center) +  (ry_ - Ry_center)*(ry_ - Ry_center));
sin_theta=(ry_-Ry_center)/r_val;
cos_theta=(rx_-Rx_center)/r_val;
val = r_val*exp(((-1.0*r_val))/(2.0*alpha));

if(abs(r_val)>0.0000001){
if(n==1 || n==3){
val *=sin_theta;
}
if(n==0 || n==2){
val *=cos_theta;
}
}


//if(n==2 || n==3){
//val=-val;}

return val;
}



void Hamiltonian::Update_Bloch_States_Using_Projection(Mat_2_Complex_doub & Psi_state_,int space_slices,double r1_min,double d_r1,double r2_min,double d_r2){


double alpha=5;
int N_bands=Psi_state_.size();
Mat_2_Complex_doub Phi_state_;
Phi_state_.resize(N_bands);
for(int b=0;b<N_bands;b++){
Phi_state_[b].resize(space_slices*space_slices);
}


//Wnr_center_x.resize(N_bands); Wnr_center_y.resize(N_bands);

//For top most 2 bands [Honeycomb lattice]
assert(N_bands==4);
/*Wnr_center_x[0]=(-1.0*Parameters_.a_moire)/(sqrt(3.0));Wnr_center_y[0]=0.0;
Wnr_center_x[1]=(-1.0*Parameters_.a_moire)/(sqrt(3.0));Wnr_center_y[1]=0.0;
Wnr_center_x[2]=(-0.5*Parameters_.a_moire)/(sqrt(3.0));Wnr_center_y[2]=0.5*Parameters_.a_moire;
Wnr_center_x[3]=(-0.5*Parameters_.a_moire)/(sqrt(3.0));Wnr_center_y[3]=0.5*Parameters_.a_moire;
*/

/*
Wnr_center_x[0]=0;Wnr_center_y[0]=0.0;
Wnr_center_x[1]=0;Wnr_center_y[1]=0.0;
Wnr_center_x[2]=0;Wnr_center_y[2]=0.0;
Wnr_center_x[3]=0;Wnr_center_y[3]=0;
*/
double rx_, ry_;
double std_dev;

Mat_2_Complex_doub Overlaps_, Bloch_overlaps, Unitary_opr;

Overlaps_.resize(N_bands);Bloch_overlaps.resize(N_bands);Unitary_opr.resize(N_bands);
for(int n=0;n<N_bands;n++){
Overlaps_[n].resize(N_bands);
Bloch_overlaps[n].resize(N_bands);
Unitary_opr[n].resize(N_bands);
}


Matrix<complex<double>> S_mat, S_inv_sqrt;
S_mat.resize(N_bands, N_bands);S_inv_sqrt.resize(N_bands, N_bands);
int r_ind;


for(int n=0;n<N_bands;n++){
for(int m=0;m<N_bands;m++){
Bloch_overlaps[m][n]=0.0;
for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                r_ind=r1_ind + (space_slices)*r2_ind;
                rx_ = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
                 ry_ = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));

        Bloch_overlaps[m][n] += d_r1*d_r2*conj(Psi_state_[m][r_ind])*Psi_state_[n][r_ind];
}
}
}
}




for(int n=0;n<N_bands;n++){
for(int m=0;m<N_bands;m++){
Overlaps_[m][n]=0.0;
for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                r_ind=r1_ind + (space_slices)*r2_ind;
		rx_ = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
		 ry_ = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));

		
                std_dev = 0.01*Parameters_.a_moire;
        //Overlaps_[m][n] += d_rx*d_ry*conj(Psi_state_[m][r_ind])*Gaussian(rx_, ry_, Wnr_center_x[n], Wnr_center_y[n], std_dev);

	Overlaps_[m][n] += d_r1*d_r2*conj(Psi_state_[m][r_ind])*P_orbital_WF(rx_, ry_, Wnr_center_x[n], Wnr_center_y[n], alpha, n);

//cout<<r1_ind<<"  "<<r2_ind<<"  "<<d_r1*d_r2*P_orbital_WF(rx_, ry_, Wnr_center_x[n], Wnr_center_y[n], alpha, n)<<endl;
	//Overlaps_[m][n] += d_r1*d_r2*conj(Psi_state_[m][r_ind]);
	//Overlaps_[m][n] += d_r1*d_r2*P_orbital_WF(rx_, ry_, Wnr_center_x[n], Wnr_center_y[n], alpha, n);
		

}
}
//assert(false);
}
}

cout<<"---Overlap matrix- <Psi_k|Porb>---"<<endl;
for(int m=0;m<N_bands;m++){
for(int n=0;n<N_bands;n++){
cout<<Overlaps_[m][n]<<"  ";
}
cout<<endl;}
cout<<"---------------------"<<endl;


cout<<"---Bloch_overlap matrix <Psi_k|Psi_k>----"<<endl;
for(int m=0;m<N_bands;m++){
for(int n=0;n<N_bands;n++){
cout<<Bloch_overlaps[m][n]<<"  ";
}
cout<<endl;}
cout<<"---------------------"<<endl;



for(int n=0;n<N_bands;n++){
for(int r1_ind=0;r1_ind<space_slices;r1_ind++){ 
        for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
	r_ind=r1_ind + (space_slices)*r2_ind;
	Phi_state_[n][r_ind] = 0.0;
for(int m=0;m<N_bands;m++){
Phi_state_[n][r_ind] += Overlaps_[m][n]*Psi_state_[m][r_ind];
}
}
}
}


for(int n=0;n<N_bands;n++){
for(int m=0;m<N_bands;m++){
S_mat(n,m)=0.0;

for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
        r_ind=r1_ind + (space_slices)*r2_ind;
	S_mat(n,m) += d_r1*d_r2*conj(Phi_state_[n][r_ind])*Phi_state_[m][r_ind];//*0.0001;
}}
}}

cout<<"---Printing S_mat-------"<<endl;
S_mat.print();
//Calculating S_inv_sqrt
//diagonalizing S_mat
    vector<double> S_eigs_;
    char jobz='V';
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int nr=S_mat.n_row();
    int lda=S_mat.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3*nr -2);
    int info;
    int lwork= -1;

    S_eigs_.resize(S_mat.n_row());
    fill(S_eigs_.begin(),S_eigs_.end(),0);
    // query:
    zheev_(&jobz,&uplo,&nr,&(S_mat(0,0)),&lda,&(S_eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&nr,&(S_mat(0,0)),&lda,&(S_eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }



for(int n=0;n<N_bands;n++){
for(int m=0;m<N_bands;m++){
S_inv_sqrt(n,m)=0.0;
for(int k=0;k<N_bands;k++){
S_inv_sqrt(n,m) += S_mat(n,k)*(1.0/sqrt(S_eigs_[k]))*conj(S_mat(m,k));
}
}
}

//cout<<"S_eigs_[0] = "<<S_eigs_[0]<<endl;
cout<<"Printing S_mat after diagonalizing----------"<<endl;
S_mat.print();

cout<<"S_eigs_[0] = "<<S_eigs_[0]<<endl;
cout<<"S_eigs_[1] = "<<S_eigs_[1]<<endl;
cout<<"S_eigs_[2] = "<<S_eigs_[2]<<endl;
cout<<"S_eigs_[3] = "<<S_eigs_[3]<<endl;
cout<<"Printing S_inv_sqrt ---------------"<<endl;
S_inv_sqrt.print();



//Create Unitary_opr
for(int band_n=0;band_n<4;band_n++){
for(int band_l=0;band_l<4;band_l++){
Unitary_opr[band_n][band_l]=0;
for(int band_m=0;band_m<4;band_m++){
Unitary_opr[band_n][band_l] += S_inv_sqrt(band_m,band_n)*Overlaps_[band_l][band_m];
}
}
}


//Update Tij matrix
 int L1_,L2_;
 L1_=Parameters_.BZ_L1;
 L2_=Parameters_.BZ_L2;
double dis_x, dis_y;
int center_=(L1_eff/2) + L1_eff*(L2_eff/2);
int center_neigh;
for(int band_alpha=0;band_alpha<4;band_alpha++){
for(int band_beta=0;band_beta<4;band_beta++){
for(int r2_=0;r2_<L2_eff;r2_++){
        for(int r1_=0;r1_<L1_eff;r1_++){
	dis_x = ((sqrt(3.0)/2.0)*(r1_-(L1_eff/2)) +  (sqrt(3.0)/2.0)*(r2_-(L2_eff/2)))*Parameters_.a_moire;
        dis_y = (-0.5*(r1_-(L1_eff/2)) + 0.5*(r2_-(L2_eff/2)))*Parameters_.a_moire;

         center_neigh = (r1_) + L1_eff*(r2_);
	for(int band_alpha_p=0;band_alpha_p<4;band_alpha_p++){

//       Tij[band_alpha][band_beta][center_][center_neigh] += (1.0/(L1_*L2_))*eigs_bands[band_alpha_p]*conj(S_inv_sqrt(band_alpha_p,band_alpha))*(S_inv_sqrt(band_alpha_p,band_beta))*exp(iota_complex*(kx_*(dis_x) +  ky_*(dis_y)));

   Tij[band_alpha][band_beta][center_][center_neigh] += (1.0/(L1_*L2_))*eigs_bands[band_alpha_p]*conj(Unitary_opr[band_alpha][band_alpha_p])*(Unitary_opr[band_beta][band_alpha_p])*exp(-1.0*iota_complex*(kx_*(dis_x) +  ky_*(dis_y)));

}

        }
    }
}
}

//--------------------


for(int n=0;n<N_bands;n++){
for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
        for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
        r_ind=r1_ind + (space_slices)*r2_ind;
        Psi_state_[n][r_ind] = 0.0;
for(int m=0;m<N_bands;m++){
Psi_state_[n][r_ind] += S_inv_sqrt(m,n)*Phi_state_[m][r_ind];
}
}
}
}




}

void Hamiltonian::Get_Wannier_function(Mat_1_int bands_){

    double b1x_, b1y_, b2x_, b2y_;
    b1x_=(2.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    b1y_=(0.0)*(2.0*PI/Parameters_.a_moire);
    b2x_=(1.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    b2y_=(1.0)*(2.0*PI/Parameters_.a_moire);

    double kx_local, ky_local;
    complex<double> temp_factor;
    complex<double> checknorm;
    int eigen_no, comp;
    double dis_x,dis_y;

    Mat_1_int band_number;
    band_number=bands_;

//---remove later-----
double d_rx, d_ry, rx_min, ry_min;
//--------------------


    double r1_min, r2_min, r1_max, r2_max, d_r1, d_r2;
    r1_min=-3.0*Parameters_.a_moire + Wnr_center_1[0];
    r2_min=-3.0*Parameters_.a_moire + Wnr_center_2[0];
    r1_max=3.0*Parameters_.a_moire + Wnr_center_1[0];
    r2_max=3.0*Parameters_.a_moire + Wnr_center_2[0];
    int r_ind;
    int space_slices=Parameters_.r_slices;//N=18n+1
    d_r1=(r1_max-r1_min)/(space_slices-1);
    d_r2=(r2_max-r2_min)/(space_slices-1);

   cout<<"d_r1, d_r2 = "<<d_r1<<", "<<d_r2<<endl;
   cout<<"r1_min,r1_max = "<<r1_min<<",  "<<r1_max<<endl;
   cout<<"r2_min,r2_max = "<<r2_min<<",  "<<r2_max<<endl;
   cout<<"Wnr_center_1[0],Wnr_center_2[0] = "<<Wnr_center_1[0]<<",  "<<Wnr_center_2[0]<<endl;  
    double qx_min, qy_min, qx_max, qy_max, d_qx, d_qy;
    qx_min=-4.0*PI;//Parameters_.a_moire;
    qy_min=-4.0*PI;//Parameters_.a_moire;
    qx_max=4.0*PI;//Parameters_.a_moire;
    qy_max=4.0*PI;//Parameters_.a_moire;
    
    int q_slices=Parameters_.q_slices;
    d_qx=(qx_max-qx_min)/(q_slices);
    d_qy=(qy_max-qy_min)/(q_slices);
    double eta_q=0.001;
    double d_sqr=360000;
    double screening=0.0;


    double q_max, d_q,  d_theta;
    q_max=0.4*PI;
    d_q=(q_max)/(q_slices);
    int theta_slices=Parameters_.theta_slices;
    d_theta = (2.0*PI)/(theta_slices);




    int N_bands;
    N_bands=bands_.size();


    int L1_,L2_;
    L1_=Parameters_.BZ_L1;
    L2_=Parameters_.BZ_L2;

    int Bottom_, Top_;
    Bottom_=0;Top_=1;
    Mat_2_Complex_doub Wnr_state_, Psi_state_;
	Wnr_state_.resize(N_bands);
	Psi_state_.resize(N_bands);
	for(int band_i=0;band_i<N_bands;band_i++){
        Wnr_state_[band_i].resize(space_slices*space_slices);
        Psi_state_[band_i].resize(space_slices*space_slices);
	}


    Mat_2_Complex_doub Mq_;
    Mq_.resize(q_slices);
    for(int q1=0;q1<q_slices;q1++){
        Mq_[q1].resize(q_slices);
    }

    Mat_2_Complex_doub Vq_;
    Vq_.resize(q_slices);
    for(int q1=0;q1<q_slices;q1++){
        Vq_[q1].resize(q_slices);
    }


    Mat_4_Complex_doub Mq_SphC;
    Mq_SphC.resize(N_bands);
    for(int bands_i=0;bands_i<N_bands;bands_i++){
    Mq_SphC[bands_i].resize(N_bands);
    for(int bands_j=0;bands_j<N_bands;bands_j++){
    Mq_SphC[bands_i][bands_j].resize(q_slices);
    for(int q_i=0;q_i<q_slices;q_i++){
        Mq_SphC[bands_i][bands_j][q_i].resize(theta_slices);
    }
    }
	}



    Mat_5_Complex_doub MqR_SphC;
    MqR_SphC.resize(N_bands);
    for(int band1=0;band1<N_bands;band1++){
    MqR_SphC[band1].resize(N_bands);
   for(int band2=0;band2<N_bands;band2++){
    MqR_SphC[band1][band2].resize(q_slices);
    for(int q_i=0;q_i<q_slices;q_i++){
        MqR_SphC[band1][band2][q_i].resize(theta_slices);
        for(int th_i=0;th_i<theta_slices;th_i++){
        MqR_SphC[band1][band2][q_i][th_i].resize(L1_eff*L2_eff);
        }
    }
    }
    }




    int center_=(L1_eff/2) + L1_eff*(L2_eff/2);
    int center_neigh;
    //Tij[center_][center_p1]=0.0;
   // eigen_no=(l1_*l2_)-1-band;

   




for(int band1=0;band1<4;band1++){
for(int band2=0;band2<4;band2++){
 for(int r2_=0;r2_<L2_eff;r2_++){
        for(int r1_=0;r1_<L1_eff;r1_++){
            center_neigh = (r1_) + L1_eff*(r2_);
           Tij[band1][band2][center_][center_neigh]=0.0;
        }
    }
}
}

 

    for(int n1=0;n1<L1_;n1++){
        for(int n2=0;n2<L2_;n2++){

            cout<<"doing "<<n1<<"  "<<n2<<endl;

            kx_=(2.0*PI/Parameters_.a_moire)*(n1*(1.0/(sqrt(3)*L1_))  +  n2*(1.0/(sqrt(3)*L2_)));
            ky_=(2.0*PI/Parameters_.a_moire)*(n1*(-1.0/(L1_))  +  n2*(1.0/(L2_)));
            HTBCreate();
            char Dflag='V';
            Diagonalize(Dflag);

            //------------------


            //Hopping-real space-------
/*      
      for(int r1_=0;r1_<L1_eff;r1_++){
                for(int r2_=0;r2_<L2_eff;r2_++){
                    center_neigh = (r1_) + L1_eff*(r2_);
                    dis_x = ((sqrt(3.0)/2.0)*(r1_-(L1_eff/2)) +  (sqrt(3.0)/2.0)*(r2_-(L2_eff/2)))*Parameters_.a_moire;
                    dis_y = (-0.5*(r1_-(L1_eff/2)) + 0.5*(r2_-(L2_eff/2)))*Parameters_.a_moire;
                    Tij[center_][center_neigh]+=(1.0/(L1_*L2_))*exp(iota_complex*( kx_*(dis_x) +  ky_*(dis_y) ))*eigs_[eigen_no];
                }
            }
*/
            //--------------------

	complex<double> phase_1=conj(Ham_(0,eigen_no))*(1.0/abs(Ham_(0,eigen_no)));
	complex<double> phase_2=conj(Ham_(0, eigen_no-1))*(1.0/abs(Ham_(0,eigen_no-1)));
 	

            for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                    r_ind=r1_ind + (space_slices)*r2_ind;
                        
			for(int band_i=0;band_i<N_bands;band_i++){
			Psi_state_[band_i][r_ind]=0.0;
			eigen_no=(l1_*l2_)-1-band_number[band_i];
			//cout<<"eigen_no"<<eigen_no<<endl;     
                   for(int i1=0;i1<l1_;i1++){
                            for(int i2=0;i2<l2_;i2++){
                                kx_local = kx_ + (-(l1_/2)+i1)*(b1x_) + (-(l2_/2)+i2)*(b2x_);
                                ky_local = ky_ + (-(l1_/2)+i1)*(b1y_) + (-(l2_/2)+i2)*(b2y_);
                                comp = Coordinates_.Nbasis(i1, i2, 0);

				dis_x = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
                                dis_y = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));

                                Psi_state_[band_i][r_ind] += ((Ham_(comp,eigen_no)))*exp(iota_complex*( kx_local*(dis_x) +  ky_local*(dis_y) ));
                            }
                        }
			
			}
                       // if( (abs(rx_min + rx_ind*d_rx)<=0.00000001 &&
                         //              abs(ry_min + ry_ind*d_ry)<=0.00000001) ){
                       //     temp_factor = Psi_state_[r_ind];
                        //}
                    

                }
            }




		for(int band_i=0;band_i<N_bands;band_i++){
	         	checknorm=0.0;
			  for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                            for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                                r_ind=r1_ind + (space_slices)*r2_ind;
                                
                                checknorm += d_r1*d_r2*Psi_state_[band_i][r_ind]*conj(Psi_state_[band_i][r_ind]);
                            }
                        }
		for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                            for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                                r_ind=r1_ind + (space_slices)*r2_ind; 
                                Psi_state_[band_i][r_ind] = (Psi_state_[band_i][r_ind])/(sqrt(abs(checknorm)));
                            }
                        }

                    }

		
	 for(int band_i=0;band_i<N_bands;band_i++){
                eigen_no=(l1_*l2_)-1-band_number[band_i];
		eigs_bands[band_i]=eigs_[eigen_no];
		}

		Update_Bloch_States_Using_Projection(Psi_state_, space_slices, r1_min, d_r1, r2_min, d_r2);


            //choosing phase
        /*    for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
                for(int ry_ind=0;ry_ind<space_slices;ry_ind++){
                    r_ind=rx_ind + (space_slices)*ry_ind;
                    Psi_state_[r_ind] *= conj(temp_factor)/abs(temp_factor);
                    
                }
            }*/

            //            for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
            //                for(int ry_ind=0;ry_ind<space_slices;ry_ind++){
            //                    r_ind=rx_ind + (space_slices)*ry_ind;
            //                    for(int orb=0;orb<2;orb++){
            //                     checknorm +=   Psi_state_[orb][r_ind]*conj(Psi_state_[orb][r_ind]);
            //                    }
            //                }
            //            }

            //-------------------


            for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                    r_ind=r1_ind + (space_slices)*r2_ind;
        for(int band_i=0;band_i<N_bands;band_i++){

			dis_x = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
                        dis_y = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));

	  //          Wnr_state_[band_i][r_ind] = P_orbital_WF(dis_x, dis_y, Wnr_center_x[band_i], Wnr_center_y[band_i], 5, band_i);
	Wnr_state_[band_i][r_ind] += (1.0/sqrt(L1_*L2_))*(Psi_state_[band_i][r_ind]);
            

		    }
		}
            }
	


        }
    }



	for(int band_i=0;band_i<N_bands;band_i++){
    checknorm=0.0;
    for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
        for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
            r_ind=r1_ind + (space_slices)*r2_ind;
            checknorm += d_r1*d_r2*Wnr_state_[band_i][r_ind]*conj(Wnr_state_[band_i][r_ind]);
        }
    }

    for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
        for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
            r_ind=r1_ind + (space_slices)*r2_ind;
            Wnr_state_[band_i][r_ind] = Wnr_state_[band_i][r_ind]*(1.0/sqrt(abs(checknorm)));
        }
    }




    checknorm=0.0;
    string file_Wnr_out="Wannier_functions_band" + to_string(band_i) +".txt";
    ofstream FileWNROut(file_Wnr_out.c_str());
    FileWNROut<<"#index r1_ind   r2_ind    r1_dis    r2_dis  rx ry W_r  ....."<<endl;
    for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
        for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
            r_ind=r1_ind + (space_slices)*r2_ind;

	    dis_x = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
            dis_y = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));
	   
            FileWNROut<<r_ind<<"  "<<r1_ind<<"  "<<r2_ind<<"   "<< (r1_min + r1_ind*d_r1)/(Parameters_.a_moire)<<"   "<<(r2_min + r2_ind*d_r2)/(Parameters_.a_moire)<<"   "<<(dis_x)/(Parameters_.a_moire)<<"   "<<(dis_y)/(Parameters_.a_moire)<<"   "<<(Wnr_state_[band_i][r_ind]).real()<<"   "<< (Wnr_state_[band_i][r_ind]).imag() <<endl;
            checknorm += d_r1*d_r2*Wnr_state_[band_i][r_ind]*conj(Wnr_state_[band_i][r_ind]);
        }
        FileWNROut<<endl;
    }

    FileWNROut<<"#a_moire = "<<Parameters_.a_moire<<endl;
    FileWNROut<<"#norm = "<<checknorm.real()<<"  "<<checknorm.imag()<<endl;


	}



  cout<<"Checking Wannier functions overlaps-----------------------"<<endl;
        for(int band1=0;band1<N_bands;band1++){
	for(int band2=0;band2<N_bands;band2++){
	complex<double> overlap_;
	overlap_=0.0;
	for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
        for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
      
	    r_ind=r1_ind + (space_slices)*r2_ind;
		overlap_ += d_r1*d_r2*Wnr_state_[band1][r_ind]*conj(Wnr_state_[band2][r_ind]);

	}
	}	
	cout<<overlap_<<"  ";
        }
	cout<<endl;
        }
   cout<<"----------------------------------------------------------"<<endl;






    //cout<<"Tij[0][0+a1] = "<<Tij[center_][center_p1]<<endl;

	cout << setprecision(6)<<endl;
	for(int band1=0;band1<4;band1++){
	for(int band2=0;band2<4;band2++){
         cout<<"--------------------------------------Tij["<<band1<< "]["<<band2<< "][center][neigh]--------------------------------------"<<endl;


     for(int r2_=0;r2_<L2_eff;r2_++){
        for(int r1_=0;r1_<L1_eff;r1_++){
            center_neigh = (r1_) + L1_eff*(r2_);
            cout<<Tij[band1][band2][center_][center_neigh]<<"  ";
        }
        cout<<endl;
    }

    cout<<"----------------------------------------------------------------------------------"<<endl;

        }
	}


    //For Uij [same orbital]--------------------------------------------------------------------//

for(int band_i=0;band_i<N_bands;band_i++){
for(int band_j=0;band_j<N_bands;band_j++){

/*if(  ((band_i==0) && (band_j==1)) 
  || ((band_i==1) && (band_j==0))
  || ((band_i==2) && (band_j==3))
  || ((band_i==3) && (band_j==2))
  || ((band_i==0) && (band_j==0))
  || ((band_i==1) && (band_j==1))
  || ((band_i==2) && (band_j==2))
  || ((band_i==3) && (band_j==3))
 ){*/

double Mq_norm, Mq_Sph_norm;
Mq_norm=0.0; Mq_Sph_norm=0.0;

double q_val, theta_val;
double qx_,qy_;
for(int q_ind=0;q_ind<q_slices;q_ind++){
 q_val = q_ind*d_q;
for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
 theta_val = theta_ind*d_theta;

qx_=q_val*cos(theta_val);
qy_=q_val*sin(theta_val);

//( kx_*((1.0/(2.0*sqrt(3.0)))*Parameters_.a_moire) + ky_*( (0.5)*Parameters_.a_moire )
      Mq_SphC[band_i][band_j][q_ind][theta_ind]=0.0;
      for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
      for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
      r_ind=r1_ind + (space_slices)*r2_ind;

            dis_x = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
            dis_y = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));
      	
	
      Mq_SphC[band_i][band_j][q_ind][theta_ind] += d_r1*d_r2*(
      conj(Wnr_state_[band_i][r_ind])*(Wnr_state_[band_j][r_ind]))
                            *exp(iota_complex*(qx_*((dis_x - 0*Wnr_center_x[band_i] )) + qy_*((dis_y  - 0*Wnr_center_y[band_i] ))));
                }
            }

        Mq_Sph_norm += sqrt((qx_*qx_) + (qy_*qy_))*d_q*d_theta*abs(Mq_SphC[band_i][band_j][q_ind][theta_ind]);

}
}
cout << "Mq_Sph_norm[" <<band_i  <<", "<<band_j <<"]= "<<Mq_Sph_norm<<endl;

//}
}
}
cout.precision(5);






//MqR_SphC[band1][band2][q_i][th_i][r1_+r2_*L1_eff];



for(int band_i=0;band_i<N_bands;band_i++){
for(int band_j=0;band_j<N_bands;band_j++){
double q_val, theta_val;
double qx_,qy_;

/* cout<<"Calculating Uij for Wannier orbital = "<<band_i<<endl;
   double V_;
    for(int r1_=0;r1_<L1_eff;r1_++){
        for(int r2_=0;r2_<L2_eff;r2_++){
            center_neigh = (r1_) + L1_eff*(r2_);
            dis_x = ((sqrt(3.0)/2.0)*(r1_-((1.0*L1_eff)/2.0)) +  (sqrt(3.0)/2.0)*(r2_-((1.0*L2_eff)/2.0)))*Parameters_.a_moire;
            dis_y = (-0.5*(r1_-((1.0*L1_eff)/2.0)) + 0.5*(r2_-((1.0*L2_eff)/2.0)))*Parameters_.a_moire;

            //dis_x = ((r1_-(L1_eff/2)))*Parameters_.a_moire;
            //dis_y = ((r2_-(L2_eff/2)))*Parameters_.a_moire;
            Uij[center_][center_neigh]=0.0;
                for(int q_ind=0;q_ind<q_slices;q_ind++){
                 q_val = q_ind*d_q;
                for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
                 theta_val = theta_ind*d_theta;
                qx_=q_val*cos(theta_val);
                qy_=q_val*sin(theta_val);


      //Mq_SphC[q_ind][theta_ind]
                    V_= (2*PI*14.399*1000)/(Parameters_.eps_DE);
                    Uij[center_][center_neigh]+= (1.0/(4.0*PI*PI))*(d_q*d_theta)*(V_*abs(Mq_SphC[band_i][q_ind][theta_ind])*abs(Mq_SphC[band_i][q_ind][theta_ind]))
                            *exp(iota_complex*( qx_*(dis_x) +  qy_*(dis_y) ));
                }
            }
        }
    }


    cout<<endl;
    cout<<"--------------------------------------Uij[center][neigh]--------------------------------------"<<endl;
     for(int r2_=L2_eff-1;r2_>=0;r2_--){
        for(int r1_=0;r1_<L1_eff;r1_++){
            center_neigh = (r1_) + L1_eff*(r2_);
            cout<<Uij[center_][center_neigh]<<"  ";
        }
        cout<<endl;
    }

    cout<<"----------------------------------------------------------------------------------"<<endl;

*/


    string MqSphC_file="MqSphC_orb" + to_string(band_i) + "_" + to_string(band_j) +".txt";
    ofstream MqSphCFILE(MqSphC_file.c_str());

    for(int q_ind=0;q_ind<q_slices;q_ind++){
        q_val = q_ind*d_q;
        for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
        theta_val = theta_ind*d_theta;
            MqSphCFILE<<q_val<<"   "<<theta_val<<"   "<<Mq_SphC[band_i][band_j][q_ind][theta_ind].real()<<"   "<<Mq_SphC[band_i][band_j][q_ind][theta_ind].imag()<<endl;
        }
        MqSphCFILE<<endl;
    }

}
}



cout<<"-----------ONSITE INTER/INTRA COULOMB REPULSION-------------"<<endl;
double q_val, theta_val;
double qx_,qy_;
U_onsite_inter.resize(N_bands);
for(int band1=0;band1<N_bands;band1++){
U_onsite_inter[band1].resize(N_bands);
for(int band2=0;band2<N_bands;band2++){
U_onsite_inter[band1][band2]=0.0;
}
}


cout<<"------band1  band2   U_onsite_inter[band1][band2]----------------------"<<endl;
for(int band1=0;band1<N_bands;band1++){

//int band2=(1-(band1%2))+2*(band1/2);
for(int band2=0;band2<N_bands;band2++){
double V_;

            //dis_x = ((r1_-(L1_eff/2)))*Parameters_.a_moire;
            //dis_y = ((r2_-(L2_eff/2)))*Parameters_.a_moire;
                for(int q_ind=0;q_ind<q_slices;q_ind++){
                 q_val = q_ind*d_q;
                for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
                 theta_val = theta_ind*d_theta;
                qx_=q_val*cos(theta_val);
                qy_=q_val*sin(theta_val);


      //Mq_SphC[q_ind][theta_ind]
                    V_= (2*PI*14.399*1000)/(Parameters_.eps_DE);
                    U_onsite_inter[band1][band2]+= (1.0/(4.0*PI*PI))*(d_q*d_theta)*(V_*conj(Mq_SphC[band1][band1][q_ind][theta_ind])*(Mq_SphC[band2][band2][q_ind][theta_ind]));
		    //U_onsite_inter[band1][band1]+= (1.0/(4.0*PI*PI))*(d_q*d_theta)*(V_*conj(Mq_SphC[band1][band1][q_ind][theta_ind])*(Mq_SphC[band1][band1][q_ind][theta_ind]));
                }
            }

cout<<U_onsite_inter[band1][band2]<<"  ";
//cout<<band1<<"  "<<band2<<"   "<<U_onsite_inter[band1][band2]<<"   "<<endl;


}
cout<<endl;
}
cout<<"----------------------------------------"<<endl;
	




cout<<"------+a1-a2   U_onsite_inter[band1][band2]----------------------"<<endl;
int r1_neigh, r2_neigh;
r1_neigh=(L1_eff/2)+1;
r2_neigh=(L2_eff/2)-1;
dis_x = ((sqrt(3.0)/2.0)*(r1_neigh-((1.0*L1_eff)/2.0)) +  (sqrt(3.0)/2.0)*(r2_neigh-((1.0*L2_eff)/2.0)))*Parameters_.a_moire;
dis_y = (-0.5*(r1_neigh-((1.0*L1_eff)/2.0)) + 0.5*(r2_neigh-((1.0*L2_eff)/2.0)))*Parameters_.a_moire;

for(int band1=0;band1<N_bands;band1++){

//int band2=(1-(band1%2))+2*(band1/2);
for(int band2=0;band2<N_bands;band2++){
double V_;
complex<double> val_;

		val_=0.0;
            //dis_x = ((r1_-(L1_eff/2)))*Parameters_.a_moire;
            //dis_y = ((r2_-(L2_eff/2)))*Parameters_.a_moire;
                for(int q_ind=0;q_ind<q_slices;q_ind++){
                 q_val = q_ind*d_q;
                for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
                 theta_val = theta_ind*d_theta;
                qx_=q_val*cos(theta_val);
                qy_=q_val*sin(theta_val);


      //Mq_SphC[q_ind][theta_ind]
                    V_= (2*PI*14.399*1000)/(Parameters_.eps_DE);
                    val_+= (1.0/(4.0*PI*PI))*(d_q*d_theta)*(V_*conj(Mq_SphC[band1][band1][q_ind][theta_ind])*(Mq_SphC[band2][band2][q_ind][theta_ind]))*exp(iota_complex* (qx_*(dis_x) + qy_*(dis_y))  );
                    //U_onsite_inter[band1][band1]+= (1.0/(4.0*PI*PI))*(d_q*d_theta)*(V_*conj(Mq_SphC[band1][band1][q_ind][theta_ind])*(Mq_SphC[band1][band1][q_ind][theta_ind]));
                }
            }

cout<<val_<<"  ";
//cout<<band1<<"  "<<band2<<"   "<<U_onsite_inter[band1][band2]<<"   "<<endl;


}
cout<<endl;
}
cout<<"----------------------------------------"<<endl;


//band2 is present at "-a2" w.r.t to site with band1
cout<<"------  '-a2'   U_onsite_inter[band1][band2]----------------------"<<endl;
r1_neigh=(L1_eff/2);
r2_neigh=(L2_eff/2)-1;
dis_x = ((sqrt(3.0)/2.0)*(r1_neigh-((1.0*L1_eff)/2.0)) +  (sqrt(3.0)/2.0)*(r2_neigh-((1.0*L2_eff)/2.0)))*Parameters_.a_moire;
dis_y = (-0.5*(r1_neigh-((1.0*L1_eff)/2.0)) + 0.5*(r2_neigh-((1.0*L2_eff)/2.0)))*Parameters_.a_moire;

for(int band1=0;band1<N_bands;band1++){

//int band2=(1-(band1%2))+2*(band1/2);
for(int band2=0;band2<N_bands;band2++){
double V_;
complex<double> val_;

                val_=0.0;
            //dis_x = ((r1_-(L1_eff/2)))*Parameters_.a_moire;
            //dis_y = ((r2_-(L2_eff/2)))*Parameters_.a_moire;
                for(int q_ind=0;q_ind<q_slices;q_ind++){
                 q_val = q_ind*d_q;
                for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
                 theta_val = theta_ind*d_theta;
                qx_=q_val*cos(theta_val);
                qy_=q_val*sin(theta_val);


      //Mq_SphC[q_ind][theta_ind]
                    V_= (2*PI*14.399*1000)/(Parameters_.eps_DE);
                    val_+= (1.0/(4.0*PI*PI))*(d_q*d_theta)*(V_*conj(Mq_SphC[band1][band1][q_ind][theta_ind])*(Mq_SphC[band2][band2][q_ind][theta_ind]))*exp(iota_complex* (qx_*(dis_x) + qy_*(dis_y))  );
                    //U_onsite_inter[band1][band1]+= (1.0/(4.0*PI*PI))*(d_q*d_theta)*(V_*conj(Mq_SphC[band1][band1][q_ind][theta_ind])*(Mq_SphC[band1][band1][q_ind][theta_ind]));
                }
            }

cout<<val_<<"  ";
//cout<<band1<<"  "<<band2<<"   "<<U_onsite_inter[band1][band2]<<"   "<<endl;


}
cout<<endl;
}
cout<<"----------------------------------------"<<endl;





cout<<"-------HUNDS COUPLING -------------------"<<endl;
for(int band1=0;band1<N_bands;band1++){
for(int band2=0;band2<N_bands;band2++){

//if(  ((band1==0) && (band2==1))
 // || ((band1==1) && (band2==0))
 // || ((band1==2) && (band2==3))
 // || ((band1==3) && (band2==2))
//)
//{

 double V_;
 complex<double> val_=0.0;
                for(int q_ind=0;q_ind<q_slices;q_ind++){
                 q_val = q_ind*d_q;
                for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
                 theta_val = theta_ind*d_theta;
                qx_=q_val*cos(theta_val);
                qy_=q_val*sin(theta_val);


      //Mq_SphC[q_ind][theta_ind]
                    V_= (2*PI*14.399*1000)/(Parameters_.eps_DE);
                    val_+= (1.0/(4.0*PI*PI))*(d_q*d_theta)*(V_*conj(Mq_SphC[band1][band2][q_ind][theta_ind])*(Mq_SphC[band1][band2][q_ind][theta_ind]));
                }
            }

cout<<val_<<"  ";

//}
}
cout<<endl;
}

cout<<"----------------------------------------"<<endl;






cout<<"-------Interaction Assisted Hopping Type-1 [band1][band2], hops from band2 to band1(with already 1 el at band1)-------------------"<<endl;
for(int band1=0;band1<N_bands;band1++){
for(int band2=0;band2<N_bands;band2++){

//if(  ((band1==0) && (band2==1))
 // || ((band1==1) && (band2==0))
 // || ((band1==2) && (band2==3))
 // || ((band1==3) && (band2==2))
//)
//{

 double V_;
 complex<double> val_=0.0;
                for(int q_ind=0;q_ind<q_slices;q_ind++){
                 q_val = q_ind*d_q;
                for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
                 theta_val = theta_ind*d_theta;
                qx_=q_val*cos(theta_val);
                qy_=q_val*sin(theta_val);


      //Mq_SphC[q_ind][theta_ind]
                    V_= (2*PI*14.399*1000)/(Parameters_.eps_DE);
                    val_+= (1.0/(4.0*PI*PI))*(d_q*d_theta)*(V_*conj(Mq_SphC[band1][band1][q_ind][theta_ind])*(Mq_SphC[band1][band2][q_ind][theta_ind]));
                }
            }

cout<<val_<<"  ";

//}
}
cout<<endl;
}

cout<<"----------------------------------------"<<endl;





for(int band_fixed=0;band_fixed<N_bands;band_fixed++){
cout<<"-------Interaction Assisted Hopping Type-2 for fixed band="<<band_fixed << ", Mat[band1][band2], hops from band2 to band1(with already 1 el at band1)-------------------"<<endl;
for(int band1=0;band1<N_bands;band1++){
for(int band2=0;band2<N_bands;band2++){

//if(  ((band1==0) && (band2==1))
 // || ((band1==1) && (band2==0))
 // || ((band1==2) && (band2==3))
 // || ((band1==3) && (band2==2))
//)
//{

 double V_;
 complex<double> val_=0.0;
                for(int q_ind=0;q_ind<q_slices;q_ind++){
                 q_val = q_ind*d_q;
                for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
                 theta_val = theta_ind*d_theta;
                qx_=q_val*cos(theta_val);
                qy_=q_val*sin(theta_val);


      //Mq_SphC[q_ind][theta_ind]
                    V_= (2*PI*14.399*1000)/(Parameters_.eps_DE);
                    val_+= (1.0/(4.0*PI*PI))*(d_q*d_theta)*(V_*conj(Mq_SphC[band_fixed][band_fixed][q_ind][theta_ind])*(Mq_SphC[band1][band2][q_ind][theta_ind]));
                }
            }

cout<<val_<<"  ";

//}
}
cout<<endl;
}

cout<<"----------------------------------------"<<endl;

}

Mat_1_int Orbital1_array,Orbital2_array,Orbital3_array,Orbital4_array;
Orbital1_array.resize(N_bands*N_bands*N_bands*N_bands);
Orbital2_array.resize(N_bands*N_bands*N_bands*N_bands);
Orbital3_array.resize(N_bands*N_bands*N_bands*N_bands);
Orbital4_array.resize(N_bands*N_bands*N_bands*N_bands);
Mat_1_Complex_doub Vals_array;Vals_array.resize(N_bands*N_bands*N_bands*N_bands);

cout<<"-------Pair Hopping Mat[band1][band2][band3][band4] pair from:bands 3,4 hop to bands 1,2   -------------------"<<endl;

int counter=0;
for(int band3=0;band3<N_bands;band3++){
for(int band4=0;band4<N_bands;band4++){
for(int band1=0;band1<N_bands;band1++){
for(int band2=0;band2<N_bands;band2++){

//if(  ((band1==0) && (band2==1))
 // || ((band1==1) && (band2==0))
 // || ((band1==2) && (band2==3))
 // || ((band1==3) && (band2==2))
//)
//{

 double V_;
 complex<double> val_=0.0;
                for(int q_ind=0;q_ind<q_slices;q_ind++){
                 q_val = q_ind*d_q;
                for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
                 theta_val = theta_ind*d_theta;
                qx_=q_val*cos(theta_val);
                qy_=q_val*sin(theta_val);


      //Mq_SphC[q_ind][theta_ind]
                 V_= (2*PI*14.399*1000)/(Parameters_.eps_DE);
                 val_+= (1.0/(4.0*PI*PI))*(d_q*d_theta)*(V_*conj(Mq_SphC[band3][band1][q_ind][theta_ind])*(Mq_SphC[band4][band2][q_ind][theta_ind]));
                }
            }

//cout<<"["<< band1<< "]"<<"["<< band2<< "]"<<"["<<band3 << "]"<<"["<< band4<< "]"<<"    "<<val_<<"  "<<endl;

Orbital1_array[counter]=band1;
Orbital2_array[counter]=band2;
Orbital3_array[counter]=band3;
Orbital4_array[counter]=band4;
Vals_array[counter]=val_;

counter++;
}
//cout<<endl;
}
}}


//Arranging in decreasing order;
assert(counter==(N_bands*N_bands*N_bands*N_bands));
complex<double> max_val=0;
Mat_1_Complex_doub Vals_array_ordered;
Vals_array_ordered=Vals_array;
Mat_1_int Orbital1_array_ordered,Orbital2_array_ordered,Orbital3_array_ordered,Orbital4_array_ordered;
Orbital1_array_ordered.resize(N_bands*N_bands*N_bands*N_bands);
Orbital2_array_ordered.resize(N_bands*N_bands*N_bands*N_bands);
Orbital3_array_ordered.resize(N_bands*N_bands*N_bands*N_bands);
Orbital4_array_ordered.resize(N_bands*N_bands*N_bands*N_bands);
int i_max;

string file_onsite_interaction_str= "Onsite_interactions.txt";
ofstream fileout_onsite_interaction(file_onsite_interaction_str.c_str());

fileout_onsite_interaction<<"#band1 band2 band3 band4  V(b1,b2,b3,b4)=<b1,b2|b3,b4>"<<endl;
for(int j=0;j<N_bands*N_bands*N_bands*N_bands;j++){

//Finding maximum val
for(int i=0;i<N_bands*N_bands*N_bands*N_bands;i++){
if(abs(Vals_array[i])>=abs(max_val)){
max_val=Vals_array[i];
i_max=i;
}
}


Vals_array_ordered[j]=Vals_array[i_max];
Vals_array[i_max]=0.0;
max_val=0;

Orbital1_array_ordered[j]=Orbital1_array[i_max];
Orbital2_array_ordered[j]=Orbital2_array[i_max];
Orbital3_array_ordered[j]=Orbital3_array[i_max];
Orbital4_array_ordered[j]=Orbital4_array[i_max];


fileout_onsite_interaction<<"["<< Orbital1_array_ordered[j]<< "]"<<"["<< Orbital2_array_ordered[j]<< "]"<<"["<<Orbital3_array_ordered[j] << "]"<<"["<< Orbital4_array_ordered[j]<< "]"<<"    "<<Vals_array_ordered[j]<<"  "<<endl;


}
cout<<" Onsite interactions printed in decreasing order"<<endl;





double dis_1, dis_2;
int r1_ind_shift, r2_ind_shift;
complex<double> val_temp1_;
for(int r1_=0;r1_<L1_eff;r1_++){
        for(int r2_=0;r2_<L2_eff;r2_++){

        if( ((r1_==int((L1_eff/2)+1)) && (r2_==int((L2_eff/2)-1))) ||
          ((r1_==int((L1_eff/2))) && (r2_==int((L2_eff/2)-1))) 
	// || ((r1_==int((L1_eff/2))) && (r2_==int((L2_eff/2))))
	
  ){

        dis_1 =  (r1_-((1.0*L1_eff)/2.0))*Parameters_.a_moire;
        dis_2 =  (r2_-((1.0*L2_eff)/2.0))*Parameters_.a_moire;


for(int band_i=0;band_i<N_bands;band_i++){
for(int band_j=0;band_j<N_bands;band_j++){

cout<<"Mqr r1,r2,band1, band2: "<<r1_<<",  "<<r2_<<",  "<<band_i<<",  "<<band_j<<endl;

/*if(  ((band_i==0) && (band_j==1)) 
  || ((band_i==1) && (band_j==0))
  || ((band_i==2) && (band_j==3))
  || ((band_i==3) && (band_j==2))
  || ((band_i==0) && (band_j==0))
  || ((band_i==1) && (band_j==1))
  || ((band_i==2) && (band_j==2))
  || ((band_i==3) && (band_j==3))
 ){*/

double Mq_norm, Mq_Sph_norm;
Mq_norm=0.0; Mq_Sph_norm=0.0;

double q_val, theta_val;
double qx_,qy_;
for(int q_ind=0;q_ind<q_slices;q_ind++){
 q_val = q_ind*d_q;
for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
 theta_val = theta_ind*d_theta;

qx_=q_val*cos(theta_val);
qy_=q_val*sin(theta_val);

//( kx_*((1.0/(2.0*sqrt(3.0)))*Parameters_.a_moire) + ky_*( (0.5)*Parameters_.a_moire )
      MqR_SphC[band_i][band_j][q_ind][theta_ind][r1_+L1_eff*r2_]=0.0;
      for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
      for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
      //r_ind=r1_ind + (space_slices)*r2_ind;

        int r1_ind_r = int((r1_ind*d_r1 + (dis_1) + d_r1*0.5)/(d_r1));
        int r2_ind_r = int((r2_ind*d_r2 + (dis_2) + d_r2*0.5)/(d_r2));


            dis_x = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
            dis_y = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));

   if( r1_ind_r >=0 && r1_ind_r<space_slices &&
            r2_ind_r >=0 && r2_ind_r<space_slices
    ) {

     int r_ind_r=r1_ind_r + (space_slices)*r2_ind_r;
     int r_ind_l=r1_ind + (space_slices)*r2_ind;


      MqR_SphC[band_i][band_j][q_ind][theta_ind][r1_+r2_*L1_eff] += d_r1*d_r2*(
      conj(Wnr_state_[band_i][r_ind_l])*(Wnr_state_[band_j][r_ind_r]))
                            *exp(iota_complex*(qx_*((dis_x - 0*Wnr_center_x[band_i] )) + qy_*((dis_y  - 0*Wnr_center_y[band_i] ))));
                }


            }
        }
}
}
}
}

}
}
}



//Here band-1 belongs to "+a1-a2" site w.r.t. band-2 site.


for(int r1_neigh_offset=0;r1_neigh_offset<2;r1_neigh_offset++){

cout<<"-------InterSite Exchange ("<<r1_neigh_offset<<"*a1, -a2) -------------------"<<endl;
r1_neigh=(L1_eff/2)+r1_neigh_offset;
r2_neigh=(L2_eff/2)-1;
dis_x = ((sqrt(3.0)/2.0)*(r1_neigh-((1.0*L1_eff)/2.0)) +  (sqrt(3.0)/2.0)*(r2_neigh-((1.0*L2_eff)/2.0)))*Parameters_.a_moire;
dis_y = (-0.5*(r1_neigh-((1.0*L1_eff)/2.0)) + 0.5*(r2_neigh-((1.0*L2_eff)/2.0)))*Parameters_.a_moire;


double dis_x_temp, dis_y_temp;


counter=0;

for(int dof1_=0;dof1_<2;dof1_++){ // dof=unitcell=sublattice (as unit cell comes with only specific sublattice here, because of n.n interaction)
for(int dof2_=0;dof2_<2;dof2_++){
for(int dof3_=0;dof3_<2;dof3_++){
for(int dof4_=0;dof4_<2;dof4_++){

if(dof2_==dof3_){
dis_x_temp=0;
dis_y_temp=0;
}
if(dof2_==1 && dof3_==0){
dis_x_temp=dis_x;
dis_y_temp=dis_y;
}
if(dof2_==0 && dof3_==1){
dis_x_temp=-1.0*dis_x;
dis_y_temp=-1.0*dis_y;
}



for(int orb1_=0;orb1_<2;orb1_++){
for(int orb2_=0;orb2_<2;orb2_++){
for(int orb3_=0;orb3_<2;orb3_++){
for(int orb4_=0;orb4_<2;orb4_++){

int band1_= orb1_ + dof1_*2;
int band2_= orb2_ + dof2_*2;
int band3_= orb3_ + dof3_*2;
int band4_= orb4_ + dof4_*2; 

//dof=0 means center site and sublattice 0
//dof=1 means neighbour (+a1,-a2) site and sublattice 1 


double V_;
 complex<double> val_=0.0;
complex<double> mat_elmt_left, mat_elmt_right;
                for(int q_ind=0;q_ind<q_slices;q_ind++){
                 q_val = q_ind*d_q;
                for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
                 theta_val = theta_ind*d_theta;
                qx_=q_val*cos(theta_val);
                qy_=q_val*sin(theta_val);


      //Mq_SphC[q_ind][theta_ind]
                    V_= (2*PI*14.399*1000)/(Parameters_.eps_DE);

	if(dof3_==dof1_ ){
	mat_elmt_left = conj(Mq_SphC[band3_][band1_][q_ind][theta_ind]);
	}
	if(dof3_==1 && dof1_==0){
	mat_elmt_left = conj(MqR_SphC[band3_][band1_][q_ind][theta_ind][r1_neigh+r2_neigh*L1_eff]);
	}
	if(dof3_==0 && dof1_==1){
	mat_elmt_left = conj(MqR_SphC[band1_][band3_][q_ind][theta_ind][r1_neigh+r2_neigh*L1_eff]*exp(iota_complex* (qx_*(dis_x) + qy_*(dis_y))) );
	}


	if(dof4_==dof2_ ){
        mat_elmt_right = (Mq_SphC[band2_][band4_][q_ind][theta_ind]);
        }
        if(dof2_==1 && dof4_==0){
        mat_elmt_right = (MqR_SphC[band2_][band4_][q_ind][theta_ind][r1_neigh+r2_neigh*L1_eff]);
        }
        if(dof2_==0 && dof4_==1){
        mat_elmt_right = (MqR_SphC[band4_][band2_][q_ind][theta_ind][r1_neigh+r2_neigh*L1_eff])*exp(iota_complex* (qx_*(dis_x) + qy_*(dis_y)));;
        }
	



        val_+= (1.0/(4.0*PI*PI))*(d_q*d_theta)*(V_*mat_elmt_left*mat_elmt_right)*exp(iota_complex* (qx_*(dis_x_temp) + qy_*(dis_y_temp)));
                }
            }

//cout<<val_<<"  ";

Orbital1_array[counter]=band1_;
Orbital2_array[counter]=band2_;
Orbital3_array[counter]=band3_;
Orbital4_array[counter]=band4_;
Vals_array[counter]=val_;
counter++;



}}}}


}}}}



//Arranging in decreasing order;
assert(counter==(N_bands*N_bands*N_bands*N_bands));
max_val=0;
int site1_, sublattice1_, orb1_;
int site2_, sublattice2_, orb2_;
int site3_, sublattice3_, orb3_;
int site4_, sublattice4_, orb4_;
Vals_array_ordered=Vals_array;
Orbital1_array_ordered.resize(N_bands*N_bands*N_bands*N_bands);
Orbital2_array_ordered.resize(N_bands*N_bands*N_bands*N_bands);
Orbital3_array_ordered.resize(N_bands*N_bands*N_bands*N_bands);
Orbital4_array_ordered.resize(N_bands*N_bands*N_bands*N_bands);
string file_NNneigh_interaction_str;
if(r1_neigh_offset==0){
file_NNneigh_interaction_str= "NearestNeighbour_ma2_interactions.txt";}
if(r1_neigh_offset==1){
file_NNneigh_interaction_str= "NearestNeighbour_pa1_ma2_interactions.txt";
}
ofstream fileout_NNneigh_interaction(file_NNneigh_interaction_str.c_str());


fileout_NNneigh_interaction<<"#site=1 is the neigh of site=0"<<endl;
fileout_NNneigh_interaction<<"#site1 sublattice1 orb1 site2 sublattice2 orb2 site3 sublattice3 orb3 site4 sublattice4 orb4  V(site1 subltc1 orb1 , -2- , -3- , -4-)"<<endl;
for(int j=0;j<N_bands*N_bands*N_bands*N_bands;j++){

//Finding maximum val
for(int i=0;i<N_bands*N_bands*N_bands*N_bands;i++){
if(abs(Vals_array[i])>=abs(max_val)){
max_val=Vals_array[i];
i_max=i;
}
}


Vals_array_ordered[j]=Vals_array[i_max];
Vals_array[i_max]=0.0;
max_val=0;

Orbital1_array_ordered[j]=Orbital1_array[i_max];
Orbital2_array_ordered[j]=Orbital2_array[i_max];
Orbital3_array_ordered[j]=Orbital3_array[i_max];
Orbital4_array_ordered[j]=Orbital4_array[i_max];

orb1_= Orbital1_array_ordered[j]%2;
orb2_= Orbital2_array_ordered[j]%2;
orb3_= Orbital3_array_ordered[j]%2;
orb4_= Orbital4_array_ordered[j]%2;

site1_ = int((1.0*Orbital1_array_ordered[j] +0.5)/(2.0));sublattice1_ = site1_;
site2_ = int((1.0*Orbital2_array_ordered[j] +0.5)/(2.0));sublattice2_ = site2_;
site3_ = int((1.0*Orbital3_array_ordered[j] +0.5)/(2.0));sublattice3_ = site3_;
site4_ = int((1.0*Orbital4_array_ordered[j] +0.5)/(2.0));sublattice4_ = site4_;


//fileout_onsite_interaction<<"["<< Orbital1_array_ordered[j]<< "]"<<"["<< Orbital2_array_ordered[j]<< "]"<<"["<<Orbital3_array_ordered[j] << "]"<<"["<< Orbital4_array_ordered[j]<< "]"<<"    "<<Vals_array_ordered[j]<<"  "<<endl;

fileout_NNneigh_interaction<<site1_<<" "<<sublattice1_<<" "<<orb1_<<" "<<site2_<<" "<<sublattice2_<<" "<<orb2_<<" "<<site3_<<" "<<sublattice3_<<" "<<orb3_<<" "<<site4_<<" "<<sublattice4_<<" "<<orb4_<<" "<<Vals_array_ordered[j]<<endl;

}
cout<<" Nearest neighbour ("<<r1_neigh_offset<<"*a1, -a2) interactions printed in decreasing order"<<endl;

cout<<"----------------------------------------"<<endl;

}//r1_neigh_offset





}

void Hamiltonian::Hoppings(){

} // ----------

void Hamiltonian::copy_eigs(int i){

    int space=2*ns_;

    if (i == 0) {
        for(int j=0;j<space;j++) {
            eigs_[j] = eigs_saved_[j];
        }
    }
    else {
        for(int j=0;j<space;j++) {
            eigs_saved_[j] = eigs_[j];
        }
    }

}


#endif
