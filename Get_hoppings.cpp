#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <complex>
#include <vector>
#include <math.h>
using namespace std;
#define PI_ 3.14159265

typedef vector< complex<double> >  Mat_1_Complex_doub;
typedef vector<Mat_1_Complex_doub> Mat_2_Complex_doub;
typedef vector<Mat_2_Complex_doub> Mat_3_Complex_doub;
typedef vector<Mat_3_Complex_doub> Mat_4_Complex_doub;

typedef vector<double>  Mat_1_doub;
typedef vector<Mat_1_doub> Mat_2_doub;


int matchstring_and_linenumber(string file_name, string hopp_type){


int line_no, line_final;
size_t found_pos;
bool found;
found=false;
string line;
ifstream readFile(file_name.c_str());

if(readFile.is_open()){
line_no=0;
while(!readFile.eof()){
getline(readFile, line);

found_pos = line.find(hopp_type);
if( found_pos != string::npos){
found=true;
line_final=line_no;
readFile.close();
}


line_no++;
}
}

if(found==false){
cout<<" \""<<hopp_type<<"\" is not present in file \""<<file_name<<"\" "<<endl;
}

return line_final;

}


int main(int argc, char *argv[]){

string file_name=argv[1];
cout<<"file used : \""<<file_name<<"\""<<endl;

ifstream in_file(file_name.c_str());

string hopp_type="Tij[0][0][center][neigh]";
int line_no=matchstring_and_linenumber(file_name, hopp_type);
//cout<<line_no<<endl;

int N1_sites=6;

 Mat_4_Complex_doub Thop;
Thop.resize(4);
for(int dof1=0;dof1<4;dof1++){
Thop[dof1].resize(4);
for(int dof2=0;dof2<4;dof2++){
Thop[dof1][dof2].resize(N1_sites);
for(int i1=0;i1<N1_sites;i1++){
Thop[dof1][dof2][i1].resize(N1_sites);
}
}
}

string line_str;
for(int line_counter=0;line_counter<line_no;line_counter++){
getline(in_file,line_str);
}
//cout<<line_str<<endl;


for(int dof1=0;dof1<4;dof1++){
for(int dof2=0;dof2<4;dof2++){
getline(in_file,line_str);
cout<<line_str<<endl;
for(int i2=0;i2<N1_sites;i2++){
for(int i1=0;i1<N1_sites;i1++){
in_file>>Thop[dof1][dof2][i1][i2];
//cout<<Thop[dof1][dof2][i1][i2]<<"  ";
}
//cout<<endl;
}
getline(in_file,line_str);
//cout<<line_str<<endl;
getline(in_file,line_str);
//cout<<line_str<<endl;
}
}




double Dis_max=1.16;
 //in units of a_m 0.58 for only HC nearest neighbour, 1.001 for upto Next-nearest neighbout on Honeycomb lattice
// 1.16 for Next to Next-nearest neighbour (NNNN)
double dis_1x, dis_1y, dis_2x, dis_2y, dis_;
complex<double> val_;
int i1_center, i2_center;
i1_center=N1_sites/2;
i2_center=i1_center;

//inter-unitcell hopping (it includes only Honeycomb NN)
string t0_hop_file="t0_mat_hop.txt";
ofstream t0_hop(t0_hop_file.c_str());

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<2;orb1++){
for(int atom1=0;atom1<2;atom1++){

for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<2;orb2++){
for(int atom2=0;atom2<2;atom2++){


dis_1x=atom1*(1.0/(sqrt(3.0)*2.0));
dis_1y=atom1*(1.0/(2.0));

dis_2x=0.0 + atom2*(1.0/(sqrt(3.0)*2.0));
dis_2y=0.0 + atom2*(1.0/(2.0));

dis_=sqrt( (dis_1x-dis_2x)*(dis_1x-dis_2x)  + (dis_1y-dis_2y)*(dis_1y-dis_2y)  );

if(dis_<Dis_max){
val_= abs(abs(1.0*spin1-1.0*spin2)-1.0)*Thop[2*atom1+orb1][2*atom2+orb2][i1_center][i2_center];
}
else{
val_=0.0;
}

if(abs(val_)<0.00001){
val_=0.0;
}

t0_hop<<val_<<"  ";

}}}
t0_hop<<endl;
}}}





//Nearest plusa1 neigbour-unitcell hopping (it includes only Honeycomb NN, NNN, NNNN)
string t1_plus_a1_hop_file="t1_plus_a1_hop.txt";
ofstream t1_plus_a1_hop(t1_plus_a1_hop_file.c_str());

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<2;orb1++){
for(int atom1=0;atom1<2;atom1++){

for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<2;orb2++){
for(int atom2=0;atom2<2;atom2++){

dis_1x=atom1*(1.0/(sqrt(3.0)*2.0));
dis_1y=atom1*(1.0/(2.0));

dis_2x=(-sqrt(3.0)/2.0) + atom2*(1.0/(sqrt(3.0)*2.0));
dis_2y= (1.0/2.0)   + atom2*(1.0/(2.0));

dis_=sqrt( (dis_1x-dis_2x)*(dis_1x-dis_2x)  + (dis_1y-dis_2y)*(dis_1y-dis_2y)  );

if(dis_<Dis_max){
val_= abs(abs(1.0*spin1-1.0*spin2)-1.0)*Thop[2*atom1+orb1][2*atom2+orb2][i1_center-1][i2_center];
}
else{
val_=0.0;
}

if(abs(val_)<0.00001){
val_=0.0;
}


t1_plus_a1_hop<<val_<<"  ";

}}}
t1_plus_a1_hop<<endl;
}}}


//Nearest plusa1 neigbour-unitcell hopping (it includes only Honeycomb NN, NNN, NNNN)
//NNNN is not included by hand
string t1_minus_a2_hop_file="t1_minus_a2_hop.txt";
ofstream t1_minus_a2_hop(t1_minus_a2_hop_file.c_str());

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<2;orb1++){
for(int atom1=0;atom1<2;atom1++){

for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<2;orb2++){
for(int atom2=0;atom2<2;atom2++){

dis_1x=atom1*(1.0/(sqrt(3.0)*2.0));
dis_1y=atom1*(1.0/(2.0));

dis_2x=((1.0*sqrt(3.0))/2.0) + atom2*(1.0/(sqrt(3.0)*2.0));
dis_2y= (1.0/2.0)   + atom2*(1.0/(2.0));

dis_=sqrt( (dis_1x-dis_2x)*(dis_1x-dis_2x)  + (dis_1y-dis_2y)*(dis_1y-dis_2y)  );

if(dis_<Dis_max){
val_= abs(abs(1.0*spin1-1.0*spin2)-1.0)*Thop[2*atom1+orb1][2*atom2+orb2][i1_center][i2_center+1];
}
else{
val_=0.0;
}

if(abs(val_)<0.00001){
val_=0.0;
}


t1_minus_a2_hop<<val_<<"  ";

}}}
t1_minus_a2_hop<<endl;
}}}


//Nearest plusa1 neigbour-unitcell hopping (it includes only Honeycomb NN, NNN, NNNN)
//NNNN is not included by hand
/*string t1_minus_a1_plus_a2_hop_file="t1_minus_a1_plus_a2_hop.txt";
ofstream t1_minus_a1_plus_a2_hop(t1_minus_a1_plus_a2_hop_file.c_str());

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<2;orb1++){
for(int atom1=0;atom1<2;atom1++){

for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<2;orb2++){
for(int atom2=0;atom2<2;atom2++){

dis_1x=atom1*(1.0/(sqrt(3.0)*2.0));
dis_1y=atom1*(1.0/(2.0));

dis_2x=0.0 + atom2*(1.0/(sqrt(3.0)*2.0));
dis_2y= -1.0  + atom2*(1.0/(2.0));

dis_=sqrt( (dis_1x-dis_2x)*(dis_1x-dis_2x)  + (dis_1y-dis_2y)*(dis_1y-dis_2y)  );

if(dis_<Dis_max){
val_= abs(abs(1.0*spin1-1.0*spin2)-1.0)*Thop[2*atom1+orb1][2*atom2+orb2][i1_center+1][i2_center-1];
}
else{
val_=0.0;
}

if(abs(val_)<0.00001){
val_=0.0;
}


t1_minus_a1_plus_a2_hop<<val_<<"  ";
//t1_minus_a1_plus_a2_hop<<dis_<<"  ";

}}}
t1_minus_a1_plus_a2_hop<<endl;
}}}

*/

string t1_plus_a1_minus_a2_hop_file="t1_plus_a1_minus_a2_hop.txt";
ofstream t1_plus_a1_minus_a2_hop(t1_plus_a1_minus_a2_hop_file.c_str());

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<2;orb1++){
for(int atom1=0;atom1<2;atom1++){

for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<2;orb2++){
for(int atom2=0;atom2<2;atom2++){

dis_1x=atom1*(1.0/(sqrt(3.0)*2.0));
dis_1y=atom1*(1.0/(2.0));

dis_2x=0.0 + atom2*(1.0/(sqrt(3.0)*2.0));
dis_2y= 1.0  + atom2*(1.0/(2.0));

dis_=sqrt( (dis_1x-dis_2x)*(dis_1x-dis_2x)  + (dis_1y-dis_2y)*(dis_1y-dis_2y)  );

if(dis_<Dis_max){
val_= abs(abs(1.0*spin1-1.0*spin2)-1.0)*Thop[2*atom1+orb1][2*atom2+orb2][i1_center-1][i2_center+1];
}
else{
val_=0.0;
}

if(abs(val_)<0.00001){
val_=0.0;
}


t1_plus_a1_minus_a2_hop<<val_<<"  ";
//t1_minus_a1_plus_a2_hop<<dis_<<"  ";

}}}
t1_plus_a1_minus_a2_hop<<endl;
}}}






string t1_plus_a1_minus_2a2_hop_file="t1_plus_a1_minus_2a2_hop.txt";
ofstream t1_plus_a1_minus_2a2_hop(t1_plus_a1_minus_2a2_hop_file.c_str());

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<2;orb1++){
for(int atom1=0;atom1<2;atom1++){

for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<2;orb2++){
for(int atom2=0;atom2<2;atom2++){

dis_1x=atom1*(1.0/(sqrt(3.0)*2.0));
dis_1y=atom1*(1.0/(2.0));

dis_2x= sqrt(3.0)/2.0 + atom2*(1.0/(sqrt(3.0)*2.0));
dis_2y= 1.5  + atom2*(1.0/(2.0));

dis_=sqrt( (dis_1x-dis_2x)*(dis_1x-dis_2x)  + (dis_1y-dis_2y)*(dis_1y-dis_2y)  );

if(dis_<Dis_max){
val_= abs(abs(1.0*spin1-1.0*spin2)-1.0)*Thop[2*atom1+orb1][2*atom2+orb2][i1_center-1][i2_center+2];
}
else{
val_=0.0;
}

if(abs(val_)<0.00001){
val_=0.0;
}


t1_plus_a1_minus_2a2_hop<<val_<<"  ";
//t1_minus_a1_plus_a2_hop<<dis_<<"  ";

}}}
t1_plus_a1_minus_2a2_hop<<endl;
}}}













}
