#include <vector>
#include <complex>
using namespace std;
#ifndef tensors
#define tensors

typedef vector< complex<double> >  Mat_1_Complex_doub;
typedef vector<Mat_1_Complex_doub> Mat_2_Complex_doub;
typedef vector<Mat_2_Complex_doub> Mat_3_Complex_doub;
typedef vector<Mat_3_Complex_doub> Mat_4_Complex_doub;




typedef vector<double> Mat_1_doub;
typedef vector<Mat_1_doub> Mat_2_doub;
typedef vector<Mat_2_doub> Mat_3_doub;
typedef vector<Mat_3_doub> Mat_4_doub;
typedef vector<Mat_4_doub> Mat_5_doub;
typedef vector<Mat_5_doub> Mat_6_doub;
typedef vector<Mat_6_doub> Mat_7_doub;






typedef vector<int> Mat_1_int;
typedef vector<Mat_1_int> Mat_2_int;
typedef vector<Mat_2_int> Mat_3_int;
typedef vector<Mat_3_int> Mat_4_int;

typedef vector<string> Mat_1_string;
typedef vector<Mat_1_string> Mat_2_string;
typedef vector<Mat_2_string> Mat_3_string;
typedef vector<Mat_3_string> Mat_4_string;


typedef vector<double> Mat_1_real;
typedef vector<Mat_1_real> Mat_2_real;
typedef vector<Mat_2_real> Mat_3_real;
typedef vector<Mat_3_real> Mat_4_real;



complex<double> const one_complex(1,0);
complex<double> const iota_complex(0,1);
complex<double> const zero_complex(0,0);


struct TUPLE_2_INT{
	
	
	int first;
	int second;
	
	
	
};
	
struct pair_doubletype_int{
    double first;
    int second;
};
typedef vector<pair_doubletype_int> Mat_1_pair_DoubInt;
typedef vector<Mat_1_pair_DoubInt> Mat_2_pair_Doubint;


struct pair_real_int{
    double first;
    int second;
};
typedef vector<pair_real_int> Mat_1_pair_realInt;
typedef vector<Mat_1_pair_realInt> Mat_2_pair_realint;



struct pair_int{
    int first;
    int second;
};
typedef vector<pair_int> Mat_1_intpair;

class INT_CSR_MAT{
	

	public:
	Mat_1_int value;
	Mat_1_int columns;
	Mat_1_int row_ind;
	
	
};


class DOUBLE_CSR_MAT{
	
	
	
	public:

	Mat_1_doub value;
	Mat_1_int columns;
	Mat_1_int row_ind;

	
	
};


class DOUBLE_COO_MAT{
	
	
	
	public:

	Mat_1_doub value;
	Mat_1_int columns;
	Mat_1_int rows;
	int nrows;
	int ncols;

	
	
};


class COMPLEXDOUBLE_CSR_MAT{
	
	
	
    public:
    Mat_1_Complex_doub value;
    Mat_1_int columns;
    Mat_1_int row_ind;
	
	
};

//template<class Mat_1_Type>
//struct TMatrix_COO{

//        Mat_1_Type value;
//        Mat_1_int  rows;
//        Mat_1_int  columns;
//        int nrows;
//        int ncols;

//    };



struct Matrix_COO{
	
	    Mat_1_doub value;
		Mat_1_int  rows;
		Mat_1_int  columns;
		int nrows;
		int ncols;
	
	};

struct Matrix_COO_Complex{

        Mat_1_Complex_doub value;
        Mat_1_int  rows;
        Mat_1_int  columns;
        int nrows;
        int ncols;

    };





typedef vector<Matrix_COO> Hamiltonian_1_COO;
typedef vector<Hamiltonian_1_COO> Hamiltonian_2_COO;
typedef vector<Hamiltonian_2_COO> Hamiltonian_3_COO;
typedef vector<Hamiltonian_3_COO> Hamiltonian_4_COO;
typedef vector<Hamiltonian_4_COO> Hamiltonian_5_COO;
typedef vector<Hamiltonian_5_COO> Hamiltonian_6_COO;
typedef vector<Hamiltonian_6_COO> Hamiltonian_7_COO;
typedef vector<Hamiltonian_7_COO> Hamiltonian_8_COO;
typedef vector<Hamiltonian_8_COO> Hamiltonian_9_COO;


#endif
