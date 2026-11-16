// DARE.cpp
// Ricatti equation implementation (DARE)
// MvM April 2019 for teaching

// Discrete-time Algebraic Ricatti Equation, recursion backwards in time
// finite horizon LQR (optimal state feedback) to choose braking at each step
// controller for lander. gives optimal K[0] K[1]... for control u=-Kx

#include <iostream>
#include <stdexcept>

#include "DARE.h"


// Q and R are weighting matrices for state and control effort
DARE::DARE(
    const Eigen::MatrixXd& A,
    const Eigen::MatrixXd& B,
    const Eigen::MatrixXd& Q,
    const Eigen::MatrixXd& R,
    int horizon
    )
    : A(A), B(B), Q(Q),R(R)
{
  int  m,k;
  state_n=A.rows();
  m=B.rows();
  k=B.cols();

  N=horizon;

  initialized=false;

 // ARRAY FOR THE  GAINS
 K_array=(Eigen::MatrixXd*)malloc(N*sizeof(Eigen::MatrixXd));
 S_array=(Eigen::MatrixXd*)malloc(N*sizeof(Eigen::MatrixXd));
}

DARE::~DARE()
{
	free(K_array);
	free(S_array);
}



void DARE::init(void)
{
	int i;
        // 0 case  (we number them backward, but note K0 is to bin x from x(n-1) to x(n) !
        K_array[N-1]=  Eigen::MatrixXd(1,state_n);
        S_array[N-1]=  Eigen::MatrixXd(state_n,state_n);

	// cost to go is the terminal penalty
	K_array[N-1].setIdentity();
        S_array[N-1].setIdentity();
	

	for(i=N-2;i>=0;i--){
		K_array[i]= Eigen::MatrixXd(state_n,1);
		S_array[i]=Eigen::MatrixXd(state_n,state_n);

		// compute optimal gain for step i
		K_array[i] = (R+ B.transpose() * S_array[i+1] *B ).inverse()  *
				  B.transpose()* S_array[i+1]* A;
		// new cost-to-go, state and control cost + future cost
 		S_array[i]= Q+  K_array[i].transpose() *R * K_array[i] +	
			      (A-B*K_array[i]).transpose()*S_array[i+1]*(A-B*K_array[i]);
	}
	initialized=true;
}



void DARE::print(void)
{

	std::cout << "A= " << A <<std::endl;
	std::cout << "B= " << B <<std::endl;
	int i;
	for(i=0;i<N;i++){
		std::cout << " i= " <<i  <<std::endl;
		std::cout << "K= " <<K_array[i] <<std::endl;
		std::cout << "S= "<< S_array[i] <<std::endl;
		std::cout <<std::endl;
	}


}


