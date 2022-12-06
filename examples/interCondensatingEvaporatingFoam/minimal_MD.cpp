#include "mpi.h"
#include <iostream>

#include "cpl.h"
#include "CPL_ndArray.h"

using namespace std;

int main() {

    bool flag;
    int MD_realm = 2;
    MPI_Comm MD_COMM, CART_COMM;
    CPL::ndArray<double> send_array, recv_array;

    MPI_Init(NULL, NULL); 
    CPL::init(MD_realm, MD_COMM);

    int npxyz[3] = {1, 1, 1}; int periods[3] = {1, 1, 1}; 
    MPI_Cart_create(MD_COMM, 3, npxyz, periods, 1, &CART_COMM);
    double xyzL[3] = {0.002, 0.005, 0.01}; 
    double xyz_orig[3] = {0.0, 0.0, 0.0};
    CPL::setup_md(CART_COMM, xyzL, xyz_orig);
    CPL::get_arrays(&recv_array, 3, &send_array, 4, MD_realm);

    std::cout << "recv_array.shape " << recv_array.shape(0) << " " 
                                     << recv_array.shape(1) << " " 
                                     << recv_array.shape(2) << " " 
                                     << recv_array.shape(3) << std::endl;
	for (int time = 0; time < 1501; time++) {
        flag = CPL::recv(&recv_array);
        for (int i=0; i<recv_array.shape(1); i++){
        for (int k=0; k<recv_array.shape(3); k++){
            std::cout << "MD recv_array " << i << " " << k 
                << " " << recv_array(0,i,0,k) 
                << " " << recv_array(1,i,0,k) 
                << " " << recv_array(2,i,0,k) << std::endl;
        }}

        for (int i=0; i<send_array.shape(1); i++){
        for (int k=0; k<send_array.shape(3); k++){
            send_array(3,i,0,k) = 1.0;
            send_array(0,i,0,k) = -1.0;
            //send_array(0,i,0,k) = i/float(send_array.shape(1)) + k/float(send_array.shape(3));
            std::cout << "MD send_array " << i << " " << k 
                   << " " << send_array(0,i,0,k)
                   << " " << send_array(1,i,0,k)
                   << " " << send_array(2,i,0,k)
                   << " " << send_array(3,i,0,k) << std::endl;
        }}
        flag = CPL::send(&send_array);
        std::cout << "MD " << time << std::endl;
    }

   // Release all coupler comms 
    CPL::finalize();
    MPI_Finalize();
   
}
