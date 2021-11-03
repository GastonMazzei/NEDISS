//
// Created by m4zz31 on 31/10/21.
//

#include "adequate_synchronization.h"

void adsync_synchronization_barrier(std::string detail, Graph &g){
    synchronize(g.process_group());
    boost::mpi::communicator().barrier();
    MPI_Barrier(MPI_COMM_WORLD);
    if (process_id(g.process_group()) == 0) {
        std::cout << "[info] " << detail << " informs: synchronization has been successful!" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void adsync_message_barrier(std::string message, Graph &g){
    boost::mpi::communicator().barrier();
    MPI_Barrier(MPI_COMM_WORLD);
    if (process_id(g.process_group()) == 0) {
        std::cout << message << " informs: barrier has been successfully applied!" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void adsync_barrier(){
    boost::mpi::communicator().barrier();
    MPI_Barrier(MPI_COMM_WORLD);
}


void adsync_message(std::string message, Graph &g){
    if (process_id(g.process_group()) == 0) {
        std::cout << message << std::endl;
    }
}


