//
// Created by m4zz31 on 31/10/21.
//

#ifndef CPPPROJCT_ADEQUATE_SYNCHRONIZATION_H
#define CPPPROJCT_ADEQUATE_SYNCHRONIZATION_H

#include "../GraphClasses/GeneralGraph.h"
#include "msleep.h"

template <int T>
void adsync_synchronization_barrier(std::string detail, Graph &g){
    mssleep(T);
    synchronize(g.process_group());
    boost::mpi::communicator().barrier();
    MPI_Barrier(MPI_COMM_WORLD);
    if (process_id(g.process_group()) == 0) {
        std::cout << "[info] " << detail << " informs: synchronization has been successful!" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    mssleep(T);
}


template <int T>
void adsync_message_barrier(std::string message, Graph &g){
    mssleep(T);
    boost::mpi::communicator().barrier();
    MPI_Barrier(MPI_COMM_WORLD);
    if (process_id(g.process_group()) == 0) {
        std::cout << message << " informs: barrier has been successfully applied!" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    mssleep(T);
}

template <int T>
void adsync_barrier(){
    mssleep(T);
    boost::mpi::communicator().barrier();
    MPI_Barrier(MPI_COMM_WORLD);
    mssleep(T);
}

template <int T>
void adsync_message(std::string message, Graph &g){
    mssleep(T);
    if (process_id(g.process_group()) == 0) {
        std::cout << message << std::endl;
    }
    mssleep(T);
}



#endif //CPPPROJCT_ADEQUATE_SYNCHRONIZATION_H
