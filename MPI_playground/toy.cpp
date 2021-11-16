

//PERFORM_REQUESTS

//-manda asincronicamente
//-waitany a ver cual llego y recibe sincronicamente


//ANSWER MSGS

//-recibe cual llego asincronicamente
//-manda sincronicamente

#include <iostream>
#include "mpi.h"
#include <omp.h>
#include <list>
#include <set>
#include <vector>
#include <queue>

#include <stdio.h>
#include <time.h>   /* Needed for struct timespec */


int mssleep(long miliseconds)
{
    struct timespec rem;
    struct timespec req= {
            (int)(miliseconds / 1000),     /* secs (Must be Non-Negative) */
            (miliseconds % 1000) * 1000000 /* nano (Must be in range of 0 to 999999999) */
    };

    return nanosleep(&req , &rem);
}


const int L = 10;


void fun2(int &ready){
	int world_size;
	int world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    printf("Hi! (fun2)\n");
	// Build
	std::queue<int> Q;
	MPI_Request R[L];
	int buffer[L];
	int response[L];
	for (int i=0; i<L; ++i){
		Q.push(i);
	}
	int i=0;
	while (!Q.empty()){	
		i = Q.front();
		buffer[i] = i;
		Q.pop();
		MPI_Ssend(&buffer[i], 1, MPI_INT,(int) ( (world_rank + 1) % world_size), 25001, MPI_COMM_WORLD);
	    printf("Sent! %d\n", buffer[i]);
		MPI_Recv(&response[i], 1, MPI_INT, (int) ( (world_rank + 1) % world_size), i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    printf("Recieved! %d\n", response[i]);
	}
	int newval=1;

	MPI_Ssend(&newval, 1, MPI_INT,(int) ( (world_rank + 1) % world_size), 9999, MPI_COMM_WORLD);
    MPI_Recv(&newval, 1, MPI_INT, (int) ( (world_rank + 1) % world_size), 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	printf("ENDED!");
#pragma omp atomic write 
	ready = 1;
	return;
}


void fun1(int &ready){
	int world_size;
	int world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	printf("Hi! (fun1)\n");
	// Build
	int V[world_size];
	int bigger_V[L];
	MPI_Request R[world_size];
	int buffer[L];
	int answers[L];
	int j = 0;
	int atomic_int;
	int statusflag = 0;
	bool first = true;

//	for (int i=0; i<L; ++i) {
//        MPI_Recv(&buffer[i], 1, MPI_INT, (int) ( (world_rank + 1) % world_size), 25001, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        answers[i] = i*5;
//        MPI_Ssend(&answers[i], 1, MPI_INT,(int) ( (world_rank + 1) % world_size), i, MPI_COMM_WORLD);
//    }

    int i = 0;
    int flag=0;
	atomic_int = 0;
	while (atomic_int != 1){
	    MPI_Status status;
	    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
	    int counter = 0;
        while ((flag != 1) && (counter < 20)){
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
            mssleep(80);
            counter ++;
        }
        if (counter <20) {
            if (status.MPI_TAG == 999){
                printf("Sending this special guy ;-)\n");
                int local = 1;
                MPI_Recv(&local, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Ssend(&local, 1, MPI_INT, status.MPI_SOURCE, 1000, MPI_COMM_WORLD);
                --i;
            } else {
                printf("Tag and source for Proc %d are %d and %d", world_rank, status.MPI_TAG, status.MPI_SOURCE);
                MPI_Recv(&buffer[i], 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                printf("Recieved! %d\n", buffer[i]);
                answers[i] = i * 5;
                MPI_Ssend(&answers[i], 1, MPI_INT, status.MPI_SOURCE, i, MPI_COMM_WORLD);
                printf("Sent! %d\n", answers[i]);
            }
        }
#pragma omp atomic read
        atomic_int = ready;
        ++i;
	}

	return;
}




int main(int argc, char** argv){
	int ready = 0;
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

#pragma omp parallel num_threads(2)
{
long MYNUM = omp_get_thread_num();	

	if (MYNUM == 0){
		fun1(ready);
	}
	else {
		fun2(ready);
	}

}
MPI_Finalize();
}
