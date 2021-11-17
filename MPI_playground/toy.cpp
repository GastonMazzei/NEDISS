

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
		MPI_Ssend(&buffer[i], 1, MPI_INT,(int) ( (world_rank + 1) % world_size), 2501, MPI_COMM_WORLD);
	    printf("Sent! now about to recieve :-) %d\n", buffer[i]);
		MPI_Recv(&response[i], 1, MPI_INT, (int) ( (world_rank + 1) % world_size), i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    printf("Recieved! %d\n", response[i]);
	}
//	int newval=1;
//	printf("Finished the loop");
//	MPI_Ssend(&newval, 1, MPI_INT,(int) ( (world_rank + 1) % world_size), 9999, MPI_COMM_WORLD);
//    	printf("finished 1 of 2 extra loops");
//	MPI_Recv(&newval, 1, MPI_INT, (int) ( (world_rank + 1) % world_size), 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//	printf("ENDED!");
//#pragma omp atomic write
//	ready = 1;
	return;
}


void fun1(int &ready){
	int world_size;
	int world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	printf("Hi! (fun1)\n");
	std::vector<int> buffer;
	std::vector<int> answers;
    int t=0, TIMEOUT = 20, DT=1;
    int flag = 0;
    int i=0;
    MPI_Request R;
    int TRIES=0;
    while (TRIES < 3){
        MPI_Status status;
        if ((flag == 1) || (i==0)) {
            R = MPI_Request();
            buffer.push_back(0);
            MPI_Irecv(&buffer[i], 1, MPI_INT,
                      MPI_ANY_SOURCE,
                      2501, MPI_COMM_WORLD, &R);
        }
        MPI_Request_get_status(R, &flag, &status);
        while ((flag != 1) && (t<TIMEOUT)){
            MPI_Request_get_status(R, &flag, &status);
            ++t;
            mssleep(DT);
            //printf("what  is happening here");
        }
        if (t>=TIMEOUT) ++TRIES;
        t=0;
        if (flag == 1){
            answers.push_back(0);
            answers[i] = i*5;
            MPI_Ssend(&answers[i], 1, MPI_INT, status.MPI_SOURCE, buffer[i], MPI_COMM_WORLD);
            ++i;
        }
    }
	return;
}




int main(int argc, char** argv){
	int ready = 0;
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

#pragma omp parallel num_threads(5)
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
