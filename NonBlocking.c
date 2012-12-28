// test programme for non-blocking communication
// (c) 2012 Carsten Urbach free under GPL
// compile on supermuc with intel compiler
//
// for normal communication (no overlapping):
// mpicc -O3 -axAVX -std=gnu99 -fopenmp NonBlocking.c
//
// for interleaved computation and communictaion:
// mpicc -O3 -axAVX -std=gnu99 -D_INTERLEAVE_ -fopenmp NonBlocking.c
// 
// for communictaion switched off:
// mpicc -O3 -axAVX -std=gnu99 -D_NOCOMM_ -fopenmp NonBlocking.c
//
// and run with 16=dims[0]*dims[1]*dims[2]*dims[3] MPI processes
// or adjust parameters dims[] accordingly

#include<stdlib.h>
#include<stdio.h>
#include <sys/types.h>
#include <stdint.h>

#include <mpi.h>
#include <omp.h>

#define SEND_BUFFER_ALIGNMENT   32
#define RECV_BUFFER_ALIGNMENT   32
#define MAX_MESSAGE_SIZE       32768

// we have four directions and forward/backward
#define NUM_DIRS               8
// number of loops to average over
#define N_LOOPS                10000
// MPI local workload
#define N_LOCAL                64

// Allocate static memory for send and receive buffers
char sendBufMemory[NUM_DIRS * MAX_MESSAGE_SIZE+ SEND_BUFFER_ALIGNMENT];
char recvBufMemory[NUM_DIRS * MAX_MESSAGE_SIZE+ SEND_BUFFER_ALIGNMENT];
// pointers to send and receive buffers
char * recvBuffers;
char * sendBuffers;

// message size and offsets in bytes
uint64_t messageSizes[NUM_DIRS];
uint64_t roffsets[NUM_DIRS], soffsets[NUM_DIRS];
uint64_t totalMessageSize;

// here come the MPI variables
int g_proc_id, g_nproc, g_cart_id, g_proc_coords[4];
int g_nproc_t, g_nproc_x, g_nproc_y, g_nproc_z;
MPI_Comm g_cart_grid;
int g_nb_up[4], g_nb_dn[4];
int g_mpi_prov = 0;

int main(int argc, char **argv) {
  int rc;
  int g_proc_id;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &g_mpi_prov);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);
  if(g_proc_id == 0){
    printf("Programme to test non-blocking communication\n");
    printf("2012 (c) Carsten Urbach free under GPL\n");
    printf("MPI implementation provided thread support = %d\n", g_mpi_prov);
    printf("We are testing with:");
#ifdef _INTERLEAVE_
    printf("Interleaving communication and computation\n");
#elif defined _NOCOMM_
    printf("Communication switched off\n");
#else
    printf("Not overlapping communication and computation\n");
#endif
  }
  int periods[] = {1,1,1,1};
  int namelen;
  int ndims = 4;
  int dims[4];
  int nalldims = 4;
  MPI_Request requests[16];
  MPI_Status hstatus[16];
  dims[0] = 2;
  dims[1] = 2;
  dims[2] = 2;
  dims[3] = 2;

  MPI_Comm_size(MPI_COMM_WORLD, &g_nproc);
  MPI_Get_processor_name(processor_name, &namelen);
  MPI_Dims_create(g_nproc, nalldims, dims);
  if(g_proc_id == 0){
    printf("# Creating the following cartesian grid for a %d dimensional parallelisation:\n# %d x %d x %d x %d\n"
            , ndims, dims[0], dims[1], dims[2], dims[3]);
  }

  g_nproc_t = dims[0];
  g_nproc_x = dims[1];
  g_nproc_y = dims[2];
  g_nproc_z = dims[3];

  // create an MPI cartesian grid
  int reorder = 1;
  MPI_Cart_create(MPI_COMM_WORLD, nalldims, dims, periods, reorder, &g_cart_grid);
  MPI_Comm_rank(g_cart_grid, &g_cart_id);
  MPI_Cart_coords(g_cart_grid, g_cart_id, nalldims, g_proc_coords);
  
  // obtain neighbouring MPI process ids
  for(int i = 0; i < ndims; i++) {
    MPI_Cart_shift(g_cart_grid, i, 1, &g_nb_dn[i], &g_nb_up[i]);
  }
  
  // determine offsets in send and receive buffers
  totalMessageSize = 0;
  for(int i = 0; i < NUM_DIRS; i ++) {
    messageSizes[i] = MAX_MESSAGE_SIZE;
    soffsets[i] = totalMessageSize;
    roffsets[i] = soffsets[i];
    totalMessageSize += messageSizes[i];
  }

  recvBuffers = (char *)(((uint64_t)recvBufMemory+RECV_BUFFER_ALIGNMENT)&~(RECV_BUFFER_ALIGNMENT-1));    
  sendBuffers = (char *)(((uint64_t)sendBufMemory+SEND_BUFFER_ALIGNMENT)&~(SEND_BUFFER_ALIGNMENT-1));

  double startTime=0, endTime = 0, tmp;
  double waitTime = 0., iTime = 0.;
  startTime = MPI_Wtime();

  MPI_Barrier(MPI_COMM_WORLD);

  // Fill send buffer
  for(int n = 0; n < totalMessageSize; n += sizeof(double)) {
    *(double*)&sendBuffers[n] = (double) g_cart_id;
  }

  double s = 0., r = 0., sum=0.;
  // perform some loops to get some averaging
  for(int l = 0; l < N_LOOPS; l++) {
    
#ifndef _NOCOMM_
    tmp = MPI_Wtime();
    // here we init the send and receive operations...
    // Irecv first
    for (int j = 0; j < ndims; j++) {
      MPI_Irecv((void*)(&recvBuffers[roffsets[2*j]]), messageSizes[2*j], MPI_CHAR, 
      	g_nb_dn[j], 80+2*j, g_cart_grid, &requests[4*j+1]);
      MPI_Irecv((void*)(&recvBuffers[roffsets[2*j+1]]), messageSizes[2*j+1], MPI_CHAR, 
      	g_nb_up[j], 80+2*j+1, g_cart_grid, &requests[4*j+3]);
    }
    // now the Isend operations
    for (int j = 0; j < ndims; j++) {
      MPI_Isend((void*)(&sendBuffers[soffsets[2*j]]), messageSizes[2*j], MPI_CHAR, 
      	g_nb_up[j], 80+2*j, g_cart_grid, &requests[4*j+0]);
      MPI_Isend((void*)(&sendBuffers[soffsets[2*j+1]]), messageSizes[2*j+1], MPI_CHAR, 
      	g_nb_dn[j], 80+2*j+1, g_cart_grid, &requests[4*j+2]);
    }
    iTime += (MPI_Wtime() - tmp);

#  ifndef _INTERLEAVE_
    // and poll for completion in case we don't overlap comm and comp
    tmp = MPI_Wtime();
    MPI_Waitall(4*ndims, requests, hstatus); 
    waitTime += (MPI_Wtime() - tmp);
#  endif
#endif

#pragma omp parallel reduction(+: s)
    {
      s=0.;
      // do some computation to hide communication
#pragma omp for
      for(int m = 0; m < N_LOCAL; m++) {
    	for(int n = 0; n < totalMessageSize; n+=sizeof(double)) {
    	  s += *(double*)&sendBuffers[n];
    	}
      }
    }
    sum += s;

#ifndef _NOCOMM_
#  ifdef _INTERLEAVE_
#    pragma omp single
    {
      // and poll for completion in case we overlap comm and comp
      tmp = MPI_Wtime();
      MPI_Waitall(4*ndims, requests, hstatus); 
      waitTime += (MPI_Wtime() - tmp);
    }
#  endif
#endif

#pragma omp parallel reduction(+: r)
    {
      r=0.;
      //do some more computation
#pragma omp for
      for(int m = 0; m < N_LOCAL; m++) {
	for(int n = 0; n < totalMessageSize; n+=sizeof(double)) {
	  r += *(double*)&recvBufMemory[n];
	}
      }
    }
    sum += r;
  }
  
  endTime = MPI_Wtime();

#ifndef _NOCOMM_
  // here we check for correctness of communications
  if(g_nb_up[0] != (int)*(double*)&recvBuffers[soffsets[0]] ||
     g_nb_dn[0] != (int)*(double*)&recvBuffers[soffsets[1]] ||
     g_nb_up[1] != (int)*(double*)&recvBuffers[soffsets[2]] ||
     g_nb_dn[1] != (int)*(double*)&recvBuffers[soffsets[3]] ||
     g_nb_up[2] != (int)*(double*)&recvBuffers[soffsets[4]] ||
     g_nb_dn[2] != (int)*(double*)&recvBuffers[soffsets[5]] ||
     g_nb_up[3] != (int)*(double*)&recvBuffers[soffsets[6]] ||
     g_nb_dn[3] != (int)*(double*)&recvBuffers[soffsets[7]]) {
    printf("neighbours wrong for %d !\n", g_cart_id);
    printf("t+ %d %d\n", g_nb_up[0], (int)*(double*)&recvBuffers[soffsets[0]]);
    printf("t- %d %d\n", g_nb_dn[0], (int)*(double*)&recvBuffers[soffsets[1]]);
    printf("x+ %d %d\n", g_nb_up[1], (int)*(double*)&recvBuffers[soffsets[2]]);
    printf("x- %d %d\n", g_nb_dn[1], (int)*(double*)&recvBuffers[soffsets[3]]);
    printf("y+ %d %d\n", g_nb_up[2], (int)*(double*)&recvBuffers[soffsets[4]]);
    printf("y- %d %d\n", g_nb_dn[2], (int)*(double*)&recvBuffers[soffsets[5]]);
    printf("z+ %d %d\n", g_nb_up[3], (int)*(double*)&recvBuffers[soffsets[6]]);
    printf("z- %d %d\n", g_nb_dn[3], (int)*(double*)&recvBuffers[soffsets[7]]);
  }
#endif
  if(g_proc_id == 0) {
    printf("total time per loop is %e seconds\n", (endTime - startTime)/(double)N_LOOPS);
    printf("total time for comm init per loop is %e seconds\n", iTime/(double)N_LOOPS);
    printf("total time for wait per loop is %e seconds\n", waitTime/(double)N_LOOPS);
  }
  if(g_proc_id == -1) printf("res for %d is %e\n",g_proc_id, sum/(double)N_LOOPS);

  MPI_Finalize();
  return(0);
}



