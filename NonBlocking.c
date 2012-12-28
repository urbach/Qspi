#include<stdlib.h>
#include<stdio.h>
#include <sys/types.h>
#include <stdint.h>

#include <mpi.h>
#include <omp.h>

#define SEND_BUFFER_ALIGNMENT   128
#define RECV_BUFFER_ALIGNMENT   128
#define MAX_MESSAGE_SIZE       32768

long long messageSizeInBytes = MAX_MESSAGE_SIZE;


// we have four directions and forward/backward
#define NUM_DIRS               8
#define N_LOOPS                10000

// Allocate static memory for send and receive buffers
char sendBufMemory[NUM_DIRS * MAX_MESSAGE_SIZE+ SEND_BUFFER_ALIGNMENT];
char recvBufMemory[NUM_DIRS * MAX_MESSAGE_SIZE+ SEND_BUFFER_ALIGNMENT];
// pointers to send and receive buffers
char * recvBuffers;
char * sendBuffers;

// in bytes
uint64_t messageSizes[NUM_DIRS];
uint64_t roffsets[NUM_DIRS], soffsets[NUM_DIRS];
uint64_t totalMessageSize;

// here come the MPI variables
int g_proc_id, g_nproc, g_cart_id, g_proc_coords[4], g_nb_list[8];
int g_nproc_t, g_nproc_x, g_nproc_y, g_nproc_z;
MPI_Comm g_cart_grid;
int g_nb_up[4], g_nb_dn[4];
int g_mpi_prov = 0;



int main(int argc, char **argv) {
  int rc;
  int g_proc_id;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &g_mpi_prov);
  //MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);
  if(g_proc_id == 0){
    printf("provided thread support = %d\n", g_mpi_prov);
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

  int reorder = 1;
  MPI_Cart_create(MPI_COMM_WORLD, nalldims, dims, periods, reorder, &g_cart_grid);
  MPI_Comm_rank(g_cart_grid, &g_cart_id);
  MPI_Cart_coords(g_cart_grid, g_cart_id, nalldims, g_proc_coords);
  
  for(int i = 0; i < ndims; i++) {
    MPI_Cart_shift(g_cart_grid, i, 1, &g_nb_dn[i], &g_nb_up[i]);
    g_nb_list[2*i  ] = g_nb_up[i];  
    g_nb_list[2*i+1] = g_nb_dn[i];
  }

  totalMessageSize = 0;
  for(int i = 0; i < NUM_DIRS; i ++) {
    messageSizes[i] = MAX_MESSAGE_SIZE;
    //if(i == 1 || i == 0) messageSizes[i] = MAX_MESSAGE_SIZE/2;
    soffsets[i] = totalMessageSize;
    roffsets[i] = soffsets[i];
    totalMessageSize += messageSizes[i];
  }
  for(int i = 0; i < 0*NUM_DIRS; i++) {
    if(i%2 == 0) {
      roffsets[i] = soffsets[i] + messageSizes[i];
    }
    else {
      roffsets[i] = soffsets[i] - messageSizes[i-1];
    }
  }

  recvBuffers = (char *)(((uint64_t)recvBufMemory+RECV_BUFFER_ALIGNMENT)&~(RECV_BUFFER_ALIGNMENT-1));    
  sendBuffers = (char *)(((uint64_t)sendBufMemory+SEND_BUFFER_ALIGNMENT)&~(SEND_BUFFER_ALIGNMENT-1));

  uint64_t totalCycles=0;
  double startTime=0, endTime = 0;
  startTime = MPI_Wtime();

  MPI_Barrier(MPI_COMM_WORLD);

  // Fill send buffer
  for(int n = 0; n < totalMessageSize; n += sizeof(double)) {
    *(double*)&sendBuffers[n] = (double) g_cart_id;
  }

  double s = 0., r = 0.;
  for(int l = 0; l < N_LOOPS; l++) {
    
    // here we init the send and receive operations...
    for (int j = 0; j < ndims; j++) {

      MPI_Isend((void*)(&sendBuffers[soffsets[2*j]]), messageSizes[2*j], MPI_CHAR, 
      	g_nb_up[j], 80+2*j, g_cart_grid, &requests[4*j+0]);
      MPI_Irecv((void*)(&recvBuffers[roffsets[2*j]]), messageSizes[2*j], MPI_CHAR, 
      	g_nb_dn[j], 80+2*j, g_cart_grid, &requests[4*j+1]);
      MPI_Isend((void*)(&sendBuffers[soffsets[2*j+1]]), messageSizes[2*j+1], MPI_CHAR, 
      	g_nb_dn[j], 80+2*j+1, g_cart_grid, &requests[4*j+2]);
      MPI_Irecv((void*)(&recvBuffers[roffsets[2*j+1]]), messageSizes[2*j+1], MPI_CHAR, 
      	g_nb_up[j], 80+2*j+1, g_cart_grid, &requests[4*j+3]);
    }

#ifndef _INTERLEAVE_
    MPI_Waitall(4*ndims, requests, hstatus); 
#endif

#pragma omp parallel reduction(+: s)
    {
      
      // do some computation to hide communication
#pragma omp for
      for(int m = 0; m < 4; m++) {
    	for(int n = 0; n < totalMessageSize; n+=sizeof(double)) {
    	  s += *(double*)&sendBuffers[n];
    	}
      }
    }

#ifdef _INTERLEAVE_
#  pragma omp single
    {
      MPI_Waitall(4*ndims, requests, hstatus); 
    }
#endif

#pragma omp parallel reduction(+: r)
    {
      //do some computation not hidding communication
#pragma omp for
      for(int m = 0; m < 4; m++) {
	for(int n = 0; n < totalMessageSize; n+=sizeof(double)) {
	  r += *(double*)&recvBufMemory[n];
	}
      }
    }
  }
  
  endTime = MPI_Wtime();

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
  if(g_proc_id == 0) {
    printf("total time per loop is %e seconds\n", (endTime - startTime)/(double)N_LOOPS);
  }
  if(g_proc_id == -1) printf("res for %d is %e\n",g_proc_id, s+r);

  MPI_Finalize();
  return(0);
}



