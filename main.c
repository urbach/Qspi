#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <stdint.h>
#include <malloc.h>
#include <assert.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>

//////////////////////////////////////////
// Basic SPI and HWI includes
//////////////////////////////////////////
#include <hwi/include/bqc/A2_core.h>
#include <hwi/include/bqc/A2_inlines.h>
#include <hwi/include/bqc/MU_PacketCommon.h>
#include <firmware/include/personality.h>
#include <spi/include/mu/Descriptor.h>
#include <spi/include/mu/Descriptor_inlines.h>
#include <spi/include/mu/InjFifo.h>
#include <spi/include/mu/RecFifo.h>
#include <spi/include/mu/Addressing.h>
#include <spi/include/mu/Addressing_inlines.h>
#include <spi/include/mu/GIBarrier.h>
#include <spi/include/kernel/MU.h>
#include <spi/include/kernel/process.h>
#include <spi/include/kernel/location.h>
#include <spi/include/kernel/collective.h>


#define ALIGN __attribute__((__aligned__(64)))

// Injection Memory FIFO Descriptor
MUHWI_Descriptor_t t_fifo[2] ALIGN;
MUHWI_Descriptor_t x_fifo[2] ALIGN;
MUHWI_Descriptor_t y_fifo[2] ALIGN;
MUHWI_Descriptor_t z_fifo[2] ALIGN;

// Injection Memory FIFO Descriptor Information Structures
MUSPI_Pt2PtMemoryFIFODescriptorInfo_t t_fifo_info[2];
MUSPI_Pt2PtMemoryFIFODescriptorInfo_t x_fifo_info[2];
MUSPI_Pt2PtMemoryFIFODescriptorInfo_t y_fifo_info[2];
MUSPI_Pt2PtMemoryFIFODescriptorInfo_t z_fifo_info[2];

uint64_t * tsend_buf, trecv_buf;
uint64_t * xsend_buf, xrecv_buf;
uint64_t * ysend_buf, yrecv_buf;
uint64_t * zsend_buf, zrecv_buf;

int g_proc_id, g_nproc, g_cart_id, g_proc_coords[4], g_nb_list[8];
MPI_Comm g_cart_grid;
int g_nb_x_up, g_nb_x_dn;
int g_nb_y_up, g_nb_y_dn;
int g_nb_t_up, g_nb_t_dn;
int g_nb_z_up, g_nb_z_dn;


int main (int argc,char *argv[]) {
  Personality_t pers;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  Kernel_GetPersonality(&pers, sizeof(pers));
  int g_proc_id;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);

  int dims[4];
  int nalldims = 4;
  dims[0] = 4;
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
  
  MPI_Cart_shift(g_cart_grid, 0, 1, &g_nb_t_dn, &g_nb_t_up);
  g_nb_list[0] = g_nb_t_up;  
  g_nb_list[1] = g_nb_t_dn;
  MPI_Cart_shift(g_cart_grid, 1, 1, &g_nb_x_dn, &g_nb_x_up);
  g_nb_list[2] = g_nb_x_up;  
  g_nb_list[3] = g_nb_x_dn;
  MPI_Cart_shift(g_cart_grid, 2, 1, &g_nb_y_dn, &g_nb_y_up);
  g_nb_list[4] = g_nb_y_up;  
  g_nb_list[5] = g_nb_y_dn;
  MPI_Cart_shift(g_cart_grid, 3, 1, &g_nb_z_dn, &g_nb_z_up);
  g_nb_list[6] = g_nb_z_up;  
  g_nb_list[7] = g_nb_z_dn;


  int my_a = 0, my_b = 0, my_c = 0, my_d = 0, my_e = 0;
  my_a = pers.Network_Config.Acoord;
  my_b = pers.Network_Config.Bcoord;
  my_c = pers.Network_Config.Ccoord;
  my_d = pers.Network_Config.Dcoord;
  my_e = pers.Network_Config.Ecoord;

  fprintf(stdout,"# Process %d of %d on %s: cart_id %d, coordinates (%d %d %d %d)\n# Process %d has personality %d %d %d %d %d\n",
	  g_proc_id, g_nproc, processor_name, g_cart_id, 
	  g_proc_coords[0], g_proc_coords[1], g_proc_coords[2], g_proc_coords[3],
	  g_proc_id, my_a, my_b, my_c, my_d, my_e);
  fflush(stdout);


  MPI_Finalize();
  return(0);
}
