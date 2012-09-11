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

double * tsend_buf, * trecv_buf;
double * xsend_buf, * xrecv_buf;
double * ysend_buf, * yrecv_buf;
double * zsend_buf, * zrecv_buf;

int g_proc_id, g_nproc, g_cart_id, g_proc_coords[4], g_nb_list[8];
int g_nproc_t, g_nproc_x, g_nproc_y, g_nproc_z;
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

  int periods[] = {1,1,1,1};
  int namelen;
  int ndims = 4;
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


  int my_c[6] = {0,0,0,0,0,0};
  int nbtplus[6], nbtminus[6];
  int nbxplus[6], nbxminus[6];
  int nbyplus[6], nbyminus[6];
  int nbzplus[6], nbzminus[6];
  my_c[0] = pers.Network_Config.Acoord;
  my_c[1] = pers.Network_Config.Bcoord;
  my_c[2] = pers.Network_Config.Ccoord;
  my_c[3] = pers.Network_Config.Dcoord;
  my_c[4] = pers.Network_Config.Ecoord;
  my_c[5] = 0;//pers.Network_config.Tcoord;

  fprintf(stdout,"# MPI Process %d of %d on %s: cart_id %d, coordinates (%d %d %d %d)\n# MPI Process %d has personality %d %d %d %d %d %d\n",
	  g_proc_id, g_nproc, processor_name, g_cart_id, 
	  g_proc_coords[0], g_proc_coords[1], g_proc_coords[2], g_proc_coords[3],
	  g_proc_id, my_c[0], my_c[1], my_c[2], my_c[3], my_c[4], my_c[5]);
  fflush(stdout);

  MPI_Status mstatus;
  MPI_Sendrecv((void*)my_c, 6, MPI_INT, g_nb_t_up, 0, 
	       (void*)nbtminus, 6, MPI_INT, g_nb_t_dn, 0,
	       g_cart_grid, &mstatus);
  MPI_Sendrecv((void*)my_c, 6, MPI_INT, g_nb_t_dn, 1, 
		 (void*)nbtplus, 6, MPI_INT, g_nb_t_up, 1, 
		 g_cart_grid, &mstatus);

  MPI_Sendrecv((void*)my_c, 6, MPI_INT, g_nb_x_up, 2, 
		 (void*)nbxminus, 6, MPI_INT, g_nb_x_dn, 2, 
		 g_cart_grid, &mstatus);
  MPI_Sendrecv((void*)my_c, 6, MPI_INT, g_nb_x_dn, 3, 
		 (void*)nbxplus, 6, MPI_INT, g_nb_x_up, 3, 
		 g_cart_grid, &mstatus);

  MPI_Sendrecv((void*)my_c, 6, MPI_INT, g_nb_y_up, 4, 
		 (void*)nbyminus, 6, MPI_INT, g_nb_y_dn, 4, 
		 g_cart_grid, &mstatus);
  MPI_Sendrecv((void*)my_c, 6, MPI_INT, g_nb_y_dn, 5, 
		 (void*)nbyplus, 6, MPI_INT, g_nb_y_up, 5, 
		 g_cart_grid, &mstatus);

  MPI_Sendrecv((void*)my_c, 6, MPI_INT, g_nb_z_up, 6, 
		 (void*)nbzminus, 6, MPI_INT, g_nb_z_dn, 6, 
		 g_cart_grid, &mstatus);
  MPI_Sendrecv((void*)my_c, 6, MPI_INT, g_nb_z_dn, 7, 
		 (void*)nbzplus, 6, MPI_INT, g_nb_z_up, 7, 
		 g_cart_grid, &mstatus);

  if(g_proc_id == 0) {
    printf("my coords (%d %d %d %d)\n", g_proc_coords[0], g_proc_coords[1], g_proc_coords[2], g_proc_coords[3]);
    printf("my kernel %d %d %d %d %d\n", my_c[0], my_c[1], my_c[2], my_c[3], my_c[4]);
    printf("my tplus  %d %d %d %d %d\n", nbtplus[0], nbtplus[1], nbtplus[2], nbtplus[3], nbtplus[4]);
    printf("my tminus  %d %d %d %d %d\n", nbtminus[0], nbtminus[1], nbtminus[2], nbtminus[3], nbtminus[4]);
    printf("my xplus  %d %d %d %d %d\n", nbxplus[0], nbxplus[1], nbxplus[2], nbxplus[3], nbxplus[4]);
    printf("my xminus  %d %d %d %d %d\n", nbxminus[0], nbxminus[1], nbxminus[2], nbxminus[3], nbxminus[4]);
    printf("my yplus  %d %d %d %d %d\n", nbyplus[0], nbyplus[1], nbyplus[2], nbyplus[3], nbyplus[4]);
    printf("my yminus  %d %d %d %d %d\n", nbyminus[0], nbyminus[1], nbyminus[2], nbyminus[3], nbyminus[4]);
    printf("my zplus  %d %d %d %d %d\n", nbzplus[0], nbzplus[1], nbzplus[2], nbzplus[3], nbzplus[4]);
    printf("my zminus  %d %d %d %d %d\n", nbzminus[0], nbzminus[1], nbzminus[2], nbzminus[3], nbzminus[4]);
  }


  MPI_Finalize();
  return(0);
}
