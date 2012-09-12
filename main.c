#include "spi.h"
#include <mpi.h>

// The ping and pong functions
int ping (int a, int b, int c, int d, int e, int t, int bytes, int my_t);
int pong (int a, int b, int c, int d, int e, int t, int bytes, int my_t);

#define NUM_LOOPS   1000

#define MAX_MESSAGE_SIZE              8192            // Multiple of 8
#define REC_MEMORY_FIFO_SIZE          0x00000FFFFULL  // 64K bytes
#define INJ_MEMORY_FIFO_SIZE          0x00000FFFFULL  // 64K bytes

// Injection Memory FIFO Descriptor
MUHWI_Descriptor_t mu_iMemoryFifoDescriptor[2] __attribute__((__aligned__(64))) ;

// Injection Memory FIFO Descriptor Information Structures
MUSPI_Pt2PtMemoryFIFODescriptorInfo_t mu_iMemoryFifoDescriptorInfo[2];

/**
 * \brief Context for book-keeping on the receiver to process incoming
 * packets.
 */
typedef struct MsgContext {
  int    done;
  int    size;
  int    bytes_recvd;
  char  *buf;
  unsigned long long HWEndTime;
} MsgContext_t;

MsgContext_t   my_recv_context[2];  /** Store receive context */

/**
 * \brief Receive packet callback handler. SPIs can poll FIFOs and
 * call a user defined dispatch callback. The packet header carries
 * the dipatch id. 
 */
int recv_packet (void                       * param,
		 MUHWI_PacketHeader_t       * hdr,
		 uint32_t                     bytes) 
{
  MsgContext_t *context   = (MsgContext_t *) param;
  context->HWEndTime      = GetTimeBase();

  SoftwareBytes_t *sw_hdr = (SoftwareBytes_t *)(&hdr->messageUnitHeader.Packet_Types.Memory_FIFO);

  char *dst = context->buf + context->bytes_recvd;
  char *src = (char *)hdr + 32;

  context->size = sw_hdr->BytesStruct.message_size_in_bytes;
  bytes -= 32; /* get the payload bytes */

  //MU can pad message with zeros
  if (context->bytes_recvd + bytes >= context->size) {
    bytes = context->size - context->bytes_recvd;
    context->done = 1;
  }
  
  context->bytes_recvd += bytes;  

  if (bytes <= 16) 
    while (bytes --) 
      *(dst++) = *(src++);
  else
    memcpy(dst, src, bytes);
  
  return 0;
}

MUSPI_InjFifoSubGroup_t   ififo_subgroup[2];
MUSPI_RecFifoSubGroup_t   rfifo_subgroup[2];

//The source of the ping pong message
int   ROOT_A = 0;
int   ROOT_B = 0;
int   ROOT_C = 0;
int   ROOT_D = 0;
int   ROOT_E = 0;
int   ROOT_T = 0;

//The address of the neighbor that bounces it back
int   NEIGHBOR_A = 0;
int   NEIGHBOR_B = 0;
int   NEIGHBOR_C = 1;
int   NEIGHBOR_D = 0;
int   NEIGHBOR_E = 0;
int   NEIGHBOR_T = 0;

// allocates area for message send buffer
uint64_t sbuf[MAX_MESSAGE_SIZE/sizeof(uint64_t)];

// allocates area for message recv buffer
uint64_t rbuf[MAX_MESSAGE_SIZE/sizeof(uint64_t)];

#define test_exit  exit

int g_proc_id, g_nproc, g_cart_id, g_proc_coords[4], g_nb_list[8];
int g_nproc_t, g_nproc_x, g_nproc_y, g_nproc_z;
MPI_Comm g_cart_grid;
int g_nb_x_up, g_nb_x_dn;
int g_nb_y_up, g_nb_y_dn;
int g_nb_t_up, g_nb_t_dn;
int g_nb_z_up, g_nb_z_dn;



int main(int argc, char **argv)
{
  int rc;
  Personality_t pers;
  //works in cnk ?
  Kernel_GetPersonality(&pers, sizeof(pers));

  int g_proc_id;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
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


  int32_t     my_a = 0;
  int32_t     my_b = 0;
  int32_t     my_c = 0;
  int32_t     my_d = 0;
  int32_t     my_e = 0;
  int32_t     my_t = 0;  //Set from kernel sys call
  uint i = 0;
  uint64_t gi_timeout = 1600000000; // about 1 sec at 16 mhz
  gi_timeout *= 30;

  char rec_memory_fifo_buffer[REC_MEMORY_FIFO_SIZE+32];
  void *rec_memory_fifo = (void*)(rec_memory_fifo_buffer);
  rec_memory_fifo = (void*)( ((uint64_t)rec_memory_fifo_buffer + 32)  & ~(31UL) );
  char inj_memory_fifo_buffer[INJ_MEMORY_FIFO_SIZE+64];
  void *inj_memory_fifo = (void*)( ((uint64_t)inj_memory_fifo_buffer + 64)  & ~(63UL) );

  my_a = pers.Network_Config.Acoord;
  my_b = pers.Network_Config.Bcoord;
  my_c = pers.Network_Config.Ccoord;
  my_d = pers.Network_Config.Dcoord;
  my_e = pers.Network_Config.Ecoord;

  //Assume atomst one process per core
  my_t = Kernel_PhysicalProcessorID ();
  
  if (my_t > 1)
    my_t = 1;

  fprintf(stdout,"# MPI Process %d of %d on %s: cart_id %d, coordinates (%d %d %d %d)\n# MPI Process %d has personality %d %d %d %d %d %d\n",
	  g_proc_id, g_nproc, processor_name, g_cart_id, 
	  g_proc_coords[0], g_proc_coords[1], g_proc_coords[2], g_proc_coords[3],
	  g_proc_id, my_a, my_b, my_c, my_d, my_e, my_t);
  fflush(stdout);

  int mypers[6];
  mypers[0] = my_a; mypers[1] = my_b;  mypers[2] = my_c;  mypers[3] = my_d;  mypers[4] = my_e; mypers[5] = my_t;
  int nbtplus[6], nbtminus[6];
  int nbxplus[6], nbxminus[6];
  int nbyplus[6], nbyminus[6];
  int nbzplus[6], nbzminus[6];
  MPI_Status mstatus;
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_t_up, 0, 
	       (void*)nbtminus, 6, MPI_INT, g_nb_t_dn, 0,
	       g_cart_grid, &mstatus);
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_t_dn, 1, 
	       (void*)nbtplus, 6, MPI_INT, g_nb_t_up, 1, 
	       g_cart_grid, &mstatus);
  
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_x_up, 2, 
	       (void*)nbxminus, 6, MPI_INT, g_nb_x_dn, 2, 
	       g_cart_grid, &mstatus);
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_x_dn, 3, 
	       (void*)nbxplus, 6, MPI_INT, g_nb_x_up, 3, 
	       g_cart_grid, &mstatus);
  
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_y_up, 4, 
	       (void*)nbyminus, 6, MPI_INT, g_nb_y_dn, 4, 
	       g_cart_grid, &mstatus);
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_y_dn, 5, 
	       (void*)nbyplus, 6, MPI_INT, g_nb_y_up, 5, 
	       g_cart_grid, &mstatus);
  
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_z_up, 6, 
	       (void*)nbzminus, 6, MPI_INT, g_nb_z_dn, 6, 
	       g_cart_grid, &mstatus);
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_z_dn, 7, 
	       (void*)nbzplus, 6, MPI_INT, g_nb_z_up, 7, 
	       g_cart_grid, &mstatus);
  
  if(g_proc_id == 0) {
    printf("my coords (%d %d %d %d)\n", g_proc_coords[0], g_proc_coords[1], g_proc_coords[2], g_proc_coords[3]);
    printf("my kernel %d %d %d %d %d\n", my_a, my_b, my_c, my_d, my_e);
    printf("my tplus  %d %d %d %d %d\n", nbtplus[0], nbtplus[1], nbtplus[2], nbtplus[3], nbtplus[4]);
    printf("my tminus  %d %d %d %d %d\n", nbtminus[0], nbtminus[1], nbtminus[2], nbtminus[3], nbtminus[4]);
    printf("my xplus  %d %d %d %d %d\n", nbxplus[0], nbxplus[1], nbxplus[2], nbxplus[3], nbxplus[4]);
    printf("my xminus  %d %d %d %d %d\n", nbxminus[0], nbxminus[1], nbxminus[2], nbxminus[3], nbxminus[4]);
    printf("my yplus  %d %d %d %d %d\n", nbyplus[0], nbyplus[1], nbyplus[2], nbyplus[3], nbyplus[4]);
    printf("my yminus  %d %d %d %d %d\n", nbyminus[0], nbyminus[1], nbyminus[2], nbyminus[3], nbyminus[4]);
    printf("my zplus  %d %d %d %d %d\n", nbzplus[0], nbzplus[1], nbzplus[2], nbzplus[3], nbzplus[4]);
    printf("my zminus  %d %d %d %d %d\n", nbzminus[0], nbzminus[1], nbzminus[2], nbzminus[3], nbzminus[4]);
  }


  int _root = 0, _neighbor = 0; 

  if (my_a == ROOT_A &&
      my_b == ROOT_B &&
      my_c == ROOT_C &&
      my_d == ROOT_D &&
      my_e == ROOT_E &&
      my_t == ROOT_T)
    {
      printf("We are the root node with coords %d %d %d %d %d %d sending to coords %d %d %d %d %d %d\n",my_a,my_b,my_c,my_d,my_e,my_t,NEIGHBOR_A,NEIGHBOR_B,NEIGHBOR_C,NEIGHBOR_D,NEIGHBOR_E,NEIGHBOR_T);
      _root = 1;
    }
  else if (my_a == NEIGHBOR_A &&
	   my_b == NEIGHBOR_B &&
	   my_c == NEIGHBOR_C &&
	   my_d == NEIGHBOR_D &&
	   my_e == NEIGHBOR_E &&
	   my_t == NEIGHBOR_T)
    {
      printf("We are the neighbor node with coords %d %d %d %d %d %d\n",my_a,my_b,my_c,my_d,my_e,my_t);      
      _neighbor = 1;
    }

  // Initializes the send buffer
  for (i=0;i<MAX_MESSAGE_SIZE/8;i++) 
    sbuf [i] = (uint64_t)(i+4) + 0x0100000000000000ull;
  
  // clears the recv buffer
  for (i=0;i<MAX_MESSAGE_SIZE/8;i++)
  {
    rbuf[i] = 0x00;
  }
  
  // Initialize Injection FIFOs
  uint32_t fifoid = 0;  
  Kernel_InjFifoAttributes_t injFifoAttrs[1];
  injFifoAttrs[0].RemoteGet = 0;
  injFifoAttrs[0].System    = 0;  
  rc = Kernel_AllocateInjFifos (my_t,
				&ififo_subgroup[my_t], 
				1,
				&fifoid, 
				injFifoAttrs);
  
  /// Map virtual address
  Kernel_MemoryRegion_t  mregion;
  Kernel_CreateMemoryRegion (&mregion,inj_memory_fifo,INJ_MEMORY_FIFO_SIZE+1);
  Kernel_InjFifoInit (&ififo_subgroup[my_t], 
		      fifoid, 
		      &mregion, 
		      (uint64_t)inj_memory_fifo - (uint64_t)mregion.BaseVa, 
		      INJ_MEMORY_FIFO_SIZE);    
  Kernel_InjFifoActivate (&ififo_subgroup[my_t], 1, &fifoid, KERNEL_INJ_FIFO_ACTIVATE);    
  
  // Initialize Reception FIFOs
  uint32_t rfifoid = 0;
  memset (&my_recv_context[my_t], 0, sizeof(my_recv_context[0]));  
  MUSPI_RegisterRecvFunction (recv_packet, &my_recv_context[my_t], my_t+1);
  my_recv_context[my_t].buf = (char *) rbuf;  
  Kernel_CreateMemoryRegion (&mregion, rec_memory_fifo, REC_MEMORY_FIFO_SIZE);
  Kernel_RecFifoAttributes_t recFifoAttrs[1];
  recFifoAttrs[0].System = 0;

  Kernel_AllocateRecFifos (my_t, 
			   &rfifo_subgroup[my_t], 
			   1,
			   &rfifoid,
			   recFifoAttrs);
  
  Kernel_RecFifoInit    (& rfifo_subgroup[my_t], rfifoid, 
			 &mregion, 
			 (uint64_t)rec_memory_fifo - (uint64_t)mregion.BaseVa, 
			 REC_MEMORY_FIFO_SIZE);

  uint64_t recFifoEnableBits=0;
  
  recFifoEnableBits |= ( 0x0000000000000001ULL << 
			 ( 15 - ( (my_t/*sgid*/*BGQ_MU_NUM_REC_FIFOS_PER_SUBGROUP) + 0/*RecFifoId*/ )) );

  Kernel_RecFifoEnable ( 0, /* Group ID */ 
			 recFifoEnableBits );
  
  // Initialize GI Barrier
  //  MUSPI_GIBarrier_t GIBarrier;  
  // Initialize the barrier, resetting the hardware.
  //rc = MUSPI_GIBarrierInit ( &GIBarrier, 0 /* classRouteId */ );
  //if (rc)
  //  {
  //    printf("MUSPI_GIBarrierInit for class route %u returned rc = %d\n",0, rc);
  //    test_exit(__LINE__);
  //  }  
  //if (rc) test_exit(__LINE__);
  
  // Enter the MU barrier
  //rc = MUSPI_GIBarrierEnter ( &GIBarrier );
  //if (rc)
  //  {
  //    printf("MUSPI_GIBarrierEnter failed on iteration = %d, returned rc = %d\n", i, rc);
  //    test_exit(1);
  //  }
  
  // Poll for completion of the barrier.
  //rc = MUSPI_GIBarrierPollWithTimeout ( &GIBarrier, gi_timeout);
  //if ( rc )
  //  {
  //    printf("MUSPI_GIBarrierPollWithTimeout failed on iteration = %d, returned rc = %d\n", i, rc);
  //    test_exit(1);
  //  }
  
  // Test ping pong
  if (_root)
    ping (NEIGHBOR_A, NEIGHBOR_B, NEIGHBOR_C, NEIGHBOR_D, NEIGHBOR_E, NEIGHBOR_T, MAX_MESSAGE_SIZE, my_t); 
  else if (_neighbor) {
    pong (ROOT_A, ROOT_B, ROOT_C, ROOT_D, ROOT_E, ROOT_T, MAX_MESSAGE_SIZE, my_t);
    
  }
  
  Delay(1000); // Make sure all processes are done
  MPI_Finalize();
  return 0;
}


int ping   ( int      a,
	     int      b,
	     int      c, 
	     int      d, 
	     int      e, 
	     int      t,
	     int      bytes,
	     int      my_t)
{
  int rc =0;
  int i=0;
  unsigned long long StartTime=0,EndTime,MeanTime;
  unsigned long long HWStartTime=0, HWTotalTime=0, HWMeanTime;
  
  rc =0;
 
  SoftwareBytes_t  SoftwareBytes;
  memset( &SoftwareBytes, 0x00, sizeof(SoftwareBytes_t) );

  Delay(500000); // Make sure receiver is ready
  
  // bit 0 ===> A- Torus FIFO
  uint64_t torusInjectionFifoMap = MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AP; 
  
  // Set up the memory fifo information structure for message
  // 1.  It is used as input to build a memory fifo descriptor.
  memset( &mu_iMemoryFifoDescriptorInfo[my_t], 0x00, 
	  sizeof(mu_iMemoryFifoDescriptorInfo[0]) );
  
  MUSPI_SetUpDestination ( &mu_iMemoryFifoDescriptorInfo[my_t].Base.Dest,
			   a, b, c, d, e );
  SoftwareBytes.BytesStruct.functionid = t+1;
  SoftwareBytes.BytesStruct.message_size_in_bytes = bytes;
  
  msg_BuildPt2PtMemoryFifoInfo ( &mu_iMemoryFifoDescriptorInfo[my_t],
				 SoftwareBytes,
				 t*4,   //use the t to isolate the recfifo
				 0, 
				 (uint64_t)sbuf,
				 bytes,
				 torusInjectionFifoMap,
				 MUHWI_PACKET_VIRTUAL_CHANNEL_DETERMINISTIC,
				 0,
				 MUHWI_PACKET_HINT_E_NONE);
  
  rc = MUSPI_CreatePt2PtMemoryFIFODescriptor 
    ( &(mu_iMemoryFifoDescriptor[my_t]),
      &(mu_iMemoryFifoDescriptorInfo[my_t]) );
  
  assert(rc ==0);
      
  // -----------------------------------------
  //    Inject this message into the fifo
  //    - For a hardware timing, start the timer after the injection.
  //      The end time is recorded in the context when the packet is received.
  // ----------------------------------------
  rc = MUSPI_InjFifoInject (MUSPI_IdToInjFifo(0, &ififo_subgroup[my_t]), 
			    &mu_iMemoryFifoDescriptor[my_t]);
  
  HWStartTime =  GetTimeBase();
  
  int hops = a + b + c + d + e;
  printf("ping:  MessageSize=%d, Iterations=%d, a=%3d, b=%3d, c=%3d, d=%3d, e=%3d, hops=%3d\n",
         bytes, 1, a, b, c, d, e, hops);
  
  return rc;
}


int pong   ( int      a,
	     int      b,
	     int      c, 
	     int      d, 
	     int      e, 
	     int      t,
	     int      bytes,
	     int      my_t)
{
  int rc =0;
  int i=0;
  unsigned long long StartTime,EndTime,MeanTime;

  SoftwareBytes_t  SoftwareBytes;
  memset( &SoftwareBytes, 0x00, sizeof(SoftwareBytes_t) );
  
  rc =0;  
  int rfifoid = 0;

  // ---------------------------------------------
  //    Reception Side
  // ---------------------------------------------
  while (!my_recv_context[my_t].done) 
    MUSPI_RecFifoPoll (MUSPI_IdToRecFifo (rfifoid, &rfifo_subgroup[my_t]), 
		       1000);
  my_recv_context[my_t].done = 0;
  my_recv_context[my_t].bytes_recvd = 0;    

    
  int hops = a+b+c+d+e;
  printf("pong:  MessageSize=%d, Iterations=%d, a=%3d, b=%3d, c=%3d, d=%3d, e=%3d, hops=%3d\n",
         bytes, 1, a, b, c, d, e,hops);
  //  assert(packets_received <= NUM_LOOPS+1);
  return rc;
}

