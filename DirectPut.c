#include "spi.h"
#include <mpi.h>

#define SEND_BUFFER_ALIGNMENT   128
#define RECV_BUFFER_ALIGNMENT   128
#define MAX_MESSAGE_SIZE       32768

// Sub message size
static int window_size  =  2048;   //  size of submessages
long long messageSizeInBytes = MAX_MESSAGE_SIZE;


// we have four directions and forward/backward
#define NUM_DIRS               8
#define NUM_INJ_FIFOS          NUM_DIRS
#define INJ_MEMORY_FIFO_SIZE  ((64*1024) -1)
#define N_LOOPS                10000

// Allocate static memory for descriptors
char muDescriptorsMemory[ NUM_INJ_FIFOS * sizeof(MUHWI_Descriptor_t) + 64 ];
// pointer to descriptor array
MUHWI_Descriptor_t *muDescriptors;

const int batsubgroupID = 0;
int do_dynamic      = 1;
// Enable different zone routing modes
uint8_t  zoneRoutingMask = 0;
unsigned zoneRoutingId   = 0;
// stay on bubble bits
uint8_t stayOnBubbleMask  = 0;
unsigned stayOnBubbleFlag = 0;

// Allocate static memory for send and receive buffers
char sendBufMemory[NUM_DIRS * MAX_MESSAGE_SIZE+ SEND_BUFFER_ALIGNMENT];
char recvBufMemory[NUM_DIRS * MAX_MESSAGE_SIZE+ SEND_BUFFER_ALIGNMENT];
// pointers to send and receive buffers
char * recvBuffers;
char * sendBuffers;

// neighbour destination cache
struct { 
  MUHWI_Destination_t dest;
  uint8_t             hintsABCD;
  uint8_t             hintsE;
} nb2dest[NUM_DIRS];

// receive counter
volatile uint64_t recvCounter;

// base addess table slot for receive buffer and counter
uint32_t recvBufBatId = 0, recvCntrBatId = 1;

// physical address of send buffers
uint64_t sendBufPAddr;

// get the destinations for all neighbours
// will be saved in nb2dest
int get_destinations(int * mypers);

// Call to create the descriptors for all eight directions
void create_descriptors();

// Call to set up the base address table id and memory regions
void setup_mregions_bats_counters();

MUSPI_GIBarrier_t GIBarrier;

void alltoall_exit(const int rc) {
  exit (rc);
  return;
}

void global_barrier() {
  int rc = 0;
  uint64_t timeoutCycles = 60UL * 1600000000UL; // about 60 sec at 1.6 ghz
  rc = MUSPI_GIBarrierEnter ( &GIBarrier );
  if (rc) {
    printf("MUSPI_GIBarrierEnter failed returned rc = %d\n", rc);
    alltoall_exit(1);
  }
  
  // Poll for completion of the barrier.
  rc = MUSPI_GIBarrierPollWithTimeout ( &GIBarrier, timeoutCycles);
  if( rc ) {
    printf("MUSPI_GIBarrierPollWithTimeout failed returned rc = %d\n", rc);
    DelayTimeBase (200000000000UL);
    alltoall_exit(1);
  }
  return;
}

// here come the MPI variables
int g_proc_id, g_nproc, g_cart_id, g_proc_coords[4], g_nb_list[8];
int g_nproc_t, g_nproc_x, g_nproc_y, g_nproc_z;
MPI_Comm g_cart_grid;
int g_nb_x_up, g_nb_x_dn;
int g_nb_y_up, g_nb_y_dn;
int g_nb_t_up, g_nb_t_dn;
int g_nb_z_up, g_nb_z_dn;
int g_mpi_prov=0;



int main(int argc, char **argv) {
  int rc;
  Personality_t pers;
  //works in cnk ?
  Kernel_GetPersonality(&pers, sizeof(pers));

  int g_proc_id;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  //MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &g_mpi_prov);
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);
  if(g_proc_id == 0){
    printf("provided thread support = %d\n", g_mpi_prov);
  }
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

  int mypers[6];
  mypers[0] = pers.Network_Config.Acoord;
  mypers[1] = pers.Network_Config.Bcoord;
  mypers[2] = pers.Network_Config.Ccoord;
  mypers[3] = pers.Network_Config.Dcoord;
  mypers[4] = pers.Network_Config.Ecoord;

  get_destinations(mypers);

  // Setup the FIFO handles
  msg_InjFifoHandle_t injFifoHandle;
  rc = msg_InjFifoInit ( &injFifoHandle,
			 0,        /* startingSubgroupId */
			 0,        /* startingFifoId     */
			 NUM_INJ_FIFOS,       /* numFifos   */
			 INJ_MEMORY_FIFO_SIZE+1, /* fifoSize */
			 NULL      /* Use default attributes */
			 );
  if(rc != 0) {
    fprintf(stderr, "msg_InjFifoInit failed with rc=%d\n",rc);
    alltoall_exit(1);
  }

  recvBuffers = (char *)(((uint64_t)recvBufMemory+RECV_BUFFER_ALIGNMENT)&~(RECV_BUFFER_ALIGNMENT-1));    
  sendBuffers = (char *)(((uint64_t)sendBufMemory+SEND_BUFFER_ALIGNMENT)&~(SEND_BUFFER_ALIGNMENT-1));

  // Set up base address table for reception counter and buffer
  setup_mregions_bats_counters();

  // Create descriptors
  // Injection Direct Put Descriptor, one for each neighbor
  muDescriptors =
    ( MUHWI_Descriptor_t *)(((uint64_t)muDescriptorsMemory+64)&~(64-1));
  create_descriptors();

  uint64_t totalCycles=0;
  uint64_t startTime=0;
  startTime = GetTimeBase();

  // Initialize the barrier, resetting the hardware.
  rc = MUSPI_GIBarrierInit ( &GIBarrier, 0 /*comm world class route */);

  if(rc) {
    printf("MUSPI_GIBarrierInit returned rc = %d\n", rc);
    alltoall_exit(__LINE__);
  }

  // Fill send buffer
  for(int n = 0; n < NUM_DIRS * MAX_MESSAGE_SIZE/sizeof(double); n+=8) {
    *(double*)&sendBufMemory[n] = (double) g_cart_id;
  }


  uint64_t descCount[NUM_INJ_FIFOS];
  double s = 0;
  for(int l = 0; l < N_LOOPS; l++) {
    // reset the recv counter 
    recvCounter = NUM_DIRS*messageSizeInBytes;
    global_barrier(); // make sure everybody is set recv counter
    
    for (uint64_t bytes = 0; bytes < messageSizeInBytes; bytes += window_size) {
      uint64_t msize = (bytes <= messageSizeInBytes - window_size) ? window_size : (messageSizeInBytes - bytes);
      for (int j = 0; j < NUM_DIRS; j++) {
	muDescriptors[j].Message_Length = msize; 
	muDescriptors[j].Pa_Payload    =  sendBufPAddr + (messageSizeInBytes * j) + bytes;
	MUSPI_SetRecPutOffset (&muDescriptors[j], (messageSizeInBytes * j) +  bytes);
	
	descCount[ j % NUM_INJ_FIFOS] =
	  msg_InjFifoInject ( injFifoHandle,
			      j % NUM_INJ_FIFOS,
			      &muDescriptors[j]);
      }
    }

    // do some computation to hide communication
    //for(int m = 0; m < 4; m++) {
    //  for(int n = 0; n < NUM_DIRS * MAX_MESSAGE_SIZE/sizeof(double); n+=8) {
    //	s += *(double*)&sendBufMemory[n];
    //  }
    //}
    
    // wait for receive completion
    while ( recvCounter > 0 );


    // wait for send completion
    unsigned sendDone;
    do {
      sendDone = 0;
      for (int j = 0; j < NUM_INJ_FIFOS; j++ )
	sendDone += msg_InjFifoCheckCompletion( injFifoHandle,
						j,
						descCount[j]);
    }
    while ( sendDone < NUM_INJ_FIFOS );
    if(g_proc_id == -1) printf("Send and receive complete... %d %d\n", g_cart_id, l);
    _bgq_msync(); // Ensure data is available to all cores.  

    // do some computation not hidding communication
    //for(int m = 0; m < 4; m++) {
    //  for(int n = 0; n < NUM_DIRS * MAX_MESSAGE_SIZE/sizeof(double); n+=8) {
    //	s += *(double*)&sendBufMemory[n];
    //  }
    //}
  }

  totalCycles = GetTimeBase() - startTime;
  if(g_proc_id == 0) {
    printf("total cycles per loop= %llu\n", (long long unsigned int) totalCycles/N_LOOPS);
  }
  printf("res for %d is %e\n",g_proc_id, s);

  msg_InjFifoTerm ( injFifoHandle );

  MPI_Finalize();
  return(0);
}



void setup_mregions_bats_counters() {
  const uint64_t buffersSize =  NUM_DIRS * messageSizeInBytes;

  // allocate bat entries for the recive buffer and the receive counter
  
  uint32_t batIds[2] = { recvBufBatId, recvCntrBatId };
  MUSPI_BaseAddressTableSubGroup_t batSubGrp;
  
  int rc =  Kernel_AllocateBaseAddressTable( batsubgroupID/*subgrpId*/,
					     &batSubGrp,
					     2,/*nbatids*/
					     batIds,
					     0 /* "User" use */);
  
  if (rc != 0) {
    fprintf(stderr, "Kernel_AllocateBaseAddressTable failed with rc=%d\n", rc);
    alltoall_exit(1);
  }
  
  // Receive buffer bat is set to the PA addr of the receive buffer
  Kernel_MemoryRegion_t memRegion;
  rc = Kernel_CreateMemoryRegion ( &memRegion,
				   recvBuffers,
				   buffersSize);
  if ( rc != 0) {
    printf("Kernel_CreateMemoryRegion failed with rc=%d\n",rc);
    alltoall_exit(1);
  }
  
  uint64_t paAddr = 
    (uint64_t)recvBuffers - 
    (uint64_t)memRegion.BaseVa + 
    (uint64_t)memRegion.BasePa;
  
  rc = MUSPI_SetBaseAddress ( &batSubGrp,
			      recvBufBatId,
			      paAddr );
  
  if(rc != 0) {
    printf("MUSPI_SetBaseAddress failed with rc=%d\n",rc);
    alltoall_exit(1);
  }
  
  // Receive counter bat is set to the MU style atomic PA addr of the receive counter
  if( (uint64_t)(&recvCounter) & 0x7 ) {
    printf("ERROR: recv counter is not 8 byte aligned\n");
    alltoall_exit(1);
  }
  
  rc = Kernel_CreateMemoryRegion ( &memRegion,
				   (void *)&recvCounter,
				   sizeof(recvCounter));
  if(rc != 0) {
    printf("Kernel_CreateMemoryRegion failed with rc=%d\n",rc);
    alltoall_exit(1);
  }
  
  paAddr = 
    (uint64_t)&recvCounter - 
    (uint64_t)memRegion.BaseVa + 
    (uint64_t)memRegion.BasePa;
  
  uint64_t paAddrAtomic =  MUSPI_GetAtomicAddress(paAddr,MUHWI_ATOMIC_OPCODE_STORE_ADD);
  
  rc = MUSPI_SetBaseAddress ( &batSubGrp,
			      recvCntrBatId,
			      paAddrAtomic );
  
  if(rc != 0) {
    printf("MUSPI_SetBaseAddress failed with rc=%d\n",rc);
    alltoall_exit(1);
  }
  
  // Get the send buffers physical address
  rc = Kernel_CreateMemoryRegion ( &memRegion,
				   sendBuffers,
				   buffersSize);
  if(rc != 0) {
    printf("Kernel_CreateMemoryRegion failed with rc=%d\n",rc);
    alltoall_exit(1);
  }
  
  sendBufPAddr = 
    (uint64_t)sendBuffers - 
    (uint64_t)memRegion.BaseVa + 
    (uint64_t)memRegion.BasePa;
  return;
}


void create_descriptors() {
  uint64_t anyFifoMap = 
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AM |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AP | 
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BM |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BP | 
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CM |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CP |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DM |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DP |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EM |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EP;
  
  uint64_t offset;
  static int  did_print =0;
 
  // loop over directions
  // CHECK offset needs to be adjusted for QCD case
  for(int i = 0, offset = 0; i < 8; i++, offset += messageSizeInBytes) {
    // Injection Direct Put Descriptor Information Structure
    MUSPI_Pt2PtDirectPutDescriptorInfo_t dinfo;
    
    memset( (void*)&dinfo, 0x00, sizeof(dinfo) );
      
    dinfo.Base.Payload_Address = sendBufPAddr + offset;
    dinfo.Base.Message_Length  = messageSizeInBytes;
    dinfo.Base.Torus_FIFO_Map  = anyFifoMap;
      
    dinfo.Base.Dest = nb2dest[i].dest;
      
    dinfo.Pt2Pt.Hints_ABCD = nb2dest[i].hintsABCD; 

    if(do_dynamic) {	  
      dinfo.Pt2Pt.Misc1 =
	nb2dest[i].hintsE |
	MUHWI_PACKET_USE_DYNAMIC_ROUTING |  
	MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE;
      
      dinfo.Pt2Pt.Misc2 = 
	MUHWI_PACKET_VIRTUAL_CHANNEL_DYNAMIC | 
	zoneRoutingMask | 
	stayOnBubbleMask;
      if ( (g_cart_id ==0) && (did_print ==0)) 
	printf(" dynamic routing  zoneRoutingMask=%d stayOnBubbleMask=%d\n",
	       zoneRoutingMask, stayOnBubbleMask);
    }
    else {	    	    
      dinfo.Pt2Pt.Misc1 =
	nb2dest[i].hintsE |
	MUHWI_PACKET_USE_DETERMINISTIC_ROUTING |  
	MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE;
	
      dinfo.Pt2Pt.Misc2 = 
	MUHWI_PACKET_VIRTUAL_CHANNEL_DETERMINISTIC | 
	zoneRoutingMask | 
	stayOnBubbleMask;
      if ( (g_cart_id ==0) && (did_print ==0)) printf(" deterministic routing\n");
    }
    did_print++;
    
    dinfo.Pt2Pt.Skip  = 8; // for checksumming, skip the header 	      
    dinfo.DirectPut.Rec_Payload_Base_Address_Id = recvBufBatId;
    dinfo.DirectPut.Rec_Payload_Offset          = offset;
    dinfo.DirectPut.Rec_Counter_Base_Address_Id = recvCntrBatId;
    dinfo.DirectPut.Rec_Counter_Offset          = 0;
      
    dinfo.DirectPut.Pacing = MUHWI_PACKET_DIRECT_PUT_IS_NOT_PACED;
      
    int rc = MUSPI_CreatePt2PtDirectPutDescriptor(&muDescriptors[i],
						  &dinfo );
    if (rc != 0) {
      fprintf(stderr, "MUSPI_CreatePt2PtDirectPutDescriptor failed with rc=%d\n",rc);
      alltoall_exit(1);
    }
      
  }
}


int get_destinations(int * mypers) {

  int tmp[6];
  MPI_Status mstatus;
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_t_up, 0, 
	       (void*)tmp, 6, MPI_INT, g_nb_t_dn, 0,
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[1].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_t_dn, 1, 
	       (void*)tmp, 6, MPI_INT, g_nb_t_up, 1, 
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[0].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );
  
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_x_up, 2, 
	       (void*)tmp, 6, MPI_INT, g_nb_x_dn, 2, 
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[3].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_x_dn, 3, 
	       (void*)tmp, 6, MPI_INT, g_nb_x_up, 3, 
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[2].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );  

  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_y_up, 4, 
	       (void*)tmp, 6, MPI_INT, g_nb_y_dn, 4, 
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[5].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );  
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_y_dn, 5, 
	       (void*)tmp, 6, MPI_INT, g_nb_y_up, 5, 
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[4].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );  
  
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_z_up, 6, 
	       (void*)tmp, 6, MPI_INT, g_nb_z_dn, 6, 
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[7].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );  
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_z_dn, 7, 
	       (void*)tmp, 6, MPI_INT, g_nb_z_up, 7, 
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[6].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );  
  
  return(0);
}
