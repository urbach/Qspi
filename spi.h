/* Init software bytes for the packet header. */
typedef union _swb
{
  struct {
    uint8_t reserved[6];
    uint8_t bytes[14];
  };
  
  struct {
    uint8_t  reserved[6];
    uint8_t  functionid;
    uint8_t  bytes[5];
    
    uint32_t message_size_in_bytes;
    uint8_t unused[4];
  } BytesStruct;
  
} SoftwareBytes_t;


static inline void msg_BuildPt2PtMemoryFifoInfo 
(MUSPI_Pt2PtMemoryFIFODescriptorInfo_t  * minfo,
 SoftwareBytes_t                          SoftwareBytes,  
 uint                                     rfifoid,
 uint64_t                                 put_offset,
 uint64_t                                 buffer,
 uint64_t                                 bytes,
 uint64_t                                 torusInjectionFifoMap,
 int                                      vc,
 uint8_t                                  hintsABCD,
 uint8_t                                  hintsE)
{
  minfo->Pt2Pt.Hints_ABCD = hintsABCD;
  if ( vc ==  MUHWI_PACKET_VIRTUAL_CHANNEL_DETERMINISTIC ) {
    minfo->Pt2Pt.Misc1      = hintsE         |
      MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE   |
      MUHWI_PACKET_USE_DETERMINISTIC_ROUTING |
      MUHWI_PACKET_DO_NOT_DEPOSIT;
    minfo->Pt2Pt.Misc2 = MUHWI_PACKET_VIRTUAL_CHANNEL_DETERMINISTIC;
  }
  else if (vc == MUHWI_PACKET_VIRTUAL_CHANNEL_DYNAMIC) {
    minfo->Pt2Pt.Misc1      = hintsE       |
      MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE |
      MUHWI_PACKET_USE_DYNAMIC_ROUTING     |
      MUHWI_PACKET_DO_NOT_DEPOSIT;
    minfo->Pt2Pt.Misc2      = MUHWI_PACKET_VIRTUAL_CHANNEL_DYNAMIC;
  }
  else if ( vc ==  MUHWI_PACKET_VIRTUAL_CHANNEL_SYSTEM ) {
    minfo->Pt2Pt.Misc1      = hintsE         |
      MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE   |
      MUHWI_PACKET_USE_DETERMINISTIC_ROUTING |
      MUHWI_PACKET_DO_NOT_DEPOSIT;
    minfo->Pt2Pt.Misc2 = MUHWI_PACKET_VIRTUAL_CHANNEL_SYSTEM;
  }
  else { MUSPI_assert(0); }
    
  minfo->Pt2Pt.Skip  = 0;
  minfo->MemFIFO.Rec_FIFO_Id    = rfifoid;
  minfo->MemFIFO.Rec_Put_Offset = put_offset;
  minfo->Base.Pre_Fetch_Only= MUHWI_DESCRIPTOR_PRE_FETCH_ONLY_NO;
  minfo->MemFIFO.Interrupt  = MUHWI_PACKET_DO_NOT_INTERRUPT_ON_PACKET_ARRIVAL;
  minfo->MemFIFO.SoftwareBit = 0;
  memcpy( minfo->MemFIFO.SoftwareBytes,
          SoftwareBytes.bytes,
          sizeof( minfo->MemFIFO.SoftwareBytes ) );

  //Warning assume V == P
  minfo->Base.Payload_Address = buffer;
  minfo->Base.Message_Length =  bytes;
  minfo->Base.Torus_FIFO_Map = torusInjectionFifoMap;
}
