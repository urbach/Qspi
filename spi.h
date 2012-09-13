/* begin_generated_IBM_copyright_prolog                             */
/*                                                                  */
/* This is an automatically generated copyright prolog.             */
/* After initializing,  DO NOT MODIFY OR MOVE                       */
/* ================================================================ */
/*                                                                  */
/* Licensed Materials - Property of IBM                             */
/*                                                                  */
/* Blue Gene/Q                                                      */
/*                                                                  */
/* (C) Copyright IBM Corp.  2008, 2012                              */
/*                                                                  */
/* US Government Users Restricted Rights -                          */
/* Use, duplication or disclosure restricted                        */
/* by GSA ADP Schedule Contract with IBM Corp.                      */
/*                                                                  */
/* This software is available to you under the                      */
/* Eclipse Public License (EPL).                                    */
/*                                                                  */
/* ================================================================ */
/*                                                                  */
/* end_generated_IBM_copyright_prolog                               */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <stdint.h>
#include <malloc.h>
#include <assert.h>
#include <string.h>

// Basic SPI and HWI includes
#include <hwi/include/bqc/A2_core.h>
#include <hwi/include/bqc/A2_inlines.h>
#include <hwi/include/bqc/MU_PacketCommon.h>
#include <firmware/include/personality.h>
#include <spi/include/mu/Descriptor.h>
#include <spi/include/mu/Descriptor_inlines.h>
#include <spi/include/mu/InjFifo.h>
#include <spi/include/mu/Addressing.h>
#include <spi/include/mu/Addressing_inlines.h>
#include <spi/include/mu/GIBarrier.h>
#include <spi/include/kernel/MU.h>
#include <spi/include/kernel/process.h>
#include <spi/include/kernel/location.h>

/**
 * \brief Injection Fifo Handle
 *
 * This is a "handle" returned from msg_InjFifoInit() and passed into subsequent
 * calls to msg_InjFifoXXXX() functions.  It is used internally within the
 * msg_InjFifoXXXX() functions to anchor resources that have been allocated.
 */
typedef struct {
  void* pOpaqueObject;
} msg_InjFifoHandle_t;


int msg_InjFifoInit ( msg_InjFifoHandle_t *injFifoHandlePtr,
                      uint32_t             startingSubgroupId,
                      uint32_t             startingFifoId,
                      uint32_t             numFifos,
                      size_t               fifoSize,
                      Kernel_InjFifoAttributes_t  *injFifoAttrs );


/**
 * \brief Terminate Injection Fifos
 * 
 * Terminate the usage of injection fifos.  This deactivates the fifos and
 * frees all of the storage associated with them (previously allocated during
 * msg_InjFifoInit()).
 * 
 * \param [in]  injFifoHandle  The handle returned from msg_InjFifoInit().
 *                             It must be passed into this function untouched
 *                             from when it was returned from msg_InjFifoInit().
 *
 * \note After this function returns, no more InjFifo functions should be called
 *       with this injFifoHandle.
 */
void msg_InjFifoTerm ( msg_InjFifoHandle_t injFifoHandle );


/**
 * \brief Inject Descriptor into Injection Fifo
 * 
 * Inject the specified descriptor into the specified injection fifo.
 * 
 * \param [in]  injFifoHandle  The handle returned from msg_InjFifoInit().
 *                             It must be passed into this function untouched
 *                             from when it was returned from msg_InjFifoInit().
 * \param [in]  relativeFifoId  The fifo number, relative to the start of
 *                              the fifos managed by this opaque object.
 *                              For example, if msg_InjFifoInit() was called
 *                              to init fifos in subgroup 2, starting with
 *                              fifo Id 3, the relativeFifoNumber of the
 *                              first fifo is 0, not 3.
 * \param [in]  descPtr         Pointer to the descriptor to be injected.
 *
 * \retval  positiveNumber  The descriptor was successfully injected.  The
 *                          returned value is the sequence number of this 
 *                          descriptor.
 * \retval  -1              The descriptor was not injected, most likely because
 *                          there is no room in the fifo.
 */
uint64_t msg_InjFifoInject ( msg_InjFifoHandle_t injFifoHandle,
                             uint32_t            relativeFifoId,
                             MUHWI_Descriptor_t *descPtr );


