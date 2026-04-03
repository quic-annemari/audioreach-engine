/* ========================================================================
  @file capi_nxp_device.c
  @brief This file contains CAPI implementation of nxp device Module

  Copyright (c) Qualcomm Innovation Center, Inc. All Rights Reserved.
  SPDX-License-Identifier: BSD-3-Clause
==============================================================================*/

/*=====================================================================
  Includes
 ======================================================================*/
#include "capi_nxp_device_i.h"
#include <zephyr/drivers/dai.h>
#include <zephyr/drivers/dma.h>
#include <zephyr/sys/sys_heap.h>

#define PLATFORM_DCACHE_ALIGN 32

extern char _end[], _heap_sentry[];
#define heapmem ((uint8_t *)ALIGN_UP((uintptr_t)_end, PLATFORM_DCACHE_ALIGN))

//Heap for dma_buffer
static struct k_heap kernel_heap;

//const struct device *dai_devs[] = { DT_FOREACH_STATUS_OKAY(DT_NODELABEL(nxp_dai_sai), DEVICE_DT_GET(DT_NODELABEL(nxp_dai_sai))) };
const struct device *dma_device = DEVICE_DT_GET(DT_NODELABEL(sdma3));
const struct device *dai_dev = DEVICE_DT_GET(DT_NODELABEL(sai3));

/*------------------------------------------------------------------------
  Function name: capi_nxp_device_common_init
  DESCRIPTION: Initialize the CAPIv2 nxp_device module and library.
  This function can allocate memory.
  -----------------------------------------------------------------------*/
capi_err_t capi_nxp_device_common_init(capi_t *_pif, capi_proplist_t *init_set_properties, uint32_t dir)
{
   capi_err_t capi_result = CAPI_EOK;
   //AR_MSG(DBG_ERROR_PRIO,
   //         "CAPI_NXP_DEVICE: Init enter");

   if (NULL == _pif || NULL == init_set_properties)
   {
      AR_MSG(DBG_ERROR_PRIO,
            "CAPI_NXP_DEVICE: Init received bad pointer, 0x%lx",
            (uint32_t)_pif);
      return CAPI_EBADPARAM;
   }

   capi_nxp_device_t *me_ptr = (capi_nxp_device_t *)_pif;
   memset((void *)me_ptr, 0, sizeof(capi_nxp_device_t));

   // Allocate vtbl
   me_ptr->vtbl.vtbl_ptr = capi_nxp_device_get_vtbl();

   // Cache direction
   me_ptr->direction = 0; //sink
   //me_ptr->i2s_media_fmt.format.data_interleaving = CAPI_DEINTERLEAVED_UNPACKED;

   capi_result = capi_nxp_device_process_set_properties(me_ptr, init_set_properties);
   capi_result ^= (capi_result & CAPI_EUNSUPPORTED); // ignore unsupported
   if (CAPI_EOK != capi_result)
   {
      AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: init set properties failed");
      return capi_result;
   }

   //me_ptr->drift_info.heap_id = (POSAL_HEAP_ID)me_ptr->heap_mem.heap_id;

   capi_result = capi_nxp_device_dai_init(me_ptr, dir);
   if (capi_result < 0) {
      AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Failed to initialize dai");
      return CAPI_EFAILED;
   }

   capi_result = capi_nxp_device_init_param(me_ptr);
   if (capi_result < 0) {
      AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Failed to initialize device parameters");
      return CAPI_EFAILED;
   }

   return capi_result;
}

/*------------------------------------------------------------------------
  Function name: capi_nxp_device_dai_init
  DESCRIPTION: Initialize the CAPIv2 nxp_device dai.
  -----------------------------------------------------------------------*/
capi_err_t capi_nxp_device_dai_init(capi_nxp_device_t *me_ptr, uint32_t dir)
{
   capi_err_t capi_result = CAPI_EOK;
   int ret;
   //struct capi_dai_data dai_data;
   struct dai_config cfg;
   struct sai_params *sai_cfg_params;
   const void* cfg_params;

   if (NULL == me_ptr)
   {
      AR_MSG(DBG_ERROR_PRIO,
            "CAPI_NXP_DEVICE: Init received null property");
      return CAPI_EBADPARAM;
   }

   ret = dai_config_get(dai_dev, &cfg, 0);
   if (ret != 0) {
      AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Failed to get dai configuration");
      capi_result = CAPI_EFAILED;
      goto end;
   }

   if (dai_probe(dai_dev)) {
      AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Failed to probe dai");
      capi_result = CAPI_EFAILED;
      goto end;
   }

   sai_cfg_params = posal_memory_malloc(sizeof(struct sai_params), (POSAL_HEAP_ID)me_ptr->heap_mem.heap_id);
   if (NULL == sai_cfg_params) {
      AR_MSG(DBG_ERROR_PRIO,
            "CAPI_NXP_DEVICE: Failed to allocate sai configuration pointer");
      capi_result = CAPI_ENOMEMORY;
      goto end;
   }

   me_ptr->dai_dir = dir;

   //Set the configuration, hardcoded from the logs
   cfg.dai_index = 3;
   cfg.format = 1;
   cfg.options = 1;
   cfg.rate = 0;
   cfg.type = DAI_IMX_SAI;
   //Fill cfg_params
   sai_cfg_params->mclk_id = 0;
   sai_cfg_params->mclk_direction = 0;
   sai_cfg_params->mclk_rate = 12288000;
   sai_cfg_params->fsync_rate = 48000;
   sai_cfg_params->bclk_rate = 3072000;
   sai_cfg_params->tdm_slots = 2;
   sai_cfg_params->rx_slots = 3;
   sai_cfg_params->tx_slots = 3;
   sai_cfg_params->tdm_slot_width = 32;
   cfg_params = sai_cfg_params;
   ret = dai_config_set(dai_dev, &cfg, cfg_params);
   if (ret) {
      return CAPI_EFAILED;
   }

end:
   if (sai_cfg_params)
      posal_memory_free(sai_cfg_params);

   return capi_result;
}

/*------------------------------------------------------------------------
  Function name: capi_nxp_device_init_param
  DESCRIPTION: Sets either a parameter value or a parameter structure containing
  multiple parameters. In the event of a failure, the appropriate error code is
  returned.
  -----------------------------------------------------------------------*/
capi_err_t capi_nxp_device_init_param(capi_nxp_device_t *me_ptr)
{
   capi_err_t capi_result = CAPI_EOK;
   uint32_t period_bytes = 0;
   uint32_t period_count = 0;

   capi_result = capi_nxp_device_init_dma_buffer(me_ptr, &period_bytes, &period_count);
   if (capi_result < 0) {
      AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Failed to initialize dma buffer");
      return CAPI_EFAILED;
   }

   capi_result = capi_nxp_device_set_dma_config(me_ptr, period_bytes, period_count); 
   if (capi_result < 0) {
      AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Failed to set dma config");
      //sleep(1);
      return CAPI_EFAILED;
   }
   //AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: dma_config finished");
   //sleep(1);

   return capi_result;
}

capi_err_t capi_nxp_device_init_dma_buffer(capi_nxp_device_t *me_ptr, uint32_t *pb, uint32_t *pc)
{
   capi_err_t capi_result = CAPI_EOK;
   //struct capi_dai_data dai_data = me_ptr->dai_data;
   struct dai_config cfg;
   int ret;
   uint32_t channels = 2;
   uint32_t sample_bytes, frame_size, period_bytes, period_count;
   uint32_t align = 0, addr_align, buf_size;
   size_t heap_size;

   ret = dai_config_get(dai_dev, &cfg, me_ptr->dai_dir);
   if (ret) {
      AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Failed to get dai configuration");
      return CAPI_EFAILED;
   }

   //Get the DMA_ATTR_BUFFER_ADDRESS_ALIGNMENT
   ret = dma_get_attribute(dma_device, DMA_ATTR_BUFFER_ADDRESS_ALIGNMENT, &addr_align);
   if (ret < 0) {
      AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Failed to get dma address align");
      return CAPI_EFAILED;
   }

   //Get the DMA_ATTR_BUFFER_SIZE_ALIGNMENT
   ret = dma_get_attribute(dma_device, DMA_ATTR_BUFFER_SIZE_ALIGNMENT, &align);
   if (ret < 0) {
      AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Failed to get dma align");
      return CAPI_EFAILED;
   }

   //Calculate frame size, period size, and period count from the values in cfg
   sample_bytes = 4;
   frame_size = sample_bytes * channels;
   period_bytes = 384;
   period_count = 2;
   buf_size = ALIGN_UP(period_count * period_bytes, align);
   *pb = period_bytes;
   *pc = period_count;
   me_ptr->dma_buf_size = buf_size;

   //If dma buffer exists, change the dma_buffer size to the period size
   if (me_ptr->dma_dest_addr) {
      //Check if the dma_buffer is the same size as the period size, otherwise resize it

      //Use sys_heap_realloc?
   } else {
      heap_size = (size_t) buf_size + 64;
      sys_heap_init(&kernel_heap.heap, heapmem, heap_size);

      me_ptr->dma_src_addr = sys_heap_aligned_alloc(&kernel_heap.heap, align, buf_size);
      if (NULL == me_ptr->dma_src_addr) {
         AR_MSG(DBG_ERROR_PRIO,
               "CAPI_NXP_DEVICE: Failed to allocate DMA destination buffer");
         return CAPI_ENOMEMORY;
      }

      //memset(me_ptr->dma__addr, 0, buf_size);
   }

   return capi_result;

}

capi_err_t capi_nxp_device_set_dma_config(capi_nxp_device_t *me_ptr, uint32_t period_bytes, uint32_t period_count)
{
   capi_err_t capi_result = CAPI_EOK;
   struct dma_config *dma_cfg = NULL;
   struct dma_block_config *dma_block_cfg = NULL;
   struct dma_block_config *prev = NULL;
   int ret, channel;
   uint32_t fifo, burst_elems, max_block_count, buf_size;

   const struct dai_properties *props = dai_get_properties(dai_dev, me_ptr->dai_dir, 0);
   uint32_t hs = props->dma_hs_id;
   burst_elems = props->fifo_depth;
   fifo = props->fifo_address;
   channel = hs & GENMASK(7, 0);

   ret = dma_get_attribute(dma_device, DMA_ATTR_MAX_BLOCK_COUNT, &max_block_count);
   if (ret < 0 || !max_block_count) {
      AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Failed to get block count");
      capi_result = CAPI_EFAILED;
      goto end;
   }

   if (max_block_count < period_count) {
      AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: unsupported period count %d", period_count);
      buf_size = period_count * period_bytes;
      do {
         if (IS_ALIGNED(buf_size, max_block_count)) {
            period_count = max_block_count;
            period_bytes = buf_size / period_count;
            break;
         } else {
            AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: alignment error for buf_size = %d, block count = %d",
                 buf_size, max_block_count);
         }
      } while (--max_block_count > 0);
   }

   dma_cfg = posal_memory_malloc(sizeof(struct dma_config), (POSAL_HEAP_ID)me_ptr->heap_mem.heap_id);
   if (NULL == dma_cfg) {
      AR_MSG(DBG_ERROR_PRIO,
            "CAPI_NXP_DEVICE: Failed to allocate dma configuration pointer");
      capi_result = CAPI_ENOMEMORY;
      goto end;
   }

   dma_cfg->channel_direction = MEMORY_TO_PERIPHERAL;
   dma_cfg->source_data_size = 4;
   dma_cfg->dest_data_size = dma_cfg->source_data_size;
   dma_cfg->source_burst_length = 8;
   dma_cfg->dest_burst_length = dma_cfg->source_burst_length;
   dma_cfg->cyclic = 1;
   dma_cfg->user_data = me_ptr;
   dma_cfg->dma_callback = &capi_nxp_device_dma_callback; //Set to dma callback function
   dma_cfg->block_count = max_block_count;
   dma_cfg->dma_slot = (hs & GENMASK(15, 8)) >> 8;

   dma_block_cfg = posal_memory_malloc(sizeof(struct dma_block_config) * dma_cfg->block_count, (POSAL_HEAP_ID)me_ptr->heap_mem.heap_id);
   if (NULL == dma_cfg) {
      AR_MSG(DBG_ERROR_PRIO,
            "CAPI_NXP_DEVICE: Failed to allocate DMA blocks");
      capi_result = CAPI_ENOMEMORY;
      goto end;
   }

   dma_cfg->head_block = dma_block_cfg;
   for (int i = 0; i < dma_cfg->block_count; i++) {
      dma_block_cfg->dest_scatter_en = 0;
      dma_block_cfg->block_size = period_bytes;
      dma_block_cfg->source_address = me_ptr->dma_src_addr;
      dma_block_cfg->dest_address = fifo;
      dma_block_cfg->source_addr_adj = DMA_ADDR_ADJ_DECREMENT;
      dma_block_cfg->dest_addr_adj = DMA_ADDR_ADJ_INCREMENT;
      prev = dma_block_cfg;
      prev->next_block = dma_block_cfg++;
   }
   if (prev)
      prev->next_block = dma_cfg->head_block;

   me_ptr->z_config = dma_cfg;

   channel = dma_request_channel(dma_device, &channel);
   if (channel < 0) {
      AR_MSG(DBG_ERROR_PRIO,
            "CAPI_NXP_DEVICE: DMA request channel failed");
      capi_result = CAPI_EFAILED;
      goto end;
   }
   me_ptr->dma_chan_index = channel;

   ret = dma_config(dma_device, me_ptr->dma_chan_index, dma_cfg);
   if (ret < 0) {
      AR_MSG(DBG_ERROR_PRIO,
            "CAPI_NXP_DEVICE: dma_config failed");
      capi_result = CAPI_EFAILED;
      goto end;
   }

   end:
   if (capi_result != 0) {
      if (dma_block_cfg)
         posal_memory_free(dma_block_cfg);
      if (dma_cfg)
         posal_memory_free(dma_cfg);

      if (me_ptr->dma_src_addr)
         sys_heap_free(&kernel_heap.heap, me_ptr->dma_src_addr);
   }

   return capi_result;

}

void capi_nxp_device_reset(capi_nxp_device_t *me_ptr)
{
   if (me_ptr->z_config->head_block)
      posal_memory_free(me_ptr->z_config->head_block);

   if (me_ptr->z_config)
      posal_memory_free(me_ptr->z_config);

   if (me_ptr->dma_src_addr)
      sys_heap_free(&kernel_heap.heap, me_ptr->dma_src_addr);
}

void capi_nxp_device_stop(capi_nxp_device_t *me_ptr)
{
   //Stop the DMA and the DAI
   dma_stop (dma_device, me_ptr->dma_chan_index);
   dai_trigger(dai_dev, me_ptr->dai_dir, DAI_TRIGGER_STOP);
}

/*---------------------------------------------------------------------
  Function name: capi_nxp_device_dma_callback
  DESCRIPTION: Will be called by the DMA driver when it needs more data
  -----------------------------------------------------------------------*/

void capi_nxp_device_dma_callback(const struct device *dev, void *user_data,
                uint32_t channel, int status) {
   capi_nxp_device_t *me_ptr;
   //AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: callback enter");
   //sleep(1);

   if (user_data)
      me_ptr = (capi_nxp_device_t *) user_data;
   else {
      AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: User data is NULL");
      return;
   }

   //Trigger the waiting thread to tell the container we need more data
   if (me_ptr->i2s_driver.i2s_intf_state.trigger_signal_ptr)
      posal_signal_send((posal_signal_t)me_ptr->i2s_driver.i2s_intf_state.trigger_signal_ptr);
   else
      AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Callback is not set yet");
}

/*---------------------------------------------------------------------
  Function name: capi_nxp_device_process_sink
  DESCRIPTION: Processes an input buffer and generates an output buffer.
  -----------------------------------------------------------------------*/
capi_err_t capi_nxp_device_process_sink(capi_t *_pif, capi_stream_data_t *input[], capi_stream_data_t *output[])
{
   capi_err_t capi_result = CAPI_EOK;
   capi_nxp_device_t *me_ptr = (capi_nxp_device_t *)_pif;
   struct dma_status stat;
   int ret;
   uint16_t port = 0;
   uint16_t i = 0;
   uint32_t free_bytes;
   uint32_t copy_bytes;

   //AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: process_sink enter");
   //sleep(1);

   //Get data sizes and read/write buffer positions from the DMA
   ret = dma_get_status(dma_device, me_ptr->dma_chan_index, &stat);

   free_bytes = stat.free; //The amount of free bytes in the output DMA buffer

   //Choose the minimum of free bytes in output buffer or available bytes in input buffer to copy
   copy_bytes = MIN(free_bytes, input[port]->buf_ptr[i].actual_data_len);

   /* return if nothing to copy */
   if (!copy_bytes) {
      //AR_MSG(DBG_ERROR_PRIO,
      //      "CAPI_NXP_DEVICE: Bytes to copy is 0");
      dma_reload(dma_device, me_ptr->dma_chan_index, 0, 0, 0);    
      return 0;
   }
  
   memcpy(me_ptr->dma_src_addr,
           input[port]->buf_ptr[i].data_ptr,
           copy_bytes);

   dma_reload(dma_device, me_ptr->dma_chan_index, 0, 0, copy_bytes);

   return capi_result;
}


capi_err_t capi_nxp_device_process_get_properties(capi_nxp_device_t *me_ptr, capi_proplist_t *proplist_ptr)
{
   capi_err_t capi_result = CAPI_EOK;
   uint32_t fwk_extn_ids[1] = { 0 };
   fwk_extn_ids[0]          = FWK_EXTN_STM;

   capi_basic_prop_t mod_prop;
   mod_prop.init_memory_req = sizeof(capi_nxp_device_t);
   mod_prop.stack_size = NXP_DEVICE_STACK_SIZE;
   mod_prop.num_fwk_extns = 1;
   mod_prop.fwk_extn_ids_arr = fwk_extn_ids;
   mod_prop.is_inplace = 0;       // NA
   mod_prop.req_data_buffering = 0; // NA
   mod_prop.max_metadata_size = 0;  // NA

   capi_result = capi_cmn_get_basic_properties(proplist_ptr, &mod_prop);
   if (CAPI_EOK != capi_result)
   {
      AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Get common basic properties failed with result %lu", capi_result);
      return capi_result;
   }

   capi_prop_t *prop_ptr = proplist_ptr->prop_ptr;
   for (uint32_t i = 0; i < proplist_ptr->props_num; i++)
   {
      capi_buf_t *payload_ptr = &prop_ptr[i].payload;
      switch (prop_ptr[i].id)
      {
         case CAPI_INIT_MEMORY_REQUIREMENT:
         case CAPI_STACK_SIZE:
         case CAPI_NUM_NEEDED_FRAMEWORK_EXTENSIONS:
         case CAPI_NEEDED_FRAMEWORK_EXTENSIONS:
         case CAPI_OUTPUT_MEDIA_FORMAT_SIZE:
         case CAPI_IS_INPLACE:
         case CAPI_REQUIRES_DATA_BUFFERING:
         case CAPI_OUTPUT_MEDIA_FORMAT_V2:
         {
            break;
         }
         case CAPI_PORT_DATA_THRESHOLD:
         {
            if (NULL == me_ptr)
            {
               AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: null ptr while querying threshold");
               return CAPI_EBADPARAM;
            }

            // based on int samples per period calculation (frame size in bytes)
            uint32_t threshold_in_bytes = 8 *
                                 16 *
                                 2;

            capi_result = capi_cmn_handle_get_port_threshold(&prop_ptr[i], threshold_in_bytes);
            break;
         }
         case CAPI_INTERFACE_EXTENSIONS:
         {
            //capi_result = capi_hw_intf_cmn_update_intf_extn_status(payload_ptr);
            break;
         }
         default:
         {
            AR_MSG(DBG_HIGH_PRIO, "CAPI_NXP_DEVICE: Skipped Get Property for 0x%x. Not supported.", prop_ptr[i].id);
            capi_result |= CAPI_EUNSUPPORTED;
            break;
         }
      }
   }
   return capi_result;
}

/*------------------------------------------------------------------------
  Function name: capi_nxp_device_process_set_properties
  DESCRIPTION: Function to set the properties for the nxp_device module
 * -----------------------------------------------------------------------*/
capi_err_t capi_nxp_device_process_set_properties(capi_nxp_device_t *me_ptr, capi_proplist_t *proplist_ptr)
{
   capi_err_t capi_result = CAPI_EOK;
   int ret;
   struct dma_status stat;
   uint32_t free_bytes;
   //AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Set properties enter");

   if (NULL == me_ptr)
   {
      AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Set property received null property.");
      return CAPI_EBADPARAM;
   }

   capi_result = capi_cmn_set_basic_properties(proplist_ptr, &me_ptr->heap_mem, &me_ptr->cb_info, FALSE);

   if (CAPI_EOK != capi_result)
   {
      AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Set basic properties failed with result %lu", capi_result);
      return capi_result;
   }

   uint32_t i = 0;
   capi_prop_t *prop_ptr = proplist_ptr->prop_ptr;

   for (i = 0; i < proplist_ptr->props_num; i++)
   {
      capi_buf_t *payload_ptr = &prop_ptr[i].payload;
      switch (prop_ptr[i].id)
      {
         case CAPI_HEAP_ID:
         case CAPI_EVENT_CALLBACK_INFO:
         case CAPI_ALGORITHMIC_RESET:
         {
            break;
         }
         case CAPI_INPUT_MEDIA_FORMAT_V2:
         {
            capi_result = CAPI_EOK;

            break;
        }
        case CAPI_PORT_NUM_INFO:
        {
            if (payload_ptr->actual_data_len >= sizeof(capi_port_num_info_t))
            {
                capi_port_num_info_t *data_ptr = (capi_port_num_info_t *)payload_ptr->data_ptr;
                if (!(data_ptr->num_input_ports == 1 && data_ptr->num_output_ports == 0) &&
                !(data_ptr->num_input_ports == 0 && data_ptr->num_output_ports == 1))
                {
                  AR_MSG(DBG_ERROR_PRIO,
                        "CAPI_NXP_DEVICE: Invalid number of input = %d, number of output ports = %d.",
                        data_ptr->num_input_ports,
                        data_ptr->num_output_ports);
                  CAPI_SET_ERROR(capi_result, CAPI_EBADPARAM);
                }
            }
            else
            {
               AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Bad param size %lu", payload_ptr->actual_data_len);
               CAPI_SET_ERROR(capi_result, CAPI_ENEEDMORE);
            }
            break;
        }
        case CAPI_CUSTOM_PROPERTY:
        {
            capi_custom_property_t *cust_prop_ptr    = (capi_custom_property_t *)payload_ptr->data_ptr;
            void *                  cust_payload_ptr = (void *)(cust_prop_ptr + 1);
  
            switch (cust_prop_ptr->secondary_prop_id)
            {
               case FWK_EXTN_PROPERTY_ID_STM_TRIGGER:
               {
                  if (payload_ptr->actual_data_len < sizeof(capi_custom_property_t) + sizeof(capi_prop_stm_trigger_t))
                  {
                     AR_MSG(DBG_ERROR_PRIO,
                            "CAPI_NXP_DEVICE: Insufficient payload size for stm trigger %d",
                            payload_ptr->actual_data_len);
                     return CAPI_EBADPARAM;
                  }

                  //AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Setting STM trigger");
                  //sleep(1);
  
                  // Get end point info
                  capi_prop_stm_trigger_t *trig_ptr                    = (capi_prop_stm_trigger_t *)cust_payload_ptr;
                  me_ptr->i2s_driver.i2s_intf_state.trigger_signal_ptr = trig_ptr->signal_ptr;
                  me_ptr->i2s_driver.i2s_intf_state.i2s_intr_cnt_ptr   = trig_ptr->raised_intr_counter_ptr;
  
                  break;
               }
               case FWK_EXTN_PROPERTY_ID_STM_CTRL:
               {
                  //Start the DMA
                   // Allocate buffer for process call
                  //AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Setting STM CTRL enter");
                  //sleep(1);
                   uint32_t buf_size = me_ptr->dma_buf_size;

                  AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Starting the DMA");
                  ret = dma_start(dma_device, me_ptr->dma_chan_index);
                  if (ret < 0) {
                     AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: DMA start failed.");
                     return CAPI_EFAILED;
                  }

                  //Start the DAI
                  dai_trigger(dai_dev, 0, DAI_TRIGGER_START);

                  AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Finished starting DMA and DAI");
                  sleep(1);

                   /*if (me_ptr->i2s_driver.i2s_intf_state.trigger_signal_ptr)
                     posal_signal_send((posal_signal_t)me_ptr->i2s_driver.i2s_intf_state.trigger_signal_ptr);
                  else
                     AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Callback is not set yet");*/
                  break;
               }
               default:
               {
                  AR_MSG(DBG_HIGH_PRIO, "CAPI_NXP_DEVICE: Unknown Custom Property[%d]", cust_prop_ptr->secondary_prop_id);
                  capi_result |= CAPI_EUNSUPPORTED;
                  break;
               }
            } // inner switch - CUSTOM Properties
            break;
        }
        case CAPI_MODULE_INSTANCE_ID:
        {
            if (payload_ptr->actual_data_len >= sizeof(capi_module_instance_id_t))
            {
               capi_module_instance_id_t *data_ptr = (capi_module_instance_id_t *)payload_ptr->data_ptr;
               me_ptr->iid = data_ptr->module_instance_id;
               AR_MSG(DBG_LOW_PRIO,
                     "CAPI_NXP_DEVICE: This module-id 0x%08lX, instance-id 0x%08lX",
                     data_ptr->module_id,
                     me_ptr->iid);
            }
            else
            {
               AR_MSG(DBG_ERROR_PRIO,
                     "CAPI_NXP_DEVICE: Set, Param id 0x%lx Bad param size %lu",
                     (uint32_t)prop_ptr[i].id,
                     payload_ptr->actual_data_len);
               CAPI_SET_ERROR(capi_result, CAPI_ENEEDMORE);
            }
            break;
        }
        case CAPI_LOGGING_INFO:
        {
           if (payload_ptr->actual_data_len >= sizeof(capi_logging_info_t))
           {
              capi_logging_info_t *data_ptr = (capi_logging_info_t *)payload_ptr->data_ptr;
              me_ptr->log_id                = data_ptr->log_id;
              me_ptr->log_id_reserved_mask  = data_ptr->log_id_mask;
              AR_MSG(DBG_LOW_PRIO,
                     "CAPI_NXP_DEVICE: log-id 0x%08lX, mask 0x%08lX",
                     me_ptr->log_id,
                     me_ptr->log_id_reserved_mask);
           }
           else
           {
              AR_MSG(DBG_ERROR_PRIO,
                     "CAPI_NXP_DEVICE: Set, Param id 0x%lx Bad param size %lu",
                     (uint32_t)prop_ptr[i].id,
                     payload_ptr->actual_data_len);
              CAPI_SET_ERROR(capi_result, CAPI_ENEEDMORE);
           }
           break;
        }
        default:
        {
            AR_MSG(DBG_ERROR_PRIO, "CAPI_NXP_DEVICE: Skipping set prop, unsupported param[%d]", prop_ptr[i].id);
            capi_result |= CAPI_EUNSUPPORTED;
            continue;
        }
      } /* Outer switch - Generic CAPI Properties */
   } /* Loop all properties */

   return capi_result;
}
