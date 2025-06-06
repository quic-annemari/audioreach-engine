/* ======================================================================== */
/**
@file capi_pcm_cnv.h

   Header file to implement the Common Audio Post Processor Interface
   for Tx/Rx Tuning PCM_CNV block
*/

/* =========================================================================
   Copyright (c) Qualcomm Innovation Center, Inc. All Rights Reserved.
   SPDX-License-Identifier: BSD-3-Clause-Clear
  ========================================================================== */

/*------------------------------------------------------------------------
 * Include files
 * -----------------------------------------------------------------------*/
#ifndef CAPI_PCM_CNV_H
#define CAPI_PCM_CNV_H

#include "capi.h"
#include "ar_defs.h"

#ifdef __cplusplus
extern "C" {
#endif /*__cplusplus*/

capi_err_t capi_pcm_cnv_get_static_properties(capi_proplist_t *init_set_properties,
                                                    capi_proplist_t *static_properties);
capi_err_t capi_pcm_cnv_init(capi_t *_pif, capi_proplist_t *init_set_properties);

#ifdef __cplusplus
}
#endif /*__cplusplus*/

#endif // CAPI_PCM_CNV_H
