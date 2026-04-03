/*
 * Copyright (c) Qualcomm Technologies, Inc. and/or its subsidiaries.
 * SPDX-License-Identifier: BSD-3-Clause
 */

#include <stdio.h>
#include <posal.h>
#include "posal.h"
#include <spf_main.h>
#include <gpr_api.h>

int main(void)
{
    int rc = 0;
    posal_init();
	
    rc = gpr_init();
    if (0 != rc) {
	    printf("gpr_init() failed with status %d\n", rc);
        return rc;
    }
	
    rc = spf_framework_pre_init();
    if (0 != rc) {
	    printf("spf_framework_pre_init() failed with status %d\n", rc);
        return rc;
    }
	
    rc = spf_framework_post_init();
    if (0 != rc) {
		printf("spf_framework_post_init() failed with status %d\n", rc);
        return rc;
    }

	printf("spf framework initialized.\n", rc);

	/*apm_test_launch();
	printf("after apm_test_launch() \n");


	spf_framework_pre_deinit();
	printf("after  spf_framework_pre_deinit call\n");

    spf_framework_post_deinit();
	printf("after  spf_framework_post_deinit call\n");
	
	printf("before gpr_deinit call\n");
	gpr_deinit();
	printf("after gpr_deinit call\n");
	
	posal_deinit();
	printf("After posal deinit call\n");*/
	    
    return 0;
}