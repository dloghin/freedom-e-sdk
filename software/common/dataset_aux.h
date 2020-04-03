/*
 * dataset_aux.h
 *
 *  Created on: 20 Mar 2020
 *      Author: dumi
 */

#ifndef SOFTWARE_COMMON_DATASET_AUX_H_
#define SOFTWARE_COMMON_DATASET_AUX_H_

#include "posit.h"
#include "dataset_posit.h"

void init_dataset() {
#if (defined WITH_POSIT_8 || defined WITH_POSIT_16 || defined WITH_POSIT_32)
	int i, j;
	for (i = 0; i < TRAINING_DATA_LENGTH; i++)
		for (j = 0; j < INPUT_LENGTH; j++)
			*((uint32_t*)&g_training_data_X[i][j]) = posit_g_training_data_X[i][j];
	for (i = 0; i < TRAINING_DATA_LENGTH; i++)
		*((uint32_t*)&g_training_data_Y[i]) = posit_g_training_data_Y[i];
	for (j = 0; j < INPUT_LENGTH; j++)
		*((uint32_t*)&g_test_data[j]) = posit_g_test_data[j];
#endif
}
#endif /* SOFTWARE_COMMON_DATASET_AUX_H_ */
