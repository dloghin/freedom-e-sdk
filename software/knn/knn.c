//KNN
// https://www.geeksforgeeks.org/k-nearest-neighbours/
// we use euclidian_distance_square instead of euclidian_distance
// Ciocirlan Stefan-Dan 14.05.2019

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "../common/posit.h"
#include "../common/dataset.h"
#include "../common/dataset_aux.h"
#include "../common/perf.h"

#define PFDEBUG

#define K_VALUE 7

element_t g_distances[K_VALUE];
element_t g_values[K_VALUE];

/* Calculate the euclidian distance between vectors X and Y with a length given
 */
element_t euclidian_distance_square (element_t X[], element_t Y[], int_t length) {
	int_t index;
	element_t sum = 0;
	for(index = 0; index < length; index++) {
		sum = sum + (X[index]-Y[index])*(X[index]-Y[index]);
	}
	return sum;
}

/* Classify a given point test_data_X with a KNN network with training_data_length
 * points with position training_data_X and values training_data_Y with given k value
 * for the bumber of neighbors. It also calculate de distances to this neighbors with
 * their values
 */
element_t classifyAPoint(element_t training_data_X[TRAINING_DATA_LENGTH][INPUT_LENGTH], element_t training_data_Y[],
		int_t training_data_length, int_t input_length,
		element_t test_data_X[], int_t k, element_t distances[],
		element_t values[]) {
	int_t index, jindex, kindex;
	element_t distance;
	int_t max_frequency;
	int_t return_value;
	int_t current_frequency;

	//initialise distances and values
	for(kindex = 0; kindex < k; kindex++) {
		distances[kindex] = 9999;
		values[kindex] = minus_one;
	}
	/*
	 * Find the closest k points from training data points to our test data point
	 * distance is a sort vector of points distance with smallest distance on position 0
	 * values is the value of the points above
	 */
	for (index = 0; index < training_data_length; index++) {
		//find the current distance
		distance = euclidian_distance_square(training_data_X[index], test_data_X, input_length);
		//find if its position in the distance vector
		for(kindex = 0; kindex < k; kindex++) {
			if(distance < distances[kindex]) {
				//move all other points on position to the right
				for(jindex = k - 1; jindex > kindex; jindex--) {
					distances[jindex] = distances[jindex-1];
					values[jindex] = values[jindex-1];
				}
				//put the point on its position
				distances[kindex] = distance;
				values[kindex] = training_data_Y[index];
				break;
			}
		}
	}

	//calculate the value with the higher frecuency in out k points
	max_frequency = 1;
	return_value = values[0];
	for(kindex = 0; kindex < k; kindex++) {
		current_frequency = 1;
		/*find other points with the same value
		going from kindex+1 because if the value was before the current
		point than was already counted
		 */
		for(jindex = kindex + 1; jindex < k; jindex++) {
			if(values[kindex] == values[jindex]) {
				current_frequency = current_frequency + 1;
			}
		}
		//if it is the the new max frequency value
		if(current_frequency > max_frequency) {
			return_value = values[kindex];
			max_frequency = current_frequency;
		}
	}

	return return_value; 
} 

int main() 
{
	init_constants();
	init_dataset();

	unsigned long long startc = read_cycles();
	int_t answer = classifyAPoint(g_training_data_X, g_training_data_Y,
			TRAINING_DATA_LENGTH, INPUT_LENGTH,
			g_test_data, K_VALUE, g_distances,
			g_values);
	unsigned long long endc = read_cycles();
	endc = endc - startc;
#ifdef PFDEBUG
	printf("Cycles %lu\n", endc);
	printf("Result: %d\n", answer);
#endif
	return 0;
}

