// CT or decision treee with gini index
// https://www.geeksforgeeks.org/decision-tree-introduction-example/
// https://medium.com/@rishabhjain_22692/decision-trees-it-begins-here-93ff54ef134
// Ciocirlan Stefan-Dan 15.05.2019

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "../common/posit.h"
#include "../common/dataset.h"
#include "../common/dataset_aux.h"
#include "../common/perf.h"

#define PFDEBUG

#define INPUT_LENGTH_SQUARE INPUT_LENGTH*INPUT_LENGTH
#define OUTPUT_GROUPS_LENGTH 3

int_t g_tree_atribute_index[INPUT_LENGTH_SQUARE];
int_t g_tree_atribute_value[INPUT_LENGTH_SQUARE];
int_t g_tree_class[INPUT_LENGTH_SQUARE];
int_t g_tree_atributes_remained[INPUT_LENGTH_SQUARE][INPUT_LENGTH];
int_t g_count_values[OUTPUT_GROUPS_LENGTH];
int_t g_training_data_Y_CT[TRAINING_DATA_LENGTH];
int_t g_training_position_tree[TRAINING_DATA_LENGTH];


/*
 * converts to a format that works with CT by grouping interval of values
 * at the same value
 */
void convert_to_CT(element_t training_data_Y[], int_t training_data_length,
		int_t training_data_Y_CT[], int_t output_groups_length) {
	int_t jindex, kindex;
	element_t min;
	element_t max;
	element_t value;


	/*
	 * find the minimum and maximum value and creates output_groups_length intervals
	 * if the value is less than the right value of the interval than the element
	 * can be in that interval
	 */
	min = training_data_Y[0];
	max = training_data_Y[0];
	for(jindex = 1; jindex < training_data_length; jindex++) {
		if(min > training_data_Y[jindex]) {
			min = training_data_Y[jindex];
		}
		if(max < training_data_Y[jindex]) {
			max = training_data_Y[jindex];
		}
	}

	for(jindex = 0; jindex < training_data_length; jindex++) {
		training_data_Y_CT[jindex] = output_groups_length - 1;
	}

	for(kindex = output_groups_length - 2; kindex > -1; kindex--) {
		value = min + ((max - min) * (kindex + 1)) / output_groups_length;
		for(jindex = 0; jindex < training_data_length; jindex++) {
			if(training_data_Y[jindex] <= value) {
				training_data_Y_CT[jindex] = kindex;
			}
		}
	}

}

/*
 * creates the tree recursive (binary tree)
 * index - the node index in the tree
 * training_data_X - the training data X values training_data_length X input_length size
 * training_data_Y - the training data Y values training_data_length size
 * tree_atrtibutes_remained - matrix with rows representin the index of the node in the tree
 *                            and the column represent the atribute of the training_data_X. The
 *                            value 1 represent that the atribute is avaible for the node, otherwise
 *                            the values is 0. The size is input_length^2 X input_length
 * output_groups_length - how many values are for training data Y
 * tree_atribute_index - the atribute chosen for the node. Size is input_length^2
 * tree_atribute_value - the atribute value of the node. Size is input_length^2
 * tree_class - the class of the node. Only for leaf nodes. Size is input_length^2
 * training_position_tree - the index position of the node in the tree where the training data is currently
 *                          Size is  training_data_length.
 * count_values - use for counting how many times a value from the training data Y is in a subpart.
 *                Size is output_groups_length
 */
void creating_tree(int_t index, element_t training_data_X[][INPUT_LENGTH], int_t training_data_length,
		int_t input_length, int_t training_data_Y[], int_t tree_atributes_remained[][INPUT_LENGTH],
		int_t output_groups_length, int_t tree_atribute_index[], int_t tree_atribute_value[],
		int_t tree_class[], int_t training_position_tree[], int_t count_values[]) {
	int_t jindex, kindex, rindex;

	/*
	 * If it is only on value on the training point for this node
	 */
	int_t sem;
	sem = 1;
	for(jindex = 0; jindex < training_data_length; jindex++) {
		if(training_position_tree[jindex] == index) {
			for(rindex = jindex + 1; rindex < training_data_length; rindex++) {
				if(training_position_tree[rindex] == index) {
					if(training_data_Y[jindex] != training_data_Y[rindex]) {
						sem = 0;
						break;
					}
				}
			}
			break;
		}
	}

	if(sem) {
		for(jindex = 0; jindex < training_data_length; jindex++) {
			if(training_position_tree[jindex] == index) {
				tree_class[index] = training_data_Y[jindex];
				return;
			}
		}
	}

	/*
	 * if there is no more input atributes to separate data
	 */
	int_t current_frequency;
	int_t max_frequency;
	sem = 1;
	for(jindex = 0; jindex < input_length; jindex++) {
		if(tree_atributes_remained[index][jindex]) {
			sem = 0;
		}
	}
	if(sem) {
		max_frequency = 0;
		for(jindex = 0; jindex < training_data_length; jindex++) {
			if(training_position_tree[jindex] == index) {
				current_frequency = 1;
				for(rindex = jindex + 1; rindex < training_data_length; rindex++) {
					if(training_position_tree[rindex] == index) {
						if(training_data_Y[jindex] == training_data_Y[index]) {
							current_frequency = current_frequency + 1;
						}
					}
				}
				if(current_frequency > max_frequency) {
					tree_class[index] = training_data_Y[jindex];
					max_frequency = current_frequency;
				}
			}
		}
		return;
	}

	/*
	 * Find the next atribute to separate the training data
	 */
	int_t current_kindex = -1;
	element_t max = zero;
	element_t mean;
	element_t total_elements;
	element_t current_count;
	element_t aux_g;
	element_t g;
	tree_class[index] = -1;
	for(rindex = 0; rindex < input_length; rindex++) {
		// only for the atributes that can be used for this node
		if(tree_atributes_remained[index][rindex]) {
			mean = zero;
			total_elements = zero;
			/*
			 * calculate the mean for that atribute with the training data
			 * for this node. Also count how many training data points are
			 */
			for(jindex = 0; jindex < training_data_length; jindex++) {
				if(training_position_tree[jindex] == index) {
					mean = mean + training_data_X[jindex][rindex];
					total_elements = total_elements + one;
				}
			}
			mean = mean / total_elements;
			/*
			 * Calculate the gini index for the atribute
			 * separete data that have the X value bigger than the mean in
			 * one group and calculate the gini index
			 * partial gini index = sum for i in output_groups_length of square of
			 *                      probability of the value
			 * probability of the value is the number of trainning points that have the value on Y
			 *                           under the total number of training points
			 *
			 * Do the same for the traiing points smaller than the mean.
			 *
			 * the final gini index =  partial gini index (bigger thean mean) * number of points (bigger than mean)
			 *                          / total number of points for the node +
			 *                          partial gini index (lower thean mean) * number of points (lower than mean)
			 *                          / total number of points for the node
			 * the real gini index is  1 minus the gini index in this program
			 */
			current_count = zero;
			for(kindex = 0; kindex < output_groups_length; kindex++) {
				count_values[kindex] = 0;
			}
			for(jindex = 0; jindex < training_data_length; jindex++) {
				if(training_position_tree[jindex] == index) {
					if(training_data_X[jindex][rindex] >= mean) {
						current_count = current_count + one;
						count_values[training_data_Y[jindex]] = count_values[training_data_Y[jindex]] + 1;
					}
				}
			}
			aux_g = zero;
			for(kindex = 0; kindex < output_groups_length; kindex++) {
				aux_g = aux_g + (count_values[kindex]* count_values[kindex]) / (current_count*current_count);
			}
			g = aux_g*current_count/total_elements;
			current_count = zero;
			for(kindex = 0; kindex < output_groups_length; kindex++) {
				count_values[kindex] = 0;
			}
			for(jindex = 0; jindex < training_data_length; jindex++) {
				if(training_position_tree[jindex] == index) {
					if(training_data_X[jindex][rindex] < mean) {
						current_count = current_count + one;
						count_values[training_data_Y[jindex]] = count_values[training_data_Y[jindex]] + 1;
					}
				}
			}
			aux_g = zero;
			for(kindex = 0; kindex < output_groups_length; kindex++) {
				aux_g = aux_g + (count_values[kindex]* count_values[kindex]) / (current_count*current_count);
			}
			g = g + aux_g*current_count/total_elements;
			/* if the gini index is bigger than max save the new atribute
			 * index and value for this node and update the max
			 *
			 */
			if(g > max) {
				max = g;
				current_kindex = rindex;
				tree_atribute_index[index] = current_kindex;
				tree_atribute_value[index] = mean;
			}
		}
	}

	// transmiting the remained atributes to sons
	for(rindex = 0; rindex < input_length; rindex++) {
		tree_atributes_remained[2*index+1][rindex] = tree_atributes_remained[index][rindex];
		tree_atributes_remained[2*index+2][rindex] = tree_atributes_remained[index][rindex];
	}
	tree_atributes_remained[2*index+1][current_kindex] = 0;
	tree_atributes_remained[2*index+2][current_kindex] = 0;

	// transmiting the data to the sons
	for(jindex = 0; jindex < training_data_length; jindex++) {
		if(training_position_tree[jindex] == index) {
			if(training_data_X[jindex][current_kindex] >= mean) {
				training_position_tree[jindex] = 2*index + 2;
			} else {
				training_position_tree[jindex] = 2*index + 1;
			}
		}
	}
	//create right son
	creating_tree(index*2+2, training_data_X, training_data_length,
			input_length, training_data_Y, tree_atributes_remained,
			output_groups_length, tree_atribute_index, tree_atribute_value,
			tree_class, training_position_tree, count_values);
	//creating left son
	creating_tree(index*2+1, training_data_X, training_data_length,
			input_length, training_data_Y, tree_atributes_remained,
			output_groups_length, tree_atribute_index, tree_atribute_value,
			tree_class, training_position_tree, count_values);
}

/*
 * Classify test_data by going through decision three
 */
int_t classify_CT(int_t index, int_t tree_atribute_index[],
		int_t tree_atribute_value[], int_t tree_class[],
		element_t test_data[]) {
	//if it is a leaf node it has class
	if(tree_class[index] != -1) {
		return tree_class[index];
	}
	/* if the value for the node atribute is higher than the mean go right,
	 * go left
	 */
	if(test_data[tree_atribute_index[index]] >= tree_atribute_value[index]) {
		return classify_CT(2*index+2, tree_atribute_index,
				tree_atribute_value, tree_class,
				test_data);
	} else {
		return classify_CT(2*index+1, tree_atribute_index,
				tree_atribute_value, tree_class,
				test_data);
	}
}


int main(void)
{   
	init_constants();
	init_dataset();

	unsigned long long startc = read_cycles();
	convert_to_CT(g_training_data_Y, TRAINING_DATA_LENGTH,
			g_training_data_Y_CT, OUTPUT_GROUPS_LENGTH);
	creating_tree(0, g_training_data_X, TRAINING_DATA_LENGTH,
			INPUT_LENGTH, g_training_data_Y_CT, g_tree_atributes_remained,
			OUTPUT_GROUPS_LENGTH, g_tree_atribute_index, g_tree_atribute_value,
			g_tree_class, g_training_position_tree, g_count_values);
	int_t answer = classify_CT(0, g_tree_atribute_index, g_tree_atribute_value,
			g_tree_class, g_test_data);
	unsigned long long endc = read_cycles();
	endc = endc - startc;
#ifdef PFDEBUG
	printf("Cycles %lu\n", endc);
	printf("Result: %d\n", answer);
#endif  // PFDEBUG
	return 0;
}
