#pragma once

#include "posit.h"

int_t linear_index(int_t i_rows, int_t i_columns, int_t i_row_index, int_t i_column_index) {
	return i_row_index * i_columns + i_column_index;
}

/**
 * The element-wise difference of two vectors.
 * Inputs:
 * 	i_vector_A: a vector
 * 	i_vector_B: a vector
 *  i_length : size of the vectors
 * Output:
 *  o_vector_C: a vector C=A-B
 */
void vector_difference (element_t *i_vector_A,
		element_t *i_vector_B,
		int_t i_length,
		element_t *o_vector_C) {

	int_t index;
	for (index = 0; index < i_length; index++) {
		o_vector_C[index] = i_vector_A[index] - i_vector_B[index];
	}
}

/**
 * Calculate the product of a matrix and a vector
 * Inputs:
 * 	i_matrix_A: a matrix i_rows x i_columns
 * 	i_rows: number of rows for A
 *  i_columns : number of columns for A
 *  i_vector_B: a vector of size i_columns x 1
 * Output:
 *  o_vector_C: a vector C=A*B size i_rows x 1
 */
void matrix_multiply_vector(element_t *i_matrix_A,
		int_t i_rows,
		int_t i_columns,
		element_t *i_vector_B,
		element_t *o_vector_C) {

	int_t index, jindex;
	for(index = 0; index < i_rows; index++) {
		o_vector_C[index] = zero;
		for(jindex = 0; jindex < i_columns; jindex++) {
			o_vector_C[index] = o_vector_C[index] +
					i_matrix_A[linear_index(i_rows, i_columns, index, jindex)] * i_vector_B[jindex];
		}
	}
}

/**
 * Product of a vector transpose and a matrix as a vector
 * Inputs:
 * 	i_vector_A: a matrix i_rows
 * 	i_matrix_B: a matrix i_rows x i_columns
 * 	i_rows: number of rows for i_matrix_A
 *  i_columns : number of columns for i_matrix_A
 * Output:
 *  o_matrix_C: a matrix C=A.T*B size 1 x i_columns
 */
void vector_T_multiply_matrix(element_t *i_vector_A,
		element_t *i_matrix_B,
		int_t i_rows,
		int_t i_columns,
		element_t o_vector_C[]) {

	int_t index, jindex;
	for(index = 0; index < i_columns; index++) {
		o_vector_C[index] = zero;
		for(jindex = 0; jindex < i_rows; jindex++) {
			o_vector_C[index] = o_vector_C[index]
										   + i_vector_A[jindex] * i_matrix_B[linear_index(i_rows, i_columns, jindex, index)];
		}
	}
}

/**
 * Product of a vector with a scalar
 * Inputs:
 * 	i_vector_A: a matrix i_length
 * 	i_length: length of i_vector_A
 * 	i_scalar: scalar
 * Output:
 *  o_vector_B: a vector B=A*scalar size i_length
 */
void vector_multiply_scalar(
		element_t *i_vector_A,
		int_t i_length,
		element_t i_scalar,
		element_t *o_vector_B) {

	int_t index;
	for (index = 0; index < i_length; index++) {
		o_vector_B[index] = i_vector_A[index] * i_scalar;
	}
}

/**
 * Product of the transpose of a matrix and a matrix
 * Inputs:
 * 	i_matrix_A: a matrix i_middle x i_rows
 * 	i_matrix_B: a matrix i_middle x i_columns
 * 	i_rows: number of rows for o_matrix_C and columns for i_matrix_A
 *  i_columns : number of columns for i_matrix_B and o_matrix_C
 *  i_middle: number of rows for i_matrix_A and i_matrix_B
 * Output:
 *  o_matrix_C: a matrix C=A.T*B size i_rows x i_columns
 */
void matrix_T_multiply_matrix(element_t *i_matrix_A,
		element_t *i_matrix_B,
		int_t i_rows,
		int_t i_middle,
		int_t i_columns,
		element_t *o_matrix_C) {

	int_t index, jindex, kindex;
	for (index = 0; index < i_rows; index++) {
		for (jindex = 0; jindex < i_columns; jindex++) {
			o_matrix_C[linear_index(i_rows, i_columns, index, jindex)] = zero;
			for (kindex = 0; kindex < i_middle; kindex++) {
				o_matrix_C[linear_index(i_rows, i_columns, index, jindex)] =
						o_matrix_C[linear_index(i_rows, i_columns, index, jindex)]
								   + i_matrix_A[linear_index(i_middle, i_rows, kindex, index)]
												* i_matrix_B[linear_index(i_middle, i_columns, kindex, jindex)];
			}
		}
	}
}

/**
 * Calculate the product of the transpose of a matrix and a vector
 * Inputs:
 * 	i_matrix_A: a matrix i_rows x i_columns
 * 	i_rows: number of rows for A
 *  i_columns : number of columns for A
 *  i_vector_B: a vector of size i_rows x 1
 * Output:
 *  o_vector_C: a vector C=A*B size i_columns x 1
 */
void matrix_T_multiply_vector(element_t *i_matrix_A,
		int_t i_rows,
		int_t i_columns,
		element_t *i_vector_B,
		element_t *o_vector_C) {

	int_t index, jindex;
	for(index = 0; index < i_columns; index++) {
		o_vector_C[index] = zero;
		for(jindex = 0; jindex < i_rows; jindex++) {
			o_vector_C[index] = o_vector_C[index]
										   + i_matrix_A[linear_index(i_rows, i_columns, jindex, index)] * i_vector_B[jindex];
		}
	}
}

/**
 * Product of transpose of a vector with a second vector resulting an element
 * Inputs:
 * 	i_vector_A: a vector of  i_length
 * 	i_vector_B: a vector of  i_length
 *  i_length : size of the vectors
 * Output:
 *  o_C: a matrix C=A.T*B size 1
 */
void vector_T_multiply_vector(element_t *i_vector_A,
		element_t i_vector_B[],
		int_t i_length,
		element_t *o_C) {

	int_t index;
	*o_C = zero;
	for(index = 0; index < i_length; index++) {
		*o_C = *o_C + i_vector_A[index] * i_vector_B[index];
	}
}

/*
* Matrix cofactor
* Inputs:
* 	i_matrix_A: a matrix i_size x i_size
* 	i_size: the size of the matrix
* 	i_row: the row to eliminate
* 	i_column: the column to eliminate
* Output:
*  o_matrix_B: a matrix (i_size-1) x (i_size-1)
* B is A without the given row and column
*/
void matrix_cofactor(element_t *i_matrix_A,
		int_t i_size,
		int_t i_row,
		int_t i_column,
		element_t *o_matrix_B) {

    int_t index, jindex, rindex, pindex;
    rindex = 0;
    pindex = 0;
    for(index = 0; index < i_size; index++) {
        for(jindex = 0; jindex < i_size; jindex++) {
            if(index != i_row && jindex != i_column) {
                o_matrix_B[linear_index(i_size - 1, i_size - 1, rindex, pindex)] = i_matrix_A[linear_index(i_size, i_size, index, jindex)];
                pindex = pindex + 1;
                if(pindex == (i_size-1)) {
                    pindex = 0;
                    rindex = rindex + 1;
                }
            }
        }
    }
}

/*
* Crout matrix decomposition
* https://en.wikipedia.org/wiki/Crout_matrix_decomposition
* Inputs:
* 	i_matrix_A: a matrix i_size x i_size
* 	i_size: the size of the matrix
* Output:
*  o_matrix_L: a matrix i_size x i_size
*  o_matrix_U: a matrix i_size x i_size
* A=L*U
*/
int_t matrix_crout_decomposition(element_t *i_matrix_A,
		int_t i_size,
		element_t *o_matrix_L,
		element_t *o_matrix_U) {

    int_t index, jindex, kindex;
    element_t sum;

    //initialize matrix
	for (index = 0; index < i_size; index++) {
	    for (jindex = 0; jindex < i_size; jindex++) {
		    o_matrix_U[linear_index(i_size, i_size, index, jindex)] = zero;
		    o_matrix_L[linear_index(i_size, i_size, index, jindex)] = zero;
        }
	}

    //initialize diagonal for U
	for (index = 0; index < i_size; index++) {
		o_matrix_U[linear_index(i_size, i_size, index, index)] = one;
	}

	for (jindex = 0; jindex < i_size; jindex++) {
		for (index = jindex; index < i_size; index++) {
			sum = zero;
			for (kindex = 0; kindex < jindex; kindex++) {
				sum = sum +
						o_matrix_L[linear_index(i_size, i_size, index, kindex)] *
						o_matrix_U[linear_index(i_size, i_size, kindex, jindex)];
			}
			o_matrix_L[linear_index(i_size, i_size, index, jindex)] = i_matrix_A[linear_index(i_size, i_size, index, jindex)] - sum;
		}

		for (index = jindex; index < i_size; index++) {
			sum = zero;
			for(kindex = 0; kindex < jindex; kindex++) {
				sum = sum +
						o_matrix_L[linear_index(i_size, i_size, jindex, kindex)] *
						o_matrix_U[linear_index(i_size, i_size, kindex, index)];
			}
			if (o_matrix_L[linear_index(i_size, i_size, jindex, jindex)] == zero) {
				return -1;
			}
			o_matrix_U[linear_index(i_size, i_size, jindex, index)] =
					(i_matrix_A[linear_index(i_size, i_size, jindex, index)] - sum) /
					o_matrix_L[linear_index(i_size, i_size, jindex, jindex)];
		}
	}

    return 0;
}

/*
* Matrix determinant
* Inputs:
* 	i_matrix_A: a matrix i_size x i_size
* 	i_size: the size of the matrix
* Output:
*  o_matrix_L: a matrix i_size x i_size
*  o_matrix_U: a matrix i_size x i_size
*  o_determinant: determinant of A  = det(A)
* A=L*U
* det(A)=det(L)*det(U), det(U) = 1 => det(A)=det(L)
*/
int_t matrix_determinant(element_t *i_matrix_A,
		int_t i_size,
		element_t *o_matrix_L,
		element_t *o_matrix_U,
		element_t *o_determinant) {

    int_t return_value;
    return_value = matrix_crout_decomposition(i_matrix_A, i_size,
                                              o_matrix_L, o_matrix_U);
    if(return_value != 0) {
        return return_value;
    }
    *o_determinant = one;
    int_t index;
	for (index = 0; index < i_size; index++) {
		*o_determinant = *o_determinant * o_matrix_L[linear_index(i_size, i_size, index, index)];
	}
    return 0;
}

/*
* Matrix co-factor
* Inputs:
* 	i_matrix_A: a matrix i_size x i_size
* 	i_size: the size of the matrix
* 	i_matrix_aux1: a matrix i_size x i_size use for co-factor
* 	i_matrix_aux2: a matrix i_size x i_size use for L calculation
* 	i_matrix_aux3: a matrix i_size x i_size use for U calculation
* Output:
*  o_matrix_B: a matrix i_size x i_size B=A^* =adj(A)
*/
int_t matrix_adjoint(element_t *i_matrix_A,
		int_t i_size,
		element_t *i_matrix_aux1,
		element_t *i_matrix_aux2,
		element_t *i_matrix_aux3,
		element_t *o_matrix_B) {

    int_t sign = 1;
    int_t index, jindex;
    int_t return_value;
    element_t determinant;

    for(index = 0; index < i_size; index++) {
        for(jindex = 0; jindex < i_size; jindex++) {
            matrix_cofactor(i_matrix_A, i_size, index, jindex, i_matrix_aux1);
            sign = ((index+jindex)%2==0)? 1: -1;
            return_value = matrix_determinant(i_matrix_aux1, i_size - 1, i_matrix_aux2,
                                              i_matrix_aux3, &determinant);
            if(return_value != 0) {
                return return_value;
            }
            o_matrix_B[linear_index(i_size, i_size, index, jindex)] = sign * determinant;
        }
    }
    return 0;
}

/*
 * Matrix inverse
 * Inputs:
 * 	i_matrix_A: a matrix i_size x i_size
 * 	i_size: the size of the matrix
 * 	i_matrix_aux1: a matrix i_size x i_size use for co-factor
 * 	i_matrix_aux2: a matrix i_size x i_size use for L calculation
 * 	i_matrix_aux3: a matrix i_size x i_size use for U calculation
 * Output:
 *  o_matrix_B: a matrix i_size x i_size B=A^-1
 */
int_t matrix_inverse(element_t *i_matrix_A,
		int_t i_size,
		element_t *i_matrix_aux1,
		element_t *i_matrix_aux2,
		element_t *i_matrix_aux3,
		element_t *o_matrix_B) {

	int_t index, jindex;
	int_t return_value;
	element_t determinant;
	return_value = matrix_determinant(i_matrix_A, i_size, i_matrix_aux2,
			i_matrix_aux3, &determinant);

	printf("det %x\n", *((uint32_t*)&determinant));
	printf("det %f\n", determinant);

	if(return_value != 0) {
		return return_value;
	}
	return_value = matrix_adjoint(i_matrix_A, i_size, i_matrix_aux1,
			i_matrix_aux2, i_matrix_aux3,
			o_matrix_B);
	if(return_value != 0) {
		return return_value;
	}

	/*
	if(determinant < eps && determinant > minus_eps) {
		return -5;
	}
	*/

	for(index = 0; index < i_size; index++) {
		for(jindex = 0; jindex < i_size; jindex++) {
			o_matrix_B[linear_index(i_size, i_size, index, jindex)] = o_matrix_B[linear_index(i_size, i_size, index, jindex)] / determinant;
		}
	}

	return 0;
}
