#ifndef _COMMON_H_
#define _COMMON_H_ 1

#define WDEBUG

#ifdef PDEBUG
#include <stdio.h>
#endif
#include <stdio.h>

#include <stdint.h>
#include <stddef.h>


#define FP32

#ifdef FP32
typedef float element_t;
#else
typedef double element_t;
#endif

// data blob
extern element_t data[];

// constants
extern element_t zero;
extern element_t one;
extern element_t minus_one;
extern element_t e;

// pointers to data in objects
/*
#if (defined WITH_POSIT_32)
extern const element_t* _binary_relu3_in_bin_posit_32_3_start;
extern const element_t* _binary_prob_sum_multiplier_bin_posit_32_3_start;
extern const element_t* _binary_ip1_weights_bin_posit_32_3_start;
extern const element_t* _binary_ip1_blob_bin_posit_32_3_start;
extern const element_t* _binary_ip1_bias_mul_bin_posit_32_3_start;
#elif (defined WITH_POSIT_16)
extern const element_t* _binary_relu3_in_bin_posit_16_2_start;
extern const element_t* _binary_prob_sum_multiplier_bin_posit_16_2_start;
extern const element_t* _binary_ip1_weights_bin_posit_16_2_start;
extern const element_t* _binary_ip1_blob_bin_posit_16_2_start;
extern const element_t* _binary_ip1_bias_mul_bin_posit_16_2_start;
#elif (defined WITH_POSIT_8)
extern const element_t* _binary_relu3_in_bin_posit_8_1_start;
extern const element_t* _binary_prob_sum_multiplier_bin_posit_8_1_start;
extern const element_t* _binary_ip1_weights_bin_posit_8_1_start;
extern const element_t* _binary_ip1_blob_bin_posit_8_1_start;
extern const element_t* _binary_ip1_bias_mul_bin_posit_8_1_start;
#else
*/
extern const element_t* _binary_relu3_in_bin_start;
extern const element_t* _binary_prob_sum_multiplier_bin_start;
extern const element_t* _binary_ip1_weights_bin_start;
extern const element_t* _binary_ip1_blob_bin_start;
extern const element_t* _binary_ip1_bias_mul_bin_start;
//#endif 		// WITH_POSIT

/*
extern const element_t* _binary_relu_conv10_in_bin_start;
extern const element_t* _binary_relu_conv10_in_bin_end;
extern const int _binary_relu_conv10_in_bin_size;

extern const element_t* _binary_res5b_relu_in_bin_start;
extern const element_t* _binary_res5b_relu_in_bin_end;
extern const int _binary_res5b_relu_in_bin_size;

extern const element_t* _binary_fc1000_in_bin_start;
extern const element_t* _binary_fc1000_in_bin_end;
extern const int _binary_fc1000_in_bin_size;

extern const element_t* _binary_fc1000_out_bin_start;
extern const element_t* _binary_fc1000_out_bin_end;
extern const int _binary_fc1000_out_bin_size;

extern const element_t* _binary_fc1000_weight_bin_start;
extern const element_t* _binary_fc1000_weight_bin_end;
extern const int _binary_fc1000_weight_bin_size;

extern const element_t* _binary_fc1000_bias_mul_bin_start;
extern const element_t* _binary_fc1000_bias_mul_bin_end;
extern const int _binary_fc1000_bias_mul_bin_size;

extern const element_t* _binary_fc1000_blob_bin_start;
extern const element_t* _binary_fc1000_blob_bin_end;
extern const int _binary_fc1000_blob_bin_size;

extern const element_t* _binary_sum_multiplier_bin_start;
extern const element_t* _binary_sum_multiplier_bin_end;
extern const int _binary_sum_multiplier_bin_size;
*/

// taken from CBLAS for compatibility
enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102 };

enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113, AtlasConj=114};


// util methods

element_t min(element_t a, element_t b);

element_t max(element_t a, element_t b);

int min_int(int a, int b);

int max_int(int a, int b);

void qsort(void *aa, size_t n, size_t es, int (*cmp)(const void *, const void *));

int comp_elements(const void* a, const void* b);

void print_top_5(const element_t* a);

void split(element_t* src, element_t* dst, size_t count);

void forward(element_t** prob, size_t* size);

void write_data(element_t* src, size_t count, char* file);

void im2col_cpu(const element_t* data_im,
		const int channels,
		const int height,
		const int width,
		const int kernel_h,
		const int kernel_w,
		const int pad_h,
		const int pad_w,
		const int stride_h,
		const int stride_w,
		const int dilation_h,
		const int dilation_w,
		element_t* data_col);


// caffe methods

void caffe_cpu_gemm(const int TransA,
		const size_t TransB,
		const size_t M,
		const size_t N,
		const size_t K,
		const element_t alpha,
		const element_t* A,
		const element_t* B,
		const element_t beta,
		element_t* C);

void caffe_cpu_gemv(const int TransA,
		const size_t M,
		const size_t N,
		const element_t alpha,
		const element_t* A,
		const element_t* x,
		const element_t beta,
		element_t* y);

void caffe_cpu_scale(const size_t n, const element_t alpha, const element_t* X, element_t* Y);

void caffe_mul(const size_t N, const element_t* X, element_t* Y);

void caffe_div(const size_t N, const element_t* A, const element_t* B, element_t* Y);

void caffe_exp(const size_t N, const element_t* A, element_t* Y);

void caffe_axpy(const size_t N, const element_t alpha, const element_t* X, element_t* Y);

void caffe_cpu_relu(element_t* bottom, element_t* top, size_t count, element_t negative_slope);

void caffe_eltwise(element_t* top_data,
		element_t* bottom1,
		element_t* bottom2,
		element_t coeff1,
		element_t coeff2,
		const size_t count,
		int op);

void caffe_set(const size_t count, element_t val, element_t* a);

void caffe_copy(const size_t count, element_t* src, element_t* dst);


#endif // _COMMON_H_
