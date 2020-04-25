// THIS IS A GENERATED FILE!
#include "common2.h"

void forward(element_t** prob, size_t* size) {
	element_t* input_ptr = (element_t*)(&_binary_relu3_in_bin_start);
	element_t* output_ptr = &data[0];
	element_t* ip1_weights_ptr = (element_t*)(&_binary_ip1_weights_bin_start);
	element_t* ip1_bias_mul_ptr = (element_t*)(&_binary_ip1_bias_mul_bin_start);
	element_t* ip1_blob_ptr = (element_t*)(&_binary_ip1_blob_bin_start);
	element_t* prob_sum_multiplier_ptr = (element_t*)(&_binary_prob_sum_multiplier_bin_start);
	element_t* aux_ptr;
	element_t* top_ptr;
	element_t* bottom_ptr;

	// relu3
	caffe_cpu_relu(input_ptr, output_ptr, 4096, zero);
#ifdef WDEBUG
	write_data(output_ptr, 4096, "relu3_out.bin");
#endif
	aux_ptr = input_ptr; input_ptr = output_ptr; output_ptr = aux_ptr;

	// pool3
#ifdef WDEBUG
	write_data(input_ptr, 4096, "pool3_in.bin");
#endif
	caffe_set(1024, zero, output_ptr);
	bottom_ptr = input_ptr; top_ptr = output_ptr;
	for (int n = 0; n < 1; ++n) {
		for (int c = 0; c < 64; ++c) {
			for (int ph = 0; ph < 4; ++ph) {
				for (int pw = 0; pw < 4; ++pw) {
					int hstart = ph * 2 - 0;
					int wstart = pw * 2 - 0;
					int hend = min_int(hstart + 3, 8);
					int wend = min_int(wstart + 3, 8);
					int pool_size = (hend - hstart) * (wend - wstart);
					hstart = max_int(hstart, 0);
					wstart = max_int(wstart, 0);
					hend = min_int(hend, 8);
					wend = min_int(wend, 8);
					for (int h = hstart; h < hend; ++h) {
						for (int w = wstart; w < wend; ++w) {
							top_ptr[ph * 4 + pw] += bottom_ptr[h * 8 + w];
						}
					}
					top_ptr[ph * 4 + pw] /= pool_size;
				}
			}
			bottom_ptr += 64;
			top_ptr += 16;
		}
	}
#ifdef WDEBUG
	write_data(output_ptr, 1024, "pool3_out.bin");
#endif
	aux_ptr = input_ptr; input_ptr = output_ptr; output_ptr = aux_ptr;

	// ip1
#ifdef WDEBUG
	write_data(input_ptr, 1024, "ip1_in.bin");
#endif
	caffe_cpu_gemm(CblasNoTrans, CblasTrans, 1, 10, 1024, one, input_ptr, ip1_weights_ptr, zero, output_ptr);
#ifdef WDEBUG
	write_data(output_ptr, 512, "ip1_out_tmp.bin");
#endif
	caffe_cpu_gemm(CblasNoTrans, CblasNoTrans, 1, 10, 1, one, ip1_bias_mul_ptr, ip1_blob_ptr, one, output_ptr);
#ifdef WDEBUG
	write_data(output_ptr, 10, "ip1_out.bin");
#endif
	aux_ptr = input_ptr; input_ptr = output_ptr; output_ptr = aux_ptr;

	// prob
	element_t scale_ptr[1];
	caffe_copy(10, input_ptr, output_ptr);
	for (int i = 0; i < 1; ++i) {
		caffe_copy(1, input_ptr + i * 10, scale_ptr);
		for (int j = 0; j < 10; j++) {
			for (int k = 0; k < 1; k++) {
				scale_ptr[k] = max(scale_ptr[k], input_ptr[i * 10 + j * 1 + k]);
			}
		}
		caffe_cpu_gemm(CblasNoTrans, CblasNoTrans, 10, 1, 1, minus_one, prob_sum_multiplier_ptr, scale_ptr, one, output_ptr);
		caffe_exp(10, output_ptr, output_ptr);
		caffe_cpu_gemv(CblasTrans, 10, 1, one, output_ptr, prob_sum_multiplier_ptr, zero, scale_ptr);
		aux_ptr = output_ptr;
		for (int j = 0; j < 10; j++) {
			caffe_div(1, aux_ptr, scale_ptr, aux_ptr);
			aux_ptr += 1;
		}
	}
#ifdef WDEBUG
	write_data(output_ptr, 10, "prob_out.bin");
#endif
	aux_ptr = input_ptr; input_ptr = output_ptr; output_ptr = aux_ptr;

	*prob = input_ptr;
	*size = 10;
}
