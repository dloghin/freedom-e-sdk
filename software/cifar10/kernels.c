#include "common2.h"

#include <math.h>

#include <stdio.h>

// utils

element_t max(element_t a, element_t b) {
	if (a > b)
		return a;
	return b;
}

element_t min(element_t a, element_t b) {
	if (a < b)
		return a;
	return b;
}

int max_int(int a, int b) {
	if (a > b)
		return a;
	return b;
}

int min_int(int a, int b) {
	if (a < b)
		return a;
	return b;
}

element_t myexp(element_t x) {
	element_t sum = one;
	element_t fact = one;
	element_t n = one;
	char sign = 0;
	if (x < zero) {
		sign = 1;
		x = minus_one * x;
	}

	for (int i = 1; i < 100; i++) {
		fact = fact * (x / n);
		sum = sum + fact;
		n = n + one;
	}
	if (sign)
		return one / sum;
	return sum;
}

// === not depending on type ===

void caffe_set(const size_t count, element_t val, element_t* a) {
	for (size_t i = 0; i < count; i++)
		a[i] = val;
}

void caffe_copy(const size_t count, element_t* src, element_t* dst) {
	for (size_t i = 0; i < count; i++)
		dst[i] = src[i];
}

void caffe_cpu_relu(element_t* bottom, element_t* top, size_t count, element_t negative_slope) {
	for (int i = 0; i < count; ++i) {
		top[i] = max(bottom[i], zero) + negative_slope * min(bottom[i], zero);
	}
}

int comp_elements(const void* a, const void* b) {
	element_t ea = *(element_t*)a;
	element_t eb = *(element_t*)b;
	if (ea == eb)
		return 0;
	if (ea > eb)
		return -1;
	return 1;
}

#ifdef PDEBUG
void print_top_5(const element_t* a) {
	printf("Top 5 predictions:\n");
	for (int i = 0; i < 5; i++)
		printf("\t%lf\n", a[i]);
}
#endif

void split(element_t* src, element_t* dst, size_t count) {
	for (int i = 0; i < count; i++)
		dst[i] = src[i];
}

enum EltwiseParameter_EltwiseOp {
	EltwiseParameter_EltwiseOp_PROD = 0,
	EltwiseParameter_EltwiseOp_SUM = 1,
	EltwiseParameter_EltwiseOp_MAX = 2
};

/**
 * top_data -> output_ptr
 */
/*
void caffe_eltwise(element_t* top_data,
		element_t* bottom1,
		element_t* bottom2,
		element_t coeff1,
		element_t coeff2,
		const size_t count,
		int op) {

	switch (op) {
	case EltwiseParameter_EltwiseOp_PROD:
		caffe_mul(count, bottom1, bottom2);
		for (int i = 0; i < count; i++)
			top_data[i] = bottom1[i];
		break;
	case EltwiseParameter_EltwiseOp_SUM:
		caffe_set(count, 0.0, top_data);
		caffe_axpy(count, coeff1, bottom1, top_data);
		caffe_axpy(count, coeff2, bottom2, top_data);
		break;
	default:
#ifdef PDEBUG
		printf("Unknown op: %d\n", op);
#endif
		return;
	}
}
*/

void CxMxM(const element_t* A, const element_t* B, element_t* C, size_t M, size_t N, size_t K, element_t alpha, char reset) {
	size_t i, j, k;
	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			if (reset)
				C[i*N+j] = zero;
			element_t sum = zero;
			for (k = 0; k < K; k++)
				sum += A[i*K+k] * B[k*N+j];
			C[i*N+j] += alpha * sum;
		}
	}
}

void CxMTxM(const element_t* A, const element_t* B, element_t* C, size_t M, size_t N, size_t K, element_t alpha, char reset) {
	size_t i, j, k;
	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			if (reset)
				C[i*N+j] = zero;
			element_t sum = zero;
			for (k = 0; k < K; k++)
				sum += A[k*M+i] * B[k*N+j];
			C[i*N+j] += alpha * sum;
		}
	}
}

void CxMxMT(const element_t* A, const element_t* B, element_t* C, size_t M, size_t N, size_t K, element_t alpha, char reset) {
	size_t i, j, k;
	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			if (reset)
				C[i*N+j] = zero;
			element_t sum = zero;
			for (k = 0; k < K; k++)
				sum += A[i*K+k] * B[j*K+k];
			C[i*N+j] += alpha * sum;
		}
	}
}

void CxMTxMT(const element_t* A, const element_t* B, element_t* C, size_t M, size_t N, size_t K, element_t alpha, char reset) {
	size_t i, j, k;
	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			if (reset)
				C[i*N+j] = zero;
			element_t sum = zero;
			for (k = 0; k < K; k++)
				sum += A[k*M+i] * B[j*K+k];
			C[i*N+j] += alpha * sum;
		}
	}
}

void MpM(element_t* A, const element_t* B, size_t M, size_t N) {
	size_t i, j;
	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			A[i*N+j] += B[i*N+j];
		}
	}
}

void CxM(element_t* A, const element_t alpha, size_t M, size_t N) {
	size_t i, j;
	for (i = 0; i < M; i++)
		for (j = 0; j < N; j++)
			A[i*N+j] *= alpha;
}

void caffe_cpu_gemm(const int TransA,
		const size_t TransB,
		const size_t M,
		const size_t N,
		const size_t K,
		const element_t alpha,
		const element_t* A,
		const element_t* B,
		const element_t beta,
		element_t* C) {

	if (TransA == CblasNoTrans) {
		if (TransB == CblasNoTrans) {
			CxM(C, beta, M, N);
			CxMxM(A, B, C, M, N, K, alpha, 0);
		}
		else {
			CxM(C, beta, M, N);
			CxMxMT(A, B, C, M, N, K, alpha, 0);
		}
	}
	else {
		if (TransB == CblasNoTrans) {
			CxM(C, beta, N, K);
			CxMTxM(A, B, C, M, N, K, alpha, 0);
		}
		else {
			CxM(C, beta, N, N);
			CxMTxMT(A, B, C, M, N, K, alpha, 0);
		}
	}
/*
	for (int j = 0; j < 10; j++)
		printf("%f ", A[j]);
	printf("\n");
	for (int j = 0; j < 10; j++)
		printf("%f ", B[j]);
	printf("\n");
	for (int j = 0; j < 10; j++)
		printf("%f ", C[j]);
	printf("\n");
*/
}

void caffe_mul(const size_t N, const element_t* X, element_t* Y) {
	size_t i;
	for (i = 0; i < N; i++)
		Y[i] = X[i] * Y[i];
}

void caffe_axpy(const size_t N, const element_t alpha, const element_t* X, element_t* Y) {
	size_t i;
	for (i = 0; i < N; i++)
		Y[i] = alpha * X[i] * Y[i];
}

void caffe_cpu_gemv(int TransA, const size_t M,
		const size_t N, const element_t alpha, const element_t* A, const element_t* x,
		const element_t beta, element_t* y) {

	size_t i, j;

	size_t sY = (TransA == CblasNoTrans) ? M : N;

	for (i = 0; i < sY; i++)
		y[i] *= beta * y[i];

	if (TransA == CblasNoTrans) {
		for (i = 0; i < M; i++) {
			element_t sum = zero;
			for (j = 0; j < N; j++)
				sum += A[i*N+j] * x[j];
			y[i] += sum * alpha;
		}
	}
	else {
		for (i = 0; i < N; i++) {
			element_t sum = zero;
			for (j = 0; j < M; j++)
				sum += A[i*M+j] * x[j];
			y[i] += sum * alpha;
		}
	}
}

#if 0
void caffe_sqr(const int n, const element_t* a, element_t* y) {
	int i;
	for (i = 0; i < n; i++)
		y[i] = a[i] * a[i];
}

element_t sqrt_asm(element_t x) {
	element_t res;
	asm("fsqrt.s %0,%1\n\t" : "=f" (res) : "f" (x) : "cc");
	return res;
}

void caffe_sqrt(const int n, const element_t* a, element_t* y) {
	int i;
	for (i = 0; i < n; i++)
		y[i] = sqrt(a[i]);
}
#endif

void caffe_exp(const size_t N, const element_t* A, element_t* Y) {
	for(size_t i = 0; i < N; i++)
		Y[i] = myexp(A[i]);
	// Y[i] = exp(A[i]);
}

void caffe_div(const size_t N, const element_t* A, const element_t* B, element_t* Y) {
	for (size_t i = 0; i < N; i++)
		Y[i] = A[i] / B[i];
}
