#include <clocale>
#include <cstdio>

// version 1
void matmul_naive(double *C, double *A, double *B, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double sum = 0.0;
            for (int k = 0; k < N; k++) {
                sum += A[i * N + k] * B[k * N + j];
            }
            C[i * N + j] = sum;
        }
    }
}

// version 2
void matmul_ikj(double *C, double *A, double *B, int N) {
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < N; k++) {
            double aik = A[i * N + k];  // use multi times
            for (int j = 0; j < N; j++) {
                C[i * N + j] += aik * B[k * N + j];  // row priority access for B
            }
        }
    }
}
// version 3
void matmul_blocked(double *C, double *A, double *B, int N) {
    const int BLOCK_SIZE = 64;  // choose the appropriate size

    for (int i_blk = 0; i_blk < N; i_blk += BLOCK_SIZE) {
        for (int j_blk = 0; j_blk < N; j_blk += BLOCK_SIZE) {
            for (int k_blk = 0; k_blk < N; k_blk += BLOCK_SIZE) {
                // BLOCK_SIZE x BLOCK_SIZE
                int i_end = min(i_blk + BLOCK_SIZE, N);
                int j_end = min(j_blk + BLOCK_SIZE, N);
                int k_end = min(k_blk + BLOCK_SIZE, N);

                for (int i = i_blk; i < i_end; i++) {
                    for (int k = k_blk; k < k_end; k++) {
                        double aik = A[i * N + k];
                        for (int j = j_blk; j < j_end; j++) {
                            C[i * N + j] += aik * B[k * N + j];
                        }
                    }
                }
            }
        }
    }
}


// version 4 Allocating large page memory on Linux
double* allocate_huge_pages(size_t size) {
    // use mmap and MAP_HUGETLB tag
    void* ptr = mmap(NULL, size, PROT_READ | PROT_WRITE,
                     MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB,
                     -1, 0);
    return (double*)ptr;
}

int main() {
    size_t matrix_size = N * N * sizeof(double);
    double* A = allocate_huge_pages(matrix_size);
    double* B = allocate_huge_pages(matrix_size);
    double* C = allocate_huge_pages(matrix_size);

    matmul_blocked(C, A, B, N);
}
