/*
 * zero_minor_map.cu
 *
 * Reads input from  input.txt
 * Writes output to  output.txt
 *
 * Input format (input.txt):
 * -------------------------
 *   N P MINOR_K MAX_POWER
 *   a[0][0] a[0][1] ... a[0][N-1]
 *   a[1][0] ...
 *   ...
 *   a[N-1][0] ... a[N-1][N-1]
 *
 * Example input.txt:
 *   3 17 2 16
 *   1  3 12
 *   2  6  8
 *   4  5  7
 *
 * Compile:  nvcc -O2 -o zero_minor_map zero_minor_map.cu
 * Run:      ./zero_minor_map
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime.h>

typedef long long ll;

// ─────────────────────────────────────────────────────────────
//  HOST MATH
// ─────────────────────────────────────────────────────────────

ll mod_pow(ll base, ll exp, ll mod) {
    ll result = 1;
    base %= mod;
    if (base < 0) base += mod;
    while (exp > 0) {
        if (exp & 1) result = result * base % mod;
        base = base * base % mod;
        exp >>= 1;
    }
    return result;
}

ll mod_inv(ll a, ll p) {
    return mod_pow(((a % p) + p) % p, p - 2, p);
}

void mat_mul_mod(ll *C, const ll *A, const ll *B, int n, ll p) {
    ll *tmp = (ll*)calloc(n * n, sizeof(ll));
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++)
            for (int j = 0; j < n; j++)
                tmp[i*n+j] = (tmp[i*n+j] + A[i*n+k] * B[k*n+j]) % p;
    memcpy(C, tmp, n * n * sizeof(ll));
    free(tmp);
}

void mat_pow_mod(ll *out, const ll *M, int exp, int n, ll p) {
    memset(out, 0, n * n * sizeof(ll));
    for (int i = 0; i < n; i++) out[i*n+i] = 1;
    ll *base = (ll*)malloc(n * n * sizeof(ll));
    memcpy(base, M, n * n * sizeof(ll));
    while (exp > 0) {
        if (exp & 1) mat_mul_mod(out, out, base, n, p);
        mat_mul_mod(base, base, base, n, p);
        exp >>= 1;
    }
    free(base);
}

int mat_inv_mod(ll *inv, const ll *M, int n, ll p) {
    ll *aug = (ll*)malloc(n * 2*n * sizeof(ll));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            aug[i*2*n+j]   = ((M[i*n+j] % p) + p) % p;
        for (int j = 0; j < n; j++)
            aug[i*2*n+n+j] = (i == j) ? 1 : 0;
    }
    for (int col = 0; col < n; col++) {
        int pivot = -1;
        for (int row = col; row < n; row++)
            if (aug[row*2*n+col] != 0) { pivot = row; break; }
        if (pivot == -1) { free(aug); return 0; }
        for (int j = 0; j < 2*n; j++) {
            ll t = aug[col*2*n+j];
            aug[col*2*n+j]   = aug[pivot*2*n+j];
            aug[pivot*2*n+j] = t;
        }
        ll inv_piv = mod_inv(aug[col*2*n+col], p);
        for (int j = 0; j < 2*n; j++)
            aug[col*2*n+j] = aug[col*2*n+j] * inv_piv % p;
        for (int row = 0; row < n; row++) {
            if (row == col) continue;
            ll factor = aug[row*2*n+col];
            for (int j = 0; j < 2*n; j++)
                aug[row*2*n+j] = ((aug[row*2*n+j] - factor * aug[col*2*n+j]) % p + p) % p;
        }
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv[i*n+j] = aug[i*2*n+n+j];
    free(aug);
    return 1;
}

// ─────────────────────────────────────────────────────────────
//  COMBINATION GENERATOR
// ─────────────────────────────────────────────────────────────

int combinations(int n, int k, int **out) {
    int count = 1;
    for (int i = 0; i < k; i++) count = count * (n - i) / (i + 1);
    *out = (int*)malloc(count * k * sizeof(int));
    int *combo = (int*)malloc(k * sizeof(int));
    for (int i = 0; i < k; i++) combo[i] = i;
    int idx = 0;
    while (1) {
        memcpy((*out) + idx * k, combo, k * sizeof(int));
        idx++;
        int i = k - 1;
        while (i >= 0 && combo[i] == n - k + i) i--;
        if (i < 0) break;
        combo[i]++;
        for (int j = i+1; j < k; j++) combo[j] = combo[j-1] + 1;
    }
    free(combo);
    return count;
}

// ─────────────────────────────────────────────────────────────
//  DEVICE FUNCTIONS
// ─────────────────────────────────────────────────────────────

__device__ long long dev_mod_pow(long long base, long long exp, long long mod) {
    long long result = 1;
    base %= mod;
    if (base < 0) base += mod;
    while (exp > 0) {
        if (exp & 1) result = result * base % mod;
        base = base * base % mod;
        exp >>= 1;
    }
    return result;
}

__device__ long long dev_mod_inv(long long a, long long p) {
    return dev_mod_pow(((a % p) + p) % p, p - 2, p);
}

__device__ long long submatrix_det_mod(
    const long long *mat, int n,
    const int *row_sel, const int *col_sel, int k,
    long long p
) {
    long long sub[8*8];
    for (int i = 0; i < k; i++)
        for (int j = 0; j < k; j++)
            sub[i*k+j] = mat[row_sel[i]*n + col_sel[j]];

    long long sign = 1;
    for (int col = 0; col < k; col++) {
        int pivot = -1;
        for (int row = col; row < k; row++)
            if (sub[row*k+col] != 0) { pivot = row; break; }
        if (pivot == -1) return 0;
        if (pivot != col) {
            for (int j = 0; j < k; j++) {
                long long t  = sub[col*k+j];
                sub[col*k+j]   = sub[pivot*k+j];
                sub[pivot*k+j] = t;
            }
            sign = (p - sign) % p;
        }
        long long inv_piv = dev_mod_inv(sub[col*k+col], p);
        for (int row = col+1; row < k; row++) {
            long long factor = sub[row*k+col] * inv_piv % p;
            for (int j = col; j < k; j++)
                sub[row*k+j] = ((sub[row*k+j] - factor * sub[col*k+j]) % p + p) % p;
        }
    }
    long long det = sign;
    for (int i = 0; i < k; i++) det = det * sub[i*k+i] % p;
    return det;
}

// ─────────────────────────────────────────────────────────────
//  CUDA KERNEL
// ─────────────────────────────────────────────────────────────

__global__ void check_minors_kernel(
    const long long *mat,
    int n,
    const int *row_combos,
    const int *col_combos,
    int num_rc,
    int num_cc,
    int k,
    long long p,
    int *result
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_rc * num_cc) return;

    int ri = tid / num_cc;
    int ci = tid % num_cc;

    long long det = submatrix_det_mod(
        mat, n,
        row_combos + ri * k,
        col_combos + ci * k,
        k, p
    );
    result[tid] = (det == 0) ? 1 : 0;
}

// ─────────────────────────────────────────────────────────────
//  PRINT MATRIX  (to file)
// ─────────────────────────────────────────────────────────────

void fprint_mat(FILE *f, const char *label, const ll *M, int n) {
    fprintf(f, "%s\n", label);
    for (int i = 0; i < n; i++) {
        fprintf(f, "  [ ");
        for (int j = 0; j < n; j++) fprintf(f, "%4lld ", M[i*n+j]);
        fprintf(f, "]\n");
    }
}

// ─────────────────────────────────────────────────────────────
//  MAIN
// ─────────────────────────────────────────────────────────────

int main(void) {

    // ── Open input file ───────────────────────────────────────
    FILE *fin = fopen("input.txt", "r");
    if (!fin) {
        fprintf(stderr, "ERROR: Cannot open input.txt\n");
        fprintf(stderr, "Expected format:\n");
        fprintf(stderr, "  Line 1: N P MINOR_K MAX_POWER\n");
        fprintf(stderr, "  Next N lines: matrix rows\n");
        return 1;
    }

    int N, P, MINOR_K, MAX_POWER;
    if (fscanf(fin, "%d %d %d %d", &N, &P, &MINOR_K, &MAX_POWER) != 4) {
        fprintf(stderr, "ERROR: Could not read N P MINOR_K MAX_POWER from input.txt\n");
        fclose(fin);
        return 1;
    }

    ll *A = (ll*)malloc(N * N * sizeof(ll));
    for (int i = 0; i < N * N; i++) {
        if (fscanf(fin, "%lld", &A[i]) != 1) {
            fprintf(stderr, "ERROR: Not enough matrix entries in input.txt\n");
            fclose(fin); free(A); return 1;
        }
    }
    fclose(fin);

    // ── Open output file ──────────────────────────────────────
    FILE *fout = fopen("output1_19.txt", "w");
    if (!fout) {
        fprintf(stderr, "ERROR: Cannot open output1_19.txt for writing\n");
        free(A); return 1;
    }

    // ── Echo config ───────────────────────────────────────────
    fprintf(fout, "=== Zero Minor Map ===\n");
    fprintf(fout, "n=%d  p=%d  minor_size=%dx%d  max_power=%d\n\n",
            N, P, MINOR_K, MINOR_K, MAX_POWER);

    // Validate MINOR_K
    if (MINOR_K > N) {
        fprintf(fout, "ERROR: MINOR_K=%d > N=%d\n", MINOR_K, N);
        fclose(fout); free(A); return 1;
    }
    if (MINOR_K > 8) {
        fprintf(fout, "ERROR: MINOR_K=%d exceeds device max of 8\n", MINOR_K);
        fclose(fout); free(A); return 1;
    }

    // ── Print input matrix ────────────────────────────────────
    fprint_mat(fout, "Input Matrix A:", A, N);

    // ── Validate row sums ─────────────────────────────────────
    fprintf(fout, "\nRow Sum Validation (expected p-1 = %d):\n", P-1);
    int row_sum_ok = 1;
    for (int i = 0; i < N; i++) {
        ll s = 0;
        for (int j = 0; j < N; j++) s += A[i*N+j];
        ll sm = ((s % P) + P) % P;
        fprintf(fout, "  Row %d: sum = %lld mod %d = %lld  %s\n",
                i+1, s, P, sm, (sm == P-1) ? "OK" : "WARNING: not p-1");
        if (sm != P-1) row_sum_ok = 0;
    }
    if (!row_sum_ok)
        fprintf(fout, "  WARNING: This may not be a valid pattern matrix.\n");

    // ── Compute inverse ───────────────────────────────────────
    ll *B = (ll*)malloc(N * N * sizeof(ll));
    if (!mat_inv_mod(B, A, N, P)) {
        fprintf(fout, "\nERROR: Matrix A is singular mod %d — cannot invert.\n", P);
        fclose(fout); free(A); free(B); return 1;
    }
    fprint_mat(fout, "\nB = A^{-1} mod p:", B, N);

    // ── Generate combinations ─────────────────────────────────
    int *row_combos, *col_combos;
    int num_rc = combinations(N, MINOR_K, &row_combos);
    int num_cc = combinations(N, MINOR_K, &col_combos);
    int num_minors = num_rc * num_cc;

    fprintf(fout, "\nC(%d,%d) = %d row combos x %d col combos = %d total %dx%d minors\n\n",
            N, MINOR_K, num_rc, num_cc, num_minors, MINOR_K, MINOR_K);

    // ── GPU allocations ───────────────────────────────────────
    int       *d_row_combos, *d_col_combos, *d_result;
    long long *d_mat;
    cudaMalloc(&d_row_combos, num_rc * MINOR_K * sizeof(int));
    cudaMalloc(&d_col_combos, num_cc * MINOR_K * sizeof(int));
    cudaMalloc(&d_result,     num_minors        * sizeof(int));
    cudaMalloc(&d_mat,        N * N             * sizeof(long long));

    cudaMemcpy(d_row_combos, row_combos,
               num_rc * MINOR_K * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_col_combos, col_combos,
               num_cc * MINOR_K * sizeof(int), cudaMemcpyHostToDevice);

    // Allocate heat map: index [k * num_minors + m]
    int *zero_map = (int*)calloc((MAX_POWER+1) * num_minors, sizeof(int));

    ll  *Bk       = (ll*)malloc(N * N * sizeof(ll));
    int *h_result = (int*)malloc(num_minors * sizeof(int));

    int threads = 256;
    int blocks  = (num_minors + threads - 1) / threads;

    // ── Per-power results ─────────────────────────────────────
    fprintf(fout, "%-8s  %-10s  %s\n",
            "Power k", "Zero count", "Zero minor positions (1-indexed)");
    fprintf(fout, "%-8s  %-10s  %s\n",
            "-------", "----------", "--------------------------------");

    for (int k = 1; k <= MAX_POWER; k += 2) {

        mat_pow_mod(Bk, B, k, N, P);
        cudaMemcpy(d_mat, Bk, N*N*sizeof(long long), cudaMemcpyHostToDevice);

        check_minors_kernel<<<blocks, threads>>>(
            d_mat, N,
            d_row_combos, d_col_combos,
            num_rc, num_cc,
            MINOR_K, P, d_result
        );
        cudaDeviceSynchronize();

        cudaMemcpy(h_result, d_result,
                   num_minors*sizeof(int), cudaMemcpyDeviceToHost);

        int count = 0;
        for (int m = 0; m < num_minors; m++) {
            zero_map[k * num_minors + m] = h_result[m];
            if (h_result[m]) count++;
        }

        fprintf(fout, "k = %-3d  |  %4d / %-4d  |  ", k, count, num_minors);
        int first = 1;
        for (int m = 0; m < num_minors; m++) {
            if (!h_result[m]) continue;
            if (!first) fprintf(fout, ", ");
            int ri = m / num_cc, ci = m % num_cc;
            fprintf(fout, "r(");
            for (int x = 0; x < MINOR_K; x++)
                fprintf(fout, "%d%s", row_combos[ri*MINOR_K+x]+1,
                        x < MINOR_K-1 ? "," : "");
            fprintf(fout, ")c(");
            for (int x = 0; x < MINOR_K; x++)
                fprintf(fout, "%d%s", col_combos[ci*MINOR_K+x]+1,
                        x < MINOR_K-1 ? "," : "");
            fprintf(fout, ")");
            first = 0;
        }
        if (count == 0) fprintf(fout, "(none)");
        fprintf(fout, "\n");

        // Print the matrix for this power
        fprintf(fout, "         B^%d =\n", k);
        for (int i = 0; i < N; i++) {
            fprintf(fout, "           [ ");
            for (int j = 0; j < N; j++) fprintf(fout, "%4lld ", Bk[i*N+j]);
            fprintf(fout, "]\n");
        }
        fprintf(fout, "\n");
    }

    // ── Heat map ──────────────────────────────────────────────
    fprintf(fout, "\n=== HEAT MAP (1=zero minor, 0=nonzero) ===\n");
    fprintf(fout, "Rows = odd powers k=1,3,...  Cols = minor index 1..%d\n\n", num_minors);

    // header
    fprintf(fout, "        |");
    for (int m = 0; m < num_minors; m++) fprintf(fout, " %3d |", m+1);
    fprintf(fout, "\n--------|");
    for (int m = 0; m < num_minors; m++) fprintf(fout, "-----|");
    fprintf(fout, "\n");
    for (int k = 1; k <= MAX_POWER; k += 2) {
        fprintf(fout, " k = %2d |", k);
        for (int m = 0; m < num_minors; m++)
            fprintf(fout, "  %d  |", zero_map[k * num_minors + m]);
        fprintf(fout, "\n");
    }

    // Minor legend
    fprintf(fout, "\nMinor Index Legend:\n");
    for (int m = 0; m < num_minors; m++) {
        int ri = m / num_cc, ci = m % num_cc;
        fprintf(fout, "  M%-3d  rows(", m+1);
        for (int x = 0; x < MINOR_K; x++)
            fprintf(fout, "%d%s", row_combos[ri*MINOR_K+x]+1,
                    x < MINOR_K-1 ? "," : "");
        fprintf(fout, ") cols(");
        for (int x = 0; x < MINOR_K; x++)
            fprintf(fout, "%d%s", col_combos[ci*MINOR_K+x]+1,
                    x < MINOR_K-1 ? "," : "");
        fprintf(fout, ")\n");
    }

    // ── Cycle detection ───────────────────────────────────────
    fprintf(fout, "\n=== CYCLE DETECTION ===\n");
    int cycle_found = 0;
    for (int k1 = 1; k1 <= MAX_POWER; k1 += 2)
        for (int k2 = k1+2; k2 <= MAX_POWER; k2 += 2) {
            int match = 1;
            for (int m = 0; m < num_minors; m++)
                if (zero_map[k1*num_minors+m] != zero_map[k2*num_minors+m])
                    { match = 0; break; }
            if (match) {
                fprintf(fout, "  k=%d matches k=%d  =>  period candidate = %d\n",
                        k1, k2, k2-k1);
                cycle_found = 1;
            }
        }
    if (!cycle_found)
        fprintf(fout, "  No exact repeat found in k=1..%d\n", MAX_POWER);

    // ── Coverage check ────────────────────────────────────────
    fprintf(fout, "\n=== COVERAGE CHECK ===\n");
    int never_zero_count = 0;
    for (int m = 0; m < num_minors; m++) {
        int ever = 0;
        for (int k = 1; k <= MAX_POWER; k += 2)
            if (zero_map[k*num_minors+m]) { ever = 1; break; }
        int ri = m / num_cc, ci = m % num_cc;
        fprintf(fout, "  M%-3d r(", m+1);
        for (int x = 0; x < MINOR_K; x++)
            fprintf(fout, "%d%s", row_combos[ri*MINOR_K+x]+1,
                    x < MINOR_K-1 ? "," : "");
        fprintf(fout, ")c(");
        for (int x = 0; x < MINOR_K; x++)
            fprintf(fout, "%d%s", col_combos[ci*MINOR_K+x]+1,
                    x < MINOR_K-1 ? "," : "");
        fprintf(fout, "): %s\n", ever ? "zero at some k" : "NEVER ZERO");
        if (!ever) never_zero_count++;
    }
    fprintf(fout, "\nSummary: %d / %d minors become zero across k=1..%d\n",
            num_minors - never_zero_count, num_minors, MAX_POWER);

    fprintf(fout, "\nDone.\n");
    fclose(fout);

    // ── Cleanup ───────────────────────────────────────────────
    free(A); free(B); free(Bk); free(h_result);
    free(row_combos); free(col_combos); free(zero_map);
    cudaFree(d_row_combos); cudaFree(d_col_combos);
    cudaFree(d_result); cudaFree(d_mat);

    printf("Done. Results written to output.txt\n");
    return 0;
}