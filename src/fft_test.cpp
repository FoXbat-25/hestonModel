#include <fftw3.h>
#include <iostream>
#include <cmath>

int main() {
    const int N = 8;

    fftw_complex *in, *out;
    fftw_plan p;

    in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    // Example input: sin wave
    for (int i = 0; i < N; ++i) {
        in[i][0] = sin(2 * M_PI * i / N);  // Real part
        in[i][1] = 0.0;                    // Imaginary part
    }

    // Create plan and execute
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    // Output FFT results
    for (int i = 0; i < N; ++i) {
        std::cout << "out[" << i << "] = " 
                  << out[i][0] << " + " << out[i][1] << "i" << std::endl;
    }

    // Clean up
    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);

    return 0;
}