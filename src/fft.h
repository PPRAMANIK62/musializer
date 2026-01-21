/**
 * =============================================================================
 * DFT & FFT - EDUCATIONAL IMPLEMENTATION
 * =============================================================================
 *
 * This file contains two implementations of the Fourier Transform:
 *   1. dft() - Naive O(N²) Discrete Fourier Transform
 *   2. fft() - Cooley-Tukey O(N log N) Fast Fourier Transform
 *
 * BACKGROUND THEORY:
 * ------------------
 * The Fourier Transform is one of the most important mathematical tools in
 * signal processing. It allows us to decompose a complex signal into its
 * constituent frequencies - essentially answering the question:
 * "What frequencies are present in this signal, and how strong are they?"
 *
 * ANALOGY:
 * Imagine you hear a chord played on a piano. Your ear perceives a single
 * sound, but it's actually made up of multiple notes (frequencies) played
 * together. The Fourier Transform is like having perfect pitch - it can
 * tell you exactly which notes are being played and how loud each one is.
 *
 * THE MATH:
 * ---------
 * Any periodic signal can be represented as a sum of sinusoids (sines and
 * cosines) at different frequencies. This is called Fourier's theorem:
 *
 *   signal(t) = Σ [A_f * sin(2π * f * t) + B_f * cos(2π * f * t)]
 *
 * Where:
 *   - f = frequency (how many cycles per unit time)
 *   - A_f, B_f = amplitudes (how strong each frequency component is)
 *   - t = time
 *
 * The DFT formula for finding frequency component 'k' is:
 *
 *   X[k] = Σ x[n] * e^(-i * 2π * k * n / N)    for n = 0 to N-1
 *
 * Using Euler's formula (e^(iθ) = cos(θ) + i*sin(θ)), this becomes:
 *
 *   X[k] = Σ x[n] * [cos(2πkn/N) - i*sin(2πkn/N)]
 *
 * This code implements the full complex DFT/FFT using C99's complex.h library,
 * applying Euler's formula directly via cexp(). The output contains both
 * real (cosine correlation) and imaginary (sine correlation) components,
 * allowing us to extract both magnitude and phase information.
 *
 * HOW CORRELATION WORKS:
 * ----------------------
 * For each target frequency f, we:
 * 1. Generate a "probe" wave (complex exponential) at that frequency
 * 2. Multiply our input signal by this probe (point by point)
 * 3. Sum up all the products
 *
 * If the input signal contains frequency f:
 *   - The probe wave and signal will be "in sync"
 *   - Positive parts multiply with positive parts, negative with negative
 *   - The sum will be large (positive or negative)
 *
 * If the input signal does NOT contain frequency f:
 *   - The waves will be out of sync
 *   - Products will sometimes be positive, sometimes negative
 *   - They cancel out, sum approaches zero
 *
 * This is called "correlation" - measuring how similar our signal is to
 * each possible frequency component.
 *
 * Mathematical basis: This is the orthogonality property of sinusoids.
 * Sine waves at different frequencies are "orthogonal" - their inner
 * product (sum of element-wise multiplication) is zero.
 *
 * =============================================================================
 */

#ifndef FFT_H
#define FFT_H

#include <assert.h>
#include <math.h>
#include <raylib.h>  /* For the PI constant - raylib defines PI as 3.14159... */
#include <stdio.h>   /* For printf() to display results */
#include <complex.h> /* For complex number support: float complex, cexp(), I, creal(), cimag() */

/**
 * =============================================================================
 * DFT - Discrete Fourier Transform (Naive Implementation)
 * =============================================================================
 *
 * COMPLEXITY: O(N²)
 * For each of N frequency bins, we iterate through all N samples.
 * This is simple but slow for large N.
 *
 * THE CORE OPERATION (using Euler's Formula):
 * -------------------------------------------
 * Euler's formula: e^(iθ) = cos(θ) + i*sin(θ)
 *
 * cexp(2 * I * PI * f * t) computes e^(i*2π*f*t) which equals:
 *   cos(2π*f*t) + i*sin(2π*f*t)
 *
 * By multiplying in[i] by this complex exponential:
 *   in[i] * e^(i*2πft) = in[i]*cos(2πft) + i*in[i]*sin(2πft)
 *
 * Summing over all samples:
 *   - Real part accumulates: Σ in[i]*cos(2πft) → cosine correlation
 *   - Imag part accumulates: Σ in[i]*sin(2πft) → sine correlation
 *
 * @param in   Input signal (real-valued time-domain samples)
 * @param out  Output spectrum (complex frequency-domain coefficients)
 * @param n    Number of samples
 */
static inline void dft(float in[], float complex out[], size_t n) {
    for (size_t f = 0; f < n; f++) {
        out[f] = 0;
        for (size_t i = 0; i < n; i++) {
            float t = (float)i / n;
            out[f] += in[i] * cexp(2 * I * PI * f * t);
        }
    }
}

/**
 * =============================================================================
 * FFT - Fast Fourier Transform (Cooley-Tukey Algorithm)
 * =============================================================================
 *
 * COMPLEXITY: O(N log N)
 * Recursively divides the problem in half, achieving massive speedup.
 *
 * THE KEY INSIGHT (Danielson-Lanczos Lemma):
 * ------------------------------------------
 * The DFT of size N can be rewritten as the sum of two DFTs of size N/2:
 *
 *   X[k] = Σ x[n] * W_N^(kn)     where W_N = e^(-i*2π/N)
 *        = Σ x[2m] * W_N^(k*2m) + Σ x[2m+1] * W_N^(k*(2m+1))
 *        = Σ x_even[m] * W_(N/2)^(km) + W_N^k * Σ x_odd[m] * W_(N/2)^(km)
 *        = E[k] + W_N^k * O[k]
 *
 * Where:
 *   E[k] = DFT of even-indexed samples: x[0], x[2], x[4], ...
 *   O[k] = DFT of odd-indexed samples:  x[1], x[3], x[5], ...
 *   W_N^k = e^(-i*2π*k/N) = "twiddle factor"
 *
 * THE BUTTERFLY OPERATION:
 * ------------------------
 * Due to periodicity of W_N^k, we can compute two outputs at once:
 *
 *   X[k]     = E[k] + W_N^k * O[k]
 *   X[k+N/2] = E[k] - W_N^k * O[k]    (because W_N^(k+N/2) = -W_N^k)
 *
 * This is called a "butterfly" because of how the data flow diagram looks:
 *
 *       E[k] ----+----> X[k] = E[k] + v
 *                 \  /
 *                  \/
 *                  /\
 *                 /  \
 *   W*O[k] = v --+----> X[k+N/2] = E[k] - v
 *
 * WHY IT'S FASTER:
 * ----------------
 * DFT:  N² operations for N samples
 *       N=1024 → 1,048,576 operations
 *
 * FFT:  N log₂(N) operations
 *       N=1024 → 10,240 operations (100x faster!)
 *
 * The recursion:
 *   T(N) = 2*T(N/2) + N    →    T(N) = N log₂(N)
 *
 * REQUIREMENT: N must be a power of 2 (for this implementation)
 *
 * THE STRIDE TRICK:
 * -----------------
 * Instead of copying data into new arrays for even/odd samples,
 * we use 'stride' to skip elements:
 *   - stride=1: access every element (full array)
 *   - stride=2: access every other element (even or odd)
 *   - stride=4: access every 4th element
 *   - etc.
 *
 * This avoids memory allocation and copying overhead.
 *
 * @param in      Input signal (real-valued samples)
 * @param stride  Step size for accessing samples (enables in-place recursion)
 * @param out     Output spectrum (complex coefficients)
 * @param n       Number of samples to process at this recursion level
 */
static inline void fft(float in[], size_t stride, float complex out[], size_t n) {
    assert(n > 0);
    
    /**
     * BASE CASE: Single element
     * The DFT of a single sample is just that sample itself.
     * X[0] = x[0] * e^0 = x[0] * 1 = x[0]
     */
    if (n == 1) {
        out[0] = in[0];
        return;
    }
    
    /**
     * RECURSIVE CASE: Divide and Conquer
     * 
     * Split into even and odd indexed samples:
     *   - Even: in[0], in[2], in[4], ... → stored in out[0..n/2-1]
     *   - Odd:  in[1], in[3], in[5], ... → stored in out[n/2..n-1]
     *
     * Double the stride to skip every other element.
     * Add 'stride' to pointer to start at odd indices.
     */
    fft(in, stride*2, out, n/2);                     /* Even indices */
    fft(in + stride, stride*2, out + n/2, n/2);      /* Odd indices */
    
    /**
     * COMBINE: Butterfly Operation
     * 
     * For each k from 0 to N/2-1:
     *   v = W_N^k * O[k]  = e^(-i*2π*k/N) * out[k + n/2]
     *   X[k]     = E[k] + v
     *   X[k+N/2] = E[k] - v
     *
     * The twiddle factor e^(-i*2π*k/N) rotates the odd contribution
     * before combining with the even contribution.
     */
    for (size_t k = 0; k < n/2; k++) {
        float t = (float)k/n;
        float complex v = cexp(-I * 2*PI*t) * out[k + n/2];  /* Twiddle × Odd */
        float complex e = out[k];                             /* Even result */
        out[k]       = e + v;  /* First half:  E[k] + W*O[k] */
        out[k + n/2] = e - v;  /* Second half: E[k] - W*O[k] */
    }
}

/**
 * =============================================================================
 * TEST - FFT Test Function
 * =============================================================================
 * Call this function to verify FFT correctness with a known test signal.
 */
static inline int fft_test(void) {
    /**
     * N = Number of samples in our signal
     *
     * This also determines our frequency resolution:
     * - With N samples, we can detect N different frequency "bins"
     * - Frequency bin f represents the frequency: f * (sample_rate / N)
     *
     * IMPORTANT: For FFT, N must be a power of 2 (2, 4, 8, 16, 32, 64, ...)
     * 
     * In real audio applications, N is typically:
     *   - 256, 512, 1024, 2048, 4096, etc.
     *   - Larger N = better frequency resolution, but more latency
     */
    #define TEST_N 16
    
    /**
     * Input buffer:  Time-domain signal (amplitude vs time)
     * Output buffer: Frequency-domain data (complex amplitude vs frequency)
     *
     * The transformation:
     *   Time Domain (in[]) ──FFT──> Frequency Domain (out[])
     *   "What does the signal look like over time?"
     *                       ──>
     *   "What frequencies make up this signal?"
     */
    float in[TEST_N];
    float complex out[TEST_N];

    /**
     * =========================================================================
     * STEP 1: SYNTHESIZE A TEST SIGNAL (Mixing Frequencies)
     * =========================================================================
     *
     * We create a test signal by adding three waves together:
     *   - cos(2π*1*t): Frequency 1 as COSINE → appears in REAL part
     *   - sin(2π*2*t): Frequency 2 as SINE   → appears in IMAGINARY part
     *   - cos(2π*3*t): Frequency 3 as COSINE → appears in REAL part
     *
     * This is like playing three notes simultaneously on a piano.
     *
     * WHY MIX COS AND SIN?
     * This demonstrates a key property: the DFT/FFT separates phase information.
     * - Cosine components → appear in the REAL part of the output
     * - Sine components   → appear in the IMAGINARY part of the output
     *
     * The formula: signal(t) = cos(2π*1*t) + sin(2π*2*t) + cos(2π*3*t)
     *
     * Where t ranges from 0 to (n-1)/n (normalized time).
     */
    for (size_t i = 0; i < TEST_N; i++) {
        /**
         * t = normalized time, ranging from 0 to (n-1)/n
         *
         * Example with n=16:
         *   i=0:  t = 0/16  = 0.0000
         *   i=1:  t = 1/16  = 0.0625
         *   i=2:  t = 2/16  = 0.1250
         *   ...
         *   i=15: t = 15/16 = 0.9375
         *
         * Note: t never reaches 1.0, which prevents sampling the same
         * point twice in a periodic signal.
         */
        float t = (float)i / TEST_N;
        
        /**
         * The '2*PI*f*t' pattern comes from the general wave formula:
         *   cos(2πft) or sin(2πft) where f is frequency and t is time
         *
         * 2π radians = 360° = one complete cycle
         *
         * - Frequency 1: completes 1 cycle over the sample window
         * - Frequency 2: completes 2 cycles over the sample window
         * - Frequency 3: completes 3 cycles over the sample window
         */
        in[i] = cosf(2 * PI * t * 1) + sinf(2 * PI * t * 2) + cosf(2 * PI * t * 3);
    }

    /**
     * =========================================================================
     * STEP 2: COMPUTE THE FOURIER TRANSFORM
     * =========================================================================
     *
     * You can use either:
     *   dft(in, out, n);     // O(N²) - slow but simple
     *   fft(in, 1, out, n);  // O(N log N) - fast!
     *
     * Both produce the same result (within floating-point precision).
     * For N=16, the difference is negligible.
     * For N=4096, FFT is ~340x faster!
     */
    fft(in, 1, out, TEST_N);

    /**
     * =========================================================================
     * STEP 3: DISPLAY RESULTS
     * =========================================================================
     *
     * Print the frequency spectrum - showing the complex coefficient of each
     * frequency component detected in our signal.
     *
     * Expected output format: "FF: +real +imagi"
     *
     * For our input signal cos(2π*1*t) + sin(2π*2*t) + cos(2π*3*t):
     *
     *   00: +0.00 +0.00i   ← Frequency 0 (DC offset) - not present
     *   01: +8.00 +0.00i   ← Frequency 1 - COSINE (real part has energy!)
     *   02: +0.00 -8.00i   ← Frequency 2 - SINE (imag part has energy!)
     *   03: +8.00 +0.00i   ← Frequency 3 - COSINE (real part has energy!)
     *   04: +0.00 +0.00i   ← Frequency 4 - not present
     *   ...
     *   13: +8.00 +0.00i   ← Mirror of frequency 3 (N-3 = 16-3 = 13)
     *   14: +0.00 +8.00i   ← Mirror of frequency 2 (N-2 = 16-2 = 14)
     *   15: +8.00 +0.00i   ← Mirror of frequency 1 (N-1 = 16-1 = 15)
     *
     * KEY INSIGHTS:
     * 1. Frequency 1 (cosine input) → energy in REAL part
     * 2. Frequency 2 (sine input)   → energy in IMAGINARY part (negative!)
     * 3. Frequency 3 (cosine input) → energy in REAL part
     *
     * NOTE ON SYMMETRY (Nyquist Theorem):
     * For real-valued input signals, the DFT output is symmetric around N/2.
     * Frequencies above N/2 are "aliases" or mirrors of lower frequencies.
     * In practice, we only use bins 0 to N/2.
     *
     * TO GET MAGNITUDE (how strong the frequency is):
     *   magnitude = cabs(out[f]) = sqrt(real² + imag²)
     *
     * TO GET PHASE (timing offset of the frequency):
     *   phase = carg(out[f]) = atan2(imag, real)
     */
    for (size_t f = 0; f < TEST_N; f++) {
        printf("%02zu: %+.2f %+.2fi\n", f, creal(out[f]), cimag(out[f]));
    }

    return 0;
}

/**
 * =============================================================================
 * COMPARISON: DFT vs FFT
 * =============================================================================
 *
 * | Aspect      | DFT              | FFT (Cooley-Tukey)        |
 * |-------------|------------------|---------------------------|
 * | Complexity  | O(N²)            | O(N log N)                |
 * | N=1024      | 1,048,576 ops    | 10,240 ops (100x faster)  |
 * | N=4096      | 16,777,216 ops   | 49,152 ops (341x faster)  |
 * | N=65536     | 4 billion ops    | 1 million ops (4000x!)    |
 * | Constraint  | Any N            | N must be power of 2      |
 * | Code        | Simple loops     | Recursive + butterfly     |
 *
 * For real-time audio (44.1kHz, 2048 samples), FFT is essential.
 *
 * =============================================================================
 * FURTHER READING & ADVANCED TOPICS
 * =============================================================================
 *
 * 1. WINDOWING:
 *    Real-world signals aren't perfectly periodic within our sample window.
 *    "Windowing" (Hamming, Hanning, Blackman) tapers the signal at edges
 *    to reduce spectral leakage - false frequencies appearing in the output.
 *
 * 2. FREQUENCY RESOLUTION:
 *    Resolution = Sample_Rate / N
 *    Example: 44100 Hz / 2048 samples = 21.5 Hz per bin
 *    More samples = finer frequency resolution, but more latency.
 *
 * 3. NYQUIST FREQUENCY:
 *    The highest frequency we can detect is Sample_Rate / 2.
 *    Frequencies above this will "alias" to lower frequencies.
 *    This is why audio is sampled at 44.1kHz - to capture up to ~22kHz.
 *
 * 4. ITERATIVE FFT:
 *    This implementation is recursive. Production FFT libraries use
 *    iterative implementations with bit-reversal permutation for better
 *    cache performance and to avoid recursion overhead.
 *
 * 5. APPLICATIONS IN AUDIO VISUALIZATION:
 *    - Music visualizers use FFT to get frequency spectrum in real-time
 *    - Typical setup: 2048 samples @ 44.1kHz = ~46ms windows, ~21Hz resolution
 *    - Low frequencies (bass)   → left side of visualization
 *    - High frequencies (treble) → right side of visualization
 *    - Magnitude values (cabs) drive bar heights or color intensity
 *
 * =============================================================================
 */

#endif /* FFT_H */
