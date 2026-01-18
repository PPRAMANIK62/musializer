/**
 * =============================================================================
 * DISCRETE FOURIER TRANSFORM (DFT) - EDUCATIONAL IMPLEMENTATION
 * =============================================================================
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
 * The DFT formula for finding frequency component 'f' is:
 *
 *   X[f] = Σ x[n] * e^(-i * 2π * f * n / N)
 *        n=0 to N-1
 *
 * Using Euler's formula (e^(iθ) = cos(θ) + i*sin(θ)), this becomes:
 *
 *   X[f] = Σ x[n] * [cos(2π * f * n/N) - i * sin(2π * f * n/N)]
 *
 * This code implements the full complex DFT using C99's complex.h library,
 * applying Euler's formula directly via cexp(). The output contains both
 * real (cosine correlation) and imaginary (sine correlation) components,
 * allowing us to extract both magnitude and phase information.
 *
 * HOW IT WORKS:
 * -------------
 * For each target frequency f, we:
 * 1. Generate a "probe" sine wave at that frequency
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
 * This process is called "correlation" - we're measuring how similar
 * our signal is to each possible frequency component.
 *
 * COMPLEXITY:
 * -----------
 * This naive DFT implementation is O(N²) - for each of N frequency bins,
 * we iterate through all N samples. The Fast Fourier Transform (FFT)
 * reduces this to O(N log N) using clever mathematical tricks.
 *
 * =============================================================================
 */

#include <math.h>
#include <raylib.h>  /* For the PI constant - raylib defines PI as 3.14159... */
#include <stdio.h>   /* For printf() to display results */
#include <complex.h> /* For complex number support: float complex, cexp(), I, creal(), cimag() */

int main() {
    /**
     * N = Number of samples in our signal
     *
     * This also determines our frequency resolution:
     * - With N samples, we can detect N different frequency "bins"
     * - Frequency bin f represents the frequency: f * (sample_rate / N)
     *
     * In this example, N=8 means we have 8 samples and can analyze
     * 8 frequency bins (0 through 7).
     *
     * Note: In real audio applications, N is typically a power of 2
     * (256, 512, 1024, 2048, etc.) for efficient FFT computation.
     */
    size_t n = 8;

    /**
     * Input buffer: Contains our time-domain signal (amplitude vs time)
     * Output buffer: Will contain frequency-domain data (amplitude vs frequency)
     *
     * The transformation is:
     *   Time Domain (in[]) ──DFT──> Frequency Domain (out[])
     *   "What does the signal look like over time?"
     *                       ──>
     *   "What frequencies make up this signal?"
     */
    float in[n];           /* Input: real-valued signal samples at each time step */
    float complex out[n];  /* Output: complex frequency components (real + imaginary) */

    /**
     * =========================================================================
     * STEP 1: SYNTHESIZE A TEST SIGNAL (Mixing Frequencies)
     * =========================================================================
     *
     * We create a test signal by adding two waves together:
     *   - One COSINE wave at frequency 1 (completes 1 cycle over our sample window)
     *   - One SINE wave at frequency 3 (completes 3 cycles over our sample window)
     *
     * This is like playing two notes simultaneously on a piano.
     * Using both cos and sin demonstrates how the DFT separates them into
     * real (cosine) and imaginary (sine) components.
     *
     * The formula: signal(t) = cos(2π * 1 * t) + sin(2π * 3 * t)
     *
     * Where t ranges from 0 to 1 (normalized time, representing one full
     * period of our sampling window).
     *
     * VISUALIZATION (8 samples of our mixed signal):
     *
     * Sample 0: t=0.000  →  cos(0) + sin(0) = 1 + 0 = 1
     * Sample 1: t=0.125  →  cos(π/4) + sin(3π/4) ≈ 0.707 + 0.707 = 1.414
     * Sample 2: t=0.250  →  cos(π/2) + sin(3π/2) = 0 + (-1) = -1
     * Sample 3: t=0.375  →  cos(3π/4) + sin(9π/4) ≈ -0.707 + 0.707 = 0
     * Sample 4: t=0.500  →  cos(π) + sin(3π) = -1 + 0 = -1
     * Sample 5: t=0.625  →  cos(5π/4) + sin(15π/4) ≈ -0.707 + (-0.707) = -1.414
     * Sample 6: t=0.750  →  cos(3π/2) + sin(9π/2) = 0 + 1 = 1
     * Sample 7: t=0.875  →  cos(7π/4) + sin(21π/4) ≈ 0.707 + (-0.707) = 0
     */
    for (size_t i = 0; i < n; i++) {
        /**
         * t = normalized time, ranging from 0 to (n-1)/n
         *
         * Example with n=8:
         *   i=0: t = 0/8 = 0.000
         *   i=1: t = 1/8 = 0.125
         *   i=2: t = 2/8 = 0.250
         *   ...
         *   i=7: t = 7/8 = 0.875
         *
         * Note: t never reaches 1.0, which prevents sampling the same
         * point twice in a periodic signal.
         */
        float t = (float)i / n;

        /**
         * Add a cosine wave and a sine wave:
         *
         * cosf(2*PI*t*1): Frequency 1 as COSINE - will appear in REAL part of DFT
         *                 Completes 1 full cycle as t goes 0→1
         *
         * sinf(2*PI*t*3): Frequency 3 as SINE - will appear in IMAGINARY part of DFT
         *                 Completes 3 full cycles as t goes 0→1
         *
         * The '2*PI*f*t' pattern comes from the general wave formula:
         *   cos(2πft) or sin(2πft) where f is frequency and t is time
         *
         * 2π radians = 360° = one complete cycle
         *
         * WHY MIX COS AND SIN?
         * This demonstrates a key property: the DFT separates phase information.
         * - Cosine components → appear in the REAL part of the output
         * - Sine components → appear in the IMAGINARY part of the output
         */
        in[i] = cosf(2 * PI * t * 1) + sinf(2 * PI * t * 3);
    }

    /**
     * =========================================================================
     * STEP 2: DISCRETE FOURIER TRANSFORM (Unmixing Frequencies)
     * =========================================================================
     *
     * Now we analyze our mixed signal to find out which frequencies are
     * present. This is the "frequency detection" or "spectral analysis" step.
     *
     * ALGORITHM:
     * For each frequency f from 0 to n-1:
     *   1. Create a complex "probe" using Euler's formula: e^(i*2π*f*t)
     *   2. Multiply each sample of our input by this complex exponential
     *   3. Sum all the products (complex accumulation)
     *   4. Store the complex result for frequency f
     *
     * WHY THIS WORKS (Correlation/Inner Product):
     *
     * When two sine waves of the SAME frequency are multiplied together,
     * the result is always positive (positive×positive = positive,
     * negative×negative = positive). The sum accumulates to a large value.
     *
     * When two sine waves of DIFFERENT frequencies are multiplied,
     * the result oscillates between positive and negative values.
     * Over a complete period, these cancel out to approximately zero.
     *
     * Mathematical basis: This is the orthogonality property of sinusoids.
     * Sine waves at different frequencies are "orthogonal" - their inner
     * product (sum of element-wise multiplication) is zero.
     *
     * EXPECTED OUTPUT for our test signal (frequencies 1 and 3):
     *   Frequency 0: ~0 (no DC component)
     *   Frequency 1: Large positive or negative value (PRESENT!)
     *   Frequency 2: ~0 (not in signal)
     *   Frequency 3: Large positive or negative value (PRESENT!)
     *   Frequency 4: ~0 (not in signal)
     *   Frequency 5: ~0 (this is actually the "mirror" of frequency 3)
     *   Frequency 6: ~0 (this is actually the "mirror" of frequency 2)
     *   Frequency 7: ~0 (this is actually the "mirror" of frequency 1)
     *
     * NOTE ON SYMMETRY (Nyquist Theorem):
     * Due to the nature of real-valued signals, the DFT output is symmetric
     * around N/2. Frequencies above N/2 are "aliases" or mirrors of lower
     * frequencies. This is why in practice we only look at bins 0 to N/2.
     */
    for (size_t f = 0; f < n; f++) {
        /**
         * Initialize accumulator for this frequency bin to zero.
         * We'll sum up all the correlation products into this variable.
         */
        out[f] = 0;

        /**
         * Inner loop: Correlate input signal with complex exponential at frequency f
         *
         * For each sample in the input signal, multiply it by the
         * corresponding complex exponential (via Euler's formula), then add to sum.
         */
        for (size_t i = 0; i < n; i++) {
            /**
             * t = normalized time position for this sample
             * Same calculation as when we generated the input signal.
             */
            float t = (float)i / n;

            /**
             * THE CORE DFT OPERATION using Euler's Formula:
             *
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
             * This is equivalent to (but more elegant than):
             *   real[f] += in[i] * cosf(2*PI*f*t);
             *   imag[f] += in[i] * sinf(2*PI*f*t);
             *
             * NOTE: Standard DFT uses e^(-i*2πft) (negative exponent).
             * Using positive exponent here gives conjugate results but
             * still correctly identifies frequencies. The sign convention
             * affects phase interpretation but not magnitude detection.
             */
            out[f] += in[i] * cexp(2 * I * PI * f * t);
        }
    }

    /**
     * =========================================================================
     * STEP 3: DISPLAY RESULTS
     * =========================================================================
     *
     * Print the frequency spectrum - showing the magnitude of each
     * frequency component detected in our signal.
     *
     * Expected output format: "FF: real, imag"
     *
     * For our input signal cos(2πt) + sin(6πt):
     *   00: ~0.00, ~0.00   ← Frequency 0 (DC offset) - not present
     *   01: ~4.00, ~0.00   ← Frequency 1 - PRESENT as COSINE (real part non-zero!)
     *   02: ~0.00, ~0.00   ← Frequency 2 - not present
     *   03: ~0.00, ~4.00   ← Frequency 3 - PRESENT as SINE (imag part non-zero!)
     *   04: ~0.00, ~0.00   ← Frequency 4 - not present
     *   05: ~0.00, ~-4.00  ← Frequency 5 - mirror of frequency 3
     *   06: ~0.00, ~0.00   ← Frequency 6 - mirror of frequency 2
     *   07: ~4.00, ~0.00   ← Frequency 7 - mirror of frequency 1
     *
     * KEY INSIGHT:
     * - Frequency 1 (cosine input) → energy in REAL part
     * - Frequency 3 (sine input) → energy in IMAGINARY part
     * This demonstrates how the complex DFT preserves phase information!
     *
     * NOTE: The actual magnitude values depend on:
     *   - The amplitude of the input sinusoids (we used 1.0)
     *   - The number of samples (n=8)
     *   - Normalization (this implementation is unnormalized)
     *
     * In normalized form, the values would be divided by N or sqrt(N).
     */
    for (size_t f = 0; f < n; f++) {
        /**
         * Print both real and imaginary components of each frequency bin.
         *
         * Format: "FF: real, imag"
         *   %02zu = zero-padded, 2-digit unsigned integer (size_t)
         *   %.2f = floating point with 2 decimal places
         *
         * creal(out[f]) = real part (cosine correlation)
         * cimag(out[f]) = imaginary part (sine correlation)
         *
         * To get magnitude (total energy at frequency f):
         *   magnitude = sqrt(creal(out[f])² + cimag(out[f])²)
         *   or use: cabs(out[f])
         *
         * To get phase (timing/offset of the frequency component):
         *   phase = atan2(cimag(out[f]), creal(out[f]))
         *   or use: carg(out[f])
         */
        printf("%02zu: %.2f, %.2f\n", f, creal(out[f]), cimag(out[f]));
    }

    return 0;
}

/**
 * =============================================================================
 * FURTHER READING & ADVANCED TOPICS
 * =============================================================================
 *
 * 1. FAST FOURIER TRANSFORM (FFT):
 *    The Cooley-Tukey algorithm reduces DFT from O(N²) to O(N log N).
 *    It works by recursively breaking down the DFT into smaller DFTs,
 *    exploiting the periodicity and symmetry of the complex exponentials.
 *
 * 2. WINDOWING:
 *    Real-world signals aren't perfectly periodic within our sample window.
 *    "Windowing" (Hamming, Hanning, Blackman) tapers the signal at edges
 *    to reduce spectral leakage - false frequencies appearing in the output.
 *
 * 3. FREQUENCY RESOLUTION:
 *    Resolution = Sample_Rate / N
 *    More samples = finer frequency resolution, but more computation
 *    and more temporal "smearing" (uncertainty principle).
 *
 * 4. NYQUIST FREQUENCY:
 *    The highest frequency we can detect is Sample_Rate / 2.
 *    Frequencies above this will "alias" to lower frequencies.
 *    This is why audio is sampled at 44.1kHz - to capture up to ~22kHz.
 *
 * 5. COMPLEX DFT:
 *    The complete DFT uses complex numbers:
 *    X[k] = Σ x[n] * e^(-i*2π*k*n/N)
 *
 *    This captures both magnitude AND phase information:
 *    - Magnitude = sqrt(real² + imag²) = "how much" of each frequency
 *    - Phase = atan2(imag, real) = "when" each frequency component starts
 *
 * 6. APPLICATIONS IN AUDIO VISUALIZATION:
 *    - Music visualizers use FFT to get frequency spectrum in real-time
 *    - Typical setup: 2048 samples @ 44.1kHz = ~46ms windows, ~21Hz resolution
 *    - Low frequencies (bass) → left side of visualization
 *    - High frequencies (treble) → right side of visualization
 *    - Magnitude values drive bar heights or color intensity
 *
 * =============================================================================
 */
