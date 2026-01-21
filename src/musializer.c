/**
 * =============================================================================
 * MUSIALIZER - Audio Frequency Spectrum Visualizer
 * =============================================================================
 *
 * A music visualizer that displays the frequency spectrum in real-time
 * using FFT (Fast Fourier Transform). Uses raylib for graphics and audio.
 *
 * CURRENT STATE:
 * - FFT-based frequency spectrum visualization
 * - N=256 frequency bins (displays 0 to ~11kHz at 44.1kHz sample rate)
 * - Left channel only
 * - Bars grow upward from center of screen
 *
 * =============================================================================
 */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <raylib.h>
#include <complex.h>
#include "fft.h"

#define ARRAY_LEN(xs) sizeof(xs)/sizeof(xs[0])

/**
 * N = FFT size (number of samples and frequency bins)
 *
 * With N=256 at 44100 Hz sample rate:
 *   - Frequency resolution: 44100/256 ≈ 172 Hz per bin
 *   - Bin 0: DC (0 Hz)
 *   - Bin 1: ~172 Hz
 *   - Bin 128: ~22050 Hz (Nyquist frequency)
 *   - We only display bins 0 to N/2 (due to symmetry)
 *
 * Must be a power of 2 for FFT algorithm
 */
#define N 256

/**
 * Frame structure for stereo audio
 * Values are floats in the range [-1.0, 1.0]
 */
typedef struct {
    float left;
    float right;
} Frame;

/**
 * FFT buffers (shared between audio thread and main thread)
 *
 * in[N]:     Time-domain samples (left channel audio)
 * out[N]:    Frequency-domain complex coefficients from FFT
 * max_amp:   Maximum amplitude across all bins (for normalization)
 *
 * NOTE: Race condition between audio callback and render loop.
 * For production code, use mutex or double-buffering.
 */
float in[N];
float complex out[N];
float max_amp;

/**
 * Fast amplitude approximation for complex numbers
 *
 * Instead of computing sqrt(real² + imag²), we use max(|real|, |imag|)
 * This is faster and sufficient for visualization purposes.
 *
 * True magnitude would be: sqrtf(crealf(z)*crealf(z) + cimagf(z)*cimagf(z))
 * Or use cabs(z) from complex.h
 */
float amp(float complex z) {
    float a = fabsf(crealf(z));
    float b = fabsf(cimagf(z));
    if (a < b) return b;
    return a;
}

/**
 * =============================================================================
 * AUDIO CALLBACK - Called by raylib's audio thread
 * =============================================================================
 *
 * Processes incoming audio data and computes FFT for visualization.
 * Runs on a separate audio thread, not the main thread.
 *
 * @param bufferData  Pointer to audio frames (stereo float pairs)
 * @param frames      Number of frames in this batch
 *
 * PROCESSING PIPELINE:
 * 1. Extract left channel samples into in[] buffer
 * 2. Apply FFT to convert time-domain → frequency-domain
 * 3. Find max amplitude for normalization in rendering
 */
void callback(void *bufferData, unsigned int frames) {
    /* Need at least N frames to fill FFT buffer */
    if (frames < N) return;

    Frame *fs = bufferData;

    /* Extract left channel samples for FFT input */
    for (size_t i = 0; i < N; i++) {
        in[i] = fs[i].left;
    }

    /* Compute FFT: time-domain → frequency-domain */
    fft(in, 1, out, N);

    /* Find maximum amplitude for normalization */
    max_amp = 0.0f;
    for (size_t i = 0; i < N; i++) {
        float a = amp(out[i]);
        if (max_amp < a) max_amp = a;
    }
}

/**
 * Argument parsing helper - shifts and returns first argument
 */
char *shift_args(int *argc, char ***argv) {
    assert(*argc > 0);
    char *result = (**argv);
    (*argv) += 1;
    (*argc) -= 1;
    return result;
}

int main(int argc, char **argv) {
    /**
     * =========================================================================
     * ARGUMENT PARSING
     * =========================================================================
     */
    const char *program = shift_args(&argc, &argv);

    if (argc == 0) {
        fprintf(stderr, "Usage: %s <input>\n", program);
        fprintf(stderr, "ERROR: no input file is provided\n");
        return 1;
    }
    const char *file_path = shift_args(&argc, &argv);

    /**
     * =========================================================================
     * INITIALIZATION
     * =========================================================================
     */
    InitWindow(800, 600, "Musializer");
    SetTargetFPS(60);

    InitAudioDevice();
    Music music = LoadMusicStream(file_path);

    /* Debug: Print audio format info */
    printf("music.frameCount: %u\n", music.frameCount);
    printf("music.stream.sampleRate: %u\n", music.stream.sampleRate);
    printf("music.stream.sampleSize: %u\n", music.stream.sampleSize);
    printf("music.stream.channels: %u\n", music.stream.channels);

    /**
     * Verify audio format assumptions
     * Frame struct requires 32-bit float stereo
     */
    assert(music.stream.sampleSize == 32);
    assert(music.stream.channels == 2);

    PlayMusicStream(music);
    SetMusicVolume(music, 0.25f);
    AttachAudioStreamProcessor(music.stream, callback);

    /**
     * =========================================================================
     * MAIN LOOP
     * =========================================================================
     */
    while (!WindowShouldClose()) {
        UpdateMusicStream(music);

        /* Spacebar: play/pause toggle */
        if (IsKeyPressed(KEY_SPACE)) {
            if (IsMusicStreamPlaying(music)) {
                PauseMusicStream(music);
            } else {
                ResumeMusicStream(music);
            }
        }

        float w = (float)GetRenderWidth();
        float h = (float)GetRenderHeight();

        /**
         * =====================================================================
         * RENDERING - Frequency Spectrum Visualization
         * =====================================================================
         *
         * Each FFT bin becomes a vertical bar:
         *   - X position: bin index mapped to screen width
         *   - Height: proportional to frequency amplitude
         *   - Bars grow upward from screen center
         *
         * Low frequencies (bass) on the left, high frequencies (treble) on right
         */
        BeginDrawing();
        ClearBackground(CLITERAL(Color) {0x1a, 0x1b, 0x26, 0xFF});

        float cell_width = (float)w/N;
        for (size_t i = 0; i < N; i++) {
            /**
             * Get normalized amplitude for this frequency bin
             * Dividing by max_amp scales all bars relative to the loudest frequency,
             * giving consistent visual output regardless of overall volume.
             * Result is in range [0, 1] where 1 = loudest frequency bin
             */
            float t = amp(out[i])/max_amp;
            
            /**
             * Draw bar from center upward
             * Y starts at (h/2 - height) so bar grows upward from h/2
             */
            DrawRectangle(i*cell_width, h/2 - h/2*t, cell_width, h/2*t, RED);
        }

        EndDrawing();
    }

    /**
     * =========================================================================
     * CLEANUP
     * =========================================================================
     * TODO: Add proper cleanup:
     *   UnloadMusicStream(music);
     *   CloseAudioDevice();
     *   CloseWindow();
     */
    return 0;
}
