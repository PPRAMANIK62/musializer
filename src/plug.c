#include "plug.h"
#include "fft.h"
#include <assert.h>
#include <math.h>
#include <raylib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

/**
 * N = FFT window size (number of samples kept in the sliding buffer)
 *
 * With N=16384 (1<<14) at 44100 Hz sample rate:
 *   - Frequency resolution: 44100/16384 ≈ 2.7 Hz per bin
 *   - Bin 0: DC (0 Hz)
 *   - Bin 1: ~2.7 Hz
 *   - Bin 8192: ~22050 Hz (Nyquist frequency)
 *   - Rendering starts at 20 Hz (bin ~7) to skip sub-bass noise
 *
 * Must be a power of 2 for the FFT algorithm. Larger N gives finer
 * frequency resolution at the cost of more memory and FFT time.
 */
#define N (1 << 14)

/*
 * why are we doing this ??
 * -> raylib audio processor callback doesn't accept user data
 * -> not hot reloadable
 */
float in[N];
float complex out[N];

/**
 * Frame structure for stereo audio
 * Values are floats in the range [-1.0, 1.0]
 */
typedef struct {
    float left;
    float right;
} Frame;

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
    if (a < b)
        return b;
    return a;
}

/**
 * Audio callback - called by raylib's audio thread, not the main thread.
 *
 * @param bufferData  Pointer to audio frames (stereo float pairs)
 * @param frames      Number of frames in this batch
 *
 * Maintains in[] as a sliding window of the N most recent left-channel
 * samples. For each new sample, the buffer is shifted left by one slot
 * (memmove) and the new sample is placed at the end. This means the FFT
 * always sees a fresh, continuous window rather than discrete chunks.
 */
void callback(void *bufferData, unsigned int frames) {
    Frame *fs = bufferData;

    for (size_t i = 0; i < frames; i++) {
        memmove(in, in + 1, (N - 1) * sizeof(in[0]));
        in[N - 1] = fs[i].left;
    }
}

void plug_hello(void) { printf("Hello from Plugin\n"); }

/**
 * Initializes plugin state. Called once at startup (or after a hot-reload).
 *
 * Loads the music file, validates the audio format (must be 32-bit float
 * stereo — matching our Frame struct), starts playback, and attaches the
 * audio callback so the FFT buffer gets filled each audio frame.
 */
void plug_init(Plug *plug, const char *file_path) {
    plug->music = LoadMusicStream(file_path);

    /* Debug: Print audio format info */
    printf("music.frameCount: %u\n", plug->music.frameCount);
    printf("music.stream.sampleRate: %u\n", plug->music.stream.sampleRate);
    printf("music.stream.sampleSize: %u\n", plug->music.stream.sampleSize);
    printf("music.stream.channels: %u\n", plug->music.stream.channels);

    /**
     * Verify audio format assumptions
     * Frame struct requires 32-bit float stereo
     */
    assert(plug->music.stream.sampleSize == 32);
    assert(plug->music.stream.channels == 2);

    SetMusicVolume(plug->music, 0.25f);
    PlayMusicStream(plug->music);
    AttachAudioStreamProcessor(plug->music.stream, callback);
}

/**
 * Called by main() just before unloading the old libplug.so.
 *
 * Detaches the audio callback so the old (soon-to-be-unloaded) callback
 * pointer is no longer called by the audio thread during the reload window.
 * Skipping this would cause a call into freed/unmapped code.
 */
void plug_pre_reload(Plug *plug) {
    DetachAudioStreamProcessor(plug->music.stream, callback);
}

/**
 * Called by main() immediately after the new libplug.so is loaded.
 *
 * Re-attaches the audio callback using the fresh function pointer from the
 * newly loaded library, resuming audio processing with the updated code.
 */
void plug_post_reload(Plug *plug) {
    AttachAudioStreamProcessor(plug->music.stream, callback);
}

/**
 * Per-frame update — called every frame from the main loop.
 *
 * Input handling:
 *   Space — play/pause toggle
 *   Q     — restart track from the beginning
 *
 * Rendering pipeline:
 *   1. fft(in[]) — transform the sliding window to frequency domain
 *   2. Find max amplitude across all bins for normalization
 *   3. Group FFT bins into logarithmic frequency bands (step=1.06, starting
 *      at 20 Hz) — this mirrors how human hearing perceives pitch, spreading
 *      bass detail rather than compressing it into a few pixels on the left
 *   4. Draw each band as a bar growing upward from the center
 *
 * Note: in[] is written by the audio thread (callback) and read here on the
 * main thread — there is a race condition, acceptable for a visualizer.
 */
void plug_update(Plug *plug) {
    UpdateMusicStream(plug->music);

    /* Spacebar: play/pause toggle */
    if (IsKeyPressed(KEY_SPACE)) {
        if (IsMusicStreamPlaying(plug->music)) {
            PauseMusicStream(plug->music);
        } else {
            ResumeMusicStream(plug->music);
        }
    }

    if (IsKeyPressed(KEY_Q)) {
        StopMusicStream(plug->music);
        PlayMusicStream(plug->music);
    }

    float w = (float)GetRenderWidth();
    float h = (float)GetRenderHeight();

    /**
     * Frequency Spectrum Visualization
     *
     * Each FFT bin becomes a vertical bar:
     *   - X position: bin index mapped to screen width
     *   - Height: proportional to frequency amplitude
     *   - Bars grow upward from screen center
     *
     * Low frequencies (bass) on the left, high frequencies (treble) on right.
     */
    BeginDrawing();
    ClearBackground(CLITERAL(Color){0x1a, 0x1b, 0x26, 0xFF});

    /* Compute FFT: time-domain → frequency-domain */
    fft(in, 1, out, N);

    /* Find maximum amplitude for normalization */
    float max_amp = 0.0f;
    for (size_t i = 0; i < N; i++) {
        float a = amp(out[i]);
        if (max_amp < a)
            max_amp = a;
    }

    float step = 1.06;
    size_t m = 0;
    for (float f = 20.0f; (size_t)f < N; f *= step) {
        m++;
    }

    float cell_width = (float)w / m;
    m = 0;
    for (float f = 20.0f; (size_t)f < N; f *= step) {
        float f1 = f * step;
        float a = 0.0f;
        for (size_t q = (size_t)f; q < N && q < (size_t)f1; q++) {
            a += amp(out[q]);
        }
        a /= (size_t)f1 - (size_t)f + 1;
        /**
         * Get normalized amplitude for this frequency bin
         * Dividing by max_amp scales all bars relative to the loudest
         * frequency, giving consistent visual output regardless of overall
         * volume. Result is in range [0, 1] where 1 = loudest frequency bin
         */
        float t = a / max_amp;

        /**
         * Draw bar from center upward
         * Y starts at (h/2 - height) so bar grows upward from h/2
         */
        DrawRectangle(m * cell_width, h / 2 - h / 2 * t, cell_width, h / 2 * t,
                      RED);
        m++;
    }

    EndDrawing();
}
