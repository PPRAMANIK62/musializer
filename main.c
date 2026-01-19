/**
 * =============================================================================
 * MUSIALIZER - Audio Waveform Visualizer
 * =============================================================================
 *
 * A simple music visualizer that displays the audio waveform in real-time.
 * Uses raylib for graphics and audio processing.
 *
 * CURRENT STATE:
 * - Displays raw audio waveform (amplitude vs time)
 * - Left channel only visualization
 * - Future: Add FFT for frequency spectrum visualization
 *
 * =============================================================================
 */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <raylib.h>

/* Macro to calculate array length at compile time */
#define ARRAY_LEN(xs) sizeof(xs)/sizeof(xs[0])

/**
 * Frame structure for stereo audio
 * Each frame contains one sample per channel (left and right)
 * Values are floats in the range [-1.0, 1.0]
 */
typedef struct {
    float left;
    float right;
} Frame;

/**
 * Ring buffer to store recent audio frames for visualization
 *
 * WHY 4410 FRAMES?
 * At 44100 Hz sample rate, 4410 frames = 0.1 seconds of audio.
 * This gives us enough data to display a smooth waveform without
 * using too much memory.
 *
 * This buffer is shared between:
 * - Audio thread (writes via callback)
 * - Main thread (reads for rendering)
 *
 * NOTE: This is a potential race condition! For production code,
 * consider using a mutex or lock-free ring buffer.
 */
Frame global_frames[4410];
size_t global_frames_count = 0;

/**
 * =============================================================================
 * AUDIO CALLBACK - Called by raylib's audio thread
 * =============================================================================
 *
 * This function is called automatically by raylib whenever new audio data
 * is available. It runs on a separate audio thread, not the main thread.
 *
 * @param bufferData  Pointer to audio frames (cast to Frame*)
 * @param frames      Number of frames in this batch
 *
 * BUFFER MANAGEMENT STRATEGY:
 * ---------------------------
 * We implement a sliding window / ring buffer approach:
 *
 * Case 1: Buffer has space → Append new frames
 *   [existing data    |  empty space  ]
 *   [existing data    | new frames    ]
 *
 * Case 2: Buffer full, new data fits → Shift left, append
 *   [old | existing data              ]
 *   [existing data        | new frames]
 *
 * Case 3: New data larger than buffer → Keep only newest
 *   [      new frames (truncated)     ]
 */
void callback(void *bufferData, unsigned int frames) {
    size_t capacity = ARRAY_LEN(global_frames);
    
    if (frames <= capacity - global_frames_count) {
        /**
         * CASE 1: Buffer has enough space
         * Simply append new frames to the end
         */
        memcpy(global_frames + global_frames_count, bufferData, sizeof(Frame)*frames);
        global_frames_count += frames;
    } else if (frames <= capacity) {
        /**
         * CASE 2: Buffer is full or nearly full, but new data fits in total capacity
         * Shift existing data left to make room, then append new frames
         * This creates a sliding window effect
         */
        memmove(global_frames, global_frames + frames, sizeof(Frame)*(capacity - frames));
        memcpy(global_frames + (capacity - frames), bufferData, sizeof(Frame)*frames);
        /* global_frames_count stays at capacity */
    } else {
        /**
         * CASE 3: New data is larger than entire buffer capacity
         * Just copy what fits (most recent samples)
         */
        memcpy(global_frames, bufferData, sizeof(Frame)*capacity);
        global_frames_count = capacity;
    }
}

int main(void) {
    /**
     * =========================================================================
     * INITIALIZATION
     * =========================================================================
     */
    
    /* Create window: 800x600 pixels */
    InitWindow(800, 600, "Musializer");
    SetTargetFPS(60);  /* Limit to 60 FPS to save CPU */

    /**
     * Initialize audio subsystem and load music file
     *
     * LoadMusicStream streams from disk (doesn't load entire file to memory)
     * Good for large audio files
     */
    InitAudioDevice();
    Music music = LoadMusicStream("Tower of Dreams.mp3");
    
    /**
     * Verify audio format assumptions
     * Our Frame struct assumes:
     * - 32-bit float samples (sampleSize = 32)
     * - Stereo (channels = 2)
     */
    assert(music.stream.sampleSize == 32);
    assert(music.stream.channels == 2);
    
    /* Debug: Print audio format info */
    printf("music.frameCount: %u\n", music.frameCount);
    printf("music.stream.sampleRate: %u\n", music.stream.sampleRate);
    printf("music.stream.sampleSize: %u\n", music.stream.sampleSize);
    printf("music.stream.channels: %u\n", music.stream.channels);
    
    PlayMusicStream(music);
    SetMusicVolume(music, 0.25f);  /* 25% volume */
    
    /**
     * Attach our callback to receive audio data
     * The callback will be called on the audio thread whenever
     * new audio frames are processed
     */
    AttachAudioStreamProcessor(music.stream, callback);

    /**
     * =========================================================================
     * MAIN LOOP
     * =========================================================================
     */
    while (!WindowShouldClose()) {
        /**
         * Update music stream - REQUIRED every frame
         * This feeds audio data to the sound card and triggers our callback
         */
        UpdateMusicStream(music);

        /**
         * Handle spacebar for play/pause toggle
         */
        if (IsKeyPressed(KEY_SPACE)) {
            if (IsMusicStreamPlaying(music)) {
                PauseMusicStream(music);
            } else {
                ResumeMusicStream(music);
            }
        }

        /* Get current window dimensions (handles resizing) */
        float w = (float)GetRenderWidth();
        float h = (float)GetRenderHeight();

        /**
         * =====================================================================
         * RENDERING - Waveform Visualization
         * =====================================================================
         */
        BeginDrawing();
        
        /* Dark background color (Tokyo Night theme: #1a1b26) */
        ClearBackground(CLITERAL(Color) {0x1a, 0x1b, 0x26, 0xFF});

        /**
         * Draw the waveform
         *
         * Each audio sample becomes a vertical bar:
         * - X position: sample index mapped to screen width
         * - Y position: centered at h/2 (middle of screen)
         * - Height: proportional to amplitude (sample value)
         *
         * Positive samples: draw upward from center
         * Negative samples: draw downward from center
         *
         * NOTE: There's a bug here! When t < 0, the height (h/2*t) is negative.
         * DrawRectangle expects positive height. This should be:
         *   DrawRectangle(x, h/2, 1, -h/2*t, RED);  // negate to get positive height
         */
        float cell_width = (float)w/global_frames_count;
        for (size_t i = 0; i < global_frames_count; i++) {
            float t = global_frames[i].left;  /* Using left channel only */
            if (t > 0) {
                /* Positive amplitude: draw bar upward from center */
                DrawRectangle(i*cell_width, h/2 - h/2*t, 1, h/2*t, RED);
            } else {
                /* Negative amplitude: draw bar downward from center */
                /* BUG: height is negative here, should use -h/2*t */
                DrawRectangle(i*cell_width, h/2, 1, h/2*t, RED);
            }
        }

        EndDrawing();
    }
    
    /**
     * =========================================================================
     * CLEANUP (TODO)
     * =========================================================================
     * Missing cleanup calls:
     *   UnloadMusicStream(music);
     *   CloseAudioDevice();
     *   CloseWindow();
     */
    return 0;
}
