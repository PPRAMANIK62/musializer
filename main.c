#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <raylib.h>

#define ARRAY_LEN(xs) sizeof(xs)/sizeof(xs[0])

typedef struct {
    float left;
    float right;
} Frame;

Frame global_frames[4410]; // Buffer to store upto 4410 frames
size_t global_frames_count = 0;

void callback(void *bufferData, unsigned int frames) {
    size_t capacity = ARRAY_LEN(global_frames);
    if (frames <= capacity - global_frames_count) {
        memcpy(global_frames + global_frames_count, bufferData, sizeof(Frame)*frames);
        global_frames_count += frames;
    } else if (frames <= capacity) {
        memmove(global_frames, global_frames + frames, sizeof(Frame)*(capacity - frames));
        memcpy(global_frames + (capacity - frames), bufferData, sizeof(Frame)*frames);
    } else {
        memcpy(global_frames, bufferData, sizeof(Frame)*capacity);
        global_frames_count = capacity;
    }
}

int main(void) {
    InitWindow(800, 600, "Musializer");
    SetTargetFPS(60);

    InitAudioDevice();
    Music music = LoadMusicStream("Tower of Dreams.mp3");
    assert(music.stream.sampleSize == 32);
    assert(music.stream.channels == 2);
    printf("music.frameCount: %u\n", music.frameCount);
    printf("music.stream.sampleRate: %u\n", music.stream.sampleRate);
    printf("music.stream.sampleSize: %u\n", music.stream.sampleSize);
    printf("music.stream.channels: %u\n", music.stream.channels);
    PlayMusicStream(music);
    SetMusicVolume(music, 0.25f);
    AttachAudioStreamProcessor(music.stream, callback);

    while (!WindowShouldClose()) {
        UpdateMusicStream(music);

        if (IsKeyPressed(KEY_SPACE)) {
            if (IsMusicStreamPlaying(music)) {
                PauseMusicStream(music);
            } else {
                ResumeMusicStream(music);
            }
        }

        float w = (float)GetRenderWidth();
        float h = (float)GetRenderHeight();

        BeginDrawing();
        ClearBackground(CLITERAL(Color) {0x1a, 0x1b, 0x26, 0xFF});

        float cell_width = (float)w/global_frames_count;
        for (size_t i = 0; i < global_frames_count; i++) {
            float t = global_frames[i].left;
            if (t > 0) {
                DrawRectangle(i*cell_width, h/2 - h/2*t, 1, h/2*t, RED);
            } else {
                DrawRectangle(i*cell_width, h/2, 1, h/2*t, RED);
            }
        }

        EndDrawing();
    }
    return 0;
}
