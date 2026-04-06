/**
 * Musializer - Audio Frequency Spectrum Visualizer
 *
 * A music visualizer that displays the frequency spectrum in real-time
 * using FFT (Fast Fourier Transform). Uses raylib for graphics and audio.
 *
 * Current state:
 * - FFT-based frequency spectrum visualization with logarithmic frequency bins
 * - N=16384 sample sliding window (~2.7 Hz resolution at 44.1kHz)
 * - Left channel only, renders 20 Hz–Nyquist range
 * - Bars grow upward from center of screen
 */

#include "plug.h"
#include <assert.h>
#include <complex.h>
#include <dlfcn.h>
#include <raylib.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define ARRAY_LEN(xs) sizeof(xs) / sizeof(xs[0])

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

const char *libplug_file_name = "libplug.so";
void *libplug = NULL;

plug_hello_t plug_hello = NULL;
plug_init_t plug_init = NULL;
plug_update_t plug_update = NULL;
plug_pre_reload_t plug_pre_reload = NULL;
plug_post_reload_t plug_post_reload = NULL;
Plug plug = {0};

/**
 * Reloads libplug.so and re-binds all function pointers.
 *
 * Called at startup and whenever the user presses R. This is what enables
 * hot-reloading: we close the old shared library, load the freshly compiled
 * one, and resolve each symbol again. The Plug state struct is preserved
 * across reloads since it lives in main(), not inside the library.
 *
 * Symbols resolved: plug_hello, plug_init, plug_pre_reload, plug_post_reload,
 * plug_update. plug_pre_reload/plug_post_reload are called by main() around
 * the reload to safely detach and re-attach the audio callback.
 *
 * dlclose() must come before dlopen() so the OS unloads the old version.
 * RTLD_NOW resolves all symbols immediately (vs RTLD_LAZY which defers until
 * first call) — we want to catch missing symbols at reload time, not later.
 */
bool reload_libplug(void) {
    if (libplug != NULL)
        dlclose(libplug);

    libplug = dlopen(libplug_file_name, RTLD_NOW);
    if (libplug == NULL) {
        fprintf(stderr, "ERROR: could not load %s: %s", libplug_file_name,
                dlerror());
        return false;
    }

    plug_hello = dlsym(libplug, "plug_hello");
    if (plug_hello == NULL) {
        fprintf(stderr, "ERROR: could not find plug_hello symbol in %s: %s",
                libplug_file_name, dlerror());
        return false;
    }

    plug_init = dlsym(libplug, "plug_init");
    if (plug_init == NULL) {
        fprintf(stderr, "ERROR: could not find plug_init symbol in %s: %s",
                libplug_file_name, dlerror());
        return false;
    }

    plug_pre_reload = dlsym(libplug, "plug_pre_reload");
    if (plug_pre_reload == NULL) {
        fprintf(stderr,
                "ERROR: could not find plug_pre_reload symbol in %s: %s",
                libplug_file_name, dlerror());
        return false;
    }

    plug_post_reload = dlsym(libplug, "plug_post_reload");
    if (plug_post_reload == NULL) {
        fprintf(stderr,
                "ERROR: could not find plug_post_reload symbol in %s: %s",
                libplug_file_name, dlerror());
        return false;
    }

    plug_update = dlsym(libplug, "plug_update");
    if (plug_update == NULL) {
        fprintf(stderr, "ERROR: could not find plug_update symbol in %s: %s",
                libplug_file_name, dlerror());
        return false;
    }

    return true;
}

int main(int argc, char **argv) {
    if (!reload_libplug())
        return 1;

    // Argument parsing
    const char *program = shift_args(&argc, &argv);

    // TODO: supply input files via drag&drop
    if (argc == 0) {
        fprintf(stderr, "Usage: %s <input>\n", program);
        fprintf(stderr, "ERROR: no input file is provided\n");
        return 1;
    }
    const char *file_path = shift_args(&argc, &argv);

    // Initialization
    InitWindow(800, 600, "Musializer");
    SetTargetFPS(60);
    InitAudioDevice();

    plug_init(&plug, file_path);

    // Main loop
    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_R)) {
            plug_pre_reload(&plug);
            if (!reload_libplug())
                return 1;
            plug_post_reload(&plug);
        }

        plug_update(&plug);
    }

    // TODO: proper cleanup (UnloadMusicStream, CloseAudioDevice, CloseWindow)
    return 0;
}
