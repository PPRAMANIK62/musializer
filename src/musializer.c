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

/*
 * Declare one variable per plug function via the X-macro.
 *
 * Hot-reload build (-DHOTRELOAD):
 *   Expands to:  plug_hello_t *plug_hello = NULL;  etc.
 *   Each is a function pointer, initially NULL, filled by dlsym() inside
 *   reload_libplug().
 *
 * Static build (no -DHOTRELOAD):
 *   Expands to:  plug_hello_t plug_hello;  etc.
 *   These are extern declarations resolved at link time from plug.c,
 *   so no dlopen/dlsym is needed and reload_libplug() is a no-op.
 */
#ifdef HOTRELOAD
#define PLUG(name) name##_t *name = NULL;
#else
#define PLUG(name) name##_t name;
#endif
LIST_OF_PLUGS
#undef PLUG

Plug plug = {0};

#ifdef HOTRELOAD
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

#define PLUG(name)                                                             \
    name = dlsym(libplug, #name);                                              \
    if (name == NULL) {                                                        \
        fprintf(stderr, "ERROR: could not find %s symbol in %s: %s", #name,    \
                libplug_file_name, dlerror());                                 \
        return false;                                                          \
    }
    LIST_OF_PLUGS
#undef PLUG

    return true;
}
#else
#define reload_libplug() true
#endif

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

    /* Detach callback before unloading so the audio thread doesn't call
     * into memory that is about to be freed. */
    plug_pre_reload(&plug);

    UnloadMusicStream(plug.music);
    CloseAudioDevice();
    CloseWindow();

#ifdef HOTRELOAD
    if (libplug != NULL)
        dlclose(libplug);
#endif

    return 0;
}
