#ifndef PLUG_H_
#define PLUG_H_

#include <complex.h>
#include <raylib.h>
#include <stddef.h>

/**
 * Hot-reloadable plugin state.
 *
 * This struct is allocated in main() and passed into every plug_* function.
 * Because it lives outside the shared library, it survives a reload — the
 * new libplug.so picks up exactly where the old one left off.
 *
 * Any state that needs to persist across reloads must live here.
 */
typedef struct {
    Music music;
} Plug;

typedef void (*plug_hello_t)(void);
typedef void (*plug_init_t)(Plug *plug, const char *file_path);
typedef void (*plug_pre_reload_t)(Plug *plug);
typedef void (*plug_post_reload_t)(Plug *plug);
typedef void (*plug_update_t)(Plug *plug);

#endif // PLUG_H_
