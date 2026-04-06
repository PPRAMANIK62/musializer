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

/*
 * Function-type typedefs (not function-pointer typedefs).
 *
 * Writing `void(foo_t)(...)` defines a function *type*, whereas
 * `void (*foo_t)(...)` would define a function-*pointer* type.
 * Using the bare function type lets the X-macro in musializer.c
 * derive both forms from the same name:
 *   - Hot-reload build:  `foo_t *foo`  → a pointer, filled by dlsym()
 *   - Static build:      `foo_t  foo`  → a direct symbol reference
 */
typedef void(plug_hello_t)(void);
typedef void(plug_init_t)(Plug *plug, const char *file_path);
typedef void(plug_pre_reload_t)(Plug *plug);
typedef void(plug_post_reload_t)(Plug *plug);
typedef void(plug_update_t)(Plug *plug);

/*
 * X-macro table of every function exported by libplug.
 *
 * Defining PLUG(name) before expanding LIST_OF_PLUGS lets one macro
 * drive three distinct uses in musializer.c without repetition:
 *
 *   1. Variable declarations  →  `name##_t *name = NULL;`  (hot-reload)
 *                                `name##_t  name;`          (static)
 *   2. dlsym() binding        →  `name = dlsym(libplug, #name);`
 *
 * To add a new plug function: declare it here and implement it in plug.c.
 * musializer.c picks it up automatically — no other changes needed.
 */
#define LIST_OF_PLUGS                                                          \
    PLUG(plug_hello)                                                           \
    PLUG(plug_init)                                                            \
    PLUG(plug_pre_reload)                                                      \
    PLUG(plug_post_reload)                                                     \
    PLUG(plug_update)

#endif // PLUG_H_
