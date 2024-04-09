include_directories(include)

# these need to be excluded
# src/alloc-override.c
# src/free.c
# src/page-queue.c
# src/static.c

add_library(
        mimalloc OBJECT
        src/alloc-aligned.c
        src/alloc-posix.c
        src/alloc.c
        src/arena.c
        src/bitmap.c
        src/bitmap.h
        src/heap.c
        src/init.c
        src/libc.c
        src/options.c
        src/os.c
        src/page.c
        src/random.c
        src/segment-map.c
        src/segment.c
        src/stats.c
        src/prim/prim.c)

set(MIMALLOC_OBJECT_FILES
        ${MIMALLOC_OBJECT_FILES} $<TARGET_OBJECTS:mimalloc>
        PARENT_SCOPE)

disable_target_warnings(mimalloc)
