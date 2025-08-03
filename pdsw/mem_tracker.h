#ifndef MEM_TRACKER_H
#define MEM_TRACKER_H

#include <stdlib.h>

// Wrapper functions
char* tracked_strdup(const char* s); // Duplicate a string and track it
void* tracked_malloc(size_t size);   // Allocate memory and track it
void tracked_free(void* ptr);        // Free memory and update tracking

// Query functions
size_t inq_malloc_use(void);         // Get current memory usage
size_t inq_max_malloc_use(void);     // Get high-water mark for memory usage

// Reset high-water mark
void clear_max_malloc(void);         // Reset the high-water mark

// Free the tracking data structure
void free_allocation_struct(void);   // Free the internal allocation tracking list


void* tracked_realloc(void* ptr, size_t new_size); // Reallocate memory and track it

#endif // MEM_TRACKER_H