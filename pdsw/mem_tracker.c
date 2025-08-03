#include "mem_tracker.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

// Static variables to track memory usage
static size_t current_memory_usage = 0;
static size_t max_memory_usage = 0;

// Hash table constants
#define HASH_TABLE_SIZE (1024 * 64)
#define CHUNK_SIZE 1024  // Array grows by this size

// Allocation metadata
typedef struct {
    void* ptr;
    size_t size;
} AllocationNode;

// Bucket for handling collisions
typedef struct {
    AllocationNode* nodes;  // Array of allocation nodes
    size_t capacity;        // Current capacity of the array
    size_t count;           // Number of entries in the array
} Bucket;

// Hash table for allocation tracking
static Bucket hash_table[HASH_TABLE_SIZE] = {0};

// Hash function for pointers
static size_t hash_ptr(void* ptr) {
    uintptr_t address = (uintptr_t)ptr;  // Convert pointer to an integer
    address = (address >> 4) ^ (address << 8);  // Mix bits by shifting
    return (address * 2654435761u) % HASH_TABLE_SIZE;
}

// Add allocation to the hash table
static void add_allocation(void* ptr, size_t size) {
    size_t index = hash_ptr(ptr);
    Bucket* bucket = &hash_table[index];

    // Allocate memory for the bucket if it is empty
    if (bucket->nodes == NULL) {
        bucket->nodes = (AllocationNode*)calloc(CHUNK_SIZE, sizeof(AllocationNode));
        if (!bucket->nodes) {
            fprintf(stderr, "Error: Memory allocation failed for bucket.\n");
            exit(EXIT_FAILURE);
        }
        bucket->capacity = CHUNK_SIZE;
        bucket->count = 0;
    }

    // Resize the bucket array if needed
    if (bucket->count >= bucket->capacity) {
        size_t new_capacity = bucket->capacity + CHUNK_SIZE;
        AllocationNode* new_nodes = (AllocationNode*)realloc(bucket->nodes, new_capacity * sizeof(AllocationNode));
        if (!new_nodes) {
            fprintf(stderr, "Error: Memory reallocation failed for bucket.\n");
            exit(EXIT_FAILURE);
        }
        bucket->nodes = new_nodes;
        bucket->capacity = new_capacity;
    }

    // Add the new allocation to the bucket
    bucket->nodes[bucket->count].ptr = ptr;
    bucket->nodes[bucket->count].size = size;
    bucket->count++;
}

// Remove allocation from the hash table and return its size
static size_t remove_allocation(void* ptr) {
    size_t index = hash_ptr(ptr);
    Bucket* bucket = &hash_table[index];

    if (bucket->nodes == NULL) {
        fprintf(stderr, "Warning: Attempt to free untracked pointer.\n");
        return 0; // Bucket is empty
    }

    // Find the allocation in the bucket
    for (size_t i = 0; i < bucket->count; i++) {
        if (bucket->nodes[i].ptr == ptr) {
            size_t size = bucket->nodes[i].size;

            // Replace the removed entry with the last entry for O(1) removal
            bucket->nodes[i] = bucket->nodes[bucket->count - 1];
            bucket->count--;

            return size;
        }
    }

    fprintf(stderr, "Warning: Attempt to free untracked pointer.\n");
    return 0; // Pointer not found
}

// Free all allocations in the hash table
void free_allocation_struct(void) {
    for (size_t i = 0; i < HASH_TABLE_SIZE; i++) {
        Bucket* bucket = &hash_table[i];
        if (bucket->nodes) {
            free(bucket->nodes); // Free the bucket's dynamic array
            bucket->nodes = NULL;
            bucket->capacity = 0;
            bucket->count = 0;
        }
    }
}

//Wrapper for strdup
char* tracked_strdup(const char* s) {
    size_t len = strlen(s) + 1;
    char* copy = (char*)tracked_malloc(len);
    if (copy) memcpy(copy, s, len);
    return copy;
}

// Wrapper for malloc
void* tracked_malloc(size_t size) {
    void* ptr = malloc(size);
    if (ptr) {
        add_allocation(ptr, size);
        current_memory_usage += size;
        if (current_memory_usage > max_memory_usage) {
            max_memory_usage = current_memory_usage;
        }
    }
    return ptr;
}

// Wrapper for free
void tracked_free(void* ptr) {
    if (ptr) {
        size_t size = remove_allocation(ptr);
        current_memory_usage -= size;
        free(ptr);
    }
}

// Returns current memory usage
size_t inq_malloc_use(void) {
    return current_memory_usage;
}

// Returns high-water mark
size_t inq_max_malloc_use(void) {
    return max_memory_usage;
}

// Reset high-water mark
void clear_max_malloc(void) {
    max_memory_usage = current_memory_usage;
}


void* tracked_realloc(void* ptr, size_t new_size) {
    if (!ptr) {
        // If ptr is NULL, behave like tracked_malloc
        return tracked_malloc(new_size);
    }

    if (new_size == 0) {
        // If new_size is 0, behave like tracked_free
        tracked_free(ptr);
        return NULL;
    }

    // Find the old size of the memory block
    size_t old_size = remove_allocation(ptr);

    // Reallocate the memory block
    void* new_ptr = realloc(ptr, new_size);
    if (!new_ptr) {
        // If realloc fails, reinsert the old allocation back into tracking
        add_allocation(ptr, old_size);
        fprintf(stderr, "Error: Memory reallocation failed.\n");
        return NULL;
    }

    // Update memory usage
    current_memory_usage += new_size - old_size;
    if (current_memory_usage > max_memory_usage) {
        max_memory_usage = current_memory_usage;
    }

    // Add the new allocation to the tracking system
    add_allocation(new_ptr, new_size);

    return new_ptr;
}