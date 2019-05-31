#include <stdlib.h>
#include <string.h>

/* strdup() is a POSIX function, not a standard C function */
char *strdup(const char *str)
{
    char *ptr;
    if (str == NULL) return NULL;

    ptr = (char*) malloc(strlen(str) + 1);
    if (ptr != NULL)
        strcpy(ptr, str);

    return ptr;
}

