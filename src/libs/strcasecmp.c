#include <ctype.h>

int strcasecmp(const char *s1, const char *s2)
{
    int c1, c2;
    do {
        c1 = tolower( (unsigned char) *s1++ );
        c2 = tolower( (unsigned char) *s2++ );
    } while (c1 == c2 && c1 != 0);
    return c1 - c2;
}
