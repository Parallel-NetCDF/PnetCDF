#include <ctype.h>

#define CHECK_NAME(A) (strlen(A) > 0 && (isalpha(A[0]) || A[0] == '_'))