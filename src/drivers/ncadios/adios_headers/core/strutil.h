#ifndef STRUTIL_H_
#define STRUTIL_H_

// copy an array of strings with allocation, return pointer
// also return the sum of string lengths in 'total_length'
char ** a2s_dup_string_array (const char ** v, int nelems, int * total_length);
void a2s_free_string_array (char ** v, int nelems);

void a2s_alloc_namelist (char ***namelist, int length);
void a2s_free_namelist (char **namelist, int length);

/* Remove leading and trailing white spaces.
 * They return a pointer inside str, they do not allocate new memory.
 */
char * a2s_trimL  (char * str);
char * a2s_trimR  (char * str);
char * a2s_trimLR (char * str);

/* Create a char** array of elements from a "n,m,z" type string used as dimension
 * specification in adios. It trims all elements for easy processing later
 */
void a2s_tokenize_dimensions (const char * str, char *** tokens, int * count);

/* Free the dimension char ** array which was created by tokenize_dimensions.
 * It expects the pointers to array and counter, so it can set them to NULL / 0.
 */
void a2s_cleanup_dimensions (char ** tokens, int count);

char * a2s_trim_spaces (const char * str);

/*******************************************************
   Processing parameter lists
**********************************************************/
/*
   Process a ;-separated and possibly multi-line text and 
   create a list of name=value pairs from each 
   item which has a "name=value" pattern. Whitespaces are removed. 
   Input is not modified. Space is allocated;
   Also, simple "name" or "name=" patterns are processed and 
   returned with value=NULL. 
*/
struct PairStruct {
    char * name;
    char * value;
    struct PairStruct * next;
};
typedef struct PairStruct PairStruct;

PairStruct * a2s_text_to_name_value_pairs (const char * text);
void a2s_free_name_value_pairs (PairStruct * pairs);


#endif
