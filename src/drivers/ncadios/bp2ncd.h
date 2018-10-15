/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#define ERR(e){if(e){printf("Error:%s\n",nc_strerror(e));return 2;}}
enum verbose_level
{
	NO_INFO = 0,
	LIST_INFO = 1,
	DEBUG_INFO = 2
};
struct var_dims_struct 
{
	uint16_t id;
	uint64_t rank;
	int	 nc_dimid;
	char	 varname[256];
};

struct time_slice_struct
{
	uint32_t from;
	uint32_t to;
};
void copy_buffer(struct adios_bp_buffer_struct_v1 *dest
                ,struct adios_bp_buffer_struct_v1 *src) {
    
    memcpy (dest, src, sizeof(struct adios_bp_buffer_struct_v1));
}
int assign_value_ulonglong(uint64_t *, enum ADIOS_DATATYPES type, void *);

int add_vardims_var(int ncid,
		    struct var_dims_struct *var_dims,
		    int *var_dims_count, 
		    struct adios_var_header_struct_v1 *var_header, 
		    struct adios_var_payload_struct_v1 *var_payload, 
		    struct adios_bp_buffer_struct_v1 *b);

int add_vardims_attribute(int ncid,
			  struct var_dims_struct* var_dims,
                      	  int *var_dims_count,
			  struct adios_attribute_struct_v1 *attribute);
int makencd(char *bp_fname, char *ncd_fname,int,int);
int parse_cmdline(int argc,
		  char **argv,
		  char **bp_fname,
		  char **ncd_fname,
		  uint32_t *from_timeindex,
		  uint32_t *to_timeindexj,
		  enum verbose_level *verb);
void print_usage(); 

int ncd_gen_name (char *fullname, char *path, char *name); 
int ncd_dataset(int ncid, 
	     	struct adios_var_header_struct_v1 *var_header, 
	     	void *var_payload, 
	     	struct var_dims_struct *var_dims, 
	     	int var_dims_count);
int ncd_scalar(int ncid, 
	       struct adios_var_header_struct_v1 *var_header,
	       void *var_payload); 

int ncd_attr(int ncid, 
             struct adios_attribute_struct_v1 * attribute,
             struct var_dims_struct * var_dims,
	     int var_dims_count); 
