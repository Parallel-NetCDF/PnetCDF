
/* on bgl the fortran compiler doesn't add underscores to symbols, and the
 * -brename trick used on seaborg and blue doesn't work with the bgl linker */

#define  ncmpi_close_file_           ncmpi_close_file 
#define  ncmpi_initialize_file_      ncmpi_initialize_file 
#define  ncmpi_write_header_info_    ncmpi_write_header_info 
#define  ncmpi_write_lrefine_        ncmpi_write_lrefine 
#define  ncmpi_write_nodetype_       ncmpi_write_nodetype 
#define  ncmpi_write_gid_            ncmpi_write_gid 
#define  ncmpi_write_coord_          ncmpi_write_coord 
#define  ncmpi_write_size_           ncmpi_write_size 
#define  ncmpi_write_bnd_box_        ncmpi_write_bnd_box 
#define  ncmpi_write_bnd_box_min_    ncmpi_write_bnd_box_min 
#define  ncmpi_write_bnd_box_max_    ncmpi_write_bnd_box_max 
#define  ncmpi_write_unknowns_       ncmpi_write_unknowns 
#define  ncmpi_write_header_info_sp_ ncmpi_write_header_info_sp 
#define  ncmpi_write_coord_sp_       ncmpi_write_coord_sp 
#define  ncmpi_write_size_sp_        ncmpi_write_size_sp 
#define  ncmpi_write_bnd_box_sp_     ncmpi_write_bnd_box_sp 
#define  ncmpi_write_bnd_box_min_sp_ ncmpi_write_bnd_box_min_sp 
#define  ncmpi_write_bnd_box_max_sp_ ncmpi_write_bnd_box_max_sp 
#define  ncmpi_write_unknowns_sp_    ncmpi_write_unknowns_sp 
#define  ncmpi_write_lrefine_sp_     ncmpi_write_lrefine_sp
#define  ncmpi_write_nodetype_sp_    ncmpi_write_nodetype_sp
#define  ncmpi_write_gid_sp_         ncmpi_write_gid_sp

