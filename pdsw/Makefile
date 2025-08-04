PNETCDF_DIR=/path/to/PnetCDF/install
PNETCDF_DIR_LIB=/path/to/PnetCDF/library-level-baseline-install
PNETCDF_DIR_FORMAT=/path/to/PnetCDF/new-format-baseline-install
CC = mpicc
H5CC = h5pcc

CFLAGS = -O0 -g
ifeq ($(MEM_TRACK),1)
    CFLAGS += -DMEM_TRACKING
endif
INCLUDES = -I$(PNETCDF_DIR)/include
INCLUDES_LIB = -I$(PNETCDF_DIR_LIB)/include
INCLUDES_FORMAT = -I$(PNETCDF_DIR_FORMAT)/include

LFLAGS = -L$(PNETCDF_DIR)/lib
LFLAGS_LIB = -L$(PNETCDF_DIR_LIB)/lib
LFLAGS_FORMAT = -L$(PNETCDF_DIR_FORMAT)/lib
LIBS = -lpnetcdf

SRCS = baseline_ncx_app.c create_ncopy_binary.c app_baseline_test_all.c app_baseline_read_test_all.c lib_baseline_test_all.c h5_baseline_test_all.c new_format_test_all.c new_format_read_test_all.c mem_tracker.c
OBJS = $(SRCS:.c=.o)

LIB_PROGRAMS = binary_to_nc lib_level_baseline_test_shared lib_level_baseline_test_read lib_level_baseline_test_dup_name lib_baseline_test_all 

all: create_ncopy_binary app_baseline_test_all app_baseline_read_test_all h5_baseline_test_all lib_baseline_test_all new_format_test_all 

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

lib_%.o: lib_%.c
	$(CC) $(CFLAGS) $(INCLUDES_LIB) -c $< -o $@

new_format_%.o: new_format_%.c
	$(CC) $(CFLAGS) $(INCLUDES_FORMAT) -c $< -o $@

app_baseline_test_all: app_baseline_test_all.o baseline_ncx_app.o $(PNETCDF_DIR)/lib/libpnetcdf.a mem_tracker.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

app_baseline_read_test_all: app_baseline_read_test_all.o baseline_ncx_app.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $@ $^ $(LIBS)

lib_baseline_test_all: lib_baseline_test_all.o baseline_ncx_app.o $(PNETCDF_DIR_LIB)/lib/libpnetcdf.a mem_tracker.o
	$(CC) $(CFLAGS) $(INCLUDES_LIB) $(LFLAGS_LIB) -o $@ $^ $(LIBS)

new_format_test_all: new_format_test_all.o baseline_ncx_app.o $(PNETCDF_DIR_FORMAT)/lib/libpnetcdf.a mem_tracker.o
	$(CC) $(CFLAGS) $(INCLUDES_FORMAT) $(LFLAGS_FORMAT) $(LIBS) -o $@ $^ 

new_format_read_test_all: new_format_read_test_all.o baseline_ncx_app.o $(PNETCDF_DIR_FORMAT)/lib/libpnetcdf.a
	$(CC) $(CFLAGS) $(INCLUDES_FORMAT) $(LFLAGS_FORMAT) $(LIBS) -o $@ $^
	
h5_baseline_test_all: h5_baseline_test_all.o baseline_ncx_app.o
	$(H5CC) $(CFLAGS) -o $@ $^ 

create_ncopy_binary: create_ncopy_binary.o baseline_ncx_app.o mem_tracker.o
	$(H5CC) $(CFLAGS) -o $@ $^

h5_%.o: h5_%.c
	$(H5CC) $(CFLAGS) $(INCLUDES) -c $< -o $@



copy:
	cp ./lib_level_baseline_old.nc ./lib_level_baseline.nc
clean:
	rm -f $(OBJS) create_ncopy_binary app_baseline_test_all app_baseline_read_test_all h5_baseline_test_all lib_baseline_test_all new_format_test_all
app: app_baseline_test_all app_baseline_read_test_all
lib: lib_baseline_test_all
h5: h5_baseline_test_all
new_format: new_format_test_all new_format_read_test_all