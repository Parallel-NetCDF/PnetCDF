/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header$
 *********************************************************************/

#include <stdio.h>		/* has getopt() under VMS */
#include <string.h>

#ifdef __hpux
#include <locale.h>
#endif
    
#include <pnetcdf.h>

#include "generic.h"
#include "ncgen.h"
#include "genlib.h"

extern int	yyparse(void);

const char *progname;			/* for error messages */
const char *cdlname;

int c_flag;
int fortran_flag;
int netcdf_flag;
char *netcdf_name = NULL;	/* name of output netCDF file to write */

extern FILE *yyin;

static const char* ubasename ( const char* av0 );
static void usage ( void );
int main ( int argc, char** argv );


/* strip off leading path */
static const char *
ubasename(
	const char *av0)
{
	const char *logident ;
#ifdef VMS
#define SEP	']'
#endif
#ifdef MSDOS
#define SEP	'\\'
#endif
#ifndef SEP
#define SEP	'/'
#endif
	if ((logident = strrchr(av0, SEP)) == NULL)
		logident = av0 ;
	else
	    logident++ ;
	return logident ;
}


static void usage(void)
{
    derror("Usage: %s [ -b ] [ -c ] [ -f ] [ -o outfile]  [ file ... ]",
	   progname);
    derror("netcdf library version %s", ncmpi_inq_libvers());
}


int
main(
	int argc,
	char *argv[])
{
    extern int optind;
    extern int opterr;
    extern char *optarg;
    int c;
    FILE *fp;

    MPI_Init(&argc, &argv);

#ifdef __hpux
    setlocale(LC_CTYPE,"");
#endif
    
#ifdef MDEBUG
	malloc_debug(2) ;	/* helps find malloc/free errors on Sun */
#endif /* MDEBUG */

    opterr = 1;			/* print error message if bad option */
    progname = ubasename(argv[0]);
    cdlname = "-";

    c_flag = 0;
    fortran_flag = 0;
    netcdf_flag = 0;

#if _CRAYMPP && 0
    /* initialize CRAY MPP parallel-I/O library */
    (void) par_io_init(32, 32);
#endif

    while ((c = getopt(argc, argv, "bcfno:")) != EOF)
      switch(c) {
	case 'c':		/* for c output */
	  c_flag = 1;
	  break;
	case 'f':		/* for fortran output */
	  fortran_flag = 1;
	  break;
	case 'b':		/* for binary netcdf output, ".nc" extension */
	  netcdf_flag = 1;
	  break;
	case 'n':		/* old version of -b, uses ".cdf" extension */
	  netcdf_flag = -1;
	  break;
	case 'o':		/* to explicitly specify output name */
	  netcdf_flag = 1;
	  netcdf_name = (char *) emalloc(strlen(optarg)+1);
	  if (! netcdf_name) {
	      derror ("%s: out of memory", progname);
	      return(1);
	  }
	  (void)strcpy(netcdf_name,optarg);
	  break;
	case '?':
	  usage();
	  return(8);
      }

    if (fortran_flag && c_flag) {
	derror("Only one of -c or -f may be specified");
	return(8);
      }
    if (fortran_flag) {
	    derror("Generating Fortran interface not supported yet");
            return(0);
    }

    argc -= optind;
    argv += optind;

    if (argc > 1) {
	derror ("%s: only one input file argument permitted",progname);
	return(6);
    }

    fp = stdin;
    if (argc > 0 && strcmp(argv[0], "-") != 0) {
	if ((fp = fopen(argv[0], "r")) == NULL) {
	    derror ("can't open file %s for reading: ", argv[0]);
	    perror("");
	    return(7);
	}
	cdlname = argv[0];
    }
    yyin = fp;
    return (yyparse());
}
