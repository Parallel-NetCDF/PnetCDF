/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0

/* Substitute the variable and function names.  */
#define yyparse ncmpiparse
#define yylex   ncmpilex
#define yyerror ncmpierror
#define yylval  ncmpilval
#define yychar  ncmpichar
#define yydebug ncmpidebug
#define yynerrs ncmpinerrs


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     NC_UNLIMITED_K = 258,
     BYTE_K = 259,
     CHAR_K = 260,
     SHORT_K = 261,
     INT_K = 262,
     FLOAT_K = 263,
     DOUBLE_K = 264,
     IDENT = 265,
     TERMSTRING = 266,
     BYTE_CONST = 267,
     CHAR_CONST = 268,
     SHORT_CONST = 269,
     INT_CONST = 270,
     FLOAT_CONST = 271,
     DOUBLE_CONST = 272,
     DIMENSIONS = 273,
     VARIABLES = 274,
     NETCDF = 275,
     DATA = 276,
     FILLVALUE = 277
   };
#endif
/* Tokens.  */
#define NC_UNLIMITED_K 258
#define BYTE_K 259
#define CHAR_K 260
#define SHORT_K 261
#define INT_K 262
#define FLOAT_K 263
#define DOUBLE_K 264
#define IDENT 265
#define TERMSTRING 266
#define BYTE_CONST 267
#define CHAR_CONST 268
#define SHORT_CONST 269
#define INT_CONST 270
#define FLOAT_CONST 271
#define DOUBLE_CONST 272
#define DIMENSIONS 273
#define VARIABLES 274
#define NETCDF 275
#define DATA 276
#define FILLVALUE 277




/* Copy the first part of user declarations.  */
#line 12 "./ncmpigen.y"

#ifndef lint
static char SccsId[] = "$Id$";
#endif
#include        <string.h>
#include	<stdlib.h>
#include	<stddef.h>  /* ptrdiff_t */
#include	<pnetcdf.h>
#include 	"generic.h"
#include        "ncmpigen.h"
#include	"genlib.h"	/* for grow_darray() et al */

typedef struct Symbol {		/* symbol table entry */
	char    	*name;
	struct Symbol   *next;
	unsigned	is_dim : 1;	/* appears as netCDF dimension */
	unsigned	is_var : 1;	/* appears as netCDF variable */
	unsigned	is_att : 1;	/* appears as netCDF attribute */
	int             dnum;	        /* handle as a dimension */
	int             vnum;	        /* handle as a variable */
	} *YYSTYPE1;

/* True if string a equals string b*/
#define	STREQ(a, b)	(*(a) == *(b) && strcmp((a), (b)) == 0)
#define NC_UNSPECIFIED ((nc_type)0)	/* unspecified (as yet) type */

#define YYSTYPE YYSTYPE1
YYSTYPE symlist;		/* symbol table: linked list */

extern int derror_count;	/* counts errors in netcdf definition */
extern int lineno;		/* line number for error messages */

static int not_a_string;	/* whether last constant read was a string */
static char termstring[MAXTRST]; /* last terminal string read */
static double double_val;	/* last double value read */
static float float_val;		/* last float value read */
static int int_val;		/* last int value read */
static short short_val;		/* last short value read */
static char char_val;		/* last char value read */
static signed char byte_val;	/* last byte value read */

static nc_type type_code;	/* holds declared type for variables */
static nc_type atype_code;	/* holds derived type for attributes */
static char *netcdfname;	/* to construct netcdf file name */
static void *att_space;		/* pointer to block for attribute values */
static nc_type valtype;		/* type code for list of attribute values  */

static char *char_valp;		/* pointers used to accumulate data values */
static signed char *byte_valp;
static short *short_valp;
static int *int_valp;
static float *float_valp;
static double *double_valp;
static void *rec_cur;		/* pointer to where next data value goes */
static void *rec_start;		/* start of space for data */

/* function prototypes */
#ifdef vms                                                                     
void                                                                           
#else                                                                          
int                                                                            
#endif                                                                         
yyerror( char *s);
int yylex(void);


/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef int YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 224 "y.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  4
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   92

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  31
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  40
/* YYNRULES -- Number of rules.  */
#define YYNRULES  77
/* YYNRULES -- Number of states.  */
#define YYNSTATES  108

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   277

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      28,    29,     2,     2,    26,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    30,    25,
       2,    27,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    23,     2,    24,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     4,     5,    14,    15,    18,    21,    25,
      27,    31,    35,    39,    43,    45,    47,    48,    51,    53,
      56,    60,    62,    64,    67,    71,    73,    76,    78,    80,
      82,    84,    86,    88,    90,    94,    95,    99,   101,   102,
     106,   108,   112,   114,   115,   120,   124,   127,   129,   131,
     133,   137,   139,   141,   143,   145,   147,   149,   151,   153,
     154,   157,   159,   162,   166,   167,   172,   174,   178,   179,
     182,   184,   186,   188,   190,   192,   194,   196
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      32,     0,    -1,    -1,    -1,    20,    23,    33,    35,    41,
      34,    63,    24,    -1,    -1,    18,    36,    -1,    37,    25,
      -1,    36,    37,    25,    -1,    38,    -1,    37,    26,    38,
      -1,    39,    27,    15,    -1,    39,    27,    17,    -1,    39,
      27,     3,    -1,    40,    -1,    10,    -1,    -1,    19,    42,
      -1,    44,    -1,    43,    25,    -1,    42,    43,    25,    -1,
      46,    -1,    55,    -1,    45,    25,    -1,    44,    45,    25,
      -1,    55,    -1,    47,    48,    -1,     4,    -1,     5,    -1,
       6,    -1,     7,    -1,     8,    -1,     9,    -1,    49,    -1,
      48,    26,    49,    -1,    -1,    51,    50,    52,    -1,    10,
      -1,    -1,    28,    53,    29,    -1,    54,    -1,    53,    26,
      54,    -1,    40,    -1,    -1,    57,    56,    27,    60,    -1,
      58,    30,    59,    -1,    30,    59,    -1,    51,    -1,    10,
      -1,    61,    -1,    60,    26,    61,    -1,    62,    -1,    13,
      -1,    11,    -1,    12,    -1,    14,    -1,    15,    -1,    16,
      -1,    17,    -1,    -1,    21,    64,    -1,    21,    -1,    65,
      25,    -1,    64,    65,    25,    -1,    -1,    58,    66,    27,
      67,    -1,    68,    -1,    67,    26,    68,    -1,    -1,    69,
      70,    -1,    13,    -1,    11,    -1,    12,    -1,    14,    -1,
      15,    -1,    16,    -1,    17,    -1,    22,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   110,   110,   113,   108,   126,   127,   129,   130,   132,
     133,   135,   141,   152,   160,   176,   178,   179,   180,   182,
     183,   185,   185,   187,   188,   190,   192,   194,   195,   196,
     197,   198,   199,   201,   202,   205,   204,   243,   245,   246,
     248,   249,   251,   273,   272,   317,   318,   324,   334,   345,
     346,   348,   357,   363,   376,   382,   388,   394,   400,   408,
     409,   410,   413,   414,   417,   416,   486,   487,   490,   490,
     544,   572,   630,   658,   686,   714,   742,   773
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "NC_UNLIMITED_K", "BYTE_K", "CHAR_K",
  "SHORT_K", "INT_K", "FLOAT_K", "DOUBLE_K", "IDENT", "TERMSTRING",
  "BYTE_CONST", "CHAR_CONST", "SHORT_CONST", "INT_CONST", "FLOAT_CONST",
  "DOUBLE_CONST", "DIMENSIONS", "VARIABLES", "NETCDF", "DATA", "FILLVALUE",
  "'{'", "'}'", "';'", "','", "'='", "'('", "')'", "':'", "$accept",
  "ncdesc", "@1", "@2", "dimsection", "dimdecls", "dimdecline", "dimdecl",
  "dimd", "dim", "vasection", "vadecls", "vadecl", "gattdecls", "gattdecl",
  "vardecl", "type", "varlist", "varspec", "@3", "var", "dimspec",
  "dimlist", "vdim", "attdecl", "@4", "att", "avar", "attr", "attvallist",
  "aconst", "attconst", "datasection", "datadecls", "datadecl", "@5",
  "constlist", "dconst", "@6", "const", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   123,   125,    59,    44,    61,    40,    41,
      58
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    31,    33,    34,    32,    35,    35,    36,    36,    37,
      37,    38,    38,    38,    39,    40,    41,    41,    41,    42,
      42,    43,    43,    44,    44,    45,    46,    47,    47,    47,
      47,    47,    47,    48,    48,    50,    49,    51,    52,    52,
      53,    53,    54,    56,    55,    57,    57,    58,    59,    60,
      60,    61,    62,    62,    62,    62,    62,    62,    62,    63,
      63,    63,    64,    64,    66,    65,    67,    67,    69,    68,
      70,    70,    70,    70,    70,    70,    70,    70
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     0,     8,     0,     2,     2,     3,     1,
       3,     3,     3,     3,     1,     1,     0,     2,     1,     2,
       3,     1,     1,     2,     3,     1,     2,     1,     1,     1,
       1,     1,     1,     1,     3,     0,     3,     1,     0,     3,
       1,     3,     1,     0,     4,     3,     2,     1,     1,     1,
       3,     1,     1,     1,     1,     1,     1,     1,     1,     0,
       2,     1,     2,     3,     0,     4,     1,     3,     0,     2,
       1,     1,     1,     1,     1,     1,     1,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,     0,     2,     1,     5,     0,    16,    15,     6,
       0,     9,     0,    14,    37,     0,     0,     3,    18,     0,
      47,    25,    43,     0,     0,     7,     0,     0,    27,    28,
      29,    30,    31,    32,    17,     0,    21,     0,    22,    48,
      46,    59,     0,    23,     0,     0,     8,    10,    13,    11,
      12,     0,    19,    26,    33,    35,    61,     0,    24,     0,
      45,    20,     0,    38,    64,    60,     0,     4,    53,    54,
      52,    55,    56,    57,    58,    44,    49,    51,    34,     0,
      36,     0,     0,    62,     0,    42,     0,    40,    68,    63,
      50,     0,    39,    65,    66,     0,    41,    68,    71,    72,
      70,    73,    74,    75,    76,    77,    69,    67
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     2,     5,    41,     7,     9,    10,    11,    12,    13,
      17,    34,    35,    18,    19,    36,    37,    53,    54,    63,
      20,    80,    86,    87,    21,    44,    22,    23,    40,    75,
      76,    77,    57,    65,    66,    81,    93,    94,    95,   106
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -69
static const yytype_int8 yypact[] =
{
       1,     4,    30,   -69,   -69,    23,    32,    -2,   -69,    32,
     -12,   -69,    13,   -69,   -69,    -4,    33,   -69,    -1,    27,
     -69,   -69,   -69,    24,    -6,   -69,    32,     7,   -69,   -69,
     -69,   -69,   -69,   -69,    -4,    28,   -69,    45,   -69,   -69,
     -69,    35,    36,   -69,    31,    33,   -69,   -69,   -69,   -69,
     -69,    37,   -69,    38,   -69,   -69,    45,    39,   -69,    34,
     -69,   -69,    45,    29,   -69,    45,    40,   -69,   -69,   -69,
     -69,   -69,   -69,   -69,   -69,    41,   -69,   -69,   -69,    32,
     -69,    42,    43,   -69,    34,   -69,   -11,   -69,   -69,   -69,
     -69,    32,   -69,    44,   -69,    22,   -69,   -69,   -69,   -69,
     -69,   -69,   -69,   -69,   -69,   -69,   -69,   -69
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -69,   -69,   -69,   -69,   -69,   -69,    50,    46,   -69,   -68,
     -69,   -69,    26,   -69,    48,   -69,   -69,   -69,     9,   -69,
     -30,   -69,   -69,   -18,    -3,   -69,   -69,   -40,    47,   -69,
     -10,   -69,   -69,   -69,    10,   -69,   -69,   -21,   -69,   -69
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
      28,    29,    30,    31,    32,    33,    14,    55,    14,    14,
      48,    85,    38,    25,    26,    91,    64,    15,    92,    46,
      26,     1,    49,    85,    50,    64,    16,     3,    16,    16,
       4,    38,    55,    98,    99,   100,   101,   102,   103,   104,
      27,     6,     8,    39,   105,    68,    69,    70,    71,    72,
      73,    74,    43,    52,    45,    14,    56,    79,    59,    24,
      51,    58,    61,    67,    62,    83,    42,    84,    89,    88,
      97,    78,    47,    96,    90,    82,   107,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    60
};

static const yytype_int8 yycheck[] =
{
       4,     5,     6,     7,     8,     9,    10,    37,    10,    10,
       3,    79,    15,    25,    26,    26,    56,    19,    29,    25,
      26,    20,    15,    91,    17,    65,    30,    23,    30,    30,
       0,    34,    62,    11,    12,    13,    14,    15,    16,    17,
      27,    18,    10,    10,    22,    11,    12,    13,    14,    15,
      16,    17,    25,    25,    30,    10,    21,    28,    27,     9,
      34,    25,    25,    24,    26,    25,    18,    26,    25,    27,
      26,    62,    26,    91,    84,    65,    97,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    45
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    20,    32,    23,     0,    33,    18,    35,    10,    36,
      37,    38,    39,    40,    10,    19,    30,    41,    44,    45,
      51,    55,    57,    58,    37,    25,    26,    27,     4,     5,
       6,     7,     8,     9,    42,    43,    46,    47,    55,    10,
      59,    34,    45,    25,    56,    30,    25,    38,     3,    15,
      17,    43,    25,    48,    49,    51,    21,    63,    25,    27,
      59,    25,    26,    50,    58,    64,    65,    24,    11,    12,
      13,    14,    15,    16,    17,    60,    61,    62,    49,    28,
      52,    66,    65,    25,    26,    40,    53,    54,    27,    25,
      61,    26,    29,    67,    68,    69,    54,    26,    11,    12,
      13,    14,    15,    16,    17,    22,    70,    68
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      fprintf (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 110 "./ncmpigen.y"
    { init_netcdf(); }
    break;

  case 3:
#line 113 "./ncmpigen.y"
    {
		       if (derror_count == 0)
			 define_netcdf(netcdfname);
		       if (derror_count > 0)
		       	exit(6);
		   }
    break;

  case 4:
#line 121 "./ncmpigen.y"
    {
		       if (derror_count == 0)
			 close_netcdf();
		   }
    break;

  case 11:
#line 136 "./ncmpigen.y"
    { if (int_val <= 0)
			 derror("dimension length must be positive");
		     dims[ndims].size = int_val;
		     ndims++;
		   }
    break;

  case 12:
#line 142 "./ncmpigen.y"
    { /* for rare case where 2^31 < dimsize < 2^32 */
		   	if (double_val <= 0)
				derror("dimension length must be positive");
			if (double_val > 4294967295.0)
				derror("dimension too large");
			if (double_val - (MPI_Offset) double_val > 0)
				derror("dimension length must be an integer");
			dims[ndims].size = (MPI_Offset) double_val;
			ndims++;
		   }
    break;

  case 13:
#line 153 "./ncmpigen.y"
    {  if (rec_dim != -1)
			 derror("only one NC_UNLIMITED dimension allowed");
		     rec_dim = ndims; /* the unlimited (record) dimension */
		     dims[ndims].size = NC_UNLIMITED;
		     ndims++;
		   }
    break;

  case 14:
#line 161 "./ncmpigen.y"
    { if ((yyvsp[(1) - (1)])->is_dim == 1) {
		        derror( "duplicate dimension declaration for %s",
		                (yyvsp[(1) - (1)])->name);
		     }
	             (yyvsp[(1) - (1)])->is_dim = 1;
		     (yyvsp[(1) - (1)])->dnum = ndims;
		     /* make sure dims array will hold dimensions */
		     grow_darray(ndims,  /* must hold ndims+1 dims */
				 &dims); /* grow as needed */
		     dims[ndims].name = (char *) emalloc(strlen((yyvsp[(1) - (1)])->name)+1);
		     (void) strcpy(dims[ndims].name, (yyvsp[(1) - (1)])->name);
		     /* name for use in generated Fortran and C variables */
		     dims[ndims].lname = decodify((yyvsp[(1) - (1)])->name);
		   }
    break;

  case 27:
#line 194 "./ncmpigen.y"
    { type_code = NC_BYTE; }
    break;

  case 28:
#line 195 "./ncmpigen.y"
    { type_code = NC_CHAR; }
    break;

  case 29:
#line 196 "./ncmpigen.y"
    { type_code = NC_SHORT; }
    break;

  case 30:
#line 197 "./ncmpigen.y"
    { type_code = NC_INT; }
    break;

  case 31:
#line 198 "./ncmpigen.y"
    { type_code = NC_FLOAT; }
    break;

  case 32:
#line 199 "./ncmpigen.y"
    { type_code = NC_DOUBLE; }
    break;

  case 35:
#line 205 "./ncmpigen.y"
    {
		    static struct vars dummyvar;

		    dummyvar.name = "dummy";
		    dummyvar.type = NC_DOUBLE;
		    dummyvar.ndims = 0;
		    dummyvar.dims = 0;
		    dummyvar.fill_value.doublev = NC_FILL_DOUBLE;
		    dummyvar.has_data = 0;

		    nvdims = 0;
		    /* make sure variable not re-declared */
		    if ((yyvsp[(1) - (1)])->is_var == 1) {
		       derror( "duplicate variable declaration for %s",
		               (yyvsp[(1) - (1)])->name);
		    }
	            (yyvsp[(1) - (1)])->is_var = 1;
		    (yyvsp[(1) - (1)])->vnum = nvars;
		    /* make sure vars array will hold variables */
		    grow_varray(nvars,  /* must hold nvars+1 vars */
				&vars); /* grow as needed */
		    vars[nvars] = dummyvar; /* to make Purify happy */
		    vars[nvars].name = (char *) emalloc(strlen((yyvsp[(1) - (1)])->name)+1);
		    (void) strcpy(vars[nvars].name, (yyvsp[(1) - (1)])->name);
		    /* name for use in generated Fortran and C variables */
		    vars[nvars].lname = decodify((yyvsp[(1) - (1)])->name);
		    vars[nvars].type = type_code;
		    /* set default fill value.  You can override this with
		     * the variable attribute "_FillValue". */
		    nc_getfill(type_code, &vars[nvars].fill_value);
		    vars[nvars].has_data = 0; /* has no data (yet) */
		   }
    break;

  case 36:
#line 238 "./ncmpigen.y"
    {
		    vars[nvars].ndims = nvdims;
		    nvars++;
		   }
    break;

  case 42:
#line 252 "./ncmpigen.y"
    {
		    if (nvdims >= NC_MAX_VAR_DIMS) {
		       derror("%s has too many dimensions",vars[nvars].name);
		    }
		    if ((yyvsp[(1) - (1)])->is_dim == 1)
		       dimnum = (yyvsp[(1) - (1)])->dnum;
		    else {
		       derror( "%s is not declared as a dimension",
			       (yyvsp[(1) - (1)])->name);
	               dimnum = ndims;
		    }
		    if (rec_dim != -1 && dimnum == rec_dim && nvdims != 0) {
		       derror("unlimited dimension must be first");
		    }
		    grow_iarray(nvdims, /* must hold nvdims+1 ints */
				&vars[nvars].dims); /* grow as needed */
		    vars[nvars].dims[nvdims] = dimnum;
                    nvdims++;
		   }
    break;

  case 43:
#line 273 "./ncmpigen.y"
    {
		       valnum = 0;
		       valtype = NC_UNSPECIFIED;
		       /* get a large block for attributes, realloc later */
		       att_space = emalloc(MAX_NC_ATTSIZE);
		       /* make all kinds of pointers point to it */
		       char_valp = (char *) att_space;
		       byte_valp = (signed char *) att_space;
		       short_valp = (short *) att_space;
		       int_valp = (int *) att_space;
		       float_valp = (float *) att_space;
		       double_valp = (double *) att_space;
		   }
    break;

  case 44:
#line 287 "./ncmpigen.y"
    {
		       {	/* check if duplicate attribute for this var */
			   int i;
			   for(i=0; i<natts; i++) { /* expensive */
			       if(atts[i].var == varnum &&
				  STREQ(atts[i].name,atts[natts].name)) {
				   derror("duplicate attribute %s:%s",
					  vars[varnum].name,atts[natts].name);
			       }
			   }
		       }
		       atts[natts].var = varnum ;
		       atts[natts].type = valtype;
		       atts[natts].len = valnum;
		       /* shrink space down to what was really needed */
		       att_space = erealloc(att_space, valnum*nctypesize(valtype));
		       atts[natts].val = att_space;
		       if (STREQ(atts[natts].name, _FillValue) &&
			   atts[natts].var != NC_GLOBAL) {
			   nc_putfill(atts[natts].type,
				       atts[natts].val,
				       &vars[atts[natts].var].fill_value);
			   if(atts[natts].type != vars[atts[natts].var].type) {
			       derror("variable %s: %s type mismatch",
				      vars[atts[natts].var].name, _FillValue);
			   }
		       }
		       natts++;
		   }
    break;

  case 46:
#line 319 "./ncmpigen.y"
    {
		    varnum = NC_GLOBAL;  /* handle of "global" attribute */
		   }
    break;

  case 47:
#line 325 "./ncmpigen.y"
    { if ((yyvsp[(1) - (1)])->is_var == 1)
		       varnum = (yyvsp[(1) - (1)])->vnum;
		    else {
		      derror("%s not declared as a variable, fatal error",
			     (yyvsp[(1) - (1)])->name);
		      YYABORT;
		      }
		   }
    break;

  case 48:
#line 335 "./ncmpigen.y"
    {
		       /* make sure atts array will hold attributes */
		       grow_aarray(natts,  /* must hold natts+1 atts */
				   &atts); /* grow as needed */
		       atts[natts].name = (char *) emalloc(strlen((yyvsp[(1) - (1)])->name)+1);
		       (void) strcpy(atts[natts].name,(yyvsp[(1) - (1)])->name);
		       /* name for use in generated Fortran and C variables */
		       atts[natts].lname = decodify((yyvsp[(1) - (1)])->name);
		   }
    break;

  case 51:
#line 349 "./ncmpigen.y"
    {
		    if (valtype == NC_UNSPECIFIED)
		      valtype = atype_code;
		    if (valtype != atype_code)
		      derror("values for attribute must be all of same type");
		   }
    break;

  case 52:
#line 358 "./ncmpigen.y"
    {
		       atype_code = NC_CHAR;
		       *char_valp++ = char_val;
		       valnum++;
		   }
    break;

  case 53:
#line 364 "./ncmpigen.y"
    {
		       atype_code = NC_CHAR;
		       {
			   /* don't null-terminate attribute strings */
			   MPI_Offset len = strlen(termstring);
			   if (len == 0) /* need null if that's only value */
			       len = 1;
			   (void)strncpy(char_valp,termstring,len);
			   valnum += len;
			   char_valp += len;
		       }
		   }
    break;

  case 54:
#line 377 "./ncmpigen.y"
    {
		       atype_code = NC_BYTE;
		       *byte_valp++ = byte_val;
		       valnum++;
		   }
    break;

  case 55:
#line 383 "./ncmpigen.y"
    {
		       atype_code = NC_SHORT;
		       *short_valp++ = short_val;
		       valnum++;
		   }
    break;

  case 56:
#line 389 "./ncmpigen.y"
    {
		       atype_code = NC_INT;
		       *int_valp++ = int_val;
		       valnum++;
		   }
    break;

  case 57:
#line 395 "./ncmpigen.y"
    {
		       atype_code = NC_FLOAT;
		       *float_valp++ = float_val;
		       valnum++;
		   }
    break;

  case 58:
#line 401 "./ncmpigen.y"
    {
		       atype_code = NC_DOUBLE;
		       *double_valp++ = double_val;
		       valnum++;
		   }
    break;

  case 64:
#line 417 "./ncmpigen.y"
    {
		       valtype = vars[varnum].type; /* variable type */
		       valnum = 0;	/* values accumulated for variable */
		       vars[varnum].has_data = 1;
		       /* compute dimensions product */
		       var_size = nctypesize(valtype);
		       if (vars[varnum].ndims == 0) { /* scalar */
			   var_len = 1;
		       } else if (vars[varnum].dims[0] == rec_dim) {
			   var_len = 1; /* one record for unlimited vars */
		       } else {
			   var_len = dims[vars[varnum].dims[0]].size;
		       }
		       for(dimnum = 1; dimnum < vars[varnum].ndims; dimnum++)
			 var_len = var_len*dims[vars[varnum].dims[dimnum]].size;
		       /* allocate memory for variable data */
		       if (var_len*var_size != (MPI_Offset)(var_len*var_size)) {
			   derror("variable %s too large for memory",
				  vars[varnum].name);
			   exit(9);
		       }
		       rec_len = var_len;
		       rec_start = malloc ((MPI_Offset)rec_len*var_size);
		       if (rec_start == 0) {
			   derror ("out of memory\n");
			   exit(3);
		       }
		       rec_cur = rec_start;
		       switch (valtype) {
			 case NC_CHAR:
			   char_valp = (char *) rec_start;
			   break;
			 case NC_BYTE:
			   byte_valp = (signed char *) rec_start;
			   break;
			 case NC_SHORT:
			   short_valp = (short *) rec_start;
			   break;
			 case NC_INT:
			   int_valp = (int *) rec_start;
			   break;
			 case NC_FLOAT:
			   float_valp = (float *) rec_start;
			   break;
			 case NC_DOUBLE:
			   double_valp = (double *) rec_start;
			   break;
			 default:
			   derror("Unhandled type %d\n", valtype);
			   break;

		       }
		 }
    break;

  case 65:
#line 471 "./ncmpigen.y"
    {
		       if (valnum < var_len) { /* leftovers */
			   nc_fill(valtype,
				    var_len - valnum,
				    rec_cur,
				    vars[varnum].fill_value);
		       }
		       /* put out var_len values */
		       /* vars[varnum].nrecs = valnum / rec_len; */
		       vars[varnum].nrecs = var_len / rec_len;
		       if (derror_count == 0)
			   put_variable(rec_start);
		       free ((char *) rec_start);
		 }
    break;

  case 68:
#line 490 "./ncmpigen.y"
    {
		       if(valnum >= var_len) {
			   if (vars[varnum].dims[0] != rec_dim) { /* not recvar */
			       derror("too many values for this variable, %d >= %d",
				      valnum, var_len);
			       exit (4);
			   } else { /* a record variable, so grow data
				      container and increment var_len by
				      multiple of record size */
			       ptrdiff_t rec_inc = (char *)rec_cur
				   - (char *)rec_start;
			       var_len = rec_len * (1 + valnum )/rec_len;
			       rec_start = erealloc(rec_start, var_len*var_size);
			       rec_cur = (char *)rec_start + rec_inc;
			       char_valp = (char *) rec_cur;
			       byte_valp = (signed char *) rec_cur;
			       short_valp = (short *) rec_cur;
			       int_valp = (int *) rec_cur;
			       float_valp = (float *) rec_cur;
			       double_valp = (double *) rec_cur;
			   }
		       }
		       not_a_string = 1;
                   }
    break;

  case 69:
#line 515 "./ncmpigen.y"
    {
		       if (not_a_string) {
			   switch (valtype) {
			     case NC_CHAR:
			       rec_cur = (void *) char_valp;
			       break;
			     case NC_BYTE:
			       rec_cur = (void *) byte_valp;
			       break;
			     case NC_SHORT:
			       rec_cur = (void *) short_valp;
			       break;
			     case NC_INT:
			       rec_cur = (void *) int_valp;
			       break;
			     case NC_FLOAT:
			       rec_cur = (void *) float_valp;
			       break;
			     case NC_DOUBLE:
			       rec_cur = (void *) double_valp;
			       break;
			     default: 
			       derror("Unhandled type %d\n", valtype);
			       break;
			   }
		       }
		   }
    break;

  case 70:
#line 545 "./ncmpigen.y"
    {
		       atype_code = NC_CHAR;
		       switch (valtype) {
			 case NC_CHAR:
			   *char_valp++ = char_val;
			   break;
			 case NC_BYTE:
			   *byte_valp++ = char_val;
			   break;
			 case NC_SHORT:
			   *short_valp++ = char_val;
			   break;
			 case NC_INT:
			   *int_valp++ = char_val;
			   break;
			 case NC_FLOAT:
			   *float_valp++ = char_val;
			   break;
			 case NC_DOUBLE:
			   *double_valp++ = char_val;
			   break; 
			 default: 
			   derror("Unhandled type %d\n", valtype);
			   break;
		       }
		       valnum++;
		   }
    break;

  case 71:
#line 573 "./ncmpigen.y"
    {
		       not_a_string = 0;
		       atype_code = NC_CHAR;
		       {
			   MPI_Offset len = strlen(termstring);

			   if(valnum + len > var_len) {
			       if (vars[varnum].dims[0] != rec_dim) {
				   derror("too many values for this variable, %d>%d", 
					  valnum+len, var_len);
				   exit (5);
			       } else {/* a record variable so grow it */
				   ptrdiff_t rec_inc = (char *)rec_cur
				       - (char *)rec_start;
				   var_len += rec_len * (len + valnum - var_len)/rec_len;
				   rec_start = erealloc(rec_start, var_len*var_size);
				   rec_cur = (char *)rec_start + rec_inc;
				   char_valp = (char *) rec_cur;
			       }
			   }
			   switch (valtype) {
			     case NC_CHAR:
			       {
				   int ld;
				   MPI_Offset i, sl;
				   (void)strncpy(char_valp,termstring,len);
				   
				   ld = vars[varnum].ndims-1;
				   if (ld > 0) {/* null-fill to size of last dim */
				       sl = dims[vars[varnum].dims[ld]].size;
				       for (i =len;i<sl;i++)
					   char_valp[i] = '\0';
				       if (sl < len)
					   sl = len;
				       valnum += sl;
				       char_valp += sl;
				   } else { /* scalar or 1D strings */
				       valnum += len;
				       char_valp += len;
				   }
				   rec_cur = (void *) char_valp;
			       }
			       break;
			     case NC_BYTE:
			     case NC_SHORT:
			     case NC_INT:
			     case NC_FLOAT:
			     case NC_DOUBLE:
			       derror("string value invalid for %s variable",
				      nctype(valtype));
			       break;
			     default: 
			       derror("Unhandled type %d\n", valtype);
			       break;
			   }
		       }
		   }
    break;

  case 72:
#line 631 "./ncmpigen.y"
    {
		       atype_code = NC_BYTE;
		       switch (valtype) {
			 case NC_CHAR:
			   *char_valp++ = byte_val;
			   break;
			 case NC_BYTE:
			   *byte_valp++ = byte_val;
			   break;
			 case NC_SHORT:
			   *short_valp++ = byte_val;
			   break;
			 case NC_INT:
			   *int_valp++ = byte_val;
			   break;
			 case NC_FLOAT:
			   *float_valp++ = byte_val;
			   break;
			 case NC_DOUBLE:
			   *double_valp++ = byte_val;
			   break;
			 default: 
			   derror("Unhandled type %d\n", valtype);
			   break;
		       }
		       valnum++;
		   }
    break;

  case 73:
#line 659 "./ncmpigen.y"
    {
		       atype_code = NC_SHORT;
		       switch (valtype) {
			 case NC_CHAR:
			   *char_valp++ = short_val;
			   break;
			 case NC_BYTE:
			   *byte_valp++ = short_val;
			   break;
			 case NC_SHORT:
			   *short_valp++ = short_val;
			   break;
			 case NC_INT:
			   *int_valp++ = short_val;
			   break;
			 case NC_FLOAT:
			   *float_valp++ = short_val;
			   break;
			 case NC_DOUBLE:
			   *double_valp++ = short_val;
			   break;
			 default: 
			   derror("Unhandled type %d\n", valtype);
			   break;
		       }
		       valnum++;
		   }
    break;

  case 74:
#line 687 "./ncmpigen.y"
    {
		       atype_code = NC_INT;
		       switch (valtype) {
			 case NC_CHAR:
			   *char_valp++ = int_val;
			   break;
			 case NC_BYTE:
			   *byte_valp++ = int_val;
			   break;
			 case NC_SHORT:
			   *short_valp++ = int_val;
			   break;
			 case NC_INT:
			   *int_valp++ = int_val;
			   break;
			 case NC_FLOAT:
			   *float_valp++ = int_val;
			   break;
			 case NC_DOUBLE:
			   *double_valp++ = int_val;
			   break;
			 default: 
			   derror("Unhandled type %d\n", valtype);
			   break;
		       }
		       valnum++;
		   }
    break;

  case 75:
#line 715 "./ncmpigen.y"
    {
		       atype_code = NC_FLOAT;
		       switch (valtype) {
			 case NC_CHAR:
			   *char_valp++ = float_val;
			   break;
			 case NC_BYTE:
			   *byte_valp++ = float_val;
			   break;
			 case NC_SHORT:
			   *short_valp++ = float_val;
			   break;
			 case NC_INT:
			   *int_valp++ = float_val;
			   break;
			 case NC_FLOAT:
			   *float_valp++ = float_val;
			   break;
			 case NC_DOUBLE:
			   *double_valp++ = float_val;
			   break;
			 default: 
			   derror("Unhandled type %d\n", valtype);
			   break;
		       }
		       valnum++;
		   }
    break;

  case 76:
#line 743 "./ncmpigen.y"
    {
		       atype_code = NC_DOUBLE;
		       switch (valtype) {
			 case NC_CHAR:
			   *char_valp++ = double_val;
			   break;
			 case NC_BYTE:
			   *byte_valp++ = double_val;
			   break;
			 case NC_SHORT:
			   *short_valp++ = double_val;
			   break;
			 case NC_INT:
			   *int_valp++ = double_val;
			   break;
			 case NC_FLOAT:
			   if (double_val == NC_FILL_DOUBLE)
			     *float_valp++ = NC_FILL_FLOAT;
			   else
			     *float_valp++ = double_val;
			   break;
			 case NC_DOUBLE:
			   *double_valp++ = double_val;
			   break;
			 default: 
			   derror("Unhandled type %d\n", valtype);
			   break;
		       }
		       valnum++;
		   }
    break;

  case 77:
#line 774 "./ncmpigen.y"
    {
		       /* store fill_value */
		       switch (valtype) {
		       case NC_CHAR:
			   nc_fill(valtype, 1, (void *)char_valp++,
				   vars[varnum].fill_value);
			   break;
		       case NC_BYTE:
			   nc_fill(valtype, 1, (void *)byte_valp++,
				   vars[varnum].fill_value);
			   break;
		       case NC_SHORT:
			   nc_fill(valtype, 1, (void *)short_valp++,
				   vars[varnum].fill_value);
			   break;
		       case NC_INT:
			   nc_fill(valtype, 1, (void *)int_valp++,
				   vars[varnum].fill_value);
			   break;
		       case NC_FLOAT:
			   nc_fill(valtype, 1, (void *)float_valp++,
				   vars[varnum].fill_value);
			   break;
		       case NC_DOUBLE:
			   nc_fill(valtype, 1, (void *)double_valp++,
				   vars[varnum].fill_value);
			   break;
		       default: 
			   derror("Unhandled type %d\n", valtype);
			   break;
		       }
		       valnum++;
		   }
    break;


/* Line 1267 of yacc.c.  */
#line 2272 "y.tab.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


#line 811 "./ncmpigen.y"


/* PROGRAMS */

#ifdef vms
void
#else
int
#endif
yyerror(	/* called for yacc syntax error */
     char *s)
{
	derror(s);
#ifndef vms
	return -1;
#endif
}

/* undefine yywrap macro, in case we are using bison instead of yacc */
#ifdef ncmpiwrap
#undef ncmpiwrap
#endif

int
ncmpiwrap(void)			/* returns 1 on EOF if no more input */
{
    return  1;
}


/* Symbol table operations for ncmpigen tool */

YYSTYPE lookup(       /* find sname in symbol table (linear search) */
	const char *sname)
{
    YYSTYPE sp;
    for (sp = symlist; sp != (YYSTYPE) 0; sp = sp -> next)
	if (STREQ(sp -> name, sname)) {
	    return sp;
	}
    return 0;			/* 0 ==> not found */
}

YYSTYPE install(  /* install sname in symbol table */
	const char *sname)
{
    YYSTYPE sp;

    sp = (YYSTYPE) emalloc (sizeof (struct Symbol));
    sp -> name = (char *) emalloc (strlen (sname) + 1);/* +1 for '\0' */
    (void) strcpy (sp -> name, sname);
    sp -> next = symlist;	/* put at front of list */
    sp -> is_dim = 0;
    sp -> is_var = 0;
    sp -> is_att = 0;
    symlist = sp;
    return sp;
}

void
clearout(void)	/* reset symbol table to empty */
{
    YYSTYPE sp, tp;
    for (sp = symlist; sp != (YYSTYPE) 0;) {
	tp = sp -> next;
	free (sp -> name);
	free ((char *) sp);
	sp = tp;
    }
    symlist = 0;
}

/* get lexical input routine generated by lex  */
#include "ncmpigenyy.c"

