## Utility Programs

### ncmpidiff
<ul>
  <li> An MPI program runs in parallel to compare the contents of the two files and
  reports the first difference to the standard output.
  </li>
  <li> <details>
  <summary>Manual page (click to expand)</summary>
  
```
 ncmpidiff(1)                   PnetCDF utilities                  ncmpidiff(1)
 
 NAME
        ncmpidiff - compares two netCDF files in parallel
 
 SYNOPSIS
        mpiexec   -n  np  ncmpidiff  [-b]  [-q]  [-h]  [-v  var1,...,varn]  [-t
               diff,ratio] file1 file2
 
 DESCRIPTION
        ncmpidiff runs in parallel on np number of MPI processes to compare the
        contents of the two files and reports the first difference to the stan-
        dard output.
 
        For variables and attributes, it reports the array indices of the first
        element  found  different when option -t is not used. When option -t is
        used, it reports the element with the largest difference that fails  to
        meet the tolerance requirements.
 
        If  neither argument -v nor -h is given besides the two file names, the
        entire files are compared.
 
        When comparing two files entirely, the difference between ncmpidiff and
        the  Unix  command  diff is that ncmpidiff skips the gaps between vari-
        ables. The gaps may occur when the alignment feature is used to  create
        a  new  file.  This alignment allows to allocate a larger space for the
        file header and align the starting file offsets of fixed-size variables
        (see  API ncmpi__enddef and PnetCDF hints). Oftentimes, the contents of
        gaps are non-zero arbitrary bytes. Thus, two netCDF files (of  same  or
        different  sizes)  can  be  reported  identical by ncmpidiff but not by
        diff.
 
 OPTIONS
        -b     Verbose mode - print results (same or different) for all  compo-
               nents (file, header, or variables) in comparison
 
        -q     Quiet mode - print nothing on the command-line output. This also
               disables verbose mode. When in quiet mode,  users  should  check
               exit status. See below in "EXIT STATUS".
 
        -h     Compare file header only
 
        -v var1,...,varn
               Compare  only  the  given  list of variables (names separated by
               comma without space).
 
        -t diff,ratio
               Compare variables element-wisely with tolerance (diff and  ratio
               separated  by  comma without space).  diff is the absolute value
               of element-wise difference of  two variables with the same  name
               but  stored  in the two input files.  ratio is the relative ele-
               ment-wise difference ratio defined as  |x  -  y|/max(|x|,  |y|),
               where  x  is  an array element from a variable in the first file
               and y is the corresponding array element of the same variable in
               the  second  file.  |x| represents the absolute value of x. Note
               when this option is used, the output reports only the first  ar-
               ray element that fails to meet both tolerance requirements.
 
 EXIT STATUS
        An  exit  status of 0 means no differences were found, and 1 means some
        differences were found.  Note on VMS-based system, the exit status val-
        ues are reversed.
 
 SEE ALSO
        ncmpidump(1), cdfdiff(1), diff(1), pnetcdf(3)
 
 DATE
        February 21, 2022
```
</details></li>
</ul>

### cdfdiff
<ul>
  <li>A sequential version of `ncmpidiff`, compares the contents of the two classic
  netCDF files and reports the first difference found to the standard output.
  The classic file formats include CDF-1, CDF-2, and CDF-5.
  </li>
  <li> <details>
  <summary>Manual page (click to expand)</summary>
  
```
 cdfdiff(1)                     PnetCDF utilities                    cdfdiff(1)
 
 NAME
        cdfdiff - compares two classic netCDF files in parallel
 
 SYNOPSIS
        cdfdiff [-b] [-q] [-h] [-v var1,...,varn] [-t diff,ratio] file1 file2
 
 DESCRIPTION
        cdfdiff,  a  sequential  version of ncmpidiff, compares the contents of
        the two classic netCDF files and reports the first difference found  to
        the standard output. The classic file formats include CDF-1, CDF-2, and
        CDF-5.
 
        For variables and attributes, it reports the array indices of the first
        element  found  different when option -t is not used. When option -t is
        used, it reports the element with the largest difference that fails  to
        meet the tolerance requirements.
 
        If  neither argument -v nor -h is given besides the two file names, the
        entire files are compared.
 
        When comparing two files entirely, the difference between  cdfdiff  and
        the Unix command diff is that cdfdiff skips the gaps between variables.
        The gaps may occur when the alignment feature is used to create  a  new
        file.  This  alignment  allows  to allocate a larger space for the file
        header and align the starting file offsets of fixed-size variables (see
        API  ncmpi__enddef and PnetCDF hints). Oftentimes, the contents of gaps
        are non-zero arbitrary bytes. Thus, two netCDF files (of same  or  dif-
        ferent sizes) can be reported identical by cdfdiff but not by diff.
 
 OPTIONS
        -b     Verbose  mode - print results (same or different) for all compo-
               nents (file, header, or variables) in comparison
 
        -q     Quiet mode - print nothing on the command-line output. This also
               disables  verbose  mode.  When in quiet mode, users should check
               exit status. See below in "EXIT STATUS".
 
        -h     Compare file header only
 
        -v var1,...,varn
               Compare only the given list of  variables  (names  separated  by
               comma without space).
 
        -t diff,ratio
               Compare  variables element-wisely with tolerance (diff and ratio
               separated by comma without space).  diff is the  absolute  value
               of  element-wise difference of  two variables with the same name
               but stored in the two input files.  ratio is the  relative  ele-
               ment-wise  difference  ratio  defined  as |x - y|/max(|x|, |y|),
               where x is an array element from a variable in  the  first  file
               and y is the corresponding array element of the same variable in
               the second file. |x| represents the absolute value  of  x.  Note
               when  this option is used, the output reports only the first ar-
               ray element that fails to meet both tolerance requirements.
 
 EXIT STATUS
        An exit status of 0 means no differences were found, and 1  means  some
        differences were found.  Note on VMS-based system, the exit status val-
        ues are reversed.
 
 SEE ALSO
        ncmpidiff(1), ncmpidump(1), diff(1), pnetcdf(3)
 
 DATE
        February 21, 2022
```
</details></li>
</ul>
  

### ncmpidump
<ul>
  <li>A utility program to generate an ASCII representation of a specified netCDF file on
  standard output.
  </li>
  <li> <details>
  <summary>Manual page (click to expand)</summary>

```
 ncmpidump(1)                   PnetCDF utilities                  ncmpidump(1)
 
 NAME
        ncmpidump - Convert netCDF files to ASCII form (CDL)
 
 SYNOPSIS
        ncmpidump  [-ch]  [-v var1,...]  [-b lang] [-f lang] [-l len] [-n name]
               [-p f_digits[,d_digits]] file
 
        ncmpidump -k file
 
 DESCRIPTION
        ncmpidump generates an ASCII representation of a specified netCDF  file
        on  standard  output.  The ASCII representation is in a form called CDL
        (``network Common Data form Language'') that can be viewed, edited,  or
        serve  as  input to ncmpigen.  ncmpigen is a companion program that can
        generate a binary netCDF file from a  CDL  file.   Hence  ncmpigen  and
        ncmpidump  can be used as inverses to transform the data representation
        between binary and ASCII representations.  See ncmpigen for a  descrip-
        tion of CDL and netCDF representations.
 
        ncmpidump  may  also  be  used to determine what kind of netCDF file is
        used (which variant of the netCDF file format) with the -k option.
 
        ncmpidump defines a default format used for each type of  netCDF  data,
        but  this  can  be  changed  if a `C_format' attribute is defined for a
        netCDF variable.  In this case, ncmpidump will use the  `C_format'  at-
        tribute  to format each value.  For example, if floating-point data for
        the netCDF variable `Z' is known to be accurate to only three  signifi-
        cant digits, it would be appropriate to use the variable attribute
 
               Z:C_format = "%.3g"
 
        ncmpidump  may  also be used as a simple browser for netCDF data files,
        to display the dimension names and sizes; variable  names,  types,  and
        shapes;  attribute names and values; and optionally, the values of data
        for all variables or selected variables in a netCDF file.
 
        ncmpidump uses `_' to represent data  values  that  are  equal  to  the
        `_FillValue'  attribute for a variable, intended to represent data that
        has not yet been written.  If a variable has no `_FillValue' attribute,
        the default fill value for the variable type is used if the variable is
        not of byte type.
 
 OPTIONS
        -c     Show the values of coordinate variables (variables that are also
               dimensions) as well as the declarations of all dimensions, vari-
               ables, and attribute  values.   Data  values  of  non-coordinate
               variables  are  not  included  in  the output.  This is the most
               suitable option to use for a brief look  at  the  structure  and
               contents of a netCDF file.
 
        -h     Show only the header information in the output, that is the dec-
               larations of dimensions, variables, and attributes but  no  data
               values  for any variables.  The output is identical to using the
               -c option except that the values of coordinate variables are not
               included.  (At most one of -c or -h options may be present.)
 
        -v var1,...,varn
               The output will include data values for the specified variables,
               in addition to the declarations of  all  dimensions,  variables,
               and attributes.  One or more variables must be specified by name
               in the comma-delimited list following  this  option.   The  list
               must  be  a single argument to the command, hence cannot contain
               blanks or other white space  characters.   The  named  variables
               must  be valid netCDF variables in the input-file.  The default,
               without this option and in the absence of the -c or -h  options,
               is to include data values for all variables in the output.
 
        -b lang
               A  brief annotation in the form of a CDL comment (text beginning
               with the characters ``//'') will be included in the data section
               of the output for each `row' of data, to help identify data val-
               ues for multidimensional variables.  If lang begins with `C'  or
               `c',  then  C  language conventions will be used (zero-based in-
               dices, last dimension varying fastest).  If lang begins with `F'
               or  `f',  then  Fortran  language conventions will be used (one-
               based indices, first  dimension  varying  fastest).   In  either
               case, the data will be presented in the same order; only the an-
               notations will differ.   This  option  is  useful  for  browsing
               through large volumes of multidimensional data.
 
        -f lang
               Full  annotations in the form of trailing CDL comments (text be-
               ginning with the characters ``//'') for every data value (except
               individual  characters  in character arrays) will be included in
               the data section.  If lang begins with `C' or `c', then  C  lan-
               guage  conventions will be used (zero-based indices, last dimen-
               sion varying fastest).  If lang begins with  `F'  or  `f',  then
               Fortran  language  conventions  will be used (one-based indices,
               first dimension varying fastest).  In either case, the data will
               be  presented  in the same order; only the annotations will dif-
               fer.  This option may be useful for piping data into other  fil-
               ters,  since  each  data value appears on a separate line, fully
               identified.
 
        -l len Changes the default maximum line length (80) used in  formatting
               lists of non-character data values.
 
        -n name
               CDL  requires  a name for a netCDF data set, for use by ncmpigen
               -b in generating  a  default  netCDF  file  name.   By  default,
               ncmpidump  constructs  this  name from the last component of the
               pathname of the input netCDF file by stripping off any extension
               it  has.   Use  the  -n option to specify a different name.  Al-
               though the output file name used by ncmpigen -b  can  be  speci-
               fied,  it  may be wise to have ncmpidump change the default name
               to avoid inadvertantly overwriting a valuable netCDF  file  when
               using  ncmpidump,  editing  the  resulting  CDL  file, and using
               ncmpigen -b to generate a new netCDF file from  the  edited  CDL
               file.
 
        -p float_digits[,double_digits]
               Specifies  default  precision  (number of significant digits) to
               use in displaying floating-point or double precision data values
               for  attributes  and  variables.  If specified, this value over-
               rides the value of the `C_format'  attribute  for  any  variable
               that  has  such  an attribute.  Floating-point data will be dis-
               played with float_digits significant digits.   If  double_digits
               is  also  specified,  double-precision  values will be displayed
               with that many significant digits.  In the  absence  of  any  -p
               specifications,  floating-point  and  double-precision  data are
               displayed with 7 and 15 significant  digits  respectively.   CDL
               files  can  be  made  smaller if less precision is required.  If
               both floating-point and double-precision precisions  are  speci-
               fied,  the  two  values  must  appear  separated  by a comma (no
               blanks) as a single argument to the command.  If you really want
               every  last bit of precision from the netCDF file represented in
               the CDL file for all possible floating-point  values,  you  will
               have  to  specify  this with -p 9,17 (according to Theorem 15 of
               the paper listed under REFERENCES).
 
        -k     Reports the kind of netCDF  file:  classic,  64-bit  offset,  or
               64-bit data.  Before netCDF version 3.6, there was only one kind
               of netCDF file, designated as `classic'  (also  know  as  format
               variant  1  or  CDF-1).   Large  file support introduced another
               variant of the format, designated as `64-bit offset'  (known  as
               format  variant  2 or CDF-2).  Large data support introduced an-
               other variant of the format, designated as `64-bit data'  (known
               as format variant 5 or CDF-5).
 
 EXAMPLES
        Look at the structure of the data in the netCDF file `foo.nc':
 
               ncmpidump -c foo.nc
 
        Produce  an  annotated  CDL  version  of  the structure and data in the
        netCDF file `foo.nc', using C-style indexing for the annotations:
 
               ncmpidump -b c foo.nc > foo.cdl
 
        Output data for only the variables `uwind' and `vwind' from the  netCDF
        file `foo.nc', and show the floating-point data with only three signif-
        icant digits of precision:
 
               ncmpidump -v uwind,vwind -p 3 foo.nc
 
        Produce a fully-annotated (one data value per line) listing of the data
        for  the  variable  `omega', using Fortran conventions for indices, and
        changing the netCDF dataset name in the resulting CDL file to `omega':
 
               ncmpidump -v omega -f fortran -n omega foo.nc > Z.cdl
 
 REFERENCES
         What Every Computer Scientist should Know About Floating-Point  Arith-
        metic, D.  Goldberg, ACM Computing Surveys, Vol. 23, No. 1, March 1991,
        pp. 5-48.
 
 SEE ALSO
        ncmpigen(1), pnetcdf(3)
 
 DATE
        February 21, 2022
 
 BUGS
        Character arrays that contain a null-byte are treated like  C  strings,
        so no characters after the null byte appear in the output.
 
        Multidimensional  character  string  arrays are not handled well, since
        the CDL syntax for breaking a long character string into several short-
        er lines is weak.
 
        There  should  be a way to specify that the data should be displayed in
        `record' order, that is with the all the values for `record'  variables
        together that have the same value of the record dimension. 
```
</details></li>
</ul>

### ncmpigen
<ul>
  <li>A utility program to generate either a netCDF file, or C or Fortran source code to
  create a netCDF file.
  </li>
  <li> <details>
  <summary>Manual page (click to expand)</summary>

```
 ncmpigen(1)                    PnetCDF utilities                   ncmpigen(1)
 
 NAME
        ncmpigen  -  From  a CDL file generate a netCDF file, a C program, or a
        Fortran program
 
 SYNOPSIS
        ncmpigen [-b] [-c] [-f]  [-n]  [-o  netcdf_filename]  [-v  file_format]
               input_file
 
 DESCRIPTION
        ncmpigen generates either a netCDF file, or C or Fortran source code to
        create a netCDF file.  The input to ncmpigen  is  a  description  of  a
        netCDF  file in a small language known as CDL (network Common Data form
        Language), described below.  If no options are  specified  in  invoking
        ncmpigen,  it merely checks the syntax of the input CDL file, producing
        error messages for any violations of CDL syntax.  Other options can  be
        used  to  create the corresponding netCDF file, to generate a C program
        that uses the netCDF C interface to create the netCDF file, or to  gen-
        erate  a Fortran program that uses the netCDF Fortran interface to cre-
        ate the same netCDF file.
 
        ncmpigen may be used with the companion program  ncmpidump  to  perform
        some  simple  operations on netCDF files.  For example, to rename a di-
        mension in a netCDF file, use ncmpidump to get a  CDL  version  of  the
        netCDF  file,  edit  the CDL file to change the name of the dimensions,
        and use ncmpigen to generate the corresponding  netCDF  file  from  the
        edited CDL file.
 
 OPTIONS
        -b     Create  a  (binary)  netCDF file.  If the -o option is absent, a
               default file name will  be  constructed  from  the  netCDF  name
               (specified  after  the netcdf keyword in the input) by appending
               the `.nc' extension.  If a file already exists with  the  speci-
               fied name, it will be overwritten.
 
        -c     Generate  C  source code that will create a netCDF file matching
               the netCDF specification.  The C source code is written to stan-
               dard output.
 
        -f     Generate  Fortran  source  code  that  will create a netCDF file
               matching the netCDF specification.  The Fortran source  code  is
               written to standard output.
 
        -o netcdf_file
               Name  for  the  binary  netCDF  file created.  If this option is
               specified, it implies the "-b" option.  (This option  is  neces-
               sary because netCDF files cannot be written directly to standard
               output, since standard output is not seekable.)
 
        -n     Like -b option, except creates netCDF  file  with  the  obsolete
               `.cdf'  extension instead of the `.nc' extension, in the absence
               of an output filename specified by the -o option.   This  option
               is only supported for backward compatibility.
 
        -v file_format
               File  format of the output netCDF file. The value of file_format
               can be: 1 or classic for CDF-1 format.  2  or  64-bit-offset  is
               CDF-2.   5  or  64-bit-variable for CDF-5.  The default (if this
               option is not given) is CDF-1, the classic format.
 
 EXAMPLES
        Check the syntax of the CDL file `foo.cdl':
 
               ncmpigen foo.cdl
 
        From the CDL file `foo.cdl', generate an equivalent binary netCDF  file
        named `x.nc':
 
               ncmpigen -o x.nc foo.cdl
 
        From the CDL file `foo.cdl', generate a C program containing the netCDF
        function invocations necessary to create an  equivalent  binary  netCDF
        file named `x.nc':
 
               ncmpigen -c -o x.nc foo.cdl
 
 USAGE
    CDL Syntax Summary
        Below is an example of CDL syntax, describing a netCDF file with sever-
        al named dimensions (lat, lon, and time), variables (Z, t, p, rh,  lat,
        lon,  time), variable attributes (units, long_name, valid_range, _Fill-
        Value), and some data.  CDL keywords are in boldface.  (This example is
        intended  to  illustrate  the syntax; a real CDL file would have a more
        complete set of attributes so that the data would  be  more  completely
        self-describing.)
 
               netcdf foo {  // an example netCDF specification in CDL
 
               dimensions:
                    lat = 10, lon = 5, time = unlimited ;
 
               variables:
                    long    lat(lat), lon(lon), time(time);
                    float   Z(time,lat,lon), t(time,lat,lon);
                    double  p(time,lat,lon);
                    long    rh(time,lat,lon);
 
                    // variable attributes
                    lat:long_name = "latitude";
                    lat:units = "degrees_north";
                    lon:long_name = "longitude";
                    lon:units = "degrees_east";
                    time:units = "seconds since 1992-1-1 00:00:00";
                    Z:units = "geopotential meters";
                    Z:valid_range = 0., 5000.;
                    p:_FillValue = -9999.;
                    rh:_FillValue = -1;
 
               data:
                    lat   = 0, 10, 20, 30, 40, 50, 60, 70, 80, 90;
                    lon   = -140, -118, -96, -84, -52;
               }
 
        All  CDL  statements  are terminated by a semicolon.  Spaces, tabs, and
        newlines can be used freely for readability.  Comments may  follow  the
        characters `//' on any line.
 
        A  CDL  description consists of three optional parts: dimensions, vari-
        ables, and data, beginning with the  keyword  dimensions:,  variables:,
        and  data, respectively.  The variable part may contain variable decla-
        rations and attribute assignments.
 
        A netCDF dimension is used to define the shape of one or  more  of  the
        multidimensional  variables contained in the netCDF file.  A netCDF di-
        mension has a name and a size.  At most one dimension in a netCDF  file
        can  have  the unlimited size, which means a variable using this dimen-
        sion can grow to any length (like a record number in a file).
 
        A variable represents a multidimensional array of values  of  the  same
        type.  A variable has a name, a data type, and a shape described by its
        list of dimensions.  Each variable may also have associated  attributes
        (see  below) as well as data values.  The name, data type, and shape of
        a variable are specified by its declaration in the variable section  of
        a  CDL  description.  A variable may have the same name as a dimension;
        by convention such a variable is one-dimensional and  contains  coordi-
        nates  of the dimension it names.  Dimensions need not have correspond-
        ing variables.
 
        A netCDF attribute contains information  about  a  netCDF  variable  or
        about  the  whole  netCDF dataset.  Attributes are used to specify such
        properties as units, special values, maximum and minimum valid  values,
        scaling  factors,  offsets,  and  parameters.  Attribute information is
        represented by single values or arrays of values.  For example, "units"
        is an attribute represented by a character array such as "celsius".  An
        attribute has an associated variable, a name, a data  type,  a  length,
        and  a value.  In contrast to variables that are intended for data, at-
        tributes are intended for metadata (data about data).
 
        In CDL, an attribute is designated by a variable  and  attribute  name,
        separated by `:'.  It is possible to assign global attributes not asso-
        ciated with any variable to the netCDF as a whole by using  `:'  before
        the  attribute  name.   The data type of an attribute in CDL is derived
        from the type of the value assigned to it.  The length of an  attribute
        is  the  number of data values assigned to it, or the number of charac-
        ters in the character string assigned to it.  Multiple values  are  as-
        signed  to  non-character attributes by separating the values with com-
        mas.  All values assigned to an attribute must be of the same type.
 
        The names for CDL dimensions, variables, and attributes must begin with
        an  alphabetic  character  or `_', and subsequent characters may be al-
        phanumeric or `_' or `-'.
 
        The optional data section of a CDL specification is where netCDF  vari-
        ables may be initialized.  The syntax of an initialization is simple: a
        variable name, an equals sign, and a comma-delimited list of  constants
        (possibly  separated  by  spaces,  tabs and newlines) terminated with a
        semicolon.  For multi-dimensional arrays,  the  last  dimension  varies
        fastest.  Thus row-order rather than column order is used for matrices.
        If fewer values are supplied than are needed to fill a variable, it  is
        extended with a type-dependent `fill value', which can be overridden by
        supplying a value for a distinguished variable attribute named  `_Fill-
        Value'.   The types of constants need not match the type declared for a
        variable; coercions are done to convert integers to floating point, for
        example.   The constant `_' can be used to designate the fill value for
        a variable.
 
    Primitive Data Types
               char characters
               byte 8-bit data
               short     16-bit signed integers
               long 32-bit signed integers
               int  (synonymous with long)
               float     IEEE single precision floating point (32 bits)
               real (synonymous with float)
               double    IEEE double precision floating point (64 bits)
 
        Except for the added data-type byte and the lack of unsigned, CDL  sup-
        ports  the same primitive data types as C.  The names for the primitive
        data types are reserved words in CDL, so the names of variables, dimen-
        sions,  and  attributes  must not be type names.  In declarations, type
        names may be specified in either upper or lower case.
 
        Bytes differ from characters in that they are intended to hold  a  full
        eight  bits  of data, and the zero byte has no special significance, as
        it does for character data.  ncmpigen  converts  byte  declarations  to
        char declarations in the output C code and to the nonstandard BYTE dec-
        laration in output Fortran code.
 
        Shorts can hold values between -32768  and  32767.   ncmpigen  converts
        short  declarations  to  short declarations in the output C code and to
        the nonstandard INTEGER*2 declaration in output Fortran code.
 
        Longs can hold values between  -2147483648  and  2147483647.   ncmpigen
        converts  long  declarations  to long declarations in the output C code
        and to INTEGER declarations in output Fortran code.   int  and  integer
        are  accepted as synonyms for long in CDL declarations.  Now that there
        are platforms with 64-bit representations for C longs, it may be better
        to use the int synonym to avoid confusion.
 
        Floats  can hold values between about -3.4+38 and 3.4+38.  Their exter-
        nal representation is as 32-bit IEEE normalized single-precision float-
        ing  point numbers.  ncmpigen converts float declarations to float dec-
        larations in the output C code and to REAL declarations in output  For-
        tran  code.   real  is  accepted as a synonym for float in CDL declara-
        tions.
 
        Doubles can hold values between about -1.7+308 and 1.7+308.  Their  ex-
        ternal representation is as 64-bit IEEE standard normalized double-pre-
        cision floating point numbers.  ncmpigen converts  double  declarations
        to  double  declarations  in  the output C code and to DOUBLE PRECISION
        declarations in output Fortran code.
 
    CDL Constants
        Constants assigned to attributes or variables may be of any of the  ba-
        sic netCDF types.  The syntax for constants is similar to C syntax, ex-
        cept that type suffixes must be appended to shorts and floats  to  dis-
        tinguish them from longs and doubles.
 
        A  byte constant is represented by a single character or multiple char-
        acter escape sequence enclosed in single quotes.  For example,
                'a'           // ASCII `a'
                '\0'          // a zero byte
                '\n'          // ASCII newline character
                '\33'         // ASCII escape character (33 octal)
                '\x2b'        // ASCII plus (2b hex)
                '\377'        // 377 octal = 255 decimal, non-ASCII
 
        Character constants are enclosed in double quotes.  A  character  array
        may  be represented as a string enclosed in double quotes.  The usual C
        string escape conventions are honored.  For example
               "a"             // ASCII `a'
               "Two\nlines\n"  // a 10-character string with two embedded newlin es
               "a bell:\007"   // a string containing an ASCII bell
        Note that the netCDF character array "a" would  fit  in  a  one-element
        variable,  since  no terminating NULL character is assumed.  However, a
        zero byte in a character array is interpreted as the end of the signif-
        icant  characters by the ncmpidump program, following the C convention.
        Therefore, a NULL byte should not be embedded in a character string un-
        less  at  the  end: use the byte data type instead for byte arrays that
        contain the zero byte.  NetCDF and CDL have no string  type,  but  only
        fixed-length character arrays, which may be multi-dimensional.
 
        short  integer  constants  are  intended for representing 16-bit signed
        quantities.  The form of a short constant is an integer  constant  with
        an `s' or `S' appended.  If a short constant begins with `0', it is in-
        terpreted as octal, except that if it begins with `0x',  it  is  inter-
        preted as a hexadecimal constant.  For example:
               -2s      // a short -2
               0123s    // octal
               0x7ffs   //hexadecimal
 
        Long  integer  constants  are  intended  for representing 32-bit signed
        quantities.  The form of a long constant is an  ordinary  integer  con-
        stant,  although it is acceptable to append an optional `l' or `L'.  If
        a long constant begins with `0', it is  interpreted  as  octal,  except
        that  if  it  begins with `0x', it is interpreted as a hexadecimal con-
        stant.  Examples of valid long constants include:
               -2
               1234567890L
               0123i         // octal
               0x7ff         // hexadecimal
 
        Floating point constants of type float are appropriate for representing
        floating  point  data with about seven significant digits of precision.
        The form of a float constant is the same as a C floating point constant
        with an `f' or `F' appended.  For example the following are all accept-
        able float constants:
                -2.0f
                3.14159265358979f    // will be truncated to less precision
                1.f
                .1f
 
        Floating point constants of type double are appropriate for  represent-
        ing floating point data with about sixteen significant digits of preci-
        sion.  The form of a double constant is the same as a C floating  point
        constant.   An  optional  `d'  or `D' may be appended.  For example the
        following are all acceptable double constants:
               -2.0
               3.141592653589793
               1.0e-20
               1.d
 
 DATE
        February 21, 2022
 
 BUGS
        The programs generated by ncmpigen when using the -c or -f use initial-
        ization statements to store data in variables, and will fail to produce
        compilable programs if you try to use them for  large  datasets,  since
        the  resulting  statements may exceed the line length or number of con-
        tinuation statements permitted by the compiler.
 
        The CDL syntax makes it easy to assign what  looks  like  an  array  of
        variable-length strings to a netCDF variable, but the strings will sim-
        ply be concatenated into a single array  of  characters,  since  netCDF
        cannot  represent  an  array  of  variable-length strings in one netCDF
        variable.
 
        NetCDF and CDL do not yet support a type corresponding to a 64-bit  in-
        teger.
```
</details></li>
</ul>

### ncoffsets
<ul>
  <li>A utility program to print the file offsets information of variables defined in a
  given netCDF file.
  </li>
  <li> <details>
  <summary>Manual page (click to expand)</summary>

```
 ncoffsets(1)                   PnetCDF utilities                  ncoffsets(1)
 
 NAME
        ncoffsets - print the starting/ending file offsets for netCDF variables
 
 SYNOPSIS
        ncoffsets [-h] | [-x] | [-sgr] [-v var1[,...]]  file
 
 DESCRIPTION
        ncoffsets prints the file offsets information of variables defined in a
        given  netCDF file. The ending offsets reported is an exclusive offset,
        i.e.  1 byte more than the last byte occupied by the variable. In other
        words, the ending offset is equal to the sum of starting offset and the
        variable size.  For record variables, only the offsets of first  record
        are printed. Add option -r to print the offsets of all records.
 
        If no argument is given, command usage information is printed.
 
 OPTIONS
        -v var1[,...]
               The output will only display  the  offset  information  for  the
               specified variables. Names of one or more variables must be pro-
               vided in the comma-delimited list which must not contain  blanks
               or  other  white  space  characters. The named variables must be
               valid netCDF variables in the input  file.  The  default,  i.e.,
               without  this  option,  is to display the offset information for
               all variables stored in the input file.
 
        -s     Print the variable size in bytes. For record variables, only the
               size of one record is printed.
 
        -g     Print the gap in bytes from the previous variable. For the first
               defined variable, print the gap from the end of file header. For
               record variables, there is no gap between records.
 
        -r     Output  the  offset  information for all records of the selected
               record variables.  Without this  option,  only  the  offsets  of
               first record are printed.
 
        -x     Check all fixed-size variable for file space gaps in between any
               two immediately adjacent variables. It prints "1" on  stdout  if
               gaps are found, "0" for otherwise. This option disables all oth-
               er options.
 
        -h     Print the available command-line options
 
 EXAMPLES
        Print the file offset information for all variables in a netCDF file.
 
        % ncoffsets -sg testfile.nc
        netcdf testfile.nc {
        //file format: CDF-1
 
        file header:
             size   = 340 bytes
             extent = 340 bytes
 
        dimensions:
             x = 100
             y = 100
             z = 100
             time = UNLIMITED // (100 currently)
 
        fixed-size variables:
             double square(x, y):
                    start file offset =         340
                    end   file offset =       80340
                    size in bytes     =       80000
                    gap from prev var =           0
             double cube(x, y, z):
                    start file offset =       80340
                    end   file offset =     8080340
                    size in bytes     =     8000000
                    gap from prev var =           0
 
        record variables:
             double time(time):
                    start file offset =     8080340    (record 0)
                    end   file offset =     8081140    (record 0)
                    size in bytes     =           8    (of one record)
                    gap from prev var =           0
             double xytime(time, x, y):
                    start file offset =     8080348    (record 0)
                    end   file offset =    16080348    (record 0)
                    size in bytes     =       80000    (of one record)
                    gap from prev var =           0
        }
 
        Check if there are gaps in between two adjacent fixed-size variables.
 
        % ncoffsets -x testfile.nc
        0
 SEE ALSO
        pnetcdf(3)
 
 DATE
        February 21, 2022
```
</details></li>
</ul>

### ncvalidator
<ul>
  <li>A utility program to check the header of a netCDF file for whether it conforms the
  classic CDF file formats.
  </li>
  <li> <details>
  <summary>Manual page (click to expand)</summary>

```
 ncvalidator(1)                 PnetCDF utilities                ncvalidator(1)
 
 NAME
        ncvalidator - validates a classic netCDF file against CDF file formats
 
 SYNOPSIS
        ncvalidator [-x] [-t] [-q] [-h] file
 
 DESCRIPTION
        ncvalidator  checks the header of a netCDF file for whether it conforms
        the classic CDF file formats. If the input file is a valid NetCDF file,
        then a message of successful validation is printed on command-line out-
        put, for example, File "testfile.nc" is a valid  NetCDF  file.   Other-
        wise, a NetCDF error message is printed.
 
 OPTIONS
        -x     Repair  the null-byte padding in file header. The null-byte pad-
               ding is required by the NetCDF  Classic  Format  Specifications.
               PnetCDF enforces this requirement, but NetCDF has never enforced
               it. This option checks the header for locations where null bytes
               are expected and replaces them with null bytes if non-null bytes
               are found. The repaired file is then conformed with the specifi-
               cation  and allows both PnetCDF and NetCDF libraries to read the
               file without reporting error code NC_ENOTNC or NC_ENULLPAD. Not-
               ed  that  this  repair  is done in place and users might want to
               backup the input file first. Once the file is repaired, one  may
               run ncmpidiff command to compare the contents of two files.
 
        -t     Turn  on  tracing  mode, printing the progress of all successful
               metadata validation. When an  error  is  detected,  the  tracing
               stops at the location of the error found.
 
        -q     Quiet  mode  - print nothing on the command-line output. When in
               quiet mode, users should check exit status. See below  in  "EXIT
               STATUS".
 
        -h     Print the available command-line options
 
 EXIT STATUS
        An exit status of 0 means the file is conform with the classic CDF file
        format, and 1 means otherwise.  Note on VMS-based system, the exit sta-
        tus values are reversed.
 
 SEE ALSO
        ncmpidump(1), pnetcdf(3)
 
 DATE
        February 21, 2022
```
</details></li>
</ul>

### pnetcdf-config
<ul>
  <li>A utility program to display the build and installation information of the
  PnetCDF library.
  </li>
  <li> <details>
  <summary>Manual page (click to expand)</summary>

```
pnetcdf-config is a utility program to display the build and installation
information of the PnetCDF library.

Usage: pnetcdf-config [OPTION]

Available values for OPTION include:

  --help                      display this help message and exit
  --all                       display all options
  --cc                        C compiler used to build PnetCDF
  --cflags                    C compiler flags used to build PnetCDF
  --cppflags                  C pre-processor flags used to build PnetCDF
  --has-c++                   whether C++ API is installed
  --c++                       C++ compiler used to build PnetCDF
  --cxxflags                  C++ compiler flags used to build PnetCDF
  --has-fortran               whether Fortran API is installed
  --f77                       Fortran 77 compiler used to build PnetCDF
  --fflags                    Fortran 77 compiler flags used to build PnetCDF
  --fppflags                  Fortran pre-processor flags used to build PnetCDF
  --fc                        Fortran 9x compiler used to build PnetCDF
  --fcflags                   Fortran 9x compiler flags used to build PnetCDF
  --ldflags                   Linker flags used to build PnetCDF
  --libs                      Extra libraries used to build PnetCDF
  --netcdf4                   Whether NetCDF-4 support is enabled or disabled
  --adios                     Whether ADIOS support is enabled or disabled
  --relax-coord-bound         Whether using a relaxed coordinate boundary check
  --in-place-swap             Whether using buffer in-place Endianness byte swap
  --erange-fill               Whether using fill values for NC_ERANGE error
  --subfiling                 Whether subfiling is enabled or disabled
  --null-byte-header-padding  Whether to check null-byte padding in header
  --burst-buffering           Whether burst buffer driver is built or not
  --profiling                 Whether internal profiling is enabled or not
  --thread-safe               Whether thread-safe capability is enabled or not
  --debug                     Whether PnetCDF is built with debug mode
  --prefix                    Installation directory
  --includedir                Installation directory containing header files
  --libdir                    Installation directory containing library files
  --version                   Library version
  --release-date              Date of PnetCDF source was released
  --config-date               Date of PnetCDF library was configured
```
</details></li>
</ul>

### pnetcdf_version
<ul>
  <li>A utility program to print the version information of PnetCDF library and the
  configure command line used to build the library
  </li>
  <li> <details>
  <summary>Manual page (click to expand)</summary>

```
 pnetcdf_version(1)             PnetCDF utilities            pnetcdf_version(1)
 
 NAME
        pnetcdf_version - print the version information of PnetCDF library
 
 SYNOPSIS
        pnetcdf_version [-v] [-d] [-c] [-b] [-h]
 
 DESCRIPTION
        pnetcdf_version  prints  the version information of PnetCDF library and
        the configure command line used to build the library
 
        If no argument is given, all information is printed.
 
 OPTIONS
        -v     Version number of this PnetCDF release.
 
        -d     Release date.
 
        -c     Configure command-line arguments used to build this PnetCDF
 
        -b     MPI compilers used to build this PnetCDF library
 
        -h     Print the available command-line options of pnetcdf_version
 
 EXAMPLES
        Print all information about the PnetCDF library by running the  command
        with no options.
 
        % pnetcdf_version
 
        PnetCDF Version:         1.12.3
        PnetCDF Release date:    February 21, 2022
        PnetCDF configure:  --with-mpi=/usr/local/bin
        MPICC:  /usr/local/bin/mpicc -g -O2
        MPICXX: /usr/local/bin/mpicxx -g -O2
        MPIF77: /usr/local/bin/mpif77 -g -O2
        MPIF90: /usr/local/bin/mpif90 -g -O2
 
 SEE ALSO
        pnetcdf(3)
 
 DATE
        February 21, 2022
```
</details></li>
</ul>

Copyright (C) 2012, Northwestern University and Argonne National Laboratory.
See COPYRIGHT notice in top-level directory.
