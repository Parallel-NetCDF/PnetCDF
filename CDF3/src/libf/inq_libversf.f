      character *80 function nfmpi_inq_libvers()
      integer ierr
      character *(80) tmpstr
C
C     Call a (C) routine that places the message into tmpstr
      ierr = nfmpi_xinq_libvers(tmpstr)
      nfmpi_inq_libvers = tmpstr
      end
