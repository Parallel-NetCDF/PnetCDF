      character *80 function nfmpi_strerror( err )
      integer err, ierr
      character *(80) tmpstr
C      
C     Call a (C) routine that places the message into tmpstr
      ierr = nfmpi_xstrerror( err, tmpstr )
      nfmpi_strerror = tmpstr
      end
