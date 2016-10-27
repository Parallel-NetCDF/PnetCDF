      character *80 function nfmpi_strerrno( err )
      integer err, ierr
      character *(80) tmpstr

      integer nfmpi_xstrerrno
      external nfmpi_xstrerrno

      ierr = nfmpi_xstrerrno( err, tmpstr )
      nfmpi_strerrno = tmpstr
      end
