      character *80 function nfmpi_strerrno( err )
      integer err, ierr
      character *(80) tmpstr

      integer nfmpi_xstrerrno
      external nfmpi_xstrerrno

      ierr = nfmpi_xstrerrno( err, tmpstr )
      if (tmpstr(1:2) .EQ. 'NC') then
          tmpstr(2:2) = 'F'
      end if
      nfmpi_strerrno = tmpstr
      end
