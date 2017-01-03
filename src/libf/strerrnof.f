      character *80 function nfmpi_strerrno( err )
      integer err
      character *(80) tmpstr

      call nfmpi_xstrerrno( err, tmpstr )
      if (tmpstr(1:2) .EQ. 'NC') then
          tmpstr(2:2) = 'F'
      end if
      nfmpi_strerrno = tmpstr
      end
