      character *80 function nfmpi_strerror( err )
      integer err
      character *(80) tmpstr

      call nfmpi_xstrerror( err, tmpstr )
      nfmpi_strerror = tmpstr
      end
