      character *80 function nfmpi_inq_libvers()
      character *(80) tmpstr

      call nfmpi_xinq_libvers(tmpstr)
      nfmpi_inq_libvers = tmpstr
      end
