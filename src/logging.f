      Subroutine impy_openlogfile(fname, opunit)
        Character*300 fname
        Integer opunit

        If (opunit==0) Then
          opunit = 66
        End If

        Open (opunit, File=fname)

      End Subroutine impy_openlogfile

      Subroutine impy_closelogfile(opunit)
        Integer opunit

        If (opunit/=0 .And. opunit/=6) Then
          Close (Unit=opunit, Iostat=ios)
          If (ios/=0) Stop 'Error closing file unit'
        Else
          write(6,*) 'Error while closing file'
        End If

      End Subroutine impy_closelogfile