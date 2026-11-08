      subroutine logo_ps(u)

      use kinparams

      implicit none

      integer :: u

      character(80):: prcsnm, prog_vrsn, nicknm, ind

      include 'process_nm.dat'

      write(u,'(a)') '#==============================================='
     &                              //'==============================#'
      write(u,'(a)') '# The multichannel phase space parameterization '
     &             //'for the process:'
      write(u,'(a,a)') '#            ',trim(prcsnm)
      write(u,'(a,a)') '# generated with PSGen v. ',trim(prog_vrsn)
     &     //' by'
      write(u,'(a)')'#                   Karol Kolodziej'
      write(u,'(a)')'#        University of Silesia, '//
     &              'Katowice, Poland'
      write(u,'(a)')'#        e-mail: karol.kolodziej@us.edu.pl, '
      write(u,'(a)')'#        www: https://www.kk.us.edu.pl/'
      write(ind,'(i10)') nkin
      write(u,'(a)') '# with '//trim(adjustl(ind))//' different '//
     & 'parameterizations, among which are:'
      write(ind,'(i10)') nkins
      write(u,'(a)') '# '//trim(adjustl(ind))//' s-channel '//
     & 'parameterizations, as defined in genps.dat'
      write(ind,'(i10)') nkins
      if (nkinst > nkins) then
       write(ind,'(i10)') nkinst-nkins
       write(u,'(a)') '# '//trim(adjustl(ind))//' t-channel '//
     & 'parameterizations'
      end if
      if (nkin > nkinst) then
       write(ind,'(i10)') nkin-nkinst
       write(u,'(a)') '# '//trim(adjustl(ind))//' photon/gluon '//
     & 'emission parameterizations'
      end if
      if (iqprec == 1) then
      write(u,'(a)') '#==============================================='
     &                              //'==============================#'
       write(u,'(a)')'# FOUR MOMENTA GENERATED WITH QUADRUPLE '//
     &               'PRECISION ARITHMETICS'
      end if
      write(u,'(a)') '#==============================================='
     &                              //'==============================#'

      end
