      subroutine print_masses(u)

      use kinparams

      implicit none

      integer :: u, i

      write(u,'(a)') '#-----------------------------------------------'
     &              //'------------------------------#'
      write(u,'(a)') '# Partcle masses:'
      do i=1,nexp
       if (iqprec == 1) then
        write(u,'(a,e39.32)') '#  m('//trim(adjustl(expc(i)))//') =',
     &                                                            m(i)
       else
        write(u,'(a,e23.16)') '#  m('//trim(adjustl(expc(i)))//') =',
     &                                                            m(i)
       end if
      end do
      write(u,'(a,e23.16)') '#  m_W =',mw
      write(u,'(a,e23.16)') '#  m_Z =',mz
      write(u,'(a,e23.16)') '#  m_H =',mh
      write(u,'(a)') '#-----------------------------------------------'
     &              //'------------------------------#'

      end
