      subroutine wrt72(nu,str)
c----------------------------------------------------------------------c
c Writes a character string 'str' in columns 7 to 72 of the file 
c associated with an opened unit 'nu'
c
c Karol Kolodziej, University of Silesia, Karol.Kolodziej@us.edu.pl
c 3 June, 2007; modified 
c----------------------------------------------------------------------c

      implicit none

      character(*) :: str
      character(6), parameter :: ccon='     &'
      character(80) :: clx
      integer :: nu,lstr,iclb,icle

      lstr=len_trim(str)
      iclb=1
      icle=72
      clx=''
      do
       if (lstr > 72) then
        write(nu,'(a)') trim(clx)//str(iclb:icle)
        iclb=icle+1
        icle=icle+66
        clx=ccon
        if (icle > lstr) then
         write(nu,'(a)') trim(clx)//str(iclb:lstr)
         exit
        end if
       else
        write(nu,'(a)') str(iclb:lstr)
        exit
       end if
      end do

      end subroutine wrt72

