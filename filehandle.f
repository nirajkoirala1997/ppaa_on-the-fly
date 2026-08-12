c ~~~~~~~~~~~  *********************     ~~~~~~~~~~~~c        
        subroutine output(run_tag,filename1)
        implicit none
        integer ierr1,ierr2
        character*100 dir_path,run_tag,dir_pathtmp,decision
     . ,filename1,filename

       ! Check if the directory exists
       dir_path ="../summary/"// trim(run_tag)
       filename ="../summary/"// trim(run_tag) //"/"// trim(filename1)
       dir_pathtmp ="../summary/temp_"// trim(run_tag)
c checking Directory
        call system("test -d "// dir_path // " && echo 1 > 
     .   command.txt || echo 0  >
     .        command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr1
       close(13)
       call system("rm command.txt")
       if( ierr1.ne. 1) then
         print*,"Directory not found " //dir_path
         print*,"Making new directory.."
         call sleep(1)
         CALL SYSTEM("mkdir -p " // dir_path )
         CALL SYSTEM("cd " // dir_path// " && echo $PWD" )
         print*,"Directory constructed successfully🤩"
         call sleep(1)
         print*," "
       endif
c checking file
       call system("test -f "// filename // " && echo 1 > 
     .   command.txt || echo 0  >
     .        command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr2
       close(13)
       call system("rm command.txt")

c Proceed with decision       

       if (ierr2 .eq. 1) then
        print*,"File  already exist  "//filename
        print*,"All previous data for this run will be lost. Overwrite ?
     .     [y/n]"
        read*,decision
        if (decision .eq. 'y' .or. decision .eq. 'Y') then
        print*,"Overwriting data in"// filename
        else
           stop
        endif
        elseif (ierr2 .eq. 0 ) then
               print*," "
               print*,"Writing data in the file  "//filename
               print*," "
        endif
        call sleep(1)
        end
c ~~~~~~~~~~~  *********************     ~~~~~~~~~~~~c        

