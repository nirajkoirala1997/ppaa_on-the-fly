      module distribs

      integer :: idis
      integer, parameter :: ndist=2,nbmax=100
      integer, dimension(ndist) :: nbs,ibin,nofr
      integer, dimension(nbmax,ndist) :: ncs
      integer, dimension(nbmax) :: ncsa
      real(kind(1.d0)), dimension(ndist) :: xmin,xmax,xd
      real(kind(1.d0)), dimension(nbmax,ndist) :: dx,di,vdi,sdi
      real(kind(1.d0)), dimension(0:nbmax,ndist) :: xdi
      real(kind(1.d0)), dimension(nbmax) :: cdmax
      real(kind(1.d0)) :: vz,gz
      character(80) desc(ndist)

      end module distribs

