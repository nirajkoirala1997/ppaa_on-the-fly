c.----------------------------------------------------------------------
c.                MRST 
c.
c.         x = parton momentum fraction ( 0 < x < 1 )
c.         scale = factorization scale in GeV
c.         iset = 0 for cteq3 lo, 1 for cteq3 nlo msbar
c. output: f filled with parton PROBABILITY distibution functions uses 
c.         Particle Data Group flavor code convention as follows:
c.           tbar,bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b,t
c.           -6  ,-5  ,-4  ,-3  ,-2  ,-1  ,0,1,2,3,4,5,6
c.
c.----------------------------------------------------------------------
      subroutine pdf(x,scale,f)
      implicit double precision (a-h,o-z)
      dimension df(-6:6),f(-6:6)
      
      y=x

       call evolvePDF(x,scale,f)
        
         do i = -6,6
       f(i) = f(i)/y
         enddo
        return
      end
