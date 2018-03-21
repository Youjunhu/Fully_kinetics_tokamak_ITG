      SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
      use precision,only:p_
      implicit none
      INTEGER m,n,NN
      REAL(p_) x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=200)
!CU    USES spline
      INTEGER j,k
      REAL(p_) y2tmp(NN),ytmp(NN)
      do 13 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
11      continue
        call spline(x2a,ytmp,n,1.d30,1.d30,y2tmp)
        do 12 k=1,n
          y2a(j,k)=y2tmp(k)
12      continue
13    continue
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.


      SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
      use precision,only:p_
      implicit none
      INTEGER m,n,NN
      REAL(p_) x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=200)
!CU    USES spline,splint
      INTEGER j,k
      REAL(p_) y2tmp(NN),ytmp(NN),yytmp(NN)
      do 12 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
          y2tmp(k)=y2a(j,k)
11      continue
        call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
12    continue
      call spline(x1a,yytmp,m,1.d30,1.d30,y2tmp)
      call splint(x1a,yytmp,y2tmp,m,x1,y)
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.
