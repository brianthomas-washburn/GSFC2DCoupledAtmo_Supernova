      subroutine stridag(a,b,c,r,u,n)
c  from Numerical Recipes
      parameter(nmax=100)
      dimension gam(nmax),a(n),b(n),c(n),r(n),u(n)
      if(b(1).eq.0) pause 'b(1).eq.0'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
          gam(j)=c(j-1)/bet
          bet=b(j)-a(j)*gam(j)
          if(bet.eq.0) pause 'bet.eq.0'
          u(j)=(r(j)-a(j)*u(j-1))/bet
   11 continue
      do 12 j=n-1,1,-1
      u(j)=u(j)-gam(j+1)*u(j+1)
   12 continue
      return
      end         
