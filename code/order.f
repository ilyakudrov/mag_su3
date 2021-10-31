      subroutine order1(am) 
      include 'paravp3'
c      include 'su2mon.cmn'

      common /var3/ z(nsite,3,3,nd)   ! MC config.
      dimension am(nsite,3,3,nd)

c array from link variables:

      do mu=1,nd
      do m=1,nsite  
        u1=z(m,1,1,mu)*dconjg(z(m,2,2,mu))
        u2=z(m,1,2,mu)*dconjg(z(m,2,1,mu))
        u3=z(m,1,1,mu)*dconjg(z(m,1,2,mu))
        u4=z(m,2,2,mu)*dconjg(z(m,2,1,mu))
        u5=z(m,1,1,mu)*dconjg(z(m,2,1,mu))
        u6=z(m,2,2,mu)*dconjg(z(m,1,2,mu))

   
        am(m,1,1,mu)=dreal(u1+u2)
        am(m,2,2,mu)=dreal(u1-u2)
        am(m,3,3,mu)=0.5*(z(m,1,1,mu)*dconjg(z(m,1,1,mu))+
     *                    z(m,2,2,mu)*dconjg(z(m,2,2,mu))-
     *                    z(m,1,2,mu)*dconjg(z(m,1,2,mu))-
     *                    z(m,2,1,mu)*dconjg(z(m,2,1,mu)))
        am(m,2,1,mu)=-dimag(u1+u2)
        am(m,1,2,mu)=dimag(u1-u2)
        am(m,3,1,mu)=dreal(u3-u4)
        am(m,1,3,mu)=dreal(u5-u6)
        am(m,3,2,mu)=dimag(u3+u4)
        am(m,2,3,mu)=-dimag(u5+u6)

      enddo 
      enddo 
        return
        end
c**********************************************************
c**********************************************************
      subroutine order2(am)

      include 'paravp3'
c      include 'su2mon.cmn'

      common /var3/ z(nsite,3,3,nd)   ! MC config.
      dimension am(nsite,3,3,nd)

c array from link variables:

      do mu=1,nd
      do m=1,nsite  
        u1=z(m,1,1,mu)*dconjg(z(m,3,3,mu))
        u2=z(m,1,3,mu)*dconjg(z(m,3,1,mu))
        u3=z(m,1,1,mu)*dconjg(z(m,1,3,mu))
        u4=z(m,3,3,mu)*dconjg(z(m,3,1,mu))
        u5=z(m,1,1,mu)*dconjg(z(m,3,1,mu))
        u6=z(m,3,3,mu)*dconjg(z(m,1,3,mu))

   
        am(m,1,1,mu)=dreal(u1+u2)
        am(m,2,2,mu)=dreal(u1-u2)
        am(m,3,3,mu)=0.5*(z(m,1,1,mu)*dconjg(z(m,1,1,mu))+
     *                    z(m,3,3,mu)*dconjg(z(m,3,3,mu))-
     *                    z(m,1,3,mu)*dconjg(z(m,1,3,mu))-
     *                    z(m,3,1,mu)*dconjg(z(m,3,1,mu)))
        am(m,2,1,mu)=-dimag(u1+u2)
        am(m,1,2,mu)=dimag(u1-u2)
        am(m,3,1,mu)=dreal(u3-u4)
        am(m,1,3,mu)=dreal(u5-u6)
        am(m,3,2,mu)=dimag(u3+u4)
        am(m,2,3,mu)=-dimag(u5+u6)

      enddo 
      enddo 
        return
        end
c**********************************************************
c**********************************************************
      subroutine order3(am)

      include 'paravp3'
c      include 'su2mon.cmn'

      common /var3/ z(nsite,3,3,nd)   ! MC config.
      dimension am(nsite,3,3,nd)

c array from link variables:

      do mu=1,nd
      do m=1,nsite  
        u1=z(m,2,2,mu)*dconjg(z(m,3,3,mu))
        u2=z(m,2,3,mu)*dconjg(z(m,3,2,mu))
        u3=z(m,2,2,mu)*dconjg(z(m,2,3,mu))
        u4=z(m,3,3,mu)*dconjg(z(m,3,2,mu))
        u5=z(m,2,2,mu)*dconjg(z(m,3,2,mu))
        u6=z(m,3,3,mu)*dconjg(z(m,2,3,mu))

   
        am(m,1,1,mu)=dreal(u1+u2)
        am(m,2,2,mu)=dreal(u1-u2)
        am(m,3,3,mu)=0.5*(z(m,2,2,mu)*dconjg(z(m,2,2,mu))+
     *                    z(m,3,3,mu)*dconjg(z(m,3,3,mu))-
     *                    z(m,2,3,mu)*dconjg(z(m,2,3,mu))-
     *                    z(m,3,2,mu)*dconjg(z(m,3,2,mu)))
        am(m,2,1,mu)=-dimag(u1+u2)
        am(m,1,2,mu)=dimag(u1-u2)
        am(m,3,1,mu)=dreal(u3-u4)
        am(m,1,3,mu)=dreal(u5-u6)
        am(m,3,2,mu)=dimag(u3+u4)
        am(m,2,3,mu)=-dimag(u5+u6)

      enddo 
      enddo 
        return
        end


