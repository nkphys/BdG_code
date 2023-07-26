 program main
       implicit none
       real(kind=8)::a,b,c,up,dn,add1,add2,add3,u,pi,omga,t
       integer::i,j,k,n,d,ix,iy,jx,jy,m,site,l
!       integer,dimension(13,13)::hop,del
       complex::cor1,img,addc,addc1,addev,add0,addod
       real,allocatable,dimension(:,:)::hop
       real,allocatable,dimension(:)::time
!       open(unit=2,file="HOP.dat")
       open(unit=2,file="thop.txt")
         pi=4.0*atan(1.0d0)
         write(*,*)"Enter the number of sites"
         read(*,*)n
         write(*,*)"Enter the VALUE of t_hop"
         read(*,*)t
         allocate(hop(n,n))
         l=(n-1)/3
          ix=0
          hop=0
         do i=1,n-1
         do j=i+1,n
         if(i.le.l)then
                 if((j-i).eq.1)hop(i,j)=t
         endif
         if(i.eq.l+1)then
                 if(j-i.eq.1)hop(i,j)=t
                 if(j-i.eq.l+1)hop(i,j)=t
         endif
         if(i.gt.(l+1).and.i.lt.(2*l+1))then
                 if(j-i.eq.1)hop(i,j)=t
         endif
         if(i.ge.(2*l+2).and.i.le.(3*l))then
                 if(j-i.eq.1)hop(i,j)=t
                 if(j-i.eq.1)write(*,*)hop(i,j),i,j
         endif
         enddo
         enddo

         do i=1,n
         do j=1,n
         hop(j,i)=hop(i,j)
         enddo
         enddo
         do i=1,n        
!         write(2,*)"[",(hop(i,j),j=1,n),"]"
         write(2,*)(int(hop(i,j)),j=1,n)
         enddo
         do i=1,n        
         do j=1,n        
         if(hop(i,j).eq.1.and.j.gt.i)then
         write(*,*)"i=",i,"j=",j
         endif
         enddo
         enddo
        end program main

