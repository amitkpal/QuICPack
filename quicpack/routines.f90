!DEFINITIONS


!
subroutine cbase(z,o) !computational basis for qubits
        implicit none
        complex*16::z(1:2,1:1),o(1:2,1:1)
        z(1,1)=dcmplx(1.0d0,0.0d0)
        z(2,1)=dcmplx(0.0d0,0.0d0)
        o(1,1)=dcmplx(0.0d0,0.0d0)
        o(2,1)=dcmplx(1.0d0,0.0d0)
end subroutine cbase
!
subroutine lproj(th,pmat) !local projectors for qubits
        implicit none
        integer::i
        double precision::th(1:2)
        complex*16::z(1:2,1:1),o(1:2,1:1),b(1:2,1:1,1:2),dcexp,pmat(1:2,1:2,1:2)
        call cbase(z,o)
        b(:,:,1)=dcos(th(1)/2.0d0)*z+dcexp(th(2))*dsin(th(1)/2.0d0)*o
        b(:,:,2)=dcexp(-th(2))*dsin(th(1)/2.0d0)*z-dcos(th(1)/2.0d0)*o
        do i=1,2
           pmat(:,:,i)=matmul(b(:,:,i),conjg(transpose(b(:,:,i))))
        end do
end subroutine lproj
!
subroutine pproj(th,b) !coefficients of basis in qubit hilbert space
        implicit none
        double precision::th(1:2)
        complex*16::b(1:2,1:2),dcexp
        b(1,1)=dcos(th(1)/2.0d0)
        b(1,2)=dcexp(th(2))*dsin(th(1)/2.0d0)
        b(2,1)=dcexp(-th(2))*dsin(th(1)/2.0d0)
        b(2,2)=-dcos(th(1)/2.0d0)
end subroutine pproj
!
subroutine paulipv(n,th) !pauli measurement parameters
        implicit none
        integer::n
        double precision::th(1:2)
        if (n==1) then
           th(1)=dacos(-1.0d0)/2.0d0
           th(2)=0.0d0
        end if
        if (n==2) th=dacos(-1.0d0)/2.0d0
        if (n==3) th=0.0d0
end subroutine paulipv
!
subroutine idmat(n,sd,id) !identity matrix in n-qudit hilbert space
        implicit none
        integer::n,i,sd
        complex*16::id(1:sd**n,1:sd**n)
        id=dcmplx(0.0d0,0.0d0)
        do i=1,sd**n
           id(i,i)=dcmplx(1.0d0,0.0d0)
        end do
end subroutine idmat
!
!RANDOM NUMBER GENERATORS


!
function ran2(idum) !uniform double precision random number between 0 and 1
        implicit none
        integer,parameter::im1=2147483563,im2=2147483399,imm1=im1-1
        integer,parameter::ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211
        integer,parameter::ir2=3791,ntab=32,ndiv=1+imm1/ntab
        double precision,parameter::am=1.0d0/im1,eps=1.2d-7,rnmx=1.0d0-eps
        integer::idum,idum2,j,k,iv(ntab),iy
        double precision::ran2
        save::iv,iy,idum2
        data idum2/123456789/,iv/ntab*0/,iy/0/
        if (idum <= 0) then
           idum=max(-idum,1)
           idum2=idum
           do j=ntab+8,1,-1
              k=idum/iq1
              idum=ia1*(idum-k*iq1)-k*ir1
              if (idum < 0) idum=idum+im1
              if (j <= ntab) iv(j)=idum
           end do
           iy=iv(1)
        end if
        k=idum/iq1
        idum=ia1*(idum-k*iq1)-k*ir1
        if (idum < 0) idum=idum+im1
        k=idum2/iq2
        idum2=ia2*(idum2-k*iq2)-k*ir2
        if (idum2 < 0) idum2=idum2+im2
        j=1+iy/ndiv
        iy=iv(j)-idum2
        iv(j)=idum
        if (iy < 1) iy=iy+imm1
        ran2=dmin1(am*iy,rnmx)
        return
        end





function gasdev(idum) !double precision gaussian random number of 0 mean and 1 standard deviation
        implicit none
        integer::idum,iset
        double precision::gasdev,fac,gset,rsq,v1,v2,ran2
        save::iset,gset
        data iset/0/
        if (idum < 0) iset=0
        if (iset == 0) then
           do
              v1=2.0d0*ran2(idum)-1.0d0
              v2=2.0d0*ran2(idum)-1.0d0
              rsq=v1**2+v2**2
              if (rsq >= 1.0d0 .or. rsq == 0.0d0) then
                 cycle
              else
                 exit
              end if
           end do
           fac=dsqrt(-2.0d0*dlog(rsq)/rsq)
           gset=v1*fac
           gasdev=v2*fac
           iset=1
        else
           gasdev=gset
           iset=0
        endif
        return
        end
!
!OPERATIONS

!pauli matrices
!swansea, 04-2017
subroutine pauli(n,s) 
        implicit none
        integer::n
        complex*16::s(1:2,1:2)
        if (n==0) then
           s(1,1)=dcmplx(1.0d0,0.0d0)
           s(1,2)=dcmplx(0.0d0,0.0d0)
           s(2,1)=dcmplx(0.0d0,0.0d0)
           s(2,2)=dcmplx(1.0d0,0.0d0)
        end if
        if (n==1) then
           s(1,1)=dcmplx(0.0d0,0.0d0)
           s(1,2)=dcmplx(1.0d0,0.0d0)
           s(2,1)=dcmplx(1.0d0,0.0d0)
           s(2,2)=dcmplx(0.0d0,0.0d0)
        end if
        if (n==2) then
           s(1,1)=dcmplx(0.0d0,0.0d0)
           s(1,2)=dcmplx(0.0d0,-1.0d0)
           s(2,1)=dcmplx(0.0d0,1.0d0)
           s(2,2)=dcmplx(0.0d0,0.0d0)
        end if
        if (n==3) then
           s(1,1)=dcmplx(1.0d0,0.0d0)
           s(1,2)=dcmplx(0.0d0,0.0d0)
           s(2,1)=dcmplx(0.0d0,0.0d0)
           s(2,2)=dcmplx(-1.0d0,0.0d0)
        end if
end subroutine pauli
!





!local noisy channels (uncorrelated) for n-qubit system
!swansea,04-2017
subroutine lchanel(lc,rho,n,p,erho)
        implicit none
        integer::n,i,j,k,l,lc,b(1:n)
        double precision::p
        complex*16::rho(1:2**n,1:2**n),erho(1:2**n,1:2**n),q(1:2,1:2)
        complex*16,dimension(:,:),allocatable::id,s,c1,c2,sp
        complex*16,dimension(:,:,:),allocatable::ss
        !global/werner
        if (lc==0) then
           allocate(id(1:2**n,1:2**n))
           call idmat(n,2,id)
           erho=(p/2**n)*id+(1.0d0-p)*rho
           deallocate(id)
        end if
        !bit-flip
        if (lc==1) then
           allocate(id(1:2,1:2))
           allocate(s(1:2,1:2))
           call pauli(0,id)
           call pauli(1,s)
           erho=dcmplx(0.0d0,0.0d0)
           do i=1,2**n
              call d2b(i-1,n,2,b)
              do j=1,n
                 if (b(j)==0) q=dsqrt(p/2.0d0)*s
                 if (b(j)==1) q=dsqrt(1.0d0-(p/2.0d0))*id
                 if (j==1) then
                    allocate(c1(1:2,1:2))
                    c1=q
                 else
                    allocate(c2(1:2**j,1:2**j))
                    call kron(c1,2**(j-1),2**(j-1),q,2,2,c2)
                    deallocate(c1)
                    allocate(c1(1:2**j,1:2**j))
                    c1=c2
                    deallocate(c2)
                 end if
              end do
              erho=erho+matmul(c1,matmul(rho,conjg(transpose(c1))))
              deallocate(c1)
           end do
           deallocate(id,s)
        end if 
        !bit-phase-flip
        if (lc==2) then
           allocate(id(1:2,1:2))
           allocate(s(1:2,1:2))
           call pauli(0,id)
           call pauli(2,s)
           erho=dcmplx(0.0d0,0.0d0)
           do i=1,2**n
              call d2b(i-1,n,2,b)
              do j=1,n
                 if (b(j)==0) q=dsqrt(p/2.0d0)*s
                 if (b(j)==1) q=dsqrt(1.0d0-(p/2.0d0))*id
                 if (j==1) then
                    allocate(c1(1:2,1:2))
                    c1=q
                 else
                    allocate(c2(1:2**j,1:2**j))
                    call kron(c1,2**(j-1),2**(j-1),q,2,2,c2)
                    deallocate(c1)
                    allocate(c1(1:2**j,1:2**j))
                    c1=c2
                    deallocate(c2)
                 end if
              end do
              erho=erho+matmul(c1,matmul(rho,conjg(transpose(c1))))
              deallocate(c1)
           end do
           deallocate(id,s)
        end if
        !phase-flip
        if (lc==3) then
           allocate(id(1:2,1:2))
           allocate(s(1:2,1:2))
           call pauli(0,id)
           call pauli(3,s)
           erho=dcmplx(0.0d0,0.0d0)
           do i=1,2**n
              call d2b(i-1,n,2,b)
              do j=1,n
                 if (b(j)==0) q=dsqrt(p/2.0d0)*s
                 if (b(j)==1) q=dsqrt(1.0d0-(p/2.0d0))*id
                 if (j==1) then
                    allocate(c1(1:2,1:2))
                    c1=q
                 else
                    allocate(c2(1:2**j,1:2**j))
                    call kron(c1,2**(j-1),2**(j-1),q,2,2,c2)
                    deallocate(c1)
                    allocate(c1(1:2**j,1:2**j))
                    c1=c2
                    deallocate(c2)
                 end if
              end do
              erho=erho+matmul(c1,matmul(rho,conjg(transpose(c1))))
              deallocate(c1)
           end do
           deallocate(id,s)
        end if
        !amplitude-damping
        if (lc==4) then
           allocate(id(1:2,1:2))
           allocate(s(1:2,1:2))
           call pauli(0,id)
           id(2,2)=dcmplx(dsqrt(1.0d0-p),0.0d0)
           s=dcmplx(0.0d0,0.0d0)
           s(1,2)=dcmplx(dsqrt(p),0.0d0)
           erho=dcmplx(0.0d0,0.0d0)
           do i=1,2**n
              call d2b(i-1,n,2,b)
              k=0
              do j=1,n
                 if (b(j)==0) q=id
                 if (b(j)==1) q=s
                 if (j==1) then
                    allocate(c1(1:2,1:2))
                    c1=q
                 else
                    allocate(c2(1:2**j,1:2**j))
                    call kron(c1,2**(j-1),2**(j-1),q,2,2,c2)
                    deallocate(c1)
                    allocate(c1(1:2**j,1:2**j))
                    c1=c2
                    deallocate(c2)
                 end if
              end do
              erho=erho+matmul(c1,matmul(rho,conjg(transpose(c1))))
              deallocate(c1)
           end do
           deallocate(id,s)
        end if
        !depolarizing
        if (lc==5) then
           allocate(id(1:2,1:2))
           allocate(ss(1:4,1:2,1:2))
           do i=1,4
              call pauli(i-1,ss(i,:,:))
           end do
           erho=dcmplx(0.0d0,0.0d0)
           do i=1,4**n
              call d2b(i-1,n,4,b)
              k=0
              do j=1,n
                 do l=1,4
                    if (b(j)==l-1.and.l==1) q=dsqrt(1.0d0-(3.0d0*p/4.0d0))*ss(l,:,:)
                    if (b(j)==l-1.and.l/=1) q=dsqrt(p/4.0d0)*ss(l,:,:)
                 end do
                 if (j==1) then
                    allocate(c1(1:2,1:2))
                    c1=q
                 else
                    allocate(c2(1:2**j,1:2**j))
                    call kron(c1,2**(j-1),2**(j-1),q,2,2,c2)
                    deallocate(c1)
                    allocate(c1(1:2**j,1:2**j))
                    c1=c2
                    deallocate(c2)
                 end if
              end do
              erho=erho+matmul(c1,matmul(rho,conjg(transpose(c1))))
              deallocate(c1)
           end do
           deallocate(id,ss)
        end if
        !phase-damping
        if (lc==6) then
           allocate(id(1:2,1:2))
           allocate(sp(1:2,1:2))
           allocate(ss(1:4,1:2,1:2))
           call pauli(0,id)
           call pauli(3,sp)
           ss(1,:,:)=dsqrt(1.0d0-p)*id
           ss(2,:,:)=(dsqrt(p)/2.0d0)*(id+sp)
           ss(3,:,:)=(dsqrt(p)/2.0d0)*(id-sp)
           erho=dcmplx(0.0d0,0.0d0)
           do i=1,3**n
              call d2b(i-1,n,3,b)
              k=0
              do j=1,n
                 do l=1,3
                    if (b(j)==l-1) q=ss(l,:,:)
                 end do
                 if (j==1) then
                    allocate(c1(1:2,1:2))
                    c1=q
                 else
                    allocate(c2(1:2**j,1:2**j))
                    call kron(c1,2**(j-1),2**(j-1),q,2,2,c2)
                    deallocate(c1)
                    allocate(c1(1:2**j,1:2**j))
                    c1=c2
                    deallocate(c2)
                 end if
              end do
              erho=erho+matmul(c1,matmul(rho,conjg(transpose(c1))))
              deallocate(c1)
           end do
           deallocate(ss,sp,id)
        end if
end subroutine lchanel
!



!unitary matrix in qubit hilbert space
!swansea, 04-2017
subroutine unimat(th,uni) 
        implicit none
        double precision::th(1:2)
        complex*16::dcexp,uni(1:2,1:2)
        uni(1,1)=dcos(th(1)/2.0d0)
        uni(1,2)=dsin(th(1)/2.0d0)*dcexp(th(2))
        uni(2,1)=dsin(th(1)/2.0d0)*dcexp(-th(2))
        uni(2,2)=-dcos(th(1)/2.0d0)
end subroutine unimat
!



!kronecker product of two states
!swansea, 04-2017
subroutine kron(a,ar,ac,b,br,bc,m) 
        implicit none
        integer::ar,ac,br,bc,i,j,k,l,temp1,temp2
        complex*16::a(1:ar,1:ac),b(1:br,1:bc),m(1:ar*br,1:ac*bc)
        temp1=0
        do i = 1,ar
           temp2=0
           do j = 1,ac
              do k = 1,br
                 do l = 1,bc
                    m(k+temp1,l+temp2)=a(i,j)*b(k,l)
                 end do
              end do
              temp2=temp2+bc
           end do
           temp1=temp1+br
        end do
end subroutine kron
!










!trace of an n-qudit state
!swansea, 04-2017
subroutine trace(rho,n,sd,tr)
        implicit none
        integer::i,n,sd
        double precision::tr
        complex*16::rho(1:sd**n,1:sd**n)
        tr=0.0d0
        do i=1,sd**n
           tr=tr+real(rho(i,i))
        end do
end subroutine trace
!



!partial trace for n-qudit mixed states
!swansea, 04-2017
subroutine mptrace(rho,n,sd,p,m,sig)
        integer::n,m,p(1:m),flag1,flag2,ind,d1,d2,i,j,k,l,b1(1:n),b2(1:n),br1(1:n-m),br2(1:n-m),sd
        complex*16::rho(1:sd**n,1:sd**n),sig(1:sd**(n-m),1:sd**(n-m))
        sig=dcmplx(0.0d0,0.0d0)
        do i=1,sd**n
           do j=1,sd**n
              if (j>=i) then
                 call d2b(i-1,n,sd,b1)
                 call d2b(j-1,n,sd,b2)
                 flag1=0
                 do k=1,m
                    if (b1(p(k))/=b2(p(k))) then
                       flag1=1
                       exit
                    end if
                 end do
                 if (flag1==0) then
                    ind=0
                    do k=1,n
                       flag2=0
                       do l=1,m
                          if (p(l)==k) then
                             flag2=1
                             exit
                          end if
                       end do
                       if (flag2==0) then
                          ind=ind+1
                          br1(ind)=b1(k)
                          br2(ind)=b2(k)
                       end if
                    end do
                    call b2d(br1,n-m,sd,d1)
                    call b2d(br2,n-m,sd,d2)
                    sig(d1,d2)=sig(d1,d2)+rho(i,j)
                    if (d1/=d2) sig(d2,d1)=sig(d2,d1)+conjg(rho(i,j))
                 end if
              end if
           end do
        end do
end subroutine mptrace
!



!partial trace for n-qudit pure states
!swansea, 04-2017
subroutine pptrace(psi,n,sd,p,nr,rho) 
        implicit none
        integer::i,j,k,l,m,n,nr,p(1:nr),b(1:n),br(1:nr),sd
        complex*16::psi(1:sd**n),rho(1:sd**(n-nr),1:sd**(n-nr))
        complex*16,dimension(:,:),allocatable::c
        rho=dcmplx(0.0d0,0.0d0)
        do i=1,sd**nr
           call d2b(i-1,nr,sd,br)
           allocate(c(1:sd**(n-nr),1:1))
           m=0
           do j=1,sd**n
              call d2b(j-1,n,sd,b)
              l=1
              do k=1,nr
                 if (br(k)/=b(p(k))) then
                    l=0
                 end if
              end do
              if (l==1) then
                 m=m+1
                 c(m,1)=psi(j)
              end if
           end do
           rho=rho+matmul(c,conjg(transpose(c)))
           deallocate(c)
        end do
end subroutine pptrace
!



!partial transpose of an n-qudit state
!swansea, 04-2017
subroutine partran(rho,n,sd,m,trho)  
        implicit none 
        integer::n,m,i,j,k,l,p,q,iflag,jflag,sd
        complex*16::rho(1:sd**n,1:sd**n),trho(1:sd**n,1:sd**n),bl(1:sd**(n-m),1:sd**(n-m),1:sd**m,1:sd**m)
        iflag=0
        k=0
        p=1
        do i=1,sd**n
           k=k+1
           if (k==sd**(n-m)) iflag=1
           if (k>sd**(n-m)) then
              if (iflag==1) then
                 p=p+1
                 k=1
                 iflag=0
              end if
           end if
           jflag=0
           l=0
           q=1
           do j=1,sd**n
              l=l+1
              if (l==sd**(n-m)) jflag=1
              if (l>sd**(n-m)) then
                 if (jflag==1) then
                    q=q+1
                    l=1
                    jflag=0
                 end if
              end if
              bl(k,l,p,q)=rho(i,j)     
           end do
        end do
        do p=1,sd**m
           do q=1,sd**m
              bl(:,:,p,q)=transpose(bl(:,:,p,q))
           end do
        end do
        iflag=0
        k=0
        p=1
        do i=1,sd**n
           k=k+1
           if (k==sd**(n-m)) iflag=1
           if (k>sd**(n-m)) then
              if (iflag==1) then
                 p=p+1
                 k=1
                 iflag=0
              end if
           end if
           jflag=0
           l=0
           q=1
           do j=1,sd**n
              l=l+1
              if (l==sd**(n-m)) jflag=1
              if (l>sd**(n-m)) then
                 if (jflag==1) then
                    q=q+1
                    l=1
                    jflag=0
                 end if
              end if
              trho(i,j)=bl(k,l,p,q)
           end do
        end do
end subroutine partran
!




!swapping two qudits in a mixed state
!swansea, 04-2017
subroutine mswap(rho,n,sd,a,b,erho)  
        implicit none 
        integer::n,a,b,i,j,k,l,base(1:n,1:2),d1,d2,sd
        complex*16::rho(1:sd**n,1:sd**n),erho(1:sd**n,1:sd**n)
        do i=1,sd**n
           do j=1,sd**n
              if (i>=j) then
                 call d2b(i-1,n,sd,base(:,1))
                 call d2b(j-1,n,sd,base(:,2))
                 do k=1,2
                    l=base(a,k)
                    base(a,k)=base(b,k)
                    base(b,k)=l
                 end do
                 call b2d(base(:,1),n,sd,d1)
                 call b2d(base(:,2),n,sd,d2)
                 erho(d1,d2)=rho(i,j)
                 erho(d2,d1)=rho(j,i)
              end if
           end do
        end do
end subroutine mswap
!



!swapping two qudits in a pure state
!swansea, 04-2017
subroutine pswap(psi,n,sd,a,b,phi)  
        implicit none 
        integer::n,a,b,i,j,k,l,base(1:n),d,sd
        complex*16::psi(1:sd**n,1:1),phi(1:sd**n,1:2)
        phi=dcmplx(0.0d0,0.0d0)
        do i=1,sd**n
           call d2b(i-1,n,sd,base)
           l=base(a)
           base(a)=base(b)
           base(b)=l
           call b2d(base,n,sd,d)
           phi(d,1)=psi(i,1)
        end do
end subroutine pswap
!




!single qubit phase gate
!swansea, 10-2016
subroutine phgate(sig)
        implicit none
        complex*16::sig(1:2,1:2)
        sig=dcmplx(0.0d0,0.0d0)
        sig(1,1)=dcmplx(1.0d0,0.0d0)
        sig(2,2)=dcmplx(0.0d0,1.0d0)
end subroutine phgate
!




!single qubit hadamard gate 
!swansea, 10-2016
subroutine hdmd(uh)
        implicit none
        complex*16::uh(1:2,1:2)
        uh(1,1)=dcmplx(1.0d0,0.0d0)
        uh(1,2)=dcmplx(1.0d0,0.0d0)
        uh(2,1)=dcmplx(1.0d0,0.0d0)
        uh(2,2)=dcmplx(-1.0d0,0.0d0)
        uh=(1.0d0/dsqrt(2.0d0))*uh
end subroutine hdmd
!



!exponential of a pauli matrix
!swansea, 10-2016
subroutine plexp(n,x,a) 
        implicit none
        integer::n
        double precision::x
        complex*16::p(1:2,1:2),a(1:2,1:2),id(1:2,1:2)
        call pauli(n,p)
        call pauli(0,id)
        a=dcos(x)*id+dcmplx(0.0d0,1.0d0)*dsin(x)*p
end subroutine plexp
!
!QUANTITIES

!maximum eigenvalue of an n-qudit mixed state
!allahabad, 07-2016
subroutine evlmax(a,n,sd,evm)
        implicit none
        integer::n,inf,sd
        double precision::hw(1:3*(sd**n)-2),ev(1:sd**n),evm
        complex*16::a(1:sd**n,1:sd**n),owork(1:2*(sd**n)+100),b(1:sd**n,1:sd**n)
        character::ch1='N',ch2='U'
        b=a
        call zheev(ch1,ch2,sd**n,b,sd**n,ev,owork,2*(sd**n)+100,hw,inf)
        evm=maxval(ev)
end subroutine evlmax
!


!minimum eigenvalue of an n-qudit mixed state
!allahabad, 07-2016
subroutine evlmin(a,n,sd,evm) 
        implicit none
        integer::n,inf,i,sd
        double precision::hw(1:3*(sd**n)-2),ev(1:sd**n),evm
        complex*16::a(1:sd**n,1:sd**n),owork(1:2*(sd**n)+100),b(1:sd**n,1:sd**n)
        character::ch1='N',ch2='U'
        b=a
        call zheev(ch1,ch2,sd**n,b,sd**n,ev,owork,2*(sd**n)+100,hw,inf)
        evm=minval(ev)
end subroutine evlmin
!


!von neumann entropy
!allahabad, 01-2016
subroutine vnent(a,n,sd,s) 
        implicit none
        integer::n,i,inf,sd
        double precision::s,hw(1:3*(sd**n)-2),ev(1:sd**n),dlog2
        complex*16::a(1:sd**n,1:sd**n),b(1:sd**n,1:sd**n),owork(1:2*(sd**n)+100)
        character::ch1='N',ch2='U'
        b=a
        call zheev(ch1,ch2,sd**n,b,sd**n,ev,owork,2*(sd**n)+100,hw,inf)
        s=0.0d0
        do i=1,sd**n
           if (ev(i)>0.0d0) s=s-ev(i)*dlog2(ev(i))
        end do
end subroutine vnent
!



!generalized geometric measure (ggm) for n-qudit pure state
!allahabad, 01-2016
subroutine ggm_full(psi,n,sd,ggm,kut) 
        implicit none 
        integer::n,nr,m,mold,i,j,k,l,u,v,kut,sd
        integer,dimension(:,:),allocatable::p,q
        integer,dimension(:),allocatable::cp
        double precision::evm,evmax,ggm
        complex*16::psi(1:sd**n)
        complex*16,dimension(:,:),allocatable::rrho
        if (n-n/2>n/2) nr=(n-1)/2
        if (n-n/2==n/2) nr=n/2
        v=0
        do k=1,nr
           call comb_l(n,k,m)
           if (k==1) then
              allocate(p(1:m,1:k))
              do i=1,m
                 do j=1,k
                    p(i,j)=i
                 end do
                 v=v+1
                 allocate(cp(1:n-k))
                 call compl(n,k,p(i,:),cp)
                 allocate(rrho(1:sd**k,1:sd**k))
                 call pptrace(psi,n,sd,cp,n-k,rrho)
                 call evlmax(rrho,k,sd,evm)
                 deallocate(rrho,cp)
                 if (v==1) then
                    evmax=evm
                    kut=k
                 else
                    if (evmax<evm) then
                       evmax=evm
                       kut=k
                    end if
                 end if
              end do
              allocate(q(1:m,1:k))
              q=p
              deallocate(p)
           else
              allocate(p(1:m,1:k))
              l=0
              do i=1,mold
                 if (q(i,k-1)<n) then
                    do u=q(i,k-1)+1,n
                       l=l+1
                       do j=1,k
                          if (j<k) then
                             p(l,j)=q(i,j)
                          else
                             p(l,k)=u
                          end if
                       end do
                       v=v+1
                       allocate(cp(1:n-k))
                       call compl(n,k,p(l,:),cp)
                       allocate(rrho(1:sd**k,1:sd**k))
                       call pptrace(psi,n,sd,cp,n-k,rrho)
                       call evlmax(rrho,k,sd,evm)
                       deallocate(rrho,cp)
                       if (v==1) then
                          evmax=evm
                          kut=k
                       else
                          if (evmax<evm) then
                             evmax=evm
                             kut=k
                          end if
                       end if
                    end do
                 end if
              end do
              deallocate(q)
              allocate(q(1:m,1:k))
              q=p
              deallocate(p)
           end if
           mold=m
        end do
        ggm=1.0d0-evmax
end subroutine ggm_full
!


!ggm (modified, upto given partition nr:rest) for n-qudit pure state
!allahabad, 01-2016
subroutine ggm_mod(psi,n,sd,nr,ggm,kut) 
        implicit none 
        integer::n,nr,m,mold,i,j,k,l,u,v,kut,sd
        integer,dimension(:,:),allocatable::p,q
        integer,dimension(:),allocatable::cp
        double precision::evm,evmax,ggm
        complex*16::psi(1:sd**n)
        complex*16,dimension(:,:),allocatable::rrho
        v=0
        do k=1,nr
           call comb_l(n,k,m)
           if (k==1) then
              allocate(p(1:m,1:k))
              do i=1,m
                 do j=1,k
                    p(i,j)=i
                 end do
                 v=v+1
                 allocate(cp(1:n-k))
                 call compl(n,k,p(i,:),cp)
                 allocate(rrho(1:sd**k,1:sd**k))
                 call pptrace(psi,n,sd,cp,n-k,rrho)
                 call evlmax(rrho,k,sd,evm)
                 deallocate(rrho,cp)
                 if (v==1) then
                    evmax=evm
                    kut=k
                 else
                    if (evmax<evm) then
                       evmax=evm
                       kut=k
                    end if
                 end if
              end do
              allocate(q(1:m,1:k))
              q=p
              deallocate(p)
           else
              allocate(p(1:m,1:k))
              l=0
              do i=1,mold
                 if (q(i,k-1)<n) then
                    do u=q(i,k-1)+1,n
                       l=l+1
                       do j=1,k
                          if (j<k) then
                             p(l,j)=q(i,j)
                          else
                             p(l,k)=u
                          end if
                       end do
                       v=v+1
                       allocate(cp(1:n-k))
                       call compl(n,k,p(l,:),cp)
                       allocate(rrho(1:sd**k,1:sd**k))
                       call pptrace(psi,n,sd,cp,n-k,rrho)
                       call evlmax(rrho,k,sd,evm)
                       deallocate(rrho,cp)
                       if (v==1) then
                          evmax=evm
                          kut=k
                       else
                          if (evmax<evm) then
                             evmax=evm
                             kut=k
                          end if
                       end if
                    end do
                 end if
              end do
              deallocate(q)
              allocate(q(1:m,1:k))
              q=p
              deallocate(p)
           end if
           mold=m
        end do
        ggm=1.0d0-evmax
end subroutine ggm_mod
!



!logarithmic negativity for n-qudit state
!swansea, 09-2016
subroutine logneg(rho,n,sd,m,neg,lneg) 
        implicit none
        integer::n,m,i,inf,j,k,sd
        complex*16::rho(1:sd**n,1:sd**n)
        complex*16::trho(1:sd**n,1:sd**n),owork(1:2*(sd**n)+100)
        complex*16::eval(1:sd**n),evl(1:sd**n,1:sd**n),evr(1:sd**n,1:sd**n)
        double precision::neg,lneg,ev(1:sd**n),hw(1:3*(sd**n)-2),dlog2
        character::ch1='N',ch2='N'
        call partran(rho,n,sd,m,trho)
        call zgeev(ch1,ch2,sd**n,trho,sd**n,eval,evl,sd**n,evr,sd**n,owork,2*(sd**n)+100,hw,inf)
        ev=real(eval)
        neg=0.0d0
        do i=1,sd**n
           if (ev(i)<0.0d0) neg=neg-ev(i)
        end do
        lneg=dlog2(2.0d0*neg+1.0d0)
end subroutine logneg
!



!unmeasured quantum mutual information for n-qudit system
!swansea, 09-2016
subroutine minfo(rho,n,sd,A,m,mi) 
        implicit none 
        integer::n,m,A(1:m),B(1:n-m),i,j,k,l,sd
        double precision::sA,sB,sAB,mi
        complex*16::rho(1:sd**n,1:sd**n),rhoA(1:sd**m,1:sd**m),rhoB(1:sd**(n-m),1:sd**(n-m))
        call compl(n,m,A,B)
        call mptrace(rho,n,sd,B,n-m,rhoA)
        call mptrace(rho,n,sd,A,m,rhoB)
        call vnent(rho,n,sd,sAB)
        call vnent(rhoA,m,sd,sA)
        call vnent(rhoB,n-m,sd,sB)
        mi=sA+sB-sAB
end subroutine minfo
!



!total correlation for n-qudit system
!swansea, 09-2016
subroutine tcorr(rho,n,sd,tc) 
        implicit none
        integer::n,A(1:1),B(1:n-1),i,j,k,l,sd
        double precision::s,tc,ss
        complex*16::rho(1:sd**n,1:sd**n),sig(1:sd,1:sd)
        call vnent(rho,n,sd,ss)
        tc=-ss
        do i=1,n
           A(1)=i
           call compl(n,1,A,B)
           call mptrace(rho,n,sd,B,n-1,sig)
           call vnent(sig,1,sd,s)
           tc=tc+s
        end do
end subroutine tcorr
!




!monogamy score of negativity raised to power 'pow'
!swansea,04-2017
subroutine enms(rho,n,sd,node,pow,ems)
        implicit none 
        integer::i,j,k,l,n,node,pow,sd,ipos(1:2),jpos(1:n-2)
        complex*16::rho(1:sd**n,1:sd**n),erho(1:sd**n,1:sd**n),sig(1:sd**2,1:sd**2)
        double precision::ems,neg,lneg
        if (node/=1) then
           call mswap(rho,n,sd,1,node,erho)
        else
           erho=rho
        end if
        call logneg(erho,n,sd,1,neg,lneg)
        ems=neg**pow
        ipos(1)=1
        do i=2,n
           ipos(2)=i
           call compl(n,2,ipos,jpos)
           call mptrace(erho,n,sd,jpos,n-2,sig)
           call logneg(sig,2,sd,1,neg,lneg)
           ems=ems-neg**pow        
        end do
end subroutine enms
!


!concurrence of two-qubit mixed state
!allahabad,01-2015
subroutine concur(a,c) 
        implicit none
        integer::i, j, inf
        double precision::rwrk(1:2*4), ev(1:4), temp, c
        complex*16::a(1:4,1:4), y(1:2,1:2), yy(1:4,1:4), mw(1:4,1:4), owork(1:100)
        complex*16::eval(1:4), evl(1:4,1:4), evr(1:4,1:4)
        character::ch1='N',ch2='N'
        call pauli(2,y)
        call kron(y,2,2,y,2,2,yy)
        mw=matmul(a,matmul(yy,matmul(conjg(a),yy)))
        call zgeev(ch1,ch2,4,mw,4,eval,evl,4,evr,4,owork,100,rwrk,inf)
        do i=1,4
           ev(i)=real(eval(i))
        end do
        do i=1,3
           do j=i+1,4
              if (ev(i)<ev(j)) then
                 temp=ev(i)
                 ev(i)=ev(j)
                 ev(j)=temp
              end if
           end do
        end do
        c=0.
        do i=1,4
           if (i==1) then
              c=c+dsqrt(dabs(ev(i)))
           else
              c=c-dsqrt(dabs(ev(i)))
           end if
        end do
        if(c<0.d0)c=0.d0
end subroutine concur
!










!conditional entropy with measurement on 1
!swansea, 04-2017
subroutine condent(rho,n,th,sc)
        implicit none
        integer::i,k,n
        double precision::th(1:2),sc,prob,sct
        complex*16::rho(1:2**n,1:2**n),mb(1:2**n,1:2**n),mop(1:2**n,1:2**n)
        complex*16::proj(1:2,1:2,1:2),s0(1:2**(n-1),1:2**(n-1))
        call lproj(th,proj)
        call idmat(n-1,2,s0)
        sc=0.0d0
        do i=1,2
           call kron(proj(:,:,i),2,2,s0,2**(n-1),2**(n-1),mop)
           mb=matmul(mop,matmul(rho,mop))
           prob=0.0d0
           do k=1,2**n
              prob=prob+real(mb(k,k))
           end do
           if (prob>0.0d0) mb=(1.0d0/prob)*mb
           call vnent(mb,n,2,sct)
           sc=sc+prob*sct
        end do 
end subroutine condent
!


!minimum conditional entropy with measurement on 1
!swansea, 04-2017
subroutine scmin(rho,n,idum,nlarge,algo,thm,scm)
        implicit none
        integer::i,j,k,nlarge,idum,algo,ires,iout,n
        integer*8::opt
        double precision::th(1:2),thm(1:2),scm,sc,ran2,lb(1:2),ub(1:2)
        complex*16::rho(1:2**n,1:2**n)
        external::scf2,scf3,scf4,scf5,scf6
        do i=1,nlarge
           do j=1,2
              th(j)=ran2(idum)
              th(j)=th(j)*dfloat(j)*dacos(-1.0d0)
           end do
           call condent(rho,n,th,sc)
           if (i==1) then
              scm=sc
              thm=th
           else
              if (sc<scm) then
                 scm=sc
                 thm=th
              end if
           end if
        end do 
        if (scm>0.0d0) then
           iout=0
           call nlo_create(opt,algo,2)
           if (n==2) call nlo_set_min_objective(ires,opt,scf2,rho)
           if (n==3) call nlo_set_min_objective(ires,opt,scf3,rho)
           if (n==4) call nlo_set_min_objective(ires,opt,scf4,rho)
           if (n==5) call nlo_set_min_objective(ires,opt,scf5,rho)
           if (n==6) call nlo_set_min_objective(ires,opt,scf6,rho)
           lb=0.0d0
           do j=1,2
              ub(j)=dfloat(j)*dacos(-1.0d0)
           end do
           call nlo_set_lower_bounds(ires,opt,lb)
           call nlo_set_upper_bounds(ires,opt,ub)
           call nlo_set_ftol_abs(ires,opt,1.0d-7)
           call nlo_optimize(ires,opt,thm,scm)
           if (ires>=0) iout=1
           call nlo_destroy(opt)
        end if
end subroutine scmin
subroutine scf2(scm,n,thm,grad,need_gradient,rho)
        implicit none
        integer,parameter::m=2
        integer::n,need_gradient
        complex*16::rho(1:2**m,1:2**m)
        double precision::thm(1:n),scm,grad(1:n) 
        call condent(rho,m,thm,scm)
end subroutine scf2
subroutine scf3(scm,n,thm,grad,need_gradient,rho)
        implicit none
        integer,parameter::m=3
        integer::n,need_gradient
        complex*16::rho(1:2**m,1:2**m)
        double precision::thm(1:n),scm,grad(1:n) 
        call condent(rho,m,thm,scm)
end subroutine scf3
subroutine scf4(scm,n,thm,grad,need_gradient,rho)
        implicit none
        integer,parameter::m=4
        integer::n,need_gradient
        complex*16::rho(1:2**m,1:2**m)
        double precision::thm(1:n),scm,grad(1:n) 
        call condent(rho,m,thm,scm)
end subroutine scf4
subroutine scf5(scm,n,thm,grad,need_gradient,rho)
        implicit none
        integer,parameter::m=5
        integer::n,need_gradient
        complex*16::rho(1:2**m,1:2**m)
        double precision::thm(1:n),scm,grad(1:n) 
        call condent(rho,m,thm,scm)
end subroutine scf5
subroutine scf6(scm,n,thm,grad,need_gradient,rho)
        implicit none
        integer,parameter::m=6
        integer::n,need_gradient
        complex*16::rho(1:2**m,1:2**m)
        double precision::thm(1:n),scm,grad(1:n) 
        call condent(rho,m,thm,scm)
end subroutine scf6
!


!classical correlation for n qubit state
!swansea, 04-2017
subroutine ccorr(rho,n,idum,nlarge,algo,thm,cc) 
        implicit none 
        integer::n,idum,algo,nlarge,mpos(1:1)
        double precision::thm(1:2),scm,sB,cc
        complex*16::rho(1:2**n,1:2**n),rhoB(1:2**(n-1),1:2**(n-1))
        mpos(1)=1
        call mptrace(rho,n,2,mpos,1,rhoB)
        call vnent(rhoB,n-1,2,sB)
        call scmin(rho,n,idum,nlarge,algo,thm,scm)
        cc=sB-scm
end subroutine ccorr
!



!quantum discord for n-qubit system
!swansea, 04-2017
subroutine qdis(rho,n,idum,nlarge,algo,thm,qd) 
        implicit none 
        integer::n,idum,algo,nlarge,mpos(1:1)
        double precision::cc,thm(1:2),mi,qd
        complex*16::rho(1:2**n,1:2**n)
        mpos(1)=1
        call minfo(rho,n,2,mpos,1,mi)
        call ccorr(rho,n,idum,nlarge,algo,thm,cc)
        qd=mi-cc
end subroutine qdis
!





!monogamy score of quantum discord raised to power 'pow'
!swansea,04-2017
subroutine qdms(rho,n,idum,nlarge,algo,node,pow,dms)
        implicit none 
        integer::i,j,k,l,n,node,pow,sd,ipos(1:2),jpos(1:n-2),idum,nlarge,algo
        complex*16::rho(1:2**n,1:2**n),erho(1:2**n,1:2**n),sig(1:2**2,1:2**2)
        double precision::dms,thm(1:2),qd
        if (node/=1) then
           call mswap(rho,n,2,1,node,erho)
        else
           erho=rho
        end if
        call qdis(erho,n,idum,nlarge,algo,thm,qd) 
        dms=qd**pow
        ipos(1)=1
        do i=2,n
           ipos(2)=i
           call compl(n,2,ipos,jpos)
           call mptrace(erho,n,2,jpos,n-2,sig)
           call qdis(sig,2,idum,nlarge,algo,thm,qd)
           dms=dms-qd**pow        
        end do
end subroutine qdms
!
!STATES



!bell vectors
!swansea, 04-2017
subroutine bellvect(k,psi)
        implicit none
        integer::i,j,k
        complex*16::psi(1:4,1:1)
        psi=dcmplx(0.0d0,0.0d0)
        if (k==1) then
           psi(1,1)=dcmplx(dsqrt(1.0d0/2.0d0),0.0d0)
           psi(4,1)=dcmplx(dsqrt(1.0d0/2.0d0),0.0d0)
        end if
        if (k==2) then
           psi(1,1)=dcmplx(dsqrt(1.0d0/2.0d0),0.0d0)
           psi(4,1)=dcmplx(-dsqrt(1.0d0/2.0d0),0.0d0)
        end if
        if (k==3) then
           psi(2,1)=dcmplx(dsqrt(1.0d0/2.0d0),0.0d0)
           psi(3,1)=dcmplx(dsqrt(1.0d0/2.0d0),0.0d0)
        end if
        if (k==4) then
           psi(2,1)=dcmplx(dsqrt(1.0d0/2.0d0),0.0d0)
           psi(3,1)=dcmplx(-dsqrt(1.0d0/2.0d0),0.0d0)
        end if
end subroutine bellvect
!





!random n qubit pure state
!swansea, 04-2017
subroutine randpsi(idum,n,sd,s) 
        implicit none
        integer::n,i,idum,sd
        double precision::c1,c2,nc,gasdev
        complex*16::s(1:sd**n,1:1)
        nc=0.0d0
        do i=1,sd**n
           c1=gasdev(idum)
           c2=gasdev(idum)
           s(i,1)=dcmplx(c1,c2)
           nc=nc+s(i,1)*conjg(s(i,1))
        end do
        s=s/dsqrt(nc)
end subroutine randpsi
!




!random 2 qubit state of rank 2
!swansea, 04-2017
subroutine q2r2(idum,rho) 
        implicit none
        integer::r(1:1),idum
        complex*16::psi(1:8,1:1),rho(1:4,1:4)
        r(1)=3
        call randpsi(idum,3,2,psi)
        call pptrace(psi,3,2,r,1,rho)
end subroutine q2r2
!





!random 2 qubit state of rank 3
!swansea, 04-2017
subroutine q2r3(idum,rho) 
        implicit none
        integer::i,j,k,idum
        double precision::c1,c2,nc,gasdev
        complex*16::a(1:12),rho(1:4,1:4),dum(1:4,1:4)
        nc=0.0d0
        do i=1,12
           c1=gasdev(idum)
           c2=gasdev(idum)
           a(i)=dcmplx(c1,c2)
           nc=nc+a(i)*conjg(a(i))
        end do
        a=a/dsqrt(nc)
        rho=dcmplx(0.0d0,0.0d0)
        j=1
        do i=1,4
           rho(i,i)=a(j)*conjg(a(j))+a(j+1)*conjg(a(j+1))+a(j+2)*conjg(a(j+2))
           j=j+3
        end do
        dum=dcmplx(0.0d0,0.0d0)
        dum(1,2)=conjg(a(1))*a(4)+conjg(a(2))*a(5)+conjg(a(3))*a(6)
        dum(1,3)=conjg(a(1))*a(7)+conjg(a(2))*a(8)+conjg(a(3))*a(9)
        dum(1,4)=conjg(a(1))*a(10)+conjg(a(2))*a(11)+conjg(a(3))*a(12)
        dum(2,3)=conjg(a(4))*a(7)+conjg(a(5))*a(8)+conjg(a(6))*a(9)
        dum(2,4)=conjg(a(4))*a(10)+conjg(a(5))*a(11)+conjg(a(6))*a(12)
        dum(3,4)=conjg(a(7))*a(10)+conjg(a(8))*a(11)+conjg(a(9))*a(12)
        rho=rho+dum+conjg(transpose(dum))
end subroutine q2r3
!




!random 2 qubit state of rank 4
!swansea, 04-2017
subroutine q2r4(idum,rho) 
        implicit none
        integer::r(1:2),idum
        complex*16::psi(1:16,1:1),rho(1:4,1:4)
        r(1)=3
        r(2)=4
        call randpsi(idum,4,2,psi)
        call pptrace(psi,4,2,r,2,rho)
end subroutine q2r4
!



!random 3 qubit state of rank 2
subroutine q3r2(idum,rho) 
        implicit none
        integer::r(1:1),idum
        complex*16::psi(1:16,1:1),rho(1:8,1:8)
        r(1)=4
        call randpsi(idum,4,2,psi)
        call pptrace(psi,4,2,r,1,rho)
end subroutine q3r2
!




!random 3 qubit state of rank 4
!swansea, 04-2017
subroutine q3r4(idum,rho) 
        implicit none
        integer::r(1:2),idum
        complex*16::psi(1:32,1:1),rho(1:8,1:8)
        r(1)=4
        r(2)=5
        call randpsi(idum,5,2,psi)
        call pptrace(psi,5,2,r,2,rho)
end subroutine q3r4
!



!n qubit ghz state
!swansea, 04-2017
subroutine ghz(n,psi) 
        implicit none 
        integer::n
        complex*16::psi(1:2**n,1:1)
        psi=dcmplx(0.0d0,0.0d0)
        psi(1,1)=dcmplx(1.0d0,0.0d0)
        psi(2**n,1)=dcmplx(1.0d0,0.0d0)
        psi=(1.0d0/dsqrt(2.0d0))*psi
end subroutine ghz
!




!n qubit random generalized ghz state
!swansea, 04-2017
subroutine genghzrand(idum,n,psi) 
        implicit none 
        integer::n,idum
        complex*16::psi(1:2**n,1:1)
        double precision::x,gasdev
        psi=dcmplx(0.0d0)
        psi(1,1)=dcmplx(gasdev(idum),0.0d0)
        psi(2**n,1)=dcmplx(gasdev(idum),0.0d0)
        x=psi(1,1)*conjg(psi(1,1))+psi(2**n,1)*conjg(psi(2**n,1))
        psi=(1.0d0/dsqrt(x))*psi
end subroutine genghzrand
!



!n qubit ghz state in x basis
!swansea, 04-2017
subroutine ghzx(n,psi) 
        implicit none 
        integer::i,j,n
        complex*16::x(1:2,1:1,1:2),z(1:2,1:1),o(1:2,1:1),psi(1:2**n,1:1)
        complex*16,dimension(:,:),allocatable::c1,c2
        call cbase(z,o)
        x(:,:,1)=(z+o)/dsqrt(2.0d0)
        x(:,:,2)=(z-o)/dsqrt(2.0d0)
        psi=dcmplx(0.0d0,0.0d0)
        do j=1,2
           do i=1,n
              if (i==1) then
                 allocate(c1(1:2,1:1))
                 c1=x(:,:,j)
              else
                 allocate(c2(1:2**i,1:1))
                 call kron(c1,2**(i-1),1,x(:,:,j),2,1,c2)
                 deallocate(c1)
                 allocate(c1(1:2**i,1:1))
                 c1=c2
                 deallocate(c2)
              end if
           end do
           psi=psi+c1/dsqrt(2.0d0)
           deallocate(c1)
        end do
end subroutine ghzx
!




!n-qubit w state
!swansea, 04-2017
subroutine ws(n,psi) 
        implicit none
        integer::n,b(1:n),d,i
        double precision::fact
        complex*16::psi(1:2**n,1:1)
        fact=1.0d0/dsqrt(dfloat(n))
        psi=dcmplx(0.0d0,0.0d0)
        b=0
        do i=1,n
           b(i)=1
           call b2d(b,n,2,d)
           psi(d,1)=dcmplx(fact,0.0d0)
           b(i)=0
        end do
end subroutine ws
!




!n qubit random generalized w state
!swansea, 04-2017
subroutine genwrand(idum,n,psi)
        implicit none 
        integer::n,idum,b(1:n),d,i
        double precision::fact,gasdev
        complex*16::psi(2**n,1:1)
        psi=dcmplx(0.0d0,0.0d0)
        fact=0.0d0
        b=0
        do i=1,n
           b(i)=1
           call b2d(b,n,2,d)
           psi(d,1)=dcmplx(gasdev(idum),0.0d0)
           fact=fact+psi(d,1)*conjg(psi(d,1))
           b(i)=0
        end do
        psi=(1.0d0/dsqrt(fact))*psi
end subroutine genwrand
!



!random w class states
!swansea, 04-2017
subroutine wclsrand(idum,s) 
        implicit none
        integer::i,idum
        double precision::a,b,nc,gasdev
        complex*16::s(1:8,1:1)
        s=dcmplx(0.d0,0.d0)
        s(1,1)=dcmplx(gasdev(idum),gasdev(idum))
        s(2,1)=dcmplx(gasdev(idum),gasdev(idum))
        s(3,1)=dcmplx(gasdev(idum),gasdev(idum))
        s(5,1)=dcmplx(gasdev(idum),gasdev(idum))
        nc=0.d0
        do i=1,8
           nc=nc+s(i,1)*conjg(s(i,1))
        end do
        s=s/dsqrt(nc)
end subroutine wclsrand
!




!nine parametric classes of four qubit states
!swansea, 04-2017
subroutine fourq(n,idum,psi)
        implicit none 
        integer::i,j,k,idum,n
        double precision::fact,gasdev,x,y,z
        complex*16::cf(1:4),psi(1:16,1:1)
        if (n/=7.or.n/=8.or.n/=9) then
           x=1.0d0/dsqrt(2.0d0)
           do i=1,4
              y=gasdev(idum)
              z=gasdev(idum)
              cf(i)=dcmplx(dabs(y),z)
           end do
           if (n==1) then
              psi=dcmplx(0.0d0,0.0d0)
              fact=0.0d0
              do i=1,16
                 if (i==1 .or. i==16) psi(i,1)=(cf(1)+cf(2))/2.0d0
                 if (i==4 .or. i==13) psi(i,1)=(cf(1)-cf(2))/2.0d0
                 if (i==6 .or. i==11) psi(i,1)=(cf(3)+cf(4))/2.0d0
                 if (i==7 .or. i==10) psi(i,1)=(cf(3)-cf(4))/2.0d0
                 fact=fact+psi(i,1)*conjg(psi(i,1))
              end do
              psi=(1.0d0/dsqrt(fact))*psi
           end if
           if (n==2) then
              psi=dcmplx(0.0d0,0.0d0)
              fact=0.0d0
              do i=1,16
                 if (i==1 .or. i==16) psi(i,1)=(cf(1)+cf(2))/2.0d0
                 if (i==4 .or. i==13) psi(i,1)=(cf(1)-cf(2))/2.0d0
                 if (i==6 .or. i==11) psi(i,1)=cf(3)
                 if (i==7) psi(i,1)=dcmplx(1.0d0,0.0d0)
                 fact=fact+psi(i,1)*conjg(psi(i,1))
              end do
              psi=(1.0d0/dsqrt(fact))*psi
           end if
           if (n==3) then
              psi=dcmplx(0.0d0,0.0d0)
              fact=0.0d0
              do i=1,16
                 if (i==1 .or. i==16) psi(i,1)=cf(1)
                 if (i==6 .or. i==11) psi(i,1)=cf(2)
                 if (i==4 .or. i==7) psi(i,1)=dcmplx(1.0d0,0.0d0)
                 fact=fact+psi(i,1)*conjg(psi(i,1))
              end do
              psi=(1.0d0/dsqrt(fact))*psi
           end if
           if (n==4) then
              psi=dcmplx(0.0d0,0.0d0)
              fact=0.0d0
              do i=1,16
                 if (i==1 .or. i==16) psi(i,1)=cf(1)
                 if (i==6 .or. i==11) psi(i,1)=(cf(1)+cf(2))/2.0d0
                 if (i==7 .or. i==10) psi(i,1)=(cf(1)-cf(2))/2.0d0
                 if (i==2 .or. i==3 .or. i==8 .or. i==12) psi(i,1)=dcmplx(0.0d0,x)
                 fact=fact+psi(i,1)*conjg(psi(i,1))
              end do
              psi=(1.0d0/dsqrt(fact))*psi
           end if
           if (n==5) then
              psi=dcmplx(0.0d0,0.0d0)
              fact=0.0d0
              do i=1,16
                 if (i==1 .or. i==6 .or. i==11 .or. i==16) psi(i,1)=cf(1)
                 if (i==2) psi(i,1)=dcmplx(0.0d0,1.0d0)
                 if (i==7) psi(i,1)=dcmplx(1.0d0,0.0d0)
                 if (i==12) psi(i,1)=dcmplx(0.0d0,-1.0d0)
                 fact=fact+psi(i,1)*conjg(psi(i,1))
              end do
              psi=(1.0d0/dsqrt(fact))*psi
           end if
           if (n==6) then
              psi=dcmplx(0.0d0,0.0d0)
              fact=0.0d0
              do i=1,16
                 if (i==1 .or. i==16) psi(i,1)=cf(1)
                 if (i==4 .or. i==6 .or. i==7) psi(i,1)=dcmplx(1.0d0,0.0d0)
                 fact=fact+psi(i,1)*conjg(psi(i,1))
              end do
              psi=(1.0d0/dsqrt(fact))*psi
           end if
        else
           if (n==7) then
              psi=dcmplx(0.0d0,0.0d0)
              do i=1,16
                 if (i==1 .or. i==6 .or. i==9 .or. i==15) psi(i,1)=dcmplx(5.0d-1,0.0d0)
              end do
           end if
           if (n==8) then
              psi=dcmplx(0.0d0,0.0d0)
              do i=1,16
                 if (i==1 .or. i==12 .or. i==14 .or. i==15) psi(i,1)=dcmplx(5.0d-1,0.0d0)
              end do           
           end if
           if (n==9) then
              x=1.0d0/dsqrt(2.0d0)
              psi=dcmplx(0.0d0,0.0d0)
              do i=1,16
                 if (i==1 .or. i==8) psi(i,1)=dcmplx(x,0.0d0)
              end do
           end if
        end if
end subroutine fourq
!
!TOOLS


!integer transform: other base to decimal
!swansea, 04-2017
subroutine b2d(b,n,m,d)
        implicit none
        integer::n,b(1:n),i,d,m
        d=1
        do i=n,1,-1
           d=d+b(i)*m**(n-i)
        end do
end subroutine b2d
!



!integer transform: decimal to other base
!swansea, 04-2017
subroutine d2b(d,n,m,b)
        implicit none
        integer::n,d,i,dum,q,r,b(1:n),m
        dum=d
        do i=n,1,-1
           q=dum/m
           r=dum-m*q
           b(i)=r
           dum=q
        end do
end subroutine d2b
!



!extension of basis (for localizable quantum correlations)
!swansea, 04-2017
subroutine extbase(n,b,pos,m,basis) 
        implicit none
        integer::n,m,pos(1:m),b(1:m),i,j,k,l,flag,basis(1:n)
        do j=1,n
           flag=0
           do k=1,m
              if (j==pos(k)) then
                 flag=k
                 exit
              end if
           end do
           if (flag==0) then
              basis(j)=-1
           else
              basis(j)=b(k)
           end if
        end do
end subroutine extbase
!



!complementary set in a set of integers
!swansea, 04-2017
subroutine compl(n,m,dpos,cpos)
        implicit none
        integer::m,n,i,j,k,l,dpos(1:m),cpos(1:n-m)
        k=0
        do i=1,n
           l=0
           do j=1,m
              if (dpos(j)==i) then
                 l=1
                 exit
              end if
           end do
           if (l==0) then
              k=k+1
              cpos(k)=i
           end if
        end do
end subroutine compl
!


!complex exponential of double-precision real variable
!swansea, 04-2017
function dcexp(x) 
        implicit none
        double precision::x
        complex*16::dcexp
        dcexp=dcmplx(dcos(x),dsin(x))
        return
        end

        
        

!logarithm of a double precision number in base 2
!swansea, 04-2017
        function dlog2(x)
        implicit none
        double precision::x,dlog2
        dlog2=dlog(x)/dlog(2.0d0)
        return
        end




!sorting array of integers in ascending order
!swansea, 04-2017
subroutine asort(a,n,b)  
        implicit none 
        integer::n,a(1:n),b(1:n),temp,i,j,k,l
        b=a
        do i=1,n
           do j=i,n
              if (b(i)>b(j)) then
                 temp=b(j)
                 b(j)=b(i)
                 b(i)=temp
              end if
           end do
        end do
end subroutine asort
!



!sorting array of integers in descending order
!swansea, 04-2017
subroutine dsort(a,n,b)  
        implicit none
        integer::n,a(1:n),b(1:n),temp,i,j,k,l
        b=a
        do i=1,n
           do j=i,n
              if (b(i)<b(j)) then
                 temp=b(j)
                 b(j)=b(i)
                 b(i)=temp
              end if
           end do
        end do
end subroutine dsort
!


!factorial
!swansea, 04-2017
subroutine fact(ni,no) 
        implicit none
        integer::ni,no,i
        if (ni==0) then
           no=1
        else
           no=1
           do i=1,ni
              no=no*i
           end do
        end if
end subroutine fact
!



!combination of large numbers
!swansea, 04-2017
subroutine comb_l(n1,n2,no)
        implicit none
        integer::i,j,k,l,n1,n2,m1,m2,no,no1,no2
        integer::pn(1:5)
        integer,dimension(:),allocatable :: num,den
        pn(1)=2
        pn(2)=3
        pn(3)=5
        pn(4)=7
        pn(5)=11
        if (n1>=n2) then
           if (n2>=n1-n2) then
             m1=n2+1
              m2=n1-n2
           else
              m1=n1-n2+1
              m2=n2
           end if
           allocate(num(n1-m1+1))
           allocate(den(m2))
           j=0
           do i=m1,n1
              j=j+1
              num(j)=i
           end do
           do i=1,m2
              den(i)=i
              k=1
              l=0
              do  
                 l=l+1
                 if(l>5) exit
                 do j=1,n1-m1+1
                    if (num(j)==pn(l)*i) then 
                       num(j)=pn(l)
                       den(i)=1
                       k=0
                    end if
                 end do
                 if (k==1) then
                    continue
                 else
                    exit
                 end if
              end do
           end do
           no1=1
           do i=1,n1-m1+1
              no1=no1*num(i)
           end do
           no2=1
           do i=1,m2
              no2=no2*den(i)
           end do
           deallocate(num,den)
           no=no1/no2
        end if
end subroutine comb_l
!



!combination of moderately large numbers
!swansea, 04-2017
subroutine comb_m(n1,n2,no) 
        implicit none
        integer::i,n1,n2,num,den,m1,m2,no
        if (n1>=n2) then
           if (n2>=n1-n2) then
              m1=n2+1
              m2=n1-n2
           else
              m1=n1-n2+1
              m2=n2
           end if
           num=1
           do i=m1,n1
              num=num*i
           end do
           den=1
           do i=1,m2
              den=den*i
           end do
           no=num/den
        end if
end subroutine comb_m
!



!combination of integers
!swansea, 04-2017
subroutine comb_s(n1,n2,no)
        implicit none
        integer::n1,n2,no,fn1,fn2,fn3
        if (n1>=n2) then
           call fact(n1,fn1)
           call fact(n2,fn2)
           call fact(n1-n2,fn3)
        end if
        no=fn1/(fn2*fn3)
end subroutine comb_s
!




!dividing 1 into n random fractions
!swansea, 04-2017
subroutine prand(idum,n,p)
        implicit none
        integer::i,j,n,idum
        double precision::x,a,sp,p(1:n),ran2
        sp=0.0d0
        do j=1,n
           if (j/=n) then
              a=1.0d0-sp
              x=ran2(idum)
              p(j)=a*x
              sp=sp+p(j)
           else
              p(j)=1.0d0-sp
           end if
        end do
end subroutine prand
!
