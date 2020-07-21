module sparse_matrix
    !sparse matrix mod
    implicit none
    !+--------------------------------------------------------------------------+!
    !|Give a suitable number of [spn],it is better for spn bigger than or equal |!
    !|to the number of matrix elements. This term is SO IMPORTANT, it can       |!
    !|influence code speed directly and significantly.                          |!
    !+--------------------------------------------------------------------------+!
    integer,parameter::spn=1002500,addn=1000
    real(8),parameter::eps=1.d-12
    type sparse
        integer::dimen
        integer::num
        integer,allocatable::row(:)
        integer,allocatable::col(:)
        real(8),allocatable::val(:)
    end type
contains
    subroutine init(a,dima)
        type(sparse),intent(inout)  ::a
        integer,intent(in)          ::dima
        a%dimen=dima
        a%num=0
        allocate(a%row(spn))
        allocate(a%col(spn))
        allocate(a%val(spn))
        a%row=0
        a%col=0
        a%val=0.d0
    end subroutine
    subroutine reinit(a)
        type(sparse),intent(inout)::a
        integer::temp_int(addn)
        real(8)::temp_real(addn)
        temp_int=0
        temp_real=0.d0
        a%row=[a%row,temp_int]
        a%col=[a%col,temp_int]
        a%val=[a%val,temp_real]
    end subroutine
    subroutine cut(a)
        type(sparse),intent(inout)::a
        a%row=[a%row(1:a%num)]
        a%col=[a%col(1:a%num)]
        a%val=[a%val(1:a%num)]
    end subroutine
    subroutine append(a,i,j,r)
        type(sparse),intent(inout)  ::a
        integer,intent(in)          ::i,j
        real(8),intent(in)          ::r
        integer                     ::idx
        logical                     ::judge
        judge=abs(r)>eps
        idx=sparse_search(a,i,j)
        if(idx<0 .and. judge)then
            !the row and col not exist
            idx=abs(idx)
            if(a%num+1>size(a%row))then
                !change the matrix dimension
                call reinit(a)
            end if
            a%row(idx+1:a%num+1)=a%row(idx:a%num)
            a%row(idx)=i
            a%col(idx+1:a%num+1)=a%col(idx:a%num)
            a%col(idx)=j
            a%val(idx+1:a%num+1)=a%val(idx:a%num)
            a%val(idx)=r
            a%num=a%num+1
        elseif(idx>0 .and. judge)then
            !col and row exist,change the value
            a%val(idx)=r
        elseif(idx>0 .and. (.not. judge))then
            !changed value to 0,pop it
            a%row(idx:a%num-1)=a%row(idx+1:a%num)
            a%row(a%num)=0
            a%col(idx:a%num-1)=a%col(idx+1:a%num)
            a%col(a%num)=0
            a%val(idx:a%num-1)=a%val(idx+1:a%num)
            a%val(a%num)=0.d0
            a%num=a%num-1
        end if
    end subroutine

    function sparse_search(a,i,j)result(idx)
        type(sparse),intent(in)  :: a
        integer,intent(in)          :: i,j
        integer                     :: idx,sizea,k1,k2,a1,b1,c1,l1,n
        sizea=a%num
        if(sizea==0)then
            idx=-1
            return
        end if
        k1=1
        k2=sizea
        n=(a%dimen)+1
        c1=(i-1)*n+j
        a1=(a%row(k1)-1)*n+a%col(k1)
        b1=(a%row(k2)-1)*n+a%col(k2)
        if(a1==c1)then
            idx=1
        elseif(b1==c1)then
            idx=sizea
        elseif(c1>b1)then
            idx=-(sizea+1)
        elseif(c1<a1)then
            idx=-1
        else
            do
                l1=(k1+k2)/2
                if(l1==k1)then
                    idx=-(l1+1)
                    exit
                end if
                a1=(a%row(l1)-1)*n+a%col(l1)
                if(a1==c1)then
                    idx=l1
                    exit
                elseif(a1>c1)then
                    k2=l1
                elseif(a1<c1)then
                    k1=l1
                end if
            end do
        end if
    end function

    function sparse_val(a,i,j)result(sv)
        type(sparse),intent(in)::a
        integer,intent(in)::i,j
        real(8)::sv
        integer::idx
        idx=sparse_search(a,i,j)
        if(idx<0)then
            sv=0.d0
        else
            sv=a%val(idx)
        end if
    end function

    function sparse_product(a,b)result(c)
        type(sparse),intent(in)     :: a
        real(8),intent(in)          :: b(:)
        integer                     :: sizea,dimb,i
        real(8),allocatable         :: c(:)
        sizea=a%num
        dimb=size(b)
        if(dimb/=a%dimen)then
            write(*,"(A,I8,A,I8)")"Error: Different shape for array,",dimb,"and",a%dimen
        else
            allocate(c(dimb))
            c=0.d0
            do i=1,sizea
                c(a%row(i))=c(a%row(i))+a%val(i)*b(a%col(i))
            end do
        end if
    end function

    function toarray(a)result(b)
        type(sparse),intent(in)     :: a
        real(8),allocatable         :: b(:,:)
        integer                     :: i
        allocate(b(a%dimen,a%dimen))
        b=0.d0
        do concurrent(i=1:size(a%row))
            b(a%row(i),a%col(i))=a%val(i)
        end do
    end function

    subroutine lanczos(a,eig)
        integer,parameter       ::ndim=400
        type(sparse),intent(in) ::a
        real(8),intent(out)     ::eig
        integer                 ::sizea,k
        real(8),external        ::dnrm2,ddot
        real(8)                 ::tri_diag(ndim),tri_off(ndim),eigv(ndim),vnorm,a0,b1
        real(8),dimension(:),pointer::pv0,pv1,ptemp
        real(8),allocatable,target     ::v0(:)
        real(8),allocatable,target     ::v1(:)
        sizea=a%dimen
        allocate(v1(sizea))
        allocate(v0(sizea))
        call random_number(v0)
        tri_diag=0.d0
        tri_off =0.d0
        vnorm=dnrm2(sizea,v0,1)
        v0=v0/vnorm
        v1=sparse_product(a,v0)
        k=1
        pv0=>v0
        pv1=>v1
        loop1:do
            write(*,*)">>>",k
            a0=ddot(sizea,pv0,1,pv1,1)
            tri_diag(k)=a0
            if(k>2)then
                call tqli(tri_diag(1:k),tri_off(1:k),k,eig,eigv(1:k))
                if(abs(eigv(k))*b1<eps)then
                    write(*,*)"Lanczos Finished"
                    exit loop1
                end if
            end if
            pv1=-a0*pv0+pv1
            b1=norm2(pv1)
            tri_off(k)=b1
            pv1=pv1/b1
            pv0=-pv0*b1

            ptemp=>pv0
            pv0=>pv1
            pv1=>ptemp

            pv1=sparse_product(a,pv0)+pv1
            k=k+1
        end do loop1
        nullify(pv0,pv1,ptemp)
    end subroutine

    subroutine tqli(d1,e1,n,ev,v)
        integer n
        real(8) ::d(n),e(n),z(n,n),v(n),d1(n),e1(n),ev
        integer ::i,iter,k,l,m
        real(8) ::b,c,dd,f,g,p,r,s,one
        d=d1
        e=e1
        v=0.d0
        z=0.d0
        do concurrent(i=1:n)
            z(i,i)=1.d0
        end do

        one=1.d0
        do l=1,n
            iter=0
            iterate : do
                do m=l,n-1
                    dd=abs(d(m))+abs(d(m+1))
                    if (abs(e(m))<eps) exit
                enddo
                if(m == l) exit iterate
                if(iter == 30) stop 'too many iterations in tqli'
                iter=iter+1
                g=(d(l+1)-d(l))/(2.d0*e(l))
                r=pythag(g,one)
                g=d(m)-d(l)+e(l)/(g+sign(r,g))
                s=1.d0
                c=1.d0
                p=0.d0
                do i=m-1,l,-1
                    f=s*e(i)
                    b=c*e(i)
                    r=pythag(f,g)
                    e(i+1)=r
                    if(r<eps) then
                        d(i+1)=d(i+1)-p
                        e(m)=0.d0
                        cycle iterate
                    endif
                    s=f/r
                    c=g/r
                    g=d(i+1)-p
                    r=(d(i)-g)*s+2.d0*c*b
                    p=s*r
                    d(i+1)=g+p
                    g=c*r-b
                    !     omit lines from here ...
                    do k=1,n
                        f=z(k,i+1)
                        z(k,i+1)=s*z(k,i)+c*f
                        z(k,i)=c*z(k,i)-s*f
                    enddo
                    !     ... to here when finding only eigenvalues.
                enddo
                d(l)=d(l)-p
                e(l)=g
                e(m)=0.d0
            end do iterate
        end do
        call eigsrt(d,z,n)
        ev=d(n)
        v=z(:,n)
    end subroutine tqli
    function pythag(a,b)
        real(8),intent(in)::a,b
        real(8)::absa,absb,pythag
        absa=abs(a)
        absb=abs(b)
        if(absa > absb) then
            pythag=absa*sqrt(1.+(absb/absa)**2)
        else
            if(absb <eps) then
                pythag=0.d0
            else
                pythag=absb*sqrt(1.d0+(absa/absb)**2)
            end if
        end if
    end function

    subroutine eigsrt(d,v,n)

        integer n
        real(8) d(n),v(n,n)
        integer i,j,k
        real(8) p
        do i=1,n-1
            k=i
            p=d(i)
            do j=i+1,n
                if(d(j).ge.p)then
                    k=j
                    p=d(j)
                end if
            end do
            if(k.ne.i)then
                d(k)=d(i)
                d(i)=p
                do j=1,n
                    p=v(j,i)
                    v(j,i)=v(j,k)
                    v(j,k)=p
                end do
            end if
        end do
    end subroutine
end module


program main
    use sparse_matrix
    implicit none
    integer,parameter::n=1000
    type(sparse)::a
    real(8)     ::eig,t1,t2,t3
    integer     ::i,j
    call init(a,n)
    call cpu_time(t1)
    do i=1,n,4
        do j=1,n,6
            call append(a,i,j,1.d0)
        end do
    end do
    do i=1,n,4
        do j=1,n,6
            call append(a,j,i,1.d0)
        end do
    end do
    do i=1,n
        call append(a,i,i,-2.d0)
    end do
    write(*,*)a%num
    call cpu_time(t3)
    write(*,*)t3-t1
    !call cut(a)
    !do i=1,a%num
    !    write(*,*)a%row(i),a%col(i),a%val(i)
    !end do
    call lanczos(a,eig)
    call cpu_time(t2)
    write(*,*)eig
    write(*,*)t2-t3
end program
