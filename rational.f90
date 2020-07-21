module frac
    implicit none
    type div
        integer(8)::num!numerator
        integer(8)::den!denominator
    end type
    interface operator(+)
        module procedure add
    end interface
    interface operator(-)
        module procedure sub
    end interface
    interface operator(*)
        module procedure mult
    end interface
    interface operator(/)
        module procedure divi
    end interface
contains
    elemental function gcd(x)result(z)
        implicit none
        type(div),intent(in)::x
        integer(8)          ::temp,small,big,sign1
        type(div)           ::z
        if(x%num==0)then
            z=div(0,1)
        else
            sign1=sign(1_8,x%num*x%den)
            z=div(abs(x%num),abs(x%den))
            small=min(z%num,z%den)
            big=max(z%num,z%den)
            do while(small>1)
                temp=mod(big,small)
                if(temp==0)exit
                big=small
                small=temp
            end do
            z=div(sign1*z%num/small,abs(z%den/small))
        end if
        return
    end function
    elemental function add(x,y)result(z)
        implicit none
        type(div),intent(in)::x,y
        type(div)           ::z
        z%num=x%num*y%den+y%num*x%den
        z%den=x%den*y%den
        z=gcd(z)
        return
    end function
    elemental function sub(x,y)result(z)
        implicit none
        type(div),intent(in)::x,y
        type(div)           ::z

        z%num=x%num*y%den-y%num*x%den
        z%den=x%den*y%den
        z=gcd(z)
        return
    end function
    elemental function mult(x,y)result(z)
        implicit none
        type(div),intent(in)::x,y
        type(div)           ::z
        z%num=x%num*y%num
        z%den=x%den*y%den
        z=gcd(z)
        return
    end function
    elemental function divi(x,y)result(z)
        implicit none
        type(div),intent(in)::x,y
        type(div)           ::z
        z%num=x%num*y%den
        z%den=x%den*y%num
        z=gcd(z)
        return
    end function

    subroutine matrixoutput(a,num)
        implicit none
        integer                    ::n,m,i,j,num1
        type(div),intent(in)       ::a(:,:)
        integer,optional,intent(in)::num
        if(present(num))then
            num1=num
        else
            num1=6
        end if
        m=size(a,1)
        n=size(a,2)
        do i=1,n
            do j=1,m
                call output(a(j,i),num1)
            end do
            write(*,*)
        end do
        write(*,*)
    end subroutine

    function matrixproduct(a,b)result(c)
        implicit none
        integer ::m,l,n
        type(div),intent(in) ::a(:,:),b(:,:)
        type(div),allocatable::c(:,:)
        integer              ::i,j,k
        m=size(a,1)
        l=size(a,2)
        n=size(b,2)
        allocate(c(m,n))
        c=div(0,1)
        do i=1,n
            do j=1,m
                do k=1,l
                    c(j,i)=c(j,i)+a(j,k)*b(k,i)
                end do
            end do
        end do
        return
    end function

    function det(a1,n)result(det1)
        implicit none
        integer,intent(in)  ::n
        type(div),intent(in)::a1(n,n)
        integer             ::k,j
        type(div)           ::a(n,n),b(n),det1,b1,c(n)
        a=a1
        det1=div(1,1)
        do j=1,n
            if(a(j,j)%num==0)then
                k=j
                do while(.true.)
                    if(k==n)then
                        det1=div(0,1)
                        return
                    end if
                    k=k+1
                    if(a(k,j)%num/=0)then
                        b=a(:,k)
                        a(:,k)=a(:,j)
                        a(:,j)=b
                        exit
                    end if

                end do
            end if
            b1=a(j,j)
            det1=det1*b1
            a(:,j)=a(:,j)/b1
            b=a(:,j)
            do k=j+1,n
                b1=a(j,k)
                c=b1*b
                a(:,k)=a(:,k)-c
            end do
        end do
        return
    end function

    function inv(a1,n)result(inv_a)
        implicit none
        integer,intent(in)  ::n
        type(div),intent(in)::a1(n,n)
        type(div)           ::b(n),invb(n),c(n),b1,inv_a(n,n),a(n,n)
        integer             ::i,j,k
        a=a1
        inv_a=div(0,1)
        do i=1,n
            inv_a(i,i)=div(1,1)
        end do
        do j=1,n
            if(a(j,j)%num==0)then
                k=j
                do while(.true.)
                    if(k==n)then
                        print*,"Singular matrix"
                        return
                    end if
                    k=k+1
                    if(a(k,j)%num/=0)then
                        b=a(:,k)
                        a(:,k)=a(:,j)
                        a(:,j)=b
                        invb=inv_a(:,k)
                        inv_a(:,k)=inv_a(:,j)
                        inv_a(:,j)=invb
                        exit
                    end if
                end do
            end if
            b1=a(j,j)
            a(:,j)=a(:,j)/b1
            inv_a(:,j)=inv_a(:,j)/b1
            b=a(:,j)
            invb=inv_a(:,j)
            do k=j+1,n
                b1=a(j,k)
                c=b1*b
                a(:,k)=a(:,k)-c
                c=b1*invb
                inv_a(:,k)=inv_a(:,k)-c
            end do
        end do
        do j=n,1,-1
            invb=inv_a(:,j)
            do k=j-1,1,-1
                b1=a(j,k)
                a(j,k)=div(0,1)
                c=b1*invb
                inv_a(:,k)=inv_a(:,k)-c
            end do
        end do
    end function

    subroutine output(x1,num)
        implicit none
        type(div),intent(in)    ::x1
        type(div)               ::x
        integer                 ::s1
        character(20)           ::form1,form2
        character(2)            ::str1
        integer,optional,intent(in)      ::num
        if(present(num))then
            write(str1,"(I2)")num
        else
            str1=" 6"
        end if
        form1="(A1,I"//trim(str1)//","//trim(str1)//"X,1X,A1)"
        form2="(A1,I"//trim(str1)//",A1,I"//trim(str1)//",A1)"
        x=x1
        if(x%num==0)then
            write(*,form1,advance="no")"(",x%num,")"
        else
            s1=int(sign(1_8,x%num*x%den),4)
            x=div(abs(x%num),abs(x%den))
            x=gcd(x)
            if(x%den==1)then
                write(*,form1,advance="no")"(",s1*x%num,")"
            else
                write(*,form2,advance="no") "(",s1*x%num,"/",x%den,")"
            end if
        end if
    end subroutine
end module

program main
    use frac
    implicit none
    integer,parameter::n=8
    type(div)        ::a(n,n),inv_a(n,n)
    integer          ::i,j
    write(*,*)"======================================================================="
    write(*,*)"|        Exact Inverse and Determinant of Rational Matrix             |"
    write(*,*)"======================================================================="
    write(*,*)"|     1.  maximal numerator or denominator is",huge(1_8),"   |"
    write(*,*)"|     2.  Format of numbers output is a number,optional,default is 6  |"
    write(*,*)"|     3.  Format of martix output is a number, optional,default is 6  |"
    write(*,*)"======================================================================="
    a=div(0,0)
    do i=1,n-1
        a(i,i+1)=div(1,1)
    enddo
    do i=2,n
        a(i,i-1)=div(1,1)
    enddo
    write(*,*)"Initial Matrix is:"
    call matrixoutput(a,1)
    write(*,*)"Determinant of the matrix is:"
    call output(det(a,n),2)
    write(*,*)

    inv_a=inv(a,n)
    write(*,*)"Inverse of the matrix is:"
    call matrixoutput(inv_a,2)
    write(*,*)"Check A*A^-1:"
    call matrixoutput(matrixproduct(a,inv_a),1)
end program
