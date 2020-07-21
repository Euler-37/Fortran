program main
    implicit none
    integer,parameter::digit=128
    character(len=digit)::a,b,temp
    character(len=digit+1)::c
    integer i,j,len_a,len_b,dab,len_temp,add1,add2,a1,b1
    do i=1,digit
        a(i:i)=" "
    end do
    do i=1,digit
        b(i:i)=" "
    end do
    do i=1,digit+1
        c(i:i)=" "
    end do
    write(*,*)"****************************************"
    write(*,*)"*              big int add             *"
    write(*,*)"****************************************"
    write(*,*)"*enter two numbers less than 128 digits*"
    write(*,*)"****************************************"
    read(*,*)a
    read(*,*)b
    write(*,*)"****************************************"
    write(*,*)"            The numbers are             "
    write(*,*)"****************************************"
    write(*,*)trim(a)
    write(*,*)trim(b)
    write(*,*)"****************************************"
    len_a=len_trim(a)
    len_b=len_trim(b)
    dab=len_a-len_b
    if(len_b>len_a)then
        temp=a
        a=b
        b=temp
        len_temp=len_a
        len_a=len_b
        len_b=len_temp
    end if
    j=0
    add2=0
    do i=len_b,1,-1
        read(a(dab+i:dab+i),"(I1)")a1
        read(b(i:i),"(I1)")b1
        add1=mod(a1+b1+add2,10)
        add2=(a1+b1+add2)/10
        write(c(digit+1-j:digit+1-j),"(I1)")add1
        j=j+1
    end do
    do i=len_a-len_b,1,-1
        read(a(i:i),"(I1)")a1
        add1=mod(add2+a1,10)
        add2=(add2+a1)/10
        write(c(digit+1-j:digit+1-j),"(I1)")add1
        j=j+1
    end do
    if(add2/=0)then
        write(c(digit+1-len_a:digit+1-len_a),"(I1)")add2
    end if
    write(*,*)"****************************************"
    print*,adjustl(c)
    write(*,*)"****************************************"
end program
