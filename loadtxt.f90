module readfile
    contains
    function loadtxt_int(filename) result(mat)
        implicit none
        integer,parameter::maxdim=1000000
        integer,parameter::initint=huge(0)
        character(1)     ::judge
        character(20)::filename
        integer::i,j,stat,row,col,matdim
        integer,allocatable::b(:),mat(:,:)
        open(10,file=filename)
        allocate(b(MAXDIM))
        b=initint
        read(10,fmt=*,iostat=stat)(b(j),j=1,maxdim)
        close(10)
        b=pack(b,b/=initint)
        matdim=size(b)
        open(10,file=filename)
        row=0
        do
            read(10,fmt=*,iostat=stat)judge
            if(stat/=0)exit
            row=row+1
        end do
        close(10)
        col=matdim/row
        allocate(mat(col,row))
        mat=reshape(b,[col,row])
    end function
    function loadtxt_real8(filename) result(mat)
        implicit none
        integer,parameter::maxdim=1000000
        real(8),parameter::initreal=huge(0.d0)
        character(1)     ::judge
        character(20)    ::filename
        integer::i,j,stat,row,col,matdim
        real(8),allocatable::b(:),mat(:,:)
        open(10,file=filename)
        allocate(b(MAXDIM))
        b=initreal
        read(10,fmt=*,iostat=stat)(b(j),j=1,maxdim)
        close(10)
        b=pack(b,b<initreal)
        matdim=size(b)
        open(10,file=filename)
        row=0
        do
            read(10,fmt=*,iostat=stat)judge
            if(stat/=0)exit
            row=row+1
        end do
        close(10)
        col=matdim/row
        allocate(mat(col,row))
        mat=reshape(b,[col,row])
    end function

end module

program main
    use readfile
    integer,allocatable::a(:,:)
    integer            ::shape1(2),i
    character(20)::filename
    filename="sample.txt"
    a=loadtxt_int(filename)
    shape1=shape(a)
    write(*,*)shape1
    do i=1,shape1(2)
        write(*,"(*(I3))")a(:,i)
    end do
end program
