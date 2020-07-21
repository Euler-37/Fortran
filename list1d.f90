module list
    implicit none
    real(8),parameter::eps=1.d-12
    interface init
        module procedure init_int
        module procedure init_real
    end interface
    interface append
        module procedure append_int
        module procedure append_real
    end interface
    interface pop
        module procedure pop_int
        module procedure pop_real
    end interface
    interface remove
        module procedure remove_int
        module procedure remove_real
    end interface
    interface reverse
        module procedure reverse_int
        module procedure reverse_real
    end interface
    interface sort
        module procedure sort_int
        module procedure sort_real
    end interface
    interface set
        module procedure set_int
        module procedure set_real
    end interface
    interface search
        module procedure search_int
        module procedure search_real
    end interface
    interface times
        module procedure times_int
        module procedure times_real
    end interface
contains

    subroutine init_int(a)
        integer,intent(inout),allocatable   :: a(:)
        a=pack([1],[1]/=1)
    end subroutine
    subroutine init_real(a)
        real(8),intent(inout),allocatable   :: a(:)
        a=pack([1.d0],[1.d0]>2.d0)
    end subroutine


    subroutine append_int(a,i,j)
        integer,intent(inout),allocatable   :: a(:)
        integer,intent(in)                  :: i
        integer,optional                    :: j
        integer                             :: sizea,k
        if(.not.allocated(a))then
            call init(a)
        end if
        sizea=size(a)
        if(.not.present(j))then
            k=sizea+1
        else
            k=j
        end if
        a=[a(1:k-1),i,a(k:sizea)]
    end subroutine
    subroutine append_real(a,i,j)
        real(8),intent(inout),allocatable   :: a(:)
        real(8),intent(in)                  :: i
        integer,optional                    :: j
        integer                             :: sizea,k
        if(.not.allocated(a))then
            call init(a)
        end if
        sizea=size(a)
        if(.not.present(j))then
            k=sizea+1
        else
            k=j
        end if
        a=[a(1:k-1),i,a(k:sizea)]
    end subroutine

    subroutine pop_int(a,i)
        integer,optional                    :: i
        integer,intent(inout),allocatable   :: a(:)
        integer                             :: sizea
        sizea=size(a)
        if(present(i))then
            if(i>sizea.or.i<1)then
                write(*,"(A,I5,A)")"Index Error,Index",i," does not exist"
            else
                a=[a(1:i-1),a(i+1:sizea)]
            end if
        else
            a=[a(1:sizea-1)]
        end if
    end subroutine
    subroutine pop_real(a,i)
        integer,optional                    :: i
        real(8),intent(inout),allocatable   :: a(:)
        integer                             :: sizea
        sizea=size(a)
        if(present(i))then
            if(i>sizea.or.i<1)then
                write(*,"(*(g0))")"Index Error,Index",i," does not exist"
            else
                a=[a(1:i-1),a(i+1:sizea)]
            end if
        else
            a=[a(1:sizea-1)]
        end if
    end subroutine

    subroutine remove_int(a,i)
        integer,intent(inout),allocatable   :: a(:)
        integer,intent(in)                  :: i
        a=pack(a,a/=i)
    end subroutine
    subroutine remove_real(a,i)
        real(8),intent(inout),allocatable   :: a(:)
        real(8),intent(in)                  :: i
        a=pack(a,abs(a-i)<eps)
    end subroutine


    subroutine sort_int(a)
        integer,intent(inout)               :: a(:)
        a=qsort_int(a)
    end subroutine
    subroutine sort_real(a)
        real(8),intent(inout)               :: a(:)
        a=qsort_real(a)
    end subroutine

    recursive function qsort_int(a)result(sort)
        integer,intent(in)                  :: a(:)
        integer                             :: sort(size(a))
        if(size(a)>1)then
            sort=[qsort_int(pack(a(2:),a(2:)<=a(1))),&
                &a(1),&
                &qsort_int(pack(a(2:),a(2:)>a(1)))]
        else
            sort=a
        end if
    end function
    recursive function qsort_real(a)result(sort)
        real(8),intent(in)                  :: a(:)
        real(8)                             :: sort(size(a))
        if(size(a)>1)then
            sort=[qsort_real(pack(a(2:),a(2:)<a(1))),&
                &a(1),&
                &qsort_real(pack(a(2:),a(2:)>a(1)))]
        else
            sort=a
        end if
    end function

    subroutine reverse_int(a)
        integer,intent(inout)               :: a(:)
        a=a(size(a):1:-1)
    end subroutine
    subroutine reverse_real(a)
        real(8),intent(inout)               :: a(:)
        a=a(size(a):1:-1)
    end subroutine

    subroutine set_int(a,seq)
        integer,intent(inout),allocatable   :: a(:)
        integer                             :: i,temp,j,sizea
        logical,optional                    :: seq
        logical                             :: judge
        if(.not.present(seq))then
            judge=.true.
        else
            judge=seq
        end if
        sizea=size(a)
        if(judge)then
            call sort(a)
            temp=a(1)
            j=0
            do i=2,sizea
                if(a(i-j)==temp)then
                    call pop(a,i-j)
                    j=j+1
                else
                    temp=a(i-j)
                end if
            end do
        else
            temp=a(1)
            i=a(1)
            do
                call remove(a,i)
                call append(a,i)
                i=a(1)
                if(temp==i)exit
            end do
        end if
    end subroutine

    subroutine set_real(a,seq)
        real(8),intent(inout),allocatable   :: a(:)
        integer                             :: i,j,sizea
        real(8)                             :: temp,temp1
        logical,optional                    :: seq
        logical                             :: judge
        if(.not.present(seq))then
            judge=.true.
        else
            judge=seq
        end if
        sizea=size(a)
        if(judge)then
            call sort(a)
            temp=a(1)
            j=0
            do i=2,sizea
                if(a(i-j)==temp)then
                    call pop(a,i-j)
                    j=j+1
                else
                    temp=a(i-j)
                end if
            end do
        else
            temp=a(1)
            temp1=a(1)
            do
                call remove(a,temp1)
                call append(a,temp1)
                temp1=a(1)
                if(abs(temp-temp1)<eps)exit
            end do
        end if
    end subroutine


    function search_int(a,i)
        integer,intent(in)                  :: i
        integer,intent(in)                  :: a(:)
        integer                             :: sizea,search_int
        search_int=1
        sizea=size(a)
        do while(search_int/=sizea+1)
            if(a(search_int)==i)exit
            search_int=search_int+1
        end do
        if(search_int==sizea+1)then
            !write(*,"(A,I4,A)")"Value ",i," does not exist"
            search_int=-1
        end if
    end function
    function search_real(a,i)
        real(8),intent(in)                  :: i
        real(8),intent(in)                  :: a(:)
        integer                             :: sizea,search_real
        search_real=1
        sizea=size(a)
        do while(search_real/=sizea+1)
            if(abs(a(search_real)-i)<eps)exit
            search_real=search_real+1
        end do
        if(search_real==sizea+1)then
            !write(*,"(A,I4,A)")"Value ",i," does not exist"
            search_real=-1
        end if
    end function

    function times_int(a,i)
        integer,intent(in)                  :: i
        integer,intent(in)                  :: a(:)
        integer,allocatable                 :: b(:)
        integer                             :: times_int
        b=a
        call remove(b,i)
        times_int=size(a)-size(b)
    end function

    function times_real(a,i)
        real(8),intent(in)                  :: i
        real(8),intent(in)                  :: a(:)
        real(8),allocatable                 :: b(:)
        integer                             :: times_real
        b=a
        call remove(b,i)
        times_real=size(a)-size(b)
    end function

end module

program main
    use list
    implicit none
    integer,allocatable                     :: a(:)
    real(8),allocatable                     :: b(:)
    integer                                 :: i
    real(8)                                 :: r1
    call random_seed()
    call init(a)
    do i=1,10
        call random_number(r1)
        call append(a,int(r1*10))
    end do
    write(*,"(10I4)")a
    write(*,*)search(a,3),times(a,3)
    call init(b)
    do i=1,10
        call random_number(r1)
        call append(b,r1*10.d0)
    end do
    write(*,"(10F7.3)")b
    call sort(b)
    write(*,"(10F7.3)")b
end program
