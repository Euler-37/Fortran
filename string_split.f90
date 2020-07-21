program main
    implicit none
    character(80)::string
    character(80),allocatable::str_array(:)
    integer::i
    read(*,"(A)")string
    call split(string,str_array)
    do i=1,size(str_array)
        write(*,"(g0)")trim(str_array(i))
    end do
contains
    subroutine split(string,str_array)
        implicit none
        character(80)::string
        character(80),allocatable::str_array(:)
        integer::i,ii,lens
        lens=len_trim(string)
        ii=1
        str_array=pack(["1"],["1"]/="1")
        do i=1,lens
            if(i/=lens)then
                if(string(i:i)==','.or.string(i:i)==" ")then
                    str_array=[str_array,string(ii:i-1)]
                    ii=i+1
                end if
            else
                str_array=[str_array,string(ii:lens)]
            end if
        end do
    end subroutine
end program


