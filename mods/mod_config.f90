module mod_config
    use iso_fortran_env, only: int32, real64
    implicit none

    private
    public :: Configuration

    type :: Configuration
        character(20) :: config_file
        real(real64) :: Re
        real(real64) :: Pe
        integer(int32) :: nk
        real(real64) :: dt
        integer(int32) :: nt
        integer(int32) :: nw
        real(real64) :: CFL
    end type Configuration

    interface Configuration
        module procedure :: constructor
    end interface Configuration
    
contains

    type(Configuration) function constructor(config_file) result(res)
        character(*), intent(in) :: config_file
        character(20) :: line
        res % config_file = config_file
        open(10, file=trim(config_file), action='read')
        read(10, *) line, res % Re
        read(10, *) line, res % Pe
        read(10, *) line, res % nk
        read(10, *) line, res % dt
        read(10, *) line, res % nt
        read(10, *) line, res % nw
        read(10, *) line, res % CFL
        close(10)
    end function constructor
    
end module mod_config
