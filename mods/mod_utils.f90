module mod_utils
    use iso_fortran_env, only: int32, real64, stdout => output_unit
    use mod_mesh, only: Mesh
    use mod_config, only: Configuration
    implicit none

    private
    public :: &
    set_field, &
    set_surface, &
    set_elapsed, &
    set_duration, &
    set_data
    
contains

    subroutine set_field(msh)
        type(Mesh), intent(in) :: msh
        integer(int32) :: i
        open(10, file='data/cxy.txt') ! field data coordinates
        do i = 1, msh % n_cells
            write(10, *) msh % cells(i) % coordinates(1), &
            msh % cells(i) % coordinates(2)
        end do
        close(10)
    end subroutine set_field

    subroutine set_surface(msh)
        type(Mesh), intent(in) :: msh
        integer(int32) :: i, ni
        open(10, file='data/nxy.txt') ! surface data coordinates
        do i = 1, size(msh % surface_nodes)
            ni = msh % surface_nodes(i)
            write(10, *) msh % nodes(ni) % coordinates(1), &
            msh % nodes(ni) % coordinates(2)
        end do
        close(10)
    end subroutine set_surface
    
    subroutine set_elapsed(ts, te, cs, ce, cr)
        real(real64), intent(in) :: ts, te
        integer(int32), intent(in) :: cs, ce, cr
        write(stdout, *) 'Elapsed time per iteration (cpu_time):', te - ts
        write(stdout, *) 'Elapsed time per iteration (system_clock):', &
        real(ce - cs, real64)/real(cr, real64)
    end subroutine set_elapsed
    
    subroutine set_duration(nt, dt, ts, te, n)
        integer(int32), intent(in) :: nt, n
        real(real64), intent(in) :: dt, ts, te
        integer(int32) :: h, m, s
        real(real64) :: t
        t = (te - ts)*(nt - n)
        h = int(t)/3600
        m = mod(int(t)/60, 60)
        s = mod(int(t), 60)
        write(stdout, *) 'Iteration', n, 'of', nt
        write(stdout, *) 'Timestep:', dt
        write(stdout, *) 'Estimated time remaining (h:m:s):', h, m, s
        write(stdout, *)
    end subroutine set_duration

    subroutine set_data(Uc, pc, msh, n)
        real(real64), dimension(:, :), intent(in) :: Uc
        real(real64), dimension(:), intent(in) :: pc
        type(Mesh), intent(in) :: msh
        integer(int32), intent(in) :: n
        character(50) :: file
        integer(int32) :: i
        write(file, '(a, i10.10, a)') 'Wc_', n, '.txt'
        open(10, file='data/'//trim(file))
        do i = 1, msh % n_cells
            write(10, *) Uc(i, :), pc(i)
        end do
        close(10)
    end subroutine set_data

end module mod_utils
