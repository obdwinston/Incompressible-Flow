program main
    use iso_fortran_env, only: int32, real64, stdout => output_unit
    use mod_mesh, only: Mesh
    use mod_config, only: Configuration
    use mod_solve, only: &
    set_initial_conditions, &
    set_timestep, &
    set_node_velocities, &
    set_face_velocities, &
    set_velocity_gradients, &
    set_diffusion_fluxes, &
    set_convection_fluxes, &
    set_intermediate_velocities, &
    set_poisson_coefficients, &
    set_poisson_constants, &
    set_cell_pressures, &
    set_face_pressures, &
    set_corrected_velocities
    use mod_utils, only: &
    set_field, &
    set_surface, &
    set_elapsed, &
    set_duration, &
    set_data
    implicit none

    type(Mesh) :: msh
    type(Configuration) :: config
    integer(int32) :: n, cs, ce, cr
    real(real64) :: ts, te
    real(real64), dimension(2) :: U0
    real(real64), dimension(:), allocatable :: pc, pf, ac, bc
    real(real64), dimension(:, :), allocatable :: Uc, Uf, Un, Ut, Utf, Fd, Fc, Du, Dv, af

    msh = Mesh()
    call msh % load_mesh('mesh/mesh.txt')
    call set_field(msh)
    call set_surface(msh)
    config = Configuration('config.txt')
    
    allocate(Uc(msh % n_cells, 2))
    allocate(Uf(msh % n_faces, 2))
    allocate(Un(msh % n_nodes, 2))
    allocate(Ut(msh % n_cells, 2))
    allocate(Utf(msh % n_faces, 2))
    allocate(Fd(msh % n_faces, 2))
    allocate(Fc(msh % n_faces, 2))
    allocate(Du(msh % n_cells, 2))
    allocate(Dv(msh % n_cells, 2))

    allocate(pc(msh % n_cells))
    allocate(pf(msh % n_faces))
    allocate(ac(msh % n_cells))
    allocate(af(msh % n_cells, 3))
    allocate(bc(msh % n_cells))

    call set_initial_conditions(U0, Uc, msh, config)

    write(stdout, *) 'Running simulation...'
    time_step: do n = 1, config % nt
        call cpu_time(ts)
        call system_clock(cs, cr)

        call set_timestep(config, Uc, msh)

        call set_node_velocities(Un, Uc, U0, msh)
        call set_face_velocities(Uf, Uc, U0, msh)
        call set_velocity_gradients(Du, Dv, Uf, msh)

        call set_diffusion_fluxes(Fd, Un, Uf, Uc, msh, config)
        call set_convection_fluxes(Fc, Un, Uf, Uc, Du, Dv, msh, config)
        call set_intermediate_velocities(Ut, Uc, Fd, Fc, msh, config)

        call set_face_velocities(Utf, Ut, U0, msh)
        call set_poisson_coefficients(ac, af, Utf, msh)
        call set_poisson_constants(bc, Utf, msh, config)

        call set_cell_pressures(pc, ac, af, bc, msh, config)
        call set_face_pressures(pf, pc, Utf, msh)
        
        call set_corrected_velocities(Uc, Ut, pf, msh, config)
        
        call cpu_time(te)
        call system_clock(ce)
        
        if (mod(n, config % nw) == 0) then
            write(stdout, *) 'umin=', minval(Uc(:, 1)), 'umax=', maxval(Uc(:, 1))
            write(stdout, *) 'vmin=', minval(Uc(:, 2)), 'vmax=', maxval(Uc(:, 2))
            call set_elapsed(ts, te, cs, ce, cr)
            call set_duration(config % nt, config % dt, ts, te, n)
            call set_data(Uc, pc, msh, n)
        end if
    end do time_step
    write(stdout, *) 'Done!'
end program main
