module mod_solve
    use iso_fortran_env, only: int32, real64
    use mod_mesh, only: Mesh
    use mod_config, only: Configuration
    implicit none

    private
    public :: &
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
    
contains

    subroutine set_initial_conditions(U0, Uc, msh, config)
        real(real64), dimension(2), intent(in out) :: U0
        real(real64), dimension(:, :), intent(in out) :: Uc
        type(Mesh), intent(in) :: msh
        type(Configuration), intent(in) :: config
        U0 = [1., 0.]
        Uc(:, 1) = U0(1)
        Uc(:, 2) = U0(2)
    end subroutine set_initial_conditions

    subroutine set_timestep(config, Uc, msh)
        type(Configuration), intent(in out) :: config
        real(real64), dimension(:, :), intent(in) :: Uc
        type(Mesh), intent(in) :: msh
        integer(int32) :: i
        real(real64) :: dt, u, x
        if (config % CFL /= 0.) then
            config % dt = 1.
            do concurrent(i = 1:msh % n_cells)
                u = sqrt(Uc(i, 1)**2 + Uc(i, 2)**2)
                x = sqrt(msh % cells(i) % cell_volume)
                dt = config % CFL*x/u
                if (config % dt > dt) config % dt = dt
            end do
        end if
    end subroutine
    
    subroutine set_node_velocities(Un, Uc, U0, msh)
        real(real64), dimension(:, :), intent(in out) :: Un
        real(real64), dimension(:, :), intent(in) :: Uc
        real(real64), dimension(2), intent(in) :: U0
        type(Mesh), intent(in) :: msh
        integer(int32) :: i, j, cj
        real(real64) :: wnj
        real(real64), dimension(2) :: Unj
        do concurrent(i = 1:msh % n_nodes)
            if (trim(msh % nodes(i) % node_type) == 'INTERIOR' .or. &
            trim(msh % nodes(i) % node_type) == 'OUTLET') then
                Unj = 0.
                do concurrent(j = 1:msh % nodes(i) % n_node_cells)
                    cj = msh % nodes(i) % node_cells(j)
                    wnj = msh % nodes(i) % node_cell_weights(j)
                    Unj = Unj + wnj*Uc(cj, :)
                end do
                Un(i, :) = Unj
            else if (trim(msh % nodes(i) % node_type) == 'INLET') then
                Un(i, :) = U0
            else if (trim(msh % nodes(i) % node_type) == 'WALL') then
                Un(i, :) = 0.
            else
                error stop 'Error: invalid node type'
            end if
        end do
    end subroutine set_node_velocities

    subroutine set_face_velocities(Uf, Uc, U0, msh)
        real(real64), dimension(:, :), intent(in out) :: Uf
        real(real64), dimension(:, :), intent(in) :: Uc
        real(real64), dimension(2), intent(in) :: U0
        type(Mesh), intent(in) :: msh
        integer(int32) :: i, c1, c2
        real(real64) :: wf
        do concurrent(i = 1:msh % n_faces)
            if (trim(msh % faces(i) % face_type) == 'INTERIOR' .or. &
            trim(msh % faces(i) % face_type) == 'OUTLET') then
                c1 = msh % faces(i) % face_cells(1)
                c2 = msh % faces(i) % face_cells(2)
                wf = msh % faces(i) % face_cell_weight
                Uf(i, :) = wf*Uc(c1, :) + (1 - wf)*Uc(c2, :)
            else if (trim(msh % faces(i) % face_type) == 'INLET') then
                Uf(i, :) = U0
            else if (trim(msh % faces(i) % face_type) == 'WALL') then
                Uf(i, :) = 0.
            else
                error stop 'Error: invalid face type'
            end if
        end do
    end subroutine set_face_velocities

    subroutine set_velocity_gradients(Du, Dv, Uf, msh)
        real(real64), dimension(:, :), intent(in out) :: Du, Dv
        real(real64), dimension(:, :), intent(in) :: Uf
        type(Mesh), intent(in) :: msh
        integer(int32) :: i, j, fj
        real(real64) :: Vc, Sf, sign
        real(real64), dimension(2) :: Duf, Dvf, nf
        do concurrent(i = 1:msh % n_cells)
            Vc = msh % cells(i) % cell_volume
            Duf = 0.
            Dvf = 0.
            do concurrent(j = 1:3)
                fj = msh % cells(i) % cell_faces(j)
                Sf = msh % faces(fj) % face_area
                nf = msh % faces(fj) % face_normal
                sign = msh % cells(i) % cell_face_signs(j)
                Duf = Duf + sign*Uf(fj, 1)*Sf*nf
                Dvf = Dvf + sign*Uf(fj, 2)*Sf*nf
            end do
            Du(i, :) = Duf/Vc
            Dv(i, :) = Dvf/Vc
        end do
    end subroutine set_velocity_gradients

    subroutine set_diffusion_fluxes(Fd, Un, Uf, Uc, msh, config)
        real(real64), dimension(:, :), intent(in out) :: Fd
        real(real64), dimension(:, :), intent(in) :: Un, Uf, Uc
        type(Mesh), intent(in) :: msh
        type(Configuration), intent(in) :: config
        integer(int32) :: i, c1, c2, na, nb
        real(real64) :: Sf, df, ef
        do concurrent(i = 1:msh % n_faces)
            c1 = msh % faces(i) % face_cells(1)
            c2 = msh % faces(i) % face_cells(2)
            na = msh % faces(i) % face_nodes(1)
            nb = msh % faces(i) % face_nodes(2)
            Sf = msh % faces(i) % face_area
            df = msh % faces(i) % face_cell_distances(1)
            ef = msh % faces(i) % face_cell_distances(2)
            if (trim(msh % faces(i) % face_type) == 'INTERIOR') then
                Fd(i, :) = (Uc(c2, :) - Uc(c1, :))/df - (Un(nb, :) - Un(na, :))/Sf/df*ef
            else ! boundary faces
                Fd(i, :) = (Uf(i, :) - Uc(c1, :))/df - (Un(nb, :) - Un(na, :))/Sf/df*ef
            end if
            Fd(i, :) = Fd(i, :)/config % Re
        end do
    end subroutine set_diffusion_fluxes

    subroutine set_convection_fluxes(Fc, Un, Uf, Uc, Du, Dv, msh, config)
        real(real64), dimension(:, :), intent(in out) :: Fc
        real(real64), dimension(:, :), intent(in) :: Un, Uf, Uc, Du, Dv
        type(Mesh), intent(in) :: msh
        type(Configuration), intent(in) :: config
        integer(int32) :: i, c1, c2
        real(real64) :: mf, Pe
        real(real64), dimension(2) :: nf, Lf, Uuu, Uff
        do concurrent(i = 1:msh % n_faces)
            nf = msh % faces(i) % face_normal
            mf = dot_product(Uf(i, :), nf)
            if (trim(msh % faces(i) % face_type) == 'INTERIOR') then
                Lf = msh % faces(i) % face_cell_vector
                Pe = abs(dot_product(Uf(i, :), Lf))*config % Re ! local peclet number
                if (Pe < config % Pe) then ! centred difference
                    Fc(i, :) = Uf(i, :)*mf
                else ! upwind tvd scheme
                    c1 = msh % faces(i) % face_cells(1)
                    c2 = msh % faces(i) % face_cells(2)
                    if (mf >= 0.) then
                        Uuu = Uc(c2, :) - 2*[dot_product(Du(c1, :), Lf), dot_product(Dv(c1, :), Lf)]
                        Uff = get_convected_quantity(Uc(c1, :), Uc(c2, :), Uuu)
                    else
                        Uuu = Uc(c1, :) + 2*[dot_product(Du(c2, :), Lf), dot_product(Dv(c2, :), Lf)]
                        Uff = get_convected_quantity(Uc(c2, :), Uc(c1, :), Uuu)
                    end if
                    Fc(i, :) = Uff*mf
                end if
            else ! boundary faces
                Fc(i, :) = Uf(i, :)*mf
            end if            
        end do
    end subroutine set_convection_fluxes

    pure elemental function get_convected_quantity(phi_u, phi_d, phi_uu) result(res)
        real(real64), intent(in) :: phi_u, phi_d, phi_uu
        real(real64) :: res
        real(real64) :: rf
        if (abs(phi_u - phi_d) > 1e-15) then ! phi_u /= phi_d
            rf = (phi_u - phi_uu)/(phi_d - phi_u)
            res = phi_u + .5*(phi_d - phi_u)*get_limiter(rf)
        else ! phi_u == phi_d
            res = phi_u
        end if
    end function get_convected_quantity

    pure elemental function get_limiter(rf) result(res)
        real(real64), intent(in) :: rf
        real(real64) :: res
        integer(int32) :: i
        res = max(0., min(1., 2.*rf), min(2., rf))
    end function get_limiter

    subroutine set_intermediate_velocities(Ut, Uc, Fd, Fc, msh, config)
        real(real64), dimension(:, :), intent(in out) :: Ut
        real(real64), dimension(:, :), intent(in) :: Uc, Fd, Fc
        type(Mesh), intent(in) :: msh
        type(Configuration), intent(in) :: config
        integer(int32) :: i, j, fj
        real(real64) :: Vc, Sf, sign
        real(real64), dimension(2) :: Ff
        do concurrent(i = 1:msh % n_cells)
            Vc = msh % cells(i) % cell_volume
            Ff = 0.
            do concurrent(j = 1:3)
                fj = msh % cells(i) % cell_faces(j)
                Sf = msh % faces(fj) % face_area
                sign = msh % cells(i) % cell_face_signs(j)
                Ff = Ff + sign*(Fd(fj, :) - Fc(fj, :))*Sf
            end do
            Ut(i, :) = Uc(i, :) + config % dt/Vc*Ff
        end do
    end subroutine set_intermediate_velocities

    subroutine set_poisson_coefficients(ac, af, msh)
        real(real64), dimension(:), intent(in out) :: ac
        real(real64), dimension(:, :), intent(in out) :: af
        type(Mesh), intent(in) :: msh
        integer(int32) :: i, j, fj
        real(real64) :: Sf, df
        do concurrent(i = 1:msh % n_cells)
            ac(i) = 0.
            do concurrent(j = 1:3)
                fj = msh % cells(i) % cell_faces(j)
                Sf = msh % faces(fj) % face_area
                df = msh % faces(fj) % face_cell_distances(1)
                if (trim(msh % faces(fj) % face_type) == 'INTERIOR') then
                    af(i, j) = Sf/df
                    ac(i) = ac(i) - Sf/df
                else if (trim(msh % faces(fj) % face_type) == 'WALL') then
                    af(i, j) = 0.
                    ac(i) = ac(i) + 0.
                else if (trim(msh % faces(fj) % face_type) == 'INLET') then
                    af(i, j) = 0.
                    ac(i) = ac(i) + 0.
                else if (trim(msh % faces(fj) % face_type) == 'OUTLET') then
                    af(i, j) = 0.
                    ac(i) = ac(i) - 2.*Sf/df ! zero pressure outlet (pcf = -pc)
                else
                    error stop 'Error: invalid face type'
                end if
            end do
        end do
    end subroutine set_poisson_coefficients
    
    subroutine set_poisson_constants(bc, Utf, msh, config)
        real(real64), dimension(:), intent(in out) :: bc
        real(real64), dimension(:, :), intent(in) :: Utf
        type(Mesh), intent(in) :: msh
        type(Configuration), intent(in) :: config
        integer(int32) :: i, j, fj
        real(real64) :: Sf, sign
        real(real64), dimension(2) :: nf
        do concurrent(i = 1:msh % n_cells)
            bc(i) = 0.
            do concurrent(j = 1:3)
                fj = msh % cells(i) % cell_faces(j)
                Sf = msh % faces(fj) % face_area
                nf = msh % faces(fj) % face_normal
                sign = msh % cells(i) % cell_face_signs(j)
                bc(i) = bc(i) + sign*dot_product(Utf(fj, :), nf)*Sf/config % dt
            end do
        end do
    end subroutine set_poisson_constants
    
    subroutine set_cell_pressures(pc, ac, af, bc, msh, config)
        real(real64), dimension(:), intent(in out) :: pc
        real(real64), dimension(:), intent(in) :: ac, bc
        real(real64), dimension(:, :), intent(in) :: af
        type(Mesh), intent(in) :: msh
        type(Configuration), intent(in) :: config
        integer(int32) :: i, j, k, cj
        real(real64) :: acf
        do concurrent(k = 1:config % nk) ! gauss-seidel
            do concurrent(i = 1:msh % n_cells)
                acf = 0.
                do concurrent(j = 1:3)
                    cj = msh % cells(i) % cell_cells(j)
                    acf = acf + af(i, j)*pc(cj)
                end do
                pc(i) = (bc(i) - acf)/ac(i)
            end do
        end do
    end subroutine set_cell_pressures

    subroutine set_face_pressures(pf, pc, msh)
        real(real64), dimension(:), intent(in out) :: pf
        real(real64), dimension(:), intent(in) :: pc
        type(Mesh), intent(in) :: msh
        integer(int32) :: i, c1, c2
        real(real64) :: wf
        do concurrent(i = 1:msh % n_faces)
            if (trim(msh % faces(i) % face_type) == 'OUTLET') then
                pf(i) = 0. ! zero pressure outlet (pf = 0)
            else
                c1 = msh % faces(i) % face_cells(1)
                c2 = msh % faces(i) % face_cells(2)
                wf = msh % faces(i) % face_cell_weight
                pf(i) = wf*pc(c1) + (1 - wf)*pc(c2)
            end if
        end do
    end subroutine set_face_pressures

    subroutine set_corrected_velocities(Uc, Ut, pf, msh, config)
        real(real64), dimension(:, :), intent(in out) :: Uc
        real(real64), dimension(:, :), intent(in) :: Ut
        real(real64), dimension(:), intent(in) :: pf
        type(Mesh), intent(in) :: msh
        type(Configuration), intent(in) :: config
        integer(int32) :: i, j, fj, c1, c2
        real(real64) :: Vc, Sf, sign
        real(real64), dimension(2) :: nf, Ff
        do concurrent(i = 1:msh % n_cells)
            Vc = msh % cells(i) % cell_volume
            Ff = 0.
            do concurrent(j = 1:3)
                fj = msh % cells(i) % cell_faces(j)
                Sf = msh % faces(fj) % face_area
                nf = msh % faces(fj) % face_normal
                sign = msh % cells(i) % cell_face_signs(j)
                Ff = Ff + sign*pf(fj)*Sf*nf
            end do
            Uc(i, :) = Ut(i, :) - config % dt/Vc*Ff
        end do
    end subroutine set_corrected_velocities

end module mod_solve
