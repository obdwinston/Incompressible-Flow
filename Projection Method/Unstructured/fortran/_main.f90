program main
    
    implicit none

    character(1000), parameter :: mesh = 'cylinder.su2'
    character(1000), parameter :: folder = 'data/'
    character(1000), parameter :: variable = 'speed'
    integer, parameter :: interval = 50

    real(8), parameter :: Ui = 5.
    real(8), parameter :: Vi = 0.
    real(8), parameter :: nu = .01
    real(8), parameter :: t = 5.
    real(8), parameter :: dt = 1e-4

    ! program start

    integer :: file_unit
    character(1000), dimension(10) :: line
    character(1000) :: file_name

    real(8) :: n1x, n1y, n2x, n2y, n3x, n3y, d1, d2
    integer :: i, j, k, n, nt, ni, n1, n2, n3, fi, fj, ci, cj, c1, c2
    integer, allocatable, dimension(:) :: l1, l2, l3

    real(8), allocatable, dimension(:, :) :: nodes
    integer, allocatable, dimension(:, :) :: faces, cells
    integer :: n_nodes, n_faces, n_cells

    integer, allocatable, dimension(:, :) :: inlet_faces, outlet_faces, wall_faces, body_faces
    integer, allocatable, dimension(:, :) :: boundary_faces, interior_faces
    integer, allocatable, dimension(:) :: i_inlet_nodes, i_outlet_nodes, i_wall_nodes, i_body_nodes
    integer, allocatable, dimension(:) :: i_inlet_faces, i_outlet_faces, i_wall_faces, i_body_faces
    integer, allocatable, dimension(:) :: i_boundary_faces, i_interior_faces
    integer :: n_inlet_nodes, n_outlet_nodes, n_wall_nodes, n_body_nodes
    integer :: n_inlet_faces, n_outlet_faces, n_wall_faces, n_body_faces
    integer :: n_boundary_faces, n_interior_faces

    integer, allocatable, dimension(:, :) :: face_cells, cell_faces, cell_cells, cell_signs
    real(8), allocatable, dimension(:) :: Vc, Sf, df, ef, wf
    real(8), allocatable, dimension(:, :) :: cc, nf, tf, cf, lf

    integer, allocatable, dimension(:, :) :: node_cells
    integer, allocatable, dimension(:) :: n_node_cells
    real(8), allocatable, dimension(:, :) :: dn, wn

    real(8) :: ufi, vfi, mfi, usfj, vsfj, msfj, val
    real(8), allocatable, dimension(:) :: u, us, un, v, vs, vn, p, pf, Jcu, Jcv, Jdu, Jdv, aC, bC
    real(8), allocatable, dimension(:, :) :: aF

    open(newunit=file_unit, file=trim(mesh), action='read')
    
    ! cells [n1, n2, n3]

    read(file_unit, *)
    read(file_unit, *) line(1), n_cells
    allocate(cells(n_cells, 3))
    do i = 1, n_cells
        read(file_unit, *) line(1), cells(i, 1), cells(i, 2), cells(i, 3)
    end do
    cells = cells + 1

    ! nodes [p1, p2]

    read(file_unit, *) line(1), n_nodes
    allocate(nodes(n_nodes, 2))
    do i = 1, n_nodes
        read(file_unit, *) nodes(i, 1), nodes(i, 2)
    end do

    ! faces [n1, n2]

    read(file_unit, *)
    read(file_unit, *)
    read(file_unit, *) line(1), n_inlet_faces
    allocate(inlet_faces(n_inlet_faces, 2))
    do i = 1, n_inlet_faces
        read(file_unit, *) line(1), inlet_faces(i, 1), inlet_faces(i, 2)
    end do
    inlet_faces = inlet_faces + 1

    read(file_unit, *)
    read(file_unit, *) line(1), n_outlet_faces
    allocate(outlet_faces(n_outlet_faces, 2))
    do i = 1, n_outlet_faces
        read(file_unit, *) line(1), outlet_faces(i, 1), outlet_faces(i, 2)
    end do
    outlet_faces = outlet_faces + 1

    read(file_unit, *)
    read(file_unit, *) line(1), n_wall_faces
    allocate(wall_faces(n_wall_faces, 2))
    do i = 1, n_wall_faces
        read(file_unit, *) line(1), wall_faces(i, 1), wall_faces(i, 2)
    end do
    wall_faces = wall_faces + 1

    read(file_unit, *)
    read(file_unit, *) line(1), n_body_faces
    allocate(body_faces(n_body_faces, 2))
    do i = 1, n_body_faces
        read(file_unit, *) line(1), body_faces(i, 2), body_faces(i, 1)  ! swap
    end do
    body_faces = body_faces + 1

    l1 = [inlet_faces(:, 1), outlet_faces(:, 1), wall_faces(:, 1), body_faces(:, 1)]
    l2 = [inlet_faces(:, 2), outlet_faces(:, 2), wall_faces(:, 2), body_faces(:, 2)]
    do i = 1, n_cells
        n1 = cells(i, 1)
        n2 = cells(i, 2)
        n3 = cells(i, 3)
        call append_array(l1, l2, n1, n2)
        call append_array(l1, l2, n2, n3)
        call append_array(l1, l2, n3, n1)
    end do
    n_faces = size(l1)
    allocate(faces(n_faces, 2))
    faces(:, 1) = l1
    faces(:, 2) = l2

    n_boundary_faces = n_wall_faces + n_body_faces + n_inlet_faces + n_outlet_faces
    allocate(boundary_faces(n_boundary_faces, 2))
    boundary_faces(:, 1) = faces(:n_boundary_faces, 1)
    boundary_faces(:, 2) = faces(:n_boundary_faces, 2)

    n_interior_faces = n_faces - n_boundary_faces
    allocate(interior_faces(n_interior_faces, 2))
    interior_faces(:, 1) = faces(n_boundary_faces + 1:, 1)
    interior_faces(:, 2) = faces(n_boundary_faces + 1:, 2)

    fi = 1
    fj = n_inlet_faces
    i_inlet_faces = [(i, i = fi, fj)]

    fi = fi + n_inlet_faces
    fj = fj + n_outlet_faces
    i_outlet_faces = [(i, i = fi, fj)]

    fi = fi + n_outlet_faces
    fj = fj + n_wall_faces
    i_wall_faces = [(i, i = fi, fj)]

    fi = fi + n_wall_faces
    fj = fj + n_body_faces
    i_body_faces = [(i, i = fi, fj)]

    fi = 1
    fj = n_boundary_faces
    i_boundary_faces = [(i, i = fi, fj)]

    fi = n_boundary_faces + 1
    fj = n_faces
    i_interior_faces = [(i, i = fi, fj)]

    i_inlet_nodes = [integer ::]
    do i = 1, n_inlet_faces
        n1 = inlet_faces(i, 1)
        n2 = inlet_faces(i, 2)
        if (.not. any(i_inlet_nodes == n1)) then
            i_inlet_nodes = [i_inlet_nodes, n1]
        end if
        if (.not. any(i_inlet_nodes == n2)) then
            i_inlet_nodes = [i_inlet_nodes, n2]
        end if
    end do
    n_inlet_nodes = size(i_inlet_nodes)

    i_outlet_nodes = [integer ::]
    do i = 1, n_outlet_faces
        n1 = outlet_faces(i, 1)
        n2 = outlet_faces(i, 2)
        if (.not. any(i_outlet_nodes == n1)) then
            i_outlet_nodes = [i_outlet_nodes, n1]
        end if
        if (.not. any(i_outlet_nodes == n2)) then
            i_outlet_nodes = [i_outlet_nodes, n2]
        end if
    end do
    n_outlet_nodes = size(i_outlet_nodes)

    i_wall_nodes = [integer ::]
    do i = 1, n_wall_faces
        n1 = wall_faces(i, 1)
        n2 = wall_faces(i, 2)
        if (.not. any(i_wall_nodes == n1)) then
            i_wall_nodes = [i_wall_nodes, n1]
        end if
        if (.not. any(i_wall_nodes == n2)) then
            i_wall_nodes = [i_wall_nodes, n2]
        end if
    end do
    n_wall_nodes = size(i_wall_nodes)

    i_body_nodes = [integer ::]
    do i = 1, n_body_faces
        n1 = body_faces(i, 1)
        n2 = body_faces(i, 2)
        if (.not. any(i_body_nodes == n1)) then
            i_body_nodes = [i_body_nodes, n1]
        end if
        if (.not. any(i_body_nodes == n2)) then
            i_body_nodes = [i_body_nodes, n2]
        end if
    end do
    n_body_nodes = size(i_body_nodes)

    ! cell faces

    allocate(cell_faces(n_cells, 3))
    do i = 1, n_cells
        n1 = cells(i, 1)
        n2 = cells(i, 2)
        n3 = cells(i, 3)
        do j = 1, n_faces
            if (compare_array(l1(j), l2(j), n1, n2)) then
                cell_faces(i, 1) = j
            end if
            if (compare_array(l1(j), l2(j), n2, n3)) then
                cell_faces(i, 2) = j
            end if
            if (compare_array(l1(j), l2(j), n3, n1)) then
                cell_faces(i, 3) = j
            end if
        end do
    end do

    ! face cells

    allocate(face_cells(n_faces, 2))
    do i = 1, n_faces
        l3 = [integer ::]
        do j = 1, n_cells
            if (any(cell_faces(j, :) == i)) then
                l3 = [l3, j]
            end if
        end do        
        if (size(l3) == 1) then
            face_cells(i, 1) = l3(1)
            face_cells(i, 2) = l3(1)
        else
            face_cells(i, 1) = l3(1)
            face_cells(i, 2) = l3(2)
        end if
    end do

    ! Vc, cc

    allocate(Vc(n_cells))
    allocate(cc(n_cells, 2))
    allocate(cell_cells(n_cells, 3))
    allocate(cell_signs(n_cells, 3))
    do i = 1, n_cells
        n1 = cells(i, 1)
        n2 = cells(i, 2)
        n3 = cells(i, 3)
        n1x = nodes(n1, 1)
        n1y = nodes(n1, 2)
        n2x = nodes(n2, 1)
        n2y = nodes(n2, 2)
        n3x = nodes(n3, 1)
        n3y = nodes(n3, 2)

        Vc(i) = .5*abs(n1x*(n2y - n3y) + n2x*(n3y - n1y) + n3x*(n1y - n2y))
        cc(i, :) = [(n1x + n2x + n3x)/3, (n1y + n2y + n3y)/3]

        do j = 1, 3
            fj = cell_faces(i, j)

            if (i == face_cells(fj, 1)) then
                cell_cells(i, j) = face_cells(fj, 2)
                cell_signs(i, j) = 1
            else
                cell_cells(i, j) = face_cells(fj, 1)
                cell_signs(i, j) = -1
            end if
        end do
    end do

    ! Sf, nf, tf, cf, df, ef, lf, wf

    allocate(Sf(n_faces))
    allocate(df(n_faces))
    allocate(ef(n_faces))
    allocate(wf(n_faces))
    allocate(nf(n_faces, 2))
    allocate(tf(n_faces, 2))
    allocate(cf(n_faces, 2))
    allocate(lf(n_faces, 2))
    do i = 1, n_faces
        n1 = faces(i, 1)
        n2 = faces(i, 2)
        n1x = nodes(n1, 1)
        n1y = nodes(n1, 2)
        n2x = nodes(n2, 1)
        n2y = nodes(n2, 2)

        Sf(i) = sqrt((n2x - n1x)**2 + (n2y - n1y)**2)
        nf(i, :) = [(n2y - n1y)/Sf(i), -(n2x - n1x)/Sf(i)]
        tf(i, :) = [(n2x - n1x)/Sf(i), (n2y - n1y)/Sf(i)]
        cf(i, :) = [.5*(n1x + n2x), .5*(n1y + n2y)]

        c1 = face_cells(i, 1)
        c2 = face_cells(i, 2)
        
        if (c1 == c2) then  ! boundary face
            lf(i, :) = cf(i, :) - cc(c1, :)
            df(i) = abs(dot_product(lf(i, :), nf(i, :)))
            ef(i) = abs(dot_product(lf(i, :), tf(i, :)))
        else  ! interior face
            lf(i, :) = cc(c2, :) - cc(c1, :)
            df(i) = abs(dot_product(lf(i, :), nf(i, :)))
            ef(i) = abs(dot_product(lf(i, :), tf(i, :)))
        end if

        d1 = 1/norm2(cc(c1, :) - cf(i, :))  ! reciprocal
        d2 = 1/norm2(cc(c2, :) - cf(i, :))  ! reciprocal
        wf(i) = d1/(d1 + d2)
    end do

    ! wn

    allocate(wn(n_nodes, 10))
    allocate(dn(n_nodes, 10))
    allocate(node_cells(n_nodes, 10))
    allocate(n_node_cells(n_nodes))
    do i = 1, n_nodes
        k = 1
        do j = 1, n_cells
            if (any(cells(j, :) == i)) then
                dn(i, k) = 1/norm2(cc(j, :) - nodes(i, :))  ! reciprocal
                node_cells(i, k) = j
                k = k + 1
            end if
        end do
        n_node_cells(i) = k - 1

        do j = 1, n_node_cells(i)
            wn(i, j) = dn(i, j)/sum(dn(i, :))
        end do
    end do

    ! save mesh

    open(100, file=trim(folder)//'cc.txt')
    do i = 1, n_cells
        write(100, *) cc(i, 1), cc(i, 2)
    end do

    ! allocate variables

    allocate(u(n_cells))
    allocate(us(n_cells))
    allocate(un(n_nodes))

    allocate(v(n_cells))
    allocate(vs(n_cells))
    allocate(vn(n_nodes))

    allocate(p(n_cells))
    allocate(pf(n_faces))

    allocate(Jcu(n_faces))
    allocate(Jcv(n_faces))
    allocate(Jdu(n_faces))
    allocate(Jdv(n_faces))

    allocate(aC(n_cells))
    allocate(aF(n_cells, 3))
    allocate(bC(n_cells))

    ! pressure coefficients

    do i = 1, n_cells
        aC(i) = 0.
        do j = 1, 3
            fj = cell_faces(i, j)

            if (any(fj == i_inlet_faces)) then
                aF(i, j) = 0.
                aC(i) = aC(i) + 0.
            else if (any(fj == i_outlet_faces)) then
                aF(i, j) = 0.  ! zero pressure outlet
                aC(i) = aC(i) - Sf(fj)/df(fj)
            else if (any(fj == i_wall_faces) .or. any(fj == i_body_faces)) then
                aF(i, j) = 0.
                aC(i) = aC(i) + 0.
            else  ! interior face
                aF(i, j) = Sf(fj)/df(fj)
                aC(i) = aC(i) - Sf(fj)/df(fj)
            end if
        end do
    end do

    ! time loop
    
    nt = int(t/dt)
    do n = 1, nt
        
        ! projection step

        do i = 1, n_nodes
            un(i) = 0.
            vn(i) = 0.
            do j = 1, n_node_cells(i)
                cj = node_cells(i, j)

                un(i) = un(i) + wn(i, j)*u(cj)
                vn(i) = vn(i) + wn(i, j)*v(cj)
            end do
        end do
        do i = 1, n_inlet_nodes
            ni = i_inlet_nodes(i)
            un(ni) = Ui
            vn(ni) = Vi
        end do
        do i = 1, n_wall_nodes
            ni = i_wall_nodes(i)
            un(ni) = 0.
            vn(ni) = 0.
        end do
        do i = 1, n_body_nodes
            ni = i_body_nodes(i)
            un(ni) = 0.
            vn(ni) = 0.
        end do
        
        do i = 1, n_interior_faces
            fi = i_interior_faces(i)
            c1 = face_cells(fi, 1)
            c2 = face_cells(fi, 2)
            n1 = faces(fi, 1)
            n2 = faces(fi, 2)

            Jdu(fi) = nu*((u(c2) - u(c1))*Sf(fi) - (un(n2) - un(n1))*ef(fi))/df(fi)
            Jdv(fi) = nu*((v(c2) - v(c1))*Sf(fi) - (vn(n2) - vn(n1))*ef(fi))/df(fi)

            ufi = wf(fi)*u(c1) + (1 - wf(fi))*u(c2)
            vfi = wf(fi)*v(c1) + (1 - wf(fi))*v(c2)
            mfi = dot_product([ufi, vfi], nf(fi, :))

            Jcu(fi) = mfi*Sf(fi)*ufi
            Jcv(fi) = mfi*Sf(fi)*vfi
        end do

        do i = 1, n_boundary_faces
            fi = i_boundary_faces(i)
            ci = face_cells(fi, 1)

            if (any(fi == i_inlet_faces)) then
                Jdu(fi) = nu*(Ui - u(ci))*Sf(fi)/df(fi)
                Jdv(fi) = nu*(Vi - v(ci))*Sf(fi)/df(fi)

                mfi = dot_product([Ui, Vi], nf(fi, :))
                Jcu(fi) = mfi*Sf(fi)*Ui
                Jcv(fi) = mfi*Sf(fi)*Vi
            end if

            if (any(fi == i_outlet_faces)) then
                Jdu(fi) = 0.
                Jdv(fi) = 0.

                mfi = dot_product([u(ci), v(ci)], nf(fi, :))
                Jcu(fi) = mfi*Sf(fi)*u(ci)
                Jcv(fi) = mfi*Sf(fi)*v(ci)
            end if

            if (any(fi == i_wall_faces) .or. any(fi == i_body_faces)) then
                Jdu(fi) = -nu*u(ci)*Sf(fi)/df(fi)
                Jdv(fi) = -nu*v(ci)*Sf(fi)/df(fi)

                Jcu(fi) = 0.
                Jcv(fi) = 0.
            end if
        end do

        do i = 1, n_cells
            us(i) = u(i)
            vs(i) = v(i)
            do j = 1, 3
                fj = cell_faces(i, j)

                us(i) = us(i) + cell_signs(i, j)*dt/Vc(i)*(Jdu(fj) - Jcu(fj))
                vs(i) = vs(i) + cell_signs(i, j)*dt/Vc(i)*(Jdv(fj) - Jcv(fj))
            end do
        end do

        ! correction step

        do i = 1, n_cells
            bC(i) = 0.
            do j = 1, 3
                fj = cell_faces(i, j)
                c1 = face_cells(fj, 1)
                c2 = face_cells(fj, 2)

                if (any(fj == i_wall_faces) .or. any(fj == i_body_faces)) then
                    bC(i) = bC(i) + 0.
                else
                    usfj = wf(fj)*us(c1) + (1 - wf(fj))*us(c2)
                    vsfj = wf(fj)*vs(c1) + (1 - wf(fj))*vs(c2)
                    msfj = dot_product([usfj, vsfj], nf(fj, :))

                    bC(i) = bC(i) + cell_signs(i, j)*msfj*Sf(fj)/dt
                end if
            end do
        end do

        do k = 1, 5  ! Gauss-Seidel
            do i = 1, n_cells
                val = 0.
                do j = 1, 3
                    cj = cell_cells(i, j)

                    val = val + aF(i, j)*p(cj)
                p(i) = (bC(i) - val)/aC(i)
                end do
            end do
        end do

        ! update solution

        do i = 1, n_faces
            c1 = face_cells(i, 1)
            c2 = face_cells(i, 2)
            pf(i) = wf(i)*p(c1) + (1 - wf(i))*p(c2)
        end do
        do i = 1, n_outlet_faces
            fi = i_outlet_faces(i)
            pf(fi) = 0.  ! zero pressure outlet
        end do

        do i = 1, n_cells
            u(i) = us(i)
            v(i) = vs(i)
            do j = 1, 3
                fj = cell_faces(i, j)

                u(i) = u(i) - dt/Vc(i)*pf(fj)*cell_signs(i, j)*nf(fj, 1)*Sf(fj)
                v(i) = v(i) - dt/Vc(i)*pf(fj)*cell_signs(i, j)*nf(fj, 2)*Sf(fj)
            end do
        end do

        write(*, '(2(i10, 5x), f10.5)') n, nt, real(n)/real(nt)*100.
        print *, maxval(u), minval(u)

        ! save data

        if (mod(n, interval) == 0) then
            write(file_name, '(a, i10.10, a)') trim(variable), n, '.txt'
            open(100, file=trim(folder)//trim(file_name))

            do i = 1, n_cells
                write(100, *) norm2([u(i), v(i)])
            end do
            print *, 'save'
        end if

    end do

contains

    subroutine append_array(a1, a2, b1, b2)
        
        integer, allocatable, intent(in out) :: a1(:), a2(:)
        integer, intent(in) :: b1, b2
        
        integer, allocatable :: e1(:), e2(:)
        logical :: in

        allocate(e1, mold=a1)
        allocate(e2, mold=a2)
        e1(:) = b1
        e2(:) = b2

        in = any((e1 == a1) .and. (e2 == a2)) .or. any((e2 == a1) .and. (e1 == a2))
        if (.not. in) then
            a1 = [a1, b1]
            a2 = [a2, b2]
        end if

    end subroutine append_array

    pure function compare_array(a1, a2, b1, b2)
    
        integer, intent(in) :: a1, a2, b1, b2
        logical :: compare_array

        if ((b1 == a1 .and. b2 == a2) .or. (b2 == a1 .and. b1 == a2)) then
            compare_array = .true.
        else
            compare_array = .false.
        end if

    end function compare_array

end program main
