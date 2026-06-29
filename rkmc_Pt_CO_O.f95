! A fortran95 program for G95
! By HY
subroutine init_random_seed()
    integer(kind=4):: i, n, clock
    integer, dimension(:), allocatable::seed

    call random_seed(size=n)
    allocate (seed(n))

    call system_clock(count=clock)

    seed = clock + 37*(/(i - 1, i=1, n)/)
    call random_seed(put=seed)

    deallocate (seed)
end subroutine

program main
    implicit none

    character(len=120)::filename1, filename2, filename3, char_seq
    character(len=120)::filename4, filename5, filename6, filename7
    character(len=120)::filename8, filename9, filename10
    character(len=2):: elem

    integer(kind=4):: i, j, k

    integer(kind=1):: restart
    integer(kind=4):: natoms
    integer(kind=8):: nLoop
    real(kind=4):: clx, cly, clz, dimx, dimy, dimz
    integer(kind=2):: nlx, nly, nlz
    real(kind=4):: bx, by, bz
    real(kind=4):: x(100000), y(100000), z(100000)
    real(kind=4):: xcent, ycent, zcent, dxa, dya, dza
    real(kind=4):: x_sum, y_sum, z_sum
    real(kind=8):: tend

    integer(kind=4):: nbulk
    real(kind=4):: nsurf10, nsurf11, nsurf12, nsurf_CO, nsurf_O
    real(kind=4):: n100, n100_CO, n100_O
    real(kind=4):: n110, n110_CO, n110_O
    real(kind=4):: n111, n111_CO, n111_O
    real(kind=4):: nedge, nedge_CO, nedge_O
    real(kind=4):: xxx(100000), yyy(100000), zzz(100000)
    real(kind=4):: xxx_100(100000), yyy_100(100000), zzz_100(100000)
    real(kind=4):: xxx_110(100000), yyy_110(100000), zzz_110(100000)
    real(kind=4):: xxx_111(100000), yyy_111(100000), zzz_111(100000)
    real(kind=4):: xxx_edge(100000), yyy_edge(100000), zzz_edge(100000)
    real(kind=4):: x3cent, y3cent, z3cent, dx3a, dy3a, dz3a
    real(kind=4):: x3_sum, y3_sum, z3_sum

    integer(kind=4):: site_type(100000), cov_type(100000), cov_type_ini(100000)

    real(kind=4):: dx, dy, dz, rx, ry, rz, dr
    integer(kind=4):: cn(100000), cn_2(100000), effcn(100000), ijpick_cn(2, 12, 100000)
    integer(kind=4):: nnsite(12, 100000), nnnsite(54, 100000)
    real(kind=8):: gcn(100000)

    real(kind=8):: r_site, rpoint, rtot, rsite(100000), revent(8, 100000), &
                   rneis(12, 8, 100000), Ea_r(2, 12, 100000), ijpick_gcn(2, 12, 100000)

    integer(kind=4):: j_cn, kpick, ni_ads, nj, ni_effcn
    integer(kind=4):: i_cn, ii_step, ipick, jpick, nn_i, nn_j, nn_k, nf_1
    integer(kind=4):: natoms_ij(108), ni_site, nj_site
    integer(kind=8):: n_step_tot, i_step, COads_step, COdes_step, &
                      Oads_step, Odes_step, COdiff_step, Odiff_step, COOx_step
    real(kind=8):: Temperature, pCO, pO2
    real(kind=8):: ctime, dtime, randum, sum_rsite, sum_revent, sum_rneis
    real(kind=8):: gcn_ni

    integer(kind=4):: neis_i, neis_j
    integer(kind=4):: nn_i_effcn, jpick_effcn

    integer(kind=8):: COOx_100, COOx_110, COOx_111, COOx_edge

    open (1, file='input')
    read (1, *) nLoop
    read (1, *) clx, cly, clz 
    read (1, *) dimx, dimy, dimz
    read (1, *) elem
    read (1, *) Temperature
    read (1, *) pCO
    read (1, *) pO2
    close (1)

    nlx = int(dimx/clx)
    nly = int(dimy/cly)
    nlz = int(dimz/clz)
    bx = (nlx - 0.5)*clx
    by = (nly - 0.5)*cly
    bz = (nlz - 0.5)*clz


    ! bulk.xyz file and get mass center
    nbulk = 0
    do i = 1, nlx
        do j = 1, nly
            do k = 1, nlz
                nbulk = nbulk + 1
                xxx(nbulk) = i*clx
                yyy(nbulk) = j*cly
                zzz(nbulk) = k*clz
                nbulk = nbulk + 1
                xxx(nbulk) = (0.5 + i)*clx
                yyy(nbulk) = (0.5 + j)*cly
                zzz(nbulk) = k*clz
                nbulk = nbulk + 1
                xxx(nbulk) = i*clx
                yyy(nbulk) = (0.5 + j)*cly
                zzz(nbulk) = (0.5 + k)*clz
                nbulk = nbulk + 1
                xxx(nbulk) = (i + 0.5)*clx
                yyy(nbulk) = j*cly
                zzz(nbulk) = (0.5 + k)*clz
            end do
        end do
    end do
    x3_sum = 0
    y3_sum = 0
    z3_sum = 0
    do i = 1, nbulk
        x3_sum = x3_sum + xxx(i)
        y3_sum = y3_sum + yyy(i)
        z3_sum = z3_sum + zzz(i)
    end do
    x3cent = real(NINT(x3_sum/nbulk))
    y3cent = real(NINT(y3_sum/nbulk))
    z3cent = real(NINT(z3_sum/nbulk))
    ! center bulk sites
    do i = 1, nbulk
        dx3a = abs(xxx(i) - x3cent)
        dy3a = abs(yyy(i) - y3cent)
        dz3a = abs(zzz(i) - z3cent)
        if (dx3a .lt. 3.0 .and. dy3a .lt. 3.0 .and. dz3a .lt. 3.0) then
            x3cent = xxx(i)
            y3cent = yyy(i)
            z3cent = zzz(i)
            exit
        end if
    end do
    ! write a bulk file
    open (1, file='bulk.xyz')
    write (1, *) nbulk
    write (1, *)
    do i = 1, nbulk
        xxx(i) = xxx(i) - x3cent
        yyy(i) = yyy(i) - y3cent
        zzz(i) = zzz(i) - z3cent
        write (1, *) elem, xxx(i), yyy(i), zzz(i)
        site_type(i) = 3  ! fake site
        cov_type(i) = 5   ! do not adsorb
    end do
    close (1)

    ! read ini.xyz file and get mass center
    x_sum = 0
    y_sum = 0
    z_sum = 0
    open (1, file='ini.xyz')
    read (1, *) natoms
    read (1, *)
    do i = 1, natoms
        read (1, *) elem, x(i), y(i), z(i)
        x_sum = x_sum + x(i)
        y_sum = y_sum + y(i)
        z_sum = z_sum + z(i)
    end do
    close (1)
    xcent = real(NINT(x_sum/natoms))
    ycent = real(NINT(y_sum/natoms))
    zcent = real(NINT(z_sum/natoms))
    ! center atom site
    do i = 1, natoms
        dxa = abs(x(i) - xcent)
        dya = abs(y(i) - ycent)
        dza = abs(z(i) - zcent)
        if (dxa .lt. 3.0 .and. dya .lt. 3.0 .and. dza .lt. 3.0) then
            xcent = x(i)
            ycent = y(i)
            zcent = z(i)
            exit
        end if
    end do
    ! transfer ini position
    open (1, file='ini_new.xyz')
    write (1, *) natoms
    write (1, *)
    do i = 1, natoms
        x(i) = x(i) - xcent
        y(i) = y(i) - ycent
        z(i) = z(i) - zcent
        write (1, *) elem, x(i), y(i), z(i)
    end do
    close (1)

    ! bulk sites and covs classification
    do i = 1, natoms
        do j = 1, nbulk
            dx = abs(x(i) - xxx(j))
            dy = abs(y(i) - yyy(j))
            dz = abs(z(i) - zzz(j))
            if (dx .le. 0.5 .and. dy .le. 0.5 .and. dz .le. 0.5) then
                site_type(j) = 1
                exit
            end if
        end do
    end do

    open (1, file='ini_new.xyz')
    write (1, *) natoms
    write (1, *)
    do i = 1, natoms
        do j = 1, nbulk
            dx = abs(x(i) - xxx(j))
            dy = abs(y(i) - yyy(j))
            dz = abs(z(i) - zzz(j))
            if (dx .le. 0.5 .and. dy .le. 0.5 .and. dz .le. 0.5) then
                write (1, *) elem, xxx(j), yyy(j), zzz(j)
                exit
            end if
        end do
    end do
    close (1)


    ! cn (nnsite) and scn (nnnsite)
    cn = 0
    cn_2 = 0
    effcn = 0
    nnsite = 0
    nnnsite = 0
    do i = 1, nbulk
        do j = i + 1, nbulk
            rx = abs(xxx(i) - xxx(j))
            ry = abs(yyy(i) - yyy(j))
            rz = abs(zzz(i) - zzz(j))
            if (rx*2 .gt. bx) rx = (bx - rx) + 0.5*clx
            if (ry*2 .gt. by) ry = (by - ry) + 0.5*cly
            if (rz*2 .gt. bz) rz = (bz - rz) + 0.5*clz
            dr = sqrt(rx**2 + ry**2 + rz**2)
            if (dr .lt. 3.0) then
                cn(i) = cn(i) + 1
                cn(j) = cn(j) + 1
                nnsite(cn(i), i) = j
                nnsite(cn(j), j) = i
                if (site_type(i) .lt. 3 .and. site_type(j) .lt. 3) then
                    effcn(i) = effcn(i) + 1
                    effcn(j) = effcn(j) + 1
                end if
            end if
            if (dr .lt. 5.8) then
                cn_2(i) = cn_2(i) + 1
                cn_2(j) = cn_2(j) + 1
                nnnsite(cn_2(i), i) = j
                nnnsite(cn_2(j), j) = i
            end if
        end do
        if (restart .eq. 0) then
            if (effcn(i) .gt. 0 .and. effcn(i) .lt. 12) then
                site_type(i) = 2
                cov_type(i) = 0
            end if
        end if
        if (restart .eq. 1) then
            if (effcn(i) .gt. 0 .and. effcn(i) .lt. 12) site_type(i) = 2
        end if
    end do

    gcn = 0
    do i = 1, nbulk
        if (site_type(i) .lt. 3) then
            do j = 1, 12
                nn_j = nnsite(j, i)
                gcn(i) = gcn(i) + effcn(nn_j)
            end do
            gcn(i) = gcn(i)/12.0
        end if
    end do

    nsurf10 = 0
    do i = 1, nbulk
        if (site_type(i) .lt. 3 .and. effcn(i) .lt. 10) then
            nsurf10 = nsurf10 + 1
        end if
    end do
    write (*, *) 'surface atoms(<10):', int(nsurf10)

    nsurf12 = 0
    do i = 1, nbulk
        if (site_type(i) .lt. 3 .and. effcn(i) .lt. 12) then
            nsurf12 = nsurf12 + 1
        end if
    end do
    write (*, *) 'surface atoms(<12):', int(nsurf12)

!open(10,file='cn.txt')
!  do i = 1, nbulk
!    if(site_type(i).lt.3) then
!        write(10,*) effcn(i),site_type(i)
!    end if
!  end do
!close(10)

    open (11, file='100.xyz')
    read (11, *) n100
    read (11, *)
    do i = 1, int(n100)
        read (11, *) elem, xxx_100(i), yyy_100(i), zzz_100(i)
    end do
    close (11)

    open (11, file='100_new.xyz')
    write (11, *) int(n100)
    write (11, *) n100
    do i = 1, nbulk
        do j = 1, int(n100)
      if (abs(xxx(i) - xxx_100(j)) .lt. 0.01 .and. abs(yyy(i) - yyy_100(j)) .lt. 0.01 &
        .and. abs(zzz(i) - zzz_100(j)) .lt. 0.01) then
            write (11, *) elem, xxx(i), yyy(i), zzz(i)
            exit
        end if
        end do
    end do
    close (11)

    open (11, file='110.xyz')
    read (11, *) n110
    read (11, *)
    do i = 1, int(n110)
        read (11, *) elem, xxx_110(i), yyy_110(i), zzz_110(i)
    end do
    close (11)

    open (11, file='110_new.xyz')
    write (11, *) int(n110)
    write (11, *) n110
    do i = 1, nbulk
        do j = 1, int(n110)
      if (abs(xxx(i) - xxx_110(j)) .lt. 0.01 .and. abs(yyy(i) - yyy_110(j)) .lt. 0.01 &
        .and. abs(zzz(i) - zzz_110(j)) .lt. 0.01) then
            xxx_110(j) = xxx(i)
            yyy_110(j) = yyy(i)
            zzz_110(j) = zzz(i)
            write (11, *) elem, xxx(i), yyy(i), zzz(i), site_type(i)
            exit
        end if
        end do
    end do
    close (11)

    open (11, file='111.xyz')
    read (11, *) n111
    read (11, *)
    do i = 1, int(n111)
        read (11, *) elem, xxx_111(i), yyy_111(i), zzz_111(i)
    end do
    close (11)

    open (11, file='111_new.xyz')
    write (11, *) int(n111)
    write (11, *) n111
    do i = 1, nbulk
        do j = 1, int(n111)
      if (abs(xxx(i) - xxx_111(j)) .lt. 0.01 .and. abs(yyy(i) - yyy_111(j)) .lt. 0.01 &
        .and. abs(zzz(i) - zzz_111(j)) .lt. 0.01) then
            write (11, *) elem, xxx(i), yyy(i), zzz(i)
            exit
        end if
        end do
    end do
    close (11)

    open (11, file='edge.xyz')
    read (11, *) nedge
    read (11, *)
    do i = 1, int(nedge)
        read (11, *) elem, xxx_edge(i), yyy_edge(i), zzz_edge(i)
    end do
    close (11)

    open (11, file='edge_new.xyz')
    write (11, *) int(nedge)
    write (11, *) nedge
    do i = 1, nbulk
        do j = 1, int(nedge)
   if (abs(xxx(i) - xxx_edge(j)) .lt. 0.01 .and. abs(yyy(i) - yyy_edge(j)) .lt. 0.01 &
    .and. abs(zzz(i) - zzz_edge(j)) .lt. 0.01) then
                write (11, *) elem, xxx(i), yyy(i), zzz(i)
                exit
            end if
        end do
    end do
    close (11)

    ! calculate rtot, rsite, rneis, revent
    rtot = 0.0
    rsite = 0.0
    revent = 0.0
    rneis = 0.0
    Ea_r = 0.0
    ijpick_cn = 0
    ijpick_gcn = 0.0
    do i = 1, nbulk
        if (site_type(i) .eq. 2) then
            do k = 1, 8
                if (k .ne. 2 .and. k .ne. 3) then
                    do j = 1, 12
                        r_site = 0
                        call rijk(r_site, i, k, j, effcn, gcn, site_type, cov_type, &
                                  nnsite, nnnsite, Ea_r, ijpick_cn, ijpick_gcn, Temperature, pCO, pO2)
                        rneis(j, k, i) = r_site
                        revent(k, i) = revent(k, i) + r_site
                    end do
                end if
                if (k .eq. 2 .or. k .eq. 3) then
                    r_site = 0
                    call rijk(r_site, i, k, 0, effcn, gcn, site_type, cov_type, &
                              nnsite, nnnsite, Ea_r, ijpick_cn, ijpick_gcn, Temperature, pCO, pO2)
                    revent(k, i) = r_site
                end if
                rsite(i) = rsite(i) + revent(k, i)
            end do
            rtot = rtot + rsite(i)
        end if
    end do

    call init_random_seed()
    n_step_tot = 0
    i_step = 0
    ii_step = 0
    COads_step = 0
    COdes_step = 0
    Oads_step = 0
    Odes_step = 0
    COdiff_step = 0
    Odiff_step = 0
    COOx_step = 0
    ctime = 0.0
    ! new files
    nf_1 = 0
    write (char_seq, '(i2.2)') nf_1
    write (filename1, *) 'atom_str_', trim(char_seq), '.xyz'  !real atoms coordination
    ! write(filename2, *) 'ijk_rec_',trim(char_seq),'.dat'   !ctime, ipick, jpick, kick of every step
    ! write(filename3, *) 'Ea_atom_',trim(char_seq),'.dat'   !energy barrier of atom jumping, cn, n_tot_step
    write (filename4, *) 'Ea_COO_', trim(char_seq), '.dat'    !energy barrier of CO oxdaition, cn, n_tot_step
    write (filename5, *) 'step_rec_', trim(char_seq), '.dat'  !steps of every event
    write (filename6, *) 'coverage_', trim(char_seq), '.dat'
    write (filename7, *) 'coverage110_', trim(char_seq), '.dat'
    write (filename8, *) 'coverage100_', trim(char_seq), '.dat'
    write (filename9, *) 'coverage111_', trim(char_seq), '.dat'
    write (filename10, *) 'coveragedge_', trim(char_seq), '.dat'

    open (10, file=filename1, status='new', action='write')
    ! open(11, file=filename2, status='new', action='write')
    ! open(12, file=filename3, status='new', action='write')
    open (13, file=filename4, status='new', action='write')
    open (14, file=filename5, status='new', action='write')
    open (15, file=filename6, status='new', action='write')
    open (16, file=filename7, status='new', action='write')
    open (17, file=filename8, status='new', action='write')
    open (18, file=filename9, status='new', action='write')
    open (19, file=filename10, status='new', action='write')

    ! record first structural coordination
    write (10, *) int(nsurf12)
    write (10, *) i_step, n_step_tot, ctime
    do i = 1, nbulk
        if (site_type(i) .lt. 3 .and. cov_type(i) .eq. 0) then
            write (10, *) 'Pt', xxx(i), yyy(i), zzz(i), cov_type(i)
        end if
    end do
    do i = 1, nbulk
        if (site_type(i) .lt. 3 .and. cov_type(i) .eq. 1) then
            write (10, *) 'Au', xxx(i), yyy(i), zzz(i), cov_type(i)
        end if
    end do
    do i = 1, nbulk
        if (site_type(i) .lt. 3 .and. cov_type(i) .eq. 2) then
            write (10, *) 'Ir', xxx(i), yyy(i), zzz(i), cov_type(i)
        end if
    end do

    COOx_100 = 0
    COOx_110 = 0
    COOx_111 = 0
    COOx_edge = 0

    ! KCM cycle
    do while (n_step_tot .le. nLoop)
        ! Circular recording
        if (n_step_tot .gt. 0 .and. mod(n_step_tot, 200000000) .eq. 0) then
            nf_1 = nf_1 + 1
            write (char_seq, '(i2.2)') nf_1
            write (filename1, *) 'atom_str_', trim(char_seq), '.xyz'  !real atoms coordination
            ! write(filename2, *) 'ijk_COO_rec_',trim(char_seq),'.dat'   !ctime, ipick, jpick, kick of every step
            ! write(filename3, *) 'Ea_atom_',trim(char_seq),'.dat'   !energy barrier of atom jumping, cn, n_tot_step
            write (filename4, *) 'Ea_COO_', trim(char_seq), '.dat'    !energy barrier of CO oxdaition, cn, n_tot_step
            write (filename5, *) 'step_rec_', trim(char_seq), '.dat'  !steps of every event
            write (filename6, *) 'coverage_', trim(char_seq), '.dat'
            write (filename7, *) 'coverage110_', trim(char_seq), '.dat'
            write (filename8, *) 'coverage100_', trim(char_seq), '.dat'
            write (filename9, *) 'coverage111_', trim(char_seq), '.dat'
            write (filename10, *) 'coveragedge_', trim(char_seq), '.dat'

            open (10, file=filename1, status='new', action='write')
            !open(11, file=filename2, status='new', action='write')
            ! open(12, file=filename3, status='new', action='write')
            open (13, file=filename4, status='new', action='write')
            open (14, file=filename5, status='new', action='write')
            open (15, file=filename6, status='new', action='write')
            open (16, file=filename7, status='new', action='write')
            open (17, file=filename8, status='new', action='write')
            open (18, file=filename9, status='new', action='write')
            open (19, file=filename10, status='new', action='write')
        end if

        ! cal time
        call random_number(randum)
        dtime = -dlog(randum)/rtot
        ctime = ctime + dtime
        ! ipick
        do while (1 .gt. 0)
            call random_number(randum)
            rpoint = rtot*randum
            sum_rsite = 0.0
            do i = 1, nbulk
                if (site_type(i) .eq. 2) then
                    sum_rsite = sum_rsite + rsite(i)
                    if (sum_rsite .ge. rpoint) goto 100
                end if
            end do
        end do
100     ipick = i

        ! event
        !do while (1 .gt. 0)
        call random_number(randum)
        rpoint = rsite(ipick)*randum
        sum_revent = 0
        do k = 1, 8
            sum_revent = sum_revent + revent(k, ipick)
            if (sum_revent .ge. rpoint) goto 200
        end do
        !end do
200     kpick = k

        ! j_cn & jpick
        if (kpick .ne. 2 .and. kpick .ne. 3) then
            !do while(1.gt.0)
            call random_number(randum)
            rpoint = revent(kpick, ipick)*randum
            sum_rneis = 0
            do j = 1, 12
                sum_rneis = sum_rneis + rneis(j, kpick, ipick)
                if (sum_rneis .ge. rpoint) goto 300
            end do
            !end do
300         j_cn = j
            jpick = nnsite(j_cn, ipick)
        end if

        n_step_tot = n_step_tot + 1

        if (kpick .eq. 8) then
            write (13, *) n_step_tot, ctime, Ea_r(2, j_cn, ipick), ipick, jpick, &
                effcn(ipick), effcn(jpick), gcn(ipick), gcn(jpick)

            do j = 1, int(n100)
                if (abs(xxx(ipick) - xxx_100(j)) .lt. 0.01 .and. abs(yyy(ipick) - yyy_100(j)) .lt. 0.01 .and. &
                    abs(zzz(ipick) - zzz_100(j)) .lt. 0.01) then
                    COOx_100 = COOx_100 + 1
                end if
                if (abs(xxx(jpick) - xxx_100(j)) .lt. 0.01 .and. abs(yyy(jpick) - yyy_100(j)) .lt. 0.01 .and. &
                    abs(zzz(jpick) - zzz_100(j)) .lt. 0.01) then
                    COOx_100 = COOx_100 + 1
                end if
            end do

            do j = 1, int(n110)
                if (abs(xxx(ipick) - xxx_110(j)) .lt. 0.01 .and. abs(yyy(ipick) - yyy_110(j)) .lt. 0.01 .and. &
                    abs(zzz(ipick) - zzz_110(j)) .lt. 0.01) then
                    COOx_110 = COOx_110 + 1
                end if
                if (abs(xxx(jpick) - xxx_110(j)) .lt. 0.01 .and. abs(yyy(jpick) - yyy_110(j)) .lt. 0.01 .and. &
                    abs(zzz(jpick) - zzz_110(j)) .lt. 0.01) then
                    COOx_110 = COOx_110 + 1
                end if
            end do

            do j = 1, int(n111)
                if (abs(xxx(ipick) - xxx_111(j)) .lt. 0.01 .and. abs(yyy(ipick) - yyy_111(j)) .lt. 0.01 .and. &
                    abs(zzz(ipick) - zzz_111(j)) .lt. 0.01) then
                    COOx_111 = COOx_111 + 1
                end if
                if (abs(xxx(jpick) - xxx_111(j)) .lt. 0.01 .and. abs(yyy(jpick) - yyy_111(j)) .lt. 0.01 .and. &
                    abs(zzz(jpick) - zzz_111(j)) .lt. 0.01) then
                    COOx_111 = COOx_111 + 1
                end if
            end do

            do j = 1, int(nedge)
                if (abs(xxx(ipick) - xxx_edge(j)) .lt. 0.01 .and. abs(yyy(ipick) - yyy_edge(j)) .lt. 0.01 .and. &
                    abs(zzz(ipick) - zzz_edge(j)) .lt. 0.01) then
                    COOx_edge = COOx_edge + 1
                end if
                if (abs(xxx(jpick) - xxx_edge(j)) .lt. 0.01 .and. abs(yyy(jpick) - yyy_edge(j)) .lt. 0.01 .and. &
                    abs(zzz(jpick) - zzz_edge(j)) .lt. 0.01) then
                    COOx_edge = COOx_edge + 1
                end if
            end do
        end if

        if (n_step_tot .gt. 0 .and. mod(n_step_tot, 100000000) .eq. 0) write (*, *) n_step_tot, ctime, &
            COOx_100, COOx_110, COOx_111, COOx_edge

        ! update site_type, cov_type, effcn, gcn  do_event
        if (kpick .eq. 1) then
            i_step = i_step + 1
            ii_step = ii_step + 1
            ni_ads = cov_type(ipick)
            site_type(jpick) = 2
            cov_type(jpick) = ni_ads
            site_type(ipick) = 3
            cov_type(ipick) = 5
            ! find all atoms of 2 neis of ipick and jpick
            ! and remove repeat
            natoms_ij = 0
            do i = 1, 54
                natoms_ij(i) = nnnsite(i, ipick)
            end do
            do j = 1, 54
                nj_site = nnnsite(j, jpick)
                nj = j + 54
                natoms_ij(nj) = nj_site
                do i = 1, 54
                    ni_site = nnnsite(i, ipick)
                    if (nj_site .eq. ni_site) then
                        natoms_ij(nj) = 0
                        exit
                    end if
                end do
            end do
            ! update effcn and gcn
            do i = 1, 108
                nn_i = natoms_ij(i)
                if (nn_i .gt. 0) then
                    effcn(nn_i) = 0
                    gcn(nn_i) = 0
                    if (site_type(nn_i) .lt. 3) then
                        ni_effcn = 0
                        gcn_ni = 0.0
                        do j = 1, 12
                            nn_j = nnsite(j, nn_i)
                            if (site_type(nn_j) .lt. 3) then
                                ni_effcn = ni_effcn + 1
                                do k = 1, 12
                                    nn_k = nnsite(k, nn_j)
                                    if (site_type(nn_k) .lt. 3) gcn_ni = gcn_ni + 1.0
                                end do
                            end if
                        end do
                        effcn(nn_i) = ni_effcn
                        gcn(nn_i) = gcn_ni/12.0
                        ! whether there is adsorbate in the Site where cn >= 10
                        if (cov_type(nn_i) .eq. 1 .and. ni_effcn .ge. 10 .and. ni_effcn .lt. 12) then
                            cov_type(nn_i) = 0
                            n_step_tot = n_step_tot + 1
                            COdes_step = COdes_step + 1
                            !write(11, *) ctime, nn_i, 0, 3
                        end if
                    end if
                end if
            end do
            ! update site_type and cov_type
            do i = 1, 12
                ! ipick, body -> surface
                neis_i = nnsite(i, ipick)
                if (site_type(neis_i) .eq. 1 .and. effcn(neis_i) .lt. 12) then
                    site_type(neis_i) = 2
                    cov_type(neis_i) = 0
                end if
                ! jpick, surface -> body
                neis_j = nnsite(i, jpick)
                if (site_type(neis_j) .eq. 2 .and. effcn(neis_j) .eq. 12) then
                    site_type(neis_j) = 1
                    cov_type(neis_j) = 5
                end if
            end do
            !site_type(jpick) = 2
            !cov_type(jpick) = ni_ads
            !site_type(ipick) = 3
            !cov_type(ipick) = 5
        end if

        ! CO adsorption
        if (kpick .eq. 2) then
            cov_type(ipick) = 1
            COads_step = COads_step + 1
        end if

        ! CO desorption
        if (kpick .eq. 3) then
            cov_type(ipick) = 0
            COdes_step = COdes_step + 1
        end if

        ! O2 adsorption
        if (kpick .eq. 4) then
            cov_type(ipick) = 2
            cov_type(jpick) = 2
            Oads_step = Oads_step + 1
        end if

        ! O2 desorption
        if (kpick .eq. 5) then
            cov_type(ipick) = 0
            cov_type(jpick) = 0
            Odes_step = Odes_step + 1
        end if

        ! CO diffusion
        if (kpick .eq. 6) then
            cov_type(ipick) = 0
            cov_type(jpick) = 1
            COdiff_step = COdiff_step + 1
            jpick_effcn = effcn(jpick)
            ! if(jpick_effcn.ge.10) then
            ! cov_type(jpick) = 0
            ! n_step_tot = n_step_tot + 1
            ! COdes_step = COdes_step + 1
            ! !write(11, *) ctime, jpick, 0, 3
            ! end if
        end if

        ! O diffusion
        if (kpick .eq. 7) then
            cov_type(ipick) = 0
            cov_type(jpick) = 2
            Odiff_step = Odiff_step + 1
        end if

        ! COOx reaction
        if (kpick .eq. 8) then
            cov_type(ipick) = 0
            cov_type(jpick) = 0
            COOx_step = COOx_step + 1
            !write(11, *) ctime, ipick, jpick
        end if

        if (mod(n_step_tot, 10000000) .eq. 0) then
            write (10, *) int(nsurf12)
            write (10, *) i_step, n_step_tot, ctime
            do i = 1, nbulk
                if (site_type(i) .lt. 3 .and. cov_type(i) .eq. 0) then
                    write (10, *) 'Pt', xxx(i), yyy(i), zzz(i), cov_type(i)
                end if
            end do
            do i = 1, nbulk
                if (site_type(i) .lt. 3 .and. cov_type(i) .eq. 1) then
                    write (10, *) 'Au', xxx(i), yyy(i), zzz(i), cov_type(i)
                end if
            end do
            do i = 1, nbulk
                if (site_type(i) .lt. 3 .and. cov_type(i) .eq. 2) then
                    write (10, *) 'Ir', xxx(i), yyy(i), zzz(i), cov_type(i)
                end if
            end do
!          do i = 1, nbulk
!                if(site_type(i).lt.3 .and. cov_type(i) .eq. 5) then
!                  write(10, *) 'Co', xxx(i), yyy(i), zzz(i), cov_type(i)
!                end if
!          end do
        end if

        if (mod(n_step_tot, 100000) .eq. 0) then
            nsurf_CO = 0
            nsurf_O = 0
            do i = 1, nbulk
                if (site_type(i) .lt. 3 .and. cov_type(i) .eq. 1) nsurf_CO = nsurf_CO + 1
                if (site_type(i) .lt. 3 .and. cov_type(i) .eq. 2) nsurf_O = nsurf_O + 1
            end do
            write (15, *) ctime, n_step_tot, nsurf_CO, nsurf_O, nsurf_CO/nsurf12, nsurf_O/nsurf12

            write (14, *) ctime, n_step_tot, i_step, COads_step, &
                COdes_step, Oads_step, Odes_step, COdiff_step, Odiff_step, COOx_step

            n110_CO = 0
            n110_O = 0
            do i = 1, nbulk
                do j = 1, int(n110)
      if (abs(xxx(i) - xxx_110(j)) .lt. 0.01 .and. abs(yyy(i) - yyy_110(j)) .lt. 0.01 &
        .and. abs(zzz(i) - zzz_110(j)) .lt. 0.01) then
                        if (cov_type(i) .eq. 1) n110_CO = n110_CO + 1
                        if (cov_type(i) .eq. 2) n110_O = n110_O + 1
                    end if
                end do
            end do
            write (16, *) ctime, n_step_tot, n110_CO, n110_O, n110_CO/n110, n110_O/n110

            n100_CO = 0
            n100_O = 0
            do i = 1, nbulk
                do j = 1, int(n100)
      if (abs(xxx(i) - xxx_100(j)) .lt. 0.01 .and. abs(yyy(i) - yyy_100(j)) .lt. 0.01 &
        .and. abs(zzz(i) - zzz_100(j)) .lt. 0.01) then
                        if (cov_type(i) .eq. 1) n100_CO = n100_CO + 1
                        if (cov_type(i) .eq. 2) n100_O = n100_O + 1
                    end if
                end do
            end do
            write (17, *) ctime, n_step_tot, n100_CO, n100_O, n100_CO/n100, n100_O/n100

            n111_CO = 0
            n111_O = 0
            do i = 1, nbulk
                do j = 1, int(n111)
      if (abs(xxx(i) - xxx_111(j)) .lt. 0.01 .and. abs(yyy(i) - yyy_111(j)) .lt. 0.01 &
        .and. abs(zzz(i) - zzz_111(j)) .lt. 0.01) then
                        if (cov_type(i) .eq. 1) n111_CO = n111_CO + 1
                        if (cov_type(i) .eq. 2) n111_O = n111_O + 1
                    end if
                end do
            end do
            write (18, *) ctime, n_step_tot, n111_CO, n111_O, n111_CO/n111, n111_O/n111

            nedge_CO = 0
            nedge_O = 0
            do i = 1, nbulk
                do j = 1, int(nedge)
   if (abs(xxx(i) - xxx_edge(j)) .lt. 0.01 .and. abs(yyy(i) - yyy_edge(j)) .lt. 0.01 &
    .and. abs(zzz(i) - zzz_edge(j)) .lt. 0.01) then
                        if (cov_type(i) .eq. 1) nedge_CO = nedge_CO + 1
                        if (cov_type(i) .eq. 2) nedge_O = nedge_O + 1
                    end if
                end do
            end do
            write (19, *) ctime, n_step_tot, nedge_CO, nedge_O, nedge_CO/nedge, nedge_O/nedge

        end if

        ! initial and cal  new relative r
        ! 2neis of ipick
        do i = 1, 54
            nn_i = nnnsite(i, ipick)
            if (nn_i .gt. 0) then
                rtot = rtot - rsite(nn_i)
                rsite(nn_i) = 0
                do k = 1, 8
                    revent(k, nn_i) = 0
                    if (k .ne. 2 .and. k .ne. 3) then
                        do j = 1, 12
                            rneis(j, k, nn_i) = 0
                        end do
                    end if
                end do
                if (site_type(nn_i) .eq. 2) then
                    do k = 1, 8
                        if (k .ne. 2 .and. k .ne. 3) then
                            do j = 1, 12
                                r_site = 0
                                call rijk(r_site, nn_i, k, j, effcn, gcn, site_type, cov_type, &
                                          nnsite, nnnsite, Ea_r, ijpick_cn, ijpick_gcn, Temperature, pCO, pO2)
                                rneis(j, k, nn_i) = r_site
                                revent(k, nn_i) = revent(k, nn_i) + r_site
                            end do
                        end if
                        if (k .eq. 2 .or. k .eq. 3) then
                            r_site = 0
                            call rijk(r_site, nn_i, k, 0, effcn, gcn, site_type, cov_type, &
                                      nnsite, nnnsite, Ea_r, ijpick_cn, ijpick_gcn, Temperature, pCO, pO2)
                            revent(k, nn_i) = revent(k, nn_i) + r_site
                        end if
                        rsite(nn_i) = rsite(nn_i) + revent(k, nn_i)
                    end do
                    rtot = rtot + rsite(nn_i)
                end if
            end if
        end do

        ! ipick site
        if (kpick .eq. 2 .or. kpick .eq. 3) then
            rtot = rtot - rsite(ipick)
            rsite(ipick) = 0
            do k = 1, 8
                ! initial
                revent(k, ipick) = 0
                if (k .ne. 2 .and. k .ne. 3) then
                    do j = 1, 12
                        rneis(j, k, ipick) = 0
                    end do
                end if
                ! cal new r
                if (site_type(ipick) .eq. 2) then
                if (k .ne. 2 .and. k .ne. 3) then
                    do j = 1, 12
                        r_site = 0
                        call rijk(r_site, ipick, k, j, effcn, gcn, site_type, cov_type, &
                                  nnsite, nnnsite, Ea_r, ijpick_cn, ijpick_gcn, Temperature, pCO, pO2)
                        rneis(j, k, ipick) = r_site
                        revent(k, ipick) = revent(k, ipick) + r_site
                    end do
                end if
                if (k .eq. 2 .or. k .eq. 3) then
                    r_site = 0
                    call rijk(r_site, ipick, k, 0, effcn, gcn, site_type, cov_type, &
                              nnsite, nnnsite, Ea_r, ijpick_cn, ijpick_gcn, Temperature, pCO, pO2)
                    revent(k, ipick) = revent(k, ipick) + r_site
                end if
                rsite(ipick) = rsite(ipick) + revent(k, ipick)
                end if
            end do
            rtot = rtot + rsite(ipick)
        end if

        ! 2neis of jpick
        if (kpick .ne. 2 .and. kpick .ne. 3) then
            do i = 1, 54
                nn_i = nnnsite(i, jpick)
                if (nn_i .gt. 0) then
                    ! initial
                    rtot = rtot - rsite(nn_i)
                    rsite(nn_i) = 0
                    do k = 1, 8
                        revent(k, nn_i) = 0
                        if (k .ne. 2 .and. k .ne. 3) then
                            do j = 1, 12
                                rneis(j, k, nn_i) = 0
                            end do
                        end if
                    end do
                    ! cal new r
                    if (site_type(nn_i) .eq. 2) then
                    do k = 1, 8
                        if (k .ne. 2 .and. k .ne. 3) then
                            do j = 1, 12
                                r_site = 0
                                call rijk(r_site, nn_i, k, j, effcn, gcn, site_type, cov_type, &
                                          nnsite, nnnsite, Ea_r, ijpick_cn, ijpick_gcn, Temperature, pCO, pO2)
                                rneis(j, k, nn_i) = r_site
                                revent(k, nn_i) = revent(k, nn_i) + r_site
                            end do
                        end if
                        if (k .eq. 2 .or. k .eq. 3) then
                            r_site = 0
                            call rijk(r_site, nn_i, k, 0, effcn, gcn, site_type, cov_type, &
                                      nnsite, nnnsite, Ea_r, ijpick_cn, ijpick_gcn, Temperature, pCO, pO2)
                            revent(k, nn_i) = revent(k, nn_i) + r_site
                        end if
                        rsite(nn_i) = rsite(nn_i) + revent(k, nn_i)
                    end do
                    rtot = rtot + rsite(nn_i)
                    end if
                end if
            end do
        end if

    end do

    close (10)
!  close(11)
!  close(12)
    close (13)
    close (14)
    close (15)
    close (16)

    write (*, *) "====== END ======="
    ! write last_one (xyz and cov_type)
    open (11, file='last_one.xyz')
    write (11, *) natoms
    write (11, *) i_step, n_step_tot, ctime
    do i = 1, nbulk
        if (site_type(i) .lt. 3 .and. cov_type(i) .eq. 0) then
            write (11, *) 'Pt', xxx(i), yyy(i), zzz(i), cov_type(i)
        end if
    end do
    do i = 1, nbulk
        if (site_type(i) .lt. 3 .and. cov_type(i) .eq. 1) then
            write (11, *) 'Au', xxx(i), yyy(i), zzz(i), cov_type(i)
        end if
    end do
    do i = 1, nbulk
        if (site_type(i) .lt. 3 .and. cov_type(i) .eq. 2) then
            write (11, *) 'Ir', xxx(i), yyy(i), zzz(i), cov_type(i)
        end if
    end do
    do i = 1, nbulk
        if (site_type(i) .lt. 3 .and. cov_type(i) .eq. 5) then
            write (11, *) 'Co', xxx(i), yyy(i), zzz(i), cov_type(i)
        end if
    end do
    close (11)

end

subroutine rijk(r_site, ipick, kpick, j_cn, effcn, gcn, &
                site_type, cov_type, nnsite, nnnsite, Ea_r, &
                ijpick_cn, ijpick_gcn, Temp, pCO, pO2)
    implicit none

    integer(kind=4):: effcn(100000), site_type(100000), cov_type(100000), nj, &
                      ni_ads, nj_ads, ijpick_cn(2, 12, 100000)
    integer(kind=4):: i, j, k
    integer(kind=4):: nnsite(12, 100000), nnnsite(54, 100000), ipick, jpick, kpick, natoms_ij(108)
    integer(kind=4):: j_cn, ni_site, nj_site, nn_i, nnei_i, nnei_j, neis_jj, effcn_ipick, effcn_jpick
    real(kind=8):: gcn(100000), Ea_r(2, 12, 100000), gcn_ipick, gcn_jpick, gcn_ni, &
                   r_site, r_jump, rco_ads, rco_des, ro_ads, ro_des, rco_diff, ro_diff, r_coo

    real(kind=8), parameter::pi = 3.141592654
    real(kind=8), parameter::kb = 8.6173324D-05
    real(kind=8), parameter::h = 4.1356676D-15
    real(kind=8), parameter::eV2J = 1.60217662D-19
    real(kind=8), parameter::Na = 6.0221409D23
    real(kind=8), parameter::Asite = (10D-10)**2
    real(kind=4), parameter::Ui = -2.80
    real(kind=4), parameter::EdiffCO = 0.68
    real(kind=4), parameter::EdiffO = 0.50
    real(kind=8)::p0, pCO, pO2, Temp, s0CO, s0O, mCO, mO, Sco_0, So2_0
    real(kind=8)::U_ini, U_mid, U_fin, Eads_i, Eads_j, Eads_ij, Eads_ini, Eads_mid, Eads_fin, &
                   Eb_ini, Eb_fin, dE, Eb, Ef, Ea
    real(kind=8)::dS_CO, dS_O2, r_K_eq
    real(kind=8):: ijpick_gcn(2, 12, 100000)

    p0 = 100000

    s0CO = 0.9
    s0O = 0.1

    mCO = 28.01D-3/Na
    mO = 15.999D-3/Na

    Sco_0 = 85.142*(Temp**0.14709)/(Na*eV2J)
    So2_0 = 89.655*(Temp**0.14489)/(Na*eV2J)

    ! CO adsorption
    if (kpick .eq. 2 .and. cov_type(ipick) .eq. 0) then
        gcn_ipick = gcn(ipick)
        if (gcn_ipick .ge. 0 .and. gcn_ipick .lt. 5.33) s0CO = 1.0
        if (effcn(ipick) .ge. 10) s0CO = 0.0
        rco_ads = (s0CO*pCO*Asite)/sqrt(2.*pi*mCO*kb*eV2J*Temp)
        r_site = rco_ads
    end if

    ! CO desorption
    if (kpick .eq. 3 .and. cov_type(ipick) .eq. 1) then
        gcn_ipick = gcn(ipick)
        ! CO ads
        if (gcn_ipick .ge. 0 .and. gcn_ipick .lt. 5.33) s0CO = 1.0
        if (effcn(ipick) .ge. 10) s0CO = 0.0
        rco_ads = (s0CO*pCO*Asite)/sqrt(2.*pi*mCO*kb*eV2J*Temp)
        ! CO equilibrium
        dS_CO = 0 - (Sco_0 - kb*log(pCO/p0))
        Eads_i = 0.0
        call Eads_site(Eads_i, ipick, gcn_ipick, effcn, site_type, cov_type, nnsite, 1)
        r_K_eq = exp(-(Eads_i - Temp*dS_CO)/(kb*Temp))
        ! r of CO des
        rco_des = rco_ads/(pCO*r_K_eq)
        r_site = rco_des
    end if

    if (kpick .ne. 2 .and. kpick .ne. 3) then
        jpick = nnsite(j_cn, ipick)
        ni_ads = cov_type(ipick)
        nj_ads = cov_type(jpick)
        gcn_ipick = gcn(ipick)
        gcn_jpick = gcn(jpick)

        ! O2 adsorption
        if (kpick .eq. 4 .and. ni_ads .eq. 0 .and. nj_ads .eq. 0) then
            if (gcn_ipick .ge. gcn_jpick) then
                if (gcn_ipick .ge. 0 .and. gcn_ipick .lt. 5.33) s0O = 1.0
                if (effcn(ipick) .ge. 10) s0O = 0.0
            else
                if (gcn_jpick .ge. 0 .and. gcn_jpick .lt. 5.33) s0O = 1.0
                if (effcn(jpick) .ge. 10) s0O = 0.0
            end if
            ro_ads = (s0O*pO2*Asite)/sqrt(2.*pi*mO*kb*eV2J*Temp)
            r_site = ro_ads
        end if

        ! O2 desorption
        if (kpick .eq. 5 .and. ni_ads .eq. 2 .and. nj_ads .eq. 2) then
            ! O2 ads
            if (gcn_ipick .ge. gcn_jpick) then
                if (gcn_ipick .ge. 0 .and. gcn_ipick .lt. 5.33) s0O = 1.0
                if (effcn(ipick) .ge. 10) s0O = 0.0
            else
                if (gcn_jpick .ge. 0 .and. gcn_jpick .lt. 5.33) s0O = 1.0
                if (effcn(jpick) .ge. 10) s0O = 0.0
            end if
            ro_ads = (s0O*pO2*Asite)/sqrt(2.*pi*mO*kb*eV2J*Temp)
            ! O2 equilibrium
            dS_O2 = 0 - (So2_0 - kb*log(pO2/p0))
            Eads_i = 0.0
            call Eads_site(Eads_i, ipick, gcn_ipick, effcn, site_type, cov_type, nnsite, 1)
            Eads_j = 0.0
            call Eads_site(Eads_j, jpick, gcn_jpick, effcn, site_type, cov_type, nnsite, 1)
            Eads_ij = Eads_i + Eads_j
            r_K_eq = exp(-(Eads_ij - Temp*dS_O2)/(kb*Temp))
            ! O2 des
            ro_des = ro_ads/(pO2*r_K_eq)
            r_site = ro_des
        end if

        ! CO diffusion
        if (kpick .eq. 6 .and. ni_ads .eq. 1 .and. nj_ads .eq. 0) then
            Eads_i = 0.0
            call Eads_site(Eads_i, ipick, gcn_ipick, effcn, site_type, cov_type, nnsite, 1)
            cov_type(ipick) = 0
            cov_type(jpick) = 1
            Eads_j = 0.0
            call Eads_site(Eads_j, jpick, gcn_jpick, effcn, site_type, cov_type, nnsite, 1)
            dE = Eads_j - Eads_i
            if (dE .le. 0) dE = 0.0
            Ea = dE + EdiffCO
            rco_diff = (kb*Temp/h)*exp(-Ea/(kb*Temp))
            r_site = rco_diff
            ! recovery
            cov_type(ipick) = 1
            cov_type(jpick) = 0
            if (effcn(jpick) .ge. 10) r_site = 0.0
        end if

        ! O diffusion
        if (kpick .eq. 7 .and. ni_ads .eq. 2 .and. nj_ads .eq. 0) then
            Eads_i = 0.0
            call Eads_site(Eads_i, ipick, gcn_ipick, effcn, site_type, cov_type, nnsite, 1)
            cov_type(ipick) = 0
            cov_type(jpick) = 2
            Eads_j = 0.0
            call Eads_site(Eads_j, jpick, gcn_jpick, effcn, site_type, cov_type, nnsite, 1)
            dE = Eads_j - Eads_i
            if (dE .le. 0) dE = 0.0
            Ea = dE + EdiffO
            ro_diff = (kb*Temp/h)*exp(-Ea/(kb*Temp))
            r_site = ro_diff
            ! recovery
            cov_type(ipick) = 2
            cov_type(jpick) = 0
        end if

        ! COOx reaction
        if (kpick .eq. 8 .and. ni_ads .eq. 1 .and. nj_ads .eq. 2) then
            Eads_i = 0.0
            call Eads_site(Eads_i, ipick, gcn_ipick, effcn, site_type, cov_type, nnsite, 1)
            Eads_j = 0.0
            call Eads_site(Eads_j, jpick, gcn_jpick, effcn, site_type, cov_type, nnsite, 1)
            Ea = (0.89 - 1)*(Eads_i + Eads_j) + 0.33
            if (Ea .le. 0) Ea = 0.0
            Ea_r(2, j_cn, ipick) = Ea
            r_coo = (kb*Temp/h)*exp(-Ea/(kb*Temp))
            r_site = r_coo
        end if
    end if
end subroutine

subroutine Eads_site(Eads_s, i_site, gcn_i, effcn, site_type, cov_type, nnsite, num)

    implicit none

    integer(kind=4):: i, j
    integer(kind=4):: site_type(100000), cov_type(100000), n_CO, n_O
    integer(kind=4):: effcn(100000), effcn_ni, effcn_nj
    integer(kind=4):: nnsite(12, 100000), i_site, nn_s, nn_j, num
    real(kind=8) Eads_s, gcn_i
    real(kind=4), parameter::Eco_co = -0.22
    real(kind=4), parameter::Eco_o = -0.18
    real(kind=4), parameter::Eo_o = -0.20

    n_CO = 0
    n_O = 0

    if (num .eq. 1) then
        effcn_ni = effcn(i_site)
        if (effcn_ni .gt. 6) then
            do i = 1, 12
                nn_s = nnsite(i, i_site)
                if (effcn(nn_s) .gt. 6) then
                    if (cov_type(nn_s) .eq. 1) n_CO = n_CO + 1
                    if (cov_type(nn_s) .eq. 2) n_O = n_O + 1
                end if
            end do
        end if
    elseif (num .eq. 2) then
        effcn_ni = 0
        do i = 1, 12
            nn_s = nnsite(i, i_site)
            if (site_type(nn_s) .lt. 3) effcn_ni = effcn_ni + 1
        end do
        if (effcn_ni .gt. 6) then
            do i = 1, 12
                nn_s = nnsite(i, i_site)
                if (cov_type(nn_s) .eq. 1 .or. cov_type(nn_s) .eq. 2) then
                    effcn_nj = 0
                    do j = 1, 12
                        nn_j = nnsite(j, nn_s)
                        if (site_type(nn_j) .lt. 3) effcn_nj = effcn_nj + 1
                    end do
                    if (effcn_nj .gt. 6) then
                        if (cov_type(nn_s) .eq. 1) n_CO = n_CO + 1
                        if (cov_type(nn_s) .eq. 2) n_O = n_O + 1
                    end if
                end if
            end do
        end if
    end if
    if (cov_type(i_site) .eq. 1) Eads_s = (0.19*gcn_i - 2.74) - (n_CO*Eco_co + n_O*Eco_o)
    if (cov_type(i_site) .eq. 2) Eads_s = (0.13*gcn_i - 1.82) - (n_CO*Eco_o + n_O*Eo_o)
end subroutine
