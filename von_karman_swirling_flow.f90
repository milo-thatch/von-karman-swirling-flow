program von_karman_swirling_flow

    use module_subroutines
    use module_algebra

    call read_input_file

    allocate(grid(n),Deltah(n))
    allocate(F(n),G(n),H(n),A(n),B(n))

    call generate_grid_uniform
    call generate_initial_guess

    call loop
        
    print*, 'Axial velocity at infinity H(inf) = ', H(n)
    print*, 'Torque exerted on the disk pi*B(0) = ', pi*B(1)
    ! print results
    filename = 'results_von_karman_swirling_flow.dat'
    open(unit=10, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        print *, "Error: Unable to open file ", trim(filename)
        stop
    end if
    write(10, '(A)') "eta   F   A   G   B   H"
    do j = 1, n
        write(10, '(9(F25.10, 3X))') grid(j), F(j), A(j), G(j), B(j), H(j)
    end do
    close(10, status='keep')


    deallocate(F,G,H,A,B)
    deallocate(grid,Deltah)

end program von_karman_swirling_flow