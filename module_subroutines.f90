module module_subroutines

    use module_algebra
    implicit none

    integer :: i, j, iter ! counters
    integer :: n, n_ext ! number of points
    real(kind=8) :: eta_max, error
    real(kind=8), parameter :: eta_ext = 7.0d0, tol = 1e-8, pi = 4*atan(1.0d0)
    real(kind=8), allocatable, dimension(:) :: grid, Deltah, F, G, H, A, B
    character(len=256) :: filename
    integer :: ios

    ! thomas algorithm
    real(kind=8), allocatable, dimension(:,:,:) :: invDeltamatrix
    real(kind=8), allocatable, dimension(:,:) :: matA, matB, matC, delta, vecw
    real(kind=8), allocatable, dimension(:) :: vecb

contains
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    !---------------------- read input file -------------------------------------------------
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    subroutine read_input_file

        character(len=20) :: key_input
        character(len=236) :: value_input
        character(len=256) :: line
        logical :: n_found, eta_max_found

        n = 0          ! Initialize output variables
        eta_max = 0.0
        n_found = .false.
        eta_max_found = .false.

        open(unit=10, file="input.in", status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error: Unable to open file ", trim(filename)
            stop
        end if

        do
            read(10, '(A)', iostat=ios) line
            if (ios /= 0) exit  ! Exit loop at end of file

            ! Remove leading and trailing spaces
            line = adjustl(line)

            ! Split key_input and value_input at the "="
            call split_key_value(line, key_input, value_input)

            if (key_input == "n") then
                read(value_input, *, iostat=ios) n
                if (ios /= 0) then
                    print *, "Error: Invalid integer value_input for 'n'."
                    stop
                end if
                n_found = .true.
            else if (key_input == "eta_max") then
                read(value_input, *, iostat=ios) eta_max
                if (ios /= 0) then
                    print *, "Error: Invalid real(8) value_input for 'eta_max'."
                    stop
                end if
                eta_max_found = .true.
            else
                print *, "Error: Unrecognized key_input '", trim(key_input), "' in input file."
                stop
            end if
        end do

        close(10)

        if (.not. n_found) then
            print *, "Error: Missing 'n' parameter in input file."
            stop
        end if

        if (.not. eta_max_found) then
            print *, "Error: Missing 'eta_max' parameter in input file."
            stop
        end if
    end subroutine read_input_file
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    !---------------------- split key value (parsing) ---------------------------------------
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    subroutine split_key_value(line, key_input, value_input)
        implicit none
        character(len=*), intent(in) :: line
        character(len=20), intent(out) :: key_input
        character(len=236), intent(out) :: value_input
        integer :: pos

        pos = index(line, "=")
        if (pos > 0) then
            key_input = adjustl(trim(line(1:pos-1)))
            value_input = adjustl(trim(line(pos+1:)))
        else
            key_input = ""
            value_input = ""
        end if
    end subroutine split_key_value
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    !---------------------- generate uniform grid -------------------------------------------
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    subroutine generate_grid_uniform

        ! initialize
        grid = 0.0d0
        Deltah = 0.0d0
        ! compute
        n_ext = int(eta_ext/eta_max*(n-1)) ! works out index n_ext to compute initial guess
        print*, 'n_ext = ', n_ext
        
        do j=1,n 
            grid(j) = real(j-1) * eta_max / real(n-1) 
        end do
        Deltah = eta_max / real(n-1) ! assign same spacing to all grid elements

    end subroutine generate_grid_uniform
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    !---------------------- initial guess base flow -----------------------------------------
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    subroutine generate_initial_guess

        F=0.0d0
        G=0.0d0
        A=0.0d0
        B=0.0d0

        ! create guesses for f, u, v, g, p
        do j=1,n_ext
            F(j) = 3.0d0/(2.0d0*eta_ext)*(grid(j)/eta_ext - (grid(j)/eta_ext)**2)
            A(j) = 0.0d0
            G(j) = 1.0d0 - 1.0d0/2.0d0*(grid(j)/eta_ext)*(3.0d0 - (grid(j)/eta_ext)**2)
            B(j) = 0.0d0
            H(j) = - 1.0d0/2.0d0*(grid(j)/eta_ext)*(3.0d0 - (grid(j)/eta_ext)**2)
        end do
        do j=n_ext+1,n
            F(j) = 0.0d0
            A(j) = 0.0d0
            G(j) = 0.0d0 
            B(j) = 0.0d0
            H(j) = - 1.0d0
        end do

        ! print initial guess
        filename = 'results_initial_guess.dat'
        open(unit=10, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            print *, "Error: Unable to open file ", trim(filename)
            stop
        end if
        write(10, '(A)') "eta   F   G   A   B   H"
        do j = 1, n
            write(10, '(9(F15.10, 3X))') grid(j), F(j), A(j), G(j), B(j), H(j)
        end do
        close(10, status='keep')

    end subroutine generate_initial_guess
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    !----------------------------- build matrix A -------------------------------------------
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    subroutine build_matA(i)

        integer :: i
        
        matA = 0.0d0 ! all the elements are set to zero
        
        ! BOTTOM BOUNDARY
        if(i==1) then
            matA(1,1) = 1.0d0
            matA(2,3) = 1.0d0
            matA(3,5) = 1.0d0
            matA(4,1) = - 1.0d0
            matA(4,2) = - Deltah(2)/2.0d0
            matA(5,3) = - 1.0d0
            matA(5,4) = - Deltah(2)/2.0d0 
            return
        end if

        ! CENTRAL POINTS (and top boundary)
        matA(1,1) = Deltah(i)
        matA(1,5) = 1.0d0
        matA(2,1) = Deltah(i)*F(i) 
        matA(2,2) = Deltah(i)*H(i)/2.0d0 - 1.0d0 
        matA(2,3) = - Deltah(i)*G(i)
        matA(2,5) = Deltah(i)/2.0d0*A(i)
        matA(3,1) = Deltah(i)*G(i)
        matA(3,3) = Deltah(i)*F(i)
        matA(3,4) = Deltah(i)/2.0d0*H(i) - 1.0d0 
        matA(3,5) = Deltah(i)/2.0d0*B(i)
        if(i==n) then ! top boundary
            matA(4,1) = 1.0d0
            matA(5,3) = 1.0d0
            return
        end if
        matA(4,1) = - 1.0d0
        matA(4,2) = - Deltah(i+1)/2.0d0
        matA(5,3) = - 1.0d0
        matA(5,4) = - Deltah(i+1)/2.0d0 
        return 

    end subroutine build_matA
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    !----------------------------- build matrix B -------------------------------------------
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    subroutine build_matB(i)

        integer :: i

        matB = 0.0d0

        ! CENTRAL POINTS
        matB(1,1) = Deltah(i)
        matB(1,5) = - 1.0d0
        matB(2,1) = Deltah(i)*F(i) 
        matB(2,2) = Deltah(i)*H(i)/2.0d0 + 1.0d0 
        matB(2,3) = - Deltah(i)*G(i)
        matB(2,5) = Deltah(i)/2.0d0*A(i)
        matB(3,1) = Deltah(i)*G(i)
        matB(3,3) = Deltah(i)*F(i)
        matB(3,4) = Deltah(i)/2.0d0*H(i) + 1.0d0 
        matB(3,5) = Deltah(i)/2.0d0*B(i)
        return

    end subroutine build_matB
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    !----------------------------- build matrix C -------------------------------------------
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    subroutine build_matC(i)
        
        integer :: i

        matC = 0.0d0

        matC(4,1) = 1.0d0
        matC(4,2) = - Deltah(i+1)/2.0d0
        matC(5,3) = 1.0d0
        matC(5,4) = - Deltah(i+1)/2.0d0 
        return 

    end subroutine build_matC
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    !----------------------------- build vector r -------------------------------------------
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    subroutine build_vecb(i)
        
        integer :: i
        
        vecb = 0.0d0

        if(i==1) then ! bottom boundary
            vecb(4) = F(1) - F(2) + Deltah(2)/2.0d0*(A(2) + A(1))
            vecb(5) = G(1) - G(2) + Deltah(2)/2.0d0*(B(2) + B(1))
            return
        else if(i==n) then ! top boundary
            vecb(1) = - Deltah(n)*(F(n) + F(n-1)) - H(n) + H(n-1)
            vecb(2) = - Deltah(n)/2.0d0*(F(n)**2.0d0 + F(n-1)**2.0d0 - G(n)**2.0d0 - G(n-1)**2.0d0) &
                        - Deltah(n)/2.0d0*(A(n)*H(n) + A(n-1)*H(n-1)) &
                        + A(n) - A(n-1)
            vecb(3) = - Deltah(n)*(F(n)*G(n) + F(n-1)*G(n-1)) &
                        - Deltah(n)/2.0d0*(B(n)*H(n) + B(n-1)*H(n-1)) &
                        + B(n) - B(n-1)    
            return
        end if
        ! internal points
        vecb(1) = - Deltah(i)*(F(i) + F(i-1)) - H(i) + H(i-1)
        vecb(2) = - Deltah(i)/2.0d0*(F(i)**2.0d0 + F(i-1)**2.0d0 - G(i)**2.0d0 - G(i-1)**2.0d0) &
                    - Deltah(i)/2.0d0*(A(i)*H(i) + A(i-1)*H(i-1)) &
                    + A(i) - A(i-1)
        vecb(3) = - Deltah(i)*(F(i)*G(i) + F(i-1)*G(i-1)) &
                    - Deltah(i)/2.0d0*(B(i)*H(i) + B(i-1)*H(i-1)) &
                    + B(i) - B(i-1)
        vecb(4) = F(i) - F(i+1) + Deltah(i+1)/2.0d0*(A(i+1) + A(i))
        vecb(5) = G(i) - G(i+1) + Deltah(i+1)/2.0d0*(B(i+1) + B(i))

    end subroutine build_vecb
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    !----------------------------- main loop ------------------------------------------------
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    subroutine loop

        allocate(matA(5,5),matB(5,5),matC(5,5),vecb(5))
        allocate(invDeltamatrix(5,5,n))
        allocate(vecw(5,n),delta(5,n))

        error = 1000.0d0
        iter = 0
        
        do while(error>=tol)

            call build_matA(1) 
            call inv_mat(matA,invDeltamatrix(:,:,1),5) 
            call build_vecb(1)
            vecw(:,1) = vecb(:)
            call build_matC(1) 
            
            ! UPWARD SWEEP
            do j=2,n 
                call build_matB(j)
                call build_matA(j)
                call build_matC(j-1)
                call inv_mat(matA - matmul(matmul(matB,invDeltamatrix(:,:,j-1)),matC),invDeltamatrix(:,:,j),5)
                call build_vecb(j)
                vecw(:,j) = vecb(:) - matmul(matmul(matB,invDeltamatrix(:,:,j-1)),vecw(:,j-1))
            end do
            delta(:,n) = matmul(invDeltamatrix(:,:,n),vecw(:,n))

            ! DOWNWARD SWEEP
            do j=n-1,1,-1
                call build_matC(j)
                delta(:,j) = matmul(invDeltamatrix(:,:,j),vecw(:,j) - matmul(matC,delta(:,j+1)))
            end do

            ! UPDATE guesses
            do j=1,n
                F(j) = F(j) + delta(1,j)
                A(j) = A(j) + delta(2,j)
                G(j) = G(j) + delta(3,j)
                B(j) = B(j) + delta(4,j)
                H(j) = H(j) + delta(5,j)
            end do
            error = abs(delta(1,1)) + abs(delta(2,1)) + abs(delta(3,1)) + abs(delta(4,1)) + abs(delta(5,1))
            print*, 'iter = ', iter,' error = ', error
            iter = iter+1
        end do

        deallocate(matA,matB,matC,vecb)
        deallocate(invDeltamatrix)
        deallocate(vecw,delta)
            
    end subroutine loop

end module module_subroutines