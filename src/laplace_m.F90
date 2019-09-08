module laplace_m

  use mpi
  use routines
  use random_m

  implicit none

  private

  public :: run_laplace


contains

  subroutine run_laplace()

    implicit none

    type(comm_t) :: comm
    integer :: ndim, comm3d, L1,L2, border(3,0:1), displ(0:1), direction, i
    integer :: n_dest, n_src, dest_temp, src_temp, src_dir(6), dest_dir(6), dest(6), src(6)
    logical :: isperiodic(3), reorder, converge
    real(8), allocatable :: v_local(:,:,:), ghost_send(:,:,:), ghost_recv(:,:,:)
    real(8) :: diff, diff_total, tol, itime, ftime, diff_old
    integer :: n_sendrecv, nstep, j, k
    real :: V1, V2
    character(4) :: file_name

    itime = mpi_wtime()


    V1 = 5.00
    V2 = 10.00
    L1 = 10
    L2 = 30

   ! Insert the number of processors for each dimension (nx,ny,nz) in comm%dims(1:3)
   ! the minimum number of processors in each direction must be 1
   ! Example: if you want to run the code in serial, you should write dims(1:3) = 1
   ! so you will have a virtual topology only with the point (0,0,0) related with one processor
   ! the number of processors request in mpiexec should be equal to nx*ny*nz, in the other way
   ! you will receive an error. 

    comm%dims(1) = 1
    comm%dims(2) = 1
    comm%dims(3) = 1 

    if(comm%dims(1).eq.0) comm%dims(1) = 1
    if(comm%dims(2).eq.0) comm%dims(2) = 1       
    if(comm%dims(3).eq.0) comm%dims(3) = 1

    isperiodic(1:3) = .false.
    reorder = .true.     ! allows the system reorganize the processors inside the topology 
    ndim = 3             ! in the most optimized way

    CALL MPI_INIT(comm%ierr)

    call MPI_CART_CREATE(MPI_COMM_WORLD, ndim, comm%dims, isperiodic, reorder, comm3d, comm%ierr)

    call MPI_COMM_SIZE(comm3d, comm%size, comm%ierr)

    call MPI_COMM_RANK(comm3d, comm%rank, comm%ierr)
    
    call MPI_CART_COORDS(comm3d, comm%rank, ndim, comm%coords, comm%ierr)

    displ(0) = -1
    displ(1) = 1
    n_dest = 0
    n_src = 0

    do direction=0, 2 
       do i=0,1
          ! obtaining the directions that the processor will send and receive data
          ! this information will be saved in shift(direcao,deslocamento)%src/dest
          ! when src/dest = -1, it means that the process will not receive/send data in the respective (dir,desl)
          call MPI_CART_SHIFT(comm3d, direction, displ(i), src_temp, dest_temp, comm%ierr) 

          if(src_temp>=0) then
            n_src = n_src + 1
            src(n_src) = src_temp
            src_dir(n_src) = (direction+1)*displ(i)  ! passing the directions from 0,1,2 to 1,2,3 and introducing the orientation
          elseif(dest_temp>=0) then
             n_dest = n_dest +1
             dest(n_dest) = dest_temp
             dest_dir(n_dest) = (direction+1)*displ(i)
          end if

       end do
    end do


    call domain_decomposition(v_local, L1,L2, V1, V2, border, comm)

    call ghost_update(v_local, border, src, dest, n_src, n_dest, src_dir, dest_dir,comm3d, comm)


    nstep = 0
    diff_total = 10
    tol = 0.00001
    converge = .true.

    do while (diff_total.ge.tol)
       
       diff = 0.d0
       nstep = nstep +1 
 
          !write(file_name,'(I4)') nstep
          !OPEN (UNIT=nstep+comm%rank,FILE='v'//trim(file_name)//'.xyz')
          !do i= border(1,0), border(1,1)
          !   do j= border(2,0), border(2,1)
                !do k= border(3,0), border(3,1)
          !         write(nstep+comm%rank,*) i,j, v_local(i,j,0)
                !end do
          !   end do
          !end do
          !close(nstep+comm%rank)
          
            
         call laplacian_calc(v_local,V1,V2, border, diff, L1, L2)

         call ghost_update(v_local, border, src, dest, n_src, n_dest, src_dir, dest_dir,comm3d, comm)

         call MPI_ALLREDUCE(diff, diff_total, 1, MPI_REAL8, MPI_SUM, comm3d, comm%ierr)

         write(*,*) diff_total
         if(abs(diff_total-diff_old).le.0.000000001) then
            converge = .false.
            EXIT
         end if
         diff_old = diff_total
         call MPI_BARRIER(comm3d,comm%ierr)

    end do

    ftime = mpi_wtime()
    if(converge.and.comm%rank.eq.0) then
       write(*,*) "Calculation converged in ", nstep," iterations."
       write(*,*) "Elapsed time:", ftime -itime
    elseif(.not.converge.and.comm%rank.eq.0) then
       write(*,*) "The calc doesn't converge. Number of iterations:", nstep
       write(*,*) "Elapsed time:", ftime -itime
    end if


    call MPI_FINALIZE(comm%ierr)

    DEALLOCATE(v_local)
  end subroutine run_laplace

end module laplace_m
