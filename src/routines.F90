module routines

  use mpi

  implicit none

  private

  public :: laplacian_calc,&
            ghost_update,&
            domain_decomposition,&
            heaviside

  type comm_t
     integer :: size, &
                rank, &
                ierr, &
                coords(3), &
                dims(3)
  end type comm_t

  contains

    subroutine laplacian_calc(v_local,V1,V2, border, diff, L1, L2)

      implicit none

      real(8), allocatable, intent(inout) :: v_local(:,:,:)
      real(8), intent(in) :: V1, V2
      real(8), allocatable :: v_new(:,:,:)
      integer, intent(in) :: L1, L2, border(1:3,0:1)
      real(8), intent(out) :: diff
      real(8) :: a
      integer :: i, j, k


      ALLOCATE(v_new(border(1,0):border(1,1), border(2,0):border(2,1), border(3,0):border(3,1) ) )
      
      a = 0.166666666667
      diff = 0.d0
    
      do i= border(1,0), border(1,1)
          do j= border(2,0), border(2,1)
              do k= border(3,0), border(3,1)
                  

               if(abs(i)>L1 .or. abs(j)> L1 .or. abs(k)> L1) then
                    v_new(i,j,k) = a*( v_local(i-1,j,k) + v_local(i+1,j,k) + &
                                       v_local(i,j-1,k) + v_local(i,j+1,k) + &
                                       v_local(i,j,k-1) + v_local(i,j,k+1) )

               else if(abs(i).le.L1 .or. abs(j).le. L1 .or. abs(k) .le. L1) then
                  ! assuring the boundary conditions
                  v_new(i,j,k) = V1
            
               end if

               ! assuring the boundary conditions
               if(i.le. -L2  ) then
                  v_new(i,j,k) = V2
               elseif(i .ge. L2 -1 ) then
                  v_new(i,j,k) = V2
               elseif(j .le.- L2 ) then
                  v_new(i,j,k) = V2
               elseif(j .ge. L2 -1 ) then
                  v_new(i,j,k) = V2
               elseif(k .le. -L2 ) then
                  v_new(i,j,k) = V2
               elseif(k .ge. L2 -1 ) then
                  v_new(i,j,k) = V2
               end if


              end do
           end do
        end do
        
        diff = 0.d0
        do i= border(1,0), border(1,1)
           do j= border(2,0), border(2,1)
              do k= border(3,0), border(3,1)
                 diff = diff +  abs((v_new(i,j,k) - v_local(i,j,k) ))/abs(v_local(i,j,k))
              end do
           end do
        end do

        v_local(border(1,0):border(1,1), border(2,0):border(2,1), border(3,0):border(3,1)) = v_new(border(1,0):border(1,1), border(2,0):border(2,1), border(3,0):border(3,1))
        
        DEALLOCATE(v_new)
    end subroutine laplacian_calc

    subroutine ghost_update(v_local, border, src, dest, n_src, n_dest, src_dir, dest_dir, comm3d, comm)

      ! this subroutine update the ghost points for each process
      ! using MPI_SENDRECV function 
     
      implicit none

      integer, intent(in) :: src(6), dest(6), n_src, n_dest, src_dir(6), dest_dir(6), comm3d
      integer, intent(in) :: border(1:3,0:1)
      real(8), allocatable, intent(inout) :: v_local(:,:,:)
      integer :: i, border_send, border_recv, send_tag, recv_tag, sendcount, recvcount
      type(comm_t) :: comm


      do i=1, n_dest

         border_recv =  border(abs(src_dir(i)),heaviside(-src_dir(i)))
         border_send =  border(abs(dest_dir(i)),heaviside(dest_dir(i)))
         send_tag = i
         recv_tag = i

         if(abs(dest_dir(i)).eq.1) then
            sendcount = size(v_local(border_send,:,:))
            recvcount = size(v_local(border_recv,:,:))
    
            call MPI_SENDRECV(v_local(border_send,:,:), sendcount, MPI_REAL8, dest(i), send_tag,&
                 v_local(border_recv-sign(1,src_dir(i)),:,:), recvcount, MPI_REAL8, src(i), recv_tag, comm3d,MPI_STATUS_IGNORE, comm%ierr )

         elseif(abs(dest_dir(i)).eq.2) then
          
            sendcount = size(v_local(:,border_send,:))
            recvcount = size(v_local(:,border_recv,:)) 

            call MPI_SENDRECV(v_local(:,border_send,:), sendcount, MPI_REAL8, dest(i), send_tag,&
                 v_local(:,border_recv-sign(1,src_dir(i)),:), recvcount, MPI_REAL8, src(i), recv_tag, comm3d,MPI_STATUS_IGNORE, comm%ierr )

         elseif(abs(dest_dir(i)).eq.3) then
        
            sendcount = size(v_local(:,:,border_send)) 
            recvcount = size(v_local(:,:,border_recv)) 

            call MPI_SENDRECV(v_local(:,:,border_send), sendcount, MPI_REAL8, dest(i), send_tag,&
                 v_local(:,:,border_recv-sign(1,src_dir(i))), recvcount, MPI_REAL8, src(i), recv_tag, comm3d,MPI_STATUS_IGNORE, comm%ierr )
         end if

      end do

    end subroutine ghost_update




    subroutine domain_decomposition(v_local, L1,L2, V1, V2, border, comm)

      ! this subroutine generate the system initial condition in each process
      ! and save the border information for later update the ghost points

      implicit none
      
      real(8), allocatable, intent(inout) :: v_local(:,:,:)
      type(comm_t), intent(in) :: comm
      real, intent(in) :: V1, V2
      integer, intent(in) :: L1, L2
      integer, intent(out) :: border(3,0:1)
      integer :: ii, is, ji, js, ki, ks, i, j, k
      integer :: gp, iseed
      

      gp = 1
      iseed = -129388384
   
         
      ii = 2*L2*(comm%coords(1))/comm%dims(1) - L2
      is = 2*L2*(comm%coords(1)+1)/comm%dims(1)- L2 -1
      
      ji = 2*L2*comm%coords(2)/comm%dims(2) - L2
      js = 2*L2*(comm%coords(2)+1)/comm%dims(2)- L2 -1
      
      ki = 2*L2*comm%coords(3)/comm%dims(3) - L2
      ks = 2*L2*(comm%coords(3)+1)/comm%dims(3) - L2 -1

      border(1,0) = ii
      border(1,1) = is
      border(2,0) = ji
      border(2,1) = js
      border(3,0) = ki
      border(3,1) = ks


      ALLOCATE(v_local(ii-gp:is+gp,ji-gp:js+gp,ki-gp:ks+gp))

      do i=ii, is
         do j=ji, js
            do k=ki,ks


               if(abs(i)>L1 .or. abs(j)> L1 .or. abs(k)> L1) then
                  v_local(i,j,k) =  9.0 !(ran2(iseed)+1)*V1
               else if(abs(i).le.L1 .or. abs(j).le. L1 .or. abs(k) .le. L1) then
                  v_local(i,j,k) = V1
               end if

               if(abs(i).eq. L2 - heaviside(i) .or. abs(j) .eq. L2  - heaviside(j) .or. abs(k) .eq. L2 - heaviside(k)) then
                    v_local(i,j,k) = V2
               end if

            end do
         end do
      end do


    end subroutine domain_decomposition

    function heaviside(x)

      ! the step function is important, because our system have their center in (0,0,0)
      ! which means that our range is from -L/2 until L/2-1, so you can write the boundary
      ! condition in a compact way using the heaviside function as: abs(x).eq.L/2 -heaviside(x)

      implicit none
      
      integer :: x, heaviside
      
      if ( x<0) then
         heaviside = 0
      else
         heaviside = 1
      end if
    end function heaviside


  end module routines
