subroutine init(a, at, itot, jtot, ktot)
    integer, intent(in) :: itot, jtot, ktot
    real*8, dimension(itot,jtot,ktot), intent(inout) :: a, at
    integer :: ijk

    do k=1,ktot
        do j=1,jtot
            do i=1,itot
                ijk = i + (j-1)*itot + (k-1)*itot*jtot
                a(i,j,k)  = (ijk-1)**2. / ijk**2.     
                at(i,j,k) = 0
            end do
        end do
    end do
end subroutine init

subroutine diff(at, a, visc, dxidxi, dyidyi, dzidzi, itot, jtot, ktot)
    integer, intent(in) :: itot, jtot, ktot
    real*8, intent(in) :: visc, dxidxi, dyidyi, dzidzi
    real*8, dimension(itot,jtot,ktot), intent(inout) :: at
    real*8, dimension(itot,jtot,ktot), intent(in) :: a

    do k=2,ktot-1
        do j=2,jtot-1
            do i=2,itot-1
                at(i,j,k) = at(i,j,k) + visc * ( &
                        + ( (a(i+1,j,k) - a(i,j,k  )) &
                          - (a(i,j,k  ) - a(i-1,j,k)) ) * dxidxi &
                        + ( (a(i,j+1,k) - a(i,j,k  )) &
                          - (a(i,j,k  ) - a(i,j-1,k)) ) * dyidyi &
                        + ( (a(1,j,k+1) - a(i,j,k  )) &
                          - (a(i,j,k  ) - a(i,j,k-1)) ) * dzidzi &
                        )
            end do
        end do
    end do
end subroutine diff

program main
    integer :: nloop = 100
    integer :: itot  = 128
    integer :: jtot  = 128
    integer :: ktot  = 128

    real :: starttime, endtime
    real*8 :: tmp = 0.1

    real*8, dimension(:,:,:), allocatable :: a, at
    allocate(a(itot,jtot,ktot), at(itot,jtot,ktot))

    call init(a, at, itot, jtot, ktot)

    ! Check results
    call diff(at, a, tmp, tmp, tmp, tmp, itot, jtot, ktot)
    print*,"at=",at(itot/2+1, 2, 2)

    ! Time performance 
    call cpu_time(starttime)

    do i=1, nloop
        call diff(at, a, tmp, tmp, tmp, tmp, itot, jtot, ktot)
    end do

    call cpu_time(endtime)
    print*, "time/iter =", (endtime-starttime)/nloop, " s (", nloop, " iters)"
end program main
