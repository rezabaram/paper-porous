program critpore_np
  !     this programm find the characteristic pore lenght lc

interface
   function distance(lam, x, y, z, ellps)
     real, intent(in) :: x, y, z, lam
     real, intent(in), dimension(10) :: ellps
   end function distance
   
   function shell_distance( x, y, z, ellps)
    real, intent(in) :: x, y, z
    real, intent(in), dimension(10) :: ellps
   end function shell_distance

   function ellip_func(x, y, z, ellps) 
     real, intent(in) :: x, y, z
     real, intent(in), dimension(10) :: ellps
   end function ellip_func

   function toBody(x, y, z, ellps) result (xrt)
     real, intent(in) :: x, y, z
     real, intent(in), dimension(10) :: ellps
     real, dimension(3) :: xrt
   end function toBody

   function slambda(x, y, z, ellps)
     real, intent(in) :: x, y, z
     real, intent(in), dimension(10) :: ellps
   end function slambda
end interface

  !     define size of porous media
  integer :: n, ngb(26,3), i, j, k, q, ii, jj, kk, qq, t
  real :: lc, lcc, lc2, dist, r, rr, lam
  logical :: pp, plc
  real, dimension(:,:,:), allocatable :: pore
  integer, dimension(:,:,:,:), allocatable :: cells
  integer, dimension(:,:,:), allocatable :: last
  integer :: i1, i2, j1, j2, k1, k2, kn, kt
  integer :: ncells ! number of ellipsoides in subdomain
  real :: dcell
  logical, dimension(:,:,:), allocatable :: paint
  real, dimension(:,:), allocatable :: ell
  integer :: l, n_ell ! number of ellipsoides\

  real, dimension(3) :: p, xrt
  real :: dlc, lc_test(10)
  logical :: plc_test(10)
 
  n = 512    ! domain size

  ! allocate domain size
  allocate(pore(n,n,n))
  allocate(paint(n,n,n))

  !   initializing neighborhood
  ngb(1,1)  = -1; ngb(1,2)  = -1; ngb(1,3)  = -1
  ngb(2,1)  =  0; ngb(2,2)  = -1; ngb(2,3)  = -1
  ngb(3,1)  =  1; ngb(3,2)  = -1; ngb(3,3)  = -1
  ngb(4,1)  = -1; ngb(4,2)  =  0; ngb(4,3)  = -1
  ngb(5,1)  =  0; ngb(5,2)  =  0; ngb(5,3)  = -1
  ngb(6,1)  =  1; ngb(6,2)  =  0; ngb(6,3)  = -1
  ngb(7,1)  = -1; ngb(7,2)  =  1; ngb(7,3)  = -1
  ngb(8,1)  =  0; ngb(8,2)  =  1; ngb(8,3)  = -1
  ngb(9,1)  =  1; ngb(9,2)  =  1; ngb(9,3)  = -1
  ngb(10,1) = -1; ngb(10,2) = -1; ngb(10,3) =  0 
  ngb(11,1) =  0; ngb(11,2) = -1; ngb(11,3) =  0 
  ngb(12,1) =  1; ngb(12,2) = -1; ngb(12,3) =  0
  ngb(13,1) = -1; ngb(13,2) =  0; ngb(13,3) =  0 

  ngb(14,1) =  1; ngb(14,2) =  0; ngb(14,3) =  0 
  ngb(15,1) = -1; ngb(15,2) =  1; ngb(15,3) =  0 
  ngb(16,1) =  0; ngb(16,2) =  1; ngb(16,3) =  0 
  ngb(17,1) =  1; ngb(17,2) =  1; ngb(17,3) =  0
  ngb(18,1) = -1; ngb(18,2) = -1; ngb(18,3) =  1
  ngb(19,1) =  0; ngb(19,2) = -1; ngb(19,3) =  1
  ngb(20,1) =  1; ngb(20,2) = -1; ngb(20,3) =  1
  ngb(21,1) = -1; ngb(21,2) =  0; ngb(21,3) =  1
  ngb(22,1) =  0; ngb(22,2) =  0; ngb(22,3) =  1
  ngb(23,1) =  1; ngb(23,2) =  0; ngb(23,3) =  1
  ngb(24,1) = -1; ngb(24,2) =  1; ngb(24,3) =  1
  ngb(25,1) =  0; ngb(25,2) =  1; ngb(25,3) =  1
  ngb(26,1) =  1; ngb(26,2) =  1; ngb(26,3) =  1


  ! reading input-file
  open (unit=99, file='#1.dat', status='old', action='read')
  read(99, *), n_ell
  allocate(ell(10,n_ell))
  read(99,*) ell
  close(unit=99)
  ! eta1.4_xi1.8.dat

  ! cut off ratio 
  r = 0.0
  do l = 1,n_ell
     r = max( max( max( ell(4,l), ell(5,l)), ell(6,l ) ), r )
  end do
  
  ! subdomain descomposition
  ncells = floor( 1/r )
  allocate(cells(ncells,ncells,ncells,100))
  allocate(last(ncells,ncells,ncells))
  last(:,:,:) = 0
  cells(:,:,:,:) = 0
  dcell = 1.0/(1.0*ncells)
 
  do l = 1,n_ell
     rr = max ( ell(4,l), max (ell(5,l), ell(6,l ) ) ) 
     
     i1 = floor((ell(1,l) - rr)/dcell) + 1
     if (i1 .lt. 1) i1 = 1
     i2 = floor((ell(1,l) + rr)/dcell) + 1
     if (i2 .gt. ncells) i2 = ncells
     j1 = floor((ell(2,l) - rr)/dcell) + 1
     if (j1 .lt. 1) j1 = 1
     j2 = floor((ell(2,l) + rr)/dcell) + 1
     if (j2.gt.ncells) j2 = ncells
     k1 = floor((ell(3,l) - rr)/dcell) + 1
     if (k1 .lt. 1) k1 = 1
     k2 = floor((ell(3,l) + rr)/dcell) + 1
     if (k2 .gt. ncells) k2 = ncells
     
     do i = i1,i2
        do j = j1,j2
           do k = k1,k2
              last(i,j,k) = last(i,j,k) + 1
              cells(i,j,k,last(i,j,k)) = l
           end do
        end do
     end do
  end do
  print *, 'subdomain descomposition created'
  print *, '# cells = ', ncells 
  print *, 'dcell = ', dcell
    

  open (unit = 17, file = "#1.txt")   
  write(17,*) n
  
  !     sphere map
  do i = 1,n
     do j = 1,n
        do k = 1,n 
           paint(i,j,k) = .false.
           pore(i,j,k) = 2.0*r
           
           ! cell coordinates
           ii = floor((i - 0.5)/(n*dcell)) + 1
           jj = floor((j - 0.5)/(n*dcell)) + 1
           kk = floor((k - 0.5)/(n*dcell)) + 1
           do l = 1,last(ii,jj,kk)
              ! translated-rotated coordinate
              xrt = toBody( (i - 0.5)/n, (j - 0.5)/n, (k - 0.5)/n, ell(:, cells(ii,jj,kk,l) ) )
              
              if ( ellip_func( xrt(1), xrt(2), xrt(3), ell(:, cells(ii,jj,kk,l) ) ) .le. 0  ) then 
                 pore(i,j,k) = 0.0
              else 
                 lam = slambda( xrt(1),xrt(2),xrt(3), ell(:, cells(ii,jj,kk,l) ) )
                 dist = distance(lam, xrt(1), xrt(2), xrt(3), ell(:, cells(ii,jj,kk,l) ) )
                 ! print *, i, j, k, l, (shell_distance( xrt(1), xrt(2), xrt(3), ell(:, cells(ii,jj,kk,l) ) ) - dist)/dist, dist
                 pore(i,j,k) = min(pore(i,j,k), dist)             
              end if
           end do
           ! write(17,*) i, j, k, pore(i,j,k)
         end do
     end do
     print *, i
  end do 
  print *, 'sphere map created'



  lc = 2.0*r
  lcc = 0.0d0
  lc2 = 2.0*r
  ! starting lc 
  do i = 1,n        
     do k = 1,n
        do j = 1,n           
           if ( pore(i,j,k) .lt. lc ) then 
              lcc = max( pore(i,j,k), lcc )
           end if
        end do
     end do
     lc2 = min(lcc, lc2)
  end do
  lc = lc2  

  dlc = lc/10.0
  kn = 1

  l = 1
  do while ( l .le. 10 )
     lc_test(l) = lc - ( l - 1 )*dlc

     ! starting from i=1
     do i = 1,n
        do k = 1,n
           if ( pore(i, 1, k) .ge. lc_test(l) ) then
              paint(i, 1, k) = .true.
           end if
        end do
     end do
     
     ! painting the path
     pp = .true.

     do while ( pp .and. kn .lt. n )
        pp = .false.
        do j = 1,kn
           do i = 1,n
              do k = 1,n
                 if ( paint(i, j, k) ) then ! only if i,j,k is true
                    do t = 1,26
                       ! calculate neighbor
                       ii = i + ngb(t, 1)
                       if (ii .lt. 1) ii = 1
                       if (ii .gt. n) ii = n
                       
                       jj = j + ngb(t, 2)
                       if (jj .lt. 1) jj = 1
                       if (jj .gt. n) jj = n
                       
                       kk = k + ngb(t, 3)
                       if (kk .lt. 1) kk = 1
                       if (kk .gt. n) kk = n
                       
                       ! checking
                       if ( (.not. paint(ii, jj, kk)) .and. pore(ii, jj, kk) .ge. lc_test(l) ) then ! this excludes the case i,j,k = ii,jj,kk
                          paint(ii, jj, kk) = .true.
                          pp = .true.

                          kn = max(kn, jj)
                       end if
                    end do
                 end if
              end do
           end do
        end do
     end do
     
     plc_test(l) = .false.
     if ( kn .eq. n ) plc_test(l) = .true.

     print *, l, lc_test(l), plc_test(l)
     l = l + 1
  end do


  do l = 1,10    
     if ( .not. plc_test(l) ) lc = lc_test(l)
  end do

  paint(:,:,:) = .false.

  plc = .false.
  kn = 1
  
  do while ( .not. plc )
     lcc = 0.0d0
     lc2 = 2.0*r
     
     ! starting lc 
     do j = 1,n        
        do k = 1,n
           do i = 1,n           
              if ( pore(i,j,k) .lt. lc ) then 
                 lcc = max( pore(i,j,k), lcc )
              end if
           end do
        end do
        lc2 = min(lcc, lc2)
     end do
     lc = lc2  
     print *, lc

     ! starting from i=1
     do i = 1,n
        do k = 1,n
           if ( pore(i, 1, k) .ge. lc ) then
              paint(i, 1, k) = .true.
           end if
        end do
     end do
     
     ! painting the path
     pp = .true.
     
     k1 = 1
     do while (pp)
        pp = .false.
        kt = n
        do j = k1,kn
           do i = 1,n
              do k = 1,n
                 if ( paint(i, j, k) ) then ! only if i,j,k is true
                    do t = 1,26
                       ! calculate neighbor
                       ii = i + ngb(t, 1)
                       if (ii .lt. 1) ii = 1
                       if (ii .gt. n) ii = n
                       
                       jj = j + ngb(t, 2)
                       if (jj .lt. 1) jj = 1
                       if (jj .gt. n) jj = n
                       
                       kk = k + ngb(t, 3)
                       if (kk .lt. 1) kk = 1
                       if (kk .gt. n) kk = n
                       
                       ! checking
                       if ( (.not. paint(ii, jj, kk)) .and. pore(ii, jj, kk) .ge. lc ) then ! this excludes the case i,j,k = ii,jj,kk
                          paint(ii, jj, kk) = .true.
                          pp = .true.
                          kt = min(kt, jj)
                          kn = max(kn, jj)
                       end if
                    end do
                 end if
              end do
           end do
        end do
        k1 = kt
     end do
     
     ! checking z = n
     plc = .false.
     if ( kn .eq. n ) plc = .true.
     
  end do
  
  write(17,*) 'lc =', lc
  close(17)

 end program critpore_np
 

! ellipsoide function
! x y z must be written and translated-rotated coordinate
 function ellip_func(x, y, z, ellps) 
   real :: ellip_func
   real, intent(in) :: x, y, z
   real, intent(in), dimension(10) :: ellps
   real :: a2, b2, c2

   a2 = ellps(4)**2
   b2 = ellps(5)**2
   c2 = ellps(6)**2
   
   ellip_func = x**2/a2 + y**2/b2 + z**2/c2 - 1.0 
 end function ellip_func

 
! translation and ratation of te coordinate system
 function toBody (x, y, z, ellps) result (xrt)
   real, intent(in) :: x, y, z
   real, intent(in), dimension(10) :: ellps
   real, dimension(3) :: xt, xrt
   real :: wp, xp, yp ,zp, wpp
   
   ! translation
   xt(1) = x - ellps(1)      
   xt(2) = y - ellps(2)      
   xt(3) = z - ellps(3)   

   ! rotation
   wp =   ( - xt(1)*ellps(8)  - xt(2)*ellps(9)  - xt(3)*ellps(10) ) 
   xp = - (   xt(1)*ellps(7)  + xt(2)*ellps(10) - xt(3)*ellps(9)  ) 
   yp = - ( - xt(1)*ellps(10) + xt(2)*ellps(7)  + xt(3)*ellps(8)  ) 
   zp = - (   xt(1)*ellps(9)  - xt(2)*ellps(8)  + xt(3)*ellps(7)  ) 

   wpp   =   ( wp*ellps(7)  - xp*ellps(8)  - yp*ellps(9)  - zp*ellps(10) ) 
   xrt(1)= - ( wp*ellps(8)  + xp*ellps(7)  + yp*ellps(10) - zp*ellps(9)  ) 
   xrt(2)= - ( wp*ellps(9)  - xp*ellps(10) + yp*ellps(7)  + zp*ellps(8)  ) 
   xrt(3)= - ( wp*ellps(10) + xp*ellps(9)  - yp*ellps(8)  + zp*ellps(7)  ) 
 end function toBody


! solve the value of lambda
! x y z must be written and translated-rotated coordinate
 function slambda(x, y, z, ellps)
   real :: slambda
   real, intent(in) :: x, y, z
   real, intent(in), dimension(10) :: ellps
   real :: a2, b2, c2, x2, y2, z2, dlambda, error, minlambda
   integer :: t
   
   a2 = ellps(4)**2
   b2 = ellps(5)**2
   c2 = ellps(6)**2
   
   x2 = x**2
   y2 = y**2
   z2 = z**2
   
   t = 0
   error = 1.0
   slambda = 0.0 

   ! Newton–Raphson method x = x - f(x)/f'(x)
   do while( error .gt. 0.01 .and. t .le. 10)
      t = t + 1
      dlambda = - ( -1.0 + a2*x2/(a2 + slambda)**2 + b2*y2/(b2 + slambda)**2 + c2*z2/(c2 + slambda)**2 ) / &
           ( -2.0*a2*x2/(a2 + slambda)**3 - 2.0*b2*y2/(b2 + slambda)**3 - 2.0*c2*z2/(c2 + slambda)**3 )
      slambda = slambda + dlambda
      error = abs(dlambda/slambda)
    end do
 
 end function slambda
 

! short distance between a point (x,y,z) and a ellipsoide
! x y z must be written and translated-rotated coordinate
 function distance(lam, x, y, z, ellps)
   real :: distance
   real, intent(in) :: x, y, z, lam
   real, intent(in), dimension(10) :: ellps
   real :: a2, b2, c2, px, py, pz

   a2 = ellps(4)**2
   b2 = ellps(5)**2
   c2 = ellps(6)**2

   px = a2*x/(a2 + lam)
   py = b2*y/(b2 + lam)
   pz = c2*z/(c2 + lam)

   distance = sqrt( (px - x)**2 + (py - y)**2 + (pz - z)**2 )
 end function distance
 

 function shell_distance( x, y, z, ellps)
   real :: shell_distance
   real, intent(in) :: x, y, z
   real, intent(in), dimension(10) :: ellps
   real :: a, b, c

   a = ellps(4)
   b = ellps(5)
   c = ellps(6)
   r = max( max( a, b), c )

   !the center of ellip is at the origin
   shell_distance = sqrt( x**2 + y**2 + z**2 ) - r
end function shell_distance

