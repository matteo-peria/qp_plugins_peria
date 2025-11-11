
 BEGIN_PROVIDER [double precision, grid3_lebedev_points, (grid3_ang_size,3) ]
&BEGIN_PROVIDER [double precision, grid3_lebedev_weights, (grid3_ang_size)]
  implicit none
  BEGIN_DOC
  ! weights and grid points for the integration on the angular variables on
  ! the unit sphere centered on (0,0,0)
  ! According to the LEBEDEV scheme
  END_DOC

  include 'constants.include.F'
  integer                        :: i
  double precision               :: accu
  double precision               :: degre_rad
  double precision               :: x(grid3_ang_size)
  double precision               :: y(grid3_ang_size)
  double precision               :: z(grid3_ang_size)
  double precision               :: w(grid3_ang_size)

  degre_rad = pi/180.d0
  accu = 0.d0

  select case (grid3_ang_size)


  case (0006)
      call LD0006(X,Y,Z,W,grid3_ang_size)
  case (0014)
      call LD0014(X,Y,Z,W,grid3_ang_size)
  case (0026)
      call LD0026(X,Y,Z,W,grid3_ang_size)
  case (0038)
      call LD0038(X,Y,Z,W,grid3_ang_size)
  case (0050)
      call LD0050(X,Y,Z,W,grid3_ang_size)
  case (0074)
      call LD0074(X,Y,Z,W,grid3_ang_size)
  case (0086)
      call LD0086(X,Y,Z,W,grid3_ang_size)
  case (0110)
      call LD0110(X,Y,Z,W,grid3_ang_size)
  case (0146)
      call LD0146(X,Y,Z,W,grid3_ang_size)
  case (0170)
      call LD0170(X,Y,Z,W,grid3_ang_size)
  case (0194)
      call LD0194(X,Y,Z,W,grid3_ang_size)
  case (0230)
      call LD0230(X,Y,Z,W,grid3_ang_size)
  case (0266)
      call LD0266(X,Y,Z,W,grid3_ang_size)
  case (0302)
      call LD0302(X,Y,Z,W,grid3_ang_size)
  case (0350)
      call LD0350(X,Y,Z,W,grid3_ang_size)
  case (0434)
      call LD0434(X,Y,Z,W,grid3_ang_size)
  case (0590)
      call LD0590(X,Y,Z,W,grid3_ang_size)
  case (0770)
      call LD0770(X,Y,Z,W,grid3_ang_size)
  case (0974)
      call LD0974(X,Y,Z,W,grid3_ang_size)
  case (1202)
      call LD1202(X,Y,Z,W,grid3_ang_size)
  case (1454)
      call LD1454(X,Y,Z,W,grid3_ang_size)
  case (1730)
      call LD1730(X,Y,Z,W,grid3_ang_size)
  case (2030)
      call LD2030(X,Y,Z,W,grid3_ang_size)
  case (2354)
      call LD2354(X,Y,Z,W,grid3_ang_size)
  case (2702)
      call LD2702(X,Y,Z,W,grid3_ang_size)
  case (3074)
      call LD3074(X,Y,Z,W,grid3_ang_size)
  case (3470)
      call LD3470(X,Y,Z,W,grid3_ang_size)
  case (3890)
      call LD3890(X,Y,Z,W,grid3_ang_size)
  case (4334)
      call LD4334(X,Y,Z,W,grid3_ang_size)
  case (4802)
      call LD4802(X,Y,Z,W,grid3_ang_size)
  case (5294)
      call LD5294(X,Y,Z,W,grid3_ang_size)
  case (5810)
      call LD5810(X,Y,Z,W,grid3_ang_size)
      case default
      print *, irp_here//': wrong grid3_ang_size. See in ${QP_ROOT}/src/becke_numerical_grid/list_angular_grid to see the possible angular grid points. Ex: '
      print *, '[ 50 | 74 | 170 | 194 | 266 | 302 | 590 | 1202 | 2030 | 5810 ]'
      stop -1
  end select

  do i = 1, grid3_ang_size
    grid3_lebedev_points(i,1) = x(i)
    grid3_lebedev_points(i,2) = y(i)
    grid3_lebedev_points(i,3) = z(i)
    grid3_lebedev_weights(i) = w(i) * 4.d0 * pi
    accu += w(i)
  enddo

END_PROVIDER
