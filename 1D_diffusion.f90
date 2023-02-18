!-----------------------------------------------------------------------------------------
! one-dimensional diffusuion program （１次元移流拡散プログラム）         by Ichiro KIMURA
!-----------------------------------------------------------------------------------------
!=========================================================================================
MODULE     parameters_m
  !======= output dat name ===============================================================
  character*100,parameter :: cwrite='result'                               !出力ファイル名
END MODULE parameters_m
!=========================================================================================
!=========================================================================================
MODULE VARIABLES_m
  implicit none
  real*8, dimension(:),allocatable:: XG, UG, AG, DG, CC, CG
  integer :: nx, i_observe
  real*8  :: clen, cl_in, cl_observe, q0, a0, e0, d0
  real*8  :: t_total, t_file, dt
  real*8  :: dx
END MODULE VARIABLES_m
!=========================================================================================
module input_para_m
  use variables_m
  implicit none
  contains
subroutine input_para
  integer :: i
  real*8  :: u0
  !
  open(10, file='input.csv', status='old')
    read(10,*) clen                !水路長 (m)
    read(10,*) cl_in               !物質投入位置(m)
    read(10,*) cl_observe          !観測地点位置(m)
    read(10,*) a0                  !水路断面積(m2)
    read(10,*) q0                  !流量(m3/s)
    read(10,*) e0                  !投入物質量(kg)
    read(10,*) d0                  !拡散係数(m2/s)
    read(10,*) t_total, t_file, dt !全計算時間, ファイル出力間隔, 計算時間ステップ (s)
    read(10,*) nx                  !(x,y)方向の格子セル数
  close(10)
  !
  dx  = clen / dble(nx)            !計算格子のx方向幅
  u0  = q0   / a0                  !流速
  !
  allocate( XG(0:nx), UG(0:nx), AG(0:nx), DG(0:nx), CC(1:nx) )
  !
  AG(0:nx) = a0
  UG(0:nx) = u0
  DG(0:nx) = d0
  !
  do i = 0, nx
    XG(i) = dble(i)*dx
  end do
  !
end subroutine input_para
end module     input_para_m
!=========================================================================================
!=========================================================================================
SUBROUTINE TIME_STEP( t_total,t_file,dt,kend,kfile )                 !時間ステップ数の計算
  implicit none
  integer,intent(out) :: kend, kfile
  real*8, intent(in ) :: t_total, t_file, dt
      KEND  = INT( t_total / DT + 0.5d0 )                  !最終時間ステップ
      KFILE = INT( t_file  / DT + 0.5d0 )                  !ファイル出力間隔ステップ
END SUBROUTINE TIME_STEP
!=========================================================================================
!=========================================================================================
module INIT_m                                                                  !変数初期化
  use parameters_m ; use variables_m
  implicit none
  contains
subroutine INIT
  integer :: i_in, i
  !
  i_in      = int( cl_in      / dx ) + 1
  i_observe = int( cl_observe / dx ) + 1
  CC( 1:nx) = 0.d0
  CC( i_in) = e0 / ( (AG(i_in)+AG(i_in-1))*0.5d0 * dx )
  !
  !! write(*,*) i_in, e0, cc(i_in) ; pause
  !
  open(20, file = trim(cwrite)//'.txt', status = 'unknown' )
  write(20,'(A22)'   ) '#     t(min)         c'
  write(20,'(F15.10,E12.5)') 0.d0, CC(i_observe)/1000.d0*(1.d+6)*(1.d+2)
end subroutine INIT
end module     INIT_m
!=========================================================================================
!=========================================================================================
module ADVEC_m                                                            !運動方程式の計算
  use parameters_m ; use variables_m
  implicit none
  contains
subroutine ADVEC
  real*8  :: ap, ai, am, adve_d, adve_u, diff_d, diff_u, C1(1:nx)
  integer :: ip, im, imm,i
  !
  do i = 1, nx
    ip  = min(i+1,nx)
    im  = max(i-1, 1)
    imm = max(i-2, 0)
    ap  = ( AG(i  ) + AG(ip ) )*0.5d0
    ai  = ( AG(i  ) + AG(i-1) )*0.5d0
    am  = ( AG(i-1) + AG(imm) )*0.5d0
    adve_d = ( (UG(i  )+abs(UG(i  )))*CC(i )*ai + (UG(i  )-abs(UG(i  )))*CC(ip)*ap )*0.5d0
    adve_u = ( (UG(i-1)+abs(UG(i-1)))*CC(im)*am + (UG(i-1)-abs(UG(i-1)))*CC(i )*ai )*0.5d0
    diff_d = ( CC(ip)-CC(i ) )/dx * AG(i  ) * DG(i  )
    diff_u = ( CC(i )-CC(im) )/dx * AG(i-1) * DG(i-1)
    C1(i) = CC(i) + dt * ( - (adve_d-adve_u)/dx + (diff_d-diff_u)/dx ) / ai
  !!C1(i) = CC(i) + dt * ( - (adve_d-adve_u)/dx )                      / ai
  !  if( c1(i)/=0.d0) then
  !    write(*,*) i,cc(i), c1(i), dt, dx, adve_d, adve_u, ug(i) ;pause
  !  endif
  end do
  !
  CC(1:nx) = C1(1:nx)
end subroutine ADVEC
end module     ADVEC_m
!=========================================================================================
!=========================================================================================
subroutine write_dat( kf, cwrite, nx, x, c )                              !結果のファイル出力
  implicit none
    integer,intent(in) :: kf, nx
    integer            :: ier, i, n
    real*8,dimension(1:nx),intent(in) :: c
    real*8,dimension(0:nx),intent(in) :: x
    character*100, intent(in) :: cwrite
    character                 :: header1*52, header2*46
    character*4               :: cnum
    real*8                    :: c_out
  write( cnum, '(i4)' ) kf
  do n = 1, 4
    if( cnum(n:n)==' ') cnum(n:n) = '0'
  end do
  !
  header1 ='#     x(km)          c'
  open( 10, file = trim(cwrite)//'.'//cnum, status = 'unknown' )
  write(10,'(A22)') trim(header1)
  do i = 1, nx
    c_out = c(i)
    if ( abs(c_out) < 1d-18 ) c_out = 0.d0
    write(10,'(F15.10, E12.5)') (x(i)+x(i-1))*0.5d0/1000.d0, c_out/1000.d0*(1.d+6)*(1.d+2)
  end do
  close(10)
end subroutine write_dat
!=========================================================================================
!=========================================================================================
PROGRAM MAIN                                                             !メインプログラム
  use parameters_m ; use variables_m
  use input_para_m
  use INIT_m
  use ADVEC_m
  implicit none
    integer :: KEND, KFILE, K
    real*8  :: time, c_out
  !==== read initial parameters =========================!input.csvから計算条件読み込み
  call input_para
  !==== preparations for computations ===================!計算の準備
  call INIT                                                !変数の初期化
  call TIME_STEP( T_TOTAL, T_FILE, dt, KEND, KFILE )       !全step数とfile出力step間隔
  !---- initial output ---------------------------------
  time = 0.d0                                              !時刻初期化
  call write_dat( 0,cwrite,nx,XG(0:nx),CC(1:nx) )          !初期値を出力
  !==== main loop start =================================!メインのループがスタート
  do K = 1, KEND                                           !main loop のdo文
    time = dble(K) * dt                                    !現在時刻の計算
    IF(mod(K,100)==0) write(*,*) 'K=',K,' T= ',time        !現在時刻とタイムステップ数表示
    call ADVEC                                             !移流拡散方程式の計算
    IF( mod(K,KFILE)==0 ) then                             !ファイル出力有無判定
      call write_dat(k/kfile,cwrite,nx,XG(0:nx),CC(1:nx))  !計算結果出力
      c_out = cc(i_observe)
      if( abs(c_out) < 1.d-18 ) c_out = 0.d0
      write(20,'(F15.10, E12.5)') time/60.d0, c_out/1000.d0*(1.d+6)*(1.d+2)
    endif
  end do
  close(20)
END PROGRAM MAIN
!=========================================================================================
