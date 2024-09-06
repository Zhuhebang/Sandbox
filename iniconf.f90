SUBROUTINE iniconf(x0,y0,z0,r0,phi0)   
USE globalvars
IMPLICIT NONE

!整数类型，表示生成初始相场的中心坐标
INTEGER :: x0,y0,z0
!r0:双精度浮点数，表示生成初始相场的半径
!phi0(nx, ny, nz): 双精度浮点数的三维数组，表示相场变量，将被初始化
!rr(nx, ny, nz): 双精度浮点数的三维数组，存储每个网格点到中心点的距离
REAL*8 :: r0,phi0(nx,ny,nz),rr(nx,ny,nz)
!ii, jj, kk: 整数类型的循环变量，用于遍历网格点
INTEGER :: ii,jj,kk

IF(ndim == 2) THEN
!PRINT*,'Two dimensionsional simulation'
DO ii=1,nx
 DO jj=1,ny
  DO kk=1,nz
    !计算每个网格点到中心 (x0, y0) 的欧氏距离，DBLE 函数将整数转换为双精度浮点数，以确保计算的准确性
   rr(ii,jj,kk) = SQRT(DBLE(ii - x0)**2 + DBLE(jj - y0)**2)
    !使用双曲正切函数 TANH 来生成平滑的相场过渡
    !TANH 函数用于模拟相场在半径 r0 处的过渡，创造一个从1到0的平滑边界，4.0 * (rr(ii, jj, kk) - r0) 控制过渡的宽度
   phi0(ii,jj,kk) = 1.0 - (0.5*TANH(4.0*(rr(ii,jj,kk) - r0)) + 0.5)
  ENDDO
 ENDDO
ENDDO

!Generate initial configuration (3D sphere)
!初始化代码尚未实现
ELSEIF (ndim == 3) THEN        
PRINT*,'Three dimensionsional simulation'
DO ii=1,nx
 DO jj=1,ny
  DO kk=1,nz
  ENDDO
 ENDDO
ENDDO
ENDIF

END SUBROUTINE iniconf
