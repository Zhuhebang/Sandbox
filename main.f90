PROGRAM main  
USE globalvars
USE invert
IMPLICIT NONE
INCLUDE 'fftw3.f' 

!---------------------------------------------------------------------------------------------------------------------
! 定义变量部分
!p: 表示 particle，即填料或颗粒。这些通常是嵌入在基体材料中的纳米颗粒或其它形式的填料，用来改变材料的整体特性
!m: 表示 matrix，即基体。基体是材料的连续相，填料或其他组分分散在其中
!s: 表示 shell ，在纳米复合材料中，shell 可以指颗粒的外壳或界面层，这一层可能有不同的材料特性（如电导率、介电常数等）
!v: 表示 void ?，即孔洞或空隙，这些空隙会影响材料的整体特性，如电导率、介电常数和机械性能等

!用于循环计数和控制流
INTEGER :: ii,jj,kk,i,j,k,l,n 
INTEGER :: kstep,nstep,np0,np1

!声明整型和双精度浮点型变量用于循环计数和存储数据，input用于存储输入数据，output用于存储输出的复数数据
REAL*8 :: input(nx,ny,nz)
COMPLEX*16 :: output(nx21,ny,nz)

! pp:从文件 given_structure.dat 中读取,赋值给phip
! px:read eta from restart.dat,赋值给eta，px,py暂未使用
REAL*8 :: px,py,pz,pp,ps,pm,pv

REAL*8 :: x11m, x12m, x13m, x21m, x22m, x23m, x31m, x32m, x33m
REAL*8 :: x11p, x12p, x13p, x21p, x22p, x23p, x31p, x32p, x33p
REAL*8 :: x11s, x12s, x13s, x21s, x22s, x23s, x31s, x32s, x33s
REAL*8 :: x11c, x12c, x13c, x21c, x22c, x23c, x31c, x32c, x33c

REAL*8 :: Eksm(3,3),Eksp(3,3),Ekss(3,3),Eksv(3,3),Eksc(3,3),Km(3,3),Kp(3,3),deta(3,3)
REAL*8 :: Em_break,Es_break,Ep_break,Ev_break

!radius:示生成的粒子的半径，该变量未使用
!R1:核的半径,R2:核和壳层在内的总半径
!E_ext(3): 这是一个三维数组，用于表示外部电场在 x,y,z 方向上的分量
!E_ext 的各个分量从输入文件input.in中读取，并随后用于模拟外部电场对材料的影响
REAL*8 :: radius,E_ext(3),R1,R2

!shell_thickness: 表示壳层的厚度
!circle(3): 长度为3的数组，设置其中心坐标，保存圆的中心坐标 x, y 和半径 R1
!circles(256,4): 这是一个二维数组，存储了多个圆的相关信息。数组的每一行表示一个圆的信息，包含圆心的 x,y 坐标
!以及内径 R1 和外径 R2。在生成多个填料或孔洞时，circles 用于记录每个圆的参数，以便后续检查是否重叠
!dd: 用于存储两个圆心之间的距离。在生成新的圆时，这个变量用于检查新生成的圆是否与已有的圆重叠
REAL*8 :: shell_thickness,circle(3),circles(256,4),dd

!ic：生成填料的数量
!nn: 用于遍历所有生成的粒子并更新相应的相场变量
!sample: 可能用于表示当前的样本编号，特别是在需要对多个样本进行模拟时，该变量仅出现在被注释的代码中
!num_samples: 表示样本的总数量
!particle_num: 表示填料（粒子）的总数量
INTEGER :: ic,nn,sample,num_samples,particle_num

REAL*8 :: dt,CA,L0,f0,g11,ifbreak(nx,ny,nz),Ebreak(nx,ny,nz),temp(nx,ny,nz)
REAL*8 :: ave_D(3)

!phi(:,:,:): 三维数组，用于表示整个模拟空间的总相场变量,可视化用
!phi1(:,:,:): 三维数组，用于表示第一种相（如核、内层填料等）在模拟空间中的分布。在生成结构时，phi1 用于初始化相场变量
!phi2(:,:,:): 三维数组，用于表示第二种相（如壳层、外层填料等）在模拟空间中的分布。通常与 phi1 一起用于描述不同相的界面
REAL*8,ALLOCATABLE :: phi(:,:,:),phi1(:,:,:),phi2(:,:,:)

!phip 是一个三维数组，表示整个模拟空间中填料相的总相场分布，matrix,shell,void?
REAL*8,ALLOCATABLE :: phip(:,:,:),phim(:,:,:),phis(:,:,:),phiv(:,:,:)

!phicp 是一个四维数组，用于存储生成的多个填料相（particle phase）在模拟空间中的相场信息
!每个 phicp(ic, :, :, :) 表示第 ic 个填料相的三维相场分布
REAL*8,ALLOCATABLE :: phicp(:,:,:,:),phicm(:,:,:,:),phics(:,:,:,:)
REAL*8,ALLOCATABLE :: el(:,:,:,:)
REAL*8,ALLOCATABLE :: Khom(:,:),Kapha(:,:,:,:,:),Eks(:,:,:,:,:),dKdeta(:,:,:,:,:)
REAL*8,ALLOCATABLE :: eta(:,:,:),force(:,:,:)
COMPLEX*16,ALLOCATABLE :: etak(:,:,:),forcek(:,:,:)
REAL*8,ALLOCATABLE :: G_elec(:,:,:),Gc(:,:,:)
!表示相分离能量、梯度能量和电能
REAL*8,ALLOCATABLE :: fsepar(:,:,:),fgrad(:,:,:),felec(:,:,:)

!声明并分配动态数组，用于存储模拟过程中的各种物理量，phi（相场变量）、el（应变张量）、eta（断裂相变量）
ALLOCATE(phi(nx,ny,nz),phi1(nx,ny,nz),phi2(nx,ny,nz))
ALLOCATE(phip(nx,ny,nz),phim(nx,ny,nz),phis(nx,ny,nz),phiv(nx,ny,nz))
ALLOCATE(phicp(256,nx,ny,nz),phicm(256,nx,ny,nz),phics(256,nx,ny,nz))
ALLOCATE(el(3,nx,ny,nz))
ALLOCATE(Khom(3,3),Kapha(3,3,nx,ny,nz),Eks(3,3,nx,ny,nz),dKdeta(3,3,nx,ny,nz))
ALLOCATE(eta(nx,ny,nz),force(nx,ny,nz))
ALLOCATE(etak(nx21,ny,nz),forcek(nx21,ny,nz))
ALLOCATE(G_elec(nx,ny,nz),Gc(nx,ny,nz))
ALLOCATE(fsepar(nx,ny,nz),fgrad(nx,ny,nz),felec(nx,ny,nz))

!--------------------------------------------------------------------------------------------------------------------
!读input.in文件部分

OPEN(UNIT = 1,FILE = 'input.in')
! 时间步长dt、步数nstep
READ(1,*)dt,nstep   
!指示是否从已有结构文件中读取结构或者创建新的结构
READ(1,*)np0,np1  
! ?
READ(1,*)g11
! 外部电场在x,y,z方向的分量
READ(1,*)E_ext(1),E_ext(2),E_ext(3) 
! 表示不同相（如基体 m，填料 p，壳层 s，和其他 c）的介电常数, c很可能表示 Composite（复合材料）的整体或 Core（核心部分）
READ(1,*)x11m,x12m,x13m   
READ(1,*)x21m,x22m,x23m
READ(1,*)x31m,x32m,x33m
READ(1,*)x11p,x12p,x13p
READ(1,*)x21p,x22p,x23p
READ(1,*)x31p,x32p,x33p
READ(1,*)x11s,x12s,x13s
READ(1,*)x21s,x22s,x23s
READ(1,*)x31s,x32s,x33s
READ(1,*)x11c,x12c,x13c
READ(1,*)x21c,x22c,x23c
READ(1,*)x31c,x32c,x33c
! 填料的数量和样本的数量
READ(1,*)particle_num,num_samples
CLOSE(1)

E_ext = E_ext*1.0E7 ! 对外部电场进行单位转换（例如，从V/m转换到更高的单位）
L0 = 1.0            ! 动力学系数的初始值
f0 = g11            ! ?可能表示与材料特性相关的系数
CA = 0.5            ! ?可能是一个系数或比例因子
Em_break = 3.7E8    ! 基体相的击穿场强
Es_break = 1.0E8    ! 壳层相的击穿场强
Ep_break = 2.0E7    ! 填料相的击穿场强
Ev_break = 3.0E6    ! ?可能是空隙或其他相的击穿场强

Eksm(1,1) = x11m
Eksm(1,2) = x12m
Eksm(1,3) = x13m
Eksm(2,1) = x21m
Eksm(2,2) = x22m
Eksm(2,3) = x23m
Eksm(3,1) = x31m
Eksm(3,2) = x32m
Eksm(3,3) = x33m

Eksp(1,1) = x11p
Eksp(1,2) = x12p
Eksp(1,3) = x13p
Eksp(2,1) = x21p
Eksp(2,2) = x22p
Eksp(2,3) = x23p
Eksp(3,1) = x31p
Eksp(3,2) = x32p
Eksp(3,3) = x33p

Ekss(1,1) = x11s
Ekss(1,2) = x12s
Ekss(1,3) = x13s
Ekss(2,1) = x21s
Ekss(2,2) = x22s
Ekss(2,3) = x23s
Ekss(3,1) = x31s
Ekss(3,2) = x32s
Ekss(3,3) = x33s

Eksv = 0.0

Eksc(1,1) = x11c
Eksc(1,2) = x12c
Eksc(1,3) = x13c
Eksc(2,1) = x21c
Eksc(2,2) = x22c
Eksc(2,3) = x23c
Eksc(3,1) = x31c
Eksc(3,2) = x32c
Eksc(3,3) = x33c

!对初始化的材料特性矩阵进行单位转换。为了确保所有的物理量具有一致的单位制，以避免在计算过程中出现数值误差
Eksm = Eksm*1.0E-9  
Eksp = Eksp*1.0E-9
Ekss = Ekss*1.0E-9
Eksv = Eksv*1.0E-9
Eksc = Eksc*1.0E-9

!FFTW plans for forward and backward transforms创建3D快速傅里叶变换的计划，用于从实数到复数和从复数到实数的转换
!（p_up 和 p_dn）定义了如何将数据从实空间转换到傅里叶空间（r2c）以及如何从傅里叶空间转换回实空间（c2r）
!nx, ny, nz: 定义了FFT的维度大小。input, output: 指定输入和输出数据数组。
!FFTW_MEASURE: 一个选项，用于指示FFTW库在生成计划时执行一些测试，以找到最快的执行路径
CALL dfftw_plan_dft_r2c_3d(p_up, nx, ny, nz, input, output, FFTW_MEASURE)
CALL dfftw_plan_dft_c2r_3d(p_dn, nx, ny, nz, output, input, FFTW_MEASURE) 
CALL periodic !设置周期性边界条件

!Define the Dirac function data定义狄拉克函数数据
!这里定义了一个三维的狄拉克函数（或更精确地说是一个小值矩阵），deta 作为一个3x3矩阵被初始化为零，然后对角线元素被设置为非常小的数（1.0E-9）。
!狄拉克函数（Dirac Delta Function） 在数值模拟中通常用于近似极小的特征，或者模拟在某些计算中几乎为零但不完全为零的情况
deta = 0.0
DO j=1,3
 DO i=1,3
  IF(i == j) deta(i,j) = 1.0E-9
 END DO
END DO

!初始化为零通常表示系统从一个均匀的初始状态开始，
!或者表明这些相场在初始化时尚未被激活（例如，所有填料、基体、壳层和空隙的初始浓度为零）
phip = 0.0  
phim = 0.0
phis = 0.0
phiv = 0.0
!GOTO 100

!------------------------------------------------------------------------------
!读结构或生成结构部分
!DO 2000 sample = 1,num_samples
!read from a given structure
!如果 np0 == 1，从文件 given_structure.dat 中读取预定义的结构数据
IF(np0 == 1) THEN 
 OPEN(UNIT = 1, FILE = 'given_structure.dat') 
  DO n=1,nx*ny*nz 
   !网格点的索引 (i, j, k) 和相应的相场值 (pp, ps, pm, pv)
   READ(1,*)i,j,k,pp,ps,pm,pv
   ii = i
   jj = j
   kk = k
   phip(ii,jj,kk) = pp
   phis(ii,jj,kk) = ps
   phim(ii,jj,kk) = pm
   phiv(ii,jj,kk) = pv
  END DO
 CLOSE(1)

ELSE
!Start to create a sample with random microstructure
!将用于生成微观结构的相场变量 phicp, phicm, phics 初始化为零，shell_thickness（壳层厚度）初始化为零
 phicp = 0.0
 phicm = 0.0
 phics = 0.0
 shell_thickness = 0.0 
 ic = 1 !初始化计数器 ic 为1，用于跟踪生成的粒子数量

!GOTO 111
!CALL create_random(circle)
!circle(1) = circle(1)*(nx - 24) + 12  !coordinate of x
!circle(2) = circle(2)*(ny - 24) + 12  !coordinate of y
!circle(3) = (circle(3) + 1.0)*16.0         !random Radius
!circle(3) = 8.0         !random Radius

!Save the first circle to circles
!初始化一个圆 circle，设置其中心坐标为网格中心 (nx/2, ny/2)，半径 R1 设置为 4， R2 设置为壳层厚度。
 circle(1) = DBLE(nx/2)  
 circle(2) = DBLE(ny/2)
 circle(3) = 4  !第一个圆的半径为4，后面随机生成的为8 ？？？
 R1 = circle(3)
 R2 = circle(3) + shell_thickness
!将初始圆的参数存储到数组 circles 中，以便后续参考
 circles(ic,1) = circle(1)
 circles(ic,2) = circle(2)
 circles(ic,3) = R1 
 circles(ic,4) = R2 

 !调用 iniconf 函数初始化相应的相场变量 phi1 和 phi2
 CALL iniconf(INT(circle(1)),INT(circle(2)),nz,R1,phi1)  
 CALL iniconf(INT(circle(1)),INT(circle(2)),nz,R2,phi2)

 !更新 phim、phis 和 phip，存储在第ic个 phicm、phics 和 phicp 中
 phim = 1.0 - phi2
 phis = phi2 - phi1
 phip = phi1
 !将三维数组 phip 赋值给四维数组 phicp 的一个特定切片（即 ic 确定的那一层）
 phicm(ic,:,:,:) = phim(:,:,:)  
 phics(ic,:,:,:) = phis(:,:,:)
 phicp(ic,:,:,:) = phip(:,:,:)
 

!检查生成的粒子数量是否达到要求 (particle_num)，如果达到则跳转到 112 
111 IF(ic >= particle_num) GOTO 112 
 !调用 create_random 函数生成随机圆
 CALL create_random(circle)  
 ! PRINT*,circle 
 !circle(1) 和 circle(2) 在调用 create_random(circle) 后，通常会返回一个范围在 0到1 之间的随机值
 !为了将这些随机值转换为实际的坐标值，需要进行缩放和偏移
 !缩放: circle(1) * (nx - 20) 将随机值放大到一个合理的坐标范围内。(nx - 20) 表示将坐标限制在模拟网格的边界内
 !减少 20 个单位的边缘空间，这有助于避免生成的圆过于靠近边界
 !偏移: 加上 10 是为了确保圆心不会在 x 或 y 坐标为 0 的位置
 !通过偏移，生成的圆心坐标范围在 10 到 nx-10 之间，这样可以避免圆心直接在边界上
  circle(1) = circle(1)*(nx - 20) + 10  !coordinate of x
  circle(2) = circle(2)*(ny - 20) + 10  !coordinate of y
  !circle(3) = (circle(3) + 1.0)*16.0         !random Radius
  circle(3) = 8.0         !random Radius?这里直接将半径设置为 8.0，可能是为了简化模型，在生成的微观结构中保持一致的填料大小
  R1 = circle(3)
  R2 = circle(3) + shell_thickness

  !Determine if the created circle is wanted
  DO i=1,ic 
   !遍历所有已经生成的圆（或球），计算新生成的圆心与每个已有圆心之间的距离 dd
   dd = SQRT( (circle(1) - circles(i,1))**2 + (circle(2) - circles(i,2))**2 )
   !如果新生成的圆与任何已有的圆重叠（距离 dd 小于等于两个圆的R2和加上一个偏置值 14），则放弃新圆并跳回 111 重新生成一个新的随机圆
   IF(dd <= (circles(i,4)+R2+14)) THEN
   !PRINT*,'Give up the unwanted and recreate'
    GOTO 111
   ENDIF
  END DO

  !如果新圆通过了重叠检查，将其信息保存到 circles 数组中。
  !ic 增加 1，表示新的粒子已生成
  ic = ic + 1
  !Save the wanted circle to circles
  !使用 iniconf 子程序初始化相应的相场变量 phi1 和 phi2，代表不同半径下的相场值
  !phim（基体相），phis（壳层相），phip（填料相）相应地被更新，表示生成的圆形区域的不同相态
  !这些相场变量被存储到四维数组 phicm, phics, phicp 中，记录每个粒子的状态
  circles(ic,1) = circle(1)
  circles(ic,2) = circle(2)
  circles(ic,3) = R1
  circles(ic,4) = R2
  CALL iniconf(INT(circle(1)),INT(circle(2)),nz,R1,phi1)  
  CALL iniconf(INT(circle(1)),INT(circle(2)),nz,R2,phi2)

  !phim 表示基体相的相场变量。这里，1.0 表示整个模拟空间的总相场值（即总和为1），phi2 是包含壳层和填料相的区域
  !1.0 - phi2 确定了基体相的存在区域。基体相存在于填料相和壳层相之外的空间，即 phi2 所定义的区域之外
  !通过这种计算，phim 得到了基体相在模拟空间中的分布，确保所有相的相场变量之和为1
  phim = 1.0 - phi2
  !phis 表示壳层相的相场变量，通过 phi2 - phi1，phis 定义了壳层相在填料相之外的区域。即在 phi2 和 phi1 之间的区域是壳层相存在的地方
  phis = phi2 - phi1
  !phi1 表示的是填料相的位置及其空间分布，在相场模型中，它的值通常在不同相的区域之间过渡（例如，中心区域值为1，边界处过渡到0）
  phip = phi1
  !存储第 ic 个粒子的相场
  phicm(ic,:,:,:) = phim(:,:,:) 
  phics(ic,:,:,:) = phis(:,:,:)
  phicp(ic,:,:,:) = phip(:,:,:)
 !PRINT*,'ic= ',ic
 GOTO 111


 !将 phi, phip, phim, phis, phiv 初始化为零，准备进行全局更新
 112 phi = 0.0 
 phip = 0.0
 phim = 0.0
 phis = 0.0
 phiv = 0.0
 DO nn=1,particle_num 
  !使用循环遍历所有粒子编号 nn
  !如果 nn 是 200 的倍数，表示当前粒子被认为是一个孔洞（或空隙），其相场被累加到 phiv
  !否则，将其相场累加到 phip（填料相）和 phis（壳层相）
   IF(MOD(nn,200) == 0) THEN   
    phiv(:,:,:) = phiv(:,:,:) + phicp(nn,:,:,:)
   ELSE
    phip(:,:,:) = phip(:,:,:) + phicp(nn,:,:,:)
    phis(:,:,:) = phis(:,:,:) + phics(nn,:,:,:)
   ENDIF
 END DO
 !使用 phim = 1.0 - phip - phis - phiv，确保所有相场变量之和为1，符合相场模型的要求
 phim(:,:,:) = 1.0 - phip(:,:,:) - phis(:,:,:) - phiv(:,:,:)  
ENDIF
!End of creating a sample with random microstructure 停止生成

!--------------------------------------------------------------------------------------------------
!写入结构数据到 structure.dat 文件
OPEN(UNIT = 2, FILE = 'structure.dat') 
  DO kk=1,nz
   DO jj=1,ny
    DO ii=1,nx
     ! 2002 是一个格式描述符，在代码末尾部分定义了格式
     !每一行的数据包括网格点坐标x,y,z和相场值（填料相 phip，壳层相 phis，基体相 phim，空隙相 phiv）
     WRITE(2,2002)ii,jj,kk,phip(ii,jj,kk),phis(ii,jj,kk),phim(ii,jj,kk),phiv(ii,jj,kk)
    END DO
   END DO
  END DO

!------------------------------------------------------------------------------------------------------
!对phi 进行赋值操作，将 phim, phis, phip, 和 phiv 乘以不同的权重（4.0，8.0，12.0，16.0）累加到 phi 中
!这些权重的选择可能用于区分不同的相，以便在后续的可视化过程中能够更容易地识别和分离不同的材料相
 phi(:,:,:) = phi(:,:,:) + 4.0*phim(:,:,:) + 8.0*phis(:,:,:) + 12.0*phip(:,:,:) + 16.0*phiv(:,:,:)

!100 DO kk=1,nz
!  DO jj=1,ny
!   DO ii=1,nx
!    IF((ABS(ii-nx/2) <= 8).AND.(ABS(jj-ny/2) <= 8)) THEN
!     phip(ii,jj,kk) = 1.0
!    ELSE
!     phim(ii,jj,kk) = 1.0
!    ENDIF
!   END DO
!  END DO
! END DO 

!将 phi 的第一个元素设置为0
!写入 phi 数组到 shape.dat 文件
 phi(1,1,1) = 0.0 
 OPEN(UNIT = 2, FILE = 'shape.dat') 
  DO ii=1,nx
   DO jj=1,ny
    !将当前网格点的坐标 (ii, jj) 和相应的 phi 值写入文件
    WRITE(2,*)ii,jj,phi(ii,jj,nz)
   END DO
   !用于换行，使得数据在文件中呈现出二维网格的结构
   WRITE(2,*)
  END DO
 CLOSE(2)

!------------------------------------------------------------------------------------------------------
!计算每个网格点的 Eks （介电常数）和 Gc（电能量密度）
 DO kk=1,nz 
  DO jj=1,ny
   DO ii=1,nx
    DO l=1,3
     DO k=1,3
      Eks(k,l,ii,jj,kk) = Eksm(k,l)*phim(ii,jj,kk) + Ekss(k,l)*phis(ii,jj,kk) + & 
                          Eksp(k,l)*phip(ii,jj,kk) + Eksv(k,l)*phiv(ii,jj,kk)
     END DO
    END DO
    Gc(ii,jj,kk) = 0.5*epson0*((Eksm(1,1)+deta(1,1))*Em_break**2*phim(ii,jj,kk) + & 
                               (Ekss(1,1)+deta(1,1))*Es_break**2*phis(ii,jj,kk) + &
                               (Eksp(1,1)+deta(1,1))*Ep_break**2*phip(ii,jj,kk) + &
                               (Eksv(1,1)+deta(1,1))*Ev_break**2*phiv(ii,jj,kk))
!   Ebreak(ii,jj,kk) = Em_break*phim(ii,jj,kk) + Es_break*phis(ii,jj,kk) + Ep_break*phip(ii,jj,kk) + Ev_break*phiv(ii,jj,kk)
!   Ebreak(ii,jj,kk) = Em_break
   END DO
  END DO
 END DO
 
 Km = Eksm + deta !初始化 Km 和 Kp
 Kp = Eksp + deta
 
!eta(ii,128,1),ii取整122-134，设这些点已被击穿
 eta = 0.0
 DO kk=1,nz
  DO jj=1,ny
   DO ii=1,nx
    IF((ABS(ii - nx/2) <= 6).AND.(ABS(jj - ny/2) < 1)) eta(ii,jj,kk) = 1.0 
   END DO
  END DO
 END DO
 
 IF(np1.EQ.1) THEN   !read eta from restart.dat
  OPEN(UNIT = 11, FILE = 'restart.dat')
   DO n=1,nx*ny*nz
    READ(11,*)i,j,k,px
    ii = i
    jj = j
    kk = k
    eta(ii,jj,kk) = px
   END DO
  CLOSE(11)
 ENDIF
 
 DO kk=1,nz
  DO jj=1,ny
   DO ii=1,nx
    DO l=1,3
     DO k=1,3 
      ! 公式(11)，空间依赖相对介电常数
      Kapha(k,l,ii,jj,kk) = (eta(ii,jj,kk)**3*(10.0 - 15.0*eta(ii,jj,kk) + 6.0*eta(ii,jj,kk)**2))* &
                             (Eksc(k,l) - Eks(k,l,ii,jj,kk)) + (Eks(k,l,ii,jj,kk) + deta(k,l))
      ! K对eta求导
      dKdeta(k,l,ii,jj,kk) = 30.0*eta(ii,jj,kk)**2*(eta(ii,jj,kk) - 1.0)**2*(Eksc(k,l) - Eks(k,l,ii,jj,kk))
     END DO
    END DO
   END DO
  END DO
 END DO

 DO jj=1,ny !将 Kapha 和 dKdeta 以及 Gc 的值写入文件
  DO ii=1,nx
   WRITE(12,*)ii,jj,Kapha(1,1,ii,jj,1),dKdeta(1,1,ii,jj,1)
    WRITE(14,*)ii,jj,Gc(ii,jj,1)
  END DO
  WRITE(12,*)
   WRITE(14,*)
 END DO
! PRINT*,'Kapha = ',Kapha(1,1,nx/2,ny/2,1)
 
 !Initilize homogenous dielectric constants初始化均匀的介电常数 Khom，王建军的一篇文章，初始的均匀状态下取的介电系数
 IF(Kp(1,1) >= Km(1,1)) Khom = Kp
 IF(Km(1,1) >= Kp(1,1)) Khom = Km
!Khom(:,:) = Kapha(:,:,nx/2,ny/2,nz)
 Khom = (Eksc + deta)*1.0

 OPEN(UNIT = 7, FILE = 'E-D.dat') !打开文件 E-D.dat
 DO 1000 kstep = 1,nstep  !进行 nstep 次循环，更新外部电场 E_ext

  E_ext(1) = 1.0E8 + DBLE(kstep*dt*1.0E6) 
!调用 separation、gradient 和 electric 子程序计算相应的相场能量
  CALL separation(eta,f0,fsepar)  
  CALL gradient(eta,g11,fgrad)  
  CALL electric(Kapha,Khom,dKdeta,E_ext,eta,felec,el,G_elec)
!计算平均电位移 ave_D 并将其写入文件
  ave_D = 0
  DO kk=1,nz
   DO jj=1,ny
    DO ii=1,nx
     DO j=1,3
      DO i=1,3
       ave_D(j) = ave_D(j) + epson0*Kapha(i,j,ii,jj,kk)*E_ext(i)  !电位移
      END DO
     END DO
    END DO
   END DO
  END DO
  IF(MOD(kstep,10) == 0) WRITE(7,*)E_ext(1),ave_D(1)
!  DO jj=1,ny
!   DO ii=1,nx
!    WRITE(13,*)ii,jj,G_elec(ii,jj,1)
!   END DO
!   WRITE(13,*)
!  END DO
!  CLOSE(13)
 
! force = fsepar + fgrad - 1.0*felec
  force = fsepar + fgrad + 1.0*felec  !计算总力 force 为 fsepar、 fgrad  felec 之和，这三个必须，其他考虑可再加
! DO kk=1,nz
!  DO jj=1,ny
!   DO ii=1,nx
!    WRITE(15,*),ii,jj,kk,force(ii,jj,kk)
!   END DO
!   WRITE(15,*)
!  END DO
! END DO
! CLOSE(15)
 
! PRINT*,'fsepar = ',fsepar(nx/2,ny/2,nz),'fgrad = ',fgrad(nx/2,ny/2,nz),'felec = ',felec(nx/2,ny/2,nz)
  
  DO kk=1,nz
   DO jj=1,ny
    DO ii=1,nx  !遍历所有网格点，比较 G_elec 和 Gc(critical) 的值
     IF(G_elec(ii,jj,kk) < Gc(ii,jj,kk)) THEN
!    IF(MAX1(ABS(el(1,ii,jj,kk)),ABS(el(2,ii,jj,kk)),ABS(el(3,ii,jj,kk))) < Ebreak(ii,jj,kk)) THEN
      force(ii,jj,kk) = 0.0 !如果 G_elec 小于 Gc，则将力 force （驱动力）设为0，并将 ifbreak 设为0
      ifbreak(ii,jj,kk) = 0.0 !是否击穿
     ELSE
      force(ii,jj,kk) = L0*force(ii,jj,kk)  !否则，将力 force 乘以系数 L0 并将 ifbreak 设为1
      ifbreak(ii,jj,kk) = 1.0
     ENDIF
!     WRITE(16,*),ii,jj,kk,force(ii,jj,kk)
!     WRITE(17,*),ii,jj,kk,G_elec(ii,jj,kk),felec(ii,jj,kk)
    END DO
!    WRITE(16,*)
!    WRITE(17,*)
   END DO
  END DO
!  CLOSE(16)
!  CLOSE(17)
 !调用 forward 函数进行傅里叶变换，将 force 和 eta 变换为频域表示 forcek 和 etak
  CALL forward(force,forcek)

  CALL forward(eta,etak)
  DO kk=1,nz
   DO jj=1,ny
    DO ii=1,nx21  !在频域中更新 etak，导空间演化的方式
     etak(ii,jj,kk) = etak(ii,jj,kk) - dt*forcek(ii,jj,kk)/(1.0 + 0.5*g11*kpow2(ii,jj,kk)*dt)
!    etak(ii,jj,kk) = etak(ii,jj,kk) - dt*forcek(ii,jj,kk)
    END DO
   END DO
  END DO
  CALL backward(etak,eta) !调用 backward 函数进行反傅里叶变换，更新 eta
! DO kk=1,nz
!  DO jj=1,ny
!   DO ii=1,nx
!    IF((ABS(ii - nx/2) <= 8).AND.(ABS(jj - ny/2) <= 2)) eta(ii,jj,kk) = 1.0
!   END DO
!  END DO
! END DO
!更新 Kapha 和 dKdet，有的地方击穿了，重新计算
  DO kk=1,nz
   DO jj=1,ny
    DO ii=1,nx
     DO l=1,3
      DO k=1,3
       Kapha(k,l,ii,jj,kk) = (eta(ii,jj,kk)**3*(10.0 - 15.0*eta(ii,jj,kk) + 6.0*eta(ii,jj,kk)**2))* &
                             (Eksc(k,l) - Eks(k,l,ii,jj,kk)) + (Eks(k,l,ii,jj,kk) + deta(k,l))
       dKdeta(k,l,ii,jj,kk) = 30.0*eta(ii,jj,kk)**2*(eta(ii,jj,kk) - 1.0)**2*(Eksc(k,l) - Eks(k,l,ii,jj,kk))
      END DO
     END DO
    END DO
   END DO
  END DO  !打印当前时间步 kstep、外部电场 E_ext 和 eta 的值
  PRINT*,'kstep = ',kstep,'E1 = ',E_ext(1), 'E2 = ',E_ext(2),'eta = ',eta(40+nx/2,ny/2,nz)
  IF((MOD(kstep,100) == 0).OR.(kstep == 1)) THEN !每100步或在第1步时，将 eta 的值用于更新 temp，并将 phi 和 el 的值写入文件
   DO ii=1,nx
    DO jj=1,ny
     IF(eta(ii,jj,1) < 0.5) temp(ii,jj,1) = 0.0 
     IF(eta(ii,jj,1) >= 0.5) temp(ii,jj,1) = 1.0 
     IF((ABS(ii - nx/2) <= 6).AND.(ABS(jj - ny/2) < 1)) temp(ii,jj,kk) = 1.0
     WRITE(kstep,2008)ii,jj,(1.0 - temp(ii,jj,1))*phi(ii,jj,1),el(1,ii,jj,1)
     WRITE(16,*),ii,jj,force(ii,jj,1),ifbreak(ii,jj,1)
     WRITE(17,*),ii,jj,G_elec(ii,jj,1),felec(ii,jj,1),el(1,ii,jj,1)
    END DO
    WRITE(kstep,*)
    WRITE(16,*)
    WRITE(17,*)
   END DO
   CLOSE(kstep)
   CLOSE(16)
   CLOSE(17)
  END IF

 1000 CONTINUE
 ! WRITE(7)

!2000 CONTINUE
 
OPEN(UNIT = 6, FILE = 'field.dat')  !打开文件 field.dat、eta.dat 和 morphology.dat
OPEN(UNIT = 7, FILE = 'eta.dat')
OPEN(UNIT = 8, FILE = 'morphology.dat')
 DO kk=1,nz
  DO ii=1,nx
   DO jj=1,ny !将 el、eta 和形态数据写入相应文件
    WRITE(6,2006)ii,jj,kk,el(1,ii,jj,kk),el(2,ii,jj,kk),el(3,ii,jj,kk)
    WRITE(7,2007)ii,jj,kk,eta(ii,jj,kk)
    WRITE(8,2007)ii,jj,kk,(1.0 - eta(ii,jj,kk))*phi(ii,jj,kk)
   END DO
   WRITE(6,*)
   WRITE(7,*)
   WRITE(8,*)
  END DO
 END DO
CLOSE(6)
CLOSE(7)

!--------------------------------------------------------------------------------------------------------
!定义文件写入格式
! 3(1x,i4):
! 3(...): 表示括号内的格式重复三次
! 1x: 表示插入一个空格（即，水平跳过一个字符位置）
! i4: 表示以整数格式输出，宽度为 4 个字符
! 3(1x,i4) 代表输出三个整数，每个整数之前有一个空格，总共输出 3 个整数，每个整数占用 4 个字符宽度
! e12.4: 表示以科学计数法（指数形式）输出一个实数，总宽度为 12 个字符，其中包含 4 位小数
! 4(1x,e12.4) 代表输出四个实数，每个实数之前有一个空格，总共输出 4 个实数，每个实数占用 12 个字符宽度
2002  FORMAT(3(1x,i4),4(1x,e12.4))
2006  FORMAT(3(1x,i4),3(1x,e12.4))
2007  FORMAT(3(1x,i4),(1x,e12.4))
2008  FORMAT(2(1x,i4),2(1x,e12.4))

!-----------------------------------------------------------------------------------------------------
!释放之前分配的所有数组内存
DEALLOCATE(phi,phi1,phi2)
DEALLOCATE(phip,phim,phis,phiv)
DEALLOCATE(phicp,phicm,phics)
DEALLOCATE(el)
DEALLOCATE(Khom,Kapha,Eks,dKdeta)
DEALLOCATE(eta,force)
DEALLOCATE(etak,forcek)
DEALLOCATE(G_elec,Gc)
DEALLOCATE(fsepar,fgrad,felec)

END PROGRAM main
