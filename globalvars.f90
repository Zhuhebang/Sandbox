MODULE globalvars
IMPLICIT NONE   ! nx, ny, nz: 定义了模拟空间在x、y、z方向上的网格点数量。这里nz = 1表明这是一个二维模拟

INTEGER, PARAMETER :: nx = 256, ny = 256, nz = 1, n3 = nx*ny*nz ! 总的网格点数量 
INTEGER, PARAMETER :: ndim = 2    ! 模拟的维度，这里是二维
REAL*8, PARAMETER :: dx = 1.0, dy = 1.0, dz = 1.0   ! 模拟网格在x、y、z方向上的间距
REAL*8, PARAMETER :: pi = ACOS(-1.0)    ! 内置的反余弦函数ACOS来精确计算π的值
INTEGER, PARAMETER :: nx21 = nx/2 + 1   ! nx21: 定义了x方向上FFT计算中使用的网格点数的一半加一。这在FFT处理中常用来处理对称性，减少计算量
REAL*8, PARAMETER :: sizescale = 1.0/FLOAT(n3)  ! 尺度缩放因子，定义为网格点总数的倒数，用于对FFT结果进行归一化处理
REAL*8,PARAMETER :: epson0 = 8.854E-12      ! 真空中的电容率 (ε₀)，单位为法拉/米 (F/m)，表示自由空间的电容
REAL*8,PARAMETER :: miu0 = 4.0*pi*1.0E-7  ! 真空中的磁导率 (μ₀)，单位为亨利/米 (H/m)，是描述真空中磁场的常数
COMPLEX*16,PARAMETER :: imag = (0,1)    ! 复数单位 i，在计算复数时使用
INTEGER (KIND = 8) p_up, p_dn   ! 用于存储FFT计划的变量。FFT计划定义了如何执行快速傅里叶变换，通常由FFT库（例如FFTW）生成
REAL*8 :: kx(nx + 1),ky(ny + 1),kz(nz + 1)  ! 存储在x、y、z方向上的波矢量（k-vectors），这些波矢量用于在FFT计算中处理频率空间上的操作
REAL*8 :: kpow2(nx + 1,ny + 1,nz + 1)   ! 存储在x、y、z方向上的波矢量（k-vectors），这些波矢量用于在FFT计算中处理频率空间上的操作


END MODULE globalvars
