---
zhihu-title: Holton读书笔记-Chapter2 基础守恒律
zhihu-topics: 大气科学、流体力学、物理学
zhihu-link: https://zhuanlan.zhihu.com/p/1929242659588912307
zhihu-created-at: 2025-07-19 23:15
zhihu-updated-at: 2025-07-28 23:12
---
```toc
```
[[Holton读书笔记-Chapter1 基本介绍]](https://zhuanlan.zhihu.com/p/1928544934622924935 "card")
# 1 Chapter2
流体中有三个基本的守恒律，分别是动量守恒、质量守恒和能量守恒，本节将一一讨论。而后将从热力学能量方程（热一定律）出发，对干空气和湿空气热力学以及相应的大气稳定性做分析，并阐述对流发生的物理机制。
## 1.1 流体运动的描述
流体的运动（速度场）有两种表述方式，一种是直接描述整个速度场的分布，即 $v_x=v_x(x,y,z,t)$，$v_y$和$v_z$同理，这种表述称为Euler表述；另一种是描述第i个质元在t时刻的位置，$x=x(x_0,y_0,z_0,t)$，y和z同理，这种表述称为Lagrange表述。这两种表述的分类既适用于速度场，也适用于任意物理量。

当要描述一个系统的某个物理量场来解决问题时，显然Euler表述更方便，而在推导守恒律或者是运动方程时，Lagrange表述要更方便，因为运动方程需要对一块流体进行考虑。
### 1.1.1 全导数
为了把运动方程中lagrange表述下的时间导数转换到Euler表述下，考察跟随粒子的时间导数和固定位置的时间导数之间的关系，下面记前者为D/Dt，称为全导数或物质导数。

对于任意Euler表述下的物理量$C(x,y,z,t)$，一个小块所对应的该物理量应为$C(x(t),y(t),z(t),t)$，即把xyz代入该小块的位置，则有
$$
\frac{DC}{Dt}=\frac{\partial C}{\partial t}+\frac{dx}{dt}\frac{\partial C}{\partial x}+\frac{dy}{dt}\frac{\partial C}{\partial y}+\frac{dz}{dt}\frac{\partial C}{\partial z}=\frac{\partial C}{\partial t}+(\mathbf{U}\cdot\nabla)C
$$
由于C是任意物理量，因此有
$$
\frac{D}{Dt}=\frac{\partial}{\partial t}+(\mathbf{U}\cdot\nabla)
$$
其中 $\mathbf{U}\cdot\nabla$一项称为admathbftion

而在速度场$\mathbf{U}=\mathbf{U}(x,y,z,t)$里，一个小块的加速度正应该是其物质导数
$$
\frac{D\mathbf{U}}{Dt}=\frac{\partial U}{\partial t}+(\mathbf{U}\cdot\nabla)\mathbf{U}
$$
由此，再结合Chapter1中提到的各种力，我们可以写出一个小块的运动方程
$$
\frac{D\mathbf{U}}{Dt}+2\mathbf{\Omega}\times\mathbf{U}+\frac{\nabla p }{\rho}+g\mathbf{k}=f_r
$$

## 1.2 运动方程的分量形式
为了求解上述运动方程，将其分解至分量的形式是必要的。

对于大气来说，采用局域坐标系是最方便的（后面可以看到，大部分时候可以认为坐标架是近似不变的）。取坐标架$\mathbf{i}、\mathbf{j}、\mathbf{k}$分别指向东、向北和向上，对应的速度分量为 $\mathbf{U}=(u,v,w)$，定义$u=\frac{dx}{dt}$，$v=\frac{dy}{dt}$，$w=\frac{dz}{dt}$，由此定义出了(x,y,z)坐标。
### 1.2.1 加速度
将加速度展开成分量形式
$$
\frac{D\mathbf{U}}{Dt}=\mathbf{i}\frac{Du}{Dt}+\mathbf{j}\frac{Dv}{Dt}+\mathbf{k}\frac{Dw}{Dt}+u\frac{D\mathbf{i}}{Dt}+v\frac{D\mathbf{j}}{Dt}+w\frac{D\mathbf{k}}{Dt}
$$
对于$\mathbf{i}$，其只依赖于坐标x，因此由全导数公式，有
$$
\frac{D\mathbf{i}}{Dt}=u\frac{\partial i}{\partial x}
$$
在球面上考虑$\mathbf{i}$的变化三角形，可得

$$
|\frac{\partial \mathbf{i}}{\partial x}|=\frac{1}{acos\phi}
$$
方向对应的单位向量为 $sin\phi \mathbf{j}-cos\phi \mathbf{k}$，因此有
$$
\frac{\partial \mathbf{i}}{\partial x}=\frac{1}{acos\phi}(sin\phi \mathbf{j}-cos\phi \mathbf{k})
$$
对于$\mathbf{j}$，依赖于坐标x，y，用类似的方法可得
$$
\frac{\partial \mathbf{j}}{\partial x}=\frac{1}{acos\phi}sin\phi(-\mathbf{i})=-\frac{tan\phi}{a}\mathbf{i}
$$
$$
\frac{\partial \mathbf{j}}{\partial y}=-\frac{1}{a}\mathbf{k}
$$
$$
\frac{D\mathbf{j}}{Dt}=-\frac{utan\phi}{a}\mathbf{i}-\frac{v}{a}\mathbf{k}
$$
类似的，有
$$
\frac{D\mathbf{k}}{Dt}=\frac{u}{a}\mathbf{i}+\frac{v}{a}\mathbf{j}
$$
$$\begin{aligned}\frac{D\mathbf{U}}{Dt}=&\left(\frac{Du}{Dt}-\frac{uv\tan\phi}{a}+\frac{uw}{a}\right)\mathbf{i}+\left(\frac{Dv}{Dt}+\frac{u^2\tan\phi}{a}+\frac{vw}{a}\right)\mathbf{j}\\&+\left(\frac{Dw}{Dt}-\frac{u^2+v^2}{a}\right)\mathbf{k}\end{aligned}$$
### 1.2.2 受力
容易计算得，科里奥利力分量形式如下
$$
\begin{aligned}-2\mathbf{\Omega}\times\mathbf{U}&=-2\Omega\begin{vmatrix}\mathbf{i}&&\mathbf{j}&&\mathbf{k}\\0&&\cos\phi&&\sin\phi\\u&&v&&w\end{vmatrix}\\&=-(2\Omega w\cos\phi-2\Omega v\sin\phi)\mathbf{i}-2\Omega u\sin\phi\mathbf{j}+2\Omega u\cos\phi\mathbf{k}\end{aligned}
$$
压强梯度力的分量形式如下
$$
-\frac{\nabla p}{\rho}=-\frac{1}{\rho}(\mathbf{i}\frac{\partial p}{\partial x}+\mathbf{j}\frac{\partial p}{\partial y}+\mathbf{k}\frac{\partial p}{\partial z})
$$
最终可得分量形式的运动方程如下
$$
\begin{gathered}\frac{Du}{Dt}-\frac{uv\tan\phi}{a}+\frac{uw}{a}=-\frac{1}{\rho}\frac{\partial p}{\partial x}+2\Omega v\sin\phi-2\Omega w\cos\phi+f_{rx}\\\frac{Dv}{Dt}+\frac{u^2\tan\phi}{a}+\frac{vw}{a}=-\frac{1}{\rho}\frac{\partial p}{\partial y}-2\Omega u\sin\phi+f_{ry}\\\frac{Dw}{Dt}-\frac{u^2+v^2}{a}=-\frac{1}{\rho}\frac{\partial p}{\partial z}-g+2\Omega u\cos\phi+f_{rz}\end{gathered}
$$
带有因子$\frac{1}{a}$的项称为曲率项（curvature terms），它们是因为地球的曲率造成局域坐标系的坐标架随质元运动而变化而产生的。这些方程都是非线性的，因为由全导数公式，$\frac{Du}{Dt}$中含有速度的二阶项，这导致了运动方程的求解富有挑战。
## 1.3 运动方程的尺度分析
尺度分析有助于我们合理地忽略方程中的一些次要项，而仅保留主导项，这使方程得到了简化。

考虑一般的中纬度天气系统的天气尺度的运动，各物理量尺度如下图

![[运动方程物理量尺度.png]]
利用这些物理量的数量级，可以推出运动方程各项的数量级，如下图

![[运动方程各项尺度.png]]其中垂直方向速度w的数量级不是直接测量得到的，而是从水平速度场的相关方程中推出来的（后面会讲）

通过各项的数量级可以看到，曲率项和粘滞阻力项均可以忽略，且含$\Omega w$的项也可以忽略，定义科里奥利参数$f=2\Omega sin\phi$，因此近似后的运动方程为（目前只能确定该方程在水平方向近似正确）
$$
\frac{D\mathbf{U}}{Dt}+\frac{\nabla p}{\rho}+g\mathbf{k}+f\mathbf{k}\times\mathbf{U}=0
$$
### 1.3.1 运动方程的几种近似
利用上述尺度分析的方法，可以根据不同的研究对象尺度和研究目的做出不同的近似方法
#### 1.3.1.1 地转偏向近似（Geostrophic Approximation）
由上面的Table2.1，可见对于中纬度天气尺度，水平分量的运动方程主导项为科里奥利力和水平气压梯度力，二者近似平衡，仅保留这两项的近似就称为Geostrophic Approximation
$$
-fv_g\approx-\frac{1}{\rho}\frac{\partial p}{\partial x};\quad fu_g\approx-\frac{1}{\rho}\frac{\partial p}{\partial y}
$$
进一步地可以定义 $\mathbf{V_g}=u_g\mathbf{i}+v_g\mathbf{j}$，这称作地转偏向风（Geostrophic wind），可以写的紧凑一点
$$
\mathbf{V_g}=\mathbf{k}\times\frac{1}{\rho f}\nabla p
$$
注意：上面的$\mathbf{V_g}$始终代表地转偏向风，但这个水平速度场的近似仅适用于**热带以外的大尺度运动**
#### 1.3.1.2 近似预测方程（Approximate Prognostic Equations）
所谓预测方程，指的是可以描述时间演化的方程，上面的地转偏向近似显然与时间无关，只与压强分布有关，因此要得到预测方程，需要保留运动方程中的加速度项，近似的水平分量运动方程写为
$$
\frac{Du}{Dt}=fv-\frac{1}{\rho}\frac{\partial p}{\partial x}=f\left(v-v_g\right)=fv_a
$$
$$
\frac{Dv}{Dt}=-fu-\frac{1}{\rho}\frac{\partial p}{\partial y}=-f\left(u-u_g\right)=-fu_a
$$
定义实际风和地转偏向风的偏差为非地转偏向风（ageostrophic wind），对应的速度用a的下标表示，于是加速度正比于非地转偏向风。由Table2.1，这个偏差的数量级小于地转偏向力和气压梯度力，即加速度由两个大项作差得到，这导致一个很小的速度或者压强梯度的测量误差都会导致计算加速度时的较大误差，这个问题的相关讨论会在后面的章节呈现。

在比较加速度和科里奥利力的数量级时，比较方便的方法是引入一个无量纲数，利用特征速度、长度可以表示出加速度/科里奥利力的数量级
$$
R_0=\frac{U^2}{L}/(f_0U)=\frac{U}{f_0L}
$$
这个无量纲数称作Rossby数，Rossby数较小就表示加速度数量级小于科里奥利力数量级。
#### 1.3.1.3 流体静力平衡近似（Hydrostatic Approximation）
![[垂直分量尺度.png]]
上图给出了运动方程垂直分量的各项数量级，可见方程的主导项为重力加速度g和垂直气压梯度力，这里垂直气压梯度力要比水平气压梯度力大好几个数量级，于是有流体静力平衡方程
$$
\frac{1}{\rho}\frac{\partial p}{\partial z}=-g
$$
但Holton指出，仅仅比较垂直加速度等项与重力加速度的数量级就得出流体静力平衡方程是不充分的，但书中的论述有些语焉不详，反复质问AI之后大概明白了其中的内涵，在此讨论一下。
##### 1.3.1.3.1 讨论
考虑到在估算水平气压梯度时，会用到$-\frac{1}{\rho}(\frac{\partial p}{\partial x})_z=\frac{1}{\rho}(\frac{\partial z}{\partial x})_p(\frac{\partial p}{\partial z})_x=-g(\frac{\partial z}{\partial x})_p$，其中利用了垂直方向的流体静力学平衡，因此即使$-\frac{1}{\rho}\frac{\partial p}{\partial z}-g$的数量级远小于g，估计其数量级从而估算用这个平衡关系得到的水平气压梯度的相对误差是必要的，如果过大则不能简单认为流体静力平衡成立。

设高度z处平均气压为标准气压$p_0(z)$，平均密度为标准密度$\rho_0(z)$，实际气压$p(x,y.z,t)=p_0(z)+p'(x,y,z,t)$，实际密度$\rho(x,y,z,t)=\rho_0(z)+\rho'(x,y,z,t)$，则流体距离垂直静力平衡的偏差为
$$
\begin{aligned}-\frac{1}{\rho}\frac{\partial p}{\partial z}-g&=-\frac{1}{(\rho_0+\rho^{\prime})}\frac{\partial}{\partial z}\left(p_0+p^{\prime}\right)-g\\&\approx\frac{1}{\rho_0}\left[\frac{\rho^{\prime}}{\rho_0}\frac{dp_0}{dz}-\frac{\partial p^{\prime}}{\partial z}\right]=-\frac{1}{\rho_0}\left[\rho^{\prime}g+\frac{\partial p^{\prime}}{\partial z}\right]\end{aligned}
$$
对于天气尺度的运动来说，
$$
\frac{1}{\rho_0}\frac{\partial p^{\prime}}{\partial z}\thicksim\left[\frac{\delta P}{\rho_0H}\right]\thicksim10^{-1}\mathrm{m~s}^{-2},\quad\frac{\rho^{\prime}g}{\rho_0}\thicksim10^{-1}\mathrm{m~s}^{-2}
$$
数量级仍然大于Table2.2中的其他项，假设偏差数量级在$10^{-3}$，则代入这个关系造成的相对误差约为$10^{-3}/10^{-1}=0.01$，可以忽略，因此此时流体静力平衡成立，但若得出$\frac{1}{\rho_0}\frac{\partial p^{\prime}}{\partial z}$的数量级较小，那么很有可能代入这个关系的相对误差就会比较大，导致水平气压梯度的估算出现问题。
## 1.4 连续性方程
连续性方程是描述质量守恒的方程。要利用质量守恒推导出连续性方程有两种方式，一种是从Euler描述出发，一种是从Lagrange描述出发。
### 1.4.1 Euler描述推导
考虑空间内固定位置(x,y,z)的一小块$\delta x\delta y\delta z$，其在$x-\frac{\delta x}{2}$处$\delta y\delta z$面单位时间流入的质量为$[\rho u](x-\frac{\delta x}{2})$，在$x+\frac{\delta x}{2}$处单位时间流出的质量为$[\rho u](x+\frac{\delta x}{2})$，因此质量净流入量为
$$
([\rho u](x-\frac{\delta x}{2})-[\rho u](x+\frac{\delta x}{2}))\delta y\delta z=-\frac{\partial (\rho u)}{\partial x}\delta x\delta y\delta z
$$
同理可知，这一小块各个面的质量净流入量总和为$-\nabla\cdot(\rho \mathbf{U})\delta x\delta y\delta z$，这等于小块内部质量的增加量$\delta x\delta y\delta z\frac{\partial \rho}{\partial t}$，于是有连续性方程
$$
\frac{\partial \rho}{\partial t}+\nabla\cdot(\rho\mathbf{U})=0
$$
利用全导数公式和矢量分析公式变形可化为
$$
\frac{1}{\rho}\frac{D\rho}{Dt}+\nabla\cdot\mathbf{U}=0
$$
### 1.4.2 Lagrange描述推导
上面的最后一个等式的物理意义由Lagrange描述解释，考虑一个随流体运动的小块$\delta x\delta y\delta z$，其质量$\delta M=\rho\delta x\delta y\delta z$应该守恒
$$
\frac{1}{\delta M}\frac{D}{Dt}(\delta M)=\frac{1}{\rho\delta V}\frac{D}{Dt}(\rho\delta V)=\frac{1}{\rho}\frac{D\rho}{Dt}+\frac{1}{\delta V}\frac{D}{Dt}(\delta V)=0
$$
对于第二项，有
$$
\frac{1}{\delta V}\frac{D}{Dt}(\delta V)=\frac{1}{\delta x}\frac{D}{Dt}(\delta x)+\frac{1}{\delta y}\frac{D}{Dt}(\delta y)+\frac{1}{\delta z}\frac{D}{Dt}(\delta z)
$$
而$\frac{D}{Dt}(\delta x)$就是小块$x$和$x+\delta x$两端的流速之差$\delta u$，因此$\frac{1}{\delta x}\frac{D}{Dt}(\delta x)=\frac{\partial u}{\partial x}$，对于y，z同理
$$
\frac{1}{\delta V}\frac{D}{Dt}(\delta V)=\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}+\frac{\partial w}{\partial z}=\nabla\cdot\mathbf{U}
$$
由此推出了一样的连续性方程
### 1.4.3 连续性方程的尺度分析
将连续性方程写开来，代入$\rho(x,y,z,t)=\rho_0(z)+\rho'(x,y,z,t)$
$$
\frac{1}{\rho_0}\left(\frac{\partial\rho^{\prime}}{\partial t}+\mathbf{U}\cdot\nabla\rho^{\prime}\right)+\frac{w}{\rho_0}\frac{d\rho_0}{dz}+\nabla\cdot\mathbf{U}=0
$$
用和前面类似的方法进行数量级分析可得，第一项的数量级要小于后面两项，因此第一项可以忽略，得连续性方程近似形式
$$
\nabla\cdot(\rho_0\mathbf{U})=0
$$
注意区分其和不可压缩流体对应的$\nabla\cdot \mathbf{U}=0$，只有在纯粹的水平流时才可以当作是不可压缩流体，一旦有了垂直流动，$\rho_0$随z的变化就必须纳入考虑。
## 1.5 热力学能量方程
热力学能量方程是描述能量守恒的方程，来自于热力学第一定律，不过由于流体是处于运动状态而非静止的，因此其是否还满足热力学第一定律的表达式值得思考。

考虑一个随流体运动的小块$\delta V=\delta x\delta y\delta z$，其能量包括动能和内能，设e为单位质量的内能，则总能量为$\rho(e+\frac{1}{2}\mathbf{U}\cdot\mathbf{U})\delta V$。对流体做功的外力主要有两种，一个是压强，一个是重力。

对于压强，以$\delta y\delta z$面为例，考虑压强对流体做功的功率，应为
$$
([pu](x-\frac{\delta x}{2})-[pu](x+\frac{\delta x}{2}))\delta y\delta z=-\frac{\partial (pu)}{\partial x}\delta V
$$
对另外两组面也是同理，于是可得压强对流体做功的功率为
$$
-\nabla\cdot(p\mathbf{U})\delta V
$$
重力对流体做功的功率容易计算，
$$
\rho\mathbf{g}\cdot\mathbf{U}\delta V
$$
设J为单位质量的吸热速率，于是由能量守恒可得
$$
\frac{D}{Dt}\left[\rho\left(e+\frac{1}{2}\mathbf{U}\cdot\mathbf{U}\right)\delta V\right]=-\nabla\cdot(p\mathbf{U})\delta V+\rho\mathbf{g}\cdot\mathbf{U}\delta V+\rho J\delta V
$$
展开
$$
\rho\delta V\frac{D}{Dt}\left(e+\frac{1}{2}\mathbf{U}\cdot\mathbf{U}\right)+\left(e+\frac{1}{2}\mathbf{U}\cdot\mathbf{U}\right)\frac{D\left(\rho\delta V\right)}{Dt}=-\mathbf{U}\cdot\nabla p\delta V-p\nabla\cdot\mathbf{U}\delta V-\rho gw\delta V+\rho J\delta V
$$
由质量守恒，$\frac{D(\rho\delta V)}{Dt}=0$，又由忽略了摩擦力的运动方程$\frac{D\mathbf{U}}{Dt}+2\mathbf{\Omega}\times\mathbf{U}+\frac{\nabla p }{\rho}+g\mathbf{k}=0$，用$\rho\mathbf{U}$在等式左边点乘，可得
$$
\rho\frac{D}{Dt}\left(\frac{1}{2}\mathbf{U}\cdot\mathbf{U}\right)=-\mathbf{U}\cdot\nabla p-\rho gw\quad(*)
$$
将上式与能量守恒的方程相结合可得
$$
\rho\frac{De}{Dt}=-p\nabla\cdot\mathbf{U}+\rho J\quad(**)
$$
上面$(*)$式如果把$gw$写成$\frac{D(gz)}{Dt}=\frac{D\Phi}{Dt}$，则$(*)$可写成
$$
\rho\frac{D}{Dt}\left(\frac{1}{2}\mathbf{U}\cdot\mathbf{U}+\Phi\right)=-\mathbf{U}\cdot\nabla p
$$
显然对应于机械能方程，可以理解为对一个质点写能量守恒方程

而对于$(**)$式，可以稍作改写，由连续性方程$\frac{1}{\rho}\frac{D\rho}{Dt}+\nabla\cdot\mathbf{U}=0$，故$\frac{1}{\rho}\nabla\cdot\mathbf{U}=-\frac{1}{\rho^2}\frac{D\rho}{Dt}=\frac{D\alpha}{Dt}$，其中$\alpha=\frac{1}{\rho}$为单位质量的体积，代入$(*)$式可得
$$
c_v\frac{DT}{Dt}+p\frac{D\alpha}{Dt}=J
$$
这就是单位质量下的热力学第一定律，可见对于流体，热一定律仍然成立。
## 1.6 干空气热力学
由理想气体状态方程，$p\alpha=RT$，两边对时间求导数得
$$
p\frac{D\alpha}{Dt}+\alpha\frac{Dp}{Dt}=R\frac{DT}{Dt}
$$
再定义单位质量等压摩尔热容$c_p=c_v+R$，则热一定律可写为
$$
c_p\frac{DT}{Dt}-\alpha\frac{Dp}{Dt}=J
$$
两边同时除以T可得
$$
c_p\frac{D\ln T}{Dt}-R\frac{D\ln p}{Dt}=\frac{J}{T}\equiv\frac{Ds}{Dt}
$$
对于一个可逆过程上式对应于单位质量熵的变化率
### 1.6.1 干空气绝热过程
即J=0，可得
$$
c_pD\ln T-RD\ln p=D\left(c_p\ln T-R\ln p\right)=0
$$
这对应于一个干空气绝热过程下的守恒量，使其具有温度的量纲，称为位温（potential temperature）
$$
\theta=T\left(p_s/p\right)^{R/c_p}
$$
对于天气尺度的运动，只要在降雨活跃区域以外，就可以近似认为$\theta$守恒

在大气分析中，位温$\theta$是很重要的物理量，我们可以用$\theta$重新表述热一定律
$$
c_p\frac{D\ln\theta}{Dt}=\frac{J}{T}=\frac{Ds}{Dt}
$$
下面考虑干空气在绝热过程下的温度衰减率，对$\theta$先取对数，再对z求导数，结合流体静力平衡方程，可得
$$
					\frac{T}{\theta} \frac{ \partial \theta }{ \partial z } =\frac{ \partial T }{ \partial z } +\frac{g}{c_{p}}
$$
在绝热情况下，$\theta$守恒，故有
$$
-\frac{dT}{dz}=\frac{g}{c_{p}}=\Gamma_{d}
$$
可见干空气绝热的温度衰减率近似为常数$\Gamma_{d}$
### 1.6.2 稳定性
考虑干空气中处于平衡的一小块气体，其温度和压强均与环境一致，设在扰动下该小块发生了沿z轴正向的$\delta z$的位移，且压强与环境始终保持一致，则z方向运动方程为
$$
\frac{Dw}{Dt}=\frac{D^2}{Dt^2}(\delta z)=-g-\frac{1}{\rho}\frac{ \partial p_{E} }{ \partial  z} =g\left(\frac{ \rho_{E}-\rho}{\rho} \right)
$$
结合理想气体方程$\rho=\frac{p}{RT}$
$$
\frac{D^2}{Dt^2}(\delta z)=g\left( \frac{T-T_{E}}{T_{E}} \right)
$$
从温度的角度来看稳定性：
$$
T(\delta z)-T_{E}(0)=-\Gamma_{d}\delta z
$$
$$
T_{E}(\delta z)-T_{E}(0)=\frac{dT_{E}}{dz}\delta z
$$
故有
$$
T-T_{E}=-\left( \frac{dT_{E}}{dz}+\Gamma_{d} \right)\delta z
$$
因此当$\frac{dT_{E}}{dz}>-\Gamma_{d}$时，稳定；当$\frac{dT_{E}}{dz}>-\Gamma_{d}$时，不稳定

也可以从位温$\theta$的角度来看稳定性：
$$
\frac{D^2}{Dt^2}(\delta z)=g\left( \frac{\theta-\theta_{E}}{\theta_{E}} \right)
$$
$$
\theta(z)=\theta_{E}(0)
$$
$$
\theta_{E}(z)-\theta_{E}(0)=\frac{d\theta_{E}}{dz}\delta z
$$
故有
$$
\frac{D^2}{Dt^2}(\delta z)=-\frac{g}{\theta_{E}}\frac{d\theta_{E}}{dz}\delta z=-g\frac{d\ln\theta_{0}}{dz}\delta z=-N^2\delta z
$$
可见当$N^2>0$时，稳定，此时发生的是buoyancy oscillation，N是buoyancy frequency；$N^2<0$时，不稳定，用$\theta$表述见下表
![[干空气稳定性.png]]
在天气尺度下，大气一般都是稳定的，因为不稳定会导致发生对流从而回归稳定
### 1.6.3 热力学能量方程的尺度分析
估算数量级时，对于波动相比于本身值要小得多的物理量，将波动项分解出来是一种常用的方法（在前面的尺度分析中都有使用），比如对于位温除了z方向的波动以外，其他方向的波动都远小于其本身，因此可以将位温分为$\theta=\theta_{0}(z)+\theta'(x,y,z,t)$，借助$|\theta/\theta_0|\ll1,|d\theta/dz|\ll d\theta_0/dz$做近似，代入热一定律得
$$
\frac{1}{\theta_0}\left(\frac{\partial\theta'}{\partial t}+u\frac{\partial\theta'}{\partial x}+v\frac{\partial\theta'}{\partial y}\right)+w\frac{d\ln\theta_0}{dz}=\frac{J}{c_pT}
$$
对于对流层，热辐射比较弱，因此$\frac{J}{c_{p}}$的数量级相比于其他项要小，可近似写成
$$
\left(\frac{\partial\theta'}{\partial t}+u\frac{\partial\theta'}{\partial x}+v\frac{\partial\theta'}{\partial y}\right)+w\frac{d\theta_0}{dz}\approx0
$$
如果把温度场也写成 $T(x,y,z,t)=T_{0}(z)+T'(x,y,z,t)$，对于相同的z，近似认为压强相等，则可得$\frac{\theta'}{\theta_{0}}\approx \frac{T'}{T_{0}}$，令$-\frac{dT_{0}}{dz}=\Gamma$，代入上面的方程可得关于温度场的近似热力学能量方程
$$
\left(\frac{\partial T'}{\partial t}+u\frac{\partial T'}{\partial x}+v\frac{\partial T'}{\partial y}\right)+w\left(\Gamma_d-\Gamma\right)\approx0
$$
## 1.7 Boussinesq 近似
当空气块垂直位移较小时，其密度的改变也就较小，这时运动方程可以近似，将$\rho$换成常数$\rho_{0}$（除了在垂直方程中），这种近似就叫Boussinesq近似，在这种近似下运动方程的分量形式如下
$$
\frac{Du}{Dt}=-\frac{1}{\rho_0}\frac{\partial p}{\partial x}+fv+F_{rx}
$$
$$
\frac{Dv}{Dt}=-\frac{1}{\rho_0}\frac{\partial p}{\partial y}-fu+F_{ry}
$$
$$
\frac{Dw}{Dt}=g\frac{\theta'}{\theta_{0}}+F_{rz}
$$
绝热方程写作
$$
c_{p} \frac{D\ln\theta}{Dt}=0\Rightarrow \frac{D\theta_{0}}{Dt}+\frac{D\theta'}{Dt}=0\Rightarrow \frac{D\theta'}{Dt}+w \frac{d\theta_{0}}{dz}=0
$$
Boussinesq近似下，连续性方程写作
$$
									\nabla\cdot\mathbf{U}=0\Rightarrow \frac{ \partial u }{ \partial x } +\frac{ \partial v }{ \partial y } +\frac{ \partial w }{ \partial z } =0
$$
## 1.8 湿空气热力学
空气中的水汽（water vapor）对于大气动力学产生的影响主要源于其可能发生相变从而带来相变潜热，改变空气团的动力学（热一定律变了）

对于干空气，理想气体状态方程为
$$
p_{d}=\rho_{d}R_{d}T
$$
对于纯的水汽，理想气体状态方程为
$$
e=\rho_{v}R_{v}T
$$
对于干空气和水汽混合形成的湿空气，理想气体状态方程可以写为
$$
p=\rho R_{d}T_{v}
$$
其中$T_{v}$是virtual temperature，指具有相同p和$\rho$的干空气需要具有的温度，计算可得
$$
T_v=\frac{T}{1-\frac{e}{p}\left(1-\frac{R_d}{R_v}\right)}>T
$$
可见如果把湿空气近似为干空气会造成温度出现偏差，而这对于稳定与否的判断很重要。

水汽在空气中的饱和蒸气压由克劳修斯-克拉伯龙方程给出
$$
\frac{1}{e_s}\frac{de_s}{dT}=\frac{L}{R_vT^2}
$$
$$
e_{s}=A\exp\left( -\frac{L}{R_{v}T} \right)
$$
相对湿度定义为
$$
RH=\frac{e}{e_{s}}\times100\%
$$
混合比（mixing ratio）定义为一块空气中水汽质量与空气质量之比
$$
q=\frac{\rho_{v}}{\rho}=\frac{\mu}{\mu_{v}} \frac{e}{p}
$$
饱和混合比 $q_{s}=\frac{\mu}{\mu_{v}} \frac{e_{s}}{p}$
### 1.8.1 相当位温
当相对湿度小于100%时，湿空气和干空气的热一定律没有差别，从而他们绝热过程中的守恒量以及大气稳定的条件（$\frac{ \partial \theta }{ \partial z }>0$）。但当湿空气上升时，温度下降，饱和蒸气压也下降，最终水汽达到饱和，从这个高度再向上，水汽就会开始凝结，此高度称为lifting condensation level（LCL）。当大气不稳定时，空气团达到一定高度再向上时，会受到正的浮力，此高度称为level of free convection（LFC）。

讨论湿空气动力学时，引入相当位温（Equivalent potential temperature）这一新的守恒量会带来方便。

相当位温$\theta_{e}$的定义是一块空气所有水汽全都凝结、相变潜热全都凝结放出之后对应的位温，即当$q=0$时，$\theta_{e}=\theta$。要想得到一块气体的相当位温，可以理解为需要先让气体上升，待所有水汽全都凝结以后，相当位温就是此时$q_{s}\approx 0$的位温。显然，只要是同一块气体，在其运动过程中$\theta_{e}$守恒。相当位温的物理意义在于，其相当于把一块湿空气里所有的潜在能量都核算在了最终的位温里。（类似的，位温的物理意义就是把干空气的所有能量都核算在位温里）

$\theta_{e}$作为新守恒量的推导如下：

先考虑空气中的水汽已经饱和的情况，此时混合比为饱和混合比$q_{s}$，单位质量空气块单位时间吸热为相变潜热
$$
J=-L_{c} \frac{Dq_{s}}{Dt}
$$
因此热一定律写作
$$
c_{p} \frac{D\ln\theta}{Dt}=- \frac{L_{c}}{T} \frac{Dq_{s}}{Dt}
$$
由于$q_{s}$的随空气块运动的变化速率远大于T的变化速率，因此可以做近似
$$
d\ln\theta\approx-d\left(L_cq_s/c_pT\right)
$$
利用$q_{s}\approx 0$时，$\theta_{e}=\theta$，假设T的变化可以忽略，即保持为最后水汽全都凝结时的温度T，则对上面的微分表达式积分可得
$$
\ln\left(\theta/\theta_e\right)\approx-L_cq_s/c_pT
$$
$$
\theta_e\approx\theta\exp\left(L_cq_s/c_pT\right)
$$
这就是一块饱和湿空气的相当位温。对于不饱和的湿空气来说，将$q_{s}$换成$q$，T换成空气块绝热达到LCL处时的温度$T_{LCL}$，则由于未发生凝结，在运动过程中$q$始终不变，且干空气绝热过程$\theta$守恒，因此这样定义的相当位温无论是在不饱和的绝热位移还是饱和的绝热位移下，都守恒，因此$\theta_{e}$是追踪湿空气块变化的利器。记J为除相变潜热外额外的吸热
$$
J=c_{p}T \frac{d\ln\theta}{dt}+L_{c} \frac{dq_{s}}{dt}=c_{p}T \frac{d\ln\theta_{e}}{dt}
$$
$$
J=c_{p} \frac{dT}{dt}-\frac{1}{\rho} \frac{dp}{dt}+L_{c} \frac{dq_{s}}{dt}=c_{p} \frac{dT}{dt}+g \frac{dz}{dt}+L_{c} \frac{dq_{s}}{dt}=\frac{dh}{dt}
$$
所以有
$$
c_{p}Td\ln\theta_{e}=dh
$$
这就是相当位温物理意义的体现。（位温物理意义的体现也对应类似的等式）
### 1.8.2 饱和绝热衰减率（Pseudoadiabatic Lapse Rate）
根据上文提到的方程，在考虑相变潜热情况下的热一定律写作
$$
\frac{dT}{dz}+\frac{g}{c_{p}}=- \frac{L_{c}}{c_{p}} \frac{dq_{s}}{dz}
$$
由于 $q_{s}=q_{s}(T,p)$，利用流体静力平衡方程
$$
									\frac{dq_{s}}{dz}=\left( \frac{ \partial q_{s} }{ \partial T }  \right)_{p} \frac{dT}{dz} - \left( \frac{ \partial q_{s} }{ \partial p }  \right)_{T}\rho g
$$
定义$\epsilon=\frac{\mu_{v}}{\mu}$，代入$e_{s}$的表达式可得
$$ 
			\left( \frac{ \partial q_{s} }{ \partial T }  \right)_{p}=\frac{L_{c}}{R_{v}T^2}q_{s}=\frac{\epsilon L_{c}}{RT^2}q_{s}
$$
$$
			\left( \frac{ \partial q_{s} }{ \partial p }  \right)_{T}=- \frac{q_{s}}{p}
$$
最终可得饱和湿空气绝热位移时的温度衰减率
$$
\Gamma_s\equiv-\frac{dT}{dz}=\Gamma_d\frac{[1+L_cq_s/(RT)]}{\left[1+\varepsilon L_c^2q_s/\left(c_pRT^2\right)\right]}<\Gamma_{d}
$$
### 1.8.3 湿空气稳定性
用和干空气稳定性类似的研究方法
$$
\frac{D^2}{Dt^2}(\delta z)=g\left( \frac{T-T_{E}}{T_{E}} \right)
$$
$$
T-T_{E}=-\left( \frac{dT_{E}}{dz}+\Gamma_{s} \right)\delta z
$$
因此从温度角度来看，当$\frac{dT_{E}}{dz}>-\Gamma_{s}$时，稳定；当$\frac{dT_{E}}{dz}<-\Gamma_{s}$时，不稳定。特别的，当$\Gamma_{s}<-\frac{dT_{E}}{dz}=\Gamma<\Gamma_{d}$时，对于干绝热位移稳定，但对于湿绝热位移不稳定，这种情况称为条件不稳定性（Conditional Instability）。

从相当位温的角度来看
$$
\frac{D^2}{Dt^2}(\delta z)=g\left( \frac{\theta-\theta_{E}}{\theta_{E}} \right)
$$
此时随着气块的上升，不仅环境的位温改变
$$
\theta_{E}(\delta z)-\theta_{E}(0)=\frac{d\theta_{E}}{dz}\delta z
$$
气块自身的位温也因为相变潜热而改变
$$
\theta(\delta z)-\theta_{E}(0)=\delta \theta
$$
$$
\frac{\delta\theta}{\theta}\approx-\delta\left(\frac{L_cq_s}{c_pT}\right)\approx-\frac{d}{d z}\left(\frac{L_cq_s}{c_pT}\right)_{E}\delta z
$$
注意上面的第二步做了个近似，将气块本身的$\frac{L_{c}q_{s}}{c_{p}T}$的变化等价于了环境的变化，严格来说这是有问题的，除非认为气块和环境保持等温（因此我觉得Holton这里的推导存在问题，有无大佬给个解释）

由于考察的气块总是饱和的，但周围的空气可能不是饱和的，此时二者就无法直接比较，因此定义一个新变量$\theta_{e}^*$为将环境大气在保持热力学结构不变（T,p不变）的情况下人为地看作饱和（即规定q取$q_{s}$）后对应的相当位温（$\theta_{e}^*$只有在大气饱和时才与$\theta_{e}$相等），此时二者才可以比较。
$$
\frac{\theta-\theta_{E}}{\theta_{E}}=-\left[\frac{1}{\theta_{E}}\frac{d\theta_{E}}{d z}+\frac{d}{d z}\left(\frac{L_cq_s}{c_pT}\right)_{E}\right]\delta z\approx-\left(\frac{d\ln\theta_e^*}{d z}\right)_{E}\delta z
$$
故当$\left(\frac{d\theta_{e}^*}{dz}\right)>0$时，稳定；当$\left(\frac{d\theta_{e}^*}{dz}\right)<0$时，不稳定

以上的讨论个人觉得略显繁琐而且存在一些问题，对于这个判据我的想法是：

在引入了$\theta_{e}^*$之后，由于$\exp\left(L_cq_s/c_pT\right)$受T的不同影响相比系数上T要小很多，因此可以直接近似认为
$$
\frac{\theta-\theta_{E}}{\theta_{E}}\approx\frac{\theta_{e}-\theta_{e}^*}{\theta_{e}^*}\quad(\#)
$$
而类似于位温，也有$\theta_{e}-\theta_{e}^*=-\frac{d\theta_{e}^*}{dz}\delta z$，因此同样可得
$$
\frac{D^2}{Dt^2}(\delta z)=-g \frac{d\ln\theta_{e}^*}{dz}\delta z
$$
这也是为什么我们常说判断湿空气稳定性时，只需要比较气块的相当位温与环境强制饱和后的相当位温的大小关系即可。(环境要强制饱和的目的是使(#)式近似成立，即两个相当位温可比)
![[相当位温截面.png]]
上图给出了一个conditionally unstable的例子（$\frac{d\theta_{e}}{dz}<0$且$\frac{d\theta}{dz}>0$），虚线代表的是一个湿空气块的相当位温随上升的变化，而虚线与$\theta_{e}^*$的交点代表的就是即将形成正浮力的临界点，即LFC。在LCL到LFC之间尽管这个饱和气块在环境中的确不稳定，但是不会发生自发的对流，需要外力将该气块抬升至LFC处以后，不稳定性释放，从而发生强对流。（因此在LFC以下的对流多为外力迫使的弱对流）

从图中我们可以发现，从地表附近抬升的气团的$\theta_{e}$会在850hPa附近与环境的$\theta_{e}^*$相交，但该气团若从850hPa开始抬升，则不会在任何地方与环境的$\theta_{e}^*$相交，因此海洋上空的对流通常需要底层气团汇聚抬升才能发生，而大陆地区的对流有时不需要底层气团汇聚，是因为地表附近的强热导致气团获得了正浮力（即LFC较低），不过持续的深度强对流仍需要底层气团汇聚。






