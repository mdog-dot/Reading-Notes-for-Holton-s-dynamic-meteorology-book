---
zhihu-title: Holton读书笔记-Chapter1 基本介绍
zhihu-topics: 大气科学
zhihu-link: https://zhuanlan.zhihu.com/p/1928544934622924935
tags:
  - 大气科学
  - 动力学
  - 流体力学
zhihu-created-at: 2025-07-17 18:00
zhihu-updated-at: 2025-07-28 23:11
---
```toc
```
# 1 前言
这是我参考大气动力学著名的教科书——Holton所写的《An Introduction to Dynamic Meteorology》读书笔记，注意这**不是**抄书式笔记，调整了部分内容的顺序，删去了书里我觉得没什么意思的一些内容，并对我觉得有意思的地方做了更详尽的解释。（在部分地方会参考Marshall的《Atmosphere，Ocean and Climate》和舒幼生的《力学》）

# 2 Chapter1
第一章主要是对于大气动力学的各种概念和大致思想做一个概括性地介绍。当考察大气中一个小块的运动时，往往采用局域坐标系（由于地球曲率很小，该坐标架往往可以看作固定的），即速度 $\vec{U}=(u,v,w)$ ，其中u为向东方向，v为向北方向，w为垂直向上方向。
## 2.1 大气中的力
考虑大气中一个小块，其受到的力分为两类：体力（body forces）和表面力（surface forces），体力包括重力以及旋转系导致的力（apprent forces），表面力包括气压梯度力、粘滞阻力，下面我们将一个一个讨论。
### 2.1.1 气压梯度力
![[Pressure gradient force.png]]
考虑小块沿x轴方向的气压造成的受力，
$$
F_{px}=p(x-\frac{\delta x}{2},y,z)\delta y\delta z-p(x+\frac{\delta x}{2},y,z)\delta y\delta z=-\frac{\partial p}{\partial x}\delta x\delta y\delta z
$$
同理，有
$$
\vec{F}_p=-\nabla p\delta x\delta y\delta z
$$
因此气压梯度造成的加速度为
$$
\frac{\vec{F}_p}{m}=-\frac{\nabla p}{\rho}
$$
### 2.1.2 重力
重力造成的加速度为
$$
\frac{F_g}{m}=\vec{g}^*=-\frac{GM}{r^2}\hat{r}
$$
如果将地球当成半径为a的球形，$r=a+z$ ，则 $\vec{g}^*=\frac{\vec{g}_0^*}{(1+z/a)^2}$ ，$\vec{g}_0^*$是海平面处的重力加速度，一般来说$z<<a$，因此往往认为$\vec{g}^*=\vec{g}_0^*$ 
### 2.1.3 旋转系导致的力
大气动力学的本质就是旋转系中的流体力学，地球自转这一性质对于分析大气受力有重要影响，因此对旋转系这一非惯性系的考虑是必要的。
考虑任一矢量 $\vec{A}$ 在旋转系和惯性系中的时间导数，设在惯性系中坐标架为 $\hat{x},\hat{y},\hat{z}$，则有
$$(\frac{d\vec{A}}{dt})_{in}=\hat{x}\frac{d\vec{A_x}}{dt}+\hat{y}\frac{d\vec{A_y}}{dt}+\hat{z}\frac{d\vec{A_z}}{dt}$$
而在旋转系中，惯性系里的坐标架会发生转动，有
$$(\frac{d\hat{x}}{dt})_{rot}=-\Omega\times\hat{x}$$
对于$\hat{y},\hat{z}$同理，故有
$$
\begin{aligned}
(\frac{d\vec{A}}{dt})_{rot}&=(\frac{d\vec{A}}{dt})_{in}+(\frac{d\hat{x}}{dt})_{rot}A_x+(\frac{d\hat{y}}{dt})_{rot}A_y+(\frac{d\hat{z}}{dt})_{rot}A_z\\&=(\frac{d\vec{A}}{dt})_{in}-\Omega\times\vec{A}
\end{aligned}
$$
取 $\vec{A}=\vec{u}$ ，即小块的速度，则有
$$
\begin{aligned}
(\frac{d\vec{u}_{in}}{dt})_{in}&=((\frac{d}{dt})_{rot}+\Omega\times)(\vec{u}_{rot}+\Omega\times r)\\&=(\frac{d\vec{u}_{rot}}{dt})_{rot}+\Omega\times\Omega\times r+2\Omega\times\vec{u}_{rot}
\end{aligned}
$$
$$\vec{F}_{rot}=\vec{F}_{in}-m\Omega\times\Omega\times r-2m\Omega\times\vec{u}_{rot}$$
由此得到了旋转系导致的力：惯性离心力 $-m\Omega\times\Omega\times r$ 和科里奥利力 $-2m\Omega\times\vec{u}$ 
#### 2.1.3.1 惯性离心力和重力再讨论
对于在地球自转系下静止的一个小块，在重力和惯性离心力的作用下保持平衡。
设 $\vec{R}$ 为从轴线指向小块位置的矢量，则惯性离心力造成的加速度为 $\Omega^2\vec{R}$ ，要使 $g^*$ 和$\Omega^2 \vec{R}$ 被支持力相平衡，则地球必须偏离球形，表面与 $\vec{g}=\vec{g}^*+\Omega^2\vec{R}$ 垂直，如下图所示
![[New gravity.png]]
因此真正的重力加速度其实是上面引出的$\vec{g}=g\hat{z}$，由引力和惯性离心力合并而成，后面的讨论中将不会再分开讨论。
### 2.1.4 粘滞阻力
先考虑向东方向的速度u（对应x方向）造成的粘滞阻力，由粘滞阻力公式，两层流体在y方向之间单位面积的摩擦力为 $f_x(y)=\mu \frac{\partial u}{\partial y}$ ，因此一块流体 $\delta x\delta y\delta z$ 在y方向之间摩擦造成的粘滞阻力为 $\mu \frac{\partial^2 u}{\partial y^2}\delta x \delta y\delta z$ ，同理在z方向之间摩擦造成的粘滞阻力为 $\mu \frac{\partial ^2 u}{\partial z^2}\delta x\delta y\delta z$ 

Holton告诉我沿着x方向（即速度方向）也存在着粘滞阻力，这是我第一次见到，但是仔细想想确实有道理。考虑粘滞阻力的本质，其实是粒子输运造成的动量交换，由于粒子的动量不同造成的力。

考虑x处的动量交换，单位时间x右侧增加的动量为 $\frac{1}{6}\rho \bar{v}(u(x-\bar{\lambda})-u(x+\bar{\lambda}))A=-\frac{1}{3}\rho\bar{v}\bar{\lambda}\frac{\partial u}{\partial x}A$ ，因此小块因为x方向粒子输运受到的力为 $\frac{1}{3}\rho\bar{v}\bar{\lambda} \frac{\partial ^2 u}{\partial x^2}\delta x\delta y\delta z$ ，类似的对y和z方向进行分析可得粘滞系数 $\mu=\frac{1}{3}\rho\bar{v}\bar{\lambda}$ (任意一本热学书都有详细的讨论)。

综上，u造成的粘滞阻力均沿x方向，其造成的加速度为
$$
F_{rx}=\frac{\mu}{\rho}(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 u}{\partial y^2}+\frac{\partial^2 u}{\partial z^2})=\frac{\mu}{\rho}\nabla^2u
$$
类似的可以求得 $F_{ry}、F_{rz}$
## 2.2 大气静态结构
### 2.2.1 流体静力学
如果垂直方向流体平衡，则有
$$
\frac{dp}{dz}=-\rho g
$$
$$
p(z)=\int_z^{\infty}\rho gdz
$$
即一块地方的压强等于其上所有空气的重力之和/面积，如果已知大气总质量为M，地球半径为a，则地球表面平均压强为
$$
\bar{p}=\frac{Mg}{4\pi a^2}
$$
利用理想气体状态方程 $p=\rho RT$ （注意在大气科学中的R往往指热学中的 $R/\mu$ ，是空气独有的常量），可得
$$
gdz=-\frac{RT}{p}dp=-RTdlnp
$$
$p_1$和$p_2$等压线之间的空气厚度为
$$
Z_T=z_2-z_1=\frac{R}{g}\int_{p_1}^{p_2}Tdlnp
$$
定义平均温度
$$
\langle T \rangle= \frac{\int_{p_{1}}^{p_{2}}Td\ln p}{\int_{p_{1}}^{p_{2}}d\ln p}
$$
则空气厚度可以写作
$$
Z_{T}=  \frac{R\left\langle T \right\rangle}{g} \ln \frac{p_{2}}{p_{1}}
$$
如果大气等温（等温大气模型），则 
$$Z=\frac{RT}{g}ln(\frac{p_0}{p})=Hln(\frac{p_0}{p})$$
$$
p(Z)=p(0)e^{-Z/H}
$$
### 2.2.2 压强坐标
在研究大气的受力方程时，气压梯度力一项 $-\frac{\nabla p}{\rho}$含有可变的 $\rho$ ，给未来微分方程的处理带来麻烦，而选用压强p代替z作为x，y之外的第三个坐标可以消去 $\rho$ 。

此时$z=z(x,y,p)$，在p坐标下，水平气压梯度力需要变成固定p情况下的偏导，因此我们重新考察 $(\frac{\partial p}{\partial x})_z$ 

借助恒等式，
$$
(\frac{\partial p}{\partial  x})_z(\frac{\partial x}{\partial z})_p(\frac{\partial z}{\partial p})_x=-1
$$
以及流体静力学平衡方程
$$
(\frac{\partial p}{\partial z})_x=-\rho g
$$
有
$$
-\frac{1}{\rho}(\frac{\partial p}{\partial x})_z=\frac{1}{\rho}(\frac{\partial z}{\partial x})_p(\frac{\partial p}{\partial z})_x=-g(\frac{\partial z}{\partial x})_p=-(\frac{\partial \Phi}{\partial x})_p
$$
故水平气压梯度力在p坐标下对应于 $-\nabla_p \Phi$ 的水平分量
### 2.2.3 一般坐标
除了选压强代替z作为坐标以外，还可以选其他物理量作为坐标，只要其为单值单调的函数，比如取$\sigma=p(x,y,z,t)/p_s(x,y,z,t)$ 为坐标，因此有必要讨论一般坐标下水平气压梯度力的表达式。
我们希望建立 $(\frac{\partial p}{\partial x})_z$ 和 $(\frac{\partial p}{\partial x})_s$ 的联系，从偏导数的定义出发，$p=p(x,y,z(x,y,s,t),t)$ 
$$
(\frac{\partial p}{\partial x})_z=\frac{p(x+\Delta x,y,z,t)-p(x,y,z,t)}{\Delta x}
$$
$$
\begin{aligned}
(\frac{\partial p}{\partial x})_s&=\frac{p(x+\Delta x,y,z(x+\Delta x,y,s,t),t)-p(x,y,z(x,y,s,t),t)}{\Delta x}\\&=\frac{p(x+\Delta x,y,z(x+\Delta x,y,s,t),t)-p(x,y,z(x+\Delta x,y,s,t),t)+p(x,y,z(x+\Delta x,y,s,t),t)-p(x,y,z(x,y,s,t),t)}{\Delta x}\\&=(\frac{\partial p}{\partial x})_z+\frac{\partial p}{\partial z}(\frac{\partial z}{\partial x})_s
\end{aligned}
$$
再利用链式法则 $\frac{\partial p}{\partial z}=\frac{\partial s}{\partial z}\frac{\partial p}{\partial s}$，可得在s坐标下的气压梯度
$$
(\frac{\partial p}{\partial x})_z=(\frac{\partial p}{\partial x})_s-\frac{\partial s}{\partial z}\frac{\partial p}{\partial s}(\frac{\partial z}{\partial x})_s
$$
## 2.3 运动学
考察流体在一点附近水平运动的结构，采用泰勒展开
$$
u(x_0+dx,y_0+dy)\approx u(x_0,y_0)+\frac{\partial u}{\partial x}dx+\frac{\partial u}{\partial y}dy
$$
$$
v(x_0+dx,y_0+dy)\approx v(x_0,y_0)+\frac{\partial v}{\partial x}dx+\frac{\partial v}{\partial y}dy
$$
流体旋转（vorticity）的度量 $\xi=\frac{\partial v}{\partial x}-\frac{\partial u}{\partial y}$ ，流体辐散（divergence）的度量 $\delta = \frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}$ ，流体的变形（deformation）的度量为 $d_1=\frac{\partial u}{\partial x}-\frac{\partial v}{\partial y}$和$\frac{\partial v}{\partial x}+\frac{\partial u}{\partial y}$，于是有
$$
u(x_0+dx,y_0+dy)\approx u(x_0,y_0)+\frac{1}{2}(\delta+d_1)dx+
\frac{1}{2}(d_2-\xi)dy
$$
$$
v(x_0+dx,y_0+dy)\approx v(x_0,y_0)+\frac{1}{2}(\xi+d_2)dx+\frac{1}{2}(\delta -d_1)dy
$$
vorticity,divergence,deformation的可视化图像如下图
![[可视化旋度散度.png]]

这三个是风场的重要属性，divergence在质量守恒方程中出现，vorticity在大气动力学中重要，deformation在创造和摧毁流体边界时重要（如前锋）。
## 2.4 尺度分析
对于大气的方程，分析各项的尺度——也就是数量级是十分重要的，因为往往存在几项是方程的主导项，而有几项可以忽略。
在估计数量级时，往往考察以下三种量的数量级：
1. 变量的大小
2. 变量波动的振幅（即变化的最大幅度）
3. 运动发生的特征长度、深度、时间

利用这三种数量级，可以估计方程中各项的数量级，从而找出可以忽略的项。比如对于大气的运动方程，设其速度和速度波动的数量级都为V，特征时间为T，特征长度L=VT，因此方程中加速度项的数量级为V/T，科里奥利力的数量级为$2\Omega V$ , 所以当运动的时间尺度较小（空间尺度较小）时，科里奥利力可以忽略，而当时间尺度（空间尺度）较大时，加速度项可以忽略，即近似为平衡。
![[空间尺度.png]]
上图给出了不同气象行为的空间尺度（水平尺度），由于大气的运动方程依赖于空间尺度，因此我们常常用空间尺度对系统分类。