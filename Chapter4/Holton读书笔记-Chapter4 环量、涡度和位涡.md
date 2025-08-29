---
zhihu-title: Holton读书笔记-Chapter4 环量、涡度和位涡
zhihu-topics: 大气科学
zhihu-link: https://zhuanlan.zhihu.com/p/1934567913379063792
zhihu-updated-at: 2025-08-06 09:22
toc: "true"
---
```toc
```
[[Holton读书笔记-Chapter1 基本介绍]](https://zhuanlan.zhihu.com/p/1928544934622924935 "card")
[[Holton读书笔记-Chapter2 基础守恒律]](https://zhuanlan.zhihu.com/p/1929242659588912307 "card")
[[Holton读书笔记-Chapter3 基本方程的初步应用]](https://zhuanlan.zhihu.com/p/1930180220087939862 "card")
# 1 Chapter4
对于经典的旋转刚体，其存在一角速度来描述其转动，相应的有角动量和角动量守恒律，而对于旋转流体来说，无法定义一个“角速度”，我们用环量（circulation）和涡度（vorticity）来描述其旋转。其中，环量是一个标量，描述的是一定区域流体旋转的宏观性质，涡度是一个矢量场，描述的是流体中一点的“微观”性质。而位涡扩展了涡度的概念，将热力学对运动的限制包含在了里面，从而构建了一个有效的解释大气动力学的框架。
## 1.1 环量定理

环量定义为，在给定闭合曲线后，
$$
C \equiv \oint \mathbf{U}\cdot d\mathbf{l}=\oint|\mathbf{U}|\cos \alpha dl
$$
$\mathbf{l}$为位置矢量，一般环路取为逆时针。

环量定理是通过对闭合曲线上的流体粒子的牛二定律沿闭合曲线环路积分得到的，先在绝对坐标系（下标a）下考虑
$$
\oint\frac{D_a\mathbf{U}_a}{Dt}\cdot d\mathbf{l}=-\oint\frac{\nabla p\cdot d\mathbf{l}}{\rho}-\oint\nabla\Phi\cdot d\mathbf{l}
$$
利用（对于标量来说不管在哪个参考系其随运动的变化速率均没有区别，因此直接写为 $\frac{D}{Dt}$）
$$
\begin{aligned}
\frac{D_{a}\mathbf{U}_{a}}{Dt}\cdot d\mathbf{l}&=\frac{D}{Dt}(\mathbf{U}_{a}\cdot d\mathbf{l})-\mathbf{U}_{a}\cdot\frac{D_{a}}{Dt}(d\mathbf{l})\\&=\frac{D}{Dt}(\mathbf{U}_{a}\cdot d\mathbf{l})-\mathbf{U_{a}}\cdot d\mathbf{U}_{a}
\end{aligned}
$$
$$
\oint \mathbf{U}_{a}\cdot d\mathbf{U}_{a}=\frac{1}{2}\oint d(\mathbf{U}_{a}\cdot \mathbf{U}_{a})=0
$$
以及
$$
\oint \nabla\Phi \cdot d\mathbf{l}=\oint d\Phi=0
$$
可得环量定理
$$
\frac{DC_a}{Dt}=\frac{D}{Dt}\oint\mathbf{U}_a\cdot d\mathbf{l}=-\oint\rho^{-1}dp
$$
方程右侧称为力管项，对于正压流体，$\rho=\rho(p)$，有
$$
\oint\rho^{-1}dp=0
$$
于是
$$
\frac{DC_{a}}{Dt}=0
$$
绝对环量守恒，称为Kelvin环量定理，对应于流体力学中的角动量守恒定律。

**讨论**：上式中所谓D/Dt指的是随什么运动的时间导数？我认为应该这样理解，初始取的闭合曲线上的每一个粒子都具有自己的速度，随时间变化这些粒子会运动到不同位置和速度，构成一个新的闭合曲线，而形成另一个位置和大小的环量$C_{a}$，这里的D/Dt就指的是这个过程中环量的变化率。

在研究气象学时，研究相对环量C更方便，绝对环量和相对环量的差来自于地球的自转，$C_{a}=C+C_{e}$
$$
C_{e}=\oint \mathbf{U}_{e}\cdot d\mathbf{l}=\int \int_{A}(\nabla\times \mathbf{U}_{e})\cdot \mathbf{n}dA
$$
$$
\mathbf{U}_{e}=\mathbf{\Omega}\times \mathbf{R}
$$
$$
\nabla \times \mathbf{U}_{e}=\nabla \times(\mathbf{\Omega} \times \mathbf{R})=\mathbf{\Omega}(\nabla \cdot \mathbf{R})=2\mathbf{\Omega}
$$
故有
$$
C_{e}=2\Omega \int \int_{A}\sin \phi dA=2\Omega A_{e}
$$
其中$A_{e}$是面积A往赤道面投影的面积，则相对环量为
$$
C=C_{a}-2\Omega A_{e}
$$
满足
$$
	\frac{DC}{Dt}=-\oint \frac{dp}{\rho}-2\Omega\frac{ DA_{e}}{Dt}
$$
这称为Bjerknes环量定理。

对于正压流体，可以将相对环量从初始积分到末态，
$$
C_{2}-C_{1}=-2\Omega(A_{2}\sin \phi_{2}-A_{1}\sin \phi_{1})
$$
流体粒子沿闭合曲线的相对环量会由于闭合曲线的面积和纬度的变化而变化。
## 1.2 涡度
绝对涡度$\mathbf{\omega}_{a}\equiv \nabla \times \mathbf{U}_{a}$，相对涡度$\mathbf{\omega}\equiv \nabla \times \mathbf{U}$，在笛卡尔坐标系里，相对涡度可以写成分量形式
$$
\mathbf{\omega}=\left(\frac{\partial w}{\partial y}-\frac{\partial v}{\partial z},\frac{\partial u}{\partial z}-\frac{\partial w}{\partial x},\frac{\partial v}{\partial x}-\frac{\partial u}{\partial y}\right)
$$
**尺度分析**：涡度水平分量尺度为$\frac{U}{H}$，垂直分量尺度为$\frac{U}{L}$，垂直分量与水平分量尺度之比为
$$
\frac{H}{L} \thicksim 10km/1000km=0.01
$$
可见涡度的主要分量为水平分量，然而后面会指出，对于大尺度气候动力学，我们主要关心的是绝对和相对涡度的垂直分量。
$$
\eta\equiv\mathbf{k}\cdot(\nabla\times\mathbf{U}_a),\quad\zeta\equiv\mathbf{k}\cdot(\nabla\times\mathbf{U})
$$
$\zeta$有其实际意义，取水平面内闭合曲线，沿逆时针方向为正，
$$
\begin{aligned}
\oint \mathbf{V}\cdot d\mathbf{l}&=\int \int_{A}(\nabla \times \mathbf{V})\cdot\mathbf{k}dA\\&=\int \int_{A}\mathbf{k}\cdot(\nabla \times \mathbf{U})dA\\&=\int \int_{A}\zeta dA
\end{aligned}
$$
所以当$\zeta>0$时，速度$\mathbf{V}$沿逆时针方向，若在北半球，则对应于气旋，在南半球气旋则对应于$\zeta<0$，因此气旋满足条件
$$
f\zeta>0
$$
绝对涡度和相对涡度之差为行星涡度（planetary vorticity），为
$$
\mathbf{k}\cdot \nabla \times \mathbf{U}_{e}=2\Omega \sin \phi=f
$$
因此有
$$
\eta=\zeta+f
$$
涡度和环量可以通过斯托克斯定理相联系
$$
\oint \mathbf{U}\cdot d\mathbf{l}=\int \int_{A}(\nabla \times \mathbf{U})\cdot \mathbf{n}dA
$$
取$\mathbf{U}$为水平速度$\mathbf{V}$，并取极限$dA\rightarrow 0$，则可得相对涡度的极限形式
$$
\zeta=\lim_{A\rightarrow 0}(\oint \mathbf{V}\cdot d\mathbf{l})A^{-1}
$$
涡度可以认为类似于局域的角速度，如果流体做刚体转动，涡度对应于两倍的角速度。
### 1.2.1 自然坐标系下的涡度
涡度垂直分量的物理直观含义在自然坐标系里看会比较方便。
![[自然坐标系下的涡度.png]]
上图为在自然坐标系下取的一个小环路，图中转角满足$d(\delta s)=\delta \beta \delta n$，则该环路的环量为
$$
\begin{aligned}
		\delta C&=V[\delta s+d(\delta s)]-\left( V+\frac{ \partial V }{ \partial n } \delta n \right)\delta s\\&=\left( V \frac{\delta \beta}{\delta s}-\frac{ \partial V }{ \partial n }  \right)\delta n\delta s
\end{aligned}
$$
取环路趋于无穷小，即$\delta n,\delta s\rightarrow 0$，则有 
$$
\begin{aligned}
				\zeta&=\lim_{ \delta n,\delta s \to 0 } \frac{\delta C}{\delta n\delta s}=V \frac{\delta \beta}{\delta s}-\frac{ \partial V }{ \partial n } \\&=\frac{V}{R_{s}}-\frac{ \partial V }{ \partial n } 
\end{aligned}
$$
其中$R_{s}$是流线的曲率半径。由上式可以看到，相对涡度的水平垂直分量包括两个部分，一是风速沿法向的变化率，称为切变涡度（shear vorticity），二是风沿流线的转向，称为曲率涡度（curvature vorticity），因此即使是直线流动的流体也可能有涡度，而曲线流动的流体也可能没有涡度（两项刚好抵消）。
## 1.3 涡度方程
上一部分讨论的是涡度的运动学性质，这一部分利用运动方程讨论涡度的动力学性质，即其时间变化率。
### 1.3.1 笛卡尔坐标系下的形式
对于天气尺度的运动，涡度方程可以用Chapter2中给出的近似运动方程推导得到。分别将x方向方程对y求偏导，y方向方程对x求偏导数，可得
$$
\frac{\partial}{\partial y}\left(\frac{\partial u}{\partial t}+u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}+w\frac{\partial u}{\partial z}-fv=-\frac{1}{\rho}\frac{\partial p}{\partial x}\right)
$$
$$
\frac{\partial}{\partial x}\left(\frac{\partial v}{\partial t}+u\frac{\partial v}{\partial x}+v\frac{\partial v}{\partial y}+w\frac{\partial v}{\partial z}+fu=-\frac{1}{\rho}\frac{\partial p}{\partial y}\right)
$$
下式减去上式，利用涡度定义$\zeta=\frac{ \partial v }{ \partial x }-\frac{ \partial u }{ \partial y }$，可得
$$
\begin{aligned}\frac{\partial\zeta}{\partial t}&+u\frac{\partial\zeta}{\partial x}+v\frac{\partial\zeta}{\partial y}+w\frac{\partial\zeta}{\partial z}+(\zeta+f)\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right)\\&+\left(\frac{\partial w}{\partial x}\frac{\partial v}{\partial z}-\frac{\partial w}{\partial y}\frac{\partial u}{\partial z}\right)+v\frac{df}{dy}=\frac{1}{\rho^2}\left(\frac{\partial\rho}{\partial x}\frac{\partial p}{\partial y}-\frac{\partial\rho}{\partial y}\frac{\partial p}{\partial x}\right)\end{aligned}
$$
由于$f$仅与y有关，故$\frac{Df}{Dt}=v \frac{df}{dy}$，于是上式可以写成更紧凑的形式
$$
\begin{aligned}\frac{D}{Dt}(\zeta+f)&=-(\zeta+f)\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right)\\&-\left(\frac{\partial w}{\partial x}\frac{\partial v}{\partial z}-\frac{\partial w}{\partial y}\frac{\partial u}{\partial z}\right)+\frac{1}{\rho^2}\left(\frac{\partial\rho}{\partial x}\frac{\partial p}{\partial y}-\frac{\partial\rho}{\partial y}\frac{\partial p}{\partial x}\right)\end{aligned}
$$
上面的方程表明，绝对涡度随运动的时间变化率由三项决定，第一项叫散度项（divergence），第二项叫倾斜项（tilting），第三项叫力管项（solenoidal）。

如果只考虑第一项，上式的意义类似于刚体中角动量守恒下转动惯量会变化的转动$I \frac{D\Omega}{Dt}+\Omega\frac{ DI}{Dt}$。当水平流是辐散时，一个由闭合曲线上的流体粒子围成的区域面积增加，如果绝对环量守恒，那么闭合曲线流体的绝对涡度会减少，与上面方程所呈现的相同。(书中将这种现象形象地称为涡度的稀释dilution，因为相当于环量随着辐散而摊开，因此涡度被稀释了，而辐合就对应于汇集concentration)

第二项从直观来看，可以理解为：水平速度沿垂直方向的变化会造成涡度，这个涡度所对应的涡线由于垂直速度沿水平方向的变化发生了倾斜，从而导致涡度的垂直分量增大。

第三项等价于环量定理中的力管项的微观（极限）形式，
$$
-\oint\alpha dp\equiv-\oint\alpha\nabla p\cdot d\mathbf{l}=-\int\int_A\nabla\times(\alpha\nabla p)\cdot\mathbf{k}dA
$$
利用$\nabla \times(\alpha \nabla p)=\nabla \alpha \times \nabla p$，可得
$$
-\oint\alpha dp=-\iint_A(\nabla\alpha\times\nabla p)\cdot\mathbf{k}dA
$$
而涡度方程中的力管项其实就是上式中的被积函数
$$
-\left(\frac{\partial\alpha}{\partial x}\frac{\partial p}{\partial y}-\frac{\partial\alpha}{\partial y}\frac{\partial p}{\partial x}\right)=-(\nabla\alpha\times\nabla p)\cdot\mathbf{k}
$$
### 1.3.2 压强坐标下的形式
涡度方程在压强坐标下形式更为简单，可以直接从运动方程的压强坐标形式出发，用矢量算符$\mathbf{k}\cdot \nabla \times$作用在方程两边即可(注意这里的$\nabla$指的是固定p后的偏导，即$\nabla_{p}$的简写)

首先利用矢量关系
$$
\begin{array}{c}(\mathbf{V}\cdot\mathbf{\nabla})\mathbf{V}=\mathbf{\nabla}\left(\frac{\mathbf{V}\cdot\mathbf{V}}{2}\right)+\mathbf{\zeta}\mathbf{k}\times\mathbf{V}\end{array}
$$
可重写压强坐标下的运动方程
$$
\frac{\partial\mathbf{V}}{\partial t}=-\nabla\left(\frac{\mathbf{V}\cdot\mathbf{V}}{2}+\Phi\right)-(\zeta+f)\mathbf{k}\times\mathbf{V}-\omega\frac{\partial\mathbf{V}}{\partial p}
$$
再用$\mathbf{k}\cdot \nabla \times$作用在等式两边
$$
\frac{\partial\zeta}{\partial t}=-\mathbf{V}\cdot\nabla\left(\zeta+f\right)-\omega\frac{\partial\zeta}{\partial p}-\left(\zeta+f\right)\nabla\cdot\mathbf{V}+\mathbf{k}\cdot\left(\frac{\partial\mathbf{V}}{\partial p}\times\nabla\omega\right)
$$
这就是压强坐标下的涡度方程，注意到该方程没有$p-\rho$的力管项，这是因为在压强坐标下，涡度的定义为$\zeta=(\partial v/\partial x-\partial u/\partial y)_{p}$，与笛卡尔坐标下固定z的偏导不同。在实际应用中影响不大，这是因为实际上力管项相当小，在天气尺度的运动中可以忽略。
### 1.3.3 涡度方程的尺度分析
各基本物理量尺度见下图
![[涡度方程基本量尺度.png]]
相对涡度量级
$$
\zeta=\frac{\partial v}{\partial x}-\frac{\partial u}{\partial y}\lesssim\frac{U}{L}\sim10^{-5}\mathrm{s}^{-1}
$$
$$
\zeta/f_0\lesssim U/\left(f_0L\right)\equiv\mathrm{Ro}\sim10^{-1}
$$
对于中纬度天气尺度系统，相对涡度相比行星涡度较小，因此有近似
$$
(\zeta+f)\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right)\approx f\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right)
$$
这个近似在强烈风暴中心不成立，那里$\zeta$和$f$大小相近。

涡度方程各项估计的尺度如下图所示
![[涡度方程各项尺度.png]]
最后三项均用了$\lesssim$符号，是因为它们都有可能出现两项相互抵消的情况，事实上散度项确实两项接近相互抵消，因为其本身的数量级大过其他项，因此天气尺度运动是准无散的。从尺度来看，散度项抵消后的尺度约为
$$
\left|\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right)\right|\lesssim10^{-6}\mathrm{~s}^{-1}
$$
仅保留涡度方程中量级在$10^{-10}s^{-2}$的项，得到天气尺度运动下的近似方程
$$
\frac{D_h(\zeta+f)}{Dt}=-f\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right)
$$
其中下标h表示水平的时间导数
$$
\frac{D_h}{Dt}\equiv\frac{\partial}{\partial t}+u\frac{\partial}{\partial x}+v\frac{\partial}{\partial y}
$$
在强烈风暴中则写为
$$
\frac{D_h(\zeta+f)}{Dt}=-(\zeta+f)\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right)
$$
上式解释了为何气旋比反气旋要强烈的多，气旋对应低压中心，散度小于0，即发生辐合，此时绝对涡度$\zeta+f$增加，于是方程的右侧项进一步增加，导致绝对涡度的增加速率越来越快；而对于反气旋，发生辐散，绝对涡旋将会减小直到变成0，即$\zeta\to -f$。

上面的近似在大气的前锋中不成立，此时涡度方程各项数量级接近。
## 1.4 位涡（Potential Vorticity）
重新回到Kelvin环量定理，这次我们将其应用的曲线做一个限制，限制曲线必须在一个等焓面上，由于$c_{p}d\ln \theta=ds$，因此等焓面也是等位温$\theta$面，对于这样的曲线，根据位温的定义有
$$
\rho=p^{c_v/c_p}(R\theta)^{-1}(p_s)^{R/c_p}
$$
可见在曲线上$\rho$只与p有关，因此力管项$\oint \frac{1}{\rho}dp=0$，因此对于绝热无摩擦的流动，初始在等位温面上的闭合曲线始终在该等位温面上，于是不管其是不是正压流体，都满足Kelvin环量定理，即
$$
\frac{DC_{a}}{Dt}=0
$$
![[位涡示意图.png]]
接下来考虑一个如上图所示的小圆柱体，其顶面和底面都限制在两个等位温面上，在绝热流动下两个面都始终保持在对应的等位温面上，且圆柱内流体粒子质量$dm$守恒。当圆柱取的极小时，将顶面取作闭合曲线，根据Kelvin环量定理以及涡度的极限定义
$$
\frac{DC_a}{Dt}\approx\frac{D}{Dt}[\mathbf{\omega}_{a}\cdot\mathbf{n}dA]=0
$$
$\mathbf{n}$指垂直于等位温面的方向，故
$$
\mathbf{n}=\frac{\nabla \theta}{\left|\nabla \theta\right|}
$$
小面积dA可以用守恒量$dm,d\theta$表示出来
$$
	dA=\frac{dm}{\rho dh}=\frac{dm}{\rho} \frac{\left|\nabla \theta\right|}{d\theta}
$$
由于$dm,d\theta$守恒，代入得
$$
\frac{D}{Dt}\left[\frac{ \mathbf{\omega}_{a}\cdot \nabla \theta}{\rho} \right]=0
$$
这就是**Ertel位涡定理**，是大气动力学中最重要的理论结果之一，它表明位涡（PV）
$$
\Pi=\frac{\boldsymbol{\omega}_a\boldsymbol{\cdot}\boldsymbol{\nabla}\theta}{\rho}
$$
在绝热无摩擦的流动中随运动守恒。位涡的重要意义在于把基础的守恒律囊括在了一个表达式中，对大气运动加上了热力学的约束。

Ertel位涡用例：若绝对涡度只有垂直分量，则PV为
$$
\frac{1}{\rho}\left(\frac{\partial v}{\partial x}-\frac{\partial u}{\partial y}+f\right)\frac{\partial\theta}{\partial z}
$$
如上面的位涡示意图所示，当圆柱往右移动时，$\frac{ \partial \theta }{ \partial z }$减小，由于PV守恒，其绝对涡度增大，又由于绝对环量守恒，因此绝对涡度增大表现为圆柱的面积减小。

### 1.4.1 尺度分析
PV的垂直分量尺度为$\frac{U\theta^*}{\rho HL}+\frac{f\theta^*}{\rho H}$，水平分量尺度为$\frac{U\theta^*}{\rho HL}$，因此垂直分量于水平分量的数量级之比为$1+R_{o}^{-1}$，对于天气或更大尺度，$R_{o}$小，故PV的主导项为其垂直分量，这就是为什么之前说大气动力学中更关心的是涡度的垂直分量，由此PV的主要成分为
$$
		\Pi \approx \frac{f}{\rho}\frac{ \partial \theta }{ \partial z } 
$$
位涡的特征数量级大概在0.5PVU左右，PVU即位涡对应的国际单位制单位的简写。
### 1.4.2 位涡应用举例
由于从对流层到平流层大气稳定性判据$\frac{ \partial \theta }{ \partial z }$会发生突变，因此PV值的突变处作为对流层顶的动力学定义是合理的。对于绝热无摩擦的流动，动力学对流层顶（等PV面）上的等位温线图对研究天气系统的情况很有用，因为位温的等值线会被风场平流输送从而直观地反映天气系统的移动和演变，而且相比于在等位温面上画等PV线，仅用一张等PV面上的位温图就能完整描述动力学对流层顶的天气状态。（在等位温面上画PV图则需要多张，因为一个位温只能在某些区域和对流层顶重合）
![[PV应用.png]]
上图a为500hPa等压面上的高度图，bcd均为在动力学对流层顶（等PV面）上画的图，分别展现了压强（高度）、位温和风速。在图a中展现出的是大气中压强的波动，有波峰——高压脊，波谷——低压槽，对应到图b可以看到槽对应于对流层顶下沉，脊对应于对流层顶抬升，这是因为槽即低压中心，空气辐合下沉，因此对流层顶被压到更低高度。在图c中，等位温线上的气团由于位涡守恒和位温守恒，不会离开等位温线，在其上形成物质涡旋（把等值线里面的空气关起来了），只有风场平流搬运才能改变其位置。在图c中，还可以看到长带状的位温强水平梯度区域， 这对应于冷暖气团交汇处，即对流层顶的峰区，这些峰区对应图d中的急流区。

### 1.4.3 完整Ertel位涡方程
更完整的Ertel位涡方程需要考虑到动量源$\mathbf{\mathcal{F}}$，如摩擦力，和熵源$\mathcal{H}$，如相变潜热
$$
\frac{D\Pi}{Dt}=\frac{\boldsymbol{\omega}_a}{\rho}\boldsymbol{\cdot}\nabla\mathcal{H}+\frac{\nabla\theta}{\rho}\boldsymbol{\cdot}\left(\nabla\times\frac{\mathbf{\mathcal{F}}}{\rho}\right)
$$
当涡度矢量指向熵源的局部最大值时，第一项为正，这多见于中低对流层处的成云降雨下方，此处为低压中心，北半球为逆时针，涡度向上指向云层中，而相变潜热的极值正在云层中，因此第一项为正。如果在云层上方，则第一项为负。

## 1.5 浅水波方程
对于各向同性的不可压缩流体，位涡守恒形式更加简单，为了得到这个，我们先考虑浅水波方程，这可以给我们未来分离大气动力学的重要部分提供帮助。

浅水波除了做了密度常数的近似以外，还做了水深度比起水平特征长度要小得多的近似，因此我们可以用流体静力平衡方程
$$
\frac{ \partial p }{ \partial z } =-\rho_{0}g
$$
设水表面高度为$h(x,y)$，表面压强$p(h)$认为是常数，从z处积分积到h处可得
$$
p(z)=\rho_{0}g(h-z)+p(h)
$$
将上式代入流体的水平运动方程中
$$
\frac{D_{h}\mathbf{V}}{Dt}=-g\nabla_{h}h-f\mathbf{k}\times \mathbf{V}\quad(*)
$$
下标h表示horizental水平。（$\frac{D_{h}}{Dt}$表示仅随水平运动变化）我们假设$\mathbf{V}$仅是x，y的函数（这可以认为源自浅水的近似），由于h也仅是x，y的函数，故$\mathbf{V}$始终保持为x，y的函数。

连续性方程在常数密度下写为
$$
\frac{ \partial u }{ \partial x } +\frac{ \partial v }{ \partial y } +\frac{ \partial w }{ \partial z }=0
$$
热力学第一定律在不可压缩流体的条件下，若流动绝热，则呈现为温度不变$dT=0$。水压强$p=p(T,\rho)$，由于$dT=0,d\rho=0$，故$dp=0$，即
$$
\frac{Dp}{Dt}=0
$$
考虑流体表面的一个粒子，其压强p应始终保持为$p(h)$，因此其始终保持在表面$z=h$，因此水面高度满足
$$
\frac{D_{h}h}{Dt}=\frac{Dz}{Dt}(h)=w(h)
$$
对连续性方程在z方向从0积到h，利用$w(0)=0$，
$$
w(h)=-h\nabla_{h}\cdot \mathbf{V}
$$
于是有
$$
\frac{D_{h}h}{Dt}=-h\nabla_{h}\cdot \mathbf{V}\quad(**)
$$
浅水波方程就由$(*)$和$(**)$式构成，描述了三个未知变量(u,v,h)的演化。对照前面涡度方程的推导过程，这里$\mathbf{V}$与z无关，且$\rho$为常数，因此第二项和第三项都没有了，简化的浅水波涡度方程写为
$$
\frac{D_h}{Dt}(\zeta+f)=-(\zeta+f)\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right)
$$
将涡度方程右边的$\nabla_{h}\cdot \mathbf{V}$用$(**)$式代入，可得守恒方程
$$
\frac{D_{h}}{Dt}\left( \frac{\zeta+f}{h} \right)=0
$$
因此浅水波位涡简化为了$\frac{\zeta+f}{h}$。
### 1.5.1 浅水波位涡的应用
利用简化后的浅水波位涡守恒，我们可以简单分析一个气体柱流过山脉的流动情况（这部分书里写的相当迷惑，很多地方感觉不太有道理，我只能按我自己的想法给出粗糙的解释，欢迎指正！）。在Chapter5中将会阐明，当气流向上靠近或向下远离山脉时，气体柱将会被垂直拉伸（h增大），当气体穿过山脉时，气体柱会被垂直压缩（h减小），后者是自然的，前者有待论证。
#### 1.5.1.1 西风流过山脉
![[西风流过山脉.png]]
如上图所示，西风流过山脉分成四个阶段。

1. 西风从西向东靠近山脉（上山），此时气体柱被垂直拉伸，h增大，由于位涡守恒，$\zeta$增大，形成气旋性流动，即向北流动，与此同时纬度增加，f增大，于是抑制了相对涡度$\zeta$的增大，形成负反馈调节。
2. 西风穿过山脉，此时由于地表的抬高，气体柱被压缩，h减小，由位涡守恒$\zeta$减小，形成反气旋性流动，即向南流动。
3. 西风向东远离山脉（下山），气体柱被垂直拉伸，h增大，$\zeta$增大，形成气旋性流动，即向北流动。
4. 气体回到原纬度时由于h的增大，$\zeta>0$，于是当气体回到原高度h时，其纬度比原纬度更高，于是气体有$\zeta<0$，进行反气旋性流动，南偏后f减小过度，又导致$\zeta>0$再次北偏，如此往复形成持续波动。
#### 1.5.1.2 东风流过山脉
前几个阶段与西风方向相反，但当气体流过山脉后的末态却与西风有显著差异，这是因为对于东风来说，$\zeta>0$时，风向南偏，f减小，于是为了保持位涡守恒$\zeta$要持续增大，形成正反馈调节，在这种调节下，风过山时会持续北偏，最后会达到稳态而不会波动，根据地转风方程，远离山后无气压梯度力，因此最终东风会平行流动。如下图所示
![[东风流过山脉.png]]
这种东风和西风的差异是由北半球$f>0$造成的，造成了很多有趣的地理现象。（以上解释甚至说服不了我自己，实在是看不懂书里在说什么）
### 1.5.2 正压位涡
在浅水波位涡的基础上，还可以对位涡进行进一步的简化，即假设流体为正压流体。这里先简单讨论一下正压流体，我认为所谓正压流体$p=p(\rho)$，指的是考虑一个流体块，其压强只与这个块的密度有关，而非考虑整个压强场（不然p总会与空间位置有关）。在这种观点下，对于浅水波假设下的一个流体块，其压强满足
$$
p(z)=p(h)+\rho_{0}g(h-z)
$$
若要求为正压流体，由于密度为常数，该流体块的压强也应保持为常数，故其h与z均为常数，故有$\frac{Dz}{Dt}=w=0$，代入连续性方程，可得正压浅水波是水平无散的
$$
\frac{ \partial u }{ \partial x } +\frac{ \partial v }{ \partial y } =0
$$
代入浅水波位涡方程可得
$$
\frac{D_{h}}{Dt}(\zeta+f)=0
$$
这就是进一步简化后的位涡（即绝对涡度），该房产称为正压涡度方程，在大尺度气象动力学理论研究中被广泛应用。

对于水平无散的运动，水平流速场可以用一个流函数$\psi(x,y)$来表征，其满足
$$
u=-\frac{ \partial \psi }{ \partial y } \quad and \quad v=\frac{ \partial \psi }{ \partial x }
$$
写成速度矢量形式为
$$
\mathbf{V}_{\psi}=\mathbf{k}\times \nabla \psi
$$
相对涡度为
$$
\zeta=\frac{ \partial v }{ \partial x } -\frac{ \partial u }{ \partial y } =\frac{ \partial^2 \psi }{ \partial x^2 } +\frac{ \partial^2 \psi }{ \partial y^2 } =\nabla_{h}^2\psi=\nabla^2\psi
$$
故速度场合涡度都可以用单一标量场$\psi(x,y)$来表示，正压涡度方程用流函数可以写为
$$
\frac{\partial}{\partial t}\nabla^2\psi=-\mathbf{V}_\psi\cdot\nabla\left(\nabla^2\psi+f\right)
$$
这个方程可以数值求解出流函数的演化，从而推算出速度场合涡度。由于对于中纬度天气尺度系统的确接近水平无散，这个预测方程对于这个情况下的短期预报竟然神奇的表现良好。
#### 1.5.2.1 正压位涡应用
与浅水位涡类似，对于正压位涡，即绝对涡度守恒，仍然存在东风和西风的不对称现象。

对于西风，向北（南）偏时对应$\zeta>0(<0)$，以及$f$增大（减小），无法保持绝对涡度守恒，因此西风不能偏转，始终保持自西向东的方向。而对于东风，向北（南）偏时对应$\zeta<0(>0)$，以及$f$增大（减小），可以保持绝对涡度守恒，因此东风可以向北也可以向南偏。
![[正压位涡举例.png]]
## 1.6 位温坐标下的Ertel位涡
这里我们详细讨论位温坐标下的Ertel位涡，包括在存在动量和熵源情况下的不守恒效应的证明。位温坐标的一大优势在于简化了热力学能量方程，将一个方程变成了一个变量，即绝热对应$\dot{\theta}=0$。
### 1.6.1 位温坐标下的运动方程
在位温坐标下垂直方向的“速度”为$\dot{\theta}\equiv \frac{D\theta}{Dt}$。对于水平运动方程，从p坐标到$\theta$坐标，除了气压梯度力需要重新推导以外，其他都只需要直接把前面推导涡度方程时用的运动方程里p改成$\theta$即可。

回顾p坐标下的水平运动方程
$$
\frac{\partial\mathbf{V}}{\partial t}=-\nabla\left(\frac{\mathbf{V}\cdot\mathbf{V}}{2}+\Phi\right)-(\zeta+f)\mathbf{k}\times\mathbf{V}-\omega\frac{\partial\mathbf{V}}{\partial p}
$$
现考虑气压梯度力$\frac{1}{\rho}\nabla p$，以x方向为例，即$\frac{1}{\rho}\left( \frac{ \partial p }{ \partial x } \right)_{z}$，其在p坐标下为$\left( \frac{ \partial \Phi }{ \partial x }  \right)_{p}$，要将其变到$\theta$坐标下
$$
\begin{aligned}
				\left( \frac{ \partial \Phi }{ \partial x }  \right)_{\theta}&=\left( \frac{ \partial \Phi }{ \partial x }  \right)_{p}+\frac{ \partial \Phi }{ \partial p } \left( \frac{ \partial p }{ \partial x }  \right)_{\theta}\\&=\left( \frac{ \partial \Phi }{ \partial x }  \right)_{p}-\frac{1}{\rho}\left( \frac{ \partial p }{ \partial x }  \right)_{\theta}
\end{aligned}
$$
在固定$\theta$的条件下，有
$$
c_{p}dT=\frac{1}{\rho}dp
$$
因此有
$$
\frac{1}{\rho}\left( \frac{ \partial p }{ \partial x }  \right)_{\theta}=\left( \frac{ \partial c_{p}T }{ \partial x }  \right)_{\theta}
$$
于是可得$\theta$坐标下的气压梯度力
$$
\left( \frac{ \partial (\Phi+c_{p}T) }{ \partial x }  \right)_{\theta}
$$
于是可以定义Montgomery流函数$\Psi \equiv c_{p}T+\Phi$，水平运动方程可以写为
$$
\frac{\partial\mathbf{V}}{\partial t}+\nabla_\theta\left(\frac{\mathbf{V}\cdot\mathbf{V}}{2}+\Psi\right)+(\zeta_\theta+f)\mathbf{k}\times\mathbf{V}=-\dot{\theta}\frac{\partial\mathbf{V}}{\partial\theta}+\mathbf{F}_r\quad(1)
$$
其中垂直涡度为$\zeta_{\theta}=\mathbf{k}\cdot (\nabla_{\theta}\times \mathbf{V})$。

对于连续性方程，考虑$\theta$坐标下的一个小流体块，根据Jacobi行列式可得质量元的变换
$$
\delta M=\rho\delta x\delta y\delta z=-\frac{1}{g}\delta x\delta y\delta p=-\frac{1}{g}\delta x\delta y\left( \frac{ \partial p }{ \partial \theta }  \right)\delta \theta=\sigma \delta x\delta y\delta \theta
$$
在$\theta$坐标下的密度定义为
$$
\sigma \equiv-\frac{1}{g}\frac{ \partial p }{ \partial \theta }\quad(2)
$$
于是质量守恒可导出$\theta$坐标下的连续性方程
$$
\frac{ \partial \sigma }{ \partial t } +\nabla_{\theta}\cdot(\sigma \mathbf{V})=-\frac{ \partial  }{ \partial \theta } (\sigma \dot{\theta})\quad(3)
$$
对于垂直运动方程，即流体静力平衡方程，在$\theta$坐标下可以转换为流函数$\Psi$和$\theta$的偏导关系，对于流函数$\Psi=c_{p}T+gz$，其T也可以看成是z的函数，而z是p的函数，因此可以把$\Psi$看成是p的函数，故有
$$
\frac{ \partial \Psi }{ \partial \theta } =\frac{ \partial \Psi }{ \partial p } \frac{ \partial p }{ \partial \theta }
$$
根据$\Psi=c_{p}T+gz$，$\theta=T\left( \frac{p_{s}}{p} \right)^{R/c_{p}}$，有
$$
\frac{\partial\Psi}{\partial\theta}=\Pi(p)\equiv c_p\left(\frac{p}{p_s}\right)^{R/c_p}=c_p\frac{T}{\theta}\quad(4)
$$
$\Pi$称作Exner函数。上面的(1)(2)(3)(4)个式子，刚好有$\mathbf{V},\sigma,\Psi,p$四个变量（如果$\dot{\theta},\mathbf{F}_{r}$已知），因此它们可以作为一组预测方程。
### 1.6.2 位涡方程
考虑笛卡尔坐标系里的Ertel位涡近似表达式$\Pi=\frac{\zeta+f}{\rho}\frac{ \partial \theta }{ \partial z }$，将其用位温坐标定义的量表达即
$$
\Pi=(\zeta_{\theta}+f)\left( -g\frac{ \partial \theta }{ \partial p }  \right)=\frac{\zeta_{\theta}+f}{\sigma}
$$

下面要想办法凑出上面的形式来，用$\mathbf{k}\cdot \nabla_{\theta}\times$作用在水平运动方程的等式两边，可得
$$
\frac{\tilde{D}}{Dt}\left(\zeta_\theta+f\right)+\left(\zeta_\theta+f\right)\nabla_\theta\cdot\mathbf{V}=\mathbf{k}\cdot\nabla_\theta\times\left(\mathbf{F}_r-\dot{\theta}\frac{\partial\mathbf{V}}{\partial\theta}\right)
$$
其中$\frac{\tilde{D}}{Dt}=\frac{\partial}{\partial t}+\mathbf{V}\cdot\nabla_\theta$表示在$\theta$坐标下随水平方向运动的物质导数，这就是位温涡度方程

要将连续性方程写成含$\sigma^{-1}$的形式，且形式与上面的运动方程类似，利用
$$
\nabla_{\theta}\cdot(\sigma \mathbf{V})=\sigma \nabla_{\theta}\cdot \mathbf{V+(\mathbf{V}\cdot \nabla_{\theta})\sigma}
$$
代入原连续性方程可得
$$
	\frac{\tilde{D}\sigma}{Dt}+\sigma \nabla_{\theta}\cdot \mathbf{V}=-\frac{ \partial  }{ \partial \theta } (\sigma \dot{\theta})
$$
两边同乘$-\sigma^{-2}$，可得
$$
\frac{\tilde{D}}{Dt}\left(\sigma^{-1}\right)-\left(\sigma^{-1}\right)\nabla_\theta\cdot\mathbf{V}=\sigma^{-2}\frac{\partial}{\partial\theta}\left(\sigma\dot{\theta}\right)
$$
给运动方程两边同乘$\sigma^{-1}$，连续性方程两边同乘$\zeta_{\theta}+f$，再相加，可得位涡方程
$$
\frac{\widetilde{D}\Pi}{Dt}=\frac{\partial\Pi}{\partial t}+\mathbf{V}\cdot\nabla_\theta\Pi=\frac{\Pi}{\sigma}\frac{\partial}{\partial\theta}\left(\sigma\dot{\theta}\right)+\sigma^{-1}\mathbf{k}\cdot\nabla_\theta\times\left(\mathbf{F}_r-\dot{\theta}\frac{\partial\mathbf{V}}{\partial\theta}\right)
$$
可见类似的，对于绝热无摩擦的流动，位温守恒即$\dot{\theta}=0$，且$\mathbf{F}_{r}=0$，故$\theta$坐标下的Ertel位涡守恒。
### 1.6.3 位温涡度的积分约束
上面得到的位温涡度方程可写为
$$
\frac{\partial\zeta_\theta}{\partial t}=-\nabla_\theta\cdot[(\zeta_\theta+f)\mathbf{V}]+\mathbf{k}\cdot\nabla_\theta\times\left(\mathbf{F}_r-\dot{\theta}\frac{\partial\mathbf{V}}{\partial\theta}\right)
$$
利用矢量分析公式
$$
\mathbf{k}\cdot(\nabla_{\theta}\times \mathbf{A})=\nabla_{\theta}\cdot(\mathbf{A}\times \mathbf{k})
$$
于是位温涡度方程可以写为
$$
\frac{\partial\zeta_\theta}{\partial t}=-\nabla_\theta\cdot\left[\left(\zeta_\theta+f\right)\mathbf{V}-\left(\mathbf{F}_r-\dot{\theta}\frac{\partial\mathbf{V}}{\partial\theta}\right)\times\mathbf{k}\right]=-\nabla_{\theta}\cdot \mathbf{J}
$$
可见位温涡度改变只能由于$\mathbf{J}$矢量的水平散度，即垂直穿过等位温线的传输不能改变位温涡度。

对上面式子两边在地球表面闭曲面积分，根据斯托克斯定理容易证明
$$
\displaystyle {\iint\kern{-20mu}{\unicode{x2B2D}}\,}_{\kern{-12mu}A}\nabla_{\theta}\cdot \mathbf{J}dA=0

$$
(具体证明需先证明$\int \int\nabla_{\theta}\cdot \mathbf{J}dA=\oint \mathbf{J}\cdot \mathbf{n}dA$，该定理可通过构造矢量$\mathbf{F}$满足$\nabla \times \mathbf{F}=\nabla \cdot \mathbf{J}$来证明)

于是有
$$
	\frac{ \partial  }{ \partial t } \displaystyle {\iint\kern{-20mu}{\unicode{x2B2D}}\,}_{\kern{-12mu}A}\zeta_{\theta}dA=0
$$
即地球表面平均位温涡度为常数。更进一步，利用$\zeta_{\theta}=\mathbf{k}\cdot(\nabla_{\theta}\times \mathbf{V})$，可得
$$
\begin{aligned}
\displaystyle {\iint\kern{-20mu}{\unicode{x2B2D}}\,}_{\kern{-12mu}A}\zeta_{\theta}dA&=\displaystyle {\iint\kern{-20mu}{\unicode{x2B2D}}\,}_{\kern{-12mu}A}(\nabla_{\theta}\times \mathbf{V})\mathbf{k}dA\\&=\oint \mathbf{V}\cdot d\mathbf{l}=0
\end{aligned}
$$
因此地球表面平均位温涡度是0，综上位温涡度此消彼长，仅由平行于等位温面的水平流$\mathbf{J}$导致稀释或浓缩（即减小或增大）。