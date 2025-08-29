---
zhihu-title: Holton读书笔记-Chapter5 大气振荡（线性摄动理论）
zhihu-topics: 大气科学
zhihu-link: https://zhuanlan.zhihu.com/p/1934580252622979224
zhihu-created-at: 2025-08-05 20:23
zhihu-updated-at: 2025-08-09 16:27
toc: "true"
tags:
  - zhihu
---
```toc
```

[[Holton读书笔记-Chapter1 基本介绍]](https://zhuanlan.zhihu.com/p/1928544934622924935 "card")
[[Holton读书笔记-Chapter2 基础守恒律]](https://zhuanlan.zhihu.com/p/1929242659588912307 "card")
[[Holton读书笔记-Chapter3 基本方程的初步应用]](https://zhuanlan.zhihu.com/p/1930180220087939862 "card")
[[Holton读书笔记-Chapter4 环量、涡度和位涡]](https://zhuanlan.zhihu.com/p/1934567913379063792 "card")

# 1 Chapter5
大气模型的固有复杂性导致我们无法用一个简单的模型就精确解释和预测期物理过程，但为了增加对大气运动的基本物理认识，研究忽略一定过程后的简单模型是很有意义的，因此这一章的主题就是对理想化情况下的大气进行分析，首先我们会介绍摄动法，一种定性分析大气中的波动的方法，然后会用该方法研究大气中几种形式的波动，在Chapter6和7，我们分别会用该方法进一步分析准地转方程和天气波动的发展。
## 1.1 摄动法
在摄动法中，所有场变量被分成基态部分和扰动部分，基态部分被设置为关于时间和经度的平均值，即其不依赖于时间和经度，而扰动量就表示局部场偏离基态的部分。

比如对于纬向速度 $u$，取$\bar{u}$为其时间和经度的平均速度，$u'$为偏离，则实际速度场为
$$
u(x,t)=\bar{u}+u'(x,t)
$$
于是其随动造成的加速度 $u\frac{ \partial u }{ \partial x }$可写为
$$
u\frac{ \partial u }{ \partial x } =(\bar{u}+u')\frac{ \partial  }{ \partial x } (\bar{u}+u')=\bar{u}\frac{ \partial u' }{ \partial x } +u'\frac{ \partial u' }{ \partial x }
$$
摄动理论的基本假设为
1. 基态变量也满足控制方程(即均取为基态量后原来的控制方程也都满足)
2. 扰动量相比于基态量极小（两个扰动量的乘积可忽略）

于是上式可以近似为
$$
u\frac{ \partial u }{ \partial x } =\bar{u}\frac{ \partial u' }{ \partial x }
$$
可见当基态给定时，非线性项被近似为了扰动量的线性项，微分方程也变成了线性微分方程，数学上容易求解。对于常系数线性微分方程，其解为正弦型或指数型，然后扰动量的解决定了波动的传播速度、垂直结构、增长或衰减等。
## 1.2 波的性质
波动运动是在时间和空间上传播的场变量的振荡。对于传播的波，其相位与空间和时间均有关，写为 $\Phi(x,t)=kx-vt-\alpha$，其相位的传播速度称为相速度（即等相位面移动的速度），按以下方式推导出了相速度c
$$
kx-vt-\alpha=k(x+dx)-v(t+dt)-\alpha\Longrightarrow c \equiv \frac{dx}{dt}=\frac{v}{k}
$$
### 1.2.1 傅里叶级数
将扰动项简单的当成正弦波看上去有些简化过度，在大气中的扰动不可能是纯粹的正弦波，因此我们需要采用傅里叶级数展开将经度的函数展开成经度平均值加上关于经度的正弦波的线性组合。

我们设一个关于经度的函数减去其经度平均值为 $f(x)$，对x采用周期性边界条件，即取$f(x)$周期为该纬度下一圈的周长L，则其傅里叶级数为
$$
f(x)=\sum_{s=1}^\infty(A_{s}\sin k_{s}x+B_{s}\cos k_{s}x)
$$
其中波数$k_{s}=\frac{2\pi s}{L}$，s为该波绕纬度一圈的波的个数，系数$A_{s},B_{s}$可以通过三角函数的正交性求出
$$
A_{s}=\frac{2}{L}\int_{0}^{L}f(x)\sin \frac{2\pi sx}{L}dx
$$
$$
B_{s}=\frac{2}{L}\int_{0}^{L}f(x)\cos \frac{2\pi sx}{L}dx
$$
我们称$f_{s}(x)=A_{s}\sin k_{s}x+B_{s}\cos k_{s}x$为s阶傅里叶分量。以等压面上的地势扰动为例，将其沿傅里叶级数展开，其中振幅最大的一个傅里叶分量，其对应的波的个数s就接近于实际观测到的高压槽或低压脊的个数。因此如果只要求做定性的分析，那可以取振幅最大的傅里叶分量做分析。

傅里叶分量也可以更紧凑地写成复数的形式
$$
f_{s}(x)=\mathrm{Re}[C_{s}\exp(ik_{s}x)]
$$
其中$C_{s}$是复数，将其转换成三角函数形式系数对应为
$$
B_{s}=\mathrm{Re}[C_{s}]\quad and \quad A_{s}=-\mathrm{Im}[C_{s}]
$$
### 1.2.2 色散和群速度
线性振子的基本性质就是其振荡频率$\nu$仅依赖于振子的物理特征（如重力加速度和绳长）而不依赖于运动本身（如振幅、波长），然而对于传播的波，其频率一般会依赖于其波数$k$，如果不满足$\nu\propto k$，则波的相速度$c=\frac{\nu}{k}$就会随k而改变，于是原本在同一个地方的许多不同波数的正弦波由于传播速度的不同就会分散开来，这就是**色散**，而$\nu=\nu(k)$就称为色散关系。

对于非色散波（如声波），其在一处由一系列傅里叶分量正弦波组成，则随着时间演化这个波群会保持形状不变地在空间按相速度传播。而对于色散波，其波群不会形状保持不变，会逐渐散开来，由于波之间的相消和相长，能量会在一些区域发生聚集，一般来说，波群的传播速度（能量）和傅里叶分量的平均相速度不同。

关于群速度的表达式，书里给出的是一种简化版的，即仅考虑两个频率差距很小的波的叠加，这里我直接对任意形式的波包的群速度做推导。注意波包需要各成分波的频率与平均频率之差很小，即准单色波。（一步到位，直接考虑三维形式的波）

考虑波包由一组波矢频率相近的平面波叠加而成
$$
\psi(\mathbf{r},t)=\int A(\mathbf{k})e^{i(\mathbf{k}\cdot\mathbf{r}-\nu(\mathbf{k})t)}d^3k
$$
设各成分波的平均频率为$\nu_{0}$，平均波数为$k_{0}$，$\nu(\mathbf{k})$在$\mathbf{k}_{0}$处可泰勒展开为
$$
\nu(\mathbf{k})\approx \nu_{0}+(\mathbf{k}-\mathbf{k}_{0})\cdot \nabla_{\mathbf{k}}\nu|_{\mathbf{k}_{0}}
$$
则波动可写为
$$
\begin{aligned}
\psi(\mathbf{r},t)&= e^{i(\mathbf{k}_0\cdot\mathbf{r}-\nu_0t)}\int A(\mathbf{k})e^{i(\mathbf{k}-\mathbf{k}_0)\cdot(\mathbf{r}-\mathbf{c}_gt)}d^3k\\&=F(\mathbf{r}-\mathbf{c}_{g}t)e^{i(\mathbf{k}_0\cdot\mathbf{r}-\nu_0t)}
\end{aligned}
$$
其中$\mathbf{c}_{g}=\nabla_{\mathbf{k}}\nu|_{\mathbf{k}_{0}}$。由于平面波波矢频率均相近，因此波动的系数$F$对应缓变项，调制波动的振幅，而后面的e指数则对应快变项（载波），直接代表波动，于是$F$对应波动的包络线，由于能量与振幅相对应，包络线对应振幅，因此能量传播的速度即群速度就对应包络线移动的速度，即F移动的速度
$$
\mathbf{c}_{g}=\nabla_{\mathbf{k}}\nu|\mathbf{k}_{0}
$$
一维下群速度则为$c_{g}=\frac{ \partial \nu }{ \partial k }|_{k_{0}}$
### 1.2.3 二维和三维下波的性质
前面仅讨论了一维传播的波，将其扩展到二维和三维是非常必要的，为了简单，我们仅讨论二维波，用类似的方法可以再扩展到三维。

二维波动表达式可以写为
$$
f(x,y,t)=\mathrm{Re}\left[Ae^{i(kx+ly-\nu t)}\right]=\mathrm{Re}\left[Ae^{i(\mathbf{k}\cdot \mathbf{r}-\nu t)}\right]=\mathrm{Re}\left[Ae^{i\phi}\right]
$$
k,l分别是波沿x和y的波数，于是波矢可以写为
$$
	\mathbf{K}=k\mathbf{i}+l\mathbf{j}=\nabla \phi
$$
总波数$\mathcal{K}=|\mathbf{K}|$，波长$\lambda=\frac{2\pi}{\mathcal{K}}$，（圆）频率为$\nu=-\frac{ \partial \phi }{ \partial t }$。将x轴取成沿$\mathbf{K}$方向可以将其变为一维，可得相速度为
$$
c=\frac{\nu}{\mathcal{K}}=-\frac{1}{|\nabla \phi|}\frac{ \partial \phi }{ \partial t }
$$
其方向沿波矢$\mathbf{K}$方向，注意相速度不能算是矢量，因为x方向相速度$c_{x}=\frac{\nu}{k}$和y方向相速度$c_{y}=\frac{\nu}{l}$的矢量叠加不等于相速度$c$。二维情况下的群速度在前面已经给过证明，为$\mathbf{c}_{g}=\nabla_{\mathbf{K}}\nu|_{\mathbf{K}_{0}}$

如果我们考虑一个与平面波略不同的波动，其波矢和频率均会随时间和空间缓慢发生变化。此时前面关于波矢和频率的公式都不会变，可得
$$
\frac{\partial\mathbf{K}}{\partial t}+\nabla\nu=0
$$
考虑色散关系$\nu=\nu(\mathbf{K})$，则$\nabla \nu$可以用链式法则写成关于波矢的偏导数
$$
\nabla \nu=(\nabla_{\mathbf{K}}\nu \cdot \nabla)\mathbf{K}=(\mathbf{c}_{g}\cdot \nabla)\mathbf{K}
$$
故有
$$
\frac{ \partial \mathbf{K} }{ \partial t } +(\mathbf{c}_{g}\cdot \nabla)\mathbf{K}=0
$$
因此在以群速度运行的参考系里看，波矢始终是守恒的。
### 1.2.4 波动解方法
这里给出用近似波动解求解问题的一般方法。
1. 选择基态量：基态是大气状态的简化表达，一般我们希望消除所有的复杂项，仅保留最至关重要、最感兴趣的项（比如取流体静力平衡的基态或者流体完全静止的基态等）（个人认为选择基态量的关键是确定该量的扰动是否远小于该量本身，比如不能取压强的基态为沿z方向的平均值，因为沿z方向p的变化很大）
2. 线性化控制方程：用摄动法可以在给定基态的情况下把微分方程转化为线性的。
3. 假设波动解：如果微分方程的系数不依赖于某个独立变量x，且沿x方向有周期性，那么可以设为关于x的波动形式的解，如果与某个变量y相关的微分方程不是常系数的（基态与y有关），或者沿y方向没有周期性，那么只能将其设成一般形式函数$F(y)$然后将偏微分方程转化为关于$F(y)$的常微分方程求解。
4. 解出色散和相位关系：将波动形式的解或者一般函数形式带回原方程中可以得到一系列常微分方程以及频率与波矢的关系，后者就是色散关系。对于波动解，振幅大小方程不会做约束（只能解出不同变量之间振幅的比例关系），但变量之间的相位关系可以解出，如同相、反相、差90度等。
## 1.3 简单的波动类型
### 1.3.1 声波
声波是纵波，其振动方向平行于传播方向，声音由交替的绝热压缩和拉伸传播，即对应于压强的波动。我们考虑沿着x方向直线传播的一维声波，因此有$v=w=0$，且$u=u(x,t)$，在这个限制条件下，流体力学三大方程：运动方程、连续性方程、热力学能量方程（绝热下）可以写为
$$
\frac{Du}{Dt}+\frac{1}{\rho}\frac{ \partial p }{ \partial x } =0
$$
$$
\frac{D\rho}{Dt}+\rho \frac{ \partial u }{ \partial x } =0
$$
$$
c_{p}dT-\frac{1}{\rho}dp=0
$$
绝热方程利用理想气体状态方程可以写为关于$p$和$\rho$的形式
$$
\frac{1}{\gamma p}\frac{Dp}{Dt}=\frac{1}{\rho} \frac{D\rho}{Dt}
$$
可以消去方程组中的$\frac{D\rho}{Dt}$，(在摄动法中，只有偏导数项会作为未知量，前面的系数由于近似只会保留基态量)
$$
\frac{1}{\gamma p} \frac{Dp}{Dt}+\frac{ \partial u }{ \partial x } =0
$$
用摄动法分解
$$
u(x,t)=\bar{u}+u'(x,t)
$$
$$
p(x,t)=\bar{p}+p'(x,t)
$$
$$
\rho(x,t)=\bar{\rho}+\rho'(x,t)
$$
代入上面的两个方程可得
$$
\left( \frac{ \partial  }{ \partial t } +\bar{u}\frac{ \partial  }{ \partial x }  \right)u'+\frac{1}{\bar{\rho}}\frac{ \partial p' }{ \partial x } =0
$$
$$
\frac{1}{\gamma \bar{p}}\left( \frac{ \partial  }{ \partial t } +\bar{u}\frac{ \partial  }{ \partial x }  \right)p'+\frac{ \partial u' }{ \partial x } =0
$$
下式用算子$\frac{ \partial  }{ \partial t }+\bar{u}\frac{ \partial  }{ \partial x }$作用，上式用算子$\frac{ \partial  }{ \partial x }$作用，两式相减可得
$$
\left( \frac{ \partial  }{ \partial t } +\bar{u}\frac{ \partial  }{ \partial x }  \right)^2p'-\frac{\gamma \bar{p}}{\bar{\rho}}\frac{ \partial ^2p' }{ \partial x^2 } =0
$$
这是比较标准的波动方程形式，可以设波动解
$$
p'=A\exp\left[ik(x-ct)\right]
$$
带回可得
$$
(-ikc+ik \bar{u})^2-\frac{\gamma \bar{p}}{\bar{\rho}}(ik)^2
$$
$$
c=\bar{u}\pm \sqrt{ \gamma \bar{p} /\bar{\rho}}=\bar{u}\pm \sqrt{ \gamma R \bar{T} }
$$
可见波速相对于纬向的平均速度流的速度为
$$
c_{s}=\sqrt{ \gamma R \bar{T} }
$$
这称为绝热声速。

纬向的平均速度流仅起到类似多普勒频移效应的作用，频率为
$$
\nu=k|c|=k|\bar{u}\pm c_{s}|
$$
以$\bar{u}>0$为例，声源静止，则在声源下游即东面听到的频率为$k(\bar{u}+c_{s})$，而在上游西面听到的频率为$k(c_{s}-\bar{u})$，因此下游听到的频率比上游更高。
### 1.3.2 浅水波
回顾一下在Chapter4中给出的浅水波方程，
$$
\frac{D\mathbf{V}}{Dt}=-g\nabla h-f\mathbf{k}\times \mathbf{V}
$$
$$
\frac{Dh}{Dt}=-h\nabla \cdot \mathbf{V}
$$
对于浅水波，我们取基态为静止状态，则速度均为0，深度为$\bar{h}$，设扰动项为$u',v',h'$，代入浅水波方程可得

$$
\begin{aligned}
&\frac{ \partial u' }{ \partial t } =-g\frac{ \partial h' }{ \partial x } +fv'\quad(1)\\&\frac{\partial v^{\prime}}{\partial t}=-g\frac{\partial h^{\prime}}{\partial y}-fu^{\prime}\quad(2)\\&\frac{\partial h^{\prime}}{\partial t}=-\bar{h}\left(\frac{\partial u^{\prime}}{\partial x}+\frac{\partial v^{\prime}}{\partial y}\right)\quad(3)\end{aligned}
$$
想办法消去$u,v$将方程写为h的偏微分方程。对(3)式关于t求偏导，再代入(1)(2)式可得
$$
\frac{ \partial ^2h' }{ \partial t^2 } =g \bar{h}\nabla_{h}^2h'-f \bar{h}\left( \frac{ \partial v' }{ \partial x } -\frac{ \partial u' }{ \partial y }  \right)
$$
再对上式关于t求偏导，代入(1)(2)式可得
$$
\frac{ \partial ^3h' }{ \partial t^3 } =g \bar{h}\nabla_{h}^2\left( \frac{ \partial h' }{ \partial t }  \right)+f^2 \bar{h}\left( \frac{ \partial u' }{ \partial x } +\frac{ \partial v' }{ \partial y }  \right)
$$
再代入(3)式可得只含h的偏微分方程
$$
\frac{\partial^3h^{\prime}}{\partial t^3}+\left(f^2-g\overline{h}\nabla_h^2\right)\frac{\partial h^{\prime}}{\partial t}=0
$$
假设侧面的边界是周期性的，可以设波动解
$$
h'=A\exp[i(kx+ly-\nu t)]
$$
代入可得
$$
\nu^3-\nu\left[f^2+g\bar{h}\left(k^2+l^2\right)\right]=0
$$
一个解为$\nu=0$，另一个为
$$
\nu^2=f^2+g \bar{h}(k^2+l^2)
$$
求得特征频率的另一种方式是假设了波动解之后直接代入原方程组中，然后利用方程组有非零解的性质，要求变量(u',v',h')的系数矩阵必须行列式等于0，从而解出特征频率$\nu$

对于第一个解$\nu=0$，波是静止的，对应方程(1)(2)(3)左边均为0，于是有
$$
\begin{aligned}
		&v'=\frac{g}{f}\frac{ \partial h }{ \partial x } \\&u'=-\frac{g}{f}\frac{ \partial h }{ \partial y } 
\end{aligned}
$$
波动对应地转平衡情况，是在限定$f$为常数下的**Rossby波**的一个例子。传播的Rossby波（$f$不是常数）的情况会在后面讨论。

对于$\nu$的非零解，这种波称作**惯性重力波**，因为粒子振荡依赖于重力和惯性力，在后面会重点讨论惯性项，在这里我们先将注意力放在重力项上。在限定$f=0$的情况下，我们得到了简单重力波
$$
\nu^2=g \bar{h}(k^2+l^2)
$$
于是可得相速度为
$$
c=\frac{\nu}{\mathcal{K}}=\frac{\nu}{\sqrt{ k^2+l^2 }}=\sqrt{ g \bar{h} }
$$
相速度为常数，因此浅水波是无色散波。（注意：浅水波方程仅适用于波长尺度远大于水深度的情况，如地震或火山喷发造成的水波）

#### 1.3.2.1 进一步讨论（换坐标）
为了进一步揭示浅水波运动的物理性质，我们将x轴设定为沿着波传播的方向，于是波矢的y分量$l=0$，（注意此时u',v'不再按照经向、纬向划分了，而是垂直于等相位面和平行于等相位面）将波动解假设带回原始方程组中可得
$$
\begin{aligned}
&-i\nu u'=-ikgh'+fv'\\&-i\nu v'=-fu'
\end{aligned}
$$
解得
$$
\begin{aligned}
&u'=\frac{\nu kg}{\nu^2-f^2}h'\\&v'=-\frac{ikgf}{\nu^2-f^2}h'
\end{aligned}
$$
根据连续性方程和上面推出的速度公式
$$
\frac{ \partial w' }{ \partial z } =-\frac{ \partial u' }{ \partial x }=\frac{1}{\bar{h}}\frac{ \partial h' }{ \partial t } =-\frac{i\nu}{\bar{h}}h'
$$
由于$h'$只与x，y有关，因此可以直接对$w'$从0到z积分，利用底部垂直速度为0，可得
$$
w'(z)=-\frac{i\nu z}{\bar{h}}h'
$$
可见$w'$和$h'$差了90度相位。

对于$\nu=0$的情况，代入上面推导的公式可得，$u'=0$，$v'=\frac{ikg}{f}h'$，$w'=0$。对于$\nu \ne 0$的重力惯性波，若限定$f=0$，则有$v'=0$，$u'=\frac{\nu }{\bar{h}k}h'=\frac{c}{\bar{h}}h'$，表明当$c>0$，即波沿x轴正向传播时，$u'$与$h'$同相，当$c<0$时，则反相。
![[浅水波示意图.png]]
上图为$c>0$，浅水波向右传播的示意图，$u'$的辐合和辐散决定了水深$h'$的升高和降低，水深的升高和水深的高位差90度相位，因此图中显示上升水流对应的高度为平均高度$\bar{h}$。又由于$u'$和$h'$同相，因此水深高位对应于水流向右，低位对应于水流向左。

#### 1.3.2.2 涡度和位涡
在把x坐标取为沿波传播方向后，有$\frac{ \partial u' }{ \partial y }=0$，相对涡度可以写为
$$
\zeta'=\frac{ \partial v' }{ \partial x } =\frac{k^2gf}{\nu^2-f^2}h'
$$
浅水波位涡方程为
$$
\frac{D_{h}}{Dt}\left( \frac{\zeta+f}{h} \right)=0
$$
在摄动法下可以近似为
$$
\frac{ \partial  }{ \partial t } \left( \frac{\zeta'+f}{\bar{h}+h'} \right)=\frac{ \partial  }{ \partial t } \left( \frac{\zeta'-\frac{f}{\bar{h}}h'}{\bar{h}} \right)=0
$$
表明摄动法线性化以后的位涡为
$$
Q'=\zeta'-\frac{f}{\bar{h}}h'
$$
其在局域守恒，即固定一点处的$Q'$不变，这对于线性动力学是非常强的约束条件。

对于浅水Rossby波（$\nu=0$），其涡度为
$$
\zeta_{RW}'=-\frac{k^2g}{f}h'
$$
与水深度场反相，根据波矢k的定义，当波长很短时，涡度很大。Rossby波的位涡除了在波长非常长的长波中会与其涡度反号，大部分时候均与其涡度同号。

对于浅水的惯性重力波，其涡度为
$$
\zeta_{GW}'=\frac{f}{\bar{h}}h'
$$
与水深度场同相，其位涡为0 
#### 1.3.2.3 总结
浅水波下的Rossby波和惯性重力波的主要不同如下
1. Rossby波在$f$为常数时静止，而惯性重力波波速非常大
2. 惯性重力波的相对涡度较小，但散度较大，而Rossby波为地转风，散度为0，相对涡度较大
3. 惯性重力波位涡为0，Rossby波有非0位涡

这些性质提供了一种有用的框架，可以从完全非线性的方程中排除惯性重力波，得到准地转方程，这在后面研究热带外的天气系统的动力学时很有用。
## 1.4 重力内（浮力）波
现在考虑在无背景旋转（即地球自转）下分层大气的重力波的传播。注意大气重力波只能发生于稳定分层的大气中，只有稳定气块才能发生如Chapter2所讨论的浮力振荡。

对于一个流体，比如海洋，其上部和下部都有边界，则传播的重力波更多是在水平平面内发生，因为垂直传播的波会在边界反射从而形成驻波。而对如大气的流体，没有上部边界，则重力波会在水平和垂直方向传播。在垂直方向传播的的波其相位是高度的函数，这种波称为内波（internal wave），尽管内波在天气尺度的天气预报中一般不重要，但它对于中尺度运动很重要，比如它导致了山的背风波，它也被认为是往中层大气传递能量动量的一种重要机制。

### 1.4.1 纯重力内波
为了简化，我们考虑无背景旋转，即$f=0$的情况，且把讨论限制在传播在x,z平面的二维重力内波，这种波的频率由Chapter2讨论的浮力模型来给出。

重力内波是横波，粒子振荡方向垂直于波传播方向，设粒子位移为$\delta s$，位移方向与垂直方向的夹角为$\alpha$，单位质量气块垂直方向浮力为$-N^2\delta z=-N^2\delta s\cos \alpha$，其中$N^2=g  \frac{d\ln\theta_{0}}{dz}$，因此平行于粒子运动路径的力为$-N^2\cos^2\alpha \delta z$，粒子振荡的运动方程为
$$
		\frac{d^2(\delta s)}{dt^2}=-(N\cos \alpha)^2\delta s
$$
解为
$$
\delta s=\exp[\pm i(N\cos \alpha)t]
$$
因此气块的简谐振动频率为$\nu=N\cos \alpha$，仅依赖于传播方向与垂直方向的夹角和大气的稳定性因子N。

上面的讨论用启发性的粒子模型给出了重力内波的频率，但是却没有给出波动的方程，现在我们考虑一个二维(x,z)重力内波的线性化方程（摄动法）。在绝热情况下，决定大气的运动方程、连续性方程、热力学能量方程如下
$$
\begin{gathered}\frac{\partial u}{\partial t}+u\frac{\partial u}{\partial x}+w\frac{\partial u}{\partial z}+\frac{1}{\rho}\frac{\partial p}{\partial x}=0\\\frac{\partial w}{\partial t}+u\frac{\partial w}{\partial x}+w\frac{\partial w}{\partial z}+\frac{1}{\rho}\frac{\partial p}{\partial z}+g=0\\ \frac{\partial u}{\partial x}+\frac{\partial w}{\partial z}=0\\\frac{\partial\theta}{\partial t}+u\frac{\partial\theta}{\partial x}+w\frac{\partial\theta}{\partial z}=0\end{gathered}
$$
位温$\theta$满足
$$
\ln \theta=\frac{1}{\gamma}\ln p-\ln \rho +C
$$
采用摄动法，取密度基态为常数密度$\rho_{0}$，压强基态为满足流体静力平衡方程的x方向平均值$\bar{p}(z)$，根据压强和密度的基态可得位温基态，水平速度的基态取为二维平均值$\bar{u}$，垂直速度的基态取为0，可得
$$
\begin{aligned}
&\rho=\rho_{0}+\rho'\\&p=\bar{p}(z)+p'\\&\theta=\bar{\theta}(z)+\theta'\\&u=\bar{u}+u'\\&w=w'
\end{aligned}
$$
其中有
$$
\frac{d\bar{p}}{dz}=-\rho_{0}g
$$
$$
\ln \bar{\theta}=\frac{1}{\gamma}\ln \bar{p}-\ln \rho_{0}+C
$$
接下来想办法消去$\rho'$（最好消去的量），第一个方程中利用近似$\rho'$可以忽略，考虑第二个方程中的最后两项
$$
\frac{1}{\rho}\frac{ \partial p }{ \partial z } +g=\frac{1}{\rho_{0}}\frac{ \partial p' }{ \partial z } +\frac{\rho'}{\rho_{0}}g
$$
想办法把$\rho'$用$\theta$和$p$表示出来，利用$\ln \theta=\frac{1}{\gamma}\ln p-\ln \rho +C$，代入扰动项可得
$$
\begin{aligned}\ln\left[\overline{\theta}\left(1+\frac{\theta^{\prime}}{\overline{\theta}}\right)\right]&=\frac{1}{\gamma}\ln\left[\overline{p}\left(1+\frac{p^{\prime}}{\overline{p}}\right)\right]\\&-\ln\left[\rho_0\left(1+\frac{\rho^{\prime}}{\rho_0}\right)\right]+\mathrm{const.}\end{aligned}
$$
由于基态量之间关系不变，因此可得近似
$$
\frac{\theta'}{\bar{\theta}}=\frac{p'}{\gamma \bar{p}}-\frac{\rho'}{\rho_{0}}
$$
写成$\rho'$的等式为
$$
\rho'=-\rho_{0} \frac{\theta'}{\theta_{0}}+\frac{p'}{c_{s}^2}
$$
其中$c_{s}^2=\frac{\gamma \bar{p}}{\rho_{0}}$为声速，对于浮力波运动，有$|\rho_{0}\theta'/\bar{\theta}|\gg|p'/c_{s}^2|$，即密度波动更多源于温度的改变而非压强的改变，于是可得近似关系式
$$
\frac{\theta'}{\bar{\theta}}=-\frac{\rho'}{\rho_{0}}
$$
于是原方程组中的$\rho'$可以消去，将其他所有摄动量均带回原方程组，可得
$$
\begin{gathered}\left(\frac{\partial}{\partial t}+\overline{u}\frac{\partial}{\partial x}\right)u^\prime+\frac{1}{\rho_0}\frac{\partial p^\prime}{\partial x}=0\quad(1)\\\begin{aligned}\left(\frac{\partial}{\partial t}+\overline{u}\frac{\partial}{\partial x}\right)w^{\prime}+\frac{1}{\rho_0}\frac{\partial p^{\prime}}{\partial z}-\frac{\theta^{\prime}}{\overline{\theta}}g=0\quad(2)\end{aligned}\\\frac{\partial u^{\prime}}{\partial x}+\frac{\partial w^{\prime}}{\partial z}=0\quad(3)\\\left(\frac{\partial}{\partial t}+\overline{u}\frac{\partial}{\partial x}\right)\theta^{\prime}+w^{\prime}\frac{d\overline{\theta}}{dz}=0\quad(4)\end{gathered}
$$
接下来消去$p'$项，用$\frac{ \partial  }{ \partial x }$作用于(2)式，减去用$\frac{ \partial  }{ \partial z }$作用于(1)式，可得
$$
\left(\frac{\partial}{\partial t}+\overline{u}\frac{\partial}{\partial x}\right)\left(\frac{\partial w^{\prime}}{\partial x}-\frac{\partial u^{\prime}}{\partial z}\right)-\frac{g}{\bar{\theta}}\frac{\partial\theta^{\prime}}{\partial x}=0
$$
再消去$u',\theta'$，用$\frac{ \partial  }{ \partial x }\left( \frac{ \partial  }{ \partial t }+\bar{u}\frac{ \partial  }{ \partial x } \right)$作用于上式，将(3)(4)式代入，可得
$$
\left(\frac{\partial}{\partial t}+\bar{u}\frac{\partial}{\partial x}\right)^2\left(\frac{\partial^2w^{\prime}}{\partial x^2}+\frac{\partial^2w^{\prime}}{\partial z^2}\right)+N^2\frac{\partial^2w^{\prime}}{\partial x^2}=0
$$
其中$N^2=g \frac{d\ln\bar{\theta}}{dz}$即之前算浮力时出现的因子，这里为了能够解析的解出波动方程，假设$N^2$为常数，这需要扰动是浅扰动。

假设波动解
$$
w'=A\exp[i(kx+mz-\nu t)]
$$
波矢为$\mathbf{\kappa}=(k,m)$，代入原方程可得
$$
(\nu-\bar{u}k)^2(k^2+m^2)-N^2k^2=0
$$
$$
\hat{\nu}=\nu-\bar{u}k=\pm \frac{Nk}{\sqrt{ k^2+m^2 }}=\pm \frac{Nk}{|\mathbf{\kappa}|}
$$
$\hat{\nu}$是本征频率，即相对于平均风场$\bar{u}$的频率，“+”号对应与k同号，等相位面相对于平均风场沿x轴正方向即向东传播，“-”号相反。

如果假设$k>0,m<0$，则等相位面$kx+mz=C$随高度z增大向东倾斜，本征频率取“+”号对应于相对于平均风场向东、向下的相速度，水平和垂直相速度分别为$c_{x}=\hat{\nu}/k,c_{z}=\hat{\nu}/m$，水平和垂直的群速度则为
$$
\begin{aligned}
		&c_{gx}=\frac{ \partial \nu }{ \partial k } =\bar{u}\pm \frac{Nm^2}{(k^2+m^2)^{3/2}}\\&c_{gz}=\frac{ \partial \nu }{ \partial m } =\bar{u}\pm\frac{(-Nkm)}{(k^2+m^2)^{3/2}}
\end{aligned}
$$
可以看到相对于平均风场的群速度的垂直分量和相速度的垂直分量有不同的符号，且相对于平均风场的群速度平行于等相位线，因此群速度垂直于相速度，这表明能量传播平行于波峰波谷，波包的移动方向垂直于相位的传播方向。

讲波动解代入连续性方程可得
$$
ku'+mw'=0
$$
这表明粒子运动的速度方向确实垂直于相速度方向，即重力内波为横波，又由于绝热运动粒子沿等位温面运动，故等位温面与等相位面重合。
![[重力内波示意图.png]]
上图为重力内波的示意图，展示了位温扰动项的波动，以及速度扰动项的波动，用与浅水波类似的方法分析可以发现二者差90度相位。还可以看到等相位面与垂直方向的夹角$\alpha$满足
$$
\cos \alpha=\frac{|k|}{|\mathbf{\kappa}|}
$$
因此本征频率$\hat{\nu}=\pm N\cos \alpha$与前面的简单粒子浮力振荡模型一致。

在大气中，重力内波可能会因为积云对流、地形上的流动以及其他可能向上传播很多高度尺度进入中层大气的过程而在对流层中产生。

## 1.5 旋转分层大气中的线性波
在不考虑旋转的纯重力内波中，粒子沿直线振荡，但在考虑了地球自转效应后，由于科里奥利力垂直于速度方向，使粒子做椭圆形振荡。
### 1.5.1 纯惯性振荡的稳定性
考虑在仅有纬向地转平均流（即仅有$u_{g}$，无$v_{g}$）的情况下的惯性运动，可以用类似于前面粒子浮力模型的方法进行分析（这个假设是合理的，因为稳定性考察的是局部的性质，完全可以将x坐标轴取成实际的地转速度$\mathbf{V}_{g}$的方向，这样就实现了仅有$u_{g}$，无$v_{g}$的条件）

在这种基态流动为纬向地转平均流$u_{g}$的情况下，运动方程为
$$
\begin{aligned}
&\frac{Du}{Dt}=fv=f \frac{Dy}{Dt}\quad(\#)\\& \frac{Dv}{Dt}=f(u_{g}-u)\quad(\#\#)
\end{aligned}
$$
考虑一个本来按基态运动的位于$y=y_{0}$的粒子，给其以沿y方向的扰动$\delta y$，由（#）式积分可得
$$
\delta u=f\delta y
$$
因此有
$$
u(y_{0}+\delta y)=u(y_{0})+f\delta y=u_{g}(y_{0})+f\delta y
$$
$$
u_{g}(y_{0}+\delta y)=u_{g}(y_{0})+\frac{ \partial u_{g} }{ \partial y } \delta y
$$
代入(##)式可得
$$
\frac{Dv}{Dt}=\frac{D^2}{Dt^2}\delta y=-f\left( f-\frac{ \partial u_{g} }{ \partial y }  \right)\delta y=-f\frac{ \partial M }{ \partial y } \delta y
$$
其中定义了绝对动量$M\equiv fy-u_{g}$

上式给出了惯性稳定性条件
$$
f\frac{\partial M}{\partial y}=f\left(f-\frac{\partial u_g}{\partial y}\right)\begin{cases}>0&\mathrm{stable}\\=0&\mathrm{neutral}\\<0&\mathrm{unstable}&\end{cases}
$$
计算可以发现，$\frac{ \partial M }{ \partial y }$是基态流动的绝对涡度。

观测表明，对于热带外的天气尺度系统，流动一般是惯性稳定的，不过在高层急流的反气旋性切变一侧经常会接近惯性中性。大区域的惯性不稳定线性会立刻激发惯性不稳定运动，像对流导致的垂直方向混合一样导致流体在横向混合，从而减少速度的切变最终使绝对涡度×f为正。
### 1.5.2 Rossby波和惯性重力波
当流体既满足惯性稳定，又满足重力稳定时，流体粒子的位移受到旋转和浮力的恢复力，这种恢复力导致的振荡称为惯性重力波，其色散关系也可以用前面的粒子位移模型来启发性地获得。

考虑y，z平面内的一个粒子，其从$(y_{0},z_{0})$位移到了$(y_{0}+\delta y,z_{0}+\delta z)$，设位移矢量与z轴夹角为$\alpha$，则沿位移方向浮力为$-N^2\delta z\cos \alpha$，如果假设$\frac{ \partial u_{g} }{ \partial y }=0$（事实上可以近似认为基态的p在水平面内均匀，因此$\frac{ \partial u_{g} }{ \partial y }$项应比f小），沿位移方向的科里奥利力为$-f^2\delta y\sin \alpha$，设粒子的位移长度为$\delta s$，则粒子的运动方程为
$$
\frac{D^2\delta s}{Dt^2}=-(f\sin \alpha)^2\delta s-(N\cos \alpha)^2\delta s
$$
本征频率为
$$
\nu^2=N^2\cos \alpha+f^2\sin^2\alpha
$$
一般来说$N^2>f^2$，这表明惯性重力波的频率范围为$f\leq|\nu|\leq N$，对于经典的中纬度对流层来说，$N/f \sim 10^2$，因此只有当$\alpha$很接近90度，即粒子位移非常接近水平时，地球自转的效应才会产生影响，此时本质频率$\nu\ll N$，即低频重力波需要地球自转的修正。

以上是启发性的粒子模型，和之前一样我们继续考虑线性化的流体力学方程，这里我们加入之前忽略的旋转项，根据之前的分析，要想旋转项起作用，需要粒子位移非常接近水平，这表明运动的水平尺度远大于垂直尺度，根据Chapter2中的讨论，可以近似认为垂直方向达到流体静力学平衡，另外，再假定基态为$\bar{u}=\bar{v}=\bar{w}=0$，$\bar{p},\bar{\theta}$的定义与前面一样，于是前面的一系列方程可以写为
$$
\begin{gathered}\frac{\partial u^{\prime}}{\partial t}-fv^{\prime}+\frac{1}{\rho_0}\frac{\partial p^{\prime}}{\partial x}=0\\\frac{\partial v^{\prime}}{\partial t}+fu^{\prime}+\frac{1}{\rho_0}\frac{\partial p^{\prime}}{\partial y}=0\\\begin{aligned}\frac{1}{\rho_0}\frac{\partial p^{\prime}}{\partial z}-\frac{\theta^{\prime}}{\overline{\theta}}g=0\end{aligned}\\\frac{\partial u^{\prime}}{\partial x}+\frac{\partial v^{\prime}}{\partial y}+\frac{\partial w^{\prime}}{\partial z}=0\\\frac{\partial\theta^{\prime}}{\partial t}+w^{\prime}\frac{d\overline{\theta}}{dz}=0\end{gathered}
$$
静力平衡方程可以代入到最后一个方程中消去$\theta'$，
$$
\frac{1}{\rho_{0}}\frac{ \partial ^2p' }{ \partial z\partial t } +N^2w'=0
$$
一个一个消元略显麻烦，直接设波动解
$$
\begin{gathered}\mathrm{u'}=\hat{u}\exp i\left(kx+ly+mz-\nu t\right)\\v^{\prime}=\hat{v}\exp i\left(kx+ly+mz-\nu t\right)\\w^{\prime}=\hat{w}\exp i\left(kx+ly+mz-\nu t\right)\\p^{\prime}/\rho_0=\hat{p}\exp i\left(kx+ly+mz-\nu t\right)\end{gathered}
$$
带回原方程组可得
$$
\begin{gathered}
-i\nu \hat{u}-f \hat{v}+ik \hat{p}=0\\-i\nu \hat{v}+f \hat{u}+il \hat{p}=0\\ik \hat{u}+il \hat{v}+im \hat{w}=0\\m\nu \hat{p}+N^2 \hat{w}=0
\end{gathered}
$$
可解得
$$
\begin{gathered}\hat{u}=\left(\nu^{2}-f^{2}\right)^{-1}\left(\nu k+ilf\right)\hat{p}\\\hat{v}=\left(\nu^2-f^2\right)^{-1}\left(\nu l-ikf\right)\hat{p}\\\hat{w}=-\left(\nu m/N^2\right)\hat{p}\end{gathered}
$$
将三者代入连续性方程（即第三个方程）中可得色散关系
$$
		m^2\nu^3-\left[N^2(k^2+l^2)+f^2m^2\right]\nu=0
$$
跟浅水波中一样，也有一个$\nu=0$的解对应静止的Rossby波，该波对应于垂直速度$w'=0$，水平速度为地转风。对于本质频率的非零解，为
$$
\nu^2=f^2+N^2(k^2+l^2)m^{-2}
$$
由于水平尺度远大于垂直尺度，即水平波长远大于垂直波长，因此有$(k^{2}+m^{2})m^{-2}\ll 1$，因此波的频率满足$|f|<|\nu|\ll N$，且如果令
$$
\sin^2\alpha \to 1,\cos^2\alpha=(k^2+l^2)m^{-2}
$$
则色散关系与前面启发性的粒子模型的结果一致。而且如果令$\bar{h}=\frac{N^{2}}{gm^{2}}$为浅水波的平均水深，则色散关系也与浅水波的色散关系一致。可见我们的解非常正确且神奇！

如果将x轴取成沿着波在水平面内的传播方向，即$l=0$（这一步等价于前面粒子模型中通过转坐标轴使$v_{g}=0$），分别计算波水平群速度和垂直群速度
$$
c_{gx}=\frac{ \partial \nu }{ \partial k } =\frac{ \partial \nu^2 }{ \partial k } /2\nu=\frac{N^2k}{\nu m^2}
$$
$$
c_{gz}=\frac{ \partial \nu^2 }{ \partial m } /2\nu=-\frac{N^2k^2}{\nu m^3}
$$
于是二者之比为
$$
|c_{gz}/c_{gx}|=|k/m|=\frac{\sqrt{ \nu^2-f^2 }}{N}
$$
可见群速度和纯重力内波一样也是平行于等相位面，即垂直于相速度方向。在$l=0$的情况下还能得到$\hat{v}=-if \hat{u}/\nu$，因此$\hat{v}$的相位比$\hat{u}$落后$\frac{\pi}{2}$，可得
$$
u'=\hat{u}\cos(kx+mz-\nu t),\quad v'=\frac{\hat{u}f}{\nu}\sin(kx+mz-\nu t)
$$
可见水平面内粒子沿椭圆轨迹运动，其水平速度矢量进行反气旋性转动（北顺南逆），当$m<0$时水平速度矢量随高度上升进行反气旋性转动，水平风场随时间和高度z的**反气旋性转动**是利用大气数据识别惯性重力波的一种**基本手段**。上述特征在下图中展现
![[惯性重力波示意图.png]]
上图展示了能量向上传播（$m<0,\quad \nu>0$）的北半球($f>0$)惯性重力波在垂直截面下的速度、压强、温度的波动，从图中可以看到这些波动量的相位关系。
### 1.5.3 涡度和位涡
惯性重力波对应的线性涡度方程可以用$\frac{ \partial  }{ \partial x }$作用于y方向运动方程减去$\frac{ \partial  }{ \partial y }$作用于x方向运动方程得到
$$
\frac{ \partial \zeta' }{ \partial t } =-f\left( \frac{ \partial u' }{ \partial x } +\frac{ \partial v' }{ \partial y }  \right)=f\frac{ \partial w' }{ \partial z }
$$
再利用绝热下的热力学能量方程
$$
w'=-\frac{1}{\left( \frac{d\bar{\theta}}{dz} \right)}\frac{ \partial \theta' }{ \partial t }
$$
于是可以令线性位涡为
$$
\Pi'=\zeta'+f\frac{ \partial  }{ \partial z } \left( \frac{\theta'}{d\bar{\theta}/dz} \right)
$$
则有位涡守恒
$$
\frac{ \partial \Pi' }{ \partial t }=0
$$
有点类似于浅水波。还可以进一步利用波动解给出涡度$\zeta'$与$\hat{p}$的关系，假设有$l=0$，对于$\nu=0$的静态Rossby波，有
$$
\zeta'=ikv'=-\frac{k^2}{f}p'
$$
对于惯性重力波有
$$
\zeta'=\frac{k^2f}{\nu^2-f^2}p'=\frac{fm^2}{N^2}p'
$$
于是当$N$为常数时，由于
$$
\frac{ \partial \theta' }{ \partial z } =\frac{\bar{\theta}}{g}\frac{ \partial ^2p' }{ \partial z^2 } =-\frac{m^2\bar{\theta}}{g}p'
$$
代入线性位涡表达式可得惯性重力波的位涡满足
$$
\Pi'=0
$$
因此线性位涡方程可以从动力学中把惯性重力波过滤掉（？）
## 1.6 向地转平衡的调节
现在我们来讨论，当初始状态没有达到地转平衡时，系统会经过怎样一个过程而达到地转平衡，这个过程称为调节（adjustment）过程

为了简化讨论，我们研究之前介绍的浅水波方程
$$
\begin{aligned}&\frac{\partial u^{\prime}}{\partial t}=-g\frac{\partial h^{\prime}}{\partial x}+fv^{\prime}\\&\frac{\partial v^{\prime}}{\partial t}=-g\frac{\partial h^\prime}{\partial y}-fu^\prime\\&\frac{\partial h^{\prime}}{\partial t}=-\bar{h}\left(\frac{\partial u^{\prime}}{\partial x}+\frac{\partial v^{\prime}}{\partial y}\right)\end{aligned}
$$
现在我们需要解决的问题是从任意初始状态开始，在长时间极限下，上面浅水波方程的解（该问题由Rossby解决，因此被称为Rossby调节问题）。解决这个问题的关键在于浅水位涡方程
$$
\frac{ \partial Q' }{ \partial t } =0
$$
其中浅水位涡为
$$
Q'(x,y,t)=\frac{\zeta'}{f}-\frac{h'}{\bar{h}}=Const
$$
其重要性在于其末态仅与初态相关$Q(x,y,t)=Q(x,y,0)$，不需要解时间依赖方程。

为了简化，我们考虑一个理想化的浅水系统，其初态为
$$
u',v'=0;\quad h'=-h_{0} sgn(x)
$$
其中$sgn(x)$即单位阶跃函数，在x正半轴为1，负半轴为-1，利用浅水位涡守恒可得
$$
Q=\frac{\zeta'}{f}-\frac{h'}{\bar{h}}=\frac{h_{0}}{\bar{h}}sgn(x)
$$
利用上式可以将涡度$\zeta'$用$h'$表示，下面回到之前浅水波关于$h'$的偏微分方程
$$
\begin{aligned}
\frac{ \partial ^2h' }{ \partial t^2 } &=g \bar{h}\nabla_{h}^2h'-f \bar{h}\left( \frac{ \partial v' }{ \partial x } -\frac{ \partial u' }{ \partial y }  \right)\\&=g \bar{h}\nabla_{h}^2h'-f \bar{h}\zeta'
\end{aligned}
$$
代入$\zeta'$可得
$$
\frac{\partial^2h^{\prime}}{\partial t^2}-c^2\left(\frac{\partial^2h^{\prime}}{\partial x^2}+\frac{\partial^2h^{\prime}}{\partial y^2}\right)+f^2h^{\prime}=-f^2h_0\mathrm{sgn}\left(x\right)
$$
其中$c=\sqrt{ g \bar{h} }$

由于初态$h'$与y无关，且上面的偏微分方程各项系数和非齐次项均与y无关，因此$h'$始终与y无关，于是在长时间极限，即稳态下，$h'$对t的偏导数为0，可得均
$$
-c^2\frac{d^2h'}{dx^2}+f^2h'=-f^2h_{0}sgn(x)
$$
利用$h'$在$x=0$处连续且$\frac{ \partial h' }{ \partial x }$在$x=0$处连续（后者是因为速度$v'$肯定在$x=0$处连续），微分方程解为
$$
\frac{h^{\prime}}{h_0}=\begin{Bmatrix}-1+\exp\left(-x/\lambda_R\right)&\mathrm{~for~}x>0\\+1-\exp\left(+x/\lambda_R\right)&\mathrm{~for~}x<0\end{Bmatrix}
$$
其中$\lambda_{R}=\frac{\sqrt{ g\bar{h} }}{f}$是Rossby变形半径，其可以当作高度场调节至地转平衡时的水平长度尺度，将上式代入浅水波方程，对于长时间极限下的稳态速度场，可得
$$
u'=0,\quad and \quad v'=\frac{g}{f}\frac{ \partial h' }{ \partial x } =-\frac{gh_{0}}{f\lambda_{R}}\exp(-|x|/\lambda_{R})
$$
该稳态速度场是无散地转风。速度场合高度场的稳态解如下图所示
![[稳态场示意图.png]]
注意，考虑长时间极限稳态时，我们其实只需要把浅水波方程中所有对时间的偏导数都取成0，就可以得到无散的地转风
$$
fu'=-g\frac{ \partial h' }{ \partial y } ,\quad fv'=g\frac{ \partial h' }{ \partial x } ,\quad \frac{ \partial u' }{ \partial x } +\frac{ \partial v' }{ \partial y } =0
$$
但无法得到末态具体的扰动高度场$h'$，这是因为我们没有充分利用已知的扰动高度场初始条件。通过守恒量浅水位涡，可以将初态和长时间极限下的稳态联系起来，从而得到具体的稳态扰动高度场。

### 1.6.1.1 地转平衡的调节过程
尽管求出末态无需解时间依赖的方程，但要求出具体的调节过程就需要解时间依赖方程，比如根据初始条件解出$h'$的偏微分方程，这有些复杂，我们仅讨论调节过程中能量的流向，这只需要计算初态和末态能量的改变量。

单位水平面积的势能为
$$
\int_{0}^{h'}\rho gzdz=\rho gh'^2/2
$$
因此y方向单位长度在调节过程中释放的势能为
$$
\begin{aligned}
	\int_{-\infty}^{+\infty} &\frac{\rho gh_{0}^2}{2}dx-\int_{-\infty}^{+\infty} \frac{\rho gh'^2}{2}dx\\&=2\int_{0}^{+\infty} \frac{\rho gh_{0}^2}{2}[1-(1-e^{-x/\lambda_{R}})^2]dx\\&=\frac{3}{2}\rho gh_{0}^2\lambda_{R}
\end{aligned}
$$
对于无地球自转的情况来说，$\lambda_{R}\to \infty$，势能释放无限，这是因为稳态时表面为平面，所有可以释放的地势能均释放掉了，由重力波携带着传播出去。

对于有地球自转的情况来说，上式的势能一部分转换为了流体运动的动能，还有一部分动能被重力波携带着传播出去，传播到了考察区域以外的地方。y方向单位长度的动能增加为
$$
2\int_{0}^{+\infty}\rho \bar{h} \frac{v'^2}{2}dx=\rho \bar{h}\left( \frac{gh_{0}}{f\lambda_{R}} \right)^2\int_{0}^{+\infty}e^{-2x/\lambda_{R}}dx=\frac{1}{2}\rho gh_{0}^2\lambda_{R}
$$
因此释放的有限势能有$\frac{1}{3}$变成了稳定地转流的动能，剩下$\frac{2}{3}$以惯性重力波的形式辐射出去。

以上的讨论揭示了下面几点：
1. 我们难以直接得出旋转流体的势能，因为$|x|\to \infty$时，$h'$有限，只能得到从初态到达到地转平衡态过程中释放的势能。
2. 位涡守恒使得我们可以具体的给出稳态下的地转调节速度以及高度场，而不用对时间积分。
3. 稳态解的长度尺度为Rossby半径$\lambda_{R}$

## 1.7 Rossby波
对于大尺度气象过程最重要的波就是Rossby波（行星波），其主要讨论的是沿着地球纬向流动的波。在正压浅水波的情况下（由于正压，故深度为常数，水平无散），Rossby波是绝对涡度守恒的运动，其产生是因为科里奥利参数$f$随纬度的变化，更一般性的，对于斜压大气，Rossby波是位涡守恒的运动，其产生是因为位涡梯度的存在。

科里奥利参数随纬度的变化能在纬度$\phi_{0}$处近似展开为
$$
f=f_{0}+\beta y
$$
其中$\beta \equiv\left( \frac{df}{dy} \right)_{\phi_{0}}=\frac{df}{d\phi}/ \frac{dy}{d\phi}=2\Omega \cos \phi/a$，这个近似被称为中纬度$\beta$-平面近似。Rossby波可以通过一个简单的模型来定性的理解，考虑一个环绕纬线一圈的闭合链条流体粒子，其绝对涡度$\eta=\zeta+f$守恒，假设初始时刻$t_{0}$有$\zeta=0$，现在假设时刻$t_{1}$时流体粒子发生了沿经向的位移$\delta y$，则根据绝对涡度守恒有
$$
(\zeta+f)_{t_{1}}=f_{t_{0}}
$$
$$
\zeta_{t_{1}}=f_{t_{0}}-f_{t_{1}}=-\beta \delta y
$$
因此初始的位移扰动会导致流体产生扰动涡度场，而扰动的涡度场又会导致经向速度场$v$，使涡度极大值西侧的流体粒子向南位移，涡度极小值西侧的流体粒子向北位移，示意图如下图
![[Rossby波定性分析示意图.png]]
其中粗波浪线表示初始的经向扰动，细波浪线表示由涡度场激发的经向速度场使流体粒子发生流动后的新经向扰动，可见流体粒子在其平衡纬度上下振动，而涡度极大极小的模式在向西传播，这种涡度场的详细传播就是Rossby波。就像位温梯度作为流体垂直位移的恢复力一样，绝对涡度（科里奥利参数$f$）的经向梯度作为经向位移的恢复力。

波动向西传播的速度$c$可以用一个简单的例子来计算，设位移满足$\delta y=a\sin k(x-ct)$，则经向速度场为$v=\frac{D}{Dt}\delta y=-kca\cos k(x-ct)$，于是可得涡度场为
$$
\zeta=\frac{ \partial v }{ \partial x } =k^2ca\sin k(x-ct)
$$
将上面的波动形式代入前面绝对涡度守恒导出的公式$\zeta=-\beta \delta y$可得
$$
k^2ca\sin k(x-ct)=-\beta a\sin k(x-ct)
$$
因此波速为
$$
c=-\frac{\beta}{k^2}
$$
相速度相对于平均流向西，反比于纬向波数的平方。
### 1.7.1 自由正压Rossby波
正压Rossby波的色散关系可以利用线性化的正压涡度方程正式地推导出来。根据Chapter4，正压涡度方程为
$$
\frac{D_{h}}{Dt}(\zeta+f)=\left( \frac{ \partial  }{ \partial t } +u\frac{ \partial  }{ \partial x } +v\frac{ \partial  }{ \partial y }  \right)\zeta+\beta v=0
$$
将其线性化，取基态为常数平均纬向速度，则有
$$
u=\overline{u}+u^{\prime},\quad v=v^{\prime},\quad\zeta=\partial v^{\prime}/\partial x-\partial u^{\prime}/\partial y=\zeta^{\prime}
$$
由于正压浅水波无散，故由Chapter4可以定义一个流函数$\psi'$，满足
$$
u^{\prime}=-\partial\psi^{\prime}/\partial y,\quad v^{\prime}=\partial\psi^{\prime}/\partial x
$$
则涡度为$\zeta'=\nabla^2\psi'$，代入正压涡度方程可将其线性化为
$$
\left(\frac{\partial}{\partial t}+\overline{u}\frac{\partial}{\partial x}\right)\nabla^2\psi^{\prime}+\beta\frac{\partial\psi^{\prime}}{\partial x}=0
$$
设波动解
$$
\psi'=\Psi\exp[i(kx+ly-\nu t)]
$$
代入可得
$$
(-\nu+k \bar{u})(-k^2-l^2)+k\beta=0
$$
$$
\nu=k \bar{u}-\frac{\beta k}{\mathbf{K}^2}
$$
其中$\mathbf{K}^2=k^2+l^2$，该波动沿x方向的相速度为
$$
c=\frac{\nu}{k}=\bar{u}-\frac{\beta}{\mathbf{K}^2}
$$
当平均风场消失且$l \to0$时上式退化至前面定性模型得到的相速度，可见Rossby波纬向相速度相对于平均风场向西，且随波长的增大而快速增大。对于经典的中纬度天气尺度下的扰动，Rossby波相对于东流的纬向平均流的速度比起平均风场要小，因此天气尺度下Rossby波通常会向东移动，但相速度小于平均风场的速度。当波长很长时，Rossby波的相速度可能会为0，此时Rossby波静止，达成这一条件需要满足
$$
\mathbf{K}^2=\frac{\beta}{\bar{u}}\equiv \mathbf{K}_{s}^2
$$
在下一子节里会讨论这个条件的重要意义。

与相速度不同，Rossby纬向的群速度相对于平均流可能向东也可能向西，取决于纬向和经向波数的大小关系
$$
c_{gx}=\frac{ \partial \nu }{ \partial k } =\bar{u}-\beta\frac{ l^2-k^2}{(k^2+l^2)^2}
$$
对于静态的Rossby波，其群速度为
$$
c_{gx}=\beta  \frac{1}{k^2+l^2}-\beta \frac{l^2-k^2}{(k^2+l^2)^2}=\frac{2\beta k^2}{(k^2+l^2)^2}
$$
相对于地面向东。天气尺度的Rossby波也倾向于有相对于地面向东的群速度，同时平均风场一般比Rossby波的相速度要大，因此相速度也相对于地面向东，但小于群速度，这表明新的扰动会在现有扰动的下游（东侧）形成，因为能量传播的比波本身要更快（能量先行，波形后到），在天气预报中这一效应的考虑非常重要。

我们也可以用更严格的分析来考察自由Rossby波，而不是简单的正压浅水波模型，但其数学过程有些复杂，其色散关系结果也与我们的浅水波模型得到的类似。

行星尺度的自由振荡尽管在仔细的观测下可以探测到，但它们的振幅往往很微小，因为其受力很小。
### 1.7.2 地形受迫Rossby波
自由传播的Rossby模式在大气中相当弱，但受迫的静止Rossby模式在理解行星尺度大气环流模式时非常重要，这种模式可能会受迫于依赖于经度的加热源或在地形上的流动，对于北半球热带以外的环流，特别重要的是由在落基山脉和喜马拉雅山脉上流过的静态Rossby模式，其为地形Rossby波，在4.3中用位涡守恒定性描述了大气越过山脉的流线。

作为一个最简单的地形Rossby波的可能的动力学模型，我们采用浅水位涡$\Pi=\frac{\zeta+f}{h}$的一个近似，假设上边界为固定高度$H$，下边界为变化的高度$h_{T}(x,y)$，其中$|h_{T}|\ll H$，再假设$\zeta$可以近似为地转风对应的涡度$\zeta_{g}$，，且假设$|\zeta_{g}|\ll f_{0}$，浅水位涡方程可以近似为
$$
\frac{D_{h}}{Dt}\left( \frac{\zeta_{g}+f}{h} \right)=\frac{1}{h}\left( \frac{ \partial  }{ \partial t } +\mathbf{V}\cdot \nabla \right)(\zeta_{g}+f)-\frac{1}{h^2}(\zeta_{g}+f) \frac{Dh}{Dt}=0
$$
化简+近似可得
$$
H\left( \frac{ \partial  }{ \partial t } +\mathbf{V}\cdot \nabla \right)(\zeta_{g}+f)=-f_{0} \frac{Dh_{T}}{Dt}
$$
用和自由Rossby波一样的基态假设可将上式线性化
$$
\left( \frac{ \partial  }{ \partial t } +\bar{u}\frac{ \partial  }{ \partial x }  \right)\zeta_{g}'+\beta v_{g}'=-\frac{f_{0}}{H}\bar{u}\frac{ \partial h_{T} }{ \partial x }
$$
现在我们在底部边界为正弦型函数的情况下考察上式的解，我们设底部为
$$
	h_{T}(x,y)=\mathrm{Re}[h_{0}\exp(ikx)]\cos ly
$$
与前面类似，用流函数来表示地转风和涡度，现考虑稳态，即场与时间无关，
$$
\psi(x,y)=\mathrm{Re}[\psi_{0}\exp(ikx)]\cos ly
$$
代入可得复振幅的关系为
$$
\psi_{0}=\frac{f_{0}h_{0}}{H(K^2-K_{s}^2)}
$$
对于长波长，$K<K_{s}$，地形涡度源主要由行星涡度$f$的经向变化即$\beta$效应来平衡；对于短波长，$K>K_{s}$，地形涡度元主要由相对涡度的纬向变化平衡。当$K=K_{s}$，即自由Rossby模式为静态Rossby波时，$\psi_{0}\to \infty$，即发生了共振现象，可以通过在上面的方程中引入一项边界层拖移导致的Ekman泵项（即给相对涡度引入线性阻尼）来消除共振发散现象，此时涡度方程写为
$$
\left( \frac{ \partial  }{ \partial t } +\bar{u}\frac{ \partial  }{ \partial x }  \right)\zeta_{g}'+\beta v_{g}'+r\zeta_{g}'=-\frac{f_{0}}{H}\bar{u}\frac{ \partial h_{T} }{ \partial x }
$$
则类似的可得其稳态解振幅关系为
$$
\psi_{0}=\frac{f_{0}h_{0}}{H(K^2-K_{s}^2-i\epsilon)}
$$
其中$\epsilon=\frac{rK^2}{k\bar{u}}$，此时当满足$K=K_{s}$时，$\psi_{0}$不会趋于$\infty$，但确实会达到最大值，此时x方向的流函数比地形受迫源领先了$\frac{\pi}{2}$的相位，因此流函数的槽（最小值）位于地形的山峰东侧四分之一个波长处（从物理上可以理解为边界层摩擦减缓了气块的响应速度，山脉脊激发的上升运动要到更下游处才能发展为低压槽）。这里之所以把流函数的最小值和低压槽等同，是因为根据地转风的受力方程可得
$$
\psi=\frac{\Phi}{f}
$$
因此流函数的槽等价于地势的槽等价于压强的槽。

上面讨论的是特殊的正弦型地形，而对于一般的地形，可以用傅里叶展开为正弦地形的组合。
## 1.8 总结
本章引入了重要的线性摄动理论，可以将非线性的复杂偏微分方程通过合理近似转化为线性的偏微分方程。利用该方法研究了一系列波动，从最简单的一维声波和浅水波开始，逐渐深入讨论了纯重力内波、考虑了地球自转的惯性重力波、由绝对涡度守恒引出的自由Rossby波以及由浅水位涡守恒引出的地形受迫Rossby波，每次都是先从一个启发性的粒子模型开始，反映了这些波的物理本质都是某种流体的稳定性（如浮力对应垂直稳定性、惯性力对应水平稳定性、f的梯度对应经向稳定性），而后再通过相应的物理学方程组加上线性摄动法加波动解假设解出相应的波动解和色散关系，接着就可以讨论它们的相速度、群速度、相位关系、涡度等等。

