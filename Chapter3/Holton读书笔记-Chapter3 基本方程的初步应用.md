---
zhihu-title: Holton读书笔记-Chapter3 基本方程的初步应用
zhihu-topics: 大气科学
zhihu-link: https://zhuanlan.zhihu.com/p/1930180220087939862
zhihu-created-at: 2025-07-28 22:48
zhihu-updated-at: 2025-08-06 09:14
---
```toc
```
[[Holton读书笔记-Chapter1 基本介绍]](https://zhuanlan.zhihu.com/p/1928544934622924935  "card")
[[Holton读书笔记-Chapter2 基础守恒律]](https://zhuanlan.zhihu.com/p/1929242659588912307 "card")
# 1 Chapter3
除了在Chapter2中提到的地转风的近似以外，还有很多其他的近似在对天气系统的分析中有用，本章对这些近似进行讨论。
## 1.1 压强坐标下的基本方程
在大气的分析中，正如Chapter1中所说，把垂直坐标从z换成p可以带来很多方便，所以我们先讨论一下在压强坐标下原来的那些基本方程应该怎么写。
### 1.1.1 水平动量方程
对于水平速度 $\mathbf{V}=(u,v)$，考虑水平方向的受力方程，近似后为
$$
\quad\frac{D\mathbf{V}}{Dt}+f\mathbf{k}\times\mathbf{V}=-\frac{1}{\rho}\nabla p
$$
由Chapter1的结果，在p坐标下写为
$$
\frac{D\mathbf{V}}{Dt}+f\mathbf{k}\times\mathbf{V}=-\nabla_p\Phi
$$
其中，原来的 $\frac{D}{Dt}$也需要变形成p坐标下的形式（这里不知为何Holton写的极为简单，我觉得这不是一件很trivial的事情）。

在z坐标下，
$$
\frac{D}{Dt}=\left(\frac{ \partial }{ \partial t }\right)_{z}+u(\frac{ \partial }{ \partial x })_{z} + v (\frac{ \partial  }{ \partial y })z+w\frac{ \partial }{ \partial z }
$$
由chapter1中所证明的（通过类似的方法即可证明一系列的等式），
$$
		\left( \frac{ \partial  }{ \partial x }  \right)_{z}=\left( \frac{ \partial  }{ \partial x }  \right)_{p}+\left(\frac{ \partial p }{ \partial x }\right) _{z}\frac{\partial}{ \partial p }  
$$
同理
$$
		\left( \frac{ \partial  }{ \partial y }  \right)_{z}=\left( \frac{ \partial  }{ \partial y }  \right)_{p}+\left(\frac{ \partial p }{ \partial y }\right) _{z}\frac{\partial}{ \partial p }  
$$
$$
						\left( \frac{ \partial  }{ \partial t }  \right)_{z}=\left( \frac{ \partial  }{ \partial t }  \right)_{p}+\left( \frac{ \partial p }{ \partial t }  \right)_{z}\frac{ \partial  }{ \partial p } 
$$

利用简单链式法则可得
$$
							\frac{ \partial  }{ \partial z } =\frac{ \partial p }{ \partial z } \frac{ \partial  }{ \partial p } 
$$
全部代入可得在p坐标下
$$
\begin{aligned}
\frac{D}{Dt}&=\left(\frac{ \partial \ }{ \partial t } \right)_{p}+u\left( \frac{ \partial  }{ \partial x }  \right)_{p}+v\left( \frac{ \partial  }{ \partial y }  \right)_{p}+\left[ \frac{ \partial p }{ \partial t } +u\frac{ \partial p }{ \partial x } +v\frac{ \partial p }{ \partial y } +w\frac{ \partial p }{ \partial z }  \right] \frac{ \partial  }{ \partial p }\\
&=\left(\frac{ \partial \ }{ \partial t } \right)_{p}+u\left( \frac{ \partial  }{ \partial x }  \right)_{p}+v\left( \frac{ \partial  }{ \partial y }  \right)_{p}+\frac{Dp}{Dt}\frac{ \partial  }{ \partial p } 
\end{aligned}
$$
如果直接令p坐标下 $w=\frac{Dp}{Dt}$（后面在p坐标下w均表示这个含义），则可以写成和z坐标下类似的形式
$$
\frac{D}{Dt}=\frac{\partial}{\partial t}+u\frac{\partial}{\partial x}+v\frac{\partial}{\partial y}+\omega\frac{\partial}{\partial p}
$$


利用p坐标下的运动方程可以得出p坐标下的地转关系：
$$
	f \mathbf{k}\times \mathbf{V_{g}}=-\nabla_{p} \Phi
$$
利用 $\mathbf{k}\times (\mathbf{k}\times \mathbf{V_{g}})=-\mathbf{V_{g}}$，可得
$$
\mathbf{V_{g}}=\frac{1}{f}\mathbf{k}\times \nabla_{p}\Phi
$$
从上式可见p坐标的一个重要优势，即$\rho$消失，因此在给定水平地势的梯度以后，同一p坐标对应同一地转风速分布，而在z坐标下，给定水平压强梯度可能对应不同地转风速分布。

当$f$是常数时，还可以推知地转风的水平散度为0：
$$
\nabla_{p}\cdot \mathbf{V}_{g}=0
$$
## 1.2 连续性方程
连续性方程直接从z坐标下的情况开始推导会遇到麻烦，更方便的做法是直接在p坐标框架里从质量守恒开始推导。

考虑z坐标框架下的体积微元 $\delta V=\delta x \delta y \delta z$ ，在p坐标框架下，由雅可比行列式
$$
				\frac{ \partial (x,y,z) }{ \partial (x,y,p) }=\frac{ \partial z }{ \partial p } =-\frac{1}{\rho g} 
$$
故p坐标下的体积微元为 $\delta V=-\frac{\delta x \delta y \delta p}{\rho g}$ ，所以质量为 $\delta M=\rho \delta V=-\frac{\delta x\delta y\delta p}{g}$ ，由质量守恒
$$
\frac{1}{\delta M}\frac{D}{Dt}(\delta M)=\frac{g}{\delta x\delta y\delta p}\frac{D}{Dt}\left(\frac{\delta x\delta y\delta p}{g}\right)=0
$$
展开可得
$$
\frac{1}{\delta x}\delta\left(\frac{Dx}{Dt}\right)+\frac{1}{\delta y}\delta\left(\frac{Dy}{Dt}\right)+\frac{1}{\delta p}\delta\left(\frac{Dp}{Dt}\right)=0
$$
取微元大小趋于0的极限可得
$$
\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right)_p+\frac{\partial\omega}{\partial p}=0
$$
在p坐标下的连续性方程不包含$\rho$和时间导数，这也是取p坐标的一大优势。
## 1.3 热力学能量方程
z坐标下热一定律写为
$$
c_{p} \frac{DT}{Dt}-\alpha\frac{ Dp}{Dt}=J
$$
在p坐标下可以写为
$$
c_p\left(\frac{\partial T}{\partial t}+u\frac{\partial T}{\partial x}+v\frac{\partial T}{\partial y}+\omega\frac{\partial T}{\partial p}\right)-\alpha\omega=J
$$
可以定义
$$
S_p\equiv\frac{RT}{c_pp}-\frac{\partial T}{\partial p}=-\frac{T}{\theta}\frac{\partial\theta}{\partial p}
$$
代入可得
$$
\left(\frac{\partial T}{\partial t}+u\frac{\partial T}{\partial x}+v\frac{\partial T}{\partial y}\right)-S_p\omega=\frac{J}{c_p}
$$
$S_{p}$的含义在于其正负可以代表稳定性：
$$
S_{p}=\frac{\Gamma_{d}-\Gamma}{\rho g}
$$
## 1.4 平衡流动
为了对大气中力的水平平衡有一个定量的认识，我们假设流动处于稳态（时间无关）且速度没有垂直分量。另外，将水平动量方程展开至自然坐标系的分量形式会带来方便。
### 1.4.1 自然坐标系
定义沿速度方向的单位向量$\mathbf{t}$、垂直于速度方向的法向单位向量$\mathbf{n}$(定义为速度左侧)、垂直于水平面的单位向量$\mathbf{k}$。

在自然坐标系下， $\mathbf{V}=V\mathbf{t}$ ，$V=\frac{Ds}{Dt}$，定义运动的曲率半径R若曲率中心在速度左侧为正，右侧为负，则加速度为
$$
\frac{D\mathbf{V}}{Dt}=\mathbf{t} \frac{DV}{Dt}+V \frac{D\mathbf{t}}{Dt}=\mathbf{t} \frac{DV}{Dt}+\mathbf{n} \frac{V^2}{R}
$$
科里奥利力写为
$$
-f\mathbf{k}\times \mathbf{V}=-fV \mathbf{n}
$$
梯度力写为(结合梯度的定义，分别考虑沿s和n方向的小位移容易得出)
$$
					-\nabla_{p}\Phi=-\left( \mathbf{t}\frac{ \partial \Phi }{ \partial s } +\mathbf{n} \frac{ \partial \Phi }{ \partial n }  \right)
$$
于是可以写出自然坐标系下的动量方程
$$
		\frac{DV}{Dt}=-\frac{ \partial \Phi }{ \partial s } 
$$
$$
		\frac{V^2}{R}+fV=-\frac{ \partial \Phi }{ \partial n } \quad(\#)
$$
当流体运动沿着等$\Phi$线（等高线）时，$\frac{ \partial \Phi }{ \partial s }=0$ ，因此水平速度大小是常数，这种流动称为梯度流（gradient flow）。如果更进一步的，$\Phi$沿运动法向的导数也为常数，则由（#）式该运动的轨迹曲率半径为一常数。在】这种情况下的流动可以根据（#）式三项的不同分成几类。
##### 1.4.1.1.1 讨论
注意这里提到的等高线的含义。由于对$\Phi$的梯度都是在等p的前提下求的，在p坐标下$z=z(x,y,p)$，因此在xyz坐标架中等p对应于一个面，而这里的等高线指的是在这个等p面上划出的一条等z线。
#### 1.4.1.2 地转流（Geostrophic Flow）
沿着一条平行于等高线的直线流动被称为地转运动（Geostrophic motion），在这种情况下有
$$
		V=-\frac{1}{f}\frac{ \partial \Phi }{ \partial n }=V_{g} 
$$
压强梯度力和科里奥利力严格平衡，V等于地转风$V_{g}$（这也许就是其名字的来历）
### 1.4.2 惯性流
如果整个等p面处处等z，则水平压强梯度力完全消失，完全由科里奥利力指导运动，
$$
\frac{V^2}{R}+fV=0
$$
此时$R=-\frac{V}{f}$为常数，气团做圆周运动，且为反气旋运动（即在北半球顺时针，南半球逆时针）周期为
$$
P=\left|\frac{2\pi R}{V}\right|=\frac{2\pi}{|f|}=\frac{\pi}{\Omega \sin\phi}
$$
由于这个运动完全由科里奥利力决定，没有外力存在，因此称为惯性振动。

在大气中，运动几乎总是由气压梯度力产生和维持，这种均匀气压的情况很罕见。但在海洋中，运动主要由瞬时的表面吹过的风造成，因此很多运动接近于惯性流运动。
#### 1.4.2.1 旋转平衡流（Cyclostrophic Flow）
当水平扰动的尺度很小，即特征长度L很小时（如龙卷风），由Rossby数$\frac{U}{fL}$很大，即科里奥利力相对于加速度很小，于是科里奥利力可忽略，n方向运动方程写为
$$
		\frac{V^2}{R}=-\frac{ \partial \Phi }{ \partial n } 
$$
得到旋转平衡风
$$
			V=\left( -R \frac{ \partial \Phi }{ \partial n }  \right)^{1/2}
$$
可见旋转平衡风可能为顺气旋也可能为逆气旋，然而观测到北半球的大多数龙卷风都是逆时针的（即顺气旋），这里面的神秘原因将在Chpater9讲述。
#### 1.4.2.2 梯度风近似
比起地转流，更贴近实际的近似是仅考虑流动沿着等高线，但不要求等p面等z，对应于V为常数的流动，称为梯度流，这时的运动由n方向方程的三项共同决定。

利用n方向运动方程可得
$$
\begin{aligned}
V&=-\frac{fR}{2}\pm\left(\frac{f^2R^2}{4}-R\frac{\partial\Phi}{\partial n}\right)^{1/2}\\&=-\frac{fR}{2}\pm\left(\frac{f^2R^2}{4}+fRV_g\right)^{1/2}
\end{aligned}
$$
上式中不是所有的根都有物理意义的，需要V为非负实数，满足该条件的情况在北半球分为四种，如下图

![[梯度风四种情况.png]]
>大量概念预警！

这里正常（regular）与异常（anomalous）的定义是根据气团相对于轴的绝对角动量的正负(包括相对运动角动量和地球自转带来的角动量），
$$
L=VR+\frac{fR^2}{2}
$$
对于北半球来说，正常对应于L>0，异常对应于L<0（南半球相反）。前面提到的气旋就是regular low，因此在北半球逆时针，反气旋就是regular high，因此在北半球顺时针。

还可以看到，对于高气压中心，无论是正常还是异常，都对应R<0且$V_{g}>0$，为保证V是实数都有不等式
$$
		0<-\frac{ \partial \Phi }{ \partial n } =fV_{g}<- \frac{f^2R}{4}
$$
因此对于高压中心，当趋于中心，即R->0时，气压梯度力会趋于0，这就是为何高压中心处的压强变化比起低压中心相对平缓、风相对柔和。

除了anomalous low的情况以外，压强梯度力均与科里奥利力反向，这种称为baric（不知道怎么翻译），而同向的anomalous low情况称为antibaric，在antibaric的情况下，地转风$V_{g}$为负，即压强梯度力指向速度右侧，此时$V_{g}$与速度反向，显然无法作为实际风的近似了。

前文所提到的气旋对应于regular low，反气旋对应于regular high。从图中可以发现，无论是南北半球，气旋要求离心力与科里奥利力同向，即Rf>0，反气旋要求二者反向，即Rf<0。

由n方向运动方程可得
$$
\frac{V_{g}}{V}=1+\frac{V}{fR}
$$
对于气旋流，$V_{g}>V$；对于反气旋流，$V_{g}<V$。在中纬度地区，梯度风V与地转风的偏差不超过10%-20%，但在热带地区，偏差较大，因此必须使用梯度风。

## 1.5 轨迹和流线
对于梯度风来说，前面运动方程中给出的曲率半径为R的路径是一个空气块的轨迹，在应用中，常常用等高线的曲率半径（容易从天气图里读出）代替轨迹的曲率半径，但事实上其等高线对应的是梯度风的流线，当等高线发生变化时，二者并不等价。

区分轨迹（trajectory）和流线（streamline）是很重要的，流线是当下**瞬时**速度场的快照，即流线上处处与该处速度平行，决定流线的方程为$\frac{dy}{dx}=\frac{v(s,y,t_{0})}{u(x,y,t_{0})}$；而轨迹是**一定时间内**一个流体块的运动轨迹，决定轨迹的方程是$\frac{Ds}{Dt}=V(x,y,t)$，后面会证明只有在稳态运动下，流线和轨迹才是重合的，而事实上天气系统的运动不是稳态，整个系统往往会以与自身风速大概同数量级的速度移动，因此用流线的曲率半径代替轨迹的曲率半径会造成错误，下面研究两个曲率半径之间的关系。

设$\beta(x,y,t)$为速度场里该点该时间下速度的角度，$R_{t}$和$R_{s}$分别为轨迹和流线的曲率半径，于是由两者的定义，有
$$
\frac{D\beta}{Ds}=\frac{1}{R_t}\quad\mathrm{and}\quad\frac{\partial\beta}{\partial s}=\frac{1}{R_s}
$$
$$
\frac{D\beta}{Dt}=\frac{D\beta}{Dt} \frac{Ds}{Dt}=\frac{V}{R_{t}}
$$
利用全导数公式可得，
$$
								\frac{D\beta}{Dt}=\frac{ \partial \beta }{ \partial t } +V\frac{ \partial \beta }{ \partial s } =\frac{ \partial \beta }{ \partial t } +\frac{V}{R_{s}}
$$
于是有
$$
				\frac{ \partial \beta }{ \partial t } =V\left( \frac{1}{R_{t}}-\frac{1}{R_{s}} \right)
$$
可见只有当风速场处处方向不随时间变化，即达到稳态，才有两者曲率半径相等，即流线与轨迹重合。

一般来说中纬度天气系统受气高空西风的影响会发生东移，在这种情况下，即使系统内部的等高线不发生变化，即系统受力保持不变，风速也会发生变化。这种情况下两个曲率半径的关系可以从一个随系统以速度$\mathbf{C}$运动的环形等高线出发。

$$
\beta(x,y,t+\delta t)=\beta(x-C_{x}\delta t,y-C_{y}\delta t,t)=\beta(x,y,t)-\mathbf{C}\cdot \nabla \beta
$$
$$
				\frac{ \partial \beta }{ \partial t } =-\mathbf{C}\cdot \nabla \beta
$$
设流线与系统运动方向之间的夹角为$\gamma$，则有
$$
						\frac{ \partial \beta }{ \partial t } =-C\frac{ \partial \beta }{ \partial s } \cos \gamma=-\frac{C}{R_{s}}\cos \gamma
$$
结合前面得出的式子，可得$R_{t}$和$R_{s}$的关系式
$$
R_t=R_s\left(1-\frac{C\cos\gamma}{V}\right)^{-1}
$$
可见当V和C相当时，若用流线的曲率半径代替轨迹的曲率半径，梯度风作为实际风场的近似不比地转风更好。
## 1.6 热成风
地转风在存在水平温度梯度的情况下风速沿垂直方向会存在切变，热成风描述的就是这一现象。

**直观理解：**
根据流体静力平衡方程，
$$
\Phi_{1}-\Phi_{0}\equiv gZ_{T}=R\langle T\rangle \ln\left( \frac{p_{0}}{p_{1}} \right)
$$
其中$Z_{T}$就是两等压面之间的厚度，可见其正比于两等压面之间的平均温度，当存在水平温度梯度时，在压强坐标下不同水平位置处的高度增长速度不同，因此水平气压梯度力沿垂直方向存在切变，进而导致地转风沿垂直方向存在切变，示意图如下

![[热成风示意图.png]]
具体来说，对于压强坐标下的地转风
$$
v_g=\frac{1}{f}\frac{\partial\Phi}{\partial x}\quad\mathrm{and}\quad u_g=-\frac{1}{f}\frac{\partial\Phi}{\partial y}
$$
利用理想气体状态方程，可以将流体静力平衡方程写为
$$
\frac{ \partial \Phi }{ \partial p } =-\frac{1}{\rho}=-\frac{RT}{p}
$$
于是有
$$
p\frac{\partial v_g}{\partial p}\equiv\frac{\partial v_g}{\partial\ln p}=-\frac{R}{f}{\left(\frac{\partial T}{\partial x}\right)}_p
$$
$$
p\frac{\partial u_g}{\partial p}\equiv\frac{\partial u_g}{\partial\ln p}=\frac{R}{f}\left(\frac{\partial T}{\partial y}\right)_p
$$
写成矢量形式为
$$
\frac{ \partial \mathbf{V_{g}} }{ \partial \ln p } =-\frac{R}{f}\mathbf{k}\times \nabla _{p}T
$$
上式称为热成风方程，热成风代表两层之间地转风的差，记为$\mathbf{V}_{T}$，
$$
\mathbf{V}_{T}=V_{g}(p_{1})-V_{g}(p_{0})=-\frac{R}{f}\int_{p_{0}}^{p_{1}}(\mathbf{k}\times \nabla_{p} T)d\ln p
$$
根据在chapter1中所说的，定义两等压面间平均温度
$$
\langle T\rangle = \frac{\int_{p_{0}}^{p_{1}}Td\ln p}{\int_{p_{0}}^{p_{1}}d\ln p}
$$
于是热成风可以写为
$$
\mathbf{V}_{T}=-\frac{R}{f}\mathbf{k}\times\nabla_{p}\langle T\rangle\ln\left( \frac{p_{1}}{p_{0}} \right)
$$
可见热成风方向平行于等温线。

利用上面的等式，我们可以仅通过风速的垂直截面分布就推算出水平温度梯度的分布，也可以在给定某个（等压）层的地转风风速的情况下利用平均温度场估算出其他（等压）层的地转风风速。
### 1.6.1 正压气体和斜压气体
正压气体指的是密度仅依赖于压强的气体，$\rho=\rho(p)$，于是等压面和等密度面重合，对于理想气体，等压面也等温，于是对于正压气体，有
$$
\nabla_{p}T=0
$$
于是热成风方程为
$$
	\frac{ \partial V_{g} }{ \partial \ln p }=0 
$$
表明对于正压气体，地转风独立于高度，这给旋转流体一个很重要的约束，其在大尺度下的运动仅依赖于水平位置与时间。

相应的，斜压气体指密度依赖于压强与温度的气体，$\rho=\rho(p,T)$，其存在速度的水平切变，尽管在大气动力学中，斜压气体更重要，但后面章节中可以看到更简单的正压气体也有它的意义。
## 1.7 垂直运动
在Chapter2里提到过，在天气尺度下垂直方向的速度分量难以直接测量，需要间接得到，这是因为一般其数量级为几cm/s，但气象的测量设备精度大约为1m/s。

为了间接得到垂直运动的速度场，有两种方法，分别为运动学法（基于连续性方程）和绝热法（基于热力学能量方程），出于方便考虑，两种方法均使用p坐标，而z坐标下的速度$w(z)$和p坐标下$\omega(p)$的关系可以由静力平衡方程得到。

$$
\begin{aligned}
	\omega \equiv \frac{Dp}{Dt}&=\frac{ \partial p }{ \partial t } +\mathbf{V}\cdot \nabla p+w\frac{ \partial p }{ \partial z } \\&=\frac{ \partial p }{ \partial t } +\mathbf{V}\cdot \nabla p-\rho gw
\end{aligned}
	
$$
对于天气尺度系统，水平速度的首阶近似为地转风，因此可以写为$\mathbf{V}=\mathbf{V}_{g}+V_{a}$，$V_{a}$即非地转风且远小于$V_{g}$，又因为$\mathbf{V}_{g}=\frac{1}{\rho f}\mathbf{k}\times \nabla p$，所以$\mathbf{V}_{g}\cdot \nabla p=0$，因此有
$$
\omega=\frac{ \partial p }{ \partial t } +\mathbf{V}_{a}\cdot \nabla p-\rho g w
$$
数量级上，
$$
\begin{gathered}\partial p/\partial t\sim10\mathrm{~hPa~d}^{-1}\\\mathbf{V}_a\boldsymbol{\cdot}\boldsymbol{\nabla}p\boldsymbol{\sim}\left(1\mathrm{~m~s}^{-1}\right)\left(1\mathrm{~Pa~km}^{-1}\right)\boldsymbol{\sim}1\mathrm{~hPa~d}^{-1}\\g\rho w\sim100\mathrm{~hPa~d}^{-1}\end{gathered}
$$
因此近似认为
$$
\omega=-\rho g w
$$
### 1.7.1 运动学方法
利用p坐标下的连续性方程
$$
			\frac{ \partial u }{ \partial x } +\frac{ \partial v }{ \partial y } +\frac{ \partial \omega }{ \partial p } =0
$$
将其从参考位置$p_{s}$积分到$p$，
$$
\begin{aligned}\omega\left(p\right)&=\omega\left(p_s\right)-\int_{p_s}^p\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right)_pdp\\&=\omega\left(p_s\right)+(p_s-p)\left(\frac{\partial\left\langle u\right\rangle}{\partial x}+\frac{\partial\left\langle v\right\rangle}{\partial y}\right)_p\end{aligned}
$$
于是z坐标下的垂直速度为
$$
w(z)=\frac{\rho\left(z_s\right)w\left(z_s\right)}{\rho\left(z\right)}-\frac{p_s-p}{\rho\left(z\right)g}\left(\frac{\partial\left\langle u\right\rangle}{\partial x}+\frac{\partial\left\langle v\right\rangle}{\partial y}\right)
$$
估算$w(z)$中的散度项时，采用差分近似
$$
\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\approx\frac{u\left(x_0+d\right)-u\left(x_0-d\right)}{2d}+\frac{v\left(y_0+d\right)-v\left(y_0-d\right)}{2d}
$$
然而对于中纬度天气尺度系统，水平风速近似为地转风，而地转风满足$\frac{ \partial u }{ \partial x }+\frac{ \partial v }{ \partial y }=0$，因此实际的散度值很小，对应于微小的非地转风项，于是测量风速时的唯象误差就会导致估算散度时的很大误差（相对），由于这个原因，不推荐用运动学方法来估算垂直速度。
### 1.7.2 绝热方法
绝热方法对水平速度测量的精度不敏感，利用热力学能量方程。如果吸热速率J相对于其他项较小(绝热方法的名称由来），则有
$$
\omega=S_p^{-1}\left(\frac{\partial T}{\partial t}+u\frac{\partial T}{\partial x}+v\frac{\partial T}{\partial y}\right)
$$
温度的水平梯度可以从地转风中较精确的推知（利用热成风方程），故只需要得到温度的相关数据就容易估算垂直速度。

这个方法的缺点是$\frac{ \partial T }{ \partial t }$这项不好测，需要我们每隔一段时间就测一下局域的温度，另外对于吸热比较强烈的情况（如暴风雨）也不成立，Chapter6将介绍一种更有效的估算方法。
## 1.8 表面压强倾向
所谓表面压强倾向，即指表面压强的时间变化率$\frac{ \partial p_{s} }{ \partial t }$，这是短期天气预报的基础，可以由连续性方程结合$\omega(p)$与$w(z)$的关系得到，取积分上限$p\rightarrow 0$，则相应的$\omega(0)=0$，对连续性方程积分可得
$$
\omega(p_{s})=-\int_{0}^{p_{s}}(\nabla_{p}\cdot \mathbf{V})dp
$$
假设表面水平，因此$w_{s}=0$，
$$
	\frac{ \partial p_{s} }{ \partial t } \approx \omega(p_{s})=-\int_{0}^{p_{s}}\nabla_{p}\cdot \mathbf{V}dp
$$
由静力平衡方程，可得
$$
	\frac{ \partial p_{s} }{ \partial t } =\nabla_{p}\cdot \int_{\infty}^{z_{s}}\rho g\mathbf{V}dz=-\nabla_{p}\cdot \int_{z_{s}}^{\infty}\rho g\mathbf{V}dz
$$
即给定位置的表面压强倾向取决于流入该位置的垂直空气质量。

尽管前面说表面压强倾向对天气预报有用，但其应用受限，一是和运动学方法估算垂直速度类似的原因，$\nabla \cdot \mathbf{V}$难以精确得到，二是因为在垂直方向上，低层大气辐散，高层大气就会辐合，因此在垂直方向上的积分对应的净辐散/辐合就相对较小。然而其也有价值，可以定性地帮助理解大气环流等现象。

假设在大气对流层中部的某处温度突然升高，于是该处上层的等压面升高，水平气压梯度力由中心往外造成气体在该处上层辐散，对应$\nabla_{p}\cdot \mathbf{V}>0$，于是表面压强减小，导致该处下部等压面下降，对应水平气压梯度力由外指向中心造成气体在该处下部辐合，又为了补偿该处上部辐散流出的空气，下部辐合进来的空气垂直向上运动，由此形成了大气环流，示意图如下

![[大气环流示意图.png]]
虚线即对应变化后的等压面。以上讨论提供了一个理解上层气体与表面相互影响动力学过程的一个视角。

另外，这个表面压强倾向方程也可以作为前面提到的一系列方程（运动方程、热力学能量方程、连续性方程等）的边界条件，在p坐标下，可以将该条件用地势来写

借助恒等式
$$
(\frac{\partial p}{\partial  t})_z(\frac{\partial t}{\partial z})_p(\frac{\partial z}{\partial p})_t=-1
$$
可得
$$
					\frac{ \partial p_{s} }{ \partial t } =-\left( \frac{ \partial z_{s} }{ \partial t }  \right)_{p}\frac{ \partial p_{s} }{ \partial z } =\rho_{s} g\left( \frac{ \partial z_{s} }{ \partial t }  \right)_{p}=\frac{p_{s}}{RT_{s}}\frac{ \partial \Phi_{s} }{ \partial t } 
$$
$$
	\frac{ \partial \Phi_{s} }{ \partial t } \approx-\frac{RT_{s}}{p_{s}}\int_{0}^{p_{s}}(\nabla_{p}\cdot V)dp
$$
在p坐标下，上面的边界条件不太好用，因为表面压强$p_{s}$会随位置和时间而变化，在实际应用中，我们会采用$\frac{p}{p_{s}}$代替$p$作为垂直坐标，这会使表面彻底变为等值面，这种方法会在Chpater10中详细讲。