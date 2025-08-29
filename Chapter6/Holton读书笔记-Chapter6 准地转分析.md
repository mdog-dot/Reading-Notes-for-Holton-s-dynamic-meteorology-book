---
zhihu-title: Holton读书笔记-Chapter6 准地转分析
zhihu-topics: 大气科学
zhihu-link: https://zhuanlan.zhihu.com/p/1940906756600489052
toc: "true"
zhihu-updated-at: 2025-08-18 22:46
---
[[Holton读书笔记-Chapter1 基本介绍]](https://zhuanlan.zhihu.com/p/1928544934622924935 "card")
[[Holton读书笔记-Chapter2 基础守恒律]](https://zhuanlan.zhihu.com/p/1929242659588912307 "card")
[[Holton读书笔记-Chapter3 基本方程的初步应用]](https://zhuanlan.zhihu.com/p/1930180220087939862 "card")
[[Holton读书笔记-Chapter4 环量、涡度和位涡]](https://zhuanlan.zhihu.com/p/1934567913379063792 "card")
[[Holton读书笔记-Chapter5 大气振荡（线性摄动理论）]](https://zhuanlan.zhihu.com/p/1934580252622979224 "card")


# 1 Chapter6
气象动力学的基本目标是用控制方程解释观察到的大尺度大气运动的结构，然而这些控制方程即使采用了流体静力平衡的近似仍然复杂，我们要寻求进一步简化的方法。对于热带外的天气尺度运动，其水平速度近似为地转风，这样的运动称为准地转运动，分析起来比热带的扰动或行星尺度的扰动要更简单。准地转运动也是传统短期天气预报感兴趣的基本系统，因此将其作为动力学分析的起点是合理的。

另外，在Chapter5分析的线性波中，性质与这种准地转运动更接近的是Rossby波（静态）而非惯性重力波。Rossby波有较慢的相速度和非0的位涡 ，而惯性重力波移动快速且几乎没有位涡，利用这个性质我们可以从完整的控制方程中过滤出Rossby波的成分，对应于准地转方程（QG方程）。在推导这个方程之前，我们先简要总结一下中纬度天气系统的观测结构和平均环流，然后我们将用PV Thinking（即基于位涡和位温的框架）推导QG方程，而后将会用另一种框架w Thinking（即基于omega方程对垂直运动的诊断）重新研究。

## 1.1 热带外环流的观测结构
>由于没有学过地理，本节第一部分内容我较为不熟悉，以复杂的地理讨论为主，读起来相当痛苦，因此我做的主要是对原书翻译

天气图上的大气环流系统很少像Chapter3那样仅仅为简单的环形旋涡，它们一般在形式上高度的非对称，其最强的风和最大温度梯度都集中在一条狭窄的带上，称为**锋**。另外，这种系统高度斜压（**正压相比于斜压的核心特征在于正压理想气体在等压面上无水平温度梯度，根据热成风方程，进而地转风就没有垂直切变。因此高度斜压表现为强烈的水平温度梯度、显著的垂直风切变，常见于冷暖气团交汇的锋区**），其地势和速度扰动的振幅和相位均随高度变化。这种复杂性部分是源于这种天气尺度系统不是在一个均匀的平均流上叠加的，而是在一个行星尺度的缓变流上叠加的，而这种缓变流本身就是高度斜压的。这种行星尺度的流动还受到地形以及陆海热对流的影响，因此它高度依赖于经度。因此尽管把天气系统扰动看作叠加在一个仅依赖于纬度和高度的平均纬向风场的扰动在理论分析中是一个很有用的一阶近似（在Chapter7中会使用这种方法），更完备的描述需要考虑到纬向的非对称性，即依赖于经度。

纬向平均下的经向截面为行星尺度环流的总体结构提供了有用的信息，下图展示了冬季和夏季经向截面中按经度和时间平均的纬向风场和温度
![[经向截面图.png]]
其中展示的包括对流层和底层平流层，本章主要讨论对流层的风场和温度场结构，平流层会在Chapter12中讨论。上图中可以看到，在北半球由极地到赤道的温度梯度冬天大于夏天，而在南半球冬天和夏天温度分布无显著区别，这主要是因为南半球很大部分面积为海洋，而海洋的比热容大于陆地，因此其热惯性比较大，温度不容易改变。由于平均风场和温度场高精度的满足热成风关系，纬向风速的季节性循环类似于经向温度梯度的循环，在冬天北半球最大的纬向风速是夏天的两倍，而在南半球冬天与夏天的差异很小。此外，无论在哪个季节，纬向风速最大值的中心（称为平均急流轴）位于对流层顶的正下方（即紧贴对流层顶，根据热成风方程对流层顶处积分上限最大，因此风度最快），纬度位于在热成风穿过对流层积分最大值处，在南北半球，冬天对应于约30°纬度处，而在夏天会向极地位移至约45°纬度处。

上面中的纬向平均下达到经向截面不对所有经度下的风场结构具有代表性，这可以通过北半球纬向风场在冬季三个月里的时间平均值在200hPa等压面上的分布看出，如下图
![[纬向风场时间平均.png]]
很明显在一些经度上时间平均纬向流与其按经度平均下的分布有较大偏差。特别的，强烈的纬向风极大值（急流，即对应等值线的中心）发生于北纬30°附近的亚洲东部、北美洲陆地上以及阿拉伯半岛的北部，明显的极小值发生于太平洋东部和大西洋东部（对应同纬度等值线的最外侧）。天气尺度扰动倾向于优先在时间平均纬向风场的极大值处发展，与太平洋西部和大西洋西部的急流联合，沿近似跟随急流轴的风暴轴向下游传播。

![[等高线图.png]]
上图为北半球1月和7月500hPa等压面上的等高线图，从图中也容易看出北半球冬天的气候性急流（等高线密集）与纬向对称的巨大偏差，这与陆地与海洋的分布有密切的联系。最明显的不对称是1月美洲和亚洲陆地东部的槽，根据图6.2可见，在北纬35°、东经140°的强烈急流是该区域半永久性的低压槽导致的，因此显然天气系统所处在的平均流动应该被认为是一个依赖于经度的时间平均流。在北半球夏天，500hPa等压面上的高度由于极地区域的变暖在高纬度地区比热带地区增加得更多，因此极地-赤道的温度差比冬天要小，且高度差和急流要更弱且位于更高纬度。

除了依赖于经度，行星尺度的波动在天与天之间也会变化，由于瞬时的天气尺度扰动的存在。事实上观测表明瞬时行星尺度波动的振幅与时间平均值可比，于是由于急流强度和位置的变化，月度的平均图会平滑掉实际的急流瞬时结构，因此在任何时候急流区域的行星尺度波动比时间平均图上展现的斜压度要更高，这一点可以由下图展现。

![[极锋示意图.png]]
上图展示了纬度-高度截面穿过一个观测到的北美洲上空的急流。（a）图为纬向风速和位温，（b）图为位温和Ertel位涡，2PVU的等位涡线近似标志着对流层顶。如图(a)中所示，急流轴位于一个狭窄倾斜、位温梯度强烈的区域的正上方，该区域称为极锋区，这个区域一般分隔了源于极地的冷空气和热带的暖空气，在这个位温梯度强烈区域的上方产生强烈急流中心是热成风平衡的结果。

上图的等位温线展示了平流层极强的静态稳定性，也体现出等位温面在急流区域附近穿过对流层顶，因此空气可以在对流层和平流层之间无需吸热或放热的移动，但对流层顶强烈的Ertel位涡梯度给沿着等位温线穿过对流层顶的流动提供了强烈的抵抗。注意，锋区的等位涡面下移，于是锋区有强烈的位涡正异常（即位涡异常大），这与急流极地一侧的强相对涡度和冷空气一侧的强静态稳定性（冷空气下沉形成稳定有关。

在流体力学中，速度切变强烈的急流相对于小扰动更可能不稳定在观测中是很常见的，任何发展成不稳定急流的小扰动会在发展过程中倾向于放大、利用急流的能量。对于中纬度天气尺度的系统，基本的不静态稳定性称为斜压不静态稳定性，因为它依赖于经向的温度梯度，或者根据热成风平衡，也可以说是依赖于垂直的风场速度切变。尽管水平温度梯度在锋区附近达到最大，斜压不静态稳定性与锋不静态稳定性不完全相同，因为大多数斜压不静态稳定性的模型描述的仅仅是地转尺度的运动，而强烈的锋区附近的扰动一定是高度非地转的，在Chapter7中会展现，斜压扰动可能会自己作用使已有的温度梯度增强从而产生锋区。

中纬度风暴轴上的扰动以斜压波和更小尺度涡旋的形式呈现。经向风场的统计分析表明，斜压波的波长约为4000km，以约10-15m/s的速度向东传播，其垂直结构为地势高度场随高度向西倾斜，且与流体静力平衡一致，温度场随高度向东倾斜。随着脊与槽之间的上升运动和槽西部的下沉运动，垂直运动在中层对流层达到峰值。这些垂直运动的有组织的模式导致了中纬度的成云和降雨，理解这些环流的动力学来源是这节的主要课题。一个理想的处于发展阶段的斜压系统的关键方面展示在下面的垂直截面图中。
![[斜压波的发展阶段.png]]
可以看到，脊轴和槽轴均随高度向西倾斜，而暖空气轴和冷空气轴往相反方向倾斜（这里的轴表示中心的意思）。在后面会展示，向西倾斜的槽和脊为了让平均流传输能量从而发展波是必要的。而斜压系统发展到了成熟阶段后，其槽在500和1000hPa处几乎同相，因此热对流和能量转换都相当弱了。

这个对斜压波统计性的描述在以热带外气旋为形式的天气扰动中是显而易见的。下图展示了经典的热带外气旋的发展阶段。在快速发展阶段，上层和表面流动之间有协同反应，可以看到强烈的冷对流发生在表面槽的西部，较弱的热对流在东部。这种热对流的模式的直接原因是500hPa的槽滞后于表面槽（槽轴随高度西倾，因此上层槽在下层槽的西部），以至于两层之间的平均地转风（沿着平均等压线吹）穿过等厚线（即等温线），在表面槽西部指向层间距厚处，在表面槽东部指向层间距薄处（根据流体静力平衡，厚度与温度成正比）。
![[热带外气旋发展阶段.png]]
## 1.2 准地转方程的推导
本节的目标是用动力学方程解释观测到的中纬度天气系统的结构。

这部分我觉得Holton在在书中用的方法在近似中非常含糊，在数学上我不太能接受，因此这里我将结合《大气动力学》（二刘）中的推导方式给出准地转方程。

### 1.2.1 小参数方法
二刘中推导准地转方程采用了小参数方法，个人认为这种方法比本书的Chapter5介绍的摄动法在近似上要更严格，因为摄动法统一地将扰动量当作同一阶的小量来处理，这有时候会存在问题（比如对于某些尺度来说，可能扰动量本身的尺度就大于某个非扰动量，或者两个扰动量之间的尺度也差了数量级）。

小参数方法的做法是：首先把方程组无量纲化，然后选择一个合适的小参数，把解各个变量都写成小参数的幂级数代入方程组求各级近似解。在准地转近似中，小参数取为Rossby数$Ro$，对于中纬度天气系统（大尺度），$Ro\sim 0.1$。

取基态为静止态，
$$
\begin{cases}u=u^{\prime},\quad v=v^{\prime},\quad w=w^{\prime},\\p=p_0(z)+p^{\prime},\quad\rho=\rho_0(z)+\rho^{\prime},\quad T=T_0(z)+T^{\prime},\quad\theta^{\prime}=\theta_0(z)+\theta^{\prime},&&\end{cases}
$$
且满足
$$
p^\prime,\rho^\prime,T^\prime,\theta^\prime\ll p_0,\rho_0,T_0,\theta_0
$$
其中$p_{0}(z)$满足
$$
\frac{dp_{0}}{dz}=-\rho_{0}g
$$
于是可以列出运动方程、连续性方程、热力学能量方程
$$
\begin{cases}
\begin{aligned}&\frac{\mathrm{d}u}{\mathrm{d}t}-fv=-\frac{1}{\rho_0}\frac{\partial p^{\prime}}{\partial x},\\&\frac{\mathrm{d}v}{\mathrm{d}t}+fu=-\frac{1}{\rho_0}\frac{\partial p^{\prime}}{\partial y},\\&\frac{\mathrm{d}w}{\mathrm{d}t}=-\frac{1}{\rho_0}\frac{\partial p^{\prime}}{\partial z}-g\frac{\rho^{\prime}}{\rho_0},\\&\frac{\mathrm{d}}{\mathrm{d}t}\left(\frac{\rho^{\prime}}{\rho_0}\right)+\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}+\frac{1}{\rho_0}\frac{\partial\rho_0w}{\partial z}=0,\\&\frac{\theta^{\prime}}{\theta_0}=\frac{1}{\gamma}\frac{p^{\prime}}{p_0}-\frac{\rho^{\prime}}{\rho_0},\\&\frac{\mathrm{d}}{\mathrm{d}t}\left(g\frac{\theta^{\prime}}{\theta_0}\right)+N^2w=0.\end{aligned}
\end{cases}
$$
然后我们对方程组进行无量纲化，即把每个变量的尺度提出来，其中以1为下标的对应无量纲量
$$
\begin{cases}(x,y)=L(x_1,y_1),\quad z=Hz_1,\quad t=\frac{L}{U}t_1,\\(u,v)=U(u_1,v_1),\quad\omega=Ro\frac{UH}{L}\omega_1,\quad f=f_0f_1,\\p^{\prime}=\rho_0f_0ULp_1^{\prime},\quad\rho^{\prime}=\rho_0\mu_0^2Ro\rho_1^{\prime},\quad\theta^{\prime}=\theta_0\mu_0^2Ro\theta_1^{\prime},&\end{cases}
$$
其中我们重点讨论几个变量的尺度
#### 1.2.1.1 垂直速度
垂直速度$w$的尺度我们可以采用下面的方法来估计，将风场分解为地转风和非地转风
$$
\mathbf{V}=\mathbf{V}_{g}+\mathbf{V}_{a}
$$
对于动力学方程
$$
\frac{D\mathbf{V}}{Dt}+f\mathbf{k}\times \mathbf{V}+\frac{1}{\rho}\nabla p'=0
$$
$$
\frac{ \partial \mathbf{V} }{ \partial t } +(\mathbf{V}\cdot \nabla)\mathbf{V}+f\mathbf{k}\times \mathbf{V}_{a}=0
$$   
由于
$$
\frac{ \partial \mathbf{V} }{ \partial t }+(\mathbf{V}\cdot \nabla)\mathbf{V}\sim \frac{U^2}{L}
$$
故有
$$
\mathbf{V}_{a}\sim \frac{U^2}{fL}=RoU
$$
再根据连续性方程有
$$
\frac{ \partial w }{ \partial z } \sim \frac{ \partial u }{ \partial x } +\frac{ \partial v }{ \partial y } =\nabla\cdot \mathbf{V}=\nabla \cdot \mathbf{V}_{a}=Ro \frac{U}{L}
$$
可得垂直速度尺度为
$$
w \sim Ro \frac{UH}{L}
$$
#### 1.2.1.2 扰动量
还是利用动力学方程，当$Ro$很小时，有
$$
\frac{p'}{\rho_{0}}\sim f_{0}UL
$$
故
$$
p'\sim \rho_{0}f_{0}UL
$$
再利用
$$
\frac{\theta^{\prime}}{\theta_0}=\frac{1}{\gamma}\frac{p^{\prime}}{p_0}-\frac{\rho^{\prime}}{\rho_0}
$$

有
$$
\frac{\theta'}{\theta_{0}}\sim \frac{p'}{p_{0}}\sim \frac{\rho'}{\rho_{0}}\sim \frac{\rho_{0}f_{0}UL}{p_{0}}
$$
再利用流体静力学平衡方程，有
$$
p_{0}\sim \rho_{0}gH
$$
所以
$$
\theta'\sim \theta_{0} \frac{f_{0}UL}{gH}
$$
$$
\rho'\sim \rho_{0} \frac{f_{0}UL}{gH}
$$
定义Obukhov参数
$$
\mu_{0}=\frac{L}{L_{R}}=\frac{f_{0}L}{\sqrt{ gH }}
$$
于是上面的数量级又可以用Obukhov参数和Rossby参数表示为
$$
\theta'\sim \theta_{0} Ro\mu_{0}^2
$$
$$
\rho'\sim \rho_{0}Ro\mu_{0}^2
$$
### 1.2.2 推导
将无量纲化下的各量代入原方程组，得到无量纲量的方程组
$$
\left\{\begin{array}{l}\frac{U^2}{L}\frac{du_1}{dt_1}+f_0U(-f_1v_1)=f_0U\left(-\frac{\partial p_1^{\prime}}{\partial x_1}\right),\\\frac{U^2}{L}\frac{dv_1}{dt_1}+f_0U(f_1u_1)=f_0U\left(-\frac{\partial p_1^{\prime}}{\partial y_1}\right),\\Ro\frac{U^2D}{L^2}\frac{dw_1}{dt_1}=\frac{f_0UL}{D}\left(-\frac{1}{\rho_0}\frac{\partial\rho_0p_1^{\prime}}{\partial z_1}\right)+g\mu_0^2Ro(-\rho_1^{\prime}),\\\mu_0^2Ro\frac{U}{L}\frac{d\rho_1^{\prime}}{dt_1}+\frac{U}{L}\left(\frac{\partial u_1}{\partial x_1}+\frac{\partial v_1}{\partial y_1}\right)+Ro\frac{U}{L}\left(\frac{1}{\rho_0}\frac{\partial\rho_0w_1}{\partial z_1}\right)=0,\\\mu_0^2Ro\theta_1^{\prime}=\frac{f_0UL}{c_s^2}p_1^{\prime}-\mu_0^2Ro\rho_1^{\prime},\\\mu_0^2Ro\frac{U}{L}\frac{d\theta_1^{\prime}}{dt_1}+Ro\frac{N^2UD}{gL}w_1=0,\end{array}\right.
$$
其中
$$
\frac{d}{dt_{1}}=\frac{ \partial  }{ \partial t_{1} }+u_{1}\frac{ \partial  }{ \partial x_{1} } +v_{1}\frac{ \partial  }{ \partial y_{1} } +Row_{1}\frac{ \partial  }{ \partial z_{1} }
$$
为了简化表达上式，下面引入一些常用的参数：

大气层结参数$\sigma_{0}\equiv-\frac{ \partial \ln \rho_{0} }{ \partial z }=\frac{N^2}{g}+\frac{g}{c_{s}^2}$，其无量纲量为
$$
\sigma_{1}\equiv-\frac{ \partial \ln \rho_{0} }{ \partial z_{1} } =H\sigma_{0}=\frac{N^2H}{g}+\frac{gH}{c_{s}^2}
$$
定义$\alpha_{0}\equiv \frac{N^2H}{g}$，$\gamma \equiv \frac{c_{p}}{c_{v}}$，则
$$
\sigma_{1}=\alpha_{0}+\frac{1}{\gamma}
$$
于是无量纲量的方程组为(最后一个方程$Ro$不消掉只是为了后面明确近似的级数)
$$
\left\{\begin{array}{l}Ro\frac{du_{1}}{dt_{1}}-f_{1}v_{1}=-\frac{\partial p_{1}^{^{\prime}}}{\partial x_{1}},\\Ro\frac{dv_{1}}{dt_{1}}+f_{1}u_{1}=-\frac{\partial p_{1}^{^{\prime}}}{\partial y_{1}},\\ \frac{H^2}{L^2}Ro^{2}\frac{dw_{1}}{dt_{1}}=-\frac{\partial p_{1}^{^{\prime}}}{\partial z_{1}}+\sigma_{1}p_{1}^{^{\prime}}-\rho_{1}^{^{\prime}},\\\mu_{0}^{2}Ro\frac{d\rho_{1}^{^{\prime}}}{dt_{1}}+\frac{\partial u_{1}}{\partial x_{1}}+\frac{\partial v_{1}}{\partial y_{1}}+Ro\frac{1}{\rho_{0}}\frac{\partial\rho_{0}w_{1}}{\partial z_{1}}=0,\\\theta_{1}^{^{\prime}}=\frac{1}{\gamma}p_{1}^{^{\prime}}-\rho_{1}^{^{\prime}},\\ Ro\left(\frac{d\theta_{1}^{^{\prime}}}{dt_{1}}+\frac{\alpha_{0}}{\mu_{0}^{2}}w_{1}\right)=0.\end{array}\right.
$$
其中第三个方程（垂直运动方程）的左边$\frac{H^2}{L^2}Ro^2\sim 10^{-4}$，因此可以忽略，即静力平衡近似，第四个方程（连续性方程），左边第一项$\mu_{0}^2Ro \sim 10^{-1}Ro \sim10^{-2}$，因此可以忽略，即非弹性近似，而最后一个方程（热力学能量方程）中的第二项系数$\frac{\alpha_{0}}{\mu_{0}^2}\sim 1$。近似过后且利用位温定义（第五个方程）消去密度扰动项$\rho_{1}'$之后得到方程组如下
$$
\left\{\begin{array}{l}Ro\frac{du_{1}}{dt_{1}}-f_{1}v_{1}=-\frac{\partial p_{1}^{\prime}}{\partial x_{1}},\\Ro\frac{dv_{1}}{dt_{1}}+f_{1}u_{1}=-\frac{\partial p_{1}^{\prime}}{\partial y_{1}},\\\frac{\partial p_{1}^{\prime}}{\partial z_{1}}-\alpha_{0}p_{1}^{\prime}-\theta_{1}^{\prime}=0,\\\frac{\partial u_{1}}{\partial x_{1}}+\frac{\partial v_{1}}{\partial y_{1}}+Ro\frac{1}{\rho_{0}}\frac{\partial\rho_{0}w_{1}}{\partial z_{1}}=0,\\Ro\left(\frac{d\theta_{1}^{\prime}}{dt_{1}}+\frac{\alpha_{0}}{\mu_{0}^2}w_{1}\right)=0.\end{array}\right.
$$

接着根据小参数方法，取$Ro$为小参数，将一系列物理量$u_{1},v_{1},w_{1},p_{1}',\theta_{1}'$均展开为$Ro$的幂级数，其中上标该级近似下各物理量对应的值，特别的，对于$w_{1}$，在零级近似中，根据连续性方程可得，$\frac{ \partial u_{1}^{(0)} }{ \partial x_{1} }+\frac{ \partial v_{1}^{(0)} }{ \partial y_{1} }=0$，于是$w_{1}^{(0)}$的零级近似满足$\frac{ \partial \rho_{0}w_{1}^{(0)} }{ \partial z_{1} }=0$，对于平面表面，有边界条件$w|_{z=0}=0$，故零级近似中的$w_{1}^{(0)}=0$，因此$w$的展开从一级近似开始。
$$
\left\{\begin{array}{l}u_{1}=u_{1}^{(0)}+Rou_{1}^{(1)}+\cdots,\\v_{1}=v_{1}^{(0)}+Rov_{1}^{(1)}+\cdots,\\w_{1}=w_{1}^{(1)}+Row_{1}^{(2)}+\cdots,\\p_{1}^{\prime}=p_{1}^{(0)}+Rop_{1}^{(1)}+\cdots,\\\theta_{1}^{\prime}=\theta_{1}^{(0)}+Ro\theta_{1}^{(1)}+\cdots,\end{array}\right.
$$
代入比较方程两端$Ro^0$系数，可得零级近似方程组
$$
\begin{cases}v_1^{(0)}=\frac{\partial p_1^{(0)}}{\partial x_1},&u_1^{(0)}=-\frac{\partial p_1^{(0)}}{\partial y_1},\\\frac{\partial p_1^{(0)}}{\partial z_1}=\theta_1^{(0)},&\frac{\partial u_1^{(0)}}{\partial x_1}+\frac{\partial v_1^{(0)}}{\partial y_1}=0.&\end{cases}
$$
可以将其还原为有量纲的形式，以第一个方程为例
$$
v_{1}^{(0)}=\frac{v^{(0)}}{U}
$$
$$
p_{1}^{(0)}=\frac{p^{(0)}}{\rho_{0}f_{0}UL}
$$
$$
x_{1}=\frac{x}{L}
$$
代入可得有量纲形式
$$
f_{0}v^{(0)}=\frac{1}{\rho_{0}}\frac{ \partial p^{(0)} }{ \partial x }
$$
类似的可得零级近似有量纲形式的方程组
$$
\begin{cases}f_0v^{(0)}=\frac{1}{\rho_0}\frac{\partial p^{(0)}}{\partial x},&f_0u^{(0)}=-\frac{1}{\rho_0}\frac{\partial p^{(0)}}{\partial y},\\\frac{\partial}{\partial z}\left(\frac{p^{(0)}}{\rho_0}\right)=g\frac{\theta^{(0)}}{\theta_0},&\frac{\partial u^{(0)}}{\partial x}+\frac{\partial v^{(0)}}{\partial y}=0.&\end{cases}
$$
可见零级近似对应于地转平衡、静力平衡和水平无散。零级近似反映了大尺度运动天气系统的主要特征，但其表征的是不随时间变化的平衡运动，要反映运动的变化，还需考虑一级近似。

比较方程两端$Ro^1$的系数可得一级近似方程组
$$
\left\{\begin{aligned}&\left(\frac{\partial}{\partial t_{1}}+u_{1}^{(0)}\frac{\partial}{\partial x_{1}}+v_{1}^{(0)}\frac{\partial}{\partial y_{1}}\right)u_{1}^{(0)}-\beta_{1}y_{1}v_{1}^{(0)}-v_{1}^{(1)}=-\frac{\partial p_{1}^{(1)}}{\partial x_{1}},\\&\left(\frac{\partial}{\partial t_{1}}+u_{1}^{(0)}\frac{\partial}{\partial x_{1}}+v_{1}^{(0)}\frac{\partial}{\partial y_{1}}\right)v_{1}^{(0)}+\beta_{1}y_{1}u_{1}^{(0)}+u_{1}^{(1)}=-\frac{\partial p_{1}^{(1)}}{\partial y_{1}},\\ &\frac{\partial p_{1}^{(1)}}{\partial z_{1}}=\theta_{1}^{(1)}+\alpha_{0}Ro^{-1}p_{1}^{(0)},\\&\frac{\partial u_{1}^{(1)}}{\partial x_{1}}+\frac{\partial v_{1}^{(1)}}{\partial y_{1}}+\frac{1}{\rho_{0}}\frac{\partial\rho_{0}w_{1}^{(1)}}{\partial z_{1}}=0,\\&\left(\frac{\partial}{\partial t_{1}}+u_{1}^{(0)}\frac{\partial}{\partial x_{1}}+v_{1}^{(0)}\frac{\partial}{\partial y_{1}}\right)\theta_{1}^{(0)}+\frac{\alpha_{0}}{\mu_{0}^2}w_{1}^{(1)}=0.\end{aligned}\right.
$$
用上面类似的方法可以还原为有量纲的形式
$$
\begin{cases}\left(\frac{\partial}{\partial t}+u^{(0)}\frac{\partial}{\partial x}+v^{(0)}\frac{\partial}{\partial y}\right)u^{(0)}-\beta_0yv^{(0)}-f_0v^{(1)}=-\frac{1}{\rho_0}\frac{\partial p^{(1)}}{\partial x},\\\left(\frac{\partial}{\partial t}+u^{(0)}\frac{\partial}{\partial x}+v^{(0)}\frac{\partial}{\partial y}\right)v^{(0)}+\beta_0yu^{(0)}+f_0u^{(1)}=-\frac{1}{\rho_0}\frac{\partial p^{(1)}}{\partial y},\\\frac{\partial}{\partial z}\left(\frac{p^{(1)}}{\rho_0}\right)=g\frac{\theta^{(1)}}{\theta_0}+\frac{N^2}{g}\frac{p^{(0)}}{\rho_0},\\\frac{\partial u^{(1)}}{\partial x}+\frac{\partial v^{(1)}}{\partial y}+\frac{1}{\rho_0}\frac{\partial\rho_0w^{(1)}}{\partial z}=0,\\\left(\frac{\partial}{\partial t}+u^{(0)}\frac{\partial}{\partial x}+v^{(0)}\frac{\partial}{\partial y}\right)\left(g\frac{\theta^{(0)}}{\theta_0}\right)+N^2w^{(1)}=0.&\end{cases}
$$
一级近似方程组不仅展示了物理量随时间的变化，还建立了零级近似与一级近似的联系，在平流项中、含$\beta$的项中水平运动均用地转风替代，因此一级近似方程组对应的就是所谓准地转运动。直接求解一级近似方程组并不方便，其中包含了$u^{(1)},v^{(1)},p^{(1)}$等一级小量，而我们讨论的目标则是利用一级近似方程组给出地转风随时间和空间的分布，因此要想办法把一级小量消去，仅保留零级量（即地转风）。

用$\frac{ \partial  }{ \partial x }$作用于第二个方程，减去$\frac{ \partial  }{ \partial y }$作用于第一个方程，利用零级地转风无散$\frac{ \partial u^{(0)} }{ \partial x }+\frac{ \partial v^{(0)} }{ \partial y }=0$，再代入连续性方程，可得准地转涡度方程
$$
		\left(\frac{\partial}{\partial t}+u^{(0)}\frac{\partial}{\partial x}+v^{(0)}\frac{\partial}{\partial y}\right)\zeta^{(0)}+\beta_0v^{(0)}=f_{0}\left( \frac{ \partial u^{(1)} }{ \partial x }+\frac{ \partial v^{(1)} }{ \partial y }   \right)=f_0\frac{1}{\rho_0}\frac{\partial \rho_0w^{(1)}}{\partial z}
$$

其中
$$
\zeta^{(0)}\equiv\frac{\partial v^{(0)}}{\partial x}-\frac{\partial u^{(0)}}{\partial y}=\frac{1}{f_{0}}\nabla_{h}^2p^{(0)}
$$
于是方程组可以改写为关于零级量的方程组
$$
\begin{cases}\left(\frac{\partial}{\partial t}+u^{(0)}\frac{\partial}{\partial x}+v^{(0)}\frac{\partial}{\partial y}\right)\zeta^{(0)}+\beta_0v^{(0)}=\left(\frac{\partial}{\partial t}+u^{(0)}\frac{\partial}{\partial x}+v^{(0)}\frac{\partial}{\partial y}\right)(\zeta^{(0)}+f)=f_0\frac{1}{\rho_0}\frac{\partial\rho_0w^{(1)}}{\partial z},\\\left(\frac{\partial}{\partial t}+u^{(0)}\frac{\partial}{\partial x}+v^{(0)}\frac{\partial}{\partial y}\right)\left(g\frac{\theta^{(0)}}{\theta_0}\right)+N^2w^{(1)}=0,&\end{cases}
$$
两个零级量由$w^{(1)}$联系在一起。上面的方程称为准地转模式，包含一个涡度方程和一个热力学能量方程，在涡度方程中所有风场均用地转风代替，除了水平散度用了非地转风。我们也可以将小参数方法应用于涡度方程而得到类似的结果。
### 1.2.3 准地转模式和准地转位涡守恒
准地转模式的方程组通过$w^{(1)}$联系起来，消去$w^{(1)}$可得完全关于零级量的偏微分方程，事实上，我们可以将所有零级量用一个函数来表示，引入准地转流函数
$$
\psi=\frac{p'}{f_{0}\rho_{0}}
$$
于是前面零级近似的方程组可以写为
$$
\begin{cases}u^{(0)}=-\frac{\partial\psi}{\partial y},\quad v^{(0)}=\frac{\partial\psi}{\partial x},\quad\frac{\partial u^{(0)}}{\partial x}+\frac{\partial v^{(0)}}{\partial y}=0,\\\frac{\partial}{\partial z}(f_0\psi)=g\frac{\theta^{\prime}}{\theta_0}.&\end{cases}
$$
由于$p',\theta'$只用到零级近似，故均用撇代替上标0，利用上面的式子可以将全部零级量均用$\psi$表示，消去$w^{(1)}$以后可得关于$\psi$的偏微分方程
$$
w=-\frac{f_0}{N^2}\left(\frac{\partial}{\partial t}+u^{(0)}\frac{\partial}{\partial x}+v^{(0)}\frac{\partial}{\partial y}\right)\left(\frac{\partial\psi}{\partial z}\right),
$$
将其代入零级量方程组第一个方程的右边项可得
$$
\begin{aligned}f_0\frac{1}{\rho_0}\frac{\partial\rho_0w}{\partial z}&=-\frac{1}{\rho_{0}}\frac{\partial}{\partial z}\left[\left(\frac{\partial}{\partial t}+u^{(0)}\frac{\partial}{\partial x}+v^{(0)}\frac{\partial}{\partial y}\right)\left(\frac{f_{0}^{2}}{N^{2}}\rho_{0}\frac{\partial\psi}{\partial z}\right)\right]\\&\mathrm{=}-\frac{1}{\rho_{0}}\left(\frac{\partial}{\partial t}+u^{(0)}\frac{\partial}{\partial x}+v^{(0)}\frac{\partial}{\partial y}\right)\frac{\partial}{\partial z}\left(\frac{f_{0}^{2}}{N^{2}}\rho_{0}\frac{\partial\psi}{\partial z}\right)\\&-\frac{f_0^2}{N^2}\left[\frac{\partial u^{(0)}}{\partial z}\frac{\partial}{\partial x}\left(\frac{\partial\psi}{\partial z}\right)+\frac{\partial v^{(0)}}{\partial z}\frac{\partial}{\partial y}\left(\frac{\partial\psi}{\partial z}\right)\right]\\&\mathrm{=}-\frac{1}{\rho_{0}}\left(\frac{\partial}{\partial t}+u^{(0)}\frac{\partial}{\partial x}+v^{(0)}\frac{\partial}{\partial y}\right)\frac{\partial}{\partial z}\left(\frac{f_{0}^{2}}{N^{2}}\rho_{0}\frac{\partial\psi}{\partial z}\right)\\&-\frac{f_0^2}{N^2}{\left[-\frac{\partial^2\psi}{\partial y\partial z}\frac{\partial^2\psi}{\partial x\partial z}+\frac{\partial^2\psi}{\partial x\partial z}\frac{\partial^2\psi}{\partial y\partial z}\right]}\\&\mathrm{=}-\frac{1}{\rho_{0}}\left(\frac{\partial}{\partial t}+u^{(0)}\frac{\partial}{\partial x}+v^{(0)}\frac{\partial}{\partial y}\right)\frac{\partial}{\partial z}\left(\frac{f_{0}^{2}}{N^{2}}\rho_{0}\frac{\partial\psi}{\partial z}\right).\end{aligned}
$$
于是可得
$$
\left.\left(\frac{\partial}{\partial t}+u^{(0)}\right.\frac{\partial}{\partial x}+v^{(0)}\frac{\partial}{\partial y}\right)q=\frac{D}{Dt_{g}}q=0,
$$
其中
$$
q\equiv f+\zeta^{(0)}+\frac{1}{\rho_0}\frac{\partial}{\partial z}\left(\frac{f_0^2}{N^2}\rho_0\frac{\partial\psi}{\partial z}\right)=f+\nabla_0^2\psi+\frac{1}{\rho_0}\frac{\partial}{\partial z}\left(\frac{f_0^2}{N^2}\rho_0\frac{\partial\psi}{\partial z}\right)
$$
称为准地转位涡，满足准地转位涡守恒定律，又称为Charney-Obukhov方程，准地转位涡守恒方程仅包含$\psi$一个未知变量，因此其有广泛的应用，可以在仅知道压强场（或流函数场）在初始时刻的状态的情况下，对系统的演变做出预测。事实上，准地转位涡守恒是大气大尺度运动特征的综合体现。
### 1.2.4 与Holton书中公式的对应
Holton书中方程表示方式均较为简洁，现在把上面推导的方程与Holton给出的方程做个对应。

首先令
$$
\frac{D}{Dt_{g}}\equiv \frac{ \partial  }{ \partial t } +u_{g}\frac{ \partial  }{ \partial x } +v_{g}\frac{ \partial  }{ \partial y }
$$
一级近似中的热力学能量方程称为准地转热力学能量方程，可以写为
$$
\frac{D\theta}{Dt_{g}}=-w \frac{d\theta_{0}}{dz}
$$
如果认为$\rho_{0}$为常数，即采用Boussinesq近似，则准地转位涡(QG PV)为
$$
q=\zeta_{g}+f+f_{0}\frac{ \partial  }{ \partial z } \left( \frac{d\theta_{0}}{dz}^{-1}\theta' \right)
$$
若做$f$-平面近似，即认为$f$近似不变，则所有$f_{0}$均可以写为$f$
$$
q=\zeta_{g}+f+f\frac{ \partial  }{ \partial z } \left( \frac{d\theta_{0}}{dz}^{-1}\theta' \right)
$$
$$
\frac{D}{Dt_{g}}q=0
$$
准地转模式中的涡度方程称为准地转涡度方程，写为
$$
\frac{D}{Dt_{g}}(\zeta_{g}+f)=f\frac{ \partial w }{ \partial z }
$$
一级近似中的连续性方程称为准地转连续性方程，其中$u,v$的一级近似对应于非地转量，用下标a表示
$$
\frac{ \partial u_{a} }{ \partial x } +\frac{ \partial v_{a} }{ \partial y } +\frac{ \partial w }{ \partial z } =0
$$
对于一级近似中的动力学方程，由于$p^{(1)}$可以取为任意的场（都对应一个特定的$\mathbf{V}^{(1)}$），我们不妨假设无非地转的压强$p^{(1)}$，即压强的扰动量$p'$就是其零级近似$p^{(0)}$。可以将该动力学方程写为矢量形式
$$
\frac{D\mathbf{V}_{g}}{Dt_{g}}=-f\mathbf{k}\times \mathbf{V}_{a}
$$

### 1.2.5 准地转近似总结
1. QG方程是在较小Rossby数和强层结的参考大气（即密度、温度、压强等物理量只与垂直坐标z有关）下的近似
2. QG PV是线性化的Ertel位涡，由绝对涡度的垂直分量以及静态稳定性的扰动组成
## 1.3 位涡思想
>在天气系统的方程中，有两种量，在$\frac{ \partial  }{ \partial t }$下的量是预测量（也可以进一步延伸为不容易直接测得的导数关系），而其他量称为诊断量。一般来说要预测系统的演化，往往可以通过给出诊断量的初始状态，再结合边界条件给出预测。下面要介绍的位涡思想就是将位涡$q$或压强场$p'$作为诊断量，后面的垂直运动思想则消去了$\frac{ \partial  }{ \partial t }$项，而将不容易测得的垂直速度$w$作为诊断量。

在Chapter4中提到过，位涡守恒用一个方程囊括了所有基础的守恒律（从上面的推导过程中可以看出，推导时同时使用了动量守恒方程、质量守恒连续性方程、热力学能量方程），其为天气系统的动力学提供了一个有用且简洁的解释基础。

首先我们考察一下QG PV的尺度，根据前面的推导，以及$f$-平面近似以及Boussinesq近似，QG PV可以写为
$$
		q-f=\frac{1}{\rho_{0}f}\nabla_{h}^2p'+\frac{1}{\rho_{0}}\frac{ \partial  }{ \partial z } \left( \frac{f}{N^2}\frac{ \partial p' }{ \partial z }  \right)
$$
利用前面的量纲分析方法可得(Holton书中习惯用上标$\hat{\#}$表示无量纲量，意义与前面的下标1相同)
$$
	q-f \sim \frac{U}{L}\left( \hat{\nabla}_{h}^2\hat{p}+\frac{f^2L^2}{N^2H^2}\frac{ \partial ^2\hat{p} }{ \partial \hat{z}^2 }  \right)=\frac{U}{L}\left( \hat{\nabla}_{h}^2\hat{p}+\frac{1}{B^2}\frac{ \partial ^2\hat{p} }{ \partial \hat{z}^2 }  \right)
$$
其中$B$为Burger数，是衡量位涡中层结（N）相对于旋转（涡度）的重要程度的基本参数，它有很多表达公式，
$$
B=\frac{NH}{fL}=\frac{L_{R}}{L}=\frac{H}{H_{R}}=\frac{Ro}{Fr}
$$
第二个等号表明Burger数为Rossby半径$L_{R}=\frac{NH}{f}$除以水平尺度L，较大的水平扰动尺度会导致较小的Burger数，此时主导运动的是旋转项（涡度），第三个等号表明Burger数为高度尺度H除以Rossby深度$H_{R}=\frac{fL}{N}$，Rossby深度衡量的是扰动的垂直影响。第四个等号表明Burger数为Rossby数$Ro$除以Froude数$Fr=\frac{U}{NH}$，两者均为无量纲数，其中Froude数是用来衡量扰动的垂直切变（对于地转风即水平温度梯度）相对于浮力频率$N$的重要性的。

对于准地转方程，由于$Ro$较小，而由于高度层结，故浮力频率$N$较大，即$Fr$较小，因此Burger数接近于1，这表明在QG PV中，旋转项（涡度项）与层结项同等重要，因此$L \sim L_{R},H \sim H_{R}$，这为我们讨论天气特征的动力学提供了指导。
### 1.3.1 位涡反演、流和分段位涡反演
位涡反演是一种运动学诊断技术，可以从位涡场推导出其他场的状态。

假设B为1，则位涡$q$的无量纲量$\hat{q}$可以写为
$$
\hat{q}-Ro^{-1}=\hat{\nabla}^2\hat{p}
$$
$Ro$虽然很大，但其仅在各处给PV加了一固定值，无动力学意义，因此可以说我们感兴趣的QG PV的特征由压强场的拉普拉斯所决定，低压区域对应较大PV，高压区域对应较小PV。所谓QG PV的反演就是从QG PV场中推导出压强场，这可以通过解
$$
\hat{p}=\hat{\nabla}^{-2}\hat{q}
$$
所给出。

回到有量纲形式的QG PV，两端同乘$f$可得
$$
f(q-f)=\frac{1}{\rho_{0}}\nabla_{h}^2p'+\frac{1}{\rho_{0}}\frac{ \partial  }{ \partial z } \frac{f^2}{N^2}\frac{ \partial p' }{ \partial z }=\mathbf{L}p'
$$
其中
$$
L=\frac{1}{\rho_{0}}\nabla_{h}^2+\frac{1}{\rho_{0}}\frac{ \partial  }{ \partial z } \frac{f^2}{N^2}\frac{ \partial  }{ \partial z }
$$
接近于三维拉普拉斯算子，但其在z轴相对于x，y轴有$\frac{f^2}{N^2}$的拉伸，在无量纲形式中忽略了这一拉伸。上面之所以同乘了$f$是因为这样可以使得位涡思想是半球无关的，即对于气旋来说无论在南北半球$f(q-f)$均为正。$\mathbf{L}$作为近似拉普拉斯算子，与前面讨论的具有相同的性质，在低压槽（p极小值）处$f(q-f)$取到极大值。类似的，我们可以定义$\mathbf{L}$的反演算符$\mathbf{L}^{-1}$，从而有
$$
p'=\mathbf{L}^{-1}(fq)
$$
其中省略了为0的项$\mathbf{L}^{-1}f^2$。在按照上面的方式推导出压强场以后，相应的地转风和位温扰动就可以利用前面推导的公式得出，这样来看，风场和位温场就是由QG PV诱导出来的。

在反演位涡时，要求解方程需要边界条件，假设水平边界条件为周期性的，表面上水平边界条件可能有Dirichlet条件（即给定表面上的压强值）、Neumann条件（即给定表面上压强沿垂直方向的导数值）、Robin条件（即给定压强和其垂直导数的线性组合值）。在没有边界条件时，其解称为自由空间解。

下面一段关于边界条件的讨论在Holton中仅给出了概括性的文字讨论，我根据deepseek的回答查阅了Hoskin的原始论文做一些简单的论述，其中我仍然觉得很多地方非常奇怪，欢迎大佬指正。

#### 1.3.1.1 边界条件讨论
表面流体静力平衡方程给出了Neumann条件
$$
\left.\frac{ \partial p' }{ \partial z } \right|_{z=0} =\frac{\rho_{0}\theta'_{surf}}{\theta_{0}}g
$$
但是由于表面位温$\theta'_{surf}$为一函数，因此边界条件为非齐次边界条件，这给我们用分离变量法求解PV反演带来了困难，一种常用的处理方法如下：

我们的目标是使方程在$z=0$处满足齐次边界条件（deepseek称其为形式边界），而在真实物理边界$z=0^+$处对应真实的边界条件。首先为了保证$z=0$处满足齐次边界条件，可以对$p'$场做偶延拓，即定义其在$z<0$处的场与$z>0$处对称，于是保证了
$$
\left.\frac{ \partial p' }{ \partial z } \right|_{z=0}=0
$$
这里其实有个问题在于如果$z=0^+$处的导数不为0，那么$z=0$处的导数实际上没有定义，deepseek说这里是广义导数，关于这里数学上的严谨性我觉得仍然存疑。

接下来，为了保证$z=0^+$处仍为真实边界，我们需要给原始方程添加一项，原始方程为
$$
\begin{cases}
q-f=\frac{1}{f}\mathbf{L}p' & z > 0 \\
\left. \frac{\partial p'}{\partial z} \right|_{z=0} = \frac{\rho_{0}\theta'_{surf}}{\theta_{0}}g & \text{(边界条件)}
\end{cases}
$$
给QG PV添加一项δ函数得到新的方程如下
$$
\begin{cases}
	q-f=\frac{1}{f}\mathbf{L}p'+\alpha \theta'_{surf}\delta(z) & z > 0 \\
\left. \frac{\partial p'}{\partial z} \right|_{z=0} = 0 & \text{(边界条件)}
\end{cases}
$$
现在可以利用偶函数性质以及$z=0^+$处边界条件为真实条件来确定待定系数$\alpha$。将新方程从$-\varepsilon$积分到$+\varepsilon$，并令$\varepsilon \to0$，可得
$$
\begin{aligned}
					\frac{f}{N^2\rho_{0}}\int_{-\varepsilon}^{+\varepsilon}\frac{ \partial ^2p' }{ \partial z^2 } dz&=\frac{f}{N^2\rho_{0}}\left.\frac{ \partial p' }{ \partial z } \right|_{-\varepsilon}^{+\varepsilon}\\&=\frac{2f}{N^2\rho_{0}}\left.\frac{ \partial p' }{ \partial z } \right|_{z=0^+}\\&=\frac{2f}{N^2\theta_{0}}\theta'_{surf}=-\alpha \theta'_{surf}
\end{aligned}
$$
可得
$$
\alpha=-\frac{2f}{N^2\theta_{0}}
$$
此时就保证了$z=0^+$边界处对应于真实边界条件。

因此可以看到，表面处的位温分布等价于给表面处QG PV添加一项正比于表面位温分布的$\delta$函数（这在Holton中被概括为峰spike），而这种将表面位温转化为表面QG PV的方法称为$PV-\theta$ perspective，对于理解对流层中的PV异常（比如对流层顶的波动或者地表气旋等）很有用。

### 1.3.2 经典案例
为了阐明PV反演概念的应用，我们讨论两个经典的扰动。

第一个是短波槽，一般发生于天气图中的对流层上部，是独立的低压槽。这种扰动是由局部的气旋性QG PV异常导致的。下图展示了最简单的一种情况，即在纬向风中嵌入一个PV异常扰动。
![[短波槽示意图.png]]
其中，图a展示的是PV异常情况下的等PV线，箭头指的是异常PV所反演得到的风场，可见PV异常可以影响到很遥远的风场。图b展示的是等压线和完整的风速场（之所以等压线整体是东西向的，是因为背景风场取为了地转风），在PV异常处可以看到明显的低压槽。图c是异常的等压线（即减去了原始的压强），可见局域PV异常的非局域影响，使得表面出现了低压区域。图d是经向风速场，以及位温异常线，可见局域PV异常对应于气旋性涡度的极大值，同时附近呈现暖空气在上冷空气在下，静态稳定性达到局部极大值。

在地球表面，暖空气对应于低压（因为密度小，上升，对应气旋性环流），冷空气对应于高压，这股环流随高度往上逐渐衰减。如果把对流层顶当作刚性边界，则冷空气对应低压（气旋性环流），暖空气对应高压，这股环流随高度向下衰减。

诱导流的概念在考虑大气的离散部分和它们的相互作用动力学时很有效，由于QG PV算子$\mathbf{L}$是线性的，大气可以拆分成一个一个部分，每个部分贡献的场的叠加就是完整的场。这称为分段位涡反演，写为
$$
p=\sum_{i=1}^Np_{i}=\sum_{i=1}^N\mathbf{L}^{-1}(fq_{i})
$$
又由于地转关系和流体静力平衡方程都是线性的，也可以先根据$p_{i}$推导出分段风速场再叠加出完整风速场。为了展示分段PV反演的概念，我们引入第二个经典的扰动：向西急流中的局域直线极大风，被称为“jet streak”。这一直线极大风是由两个QG PV异常（下图中的椭圆形）导致的，其中气旋性的QG PV异常比反气旋性的要更靠极地。
![[急流中的极大风.png]]
上图可以看出，这一上一下一气旋一反气旋的QG PV异常导致两者中间处的诱导流极强，因为两个涡旋都对中间的流动起到增强作用。从c，d图可以看到，这里的整体场就可以用前面的分段位涡反演得到，即将c、d叠加起来得到的就是a图。

从c、d中也可以看到，自身的QG PV异常对自身异常（椭圆）的平流作用很小，主要是由另一个QG PV异常起到平流作用，这种平流作用会使得异常以一个比就留风速小的速度对称地向下游移动。
### 1.3.3 位涡守恒和准地转压强倾向方程
前面讨论的是准地转位涡的反演，利用反演的功能，我们还可以通过位涡思想进行预测。

位涡守恒可以写为
$$
\frac{ \partial q }{ \partial t } =-\mathbf{V}_{g}\cdot \nabla_{h}q
$$
如果给定了初始QG PV的分布，那么根据位涡反演我们就可以得到压强场分布，从而得到地转风分布$\mathbf{V}_{g}$，代入可得到$q$随时间的变化，而后再可利用位涡反演推得新的地转风，无限循环下去。

以上过程的本质其实就是个关于压强场$p'$的偏微分方程，因此边界条件也需要考虑。在地球表面，假设其可以近似为刚性水平面，则$\left.w\right|_{z=0}=0$，因此前面的准地转热力学能量方程可以写为
$$
\frac{ \partial \theta }{ \partial t } =-\mathbf{V}_{g}\cdot \nabla_{h}\theta
$$
这方程形式与前面的位涡守恒方程一致，因此推导QG PV的同时，也可以推导位温，这给前面的偏微分方程提供了边界条件。除了地球表面以外，上部边界条件也需要考虑，可以假设为另一个刚性水平边界（粗略地代替对流层顶），然后也可以用上式的热力学能量方程表达。

然而，一般我们不会直接测量QG PV，而是采取更直观的方式来预测压强场，正如前面所说，以上无限循环下去的过程其实本质就是求解一个关于压强场$p'$的偏微分方程，我们可以通过前面$q$的定义式得到这个方程
$$
\mathbf{L}\frac{ \partial p' }{ \partial t } =-\mathbf{V}_{g}\cdot \nabla_{h}(\mathbf{L}p')
$$
将右侧代入$\mathbf{L}$的定义可得
$$
\begin{aligned}
		\mathbf{L}\frac{ \partial p' }{ \partial t } =-\mathbf{V}_{g}\cdot \nabla_{h}(f\zeta_{g})-f^2\mathbf{V}_{g}\cdot \nabla_{h}\left( \frac{ \partial  }{ \partial z } \left( \frac{d\theta_{0}}{dz}^{-1}\theta' \right) \right)
\end{aligned}
$$
利用热成风方程（z坐标下的）
$$
\frac{ \partial \mathbf{V}_{g} }{ \partial z } =\frac{g}{f\theta_{0}}\mathbf{k}\times \nabla_{h}\theta'
$$
由于
$$
(\mathbf{k}\times \nabla_{h}\theta')\cdot \nabla_{h}\theta'=0
$$
可得
$$
	\mathbf{L}\frac{ \partial p' }{ \partial t } =-\mathbf{V}_{g}\cdot \nabla_{h}(f\zeta_{g})-f^2\frac{ \partial  }{ \partial z } \left[ \frac{d\theta_{0}}{dz}^{-1}\left(\mathbf{V}_{g}\cdot \nabla_{h}\theta' \right) \right]
$$
这就是准地转压强倾向方程。其边界条件可以利用流体静力平衡方程推知
$$
			\frac{ \partial  }{ \partial z } \frac{ \partial p' }{ \partial t } =\frac{\rho_{0}g}{\theta_{0}}(-\mathbf{V}_{g}\cdot \nabla_{h}\theta')
$$
可见表面暖空气的平流（即从暖移到冷）会带来压强倾向从表面向上增大。

由于$\mathbf{L}$是拉伸的拉普拉斯算子，因此也具有类似的负定性，可以定性地理解为，当$\mathbf{L}\psi=S>0$时，对应于$\psi$的局域最低点，而当$\mathbf{L}\psi=S<0$时，对应于$\psi$的局域最高点，因此近似可以认为算符$\mathbf{L}$具有反号的功能。

根据压强倾向方程以及$\mathbf{L}$的反号功能，压强场$p'$可能会在下面两种情况中随时间衰减，分别对应方程的两项：
1. 气旋性的地转涡度从异常中心随地转风平流（对应$\mathbf{V}_{g}\cdot \nabla_{h}(f\zeta_{g})<0$）
2. 位温随地转风的平流随高度上升而减小
>不知为何Holton书中所有说的平流都指的是$-\mathbf{V}_{g}\cdot \nabla_{h}$对应的平流，但我认为这并不是真正的地转运动造成的平流，因此这里我关于平流的描述均与书中恰好相反

第一种情况的例子：发生于低压槽的下游，由于气旋性涡度在低压槽中心处达到极大值。第二种情况的例子：在低对流层处，冷空气平流随高度上升而减小。注意如果大气层结非常稳定（即$\frac{d\theta_{0}}{dz}$很大），则主导压强倾向的就是涡度平流项。
![[压强倾向方程应用实例.png]]还是用前面的均匀西风（背景风为地转风$\mathbf{V}_{g}$）中嵌入的QG PV异常为例，上图为此例子下的压强倾向（$\frac{ \partial p }{ \partial t }$）等值线图。根据压强倾向方程，在异常上游（西侧），气旋性地转涡度随地转风的平流大于0，因此压强随时间增加，在异常下游（东侧），则相反，压强随时间减少，因此如上图a所示，在异常西侧为+，东侧为- 。而在垂直方向，如图b所示，在QG PV异常处东侧，压强随时间下降，而且上层的这一压强随时间下降的特征一直影响到了表面，形成很深的一层压强下降区。根据上面的边界条件，PV异常东侧的暖空气平流导致了表面处$\frac{ \partial p }{ \partial t }$随z而增加（位温的地转平流小于0）。由于前面提到的，在表面暖空气对应气旋性QG PV，而在对流层顶的边界处暖空气对应反气旋性QG PV，因此表面的温度平流（暖空气平流）会增强对流层中PV异常向表面传播的正PV的影响，而对流层顶的温度平流（暖空气平流）会减弱PV异常向对流层顶传播的正PV的影响。

上图中的c，d展示的是涡度平流项和位温平流随高度变化项单独对压强倾向的贡献，可以看到由于西风随高度增加，涡度平流造成的压强倾向的极值在PV异常点的上方一些（即极值处还不是平流最大的地方，可能是因为上方地转风更快），而温度平流造成的压强倾向在表面与对流层之间会发生反号，这部分是由于边界条件的影响。
## 1.4 垂直运动思想（w-Thinking）
在PV Thinking中垂直运动不起作用，因为位涡的平流完全由水平的地转风所决定。而在垂直运动思想的框架下，准地转动力学由前面推出的准地转模式方程组（热力学能量方程&涡度方程）开始讨论
$$
\begin{cases}
		 &\frac{ \partial \theta' }{ \partial t  }=-\mathbf{V} _{g}\cdot \nabla_{h}\theta'-w \frac{d\theta_{0}}{dz} \\&\frac{ \partial  }{ \partial t } (\zeta_{g}+f)=-\mathbf{V}_{g}\cdot \nabla_{h}(\zeta_{g}+f)+f\frac{ \partial w }{ \partial z } 
\end{cases}
$$
可以看到两个零级量之间通过垂直速度$w$建立了联系，在前面推导位涡守恒的时候，我们是消去了$w$，当然现在我们也可以选择消去别的量，比如可以将含$\frac{ \partial  }{ \partial t }$的预测量消去。利用
$$
\frac{ \partial \theta' }{ \partial t } =\frac{\theta_{0}}{\rho_{0}g}\frac{ \partial ^2p' }{ \partial z\partial t } ,\quad \frac{ \partial \zeta_{g} }{ \partial t } =\frac{1}{\rho_{0}f}\nabla_{h}^2\frac{ \partial p' }{ \partial t }
$$
对第一个方程作用$\frac{g}{f\theta_{0}}\nabla_{h}^2$，对第二个方程作用$\frac{ \partial  }{ \partial z }$，取为$f$-平面近似，两式相减可得
$$
\nabla_{h}^2w+\frac{f^2}{N^2}\frac{ \partial ^2w }{ \partial z^2 } =-\frac{g}{\theta_{0}N^2}\nabla_{h}^2(\mathbf{V}_{g}\cdot \nabla_{h}\theta')+\frac{f}{N^2}\frac{ \partial  }{ \partial z } (\mathbf{V}_{g}\cdot \nabla_{h}\zeta_{g})
$$
这就是传统形式的准地转垂直运动方程，又称为$\omega$-方程，方程左边作用于$w$的算符是准拉普拉斯算符，和前面的$\mathbf{L}$只差了一个系数
$$
\tilde{\mathbf{L}}=\nabla_{h}^2+\frac{f^2}{N^2}\frac{ \partial ^2 }{ \partial z^2 }
$$
可以看到，当$f$较小而静态稳定性较大时，对于给定的诊断量，$w$的反应主要发生在水平方向上。

从上面的准地转垂直运动方程中可以看到，利用$\tilde{\mathbf{L}}$的负定性（即反号），（和前面类似，这里同样与Holton原书中讨论相反）可知向下的垂直运动（$w<0$）总是与位温局部极大值处的暖空气平流和地转涡度平流随高度的增大相联系。

传统形式的$\omega$-方程有个问题在于方程右边两项其实有一部分可以相消，去除这部分需要采用爱因斯坦求和约定。（Holton书里做了更详细的介绍，我就默认这是很基础的数学知识了）
$$
\begin{aligned}
\nabla_{h}^2(\mathbf{V}_{g}\cdot \nabla_{h}\theta')&=\partial_{i}\partial_{i}(v_{gj}\partial_{j}\theta')\\&=\partial_{i}(v_{gj}\partial_{i}\partial_{j}\theta'+\partial_{i}v_{gj}\partial_{j}\theta')\\&=v_{gj}\partial_{j}\partial_{i}^2\theta'+\partial_{i}^2v_{gj}\partial_{j}\theta'+2\partial_{i}v_{gj}\partial_{j}\partial_{i}\theta'\\&=\mathbf{V}_{g}\cdot \nabla_{h}\nabla_{h}^2\theta'+\nabla_{h}^2\mathbf{V}_{g}\cdot \nabla_{h}\theta'+2\nabla_{h}\mathbf{V}_{g}:\nabla_{h}\nabla_{h}\theta'
\end{aligned}
$$
其中双点乘使两个并矢相互作用变成标量，双点乘的计算规则可以理解为先计算符号左右靠近的两个矢量点积（即下标保持相同），再计算较远的两个矢量的点积。
$$
\begin{aligned}
			\frac{ \partial  }{ \partial z } (\mathbf{V}_{g}\cdot \nabla_{h}\zeta_{g})&=\frac{ \partial \mathbf{V}_{g} }{ \partial z } \cdot \nabla_{h}\zeta_{g}+\mathbf{V}_{g}\cdot \nabla_{h}\frac{ \partial \zeta_{g} }{ \partial z } 
\end{aligned}
$$
利用
$$
\theta'=\frac{\theta_{0}}{\rho_{0}g}\frac{ \partial p' }{ \partial z } ,\quad \zeta_{g}=\frac{1}{f\rho_{0}}\nabla_{h}^2p'
$$
可得
$$
\begin{aligned}
			\frac{g}{\theta_{0}N^2}\mathbf{V}_{g}\cdot\nabla_{h}\nabla_{h}^2\theta'&=\frac{1}{\rho_{0}N^2}\mathbf{V}_{g}\cdot \nabla_{h}\frac{ \partial  }{ \partial z } (\nabla_{h}^2p')\\&=\frac{f}{N^2}\mathbf{V}_{g}\cdot\nabla_{h}\frac{ \partial  }{ \partial z } \left( \frac{1}{f\rho_{0}}\nabla_{h}^2p' \right)\\&=\frac{f}{N^2}\mathbf{V}_{g}\cdot \nabla_{h}\frac{ \partial  }{ \partial z } \zeta_{g}
\end{aligned}
$$
这就是相消两项中相消的一部分，还有一部分在于
$$
\nabla_{h}^2\mathbf{V}_{g}\cdot \nabla_{h}\theta'=\frac{1}{f\rho_{0}}\nabla_{h}^2(\mathbf{k}\times \nabla_{h}p')\cdot \frac{\theta_{0}}{\rho_{0}g}\frac{ \partial  }{ \partial z } \nabla_{h}p'
$$
利用
$$
\frac{ \partial \mathbf{V}_{g} }{ \partial z } =\frac{1}{f\rho_{0}}\mathbf{k}\times \nabla_{h}p'
$$
即
$$
\nabla_{h}p'=-f\rho_{0}\mathbf{k}\times \frac{ \partial \mathbf{V}_{g} }{ \partial z }
$$
因此上式可化为
$$
\begin{aligned}
					\nabla_{h}^2\mathbf{V}_{g}\cdot \nabla_{h}\theta'&=-\frac{\theta_{0}}{\rho_{0}g}\nabla_{h}^2(\mathbf{k}\times \nabla_{h}p')\cdot \left( \mathbf{k}\times \frac{ \partial \mathbf{V}_{g} }{ \partial z }  \right)\\&=-\frac{\theta_{0}}{\rho_{0}g}\left[ \frac{ \partial v_{g} }{ \partial z }\nabla_{h}^2\frac{ \partial p' }{ \partial y } +\frac{ \partial u_{g} }{ \partial z } \nabla_{h}^2\frac{ \partial p' }{ \partial x }   \right]\\&=-\frac{\theta_{0}}{\rho_{0}g}\frac{ \partial \mathbf{V}_{g} }{ \partial z } \cdot \nabla_{h}\nabla_{h}^2p'\\&=-\frac{f\theta_{0}}{g}\frac{ \partial \mathbf{V}_{g} }{ \partial z } \cdot \nabla_{h}\zeta_{g}
\end{aligned}
$$
可以看到这和原方程中第二项的一部分形式一模一样，于是原方程可以写为
$$
\begin{aligned}
\tilde{\mathbf{L}}w&=-\frac{2g}{\theta_{0}N^2}\nabla_{h}\mathbf{V}_{g}:\nabla_{h}\nabla_{h}\theta'+\frac{2f}{N^2}\frac{ \partial \mathbf{V}_{g} }{ \partial z } \cdot \nabla_{h}\zeta_{g}\\&=-\frac{2g}{\theta_{0}N^2}\nabla_{h}\mathbf{V}_{g}:\nabla_{h}\nabla_{h}\theta'-\frac{2g}{\theta_{0}N^2}\nabla_{h}^2\mathbf{V}_{g}\cdot \nabla_{h}\theta'
\end{aligned}
$$
神奇的是，这里的右式两项也可以融合，写成一项，用Einstein求和法则可得
$$
\begin{aligned}
				RHS&=-\frac{2}{\rho_{0}N^2}\left(\partial_{i}v_{gj}\partial_{j}\partial_{i}\frac{ \partial p' }{ \partial z } +\partial_{i}^2v_{gj}\partial _{j}\frac{ \partial p' }{ \partial z } \right)\\&=-\frac{2}{\rho_{0}N^2}\partial_{i}\left( \partial_{i}v_{gj}\partial_{j}\frac{ \partial p' }{ \partial z }  \right)\\&=-\frac{2}{\rho_{0}N^2}\nabla_{h} \cdot\left( \nabla_{h} \mathbf{V}_{g}\cdot \nabla_{h}\frac{ \partial p' }{ \partial z }  \right)
\end{aligned}
$$
我们可以将其用位温扰动场来表示出来
$$
RHS=-\frac{2g}{\theta_{0}N^2}\nabla_{h} \cdot(\nabla_{h}\mathbf{V}_{g}\cdot \nabla_{h}\theta')
$$
定义矢量$\mathbf{Q}$为（Holton书中用的是按Einstein求和法则展开的形式，二者是等价的）
$$
\mathbf{Q}\equiv-\frac{g}{\theta_{0}N^2}\nabla_{h}\mathbf{V}_{g}\cdot \nabla_{h}\theta'
$$
则通过上面一番恐怖的计算可得垂直运动的方程为
$$
\tilde{\mathbf{L}}w=2\nabla_{h}\cdot \mathbf{Q}
$$
上升运动($w>0$)会发生在$\mathbf{Q}$辐合的地方

天气图上的$\mathbf{Q}$矢量可以通过下面的方法来估计。先将其展开为分量形式
$$
\mathbf{Q}=-\frac{g}{\theta_{0}N^2}\left( \frac{ \partial \mathbf{V}_{g} }{ \partial x } \cdot \nabla_{h}\theta',\frac{ \partial \mathbf{V}_{g} }{ \partial y } \cdot \nabla_{h}\theta' \right)
$$
我们可以将x轴取为沿着等位温线的，且冷空气在x轴左侧，于是有
$$
\frac{ \partial \theta' }{ \partial x } =0,\quad \frac{ \partial \theta' }{ \partial y } <0
$$
于是$\mathbf{Q}$可以写为
$$
\mathbf{Q}=-\frac{g}{\theta_{0}N^2}\left( \frac{ \partial v_{g} }{ \partial x } \frac{ \partial \theta' }{ \partial y } ,\frac{ \partial v_{g} }{ \partial y } \frac{ \partial \theta' }{ \partial y }  \right)
$$
利用地转风水平无散的性质可得
$$
\begin{aligned}
		\mathbf{Q}&=-\frac{g}{\theta_{0}N^2}\left( \frac{ \partial v_{g} }{ \partial x } \frac{ \partial \theta' }{ \partial y } ,-\frac{ \partial u_{g} }{ \partial  x} \frac{ \partial \theta' }{ \partial y }  \right)\\&=-\frac{g}{\theta_{0}N^2}\left|\frac{ \partial \theta' }{ \partial y } \right|\mathbf{k}\times \frac{ \partial \mathbf{V}_{g} }{ \partial x } 
\end{aligned}
$$
因此$\mathbf{Q}$矢量可以通过估计$\mathbf{V}_{g}$沿等位温线的变化（方向保持冷空气在左侧），然后顺时针旋转90度，并乘以系数$\frac{g}{\theta_{0}N^2}\left|\frac{ \partial \theta' }{ \partial y }\right|$得到，下图为此方法的一个例子。
![[估计Q矢量.png]]
这是西风中嵌入的气旋和反气旋扰动，形成压强场中HLH的排列，于是在左边的反气旋气旋间隙中为向南运动，右边则为向北运动。取x轴由西向东沿着等位温线，冷空气在x轴左侧（即北边），在低压中心附近，地转风由南变北，$\frac{ \partial \mathbf{V}_{g} }{ \partial x }$向北，顺时针旋转后指向东。由于冷空气在北边，因此$\nabla_{h}\theta'$指向南边，对应热成风$\frac{g}{g\theta_{0}}\mathbf{k}\times \nabla_{h}\theta'$指向东，因此低压中心附近$\mathbf{Q}$矢量平行于热成风，而在高压中心附近则反过来，$\mathbf{Q}$矢量反平行于热成风。

根据上面的分析，可以看到在冷空气在低压槽西部平流（即左边的HL间隙）处，$\mathbf{Q}$矢量辐散，$\nabla \cdot \mathbf{Q}>0$，因此$w<0$，发生下降运动；在热空气在低压槽东部平流处，$\mathbf{Q}$矢量辐合，$\nabla \cdot \mathbf{Q}<0$，因此$w>0$，发生上升运动。

$\omega$-方程还有另一种形式，也就是在最后一步中不融合两项所得到的
$$
\tilde{\mathbf{L}}w=-\frac{2g}{\theta_{0}N^2}\nabla_{h}\mathbf{V}_{g}:\nabla_{h}\nabla_{h}\theta'+\frac{2f}{N^2}\frac{ \partial \mathbf{V}_{g} }{ \partial z } \cdot \nabla_{h}\zeta_{g}
$$
记$D=-\frac{2g}{\theta_{0}N^2}\nabla_{h}\mathbf{V}_{g}:\nabla_{h}\nabla_{h}\theta'$称为变形项（deformation term），并代入热成风
$$
\mathbf{V}_{T}=\frac{ \partial \mathbf{V}_{g} }{ \partial z }
$$
可得
$$
\tilde{\mathbf{L}}w=\frac{2}{N^2}\mathbf{V}_{T}\cdot \nabla_{h}(f\zeta_{g})+D
$$
这称为QG $\omega$-方程的Sutcliffe形式。如果忽略变形项，可以看到上升的空气来自于气旋性涡度异常中心处的热成风平流，尽管这个近似忽略了有可能很重要的变形项，但涡度的热成风平流在天气图中是容易看出的。比如还是对于速度随高度增大的西风急流，其热成风指向东，于是在PV异常下游（右侧），热成风的地转涡度平流小于0，因此$w>0$，发生上升运动，这结果与通过$\mathbf{Q}$矢量分析的一致。这种方法可以快速地通过天气图判断出上升和下降气流所在的位置。

我们继续以西风场中的PV异常为例子研究$w$- Thinking分析动力学的方法。(这边书里写的非常含糊，很多地方我没理解他的逻辑在哪，只能按自己粗糙的理解写了)
![[垂直运动.png]]
在PV异常下游，$\mathbf{Q}$矢量辐合，于是发生向上的垂直运动，上游则反之。出于连续性方程的要求，这种垂直运动会和非地转风的环流联系起来。现在考虑QG PV极大值附近的位温倾向
$$
\frac{ \partial \theta' }{ \partial t } =-\mathbf{V}_{g}\cdot \nabla_{h}\theta'-w \frac{d\theta_{0}}{dz}
$$
在PV极大值右侧，由于是气旋性涡度，故$\mathbf{V}_{g}$指向北，而根据前面所规定的，冷空气在x轴左侧，即北侧，因此$\frac{ \partial \theta' }{ \partial y }<0$，故位温平流项>0，位温在增加。在位温变化和垂直运动的综合影响下，地转平衡与流体静力平衡逐渐被打破，而后又逐渐被非地转的环流所修复。

我们以流体静力平衡打破后的调节为例。流体静力平衡可以等价为热成风关系
$$
\frac{ \partial \mathbf{V}_{g} }{ \partial z } =\frac{g}{f\theta_{0}}\mathbf{k}\times \nabla_{h}\theta'
$$
由于上面所说的，在PV异常下游，位温的平流项导致该处$\theta'$增大，于是$\frac{ \partial \theta' }{ \partial x }$增大，打破了热成风关系。而其恢复机制可以从两个角度看出来，一是PV异常下游的上升运动$w>0$对平流项造成的位温增加起到了抑制作用，这很容易从上面的位温倾向方程看出来。二是非地转风队地转风的调节作用，我们考虑前面一级近似下的动力学方程
$$
\left(\frac{\partial}{\partial t}+u^{(0)}\frac{\partial}{\partial x}+v^{(0)}\frac{\partial}{\partial y}\right)v^{(0)}+\beta_0yu^{(0)}+f_0u^{(1)}=-\frac{1}{\rho_0}\frac{\partial p^{(1)}}{\partial y}
$$
按$f$-平面近似，以及取$p^{(1)}=0$，可得
$$
\frac{Dv_{g}}{Dt_{g}}+fu_{a}=0
$$
在PV异常上部，$u_{a}<0$，于是$v_{g}$随时间，从而使热成风关系$\frac{ \partial v_{g} }{ \partial z }=\frac{g}{f\theta_{0}}\frac{ \partial \theta' }{ \partial x }$恢复。

在准地转方程中，调节过程是瞬时的，因此大气总是准地转平衡和准流体静力平衡的。回到位涡视角下的动力学来看，可以说非地转环流对于维持地转平衡和流体静力平衡是必要的。

下面以西风急流带（即上气旋下反气旋，中间夹的就是西风急流带）为例展示垂直运动以及连续性方程引发的非地转环流。如下图
![[西风急流带-非地转环流.png]]
气旋和反气旋一上一下分割出四个区域，如果定义入口是风速开始急剧变大的地方（$\frac{ \partial u_{g} }{ \partial x }>0$），出口是风速开始急剧变小的地方，那么可以利用$\mathbf{Q}$矢量推知上升运动发生的两个区域分别对应于“右入口”和“左出口”，下降运动对应于“左入口”和“右出口”。图6.17中左边对应的是y=-700km处（即出口）垂直截面上的非地转环流，右边则是y=700km（即入口）处垂直截面上的非地转环流，这些环流使得流入急流的风（即入口）加速，流出的则相反，同时也增加了纬向风在急流入口区域的垂直切变（这部分不知道是为什么？）。注意到，在急流入口区域，由于急流左侧对应冷空气平流，右侧对应暖空气平流，因此暖空气上升，冷空气下降，由于暖空气密度偏小，势能下降转变为急流带中的强动能，这从能量角度对急流带的产生与入口出口区分做了解释。

## 1.5 斜压扰动的理想模型
我们有了两个解释准地转动力学的框架，现在我们用这个框架来理解增强的热带外气旋，这一过程称为cyclogenesis。下图展示了这一发展过程中的基本元素。（这部分可以与第一部分中定性描述的现象对应起来）
>这部分内容有点晦涩，暂时只能粗略地翻译一下，很多地方都没有理解书里说的是什么意思

![[热带外气旋的增强.png]]
几乎所有的气旋发展都是由对流层顶波动向下传播造成的上层大气扰动导致的，这种向下传播与气旋性的PV异常相联系，在上图中用“+”号表示。前面讨论过，气旋性的PV异常对应于低压、气旋性风场、上暖下冷。在地表，发展中的气旋产生于上空PV异常下游（东侧）的较暖区域，这种地表的暖异常区域也对应气旋性的PV异常、低压，且随高度升高而衰减。总体来说，这种PV异常产生的压强场随高度西倾，而位温场则随高度东倾。

从PV-$\theta$视角来看，由于上空PV异常的诱导流通过暖空气平流使表面暖异常幅度增大，因此表面低压中心发展于上空PV异常东侧。注意这种机制需要背景位温场随往极地方向而减小，而根据热成风西风向上速度增加。类似的，在上空PV异常处给定一个沿极地方向增加的PV，则PV异常的幅度会因为地表暖异常诱导的赤道方向的PV平流而增大。（看不懂？？）

从高度倾向（即压强倾向方程）的视角来看，我们知道由于背景西风的QG PV平流使上空PV异常东侧压强下降。此外，在地表的发展中低压中心西侧的向上衰减的冷空气对流导致了上层低压的加剧。通过类似的原因，可以理解表面高压在PV异常的上游的发展。

从$w$视角，由于与上升运动相联系的涡旋拉伸发生在上空PV异常的东侧，通过前面的分析可以知道背景切变西风中的PV异常导致了非地转环流的产生，而该处附近发生辐合，因此低层涡度在该处增加。通过热力学能量方程可知，上升运动倾向于降低附近空气的温度，从而部分地抵消暖空气平流造成的升温。类似的，可以通过空气下降导致的涡度压缩来理解表面高压在PV异常上游的发展。

现在我们可以往干的cyclogenesis的描述中加入云与降雨的贡献。由于相变潜热的释放，在上升区域会成云。从$w$视角来看，这种释放的热量通过加热气旋附近的空气进一步增强了气旋。从PV-$\theta$视角来看，这种热量的释放意味着PV随空气块的运动是不守恒的，Chapter4中讨论有热源的Ertel位涡方程时就指出，气旋性位涡发生在热量最大值下方，这个新的气旋性PV位于地表气旋的上方，通过PV反演和叠加原理可知，这导致了压强的减少和气旋性环量的增加。

最后，cyclogenesis过程当背景静态稳定性降低时会被加速，这可以通过在PV-$\theta$视角下Rossby深度的增大来理解，当Rossby深度增加时，地表的诱导流会被增强。从$w$视角来看，更小的静态稳定性也导致了更强的非地转环流。相反，当背景静态稳定性增大时，cyclogenesis过程减慢，这就是为什么主要的热带外气旋风暴轨迹在太平洋和大西洋的西部被发现，海洋上空的高热量使得在冬天水温仍较暖，从而可以加热海水上面的低层大气，减弱对流层的静态稳定性，从而加速cyclogenesis的发展。
## 1.6 压强坐标形式的准地转方程
天气数据和天气图往往是用压强坐标来表达的，这一方面是历史因素，另一方面也是因为在压强坐标下密度$\rho$消失了，简化了方程的未知量。压强形式的准地转方程当然可以用前面的小参数方法来推导，这里我们介绍一种新的方式，即通过更为纯粹的矢量分析方法来推导，这种方法可以更清楚地看出小参数方法的本质。(这部分内容借鉴自@mosekyo的文章)

考虑大气的基本方程
$$
\begin{cases}
					&\frac{D\mathbf{V}}{Dt}+f\mathbf{k}\times \mathbf{V}+\nabla \phi=0\\&\nabla \cdot \mathbf{V}+\frac{ \partial \omega }{ \partial p } =0\\&\frac{ \partial \phi }{ \partial p } =-\frac{RT}{p}\\&\left( \frac{ \partial  }{ \partial t } +\mathbf{V}\cdot \nabla  \right)T-\frac{\sigma p}{R}\omega=\frac{J}{c_{p}}
\end{cases}
$$
其中大气静态稳定性参数$\sigma=-\frac{RT}{p}\frac{ \partial \ln \theta }{ \partial p }$，四个变量$\mathbf{V},w,\phi,T$，四个方程。

现在做准地转近似，即$\mathbf{V}=\mathbf{V}_{g}+V_{a}$，$\frac{V_{a}}{V_{g}}\sim Ro \sim0.1$，其中$\mathbf{V}_{g}=\frac{1}{f_{0}}\mathbf{k}\times \nabla \phi$，再做$\beta$-平面近似，代入可得
$$
\begin{cases}
&(\frac{ \partial  }{ \partial t } +\mathbf{V _{g}\cdot \nabla })\mathbf{V _{g}}+f_{0}\mathbf{k}\times \mathbf{V}_{a}+\beta y\mathbf{k}\times \mathbf{V}_{g}=0\\&\nabla \cdot \mathbf{V}_{a}+\frac{ \partial \omega }{ \partial p } =0\\&\frac{ \partial \phi }{ \partial p } =-\frac{RT}{p}\\&\left( \frac{ \partial  }{ \partial t } +\mathbf{V}_{g}\cdot \nabla  \right)T-\frac{\sigma p}{R}\omega=\frac{J}{c_{p}}
\end{cases}
$$
这和前面小参数方法推导得到的一级近似方程完全等价，和前面处理方法类似，消去不感兴趣的非地转量，用$\mathbf{k}\times \nabla$作用于第一个方程可得
$$
\left( \frac{ \partial  }{ \partial t } +\mathbf{V}_{g}\cdot \nabla \right)(\zeta_{g}+f)=-f_{0}\nabla \cdot \mathbf{V}_{a}=f_{0}\frac{ \partial \omega }{ \partial p }
$$
将流体静力平衡方程代入热力学能量方程消去T可得
$$
\left( \frac{ \partial  }{ \partial t } +\mathbf{V}_{g}\cdot \nabla \right)\left( \frac{ \partial \phi }{ \partial p }  \right)+\sigma \omega=-\frac{\kappa J}{p}
$$
这两个方程合起来就是准地转模式方程组，其中$\zeta_{g}=\frac{1}{f_{0}}\nabla^2\phi$，定义$\chi=\frac{ \partial \phi }{ \partial t }$，消去$\omega$可得
$$
\begin{aligned}\left[\nabla^{2}+\frac{\partial}{\partial p}{\left(\frac{f_{0}^{2}}{\sigma}\frac{\partial}{\partial p}\right)}\right]\chi=&-\left.f_0\boldsymbol{V}_g\cdot\nabla\left[\frac{1}{f_0}\nabla^2\Phi+f\right]\right.\\&-\frac{\partial}{\partial p}\left[\frac{f_0^2}{\sigma}\left(\boldsymbol{V}_g\cdot\nabla\left(\frac{\partial\Phi}{\partial p}\right)+\frac{\kappa J}{p}\right)\right]\end{aligned}
$$
这称为准地转-$\chi$方程（准地转位势倾向方程），右侧为通过$\phi$瞬时计算得到的诊断量，坐标$\phi$是预测量，这就是$\phi$场的预报方程。

上面方程两边同除$f_{0}$后，左边可以写为
$$
\left[\frac{1}{f_0}\nabla^2+\frac{\partial}{\partial p}\left(\frac{f_0}{\sigma}\frac{\partial}{\partial p}\right)\right]\chi=\frac{\partial}{\partial t}\left[\frac{1}{f_0}\nabla^2\Phi+f+\frac{\partial}{\partial p}\left(\frac{f_0}{\sigma}\frac{\partial\Phi}{\partial p}\right)\right]+\frac{\partial}{\partial p}\left(\frac{f_0}{\sigma^2}\frac{\partial\sigma}{\partial t}\frac{\partial\Phi}{\partial p}\right)
$$
右边第二项可以写为（这里认为$\sigma$与水平方向x，y无关）
$$
\begin{aligned}\frac{\partial}{\partial p}{\left[\frac{f_{0}}{\sigma}\boldsymbol{V}_{\boldsymbol{g}}\cdot\nabla\left(\frac{\partial\Phi}{\partial p}\right)\right]}&=\frac{f_{0}}{\sigma}\frac{\partial\boldsymbol{V}_{g}}{\partial p}\cdot\nabla\frac{\partial\Phi}{\partial p}+\boldsymbol{V}_{g}\cdot\nabla\left[\frac{\partial}{\partial p}{\left(\frac{f_{0}}{\sigma}\frac{\partial\Phi}{\partial p}\right)}\right]\\&=\frac{1}{\sigma}\left(\boldsymbol{k}\times\nabla\frac{\partial\Phi}{\partial p}\right)\cdot\nabla\frac{\partial\Phi}{\partial p}+\boldsymbol{V}_g\cdot\nabla\left[\frac{\partial}{\partial p}\left(\frac{f_0}{\sigma}\frac{\partial\Phi}{\partial p}\right)\right]\\&=\boldsymbol{V}_g\cdot\nabla\left[\frac{\partial}{\partial p}{\left(\frac{f_0}{\sigma}\frac{\partial\Phi}{\partial p}\right)}\right]\end{aligned}
$$
>这些公式包括后面的都可以用前面使用的Einstein求和法则轻松地推导出来

于是准地转-$\chi$方程可以写为
$$
\left(\frac{\partial}{\partial t}+\boldsymbol{V}_g\cdot\nabla\right)q=-\frac{\partial}{\partial p}\left(\frac{f_0\kappa J}{\sigma p}+\frac{f_0}{\sigma^2}\frac{\partial\sigma}{\partial t}\frac{\partial\Phi}{\partial p}\right)
$$
其中准地转位涡$q$为
$$
q=\frac{1}{f_{0}}\nabla^2\phi+f+\frac{ \partial  }{ \partial p } \left( \frac{f_{0}}{\sigma}\frac{ \partial \phi }{ \partial p }  \right)
$$
当绝热$J=0$且大气层结不随时间变化$\frac{ \partial \sigma }{ \partial t }=0$时，可得准地转位涡守恒。

在推导准地转-$\chi$方程时消去的是$\omega$，我们也可以消去预测量$\chi$，得到的就是$\omega$-方程
$$
\begin{aligned}\left[\nabla^{2}(\sigma\cdot)+f_{0}^{2}\frac{\partial^{2}}{\partial p^{2}}\right]\omega=&-\frac{\partial}{\partial p}\left[-f_0\boldsymbol{V}_g\cdot\nabla\left(\frac{1}{f_0}\nabla^2\Phi+f\right)\right]\\&-\nabla^2\left[-\boldsymbol{V}_g\cdot\nabla\left(-\frac{\partial\Phi}{\partial p}\right)+\frac{\kappa J}{p}\right]\end{aligned}
$$
用和前面类似的分析方法可以将$\omega$-方程右侧部分抵消的两项消掉，得到
$$
\begin{aligned}\left[\nabla^{2}(\sigma\cdot)+f_{0}^{2}\frac{\partial^{2}}{\partial p^{2}}\right]\omega=&2\frac{\partial\boldsymbol{V}_{g}}{\partial p}\cdot\nabla\nabla^{2}\Phi-2\nabla\boldsymbol{V}_{g}:\nabla\nabla\frac{\partial\Phi}{\partial p}\\&+f_0\frac{\partial\boldsymbol{V}_g}{\partial p}\cdot\nabla f-\nabla^2\frac{\kappa J}{p}\end{aligned}
$$
这里多出了后面两项是因为没有采用$f$-平面近似以及绝热近似，如果采取这两个近似就可以得到stucliffe形式的$\omega$-方程，更进一步的，还可以类似的得到$\mathbf{Q}$矢量
$$
\mathbf{Q}=\frac{1}{\sigma}\nabla \mathbf{V}_{g}\cdot \nabla \frac{ \partial \phi }{ \partial p } =-\frac{R}{\sigma p}\nabla \mathbf{V}_{g}\cdot \nabla T
$$
$$
\left( \nabla^2+\frac{f_{0}^2}{\sigma}\frac{ \partial ^2 }{ \partial p^2 }  \right)\omega=-2\nabla \cdot \mathbf{Q}
$$
可以与z坐标下的各种方程对应上。