# 0. 前置知识

## 0.1 随机过程与布朗运动

**随机过程**就是“随时间演化的随机变量族”。  
也就是说，对每个时刻 $t$，系统状态 $X_t$ 都是一个随机变量；把所有时刻连起来，就得到随机过程：
$$
\{X_t\}_{t\ge 0}.
$$

它和普通随机变量的区别在于：

- 普通随机变量只描述“一次抽样结果”
- 随机过程描述“随时间变化的一串随机结果”

---

**布朗运动**是最经典的随机过程之一。  
它最初来自对悬浮微粒在液体中不规则运动的观察。微粒不断受到周围大量分子的随机碰撞，于是轨迹表现为杂乱无规则的抖动。

从建模角度看：

- 单个粒子的轨迹是随机的
- 但大量重复实验的统计规律是可描述的

这就引出两种互补语言：

1. **轨迹语言**：研究单条随机轨道如何演化  
2. **分布语言**：研究许多轨道组成的概率密度如何演化  

Langevin 方程属于第一类，Fokker–Planck 方程属于第二类。

---

## 0.2 Wiener 过程与白噪声

### 0.2.1 Wiener 过程

Wiener 过程通常记作 $W_t$，它是理想布朗运动的数学模型。它满足以下四个基本性质：

#### (1) 初值
$$
W_0 = 0
$$

#### (2) 独立增量
对任意不重叠时间区间，增量彼此独立。  
例如：
$$
W_{t_2}-W_{t_1}, \qquad W_{t_4}-W_{t_3}
$$
若区间 $[t_1,t_2]$ 与 $[t_3,t_4]$ 不重叠，则这两个增量独立。

#### (3) 平稳高斯增量
任意小时间步 $\Delta t$ 上，
$$
W_{t+\Delta t}-W_t \sim \mathcal N(0,\Delta t)
$$
也就是说：

- 均值为 $0$
- 方差为 $\Delta t$

#### (4) 路径连续但几乎处处不可导
$W_t$ 关于 $t$ 连续，但几乎处处不可导。  
所以它可以连续变化，却“抖得太厉害”，没有通常意义下的切线。

---

在随机微积分中，常把它写成微分形式：
$$
dW_t \sim \mathcal N(0,dt)
$$
并满足
$$
\langle dW_t\rangle = 0, \qquad \langle dW_t^2\rangle= dt.
$$

> 注：$\langle dW_t^2\rangle=\mathbb{E}[X^2] = \mathrm{Var}(X) + (\mathbb{E}[X])^2 = dt + 0  = dt$

这里第二个式子非常关键，因为它体现了Itô 微积分（随机微积分）与普通微积分最根本的区别。

---

### 0.2.2 白噪声

白噪声常记作 $\xi(t)$，形式上常写成
$$
\xi(t)=\frac{dW_t}{dt}.
$$

但要特别强调：

> 这不是普通函数意义下的导数，而是**广义函数 / 分布意义下**的写法。

Wiener 过程是白噪声的时间积分：
$$
W_t = \int_0^t \xi(s)\,ds.
$$

所以可以把它们理解成：

- $\xi(t)$：随机“速度”或随机驱动
- $W_t$：其积分后的随机位移

---

此外，白噪声满足：

$$
\langle \xi(t)\rangle = 0,
$$
$$
\langle \xi(t)\xi(t')\rangle = \delta(t-t').
$$

若带噪声强度，也常写作
$$
\langle \xi(t)\xi(t')\rangle = 2D\,\delta(t-t').
$$

这里 $\delta(t-t')$ 的含义是：

- 只有在“同一时刻”才相关
- 不同时刻之间不相关
- 严格说只能放在积分里理解（作为核）

例如
$$
\int f(t')\,\delta(t-t')\,dt' = f(t).
$$

> 随机过程的自相关函数（物理）是协方差核（数学），白噪声的自相关函数/协方差核是dirac函数（作用之后返回过程本身）

---

**总结：白噪声两种等价的刻画方式**

- **从 Wiener 过程出发**，把白噪声看成 $W_{t}$ 的广义导数
- **从协方差结构出发**，把白噪声定义为满足 $\langle \xi(t) \xi(t') \rangle = \delta(t - t')$ 的零均值高斯广义过程

---

**白噪声为什么叫“白”**？

这是和白光类比：

- 白光包含各种频率成分
- 白噪声对应近似平坦的功率谱

所以“白”表示**频谱平坦**，不是颜色。

## 0.3 Langevin 方程的基本形式

Langevin 方程描述单个随机系统轨迹的演化。最常见的一维形式是
$$
\frac{dx}{dt}=A(x)+B(x)\,\xi(t),
$$
或者更严格地写成随机微分方程：
$$
dx=A(x)\,dt+B(x)\,dW_t.
$$

其中：

- $A(x)$：漂移项，表示确定性趋势
- $B(x)$：噪声振幅
- $B(x)\xi(t)$：随机扰动

如果 $B(x)$ 为常数，叫**加性噪声**；  
如果 $B(x)$ 依赖于状态 $x$，叫**乘性噪声**。

---

Langevin 方程的物理意义是：

- $A(x)$ 负责“平均往哪里走”
- 噪声项负责“在平均趋势之上不断抖动”

所以它给出的是**单条随机轨迹**的动力学描述。

---

## 0.4 概率密度、概率流与连续性方程

设 $P(x,t)$ 表示系统在时刻 $t$ 出现在位置 $x$ 附近的概率密度，则有归一化条件：
$$
\int P(x,t)\,dx=1.
$$

如果概率守恒，那么概率密度应满足连续性方程：
$$
\frac{\partial P}{\partial t}=-\frac{\partial J}{\partial x},
$$
其中 $J(x,t)$ 称为**概率流**。

它和流体中的质量守恒方程完全类比：
$$
\frac{\partial \rho}{\partial t} + \nabla \cdot \mathbf{j} = 0  
$$
- 流体密度 $\rho \leftrightarrow P$ 概率密度
- 质量流 $\mathbf j \leftrightarrow J$ 概率流

物理上就是：  
某点概率密度的变化，来源于概率从那里流入或流出。

---

# 1. Fokker–Planck 方程的提出

## 1.1 为什么要从单条轨迹转向概率分布

Langevin 方程描述单条随机轨迹，但随机系统的一个根本特点是：

> 你通常无法精确预测某一条轨迹，只能预测许多重复实验的统计规律。

因此，除了研究单个粒子怎么走，还需要研究：

- 某时刻粒子更可能出现在哪里
- 概率分布如何扩散、漂移、趋于稳态
- 噪声如何改变整体统计行为

这就从“轨迹语言”转向了“分布语言”。

---

## 1.2 Fokker–Planck 方程的标准形式

对于一维过程，Fokker–Planck 方程通常写作
$$\frac{\partial P(x,t)}{\partial t}=-\frac{\partial}{\partial x}\bigl[A(x)P(x,t)\bigr]+\frac{\partial^2}{\partial x^2}\bigl[D(x)P(x,t)\bigr].$$
其中：

- 漂移项 $-\frac{\partial}{\partial x}(A P)$：表示概率分布被确定性动力学推动着整体移动。
- 扩散项 $\frac{\partial^2}{\partial x^2}(D P)$：表示随机噪声导致概率分布在空间中摊开。

所以它描述的是：

> **确定性趋势 + 随机扩散** 共同决定概率分布的演化。

---
**如果只有漂移项**：

如果 $D(x) = 0$，方程变成  
$$\frac{\partial P}{\partial t} = -\frac{\partial}{\partial x}\big(A(x)P\big)$$
退化为连续性方程

它实际上就是  
$$\frac{\partial P}{\partial t} + \frac{\partial J}{\partial x} = 0$$
其中概率流
$$J(x,t) = A(x)P(x,t)$$

这意味着：  
- 没有随机性
- 所有粒子完全按确定性动力学运动
- 概率密度在相空间中只进行输运

对应的粒子动力学就是  
$$
\frac{dx}{dt} = A(x)
$$
这是一条普通微分方程。  

---
**如果只有扩散项**：

如果 $A(x) = 0$，方程变成
$$\frac{\partial P}{\partial t} = \frac{\partial^2}{\partial x^2}(D(x)P)$$
如果扩散系数是常数
$$D(x) = D$$
就得到经典的扩散方程：
$$\frac{\partial P}{\partial t} = D \frac{\partial^2 P}{\partial x^2}$$
这意味着：
- 没有确定性力
- 只有随机噪声
- 概率分布为一个不断变宽的高斯

粒子动力学是
$$dx = \sqrt{2D}\, dW_t$$
也就是纯布朗运动

# 2. 从 Langevin 方程到 Fokker–Planck 方程

## 2.1 随机增量的一阶矩与二阶矩

考虑 Itô 形式随机微分方程：
$$
dx=A(x)\,dt+B(x)\,dW_t.
$$
在一个很小时间步 $dt$ 内，粒子增量为
$$
\Delta x = A(x)\,dt + B(x)\,dW_t.
$$
因为
$$
\langle dW_t\rangle=0,\qquad \langle dW_t^2\rangle=dt,
$$
所以有：

- **一阶矩**
$$
\langle \Delta x\rangle = A(x)\,dt
$$
- **二阶矩**
$$
\langle (\Delta x)^2\rangle = B^2(x)\,dt
$$

> $(\Delta x)^2 = A^2 dt^2 + 2AB\, dt\, dW_t + B^2 (dW_t)^2$

---
## 2.2 漂移系数与扩散系数

由上面两式可见：

- **一阶矩**给出平均漂移
- **二阶矩**给出扩散强度

定义
$$
M^{(1)}(x)=\lim_{dt\to 0}\frac{\langle \Delta x\rangle}{dt}=A(x),
$$
$$
M^{(2)}(x)=\lim_{dt\to 0}\frac{\langle (\Delta x)^2\rangle}{dt}=B^2(x).
$$
于是扩散系数常定义为
$$
D(x)=\frac{1}{2}M^{(2)}(x)=\frac{B^2(x)}{2}.
$$
则 Planck 方程可重新写作：
$$  
\frac{\partial P}{\partial t} = -\frac{\partial}{\partial x} \left( M^{(1)} P \right) + \frac{1}{2} \frac{\partial^2}{\partial x^2} \left( M^{(2)} P \right) 
$$
> 记：$D(x) = \frac{1}{2} M^{(2)}(x)$

这里的 $M$ 指的是 Kramers–Moyal 矩/跃迁矩。
 
---

**两个方程联系的桥梁**：

Langevin 方程描述的是单条随机轨迹：
$$
dx = A(x)\,dt + B(x)\,dW_t.
$$
Fokker–Planck 方程描述的是概率分布 $P(x,t)$ 的演化。

如果现在已知单条轨迹在很短时间内怎么走，怎样推出大量粒子的概率分布怎么变？

答案就是看短时间增量 $\Delta x$ 的统计性质。

---

## 2.3 Fokker–Planck 方程的推导思路

推导的核心思想是：

1. **局域短时增量** $\Delta x$，由 Langevin 方程给出  
2. 一阶矩 $\langle \Delta x\rangle = A(x)\,dt$，决定概率的平均流动  
3. 二阶矩 $\langle (\Delta x)^2\rangle = B^2(x)\,dt$，决定概率的扩散  
4. 再结合概率守恒，就得到 Fokker–Planck 方程

最终：

$$\frac{\partial P}{\partial t}=-\frac{\partial}{\partial x}\bigl[A(x)P\bigr]+\frac{1}{2}\frac{\partial^2}{\partial x^2}\bigl[B^2(x)P\bigr].$$

换句话说：

> Langevin 方程描述微观轨迹的局域统计；  
> Fokker–Planck 方程把这些局域统计转写成整体概率密度的演化方程。

---

# 3. 从 Chapman–Kolmogorov 到 Kramers–Moyal 展开

## 3.1 Chapman–Kolmogorov 方程

Chapman–Kolmogorov 方程是马尔可夫（Markov）过程的基本拼接关系。  

设系统从 $x_0,t_0$ 演化到 $x,t+\tau$，可以先经过中间状态 $x',t$，则

$$P(x,t+\tau|x_0,t_0)=\int P(x,t+\tau|x',t)\,P(x',t|x_0,t_0)\,dx'.$$

若写成短时转移核形式，则

$$P(x,t+\tau)=\int W(x|x',\tau)\,P(x',t)\,dx'.$$

其中 $W(x|x',\tau)$ 是短时转移概率核，表示在一个很短时间 $\tau$ 内，系统若现在在 $x'$，之后出现在 $x$ 的转移概率密度。

这就是后续展开的起点。

---
## 3.2 Kramers–Moyal 展开

对短时转移概率做小时间展开，可以定义跃迁矩：

$$M^{(n)}(x,\tau)=\int r^n W(x+r|x,\tau)\,dr.$$

进一步定义 Kramers–Moyal 系数：

$$D^{(n)}(x)=\frac{1}{n!}\lim_{\tau\to 0}\frac{M^{(n)}(x,\tau)}{\tau}.$$

则概率密度满足 Kramers–Moyal 级数：

$$\frac{\partial P(x,t)}{\partial t}=\sum_{n=1}^{\infty}\left(-\frac{\partial}{\partial x}\right)^n\left[D^{(n)}(x)P(x,t)\right].$$

这是一种非常一般的连续马尔可夫过程描述。

---
## 3.3 具体推导

设概率密度满足短时传播公式
$$
P(x,t+\tau)=\int W(x|x',\tau)\,P(x',t)\,dx',
$$
其中 $W(x|x',\tau)$ 是从 $x'$ 在时间 $\tau$ 后到达 $x$ 的短时转移概率密度。

### 1. 引入短时跃迁变量

令
$$
r=x-x', \qquad x'=x-r, \qquad dx'=dr,
$$
则上式写成
$$
P(x,t+\tau)=\int W(x|x-r,\tau)\,P(x-r,t)\,dr.
$$
这表示：系统先在 $x-r$，再在短时间 $\tau$ 内跳跃 $r$ 到达 $x$。

### 2. 对被积函数作小跃迁展开

当 $\tau\to 0$ 时，主要贡献来自小跃迁 $r$，于是对 $P(x-r,t)$ 在 $x$ 处按位移 $-r$ 作 Taylor 展开：

$$P(x-r,t)=\sum_{n=0}^\infty \frac{(-r)^n}{n!}\,\partial_x^n P(x,t).$$

代入得

$$P(x,t+\tau)=\int W(x|x-r,\tau)\sum_{n=0}^\infty \frac{(-r)^n}{n!}\,\partial_x^n P(x,t)\,dr.$$

交换积分与求和：

$$P(x,t+\tau)=\sum_{n=0}^\infty\frac{(-1)^n}{n!}\partial_x^n P(x,t)\int r^n W(x|x-r,\tau)\,dr.$$

### 3. 定义短时跃迁矩

定义第 $n$ 阶跃迁矩

$$M^{(n)}(x,\tau)=\int r^n\,W(x+r|x,\tau)\,dr.$$

于是可将传播式整理为标准形式

$$P(x,t+\tau)=\sum_{n=0}^\infty\left(-\frac{\partial}{\partial x}\right)^n\left[\frac{M^{(n)}(x,\tau)}{n!}P(x,t)\right].$$

这个式子表示：时刻 $t + \tau$ 的概率分布，是时刻 $t$ 的分布经过“一次短时间随机跳跃”后得到的结果；级数表示考虑短时间内所有可能发生的的小跃迁，把这些跃迁对分布的影响，按跃迁长度 $r$ 的各阶矩展开，最后加权得到新的分布

> **矩**本质上表示对随机变量的幂 $X^n$ 求平均，用来刻画分布的形状。一阶矩是平均值，二阶矩反应分布的展宽程度，三阶矩、四阶矩则和偏斜、尖峭等高阶形状有关

把式子展开到前三项：

$$P(x, t + \tau) = P(x, t) - \frac{\partial}{\partial x} \left[ M^{(1)}(x, \tau) P(x, t) \right] + \frac{1}{2} \frac{\partial^2}{\partial x^2} \left[ M^{(2)}(x, \tau) P(x, t) \right] + \cdots$$

#### 1. $n = 0$：什么也不发生（概率守恒）

这一项是
$$
M^{(0)}(x,\tau) P(x,t)
$$
而
$$
M^{(0)}(x,\tau) = \int W(x+r|x,\tau) dr = 1
$$
所以
$$
n = 0 \Rightarrow P(x,t)
$$

这一项表示：

- 如果暂时忽略所有跃迁效应，概率分布保持原样。

或者说：

- 这只是保证总概率为 $1$ 的归一化项。

它不描述任何动力学变化，只是传播公式里的“基底”。

#### 2. $n = 1$：平均位移 $\to$ 漂移

这一项是

$$\frac{\partial}{\partial x} \left[ M^{(1)}(x, \tau) P(x, t) \right]$$

其中

$$M^{(1)}(x, \tau) = \int r \, W(x + r | x, \tau) \, dr$$

是位移 $r$ 的平均值。

- 如果 $M^{(1)} > 0$，说明粒子平均会向 $x$ 正方向移动。
- 如果 $M^{(1)} < 0$，说明粒子平均会向 $x$ 负方向移动。

所以，一阶矩控制概率分布的整体漂移。

这对应 Fokker–Planck 方程里的漂移项：

$$\partial_x (A(x) P)$$

其中

$$A(x) = \lim_{\tau \to 0} \frac{M^{(1)}(x, \tau)}{\tau}$$

#### 3. $n = 2$：随机展宽 $\to$ 扩散

这一项是

$$\frac{1}{2} \frac{\partial^2}{\partial x^2} \left[ M^{(2)}(x,\tau) P(x,t) \right]$$

其中

$$M^{(2)}(x,\tau) = \int r^2 W(x+r|x,\tau) dr$$

是均方位移。

即使平均位移为 0，粒子仍然会向左右随机跳动。

于是：

- 有的粒子向左  
- 有的粒子向右  

概率分布就会变宽。

这就是扩散。

所以，二阶矩控制概率分布的展宽（扩散）。

在连续极限里：

$$D(x) = \frac{1}{2} \lim_{\tau \to 0} \frac{M^{(2)}(x,\tau)}{\tau}$$

得到 Fokker–Planck 方程的扩散项：

$$\partial_x^2 (D P)$$

#### 4. 更高阶 $n \geq 3$

这些项表示更复杂的跃迁统计：

例如

- 三阶矩 $\to$ 跃迁分布不对称（偏斜）  
- 四阶矩 $\to$ 尾部更重（尖峭）

如果过程是连续扩散型过程（如布朗运动），这些高阶项在极限里**通常消失**。

于是只剩前两项：  

$$\frac{\partial P}{\partial t} = -\partial_x (A P) + \partial^2_x (D P)$$

这就是 Fokker–Planck 方程。

#### 5. Pawula 定理

若一个 Kramers–Moyal 展开对应某个马尔可夫过程的概率密度演化，并且要求演化保持概率密度的正性。

如果所有高阶系数从某一阶开始严格为零，那么除了二阶截断外，不允许在更高有限阶截断。

也就是说：

- 可以是二阶：  $D^{(n)} = 0, \quad n \geq 3$
- 也可以是无穷阶：$D^{(n)} \text{ 有无限多个非零}$
- 但不能是：$D^{(3)} \neq 0, \quad D^{(n)} = 0 \text{ for } n \geq 4$，或类似的有限高阶截断

因为这种有限高阶截断一般不保持正性。

### 4. 对左边作小时间展开

回到：

$$P(x,t+\tau)=\sum_{n=0}^\infty\left(-\frac{\partial}{\partial x}\right)^n
\left[\frac{M^{(n)}(x,\tau)}{n!}P(x,t)\right].$$

对左边作 Taylor 展开：

$$P(x,t+\tau)=P(x,t)+\tau\,\frac{\partial P(x,t)}{\partial t}+o(\tau).$$

代入并消去两边共同的 $P(x,t)$，得到

$$\tau\,\frac{\partial P}{\partial t}=\sum_{n=1}^\infty\left(-\frac{\partial}{\partial x}\right)^n\left[\frac{M^{(n)}(x,\tau)}{n!}P(x,t)\right]+o(\tau).$$

两边除以 $\tau$，再取极限 $\tau\to 0$。

### 5. 定义 Kramers–Moyal 系数

定义

$$D^{(n)}(x)=\frac{1}{n!}\lim_{\tau\to 0}\frac{M^{(n)}(x,\tau)}{\tau}.$$

于是得到 Kramers–Moyal 展开：

$$
\frac{\partial P(x,t)}{\partial t}=\sum_{n=1}^{\infty}\left(-\frac{\partial}{\partial x}\right)^n\left[D^{(n)}(x)P(x,t)\right].
$$

### 6. 前两项就是 Fokker–Planck 方程

若只保留前两阶，

$$
\frac{\partial P}{\partial t}=-\frac{\partial}{\partial x}\bigl[D^{(1)}(x)P\bigr]+\frac{\partial^2}{\partial x^2}\bigl[D^{(2)}(x)P\bigr].
$$

记

$$
D^{(1)}(x)=A(x), \qquad D^{(2)}(x)=D(x),
$$

就得到标准 Fokker–Planck 方程：

$$
\frac{\partial P}{\partial t}
=-\frac{\partial}{\partial x}[A(x)P]+\frac{\partial^2}{\partial x^2}[D(x)P].
$$

### 7. 一句话总结

1. Chapman–Kolmogorov 给出短时传播；
2. 对“小跃迁”按跃迁长度 $r$ 做矩展开；
3. 再取连续时间极限，就得到 Kramers–Moyal 级数；
4. 若高阶矩消失，只剩前两项，就得到 Fokker–Planck 方程。

## 3.3 为什么 Fokker–Planck 方程只保留到二阶

若过程来自扩散型随机微分方程

$$
\Delta x = A(x)\tau + B(x)\Delta W,
\qquad \Delta W\sim \mathcal N(0,\tau),
$$

则

$$
\langle \Delta x\rangle = A(x)\tau,
$$

$$
\langle (\Delta x)^2\rangle = B^2(x)\tau + O(\tau^2),
$$

而更高阶矩满足

$$
\langle (\Delta x)^n\rangle = O(\tau^{n/2}), \qquad n\ge 3.
$$

除以 $\tau$ 后令 $\tau\to 0$，高阶项消失，所以只剩：

$$
D^{(1)}(x)=A(x),
$$

$$
D^{(2)}(x)=\frac{B^2(x)}{2}.
$$

于是 Kramers–Moyal 级数截断为：

$$
\frac{\partial P}{\partial t}
=-\frac{\partial}{\partial x}\left[D^{(1)}P\right]+\frac{\partial^2}{\partial x^2}\left[D^{(2)}P\right],
$$

即 Fokker–Planck 方程。

这说明，对于 Itô 扩散过程，Kramers–Moyal 展开天然只到二阶。

而 Chapman–Kolmogorov 方程出发的 Kramers–Moyal 展开，则从更一般的短时转移概率出发，为连续马尔可夫过程提供了统一描述：Fokker–Planck 方程只是该展开在高阶项消失时的二阶截断形式。

---

## 3.4 概率流形式与连续性方程的统一写法

对
$$
\frac{\partial P}{\partial t} = -\frac{\partial}{\partial x}[A(x)P] + \frac{\partial^2}{\partial x^2}[D(x)P]
$$
可写成
$$
\frac{\partial P}{\partial t} = -\frac{\partial J}{\partial x}
$$
其中
$$
J(x,t) = A(x)P(x,t) - \frac{\partial}{\partial x}[D(x)P(x,t)].
$$
- 第一项 $A P$ 是漂移流  
- 第二项 $-\partial_x (D P)$ 是扩散流  
- 稳态时不一定要求 $P$ 不动，而是要求 $\partial_t P = 0\Rightarrow\partial_x J = 0$，即概率流为常数  
- 平衡态常常对应 $J=0$  
- 非平衡稳态则可能是 $J \neq 0$

# 4. 稳态解与平衡分布
## 4.1 稳态条件

在许多物理问题中，我们关心系统在长时间之后会趋于什么状态。

如果概率分布不再随时间变化，即
$$
\frac{\partial P(x,t)}{\partial t}=0
$$
则称系统达到 **稳态（stationary state）**。

Fokker–Planck 方程为

$$
\frac{\partial P}{\partial t}=-\frac{\partial}{\partial x}[A(x)P]+\frac{\partial^2}{\partial x^2}[D(x)P].
$$

代入稳态条件为

$$0=-\frac{\partial}{\partial x}[A(x)P]+\frac{\partial^2}{\partial x^2}[D(x)P].$$

考虑概率流形式：

$$
\frac{\partial P}{\partial t}=-\frac{\partial J}{\partial x}
$$

因此稳态满足

$$
\frac{\partial J}{\partial x}=0
$$

即

$$
J(x)=\text{常数}.
$$

---

## 4.2 平衡态

- **稳态**：分布不再随时间变化
- **平衡态**：不仅分布不变，而且没有净概率流、没有持续循环流动

Fokker–Planck 方程为
$$
\frac{\partial P}{\partial t} = -\frac{\partial}{\partial x}[A(x)P] + \frac{\partial^2}{\partial x^2}[D(x)P]
$$
考虑概率流形式：
$$
\frac{\partial P}{\partial t} = -\frac{\partial J}{\partial x}
$$
得
$$
J(x) = A(x)P(x) - \frac{\partial}{\partial x}[D(x)P(x)].
$$
在许多热平衡系统中，概率流为零：
$$
J(x)=0.
$$
因此
$$
A(x)P(x)-\frac{d}{dx}[D(x)P(x)]=0.
$$
整理得到
$$
\frac{d}{dx}[D(x)P(x)]=A(x)P(x).
$$
若扩散系数为常数
$$
D(x)=D
$$
则
$$
D\frac{dP}{dx}=A(x)P.
$$
于是
$$
\frac{dP}{P}=\frac{A(x)}{D}dx.
$$
积分得到
$$
P(x)\propto
\exp\!\left(
\int^x \frac{A(x')}{D}\,dx'
\right).
$$
只要知道漂移 $A(x)$ 和扩散系数 $D$，就能直接写出零概率流平衡态下的稳态分布 $P(x)$。

---

## 4.3 势场中的布朗粒子

若系统受到势能 $U(x)$ 的作用，

确定性力为
$$
F(x)=-U'(x).
$$
在阻尼介质中
$$
A(x)=\mu F(x)
$$
其中 $\mu$ 为迁移率。

若满足爱因斯坦关系
$$
D=\mu k_B T
$$
则稳态分布为
$$
P(x)\propto e^{-U(x)/k_B T}.
$$
这正是 **Boltzmann 分布**。

这说明：

- Langevin 动力学  
- Fokker–Planck 方程  
- Boltzmann 分布  

三者之间是完全一致的。

或者说，带热噪声的随机动力学自然会把系统带到统计物理的平衡分布。

---
# 5. 总结

### 5.1 随机过程

系统状态随时间演化：
$$
X_t
$$
典型例子：布朗运动。

---
### 5.2 微观轨迹描述（Langevin 方程）

描述 **单条随机轨迹**：
$$
dx=A(x)dt+B(x)dW_t
$$
- $A(x)$：漂移
- $B(x)$：噪声强度

---
### 5.3 概率分布描述（Fokker–Planck 方程）

描述 **概率密度的演化**：

$$\frac{\partial P}{\partial t}=-\partial_x(A P)+\partial_x^2(D P)$$

其中

$$
D=\frac{B^2}{2}.
$$

---
### 5.4 概率守恒结构

Fokker–Planck 方程本质是连续性方程：
$$
\frac{\partial P}{\partial t}=-\partial_x J
$$
概率流
$$
J=A P-\partial_x(DP).
$$

---
### 5.5 Chapman–Kolmogorov → Kramers–Moyal

从马尔可夫过程出发：
$$
P(x,t+\tau)=\int W(x|x',\tau)P(x',t)dx'
$$

做短时跃迁展开得到

$$
\frac{\partial P}{\partial t}=\sum_{n=1}^\infty(-\partial_x)^n[D^{(n)}P].
$$

---
### 5.6 二阶截断 → Fokker–Planck

若过程为扩散型随机过程：
$$
D^{(n)}=0 \quad (n\ge3)
$$

就得到 Fokker–Planck 方程。

---
### 5.7 长时间极限 → 平衡分布

在势场中：
$$
P(x)\propto e^{-U(x)/k_B T}.
$$

---

**一句话：**

> Langevin 方程描述随机轨迹；  
> Fokker–Planck 方程描述概率分布；  
> Chapman–Kolmogorov 与 Kramers–Moyal 提供它们之间的数学桥梁。

---

# 6. 补充与拓展

## 6.1 Kramers 方程与 Smoluchowski 极限

如果粒子具有惯性，需要同时考虑位置与速度：

此时 Langevin 方程写作：

$$
m\frac{dv}{dt}=-\gamma v+F(x)+\sqrt{2\gamma k_B T}\,\xi(t)
$$

以及

$$
\frac{dx}{dt}=v.
$$

对应的概率分布为

$$
P(x,v,t)
$$

其演化满足 **Kramers 方程**。

---

当阻尼很强

$$
m\to0
$$

可以消去速度变量，得到 **Smoluchowski 方程**：

$$
\frac{\partial P(x,t)}{\partial t}=-\partial_x(A P)+\partial_x^2(D P).
$$

形式和**位置空间**中的 Fokker–Planck 方程相同。

---

## 6.2 Fokker–Planck 方程的一般形式

若状态变量记为 $\mathbf{y} = (y_1, \dots, y_n)$，那么概率密度是 $P(\mathbf{y}, t)$

对应的 Fokker–Planck 方程就是

$$
\frac{\partial P(\mathbf{y}, t)}{\partial t} = -\sum_i \frac{\partial}{\partial y_i} \left[ A_i(\mathbf{y}, t) P(\mathbf{y}, t) \right] + \sum_{ij} \frac{\partial^2}{\partial y_i \partial y_j} \left[ D_{ij}(\mathbf{y}, t) P(\mathbf{y}, t) \right].
$$

这里 $\mathbf{y}$ 可以是任何连续状态变量。

- 当状态空间 $\mathbf{y}$ 为 $(x, v)$ 时，$P = P(x, v, t)$，对应 **Kramers 方程**

---

或定义扩散张量

$$
D_{ij}=\frac{1}{2}(BB^\top)_{ij},
$$

则可写成

$$
\frac{\partial P}{\partial t}=-\nabla\cdot(\mathbf A P)+\nabla\nabla:(\mathbf D P).
$$

## 6.3 更一般的随机过程

Fokker–Planck 方程只适用于 **连续扩散过程**。

若跃迁是离散的，例如：

- 化学反应
- 生灭过程
- 跳跃过程

则需要使用 **Master equation**：

$$
\frac{dP_n}{dt}=\sum_m[W_{nm}P_m-W_{mn}P_n].
$$

Kramers–Moyal 展开可以看作 Master equation 的连续极限。

## 6.4 噪声理论

### 6.4.1 什么是噪声

最朴素地说，**噪声就是系统中那些你不显式跟踪、但又会真实影响动力学的随机扰动**。

比如：

- 布朗粒子被大量流体分子不停撞击
- 电路里电子热运动引起电压波动
- 激光、探测器、传感器里存在涨落
- 生物系统里分子数有限，反应有随机性
- 金融模型里价格受大量不可控因素影响

这些微观扰动太复杂，没法一一追踪，于是把它们合并成一个“随机驱动项”，这就是噪声。

### 6.4.2 为什么需要噪声理论

因为现实系统往往不是纯确定性的。  

如果只写  
$$  
\dot{x} = f(x)  
$$
这表示：给定初值，以后轨迹完全确定。  
  
但很多系统更像：  
$$  
\dot{x} = f(x) + \text{噪声}  
$$
也就是：  
- 一部分是确定性动力学  
- 一部分是环境、热浴、未建模自由度造成的随机扰动  

噪声理论研究的就是：  
1. 噪声如何建模  
2. 噪声如何改变轨迹  
3. 噪声如何改变概率分布  
4. 噪声如何导致扩散、跃迁、退相干、线宽展宽、逃逸、随机共振等现象

### 6.4.3 噪声的两个基本视角

**轨迹视角**

研究单个系统随时间怎么走：  
$$dx = a(x, t)\, dt + b(x, t)\, dW_t$$  
这对应 Langevin / 随机微分方程。  

---

**分布视角**

研究很多重复实验的概率密度怎么演化：  
$$\frac{\partial P}{\partial t} = \cdots$$  
这对应 Fokker–Planck / Master equation / Kramers equation 等。  
  
所以噪声理论本质上总在两种语言间切换：  
- 随机轨迹  
- 概率分布

### 6.4.4 分类

#### 1. 按均值分

**零均值噪声**  
$$\langle \xi(t) \rangle = 0$$
最常见。表示没有整体偏置，纯随机扰动。  
  
**非零均值噪声**  
$$\langle \xi(t) \rangle \neq 0$$
通常可拆成“平均部分 + 零均值波动部分”。

---
#### 2. 按时间相关性分

**白噪声**  

相关时间极短，理想化写作  
$$\langle \xi(t) \xi(t') \rangle = 2D \, \delta(t - t')$$不同时刻之间几乎不相关。  

**有色噪声**  

存在非零相关时间，例如  
$$\langle \xi(t) \xi(t') \rangle \sim e^{-|t - t'|/\tau_c}$$  
这里 $\tau_c$ 是相关时间。  

若 $\tau_c \neq 0$，则现在的噪声与稍后的一段时间还有记忆。

---
#### 3. 按分布分

**高斯噪声**

任意时刻或任意线性组合都服从高斯分布。

最重要，因为：

- 中心极限定理支持它常出现
- 数学上最容易处理
- 布朗运动、热噪声常近似成高斯噪声

**非高斯噪声**

如 Lévy 噪声、跳跃噪声、脉冲噪声。

它们可能有厚尾、大跳跃、非平滑行为。

---
#### 4. 按进入方程的方式分

**加性噪声**  
$$\dot{x} = f(x) + \xi(t)$$
噪声强度与状态无关。  

**乘性噪声**  
$$\dot{x} = f(x) + g(x)\xi(t)$$
噪声强度依赖状态。

---
#### 5. 按是否有记忆分

**马尔可夫噪声**

未来只和当前有关，不依赖更久远历史。

**非马尔可夫噪声**

系统存在记忆核、历史依赖。  
这时简单 Fokker–Planck 往往不够，需要更一般方法。

## 6.5 Itô 和 Stratonovich

先看随机微分方程：  
$$ dx = a(x, t)\, dt + b(x, t)\, dW_t $$
这里 $dW_t$ 是 Wiener 过程增量。  

问题在于：  

因为 $W_t$ 的路径不可导、非常粗糙，普通微积分里的积分定义不再直接适用。  
所以像  
$$ \int b(x, t)\, dW_t $$  
这种积分，必须说明你怎么取采样点。

---

**Itô 积分**

Itô 解释里，在每个小区间 $[t_n, t_{n+1}]$ 上，系数取左端点：  
$$  
\int_0^t b(X_s, s) \, dW_s \approx \sum_n b(X_{t_n}, t_n) \left( W_{t_{n+1}} - W_{t_n} \right)  
$$
也就是“先看当前状态，再乘上这一步噪声”。  

这在概率论里特别自然，因为它是非预见的：  
  
- 你在时刻 $t_n$ 只能用当前信息  
- 不能偷看未来的噪声增量  

所以 Itô 很适合金融、马尔可夫过程、严格概率论。

---

**Stratonovich 积分**

Stratonovich 解释里，系数取中点（更准确说是对称取样）：  
$$  
\int_0^t b(X_s, s) \circ dW_s \approx \sum_n b\left(\frac{X_{t_n} + X_{t_{n+1}}}{2}, \frac{t_n + t_{n+1}}{2}\right)(W_{t_{n+1}} - W_{t_n})  
$$
常用符号是小圆圈：  
$$  
dx = a(x, t)\, dt + b(x, t) \circ dW_t  
$$
这里的 $\circ\  dW_t$ 就表示 Stratonovich。  
  
它更像普通微积分，保留通常的链式法则，所以在物理里很常见，尤其在：  
  
- 噪声有很短但非零相关时间  
- 白噪声是某个光滑有色噪声极限  
  
这种情况下，Stratonovich 常更自然。

---
## 6.5.1 为什么会有区别

核心原因只有一句：Wiener 过程太粗糙，$(dW)^2 \sim dt$ 不能忽略。  

普通微积分里，小量平方往往比一阶小得多，可以丢掉。  

但随机微积分里：  
$$dW \sim \sqrt{dt}$$
所以  
$$(dW)^2 \sim dt$$
这和 $dt$ 是同阶的，不能随便扔。  

于是你在展开 $b(x)$ 时，采样点选左端还是中点，就会在平均漂移里留下额外项。

---

**什么时候它们没区别**

如果噪声强度 $b$ 与状态无关，即  
$$  
dx = a(x)\,dt + \sigma\,dW_t  
$$

这里 $\sigma$ 是常数，那么 Itô 和 Stratonovich 给出同一个结果。  

因为这时 $b'(x) = 0$，不会出现额外漂移修正。  

这叫**加性噪声**。  

---

**什么时候它们有区别**

如果噪声强度依赖于 $x$：  
$$  
dx = a(x)\,dt + b(x)\,dW_t  
$$  
这叫**乘性噪声**。  

这时 $b(x)$ 在不同位置大小不同，噪声本身“带结构”，于是不同积分解释会改变有效漂移。

## 6.5.2 Itô 和 Stratonovich 的转换公式

这是最关健的公式。  

如果 Stratonovich 方程写成  
$$  
dx = a_S(x, t)\, dt + b(x, t) \circ dW_t  
$$
那么与它等价的 Itô 形式是  
$$  
dx = a_I(x, t)\, dt + b(x, t)\, dW_t  
$$
其中  
$$  
a_I(x, t) = a_S(x, t) + \frac{1}{2} b(x, t)\, \frac{\partial b(x, t)}{\partial x}  
$$
一维时就是这样。  
  
也就是说，Stratonovich 漂移转成 Itô 漂移，要多出一个“噪声诱导漂移” $\frac{1}{2}\, b\, b'$

---

**为什么会多出这个项**

直觉上：

如果 $b(x)$ 在不同位置不一样大，那么粒子在噪声强的地方会被“抖”得更厉害。

这种不对称抖动，平均起来就会产生一个净漂移。

所以虽然随机项均值是 0，  

但“噪声强度随位置变化” + “随机路径粗糙”会导致额外平均效应。

这就是为什么 Stratonovich 和 Itô 不再只是“写法不同”，而是漂移项真的不同。

---

**它为什么影响 Fokker–Planck 形式**  

因为 Fokker–Planck 方程中的漂移项，针对的是 Itô 形式的漂移系数。  
  
对于一维 Itô SDE：  
$$  
dx = a_I(x)\,dt + b(x)\,dW_t  
$$  
对应 Fokker–Planck 方程是  
$$  
\frac{\partial P}{\partial t} = -\frac{\partial}{\partial x}\big[a_I(x)P\big] + \frac{1}{2}\frac{\partial^2}{\partial x^2}\big[b^2(x)P\big]  
$$
如果你原来写的是 Stratonovich 形式：  
$$  
dx = a_S(x)\,dt + b(x)\circ dW_t  
$$
那就要先转成 Itô：  
$$  
a_I(x) = a_S(x) + \frac{1}{2}b(x)b'(x)  
$$
代进去得到  
$$  
\frac{\partial P}{\partial t} = -\frac{\partial}{\partial x}\left[\left(a_S + \frac{1}{2}bb'\right)P\right] + \frac{1}{2}\frac{\partial^2}{\partial x^2}[b^2 P]  
$$
所以你看见没有：

- 扩散项仍由 $b^2$ 决定
- 但漂移项会因为 Itô/Stratonovich 解释不同而不同

这就是“为什么它会影响 Fokker–Planck 形式”。

---

**一个很简单的例子**  

看这个方程：  
$$dx = x \circ dW_t$$
这是 Stratonovich 形式，没有显式漂移项，即  
$$a_S(x) = 0, \quad b(x) = x$$
转成 Itô：  
$$a_I(x) = 0 + \frac{1}{2}x \cdot \frac{d}{dx}x = \frac{1}{2}x$$
所以等价的 Itô 方程是  
$$dx = \frac{1}{2}x\,dt + x\,dW_t$$
你会发现：  

- Stratonovich 写起来“只有噪声”  
- Itô 写起来却“多了一个确定性漂移 $\frac{1}{2}x$”  

于是对应的 Fokker–Planck 方程不是纯扩散，而是  
$$\frac{\partial P}{\partial t} = -\frac{\partial}{\partial x}\left(\frac{1}{2}xP\right) + \frac{1}{2}\frac{\partial^2}{\partial x^2}(x^2 P)$$
所以“同一个物理过程”的不同解释，本质上不是矛盾，而是：

> **你必须先说明随机积分的定义，才能正确写出对应的 Fokker–Planck 方程。**

---

**它们在链式法则上也不同**

**Stratonovich 保留普通链式法则**  
  
如果  
$$dx = a_S\,dt + b \circ dW$$
对一个函数 $f(x)$，Stratonovich 下你几乎可以像普通微积分一样写：  
$$df = f'(x)\,dx$$
当然展开后再代入 $dx$。  
  
这也是很多物理学家喜欢它的原因：形式上更接近经典微积分。

**Itô 要用 Itô 引理**
  
Itô 下不能直接用普通链式法则，必须加二阶项：  
$$
df(x,t) = \left( \frac{\partial f}{\partial t} + a_I \frac{\partial f}{\partial x} + \frac{1}{2} b^2 \frac{\partial^2 f}{\partial x^2} \right) dt + b \frac{\partial f}{\partial x} dW_t  
$$
这个额外的  
$$
\frac{1}{2} b^2 f_{xx}  
$$
就是因为 $(dW)^2 \sim dt$。 所以： 
- Stratonovich：链式法则像普通微积分 
- Itô：链式法则带额外二阶修正

## 6.5.3 小结

**物理上该选哪个**

这不是“谁对谁错”，而是**模型来源不同，解释不同**。

**Itô 更自然于：**

- 纯数学概率论
- 马尔可夫过程生成元
- 金融建模
- 直接写 Fokker–Planck 对应关系时

**Stratonovich 更自然于：**

- 物理中的白噪声极限
- 把真实短相关时间噪声压缩成白噪声时
- 希望保留普通微积分形式时

很多物理系统从“有色噪声极限”推到白噪声时，最后自然出现 Stratonovich。

---

最稳妥的记法是：  
  
**如果 SDE 是 Itô 形式**  
$$  
dx = a_I(x)\,dt + b(x)\,dW_t  
$$
那么 Fokker–Planck 直接是  
$$  
\frac{\partial P}{\partial t} = -\partial_x (a_I P) + \frac{1}{2}\partial_x^2 (b^2 P)  
$$

**如果 SDE 是 Stratonovich 形式**  
$$  
dx = a_S(x)\,dt + b(x) \circ dW_t  
$$
先转成 Itô：  
$$  
a_I = a_S + \frac{1}{2}bb'  
$$
再写 Fokker–Planck：  
$$  
\frac{\partial P}{\partial t} = -\partial_x \left[ \left( a_S + \frac{1}{2}bb' \right) P \right] + \frac{1}{2}\partial_x^2 (b^2 P)  
$$
所以真正受影响的是漂移项的定义。


