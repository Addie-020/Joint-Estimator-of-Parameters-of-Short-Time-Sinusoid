# Instruction about Correlation Computation

## 1	Definition

The Pearson correlation coefficient is defined as
$$
\rho(A, B) = \frac{1}{N-1} \sum\limits_{i=1}^N \left( \frac{A_i-\mu_A}{\sigma_A} \right) \left( \frac{B_i-\mu_B}{\sigma_B} \right)
$$
where $N$ is the number of signal samples, $\mu_A,\mu_B$ is the signal's mean and $\sigma_A,\sigma_B$ is the signal's standard variance.

## 2	Computation with Matrix

Since we have already known all the information of the given signal $x(n)$ before computation, we can do some necessary computations in advance.

Let $\boldsymbol{Ct} = \left( \frac{A_1-\mu_A}{\sigma_A}, \frac{A_2-\mu_A}{\sigma_A}, \dots, \frac{A_N-\mu_A}{\sigma_A} \right)$ and $\boldsymbol{C_e} = \left( \frac{B_1-\mu_B}{\sigma_B}, \frac{B_2-\mu_B}{\sigma_B}, \dots, \frac{B_N-\mu_B}{\sigma_B} \right)$. Then $\rho(A, B) = \frac{1}{N-1} \boldsymbol{Ct} \boldsymbol{C_e}^\top$ and $\boldsymbol{Ct}$ can be computed before optimization process to reduce computation complex in the next steps.

### 2.1	Compute $\boldsymbol{Ct}$

We can easily calculate the exact value of $\boldsymbol{Ct}$, the process is clearly presented in the code.

### 2.2	Compute $\boldsymbol{C_e}$

Since the object function should be adapt to vectorized inputs and outputs, the dimension of matrix during the computation process including $\boldsymbol{C_e}$ is different from $\boldsymbol{Ct}$.

Assume that the dimension of variable is $D$ and the number of input variable vectors is $R$. That is to say, the object function program have to deal with $R$ different sequences at a time.

To show the calculation process more clearly, we can start from looking at the result. The output of the calculation should be $\boldsymbol{Y}_{1 \times R} = (y_1, y_2,\dots, y_R)$ and each correlation coefficient $y_i = \frac{1}{N-1} \boldsymbol{Ct} \boldsymbol{C_{ei}}^\top$. Thus, the output vector can be expressed as
$$
\begin{align}
\boldsymbol{Y} &= \frac{1}{N-1} \boldsymbol{Ct} (\boldsymbol{C_{e1}}^\top, \boldsymbol{C_{e2}}^\top, \dots, \boldsymbol{C_{eR}}^\top) \\ \\
&= \frac{1}{N-1} \boldsymbol{C_t} \boldsymbol{C_e}^\top
\end{align}
$$
where $\boldsymbol{C_e} = (\boldsymbol{C_{e1}}; \boldsymbol{C_{e2}}; \dots; \boldsymbol{C_{eR}})$.

On this basis, we then consider calculate $\boldsymbol{C_{2i}}$ and finally get the value of $\boldsymbol{C_e}$.
$$
\begin{align}
\boldsymbol{C_{2i}} &= \left( \frac{B_{i1}-\mu_{B_i}}{\sigma_{B_i}}, \frac{B_{i2}-\mu_{B_i}}{\sigma_{B_i}}, \dots, \frac{B_{iN}-\mu_{B_i}}{\sigma_{B_i}} \right) \\ \\
&= \frac{1}{\sigma_{B_i}} \left[ (B_{i1}, B_{i2}, \dots, B_{iN}) - (\mu_{B_i}, \mu_{B_i}, \dots, \mu_{B_i}) \right] \\ \\
& = \frac{1}{\sigma_{B_i}} (\boldsymbol{B_i} - \boldsymbol{\Mu_{B_i}})
\end{align}
$$
where $\boldsymbol{B_i} = a_i\mathrm{sin}(2\pi f_i \boldsymbol{X_t} + p_i)$ and $\boldsymbol{X_t} = \frac{1}{Fs} (0, 1, \dots, N)$.

To calculate $\boldsymbol{C_e}$, we combine $\boldsymbol{C_{2i}}$s together.
$$
\boldsymbol{C_e} = 
\begin{bmatrix}
\frac{1}{\sigma_{B_1}} (\boldsymbol{B_1} - \boldsymbol{\Mu_{B_1}}) \\
\frac{1}{\sigma_{B_2}} (\boldsymbol{B_2} - \boldsymbol{\Mu_{B_2}}) \\
\cdots \\
\frac{1}{\sigma_{B_R}} (\boldsymbol{B_R} - \boldsymbol{\Mu_{B_R}})
\end{bmatrix} = 
\begin{bmatrix}
\boldsymbol{B_1} - \boldsymbol{\Mu_{B_1}} \\
\boldsymbol{B_2} - \boldsymbol{\Mu_{B_2}} \\
\cdots \\
\boldsymbol{B_R} - \boldsymbol{\Mu_{B_R}}
\end{bmatrix} ./
\begin{bmatrix}
\sigma_{B_1} & \sigma_{B_1} & \cdots & \sigma_{B_1} \\
\sigma_{B_2} & \sigma_{B_2} & \cdots & \sigma_{B_2} \\
\vdots & \vdots & \ddots & \vdots \\
\sigma_{B_R} & \sigma_{B_R} & \cdots & \sigma_{B_R}
\end{bmatrix} = 
(\boldsymbol{B} - \boldsymbol{\Mu_B}) ./ \boldsymbol{\Sigma_B}
$$
Then we look into the value of $\boldsymbol{B}$.
$$
\begin{align}
\boldsymbol{B} &=
\begin{bmatrix}
a_1\mathrm{sin}(2\pi f_1 \boldsymbol{X_t} + p_1) \\
a_2\mathrm{sin}(2\pi f_2 \boldsymbol{X_t} + p_2) \\
\cdots \\
a_R\mathrm{sin}(2\pi f_R \boldsymbol{X_t} + p_R)
\end{bmatrix} =
\begin{bmatrix}
a_1 & a_1 & \cdots & a_1 \\ a_2 & a_2 & \cdots & a_2 \\
\vdots & \vdots & \ddots & \vdots \\ a_R & a_R & \cdots & a_R
\end{bmatrix}\ .*\ \mathrm{sin} \left(2\pi
\begin{bmatrix}
f_1 \\ f_2 \\ \cdots \\ f_R
\end{bmatrix}
\begin{bmatrix}
x_1 & x_2 & \cdots & x_N
\end{bmatrix} +
\begin{bmatrix}
p_1 & p_1 & \cdots & p_1 \\ p_2 & p_2 & \cdots & p_2 \\
\vdots & \vdots & \ddots & \vdots \\ p_R & p_R & \cdots & p_R
\end{bmatrix} \right) \\ \\ &= \boldsymbol{A_{rep}} \odot \mathrm{sin} \left( 2\pi \boldsymbol{F} \boldsymbol{X_t} + \boldsymbol{P_{rep}} \right)
\end{align}
$$
$\boldsymbol{\Mu_B}, \boldsymbol{\Sigma_B}$ then can be calculated easily.