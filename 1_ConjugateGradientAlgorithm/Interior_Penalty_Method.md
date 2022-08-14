# Interior Penalty Method

## 1	Problem Description

$$
\begin{align}
&\min \quad f(x) \\
&s.t. \quad g_i(x) \le 0, \quad i = 1,2,\dots,m
\end{align}
$$

For the joint estimation problem, the inequation constraints are bound constraints:
$$
f_L \le f \le f_U \\
p_L \le p \le p_U
$$
which is equal to four inequations:
$$
\left\{
\begin{align}
f - f_U \le 0 \\
f_L - f \le 0 \\
p - p_U \le 0 \\
p_L - p \le 0
\end{align}
\right.
$$

## 2	Theory of The Algorithm

First, define a penalty function:
$$
B(x) = \sum\limits_{i=1}^L g_i^+(x)
$$
where
$$
g_i^+(x) = -\mathrm{ln}(-g_i(x)) \quad or \quad g_i^+(x) = -\frac{1}{g_i(x)}
$$
On this basis, a enhanced object function is defined:
$$
F(x,r_k) = f(x) + r_k B(x)
$$
where $r_k > 0$ and $r_1>r_2>\cdots>r_k>\cdots$, $\lim\limits_{k \to \infin} r_k = 0$.

## 3	Algorithm Implementation

Step 1:	Select initial penalty factor $r_1$ and make sure $r_1 > 0$. Set search accuracy $\epsilon > 0$ and penalty factors iteration factor $c \ge 2$ (normally $c \in [4,10]$).

Step 2:	Find an interior point in the feasible domain $S$ as the starting point $x_0$ and set iteration time $k=1$.

Step 3:	Iteration loop

(1)	Solve an equivalent unconstrained optimization problem
$$
\begin{align}
&\min \quad F(x,r_k) = f(x) + r_k B(x) \\
&s.t. \quad x \in \mathrm{int} S
\end{align}
$$
The solution is $x_k = x(r_k)$.

(2) If $x_k$ satisfies the end criterion, the iteration ends. Otherwise, let $r_{k+1} = r_k / c$ and $k = k+1$, the iteration continues.

## 4	End Criterion

$$
r_k B(x_k) < \epsilon
$$

This can be interpreted as the boundary constraint is small and has limited effect on the optimal solution.
$$
\Vert x_k - x_{k-1} \Vert < \epsilon
$$
