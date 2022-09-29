# Particle Swarm Optimization

## 1	Parameters

(1)	$D$ :	Searching dimension (number of variables)

(2)	$P$ :	Number of populations

(3)	$G$ :	Maximum number of generations

(4)	$x_{ij}^t$ :	Position of $i$th particle of $j$th population in $t$th generation

(5)	$v_{ij}^t$ :	Velocity of $i$th particle of $j$th population in $t$th generation



## 2	Optimization Process

### 2.1	Initialization

#### 2.1.1	Initialize particle position

Randomly generate two set of particles $X1,X2$ (Initial guess) $\to$ Move initial points towards initial guess, by convex combination $\to$ Initialize particle position and velocity
$$
X^0 = (\vec{x}_1^0,\ \vec{x}_2^0,\ \dots\ ,\ \vec{x}_P^0) \\
V^0 = (\vec{v}_1^0,\ \vec{v}_2^0,\ \dots\ ,\ \vec{v}_P^0)
$$
where
$$
\vec{x}_j^t = (x_{1j}^t,\ x_{2j}^t,\ \dots\ ,\ x_{Dj}^t)^{\mathrm{T}} \\
\vec{v}_j^t = (v_{1j}^t,\ v_{2j}^t,\ \dots\ ,\ v_{Dj}^t)^{\mathrm{T}}
$$
If defined "warm start", particles of the first population are just the initial particles input.

#### 2.1.2	Initialize fitness function

$$
F_{1 \times P}^0 = f(X)
$$

#### 2.1.3	Initialize best points

$$
P_{best}^0 = X^0 \quad (D \times P) \\
F_{best}^0 = F^0 \quad (1 \times P) \\
F_{global}^0 = \mathop{max}_j\{F_{best,i}\} \\
G_{best}^0 = \mathop{max}_{\vec{x}_j^0}\{F_j^0\} \quad (D \times 1)
$$

#### 2.1.4	Allocate memory for info and data log

### 2.2	Iteration

#### 2.2.1	Generate random factors $r_{1_{ij}}^{t+1}$ and $r_{2_{ij}}^{t+1}$

$r_{1_{ij}}^{t+1}$ and $r_{2_{ij}}^{t+1}$ are two independent random numbers at $(t+1)$th generation, which are uniformly distributed within the range $[0,1]$, and sampled independently for every particle $j$ at each dimension $i$.

#### 2.2.2	Update particle velocity

$$
v_{ij}^{t+1} = \omega v_{ij}^t + c_1 r_{1_{ij}}^{t+1}(p_{best_{ij}}^t-x_{ij}^t) + c_2 r_{2_{ij}}^{t+1}(g_{best_i}^t-x_{ij}^t)
$$

where $\omega$ is inertia coefficient,  $c_1$ is cognitive coefficient, $c_2$ is social coefficient

#### 2.2.3	Update particle position

$$
x_{ij}^{t+1}=x_{ij}^t+v_{ij}^{t+1}
$$

Check if new position is in the bound after that.

#### 2.2.4	Update key indicators

(1)	Update fitness value of all particles $F$

(2)	Update best fitness value of each population $F_{best}$

(3)	Update best particle of each population $P_{best}$

(4)	Update best fitness value of all particles $F_{global}$

(5)	Update particle with best fitness value of all populations $G_{best}$