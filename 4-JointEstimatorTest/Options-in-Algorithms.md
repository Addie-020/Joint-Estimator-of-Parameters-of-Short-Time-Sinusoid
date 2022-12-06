# Options in Algorithms

## 1 Options in Particle Swarm Optimization

|    Option Name     |         Default Value          |
| :----------------: | :----------------------------: |
|    inertiaRange    |          $[0.1,1.1]$           |
|   selfAdjustment   |              1.49              |
|  socialAdjustment  |              1.49              |
|     swarmSize      | $\min(100,10*\mathrm{nvars})$  |
|    initialSwarm    |             $[\ ]$             |
|  initialSwarmSpan  |              2000              |
|   maxGeneration    | $\min(100,200*\mathrm{nvars})$ |
|      maxTime       |           $+\infin$            |
|  minNeighborFrac   |              0.25              |
|   objectiveLimit   |           $-\infin$            |
|    tolFunValue     |           $10^{-9}$            |
|   stallIterLimit   |               20               |
|   stallTimeLimit   |           $+\infin$            |
|      display       |             'none'             |
|  displayInterval   |               1                |
| displaySectionSize |               20               |

## 2 Options in Gradient Optimization

1. **Algorithm**: 'interior-point', 'trust-region-reflective', 'sqp', 'sqp-legacy', 'active-set'

2. **CheckGradients**:  Compare user-supplied derivatives (gradients of objective or constraints) to finite-differencing derivatives

3. **ConstraintTolerance**: FunctionTolerance & StepTolerance

4. FiniteDifferenceType:  Finite differences, used to estimate gradients, are either `'forward'` (default), or `'central'` (centered). `'central'` takes twice as many function evaluations but should be more accurate. The trust-region-reflective algorithm uses `FiniteDifferenceType` only when `CheckGradients` is set to `true`.

   `fmincon` is careful to obey bounds when estimating both types of finite differences. So, for example, it could take a backward, rather than a forward, difference to avoid evaluating at a point outside bounds. However, for the `interior-point` algorithm, `'central'` differences might violate bounds during their evaluation if the `HonorBounds` option is set to `false`.

5. FiniteDifferenceStepSize: The default is `sqrt(eps)` for forward finite differences, and `eps^(1/3)` for central finite differences.

|      Option Name       |        Default Value        |
| :--------------------: | :-------------------------: |
|     DiffMaxChange      |          $+\infin$          |
|     DiffMinChange      |              0              |
|        Display         |           'final'           |
|      FunValCheck       |             off             |
| MaxFunctionEvaluations |            3000             |
|     MaxIterations      |            1000             |
|     StepTolerance      |         $10^{-10}$          |
|        TypicalX        | `ones(numberofvariables,1)` |





















