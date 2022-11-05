# Particle Swarm Optimization (Improved)

## 1 Preparation

**Note:**	Variables are expressed as line vectors: $x = [f, \phi]^\top$.

**Step1**	Judge if input bounds are line vectors

**Step2**	Merge user options with default options with function "SetOptions"

Including: Inertia range, Self adjustment factor, Social adjustment factor, Swarm size, Minimum neighborhood size fraction, Maximum generation, Initial swarm, Initial swarm span, Parameter necessary for objective function, Display option

**Step 3**	Assign some parameters according to options

## 2	Initialization

**Step1**	Create initial state

Including: iteration, startTime, stopFlag, lastImprovement, lastImprovementTime, funEval, velocities, fvals, individualBestFvals, individualBestPositions

**Step 2**	Set adaptive inertia and its counter

**Step 3**	Set display header

Including: iteration, f-count, best f(x), mean f(x), stall iterations

## 3	Iteration

**Step 1**	Update iteration

**Step 2**	Generate best neighbor's index

**Step 3**	Update particle velocities

**Step 4**	Update particle positions

**Step 5**	Update particle velocities on the boundary

**Step 6**	Update current objective function value

**Step 6**	Update optimization state

Including: funEval, individualBestFvals, individualBestPositions

**Step 7**	Update adaptive inertia