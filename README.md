-   [<span class="toc-section-number">1</span> System of Systems of
    ODE’s](#system-of-systems-of-odes)
-   [<span class="toc-section-number">2</span> Each piece of the
    System](#each-piece-of-the-system)
-   [<span class="toc-section-number">3</span> Examples](#examples)
    -   [<span class="toc-section-number">3.1</span> Example
        1](#example-1)
    -   [<span class="toc-section-number">3.2</span> Example 2 (with
        integration)](#example-2-with-integration)
    -   [<span class="toc-section-number">3.3</span> Example 3 (using
        integration object)](#example-3-using-integration-object)

# System of Systems of ODE's

The object `kSosode()` stands for a system of several systems of ordinary
differential equations. The object concatenates several systems of ODE's in a
unique callable function that can be solved by any numerical integration
function or object. The identification of the variables of each equation/system
is symbolic, and the object `kSosode()` builds the system based on the symbols.

When a system depends on parameters calculated by different functions, or a
cascade of functions, the object determines the sequence of calls to solve the
parameters.

The goal is to generate a function as:

$$ dy(t)/dt = f(y,t) $$

where $f(y,t)$ is callable by the integrator, and $y$ is the global state
$[y_1, \dots, y_L] $.  Internally, hidden to the Integrator, $f(y,t)$ can be a
complicated mess of ODE's and parameters.

The constituent functions to generate $f(y,t)$ are $f_i(y_i,t,P_i)$, where
$y_i$ is a subset of $y$, and $P_i$ is a subset of all available parameters.

The function $f_i$ is registered with four descriptors:

| descriptor | description |
| ---------- | ----------- |
| `i_state` | list with labels of input states of $f_i$ |
| `i_param` | list with labels of input params of $f_i$ |
| `o_state` | list with labels of y states by which derivatives are calculated by $f_i$ |
| `o_param` | list with labels of output params calculated by $f_i$ |

The order of parameters and returns is relevant when creating the lists. For
example, the function `[a,b] = f([c,d],t,[e,f])` needs `i_state=['c','d']` with
this order, as much as `i_param=['e','f']`. The return with also be
`o_state=['a','b']`.

A sanity check is provided to certify that all necessary inputs and params are
available to calculate the global state $y$.

# Each piece of the System

Each piece of the system is created by the object `kSosodeFunction()`, and
defined by the input states and parameters, and output d(states)/dt or
parameters.

# Examples

## Example 1

Consider the following differential system of equations:

$$ \[ \dot{y}_1, \dot{y}_3 \]^T = f_1(y,t,\theta_1,\theta_2) $$

and

$$ \[ \dot{y}_2, \dot{y}_4 \]^T = f_2(y,t) $$

with:

$$ \[ \theta_1, \theta_2 \]^T = g(y,t) $$

The single system of systems is build with:
```
    fn1 = kSosodeFunction(f1)
    fn1.set_i_state([ "y 1", "y 2", "y 3", "y 4" ])
    fn1.set_i_param([ "theta 1", "theta 2" ])
    fn1.set_o_state([ "y 1", "y 3" ])

    fn2 = kSosodeFunction(f2)
    fn2.set_i_state([ "y 1", "y 2", "y 3", "y 4" ])
    fn2.set_o_state([ "y 2", "y 4" ])

    fn3 = kSosodeFunction(g1)
    fn3.set_o_param([ "theta 1", "theta 2" ])

    # : reverse      : True if input order will be (y,t) instead of (t,y)
    # : order_states : list with the order of states in the state vector.

    a = kSosode( fn1, fn2, fn3, reverse=True, order_states= [ "y 2", "y 3", "y 4", "y 1" ] )
    a.create_nets()
    a.showregisteredfunctions()
```

## Example 2 (with integration)

```
    def tta(t):
        return 8 if t < 10 else -2

    def eq0(t, y, *args):
        tta = args[0]
        return (-8.*y[0]) + (tta*y[1])

    def eq1(t, y):
        return (30*y[0]) + y[1] - (y[0]*y[2])

    def eq2(t, y):
        return (y[0]*y[1]) - (8.*y[2]/3)

    fn0 = kSosodeFunction(eq0)
    fn0.set_i_state([ 'y0', 'y1' ])
    fn0.set_i_param([ 'tta' ])
    fn0.set_o_state([ 'y0' ])

    fn1 = kSosodeFunction(eq1)
    fn1.set_i_state([ 'y0', 'y1', 'y2' ])
    fn1.set_o_state([ 'y1' ])

    fn2 = kSosodeFunction(eq2)
    fn2.set_i_state([ 'y0', 'y1', 'y2' ])
    fn2.set_o_state([ 'y2' ])

    g0 = kSosodeFunction(tta)
    g0.set_o_param([ 'tta' ])

    # : reverse      : True if input order will be (y,t) instead of (t,y)
    # : order_states : list with the order of states in the state vector.

    b = kSosode( fn0, fn1, fn2, g0, reverse=True, order_states= [ 'y0', 'y1', 'y2' ] )
    b.create_nets()
    b.showregisteredfunctions()

    T = np.linspace(0,12.0,1000)
    R = Int.odeint(b, [-1, 0, 1], T, args=() )
```

## Example 3 (using integration object)

(See example in `__main__` in source files.)

```
class kExample_RC_discharge(kSosodeUtils, kSosodeIntegrator, kExample_RC_discharge_system):
    def __init__(self, sample_freq_Hz, **kargs):
        # initialize the System of Systems model:
        super().__init__(**kargs)

        assert sample_freq_Hz is not None
        assert sample_freq_Hz > 0

        self.dt        = 1./sample_freq_Hz

    def get_V(self):
        return self.pick_from_state( 'V' )

for rc in [kExample_RC_discharge(100, R=1e6), kExample_RC_discharge(100, V=3, R=100e3, C=4.7e-6)]:
    for t in T:
        # update current state vector:
        rc.update()
        print(rc.get_V())
```
