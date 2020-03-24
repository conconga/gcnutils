# System of Systems of ODE's

The object cSOSODE() stands for a system of several systems of ordinary differential equations. The object concatenates several systems of ODE's in a unique callable function that can be solved by any numerical integration function or object. The identification of the variables of each equation/system is symbolic, and the object cSOSODE() builds the system based on the symbols.

When a system depends on parameters calculated by different functions, or a cascade of functions, the object determines the sequence of calls to solve the parameters.

# Each System

Each system is created by the object cSOSODE_FUNCTION(), and defined by the input states and parameters, and output d(states)/dt or parameters.

# Examples

The namespace __ main __ in the code provides several standalone examples of usage of the object.

