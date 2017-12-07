# Function output convention for exact rational objects

For speed, we wrote all our routines using the default `double` MATLAB floating-point type. In some cases, however, computations should be done in exact arithmetic to avoid numerical problems.

Most of the functions/methods in our code allow the use of exact arithmetic with the following trick.

When a function/method is called with a single output parameter, a floating-point approximation is returned. However, when two output parameters are required, the method provides a (numerator, denominator) pair, and the exact solution is represented by numerator/denominator.

Example:

```matlab

>> S = qvx.di.Scenario.homogenous(1,2,2)

S = 

{[2 2]}
>> S.conversionMatrix('SG', 'P')

ans =

    0.5000    0.5000    0.5000    0.5000
    0.7500   -0.2500    0.2500    0.2500
    0.2500    0.2500    0.7500   -0.2500
    1.0000    1.0000   -1.0000   -1.0000

>> [num den] = S.conversionMatrix('SG', 'P')

num =

     2     2     2     2
     3    -1     1     1
     1     1     3    -1
     4     4    -4    -4


den =

     4
```
