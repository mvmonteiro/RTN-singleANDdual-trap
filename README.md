# Random telegraph noise trap simulation

This repository has the aim of build a struct that can simulate a single and a dual trap behavior. 

## Context

The downscaling of channel area of transistor MOS leads to problems that impact in the clarity of the signal inside the circuits. So we can
increase the reability on this kind of system by some simulation tests with charge trapping.

![flux](https://user-images.githubusercontent.com/48728730/217322775-eeea2e74-48c9-4970-83b0-d1a5ce2c9dec.png)

The target of this simulations is based on the RTN behavior with the trap activity probabilistic calculations following this kind of struct of noise.

## Flux

The code follows the code flux above for each simulator:

![flux2](https://user-images.githubusercontent.com/48728730/217324507-dc9c5c2c-b6b3-4bca-a365-36ddb62c5472.png)

The difference between them both is in the "rtn_calc" function. When using the dual trap simulator, this function has two differents traps 
working simultaneously, enabling the correlated activity of them or not, it depends on the probabilistic of the emission time and capture time (tc[0] and tc[1] respectivily) of the traps having some relation.

## Execution

Single trap: gcc rtn_simple.c -o rtn -lm -std=c99 -lfftw3l
Dual trap: gcc rtn_dual.c -o rtn -lm -std=c99 -lfftw3l

**For both simulations you will need the FFT libery for C. You can find all you need about this here: https://www.fftw.org**

## Tools

- MATLAB
- C
- Python
- ShellScript
