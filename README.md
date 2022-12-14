# D12 Rice Gray

Calculates the tracer diffusivity of a given solute in a solvent. If used please cite: "Zêzere, B.; Portugal, I.; Gomes, J.R.B.; Silva, C.M. Modeling Tracer Diffusion Coefficients of Any Type of Solutes in Polar and Non-Polar Dense Solvents. Materials 2022, 15, 6416. [https://doi.org/10.3390/ma1518641](https://doi.org/10.3390/ma15186416)]"

## Requiered Data:

* `T` - Absolute temperature in K in column format;
* `Density` - Density of the solvent in g/cucm  in column format;
* `M` - Molecular mass of the compounds by the following order solute, solvent;
* `Tc` - Critical temperature in K of the compounds by the following order solute, solvent;
* `Vc` - Critical volume in cucm/gmol of the compounds by the following order solute, solvent;
* `k12` - Binary binary interaction parameter introduced in the LJ diameter combination rule;
* `B12` - Binary binary interaction parameter.

`M`, `Tc`, `Vc`, `k12` and `B12` for known systems can be found in: doi: [https://doi.org/10.3390/ma1518641](https://doi.org/10.3390/ma15186416)

## Code examples

Copy paste on Matlab to run, tested in Matlab 2021b

### Example 1 CO2/ehanol 1 data point
```matlab
T = 313.21; %K

Density = 0.74364; % g/cum

M = [46.069 44.01];
Tc = [513.9 304.1];
Vc = [167.1 93.9];
k12 = 0.240650;
B12 = 1.4688970;

D12calc_P2 = Rice_and_Gray(T, Density, M, Tc, Vc, k12, B12)
```


### Example 2 CO2/ibuprofen Multiple Data Points
```matlab
T = [313.15; 313.15; 318.15; 318.15]; % K

Density = [0.93481; 0.95607; 0.78462; 0.78742]; % g/cum

M = [206.29 44.01];
Tc = [769.6305 304.1];
Vc = [686.35 93.9];
k12 = 0.19035;
B12 = 0.700406;

D12calc_P2 = Rice_and_Gray(T, Density, M, Tc, Vc, k12, B12)
```
