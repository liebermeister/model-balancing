Model balancing
===============

[Model balancing](https://www.metabolic-economics.de/model-balancing/index.html) is a computational method to determine 
consistent states and kinetic constants of metabolic models.

## Dependencies
### Matlab
- [SBML toolbox](http://sbml.org/Software/SBMLToolbox)
- Clone the following [GitHub](https://github.com/liebermeister) repositories
    - [`matlab-utils`](https://github.com/liebermeister/matlab-utils) - utility functions
    - [`metabolic-network-toolbox`](https://github.com/liebermeister/metabolic-network-toolbox) - metabolic network toolbox
    - [`sbtab-matlab`](https://github.com/liebermeister/sbtab-matlab) - SBtab toolbox
    - [`enzyme-cost-minimization`](https://github.com/liebermeister/enzyme-cost-minimization) - enzyme cost minimization toolbox
-  Make sure all the directories and subdirectories are included in your Matlab path
### Python
- Install using `pip`:
    - `cvxpy` ~= 1.1
    - `pint` ~= 0.16
    - `numpy` ~= 1.20
    - `scipy` ~= 1.6
    - `pandas` ~= 1.2
    - `Mosek` ~= 9.2

## License
This package is released under the [GNU General Public License](LICENSE).

## Contact
Please contact [Wolfram Liebermeister](mailto:wolfram.liebermeister@gmail.com) and [Elad Noor](mailto:elad.noor@weizmann.ac.il) with any questions or comments.

## References
Liebermeister W. (2019), *Model balancing: consistent in-vivo kinetic constants and metabolic states obtained by convex optimisation*
[bioRxiv doi:10.1101/2019.12.23.887166v1](https://www.biorxiv.org/content/10.1101/2019.12.23.887166v1)
