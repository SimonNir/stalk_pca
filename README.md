# Surrogate Theory Accelerated Line-search Kit (STALK)

**Note!** *This is NOT the official release but rather an experimental version to play around with
incorporating automatic parameterization and transition pathway functionalities.

For a well-tested implementation of the automated transition pathway workflow, please refer to https://github.com/SimonNir/legacy_surrogate_hessian_pca_workflow.git
*

Surrogate Theory Accelerated Line-search Kit (STALK) is a Python implementation of The
Surrogate Hessian Accelerated Parallel Line-search method. The method is intended for
optimizing and performing energy minimization of atomic structures in the presence of
statistical noise.

**Note!** *The repository was renamed 'surrogate_hessian_relax'->'stalk' on 19 Dec 2025, and
the code usage has changed substantially upon python packaging. To complete projects in the
old code base, keep using [v0.1](https://github.com/QMCPACK/stalk/releases/tag/v0.1) or
reach out for help in migration. We regret any inconvenience and anticipate fewer breaking
changes in the future.*

## CITING

Upon publishing results based on the method, we kindly ask you to cite [The original
work](https://doi.org/10.1063/5.0079046).

> Juha Tiihonen, Paul R. C. Kent, and Jaron T. Krogel \
The Journal of Chemical Physics \
156, 054104 (2022)

For results based on the transition pathway or saddle search methods, we kindly ask you to cite [the work introducing these methods](https://pubs.acs.org/doi/10.1021/acs.jctc.4c00214). 
> Gopal R. Iyer, Noah Whelpley, Juha Tiihonen, Paul R. C. Kent, Jaron T. Krogel, and Brenda Rubenstein \
The Journal of Chemical Theory and Computation \
2024, 20, 17, 7416–7429

<!-- And, if utilizing the automatic parameterization pipeline for transition paths, we kindly ask you to cite [the work introducing it]().
> TBD
-->


## SUPPORT

Installation instructions, examples and other information can be found in the online
[Documentation](https://stalk.readthedocs.io/en/latest/).

The software and its documentation are under development with no warranties. Support may be
inquired by [contacting the authors](mailto:tiihonen@iki.fi).

## ACKNOWLEDGEMENTS

The authors of this method are Juha Tiihonen, Paul R. C. Kent and Jaron T.
Krogel, working in the Center for Predictive Simulation of Functional Materials
(https://cpsfm.ornl.gov/)

This work has been authored in part by UT-Battelle, LLC, under contract
DE-AC05-00OR22725 with the US Department of Energy (DOE).
