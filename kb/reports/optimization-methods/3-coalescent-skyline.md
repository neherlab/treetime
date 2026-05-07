# Chapter 3: Coalescent skyline optimization

[Back to index](README.md) | Previous: [Chapter 2: Branch length optimization](2-branch-length-optimization.md) | Next: [Chapter 4: Outer loops](4-outer-loops.md)

The coalescent skyline estimates how the effective population size Ne(t) (or its inverse, the coalescent time scale Tc(t)) varied over time. TreeTime parameterizes this as `log(Tc)` on a grid and minimizes a penalized negative log-likelihood. This is a multi-dimensional optimization problem.

## Model evolution

| Model               | Year | Parameterization                             | Optimization                 | Reference                                                                                                         |
| :------------------ | :--- | :------------------------------------------- | :--------------------------- | :---------------------------------------------------------------------------------------------------------------- |
| Classic skyline     | 2000 | 1 param per coalescent event                 | MLE per interval             | <a id="cite-1"></a>[Pybus et al. 2000](https://doi.org/10.1093/genetics/155.3.1429) [[1](#ref-1)]                 |
| Generalized skyline | 2001 | AIC-grouped intervals                        | AIC selection                | <a id="cite-2"></a>[Strimmer and Pybus 2001](https://doi.org/10.1093/oxfordjournals.molbev.a003776) [[2](#ref-2)] |
| BSP                 | 2005 | Multiple change-point                        | MCMC                         | <a id="cite-3"></a>[Drummond et al. 2005](https://doi.org/10.1093/molbev/msi103) [[3](#ref-3)]                    |
| Skyride             | 2008 | GMRF on log(Ne) at coalescent times          | MCMC + block update          | <a id="cite-4"></a>[Minin et al. 2008](https://doi.org/10.1093/molbev/msn090) [[4](#ref-4)]                       |
| Skygrid             | 2013 | GMRF on regular grid, multi-locus            | MCMC                         | <a id="cite-5"></a>[Gill et al. 2013](https://doi.org/10.1093/molbev/mss265) [[5](#ref-5)]                        |
| TreeTime            | 2018 | Penalized likelihood, piecewise-linear Tc(t) | SLSQP (v0) / NelderMead (v1) | <a id="cite-6"></a>[Sagulenko et al. 2018](https://doi.org/10.1093/ve/vex042) [[6](#ref-6)]                       |
| mlesky              | 2023 | Cross-validated ML                           | ML optimization              | <a id="cite-7"></a>[Didelot et al. 2023](https://doi.org/10.1093/ve/vead028) [[7](#ref-7)]                        |
| Horseshoe MRF       | 2020 | Locally adaptive smoothing                   | MCMC                         | <a id="cite-8"></a>[Faulkner et al. 2020](https://doi.org/10.1111/biom.13276) [[8](#ref-8)]                       |

## GMRF prior

All modern implementations share a Gaussian Markov Random Field prior on `theta = log(Ne)` (denoted `gamma` in the algorithm inventory, following the Skyride convention; this report uses `theta` throughout):

```
p(theta | tau) ~ tau^((m-1)/2) * exp(-tau/2 * theta^T Q theta)
```

where Q is a tridiagonal precision matrix. Two constructions:

- **Time-weighted (Skyride)**: `Q[i,i+1] = -1/dt_i`. Penalizes rate of change proportionally to time interval.
- **Uniform (Skygrid, TreeTime)**: `Q[i,i+1] = -1`. Equivalent to time-weighted with uniform spacing.

TreeTime v1's stiffness penalty `stiffness * sum(diff(log_Tc)^2)` is the GMRF quadratic form with uniform weights. The boundary penalty (for log_Tc outside [-100, 0]) is TreeTime-specific.

## Gradient availability

The coalescent + GMRF gradient is analytically tractable:

- **GMRF penalty**: `d/d(theta_k) [theta^T Q theta] = 2*(Q*theta)_k` = tridiagonal Laplacian, O(n)
- **Coalescent likelihood**: `d/d(theta_k) = n_coal_k - C_k * exp(-theta_k)` where n_coal_k and C_k are precomputable sufficient statistics
- **Boundary penalty**: piecewise constant (subgradient)

<a id="cite-9"></a>[Lan et al. 2015](https://doi.org/10.1093/bioinformatics/btv378) [[9](#ref-9)] demonstrated that gradient-based methods (HMC) outperform gradient-free methods for coalescent parameter estimation. <a id="cite-10"></a>[Fourment et al. 2023](https://doi.org/10.1093/gbe/evad099) [[10](#ref-10)] showed hand-coded gradients are faster than autodiff for phylogenetics.

## NelderMead vs LBFGS

v1 uses NelderMead (derivative-free) for the skyline. v0 uses SLSQP (gradient-based via scipy). This is a known difference (`M-timetree-skyline-nelder-mead-optimizer`).

| Property              | Nelder-Mead          | L-BFGS                    | SLSQP (v0)             |
| :-------------------- | :------------------- | :------------------------ | :--------------------- |
| Convergence guarantee | None above 2D        | Global (uniformly convex) | Local superlinear      |
| Gradients needed      | No                   | Yes                       | Yes (FD or analytical) |
| Constraint handling   | Penalty terms        | L-BFGS-B: box constraints | Full (eq + ineq + box) |
| Dimension scaling     | Poor (>10D degrades) | O(mn), scales to millions | O(n^2), small-medium   |

Key results:

- <a id="cite-11"></a>[McKinnon 1998](https://doi.org/10.1137/S1052623496303482) [[11](#ref-11)]: NelderMead converges to a nonstationary point in 2D (explicit counterexample). No convergence guarantee above 2D.
- <a id="cite-12"></a>[Lagarias et al. 1998](https://doi.org/10.1137/S1052623496303470) [[12](#ref-12)]: NelderMead convergence proved only for 1D strictly convex.
- <a id="cite-13"></a>[Rios and Sahinidis 2013](https://doi.org/10.1007/s10898-012-9951-y) [[13](#ref-13)]: "The ability of all [derivative-free] solvers to obtain good solutions diminishes with increasing problem size."

v1's default 10 grid points means 10D optimization without convergence guarantees. Replacing NelderMead with LBFGS is well-supported (see [audit proposal P8](7-audit.md)).

## Code locations

**BEAST** ([beast-mcmc](https://github.com/beast-dev/beast-mcmc)):

- Skyride likelihood: `src/dr/evomodel/coalescent/GMRFSkyrideLikelihood.java`
- Block update operator (Newton+Cholesky): `src/dr/evomodel/coalescent/operators/GMRFSkyrideBlockUpdateOperator.java`
- Gradient for HMC: `src/dr/evomodel/coalescent/GMRFSkyrideGradient.java`

**skygrowth** ([skygrowth](https://github.com/mrc-ide/skygrowth)):

- MAP estimation via BFGS: `R/skygrowth.R`

**TreeTime v0** ([treetime](https://github.com/neherlab/treetime)):

- Skyline: [treetime/merger_models.py#L281-L338](https://github.com/neherlab/treetime/blob/master/treetime/merger_models.py#L281-L338) `optimize_skyline()` via `scipy.optimize.minimize(method='SLSQP')`

**TreeTime v1**:

- Skyline: [packages/treetime/src/commands/timetree/coalescent/skyline.rs#L68-L157](../../../packages/treetime/src/commands/timetree/coalescent/skyline.rs#L68-L157) `optimize_skyline()` via argmin `NelderMead`

## References

1. <a id="ref-1"></a> Pybus, Oliver G., Andrew Rambaut, and Paul H. Harvey. 2000. "An Integrated Framework for the Inference of Viral Population History from Reconstructed Genealogies." _Genetics_ 155(3):1429-1437. https://doi.org/10.1093/genetics/155.3.1429 [↩](#cite-1)
2. <a id="ref-2"></a> Strimmer, Korbinian, and Oliver G. Pybus. 2001. "Exploring the Demographic History of DNA Sequences Using the Generalized Skyline Plot." _Mol. Biol. Evol._ 18(12):2298-2305. https://doi.org/10.1093/oxfordjournals.molbev.a003776 [↩](#cite-2)
3. <a id="ref-3"></a> Drummond, Alexei J., Andrew Rambaut, Beth Shapiro, and Oliver G. Pybus. 2005. "Bayesian Coalescent Inference of Past Population Dynamics from Molecular Sequences." _Mol. Biol. Evol._ 22(5):1185-1192. https://doi.org/10.1093/molbev/msi103 [↩](#cite-3)
4. <a id="ref-4"></a> Minin, Vladimir N., Erik W. Bloomquist, and Marc A. Suchard. 2008. "Smooth Skyride through a Rough Skyline: Bayesian Coalescent-Based Inference of Population Dynamics." _Mol. Biol. Evol._ 25(7):1459-1471. https://doi.org/10.1093/molbev/msn090 [↩](#cite-4)
5. <a id="ref-5"></a> Gill, Mandev S., Philippe Lemey, Nuno R. Faria, Andrew Rambaut, Beth Shapiro, and Marc A. Suchard. 2013. "Improving Bayesian Population Dynamics Inference: A Coalescent-Based Model for Multiple Loci." _Mol. Biol. Evol._ 30(3):713-724. https://doi.org/10.1093/molbev/mss265 [↩](#cite-5)
6. <a id="ref-6"></a> Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-Likelihood Phylodynamic Analysis." _Virus Evol._ 4(1):vex042. https://doi.org/10.1093/ve/vex042 [↩](#cite-6)
7. <a id="ref-7"></a> Didelot, Xavier, Vinicius Franceschi, Simon D. W. Frost, Ann Dennis, and Erik M. Volz. 2023. "Model Design for Nonparametric Phylodynamic Inference and Applications to Pathogen Surveillance." _Virus Evol._ 9(1):vead028. https://doi.org/10.1093/ve/vead028 [↩](#cite-7)
8. <a id="ref-8"></a> Faulkner, James R., Andrew F. Magee, Beth Shapiro, and Vladimir N. Minin. 2020. "Horseshoe-Based Bayesian Nonparametric Estimation of Effective Population Size Trajectories." _Biometrics_ 76(3):677-690. https://doi.org/10.1111/biom.13276 [↩](#cite-8)
9. <a id="ref-9"></a> Lan, Shiwei, Julia A. Palacios, Michael D. Karcher, Vladimir N. Minin, and Babak Shahbaba. 2015. "An Efficient Bayesian Inference Framework for Coalescent-Based Nonparametric Phylodynamics." _Bioinformatics_ 31(20):3282-3289. https://doi.org/10.1093/bioinformatics/btv378 [↩](#cite-9)
10. <a id="ref-10"></a> Fourment, Mathieu, Christiaan J. Swanepoel, Jared G. Galloway, Xiang Ji, Karthik Gangavarapu, Marc A. Suchard, and Frederick A. Matsen IV. 2023. "Automatic Differentiation is no Panacea for Phylogenetic Gradient Computation." _Genome Biol. Evol._ 15(6):evad099. https://doi.org/10.1093/gbe/evad099 [↩](#cite-10)
11. <a id="ref-11"></a> McKinnon, Ken I. M. 1998. "Convergence of the Nelder-Mead Simplex Method to a Nonstationary Point." _SIAM J. Optim._ 9(1):148-158. https://doi.org/10.1137/S1052623496303482 [↩](#cite-11)
12. <a id="ref-12"></a> Lagarias, Jeffrey C., James A. Reeds, Margaret H. Wright, and Paul E. Wright. 1998. "Convergence Properties of the Nelder-Mead Simplex Method in Low Dimensions." _SIAM J. Optim._ 9(1):112-147. https://doi.org/10.1137/S1052623496303470 [↩](#cite-12)
13. <a id="ref-13"></a> Rios, Luis Miguel, and Nikolaos V. Sahinidis. 2013. "Derivative-Free Optimization: A Review of Algorithms and Comparison of Software Implementations." _J. Global Optim._ 56(3):1247-1293. https://doi.org/10.1007/s10898-012-9951-y [↩](#cite-13)
