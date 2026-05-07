# Rate susceptibility citation references wrong paper section

## Summary

Four code locations cite "Sagulenko, Puller & Neher 2018, Section 2.5" for rate susceptibility analysis, but Section 2.5 covers the coalescent model, not rate susceptibility. The correct section appears to be Section 2.4. Requires paper verification.

## Details

The following locations reference Section 2.5:

- `representation/payload/timetree.rs:29:`
- `commands/timetree/output/confidence.rs:38:`
- `commands/timetree/output/confidence.rs:161:`
- `commands/timetree/output/confidence.rs:165:`

Two independent code reviews flagged this discrepancy. The rate susceptibility analysis (perturbing the clock rate and measuring date shifts) is described in Section 2.4 of Sagulenko et al. 2018, while Section 2.5 describes the coalescent tree prior.

## Required investigation

- Verify which section of Sagulenko, Puller & Neher (TreeTime: Maximum-likelihood phylodynamic analysis, Virus Evolution, 2018) describes the rate susceptibility / confidence interval method
- Update all four citations to reference the correct section
- Verify that the implementation matches the cited method

## Impact

- Incorrect citations make it harder to verify implementation correctness against the paper
- If the implementation was written following Section 2.5 instead of 2.4, the confidence interval method itself could be wrong
