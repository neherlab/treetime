import numpy as np
from treetime import GTR

# generate a general GTR rate matrix
## equilibrium distributiontest
p = np.array([0.1, 0.15, 0.3, 0.4, 0.05])

## symmetric rate matrix
W = np.array([[0.0, 0.2, 0.5, 0.2, 0.1],
              [0.0, 0.0, 0.3, 0.5, 0.1],
              [0.0, 0.0, 0.0, 0.1, 0.1],
              [0.0, 0.0, 0.0, 0.0, 0.1],
              [0.0, 0.0, 0.0, 0.0, 0.0]])
W = W + W.T
testGTR = GTR.custom(W=W, pi=p)

## tests for the rate matrices and eigenvalues (set up of GTR model)
print("Rate Matrix:\n")
print(np.round(testGTR.Q,10))

print("\nEigenvalues:\n")
print(np.round(testGTR.eigenvals,10))

# tests for the eigenvectors
print("\nEigenvectors:\n")
print(np.round(testGTR.v,10))

print("\nInverse Eigenvectors:\n")
print(np.round(testGTR.v_inv,10))

# test consistency of eigenvectors and equilibrium distribution
print("\nOrthogonality of Eigenvectors:\n")
print(np.round(np.dot(testGTR.v_inv,testGTR.v),10))
print(np.abs(np.dot(testGTR.v_inv,testGTR.v) - np.eye(len(p))).sum()<1e-10)

print("\nEquilibrium distribution:")
eq_dis = testGTR.v[:,-1]
print(eq_dis)
print(np.abs(eq_dis - p).sum()<1e-10)

# test dynamics
print("\nexpQT:")
print(np.round(testGTR.expQt(0.1),10))

test_profile = np.array([[1, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0],
                         [0, 0, 1, 0, 0],
                         [0, 0, 0, 1, 0],
                         [0, 0, 0, 1, 1],
                         [0, 1, 1, 0, 0],
                         [1, 1, 1, 1, 1]], dtype=float)

# propagation tests -- this is backwards in time
# propagate for short time
propagated = testGTR.propagate_profile(test_profile, 0.1)
print("\nPropagated for short time:")
print(propagated)

# propagate for long time
propagated = testGTR.propagate_profile(test_profile, 1000)
print("\nPropagated for long time:")
print(propagated)

# evolution tests -- this is forward in time
test_profile = np.array([[1, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0],
                         [0, 0, 1, 0, 0],
                         [0, 0, 0, 1, 0],
                         [0, 0, 0, 0.5, 0.5],
                         [0, 0.8, 0.2, 0, 0],
                         p], dtype=float)

# evolve for short time
evolved = testGTR.evolve(test_profile, 0.1)
print("\nEvolved for short time:")
print(evolved)

# evolve for short time
evolved = testGTR.evolve(test_profile, 1000)
print("\nEvolved for long time:")
print(evolved)
print([np.abs(v - p).sum()<1e-10 for v in evolved])

