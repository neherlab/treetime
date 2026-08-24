"""Multiplying and dividing log-densities that live on regular grids.

Both operations are the same thing: evaluate every term on one common grid and
add its ``y = -log p`` values with a sign (+1 to multiply, -1 to divide), then
refit the two extrapolation policies to the result.

Boundaries, for products and quotients alike:

  hard left  xmin -> the most stringent (largest) of the terms.  Below it at
                     least one term has no support, so the product is zero and
                     the quotient undefined.  A quotient therefore cannot
                     recover a factor's support that the dividend has already
                     truncated -- see the second half of the demo.
  soft right xmax -> the most lenient (largest) of the terms, so no explicitly
                     sampled value is discarded; the shorter terms are carried
                     out there by their own exponential tail policy.

Spacing is the finest of the terms.  Because y[0] sits on the hard boundary
where p = 0, the combination there is inf + inf or inf - inf; it is set back to
inf, which is what the `LogDensity` convention already says y[0] means.
"""

from __future__ import annotations

import math
from functools import lru_cache
from types import SimpleNamespace

import numpy as np

from regrid_density_simple import LogDensity, diagnostics


def combine(terms, signs=None, edge_window=4, tail_window=4):
    """Add the -log p of `terms` with `signs` (+1 / -1) on a common grid."""
    terms = list(terms)
    signs = [1.0] * len(terms) if signs is None else list(signs)
    if len(terms) != len(signs) or not terms:
        raise ValueError("need one sign per term")

    lo = max(d.xmin for d in terms)          # most stringent hard boundary
    hi = max(d.xmax for d in terms)          # most lenient soft boundary
    dx = min(d.dx for d in terms)            # finest spacing
    if not hi > lo:
        raise ValueError("terms do not overlap")

    x = lo + dx * np.arange(int(np.ceil((hi - lo) / dx)) + 1)
    with np.errstate(invalid="ignore"):      # x[0] gives inf +/- inf
        y = sum(s * d(x) for s, d in zip(signs, terms))
    y[0] = np.inf                            # the hard boundary: p = 0
    # Both policies are refit from the combined values.  (The alternative is
    # exact policy algebra: the edge exponents of the terms sitting on the
    # hard boundary add with their signs, and so do the tail slopes -- but only
    # while every term's grid reaches the same xmax.)
    return LogDensity.from_grid(x, y, edge_window, tail_window)


def multiply(*terms, **kw):
    """Density proportional to the product of `terms`."""
    return combine(terms, **kw)


def divide(dividend, divisor, **kw):
    """Density proportional to `dividend` / `divisor`."""
    return combine([dividend, divisor], [1.0, -1.0], **kw)


# ----------------------------------------------------------------------------
def shifted_gamma(k, theta, x0):
    """-log p of a Gamma(k, theta) shifted to start at x0, as a callable."""
    const = math.lgamma(k) + k * math.log(theta)

    def neg_log_p(xv):
        u = np.asarray(xv, dtype=float) - x0
        with np.errstate(divide="ignore", invalid="ignore"):
            out = -(k - 1.0) * np.log(u) + u / theta + const
        return np.where(u > 0.0, out, np.inf)

    neg_log_p.k, neg_log_p.theta, neg_log_p.x0 = k, theta, x0
    return neg_log_p


def sample(neg_log_p, dx, xmax):
    """A LogDensity sampled from `neg_log_p` on [x0, xmax] with spacing dx."""
    x = neg_log_p.x0 + dx * np.arange(int(round((xmax - neg_log_p.x0) / dx)) + 1)
    return LogDensity.from_grid(x, neg_log_p(x))


def max_dev(a, b, lo, hi, n=20001):
    """max |Δ(-log p)| between two callables on [lo, hi], and where it occurs.

    Keep `lo` at or above both operands' second grid point: inside an edge cell
    the two power-law policies differ by Δn*|log u|, which diverges as u -> 0
    however small Δn is, and would drown out everything else.
    """
    x = np.linspace(lo, hi, n)
    e = np.abs(np.asarray(a(x), float) - np.asarray(b(x), float))
    return float(e.max()), float(x[np.argmax(e)])


# ----------------------------------------------------------------------------
# Tests -- multiplication and division only.  Properties of the LogDensity
# representation itself (policy fits, mass model, cdf/quantile) are tested in
# test_log_density.py.  Tolerances are the measured values with a little
# headroom, and are commented where the size of the number is itself the point.


@lru_cache(maxsize=None)
def _case():
    """Three shifted Gammas with different supports, spacings and ranges."""
    A = shifted_gamma(k=3.0, theta=1.2, x0=0.00)
    B = shifted_gamma(k=2.5, theta=0.8, x0=0.30)    # most stringent boundary
    C = shifted_gamma(k=4.0, theta=2.0, x0=0.15)
    fa = sample(A, 0.01, 25.0)
    fb = sample(B, 0.02, 18.0)
    fc = sample(C, 0.005, 30.0)
    prod = multiply(fa, fb, fc)                     # A*B*C
    return SimpleNamespace(
        A=A, B=B, C=C, fa=fa, fb=fb, fc=fc,
        prod=prod,
        quot=divide(prod, fc),                      # (A*B*C)/C  ->  A*B
        ref=multiply(fa, fb),                       # A*B, built directly
        truth=lambda xv: A(xv) + B(xv),             # analytic A*B
    )


def test_product_takes_stringent_left_lenient_right_finest_dx():
    c, terms = _case(), (_case().fa, _case().fb, _case().fc)
    assert c.prod.xmin == max(d.xmin for d in terms) == c.fb.xmin
    hi = max(d.xmax for d in terms)
    assert hi <= c.prod.xmax < hi + c.prod.dx      # ceil may overshoot by <dx
    assert np.isclose(c.prod.dx, min(d.dx for d in terms), rtol=1e-12)


def test_quotient_takes_stringent_left_lenient_right():
    c = _case()
    assert c.quot.xmin == max(c.prod.xmin, c.fc.xmin) == c.prod.xmin
    assert c.quot.xmax == max(c.prod.xmax, c.fc.xmax)


def test_combine_writes_inf_on_the_hard_boundary():
    """Multiplying gives inf + inf at x[0] and dividing inf - inf = nan, so
    combine() writes it back to inf.  That nothing ever *reads* y[0] is a
    property of the representation -- see test_log_density.py."""
    c = _case()
    for d in (c.prod, c.quot, c.ref):
        assert np.isposinf(d.y[0]), d.y[0]
        assert np.all(np.isfinite(d.y[1:]))


def test_quotient_equals_the_explicit_product_beyond_the_second_cell():
    """Away from the edge the two agree to machine precision everywhere, not
    just on the nodes."""
    c = _case()
    dev, _ = max_dev(c.quot, c.ref, c.ref.x[2], c.ref.xmax)
    assert dev < 1e-11, dev


def test_multiply_is_order_independent():
    c = _case()
    ab, ba = multiply(c.fa, c.fb), multiply(c.fb, c.fa)
    assert np.array_equal(ab.x, ba.x) and np.array_equal(ab.y, ba.y)


def test_multiply_then_divide_is_idempotent_after_the_first_cycle():
    """Re-gridding onto the finest spacing costs one interpolation; repeating
    the cycle changes nothing further."""
    c = _case()
    fd = sample(shifted_gamma(k=2.0, theta=1.5, x0=0.0), 0.004, 20.0)
    q1 = divide(multiply(c.fa, fd), fd)
    q2 = divide(multiply(q1, fd), fd)
    assert np.array_equal(q1.x, q2.x) and np.array_equal(q1.y, q2.y)
    assert (q1.left.n, q1.right.slope) == (q2.left.n, q2.right.slope)
    assert abs(q1.right.slope - c.fa.right.slope) < 1e-4


def test_dividing_out_the_boundary_factor_cannot_restore_support():
    """(A*B*C)/B keeps B's hard boundary: below it the dividend is already
    zero, so A*C's support there is gone for good."""
    c = _case()
    quot_b, ref_b = divide(c.prod, c.fb), multiply(c.fa, c.fc)
    assert quot_b.xmin == c.prod.xmin > ref_b.xmin
    xs = np.array([0.20, 0.25, 0.29])            # inside A*C, outside A*B*C
    assert np.all(np.isinf(quot_b(xs)))
    assert np.all(np.isfinite(c.A(xs) + c.C(xs)))
    # inside the surviving range it is still the right function
    dev, _ = max_dev(quot_b, lambda v: c.A(v) + c.C(v), quot_b.x[1],
                     ref_b.xmax)
    assert dev < 0.05, dev


# ----------------------------------------------------------------------------
def run_tests(namespace):
    """Run every test_* in `namespace`, in definition order.

    A stand-in for pytest, which collects the same functions unchanged.
    """
    cases = sorted(((f.__code__.co_firstlineno, n, f)
                    for n, f in list(namespace.items())
                    if n.startswith("test_") and callable(f)))
    failed = 0
    for _, name, fn in cases:
        try:
            fn()
            print(f"  PASS  {name}")
        except Exception as err:                                  # noqa: BLE001
            failed += 1
            print(f"  FAIL  {name}\n          {type(err).__name__}: {err}")
    print(f"\n{len(cases) - failed}/{len(cases)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(run_tests(globals()))
