"""Re-grid a probability density stored as ``y = -log p(x)`` on a regular grid.

Fixed-xmin variant of regrid_density.py: the support edge xmin is never moved,
so only the exponential tail is truncated and there is no `epsilon_left`.

A density is a `LogDensity`: a regular grid ``x``, values ``y = -log p``, and
two extrapolation policies that define the density off the interior:

  left  `PowerLawEdge(n)`     y(x) = y[1] - n*log((x - xmin)/dx)   xmin < x < x[1]
                              i.e. p ∝ (x - xmin)**n; n > 0 is a polynomial
                              approach to zero at the support edge xmin.
                              p == 0 for x <= xmin, so y[0] is formally +inf
                              and is never used.
  right `ExponentialTail(s)`  y(x) = y[-1] + s*(x - xmax)          x > xmax

The interior (x[1] .. xmax) is the log-linear interpolant of y, i.e. piecewise
exponential in p.  All masses below are exact for that model.

`regrid` returns the same kind of object on a coarse grid: xmin stays put --
the density vanishes there, so no mass is lost -- and only the exponential tail
is cut, at the 1-epsilon quantile.  The y values are copied from
the input, unshifted, so the density is preserved as given (no normalization).
The two policy parameters of the output are chosen to reproduce exactly the
input mass of the first cell and of the truncated tail.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def _phi(z):
    """(exp(z) - 1) / z, evaluated stably (limit 1 at z = 0)."""
    z = np.asarray(z, dtype=float)
    tiny = np.abs(z) < 1e-8
    zs = np.where(tiny, 1.0, z)
    return np.where(tiny, 1.0 + 0.5 * z, np.expm1(zs) / zs)


@dataclass
class PowerLawEdge:
    """Left policy: p ∝ (x - xmin)**n between xmin and the second point."""

    n: float

    def __post_init__(self):
        if not self.n > -1.0:
            raise ValueError(f"edge exponent n={self.n:.4g} must be > -1 "
                             "for the edge cell to have finite mass")


@dataclass
class ExponentialTail:
    """Right policy: p ∝ exp(-slope*(x - xmax)) beyond the last point."""

    slope: float

    def __post_init__(self):
        if not self.slope > 0.0:
            raise ValueError(f"tail slope={self.slope:.4g} must be > 0 "
                             "for the tail to have finite mass")


@dataclass
class LogDensity:
    x: np.ndarray
    y: np.ndarray          # -log p, y[0] unused as hard p=0 boundary at xmin
    left: PowerLawEdge
    right: ExponentialTail

    # ---------------------------------------------------------------- setup
    def __post_init__(self):
        self.x = np.asarray(self.x, dtype=float)
        self.y = np.asarray(self.y, dtype=float)
        if self.x.shape != self.y.shape or self.x.ndim != 1 or self.x.size < 3:
            raise ValueError("x and y must be 1-D of equal length >= 3")
        d = np.diff(self.x)
        if np.ptp(d) > 1e-9 * abs(self.dx):
            raise ValueError("x must be a regular grid")
        self._build()

    def _build(self):
        # Cell masses in units of exp(-ref), so that exp(-y) stays at O(1).
        ref = float(np.min(self.y[1:]))
        ys = self.y - ref
        h, n = self.dx, self.left.n
        m_edge = np.exp(-ys[1]) * h / (1.0 + n)          # power-law edge cell
        s = np.diff(ys[1:]) / h
        m_int = np.exp(-ys[1:-1]) * h * _phi(-s * h)     # log-linear cells
        m_tail = np.exp(-ys[-1]) / self.right.slope      # exponential tail
        cum = np.cumsum(np.concatenate(([0.0, m_edge], m_int, [m_tail])))
        self._log_mass = float(np.log(cum[-1]) - ref)
        # The CDF at the grid points: _cdf[i] is the *fraction* of the mass
        # below x[i] for i = 0 .. size-1, and the extra last entry is
        # cdf(+inf) = 1.  So _cdf[1] is the mass of the edge cell and
        # 1 - _cdf[-2] the tail's.
        self._cdf = cum / cum[-1]

    @property
    def dx(self):
        return float(self.x[1] - self.x[0])

    @property
    def xmin(self):
        return float(self.x[0])

    @property
    def xmax(self):
        return float(self.x[-1])

    @classmethod
    def from_grid(cls, x, y, edge_window=4, tail_window=4):
        """Build a LogDensity, fitting both policies to the outermost points."""
        x, y = np.asarray(x, float), np.asarray(y, float)
        k = int(np.clip(edge_window, 2, y.size - 1))
        lg = np.log(np.arange(1, k + 1))            # log((x_i - xmin)/dx)
        n = -float(lg @ (y[1:1 + k] - y[1]) / (lg @ lg))
        m = int(np.clip(tail_window, 2, y.size - 1))
        slope = float(np.polyfit(x[-m:], y[-m:], 1)[0])
        return cls(x, y, PowerLawEdge(n), ExponentialTail(slope))

    # ----------------------------------------------------------- evaluation
    def __call__(self, xq):
        """-log p at arbitrary points, using both policies outside the grid."""
        xq = np.atleast_1d(np.asarray(xq, dtype=float))
        out = np.empty_like(xq)
        out[:] = np.interp(xq, self.x, self.y)
        edge = (xq > self.xmin) & (xq < self.x[1])
        out[edge] = self.y[1] - self.left.n * np.log((xq[edge] - self.xmin)
                                                    / self.dx)
        out[xq <= self.xmin] = np.inf                       # p = 0
        hi = xq > self.xmax
        out[hi] = self.y[-1] + self.right.slope * (xq[hi] - self.xmax)
        return out

    # ------------------------------------------------------------- integrals
    @property
    def log_mass(self):
        """log of the total mass (the density is kept unnormalized)."""
        return self._log_mass

    def cdf(self, xq):
        """Fraction of the total mass below xq."""
        h, n = self.dx, self.left.n
        if xq <= self.xmin:
            return 0.0
        if xq < self.x[1]:                                    # edge cell
            return self._cdf[1] * ((xq - self.xmin) / h) ** (1.0 + n)
        if xq >= self.xmax:                                   # exponential tail
            return 1.0 - (1.0 - self._cdf[-2]) * np.exp(
                -self.right.slope * (xq - self.xmax))
        # Fraction of cell i below xq, from the exact log-linear cell mass.
        i = min(int((xq - self.xmin) / h), self.y.size - 2)
        s = (self.y[i + 1] - self.y[i]) / h
        u = xq - self.x[i]
        return self._cdf[i] + (self._cdf[i + 1] - self._cdf[i]) * (
            u * float(_phi(-s * u)) / (h * float(_phi(-s * h))))

    def quantile(self, q):
        """Inverse CDF (analytic in every regime)."""
        if not 0.0 < q < 1.0:
            raise ValueError("q must be in (0, 1)")
        cdf, h, n = self._cdf, self.dx, self.left.n
        if q <= cdf[1]:                                       # edge cell with hard boundary at xmin
            return self.xmin + h * (q / cdf[1]) ** (1.0 / (1.0 + n))
        if q >= cdf[-2]:                                      # exponential tail past the soft xmax
            return self.xmax - np.log((1.0 - q) / (1.0 - cdf[-2])) \
                / self.right.slope
        # Find the interior cell i such that cdf[i] <= q < cdf[i+1].
        i = min(np.searchsorted(cdf, q, side="right") - 1, self.y.size - 2)
        s = (self.y[i + 1] - self.y[i]) / h
        remainder = q - cdf[i]
        # p(x[i]) as a remainder of the total, again from the exact cell mass
        a = (cdf[i + 1] - cdf[i]) / (h * float(_phi(-s * h)))
        if abs(s * h) < 1e-8:                                 # flat cell
            return self.x[i] + remainder / a
        return self.x[i] - np.log1p(-remainder * s / a) / s


# ----------------------------------------------------------------------------
def regrid(f: LogDensity, n_points, epsilon=1e-4):
    """Resample `f` onto `n_points` regular points, keeping its policies.

    The new range is [f.xmin, f.quantile(1 - epsilon)]: the support edge is
    kept, so the only mass dropped is the epsilon in the exponential tail.

    Returns the coarse LogDensity; pass it to `diagnostics` to see how well
    it tracks `f`.
    """
    if n_points < 3:
        raise ValueError("n_points must be >= 3")
    if not 0.0 < epsilon < 0.5:
        raise ValueError("epsilon must be in (0, 0.5)")

    x_new = np.linspace(f.xmin, f.quantile(1.0 - epsilon), n_points)
    y_new = f(x_new)                       # y_new[0] = inf: the support edge
    h_new = float(x_new[1] - x_new[0])

    # Policy parameters that reproduce the input mass of the edge cell and of
    # the discarded tail exactly.  Everything is a fraction of f's total mass,
    # so the densities are taken as exp(-y)/mass too; the tail fraction beyond
    # the last point is epsilon by construction, and f.cdf(f.xmin) is 0.
    m_edge = f.cdf(x_new[1])
    if not m_edge > 0.0:                   # only reachable by underflow
        raise ValueError("the first cell of the new grid carries no mass")
    # Edge cell, with u = (x - xmin)/h and p1 = p(x[1]):
    #   p(x) = p1 * u**n  ->  m_edge = p1*h * ∫_0^1 u**n du = p1*h/(1 + n)
    # hence n = p1*h/m_edge - 1  (n > 0 for a density vanishing at xmin).
    n_edge = np.exp(-(y_new[1] + f.log_mass)) * h_new / m_edge - 1.0
    slope = np.exp(-(y_new[-1] + f.log_mass)) / epsilon

    return LogDensity(x_new, y_new, PowerLawEdge(n_edge),
                      ExponentialTail(slope))


# ----------------------------------------------------------------------------
def diagnostics(before: LogDensity, after: LogDensity):
    """How faithfully `after` represents `before`.

    Only the two objects are needed, so this also compares any other pair of
    densities covering the same range.  The mass and quantile entries converge
    as `after` is refined; `max_abs_dlogp` does not -- see below.
    """
    f, g = before, after

    # Deviation of -log p on f's own nodes.  `..._away_from_edge` skips the two
    # cells next to xmin: there the interpolation error of n*log(x - xmin) is
    # independent of the spacing (~|n|/8 per cell), so the plain maximum
    # saturates while everything else converges like dx**2.
    def dev(mask):
        if not mask.any():
            return 0.0, float("nan")
        e = np.abs(g(f.x[mask]) - f.y[mask])
        return float(e.max()), float(f.x[mask][np.argmax(e)])

    interior = (f.x >= g.x[1]) & (f.x <= g.xmax)
    d_all, at_all = dev(interior)
    d_smooth, at_smooth = dev(interior & (f.x >= g.x[min(3, g.x.size - 1)]))
    return {
        "mass_covered": f.cdf(g.xmax),          # f.cdf(g.xmin) == 0
        "mass_ratio": float(np.exp(g.log_mass - f.log_mass)),
        "compression": float(f.x.size / g.x.size),
        "max_abs_dlogp": d_all,
        "argmax_abs_dlogp": at_all,
        "max_abs_dlogp_away_from_edge": d_smooth,
        "argmax_abs_dlogp_away_from_edge": at_smooth,
        "max_dq": max(abs(g.quantile(q) - f.quantile(q))
                      for q in (0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)),
    }


# ----------------------------------------------------------------------------
# Categorical slots 1-3 of a palette validated for all-pairs CVD separation.
_C_TRUE, _C_FINE, _C_NEW = "#2a78d6", "#eb6834", "#1baf7a"
_SURFACE, _INK, _INK_2, _GRIDLINE = "#fcfcfb", "#0b0b0b", "#52514e", "#d5d5d0"


def plot(fine: LogDensity, coarse: LogDensity, truth=None,
         path="regrid_density.png"):
    """Plot -log p of the true function, the fine grid and the re-gridded one.

    `truth` is an optional callable x -> -log p on the same (unnormalized)
    scale.  Both policy regions of the coarse density -- its power-law edge
    cell and its exponential tail -- are drawn dashed, so what is grid and
    what is extrapolation stays visible.  Writes `path` and returns the figure.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # The coarse grid can reach past the fine one (small epsilon), so plot the
    # union.  Sample geometrically away from xmin: there -log p diverges.
    span = max(fine.xmax, coarse.xmax) - fine.xmin
    xd = fine.xmin + span * np.geomspace(1e-7, 1.0, 4000)
    y_fine, y_new = fine(xd), coarse(xd)
    y_true = None if truth is None else np.asarray(truth(xd), dtype=float)

    fig, axes = plt.subplots(3, 1, figsize=(8.0, 9.0), constrained_layout=True,
                             gridspec_kw={"height_ratios": [3.0, 1.5, 2.2]})
    fig.patch.set_facecolor(_SURFACE)
    for ax in axes:
        ax.set_facecolor(_SURFACE)
        ax.grid(color=_GRIDLINE, lw=0.7, alpha=0.7)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(_GRIDLINE)
        ax.tick_params(colors=_INK_2, labelsize=9)

    def draw(ax, x, y, color, label, lw=2.0, alpha=1.0, dash_outside=None):
        """One series; the policy regions of `dash_outside` are drawn dashed."""
        if dash_outside is None:
            return ax.plot(x, y, color=color, lw=lw, alpha=alpha,
                           label=label, solid_capstyle="round")
        lo, hi = dash_outside.x[1], dash_outside.xmax
        body = (x >= lo) & (x <= hi)
        ax.plot(x[body], y[body], color=color, lw=lw, label=label,
                solid_capstyle="round")
        for part in (x < lo, x > hi):        # edge cell / exponential tail
            ax.plot(np.where(part, x, np.nan), y, color=color, lw=lw,
                    ls=(0, (4, 2.5)))

    # -- 1: the three functions ------------------------------------------
    ax = axes[0]
    if y_true is not None:
        draw(ax, xd, y_true, _C_TRUE, "true function", lw=5.0, alpha=0.30)
    draw(ax, xd, y_fine, _C_FINE, f"fine grid ({fine.x.size} points)",
         lw=1.7, dash_outside=fine)
    draw(ax, xd, y_new, _C_NEW, f"re-gridded ({coarse.x.size} points)",
         dash_outside=coarse)
    ax.plot(coarse.x[1:], coarse.y[1:], "o", ms=5.5, color=_C_NEW,
            mec=_SURFACE, mew=1.4, zorder=5)
    ax.axvline(coarse.xmax, color=_INK_2, lw=0.8, ls=":", zorder=1)
    ax.annotate(f"coarse xmax = q(1-$\\epsilon$) = {coarse.xmax:.2f}  ",
                (coarse.xmax, 0.04), xycoords=("data", "axes fraction"),
                ha="right", va="bottom", fontsize=8.5, color=_INK_2)
    bulk = xd > fine.xmin + 0.02 * span
    ax.set_ylim(np.min(y_fine[bulk]) - 1.0, np.max(y_fine[bulk]) + 2.0)
    ax.set_xlim(fine.xmin - 0.02 * span, fine.xmin + 1.02 * span)
    ax.set_ylabel("$-\\log p$", color=_INK, fontsize=10)
    ax.set_title("Re-gridding a log-density: dashed = extrapolation policy",
                 loc="left", color=_INK, fontsize=11.5)
    leg = ax.legend(loc="upper center", frameon=False, fontsize=9.5,
                    labelcolor=_INK_2)
    leg.set_zorder(6)

    # -- 2: residual against the reference -------------------------------
    ax = axes[1]
    ref, ref_name = (y_true, "true") if y_true is not None else (y_fine, "fine")
    ax.axhline(0.0, color=_INK_2, lw=0.9)
    if y_true is not None:
        draw(ax, xd, y_fine - ref, _C_FINE, "fine", lw=1.7,
             dash_outside=fine)
    draw(ax, xd, y_new - ref, _C_NEW, "re-gridded", dash_outside=coarse)
    inner = xd >= coarse.x[min(3, coarse.x.size - 1)]
    lim = max(0.05, 1.6 * np.max(np.abs((y_new - ref)[inner])))
    ax.set_ylim(-lim, lim)
    ax.set_xlim(*axes[0].get_xlim())
    ax.set_ylabel(f"$\\Delta(-\\log p)$ vs {ref_name}", color=_INK, fontsize=10)
    ax.set_xlabel("x", color=_INK, fontsize=10)
    ax.set_title(f"residual (clipped at $\\pm${lim:.2f}; the edge cell runs "
                 "off scale)", loc="left", color=_INK_2, fontsize=9.5)

    # -- 3: the edge, where the power-law policy lives --------------------
    ax = axes[2]
    if y_true is not None:
        draw(ax, xd, y_true, _C_TRUE, None, lw=5.0, alpha=0.30)
    draw(ax, xd, y_fine, _C_FINE, None, lw=1.7, dash_outside=fine)
    draw(ax, xd, y_new, _C_NEW, None, dash_outside=coarse)
    ax.plot(coarse.x[1:3], coarse.y[1:3], "o", ms=5.5, color=_C_NEW,
            mec=_SURFACE, mew=1.4, zorder=5)
    ax.set_xscale("log")
    ax.set_xlim(1e-5 * span, 0.35 * span)
    lo = float(np.min(y_fine[bulk]))
    ax.set_ylim(lo - 1.0, lo + 20.0)
    ax.set_xlabel("$x - x_{min}$   (log scale)", color=_INK, fontsize=10)
    ax.set_ylabel("$-\\log p$", color=_INK, fontsize=10)
    ax.set_title(f"edge cell: straight lines of slope $-n$, "
                 f"n = {fine.left.n:.2f} fine vs {coarse.left.n:.2f} re-gridded",
                 loc="left", color=_INK_2, fontsize=9.5)
    # Direct labels, stacked in the empty lower-left corner of this panel.
    for j, (txt, col) in enumerate((("true", _C_TRUE), ("fine", _C_FINE),
                                    ("re-gridded", _C_NEW))):
        ya = lo + 4.4 - 1.7 * j
        ax.plot([4e-5 * span], [ya], "s", ms=6, color=col, mec=_SURFACE,
                mew=1.0, clip_on=False)
        ax.annotate(f"  {txt}", (4e-5 * span, ya), fontsize=9, color=_INK_2,
                    va="center")

    fig.savefig(path, dpi=140, facecolor=_SURFACE)
    return fig


# ----------------------------------------------------------------------------
if __name__ == "__main__":
    from math import gamma as _gammafn

    # Bimodal: Gamma(k=3, theta=0.8) + Gaussian(5, 0.7), on a needlessly fine
    # grid.  p ~ x**2 at the left edge, so the true edge exponent is n = 2.
    k, theta, w = 3.0, 0.8, 0.65
    xmin, dx, npts = 0.0, 0.005, 2601

    def neg_log_p(xv):
        """-log p of the mixture, with an arbitrary offset to be preserved."""
        xv = np.asarray(xv, dtype=float)
        p = (w * xv ** (k - 1) * np.exp(-xv / theta)
             / (_gammafn(k) * theta ** k)
             + (1 - w) * np.exp(-0.5 * ((xv - 5.0) / 0.7) ** 2)
             / (0.7 * np.sqrt(2 * np.pi)))
        with np.errstate(divide="ignore"):
            return -np.log(p) + 4.2

    x = xmin + dx * np.arange(npts)
    f = LogDensity.from_grid(x, neg_log_p(x))
    local_slope = 1 / theta - (k - 1) / f.xmax     # d(-log p)/dx at xmax
    print(f"input : {npts} points on [{f.xmin}, {f.xmax}]\n"
          f"        edge n={f.left.n:.4f} (exact {k - 1:.4f}), "
          f"tail slope={f.right.slope:.4f} (exact {local_slope:.4f})\n"
          f"        log mass = {f.log_mass:.6f} (exact {-4.2:.6f})\n")

    for n_new in (9, 17, 33, 65):
        g = regrid(f, n_new, epsilon=1e-6)
        d = diagnostics(f, g)
        print(f"n={n_new:3d} x in [{g.xmin:.3f},{g.xmax:6.3f}] dx={g.dx:6.4f} "
              f"n_edge={g.left.n:6.3f} slope={g.right.slope:.3f} "
              f"mass_ratio={d['mass_ratio']:.6f} max|dq|={d['max_dq']:.2e}\n"
              f"       max|dy| = {d['max_abs_dlogp']:.2e} at "
              f"x={d['argmax_abs_dlogp']:5.2f} (edge cells), "
              f"{d['max_abs_dlogp_away_from_edge']:.2e} at "
              f"x={d['argmax_abs_dlogp_away_from_edge']:5.2f} beyond them")

    plot(f, regrid(f, 100, epsilon=1e-6), path="regrid_density_simple.png",
         truth=neg_log_p)
    print("\nwrote regrid_density_simple.png")
