from __future__ import annotations

import hashlib
import random
from itertools import product

import sympy as sp


N_DIMENSIONS = 4
RANDOM_SEED = 20260708


def simplify_expr(expr):
    return sp.cancel(sp.together(sp.expand(expr)))


class DiagMetric:
    """Curvature machinery for a diagonal metric over given coordinates."""

    def __init__(self, gdiag, coords):
        self.coords = coords
        self.g = {(i, i): gdiag[i] for i in range(N_DIMENSIONS)}
        self.gi = {(i, i): simplify_expr(1 / gdiag[i]) for i in range(N_DIMENSIONS)}
        self._build()

    def gd(self, a, b):
        return self.g.get((a, b), sp.Integer(0))

    def gu(self, a, b):
        return self.gi.get((a, b), sp.Integer(0))

    def _build(self):
        c = self.coords
        self.Gam = {}
        for a, b, d in product(range(N_DIMENSIONS), repeat=3):
            if b > d:
                continue
            expr = sum(self.gu(a, f) * (sp.diff(self.gd(f, b), c[d])
                                        + sp.diff(self.gd(f, d), c[b])
                                        - sp.diff(self.gd(b, d), c[f]))
                       for f in range(N_DIMENSIONS)) / 2
            expr = simplify_expr(expr)
            if expr != 0:
                self.Gam[(a, b, d)] = expr
                self.Gam[(a, d, b)] = expr

        G = lambda a, b, d: self.Gam.get((a, b, d), sp.Integer(0))
        self.G = G

        self.Riem = {}
        for a, b, cc, d in product(range(N_DIMENSIONS), repeat=4):
            if cc >= d:
                continue
            expr = (sp.diff(G(a, b, d), c[cc]) - sp.diff(G(a, b, cc), c[d])
                    + sum(G(a, cc, f) * G(f, b, d) - G(a, d, f) * G(f, b, cc)
                          for f in range(N_DIMENSIONS)))
            expr = simplify_expr(expr)
            if expr != 0:
                self.Riem[(a, b, cc, d)] = expr
                self.Riem[(a, b, d, cc)] = -expr
        Rm = lambda a, b, cc, d: self.Riem.get((a, b, cc, d), sp.Integer(0))

        self.Ric = {}
        for b, d in product(range(N_DIMENSIONS), repeat=2):
            expr = simplify_expr(sum(Rm(a, b, a, d) for a in range(N_DIMENSIONS)))
            if expr != 0:
                self.Ric[(b, d)] = expr
        self.Rc = lambda b, d: self.Ric.get((b, d), sp.Integer(0))
        self.Rs = simplify_expr(sum(self.gu(a, a) * self.Rc(a, a) for a in range(N_DIMENSIONS)))

        self.Weyl = {}
        for a, b, cc, d in product(range(N_DIMENSIONS), repeat=4):
            if a >= b or cc >= d or (a, b) > (cc, d):
                continue
            rabcd = sum(self.gd(a, f) * Rm(f, b, cc, d) for f in range(N_DIMENSIONS))
            expr = (rabcd
                    - sp.Rational(1, 2) * (self.gd(a, cc) * self.Rc(b, d)
                                           - self.gd(a, d) * self.Rc(b, cc)
                                           + self.gd(b, d) * self.Rc(a, cc)
                                           - self.gd(b, cc) * self.Rc(a, d))
                    + self.Rs / 6 * (self.gd(a, cc) * self.gd(b, d)
                                     - self.gd(a, d) * self.gd(b, cc)))
            expr = simplify_expr(expr)
            if expr != 0:
                self.Weyl[(a, b, cc, d)] = expr

    def C(self, a, b, cc, d):
        sgn = 1
        if a > b:
            a, b, sgn = b, a, -sgn
        if cc > d:
            cc, d, sgn = d, cc, -sgn
        if a == b or cc == d:
            return sp.Integer(0)
        if (a, b) > (cc, d):
            a, b, cc, d = cc, d, a, b
        return sgn * self.Weyl.get((a, b, cc, d), sp.Integer(0))

    def bach(self, a, b, memo):
        """B_ab = 2 D^m D^n C_{mabn} + C_{mabn} R^{mn}."""
        c, G, C = self.coords, self.G, self.C

        def DC(e, i, j, kk, l):
            key = (e, i, j, kk, l)
            if key in memo:
                return memo[key]
            expr = (sp.diff(C(i, j, kk, l), c[e])
                    - sum(G(f, e, i) * C(f, j, kk, l) for f in range(N_DIMENSIONS))
                    - sum(G(f, e, j) * C(i, f, kk, l) for f in range(N_DIMENSIONS))
                    - sum(G(f, e, kk) * C(i, j, f, l) for f in range(N_DIMENSIONS))
                    - sum(G(f, e, l) * C(i, j, kk, f) for f in range(N_DIMENSIONS)))
            expr = simplify_expr(expr)
            memo[key] = expr
            return expr

        def DDC(f, e, i, j, kk, l):
            expr = (sp.diff(DC(e, i, j, kk, l), c[f])
                    - sum(G(h, f, e) * DC(h, i, j, kk, l) for h in range(N_DIMENSIONS))
                    - sum(G(h, f, i) * DC(e, h, j, kk, l) for h in range(N_DIMENSIONS))
                    - sum(G(h, f, j) * DC(e, i, h, kk, l) for h in range(N_DIMENSIONS))
                    - sum(G(h, f, kk) * DC(e, i, j, h, l) for h in range(N_DIMENSIONS))
                    - sum(G(h, f, l) * DC(e, i, j, kk, h) for h in range(N_DIMENSIONS)))
            return simplify_expr(expr)

        term1 = 2 * sum(self.gu(m, m) * self.gu(nn, nn) * DDC(m, nn, m, a, b, nn)
                        for m in range(N_DIMENSIONS) for nn in range(N_DIMENSIONS))
        term2 = sum(self.gu(m, m) * self.gu(nn, nn) * self.Rc(m, nn) * C(m, a, b, nn)
                    for m in range(N_DIMENSIONS) for nn in range(N_DIMENSIONS))
        raw = term1 + term2
        return raw, sp.simplify(simplify_expr(raw))


def component_rng(case_id: str, component: tuple[int, int]) -> random.Random:
    material = f"{RANDOM_SEED}:{case_id}:{component[0]}:{component[1]}"
    seed = int.from_bytes(hashlib.sha256(material.encode("utf-8")).digest()[:8], "big")
    return random.Random(seed)


def sample_points(symbols_in_expr, radial_symbol, polar_symbol, beta_symbol, gamma_symbol, kappa_symbol, lambda_symbol, *, case_id: str, component: tuple[int, int]):
    """Three exact rational sample points for numerical cross-checks.

    Sampling is deterministic, component-keyed, and call-order independent.
    """
    rng = component_rng(case_id, component)
    points = []
    for _ in range(3):
        substitutions = {
            radial_symbol: sp.Rational(rng.randint(3, 40), rng.randint(1, 4)),
            polar_symbol: sp.pi / rng.choice([3, 4, 6]),
            beta_symbol: sp.Rational(rng.randint(1, 9), 100),
            gamma_symbol: sp.Rational(rng.randint(1, 9), 1000),
            kappa_symbol: sp.Rational(rng.randint(1, 9), 10000),
            lambda_symbol: sp.Rational(rng.randint(1, 9), 10000),
        }
        point = {symbol: value for symbol, value in substitutions.items() if symbol in symbols_in_expr or symbol in (radial_symbol, polar_symbol)}
        if point[radial_symbol] == 0:
            continue
        points.append(point)
    return points
