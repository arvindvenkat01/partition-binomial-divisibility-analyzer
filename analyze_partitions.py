#!/usr/bin/env python3

"""
================================================================================
A script to find and analyze solutions to the divisibility condition
p(n) | C(2n,n), where p(n) is the integer partition function and C(2n,n)
is the central binomial coefficient.

This script uses a partial factorization method to efficiently test for
divisibility up to a large N and can perform congruence analysis on the results.
================================================================================

__author__ = "Arvind N. Venkat"
__email__ = "arvind.venkat01@gmail.com"
__date__ = "September 18, 2025"

# ==============================================================================
# Copyright 2025 Arvind N. Venkat

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE

# SOFTWARE.
# ==============================================================================

# Required libraries: pandas, sympy, numba
# Install with: pip install pandas sympy numba

# Usage: python <script_name>.py

"""

import time
from pathlib import Path
import pandas as pd
import sympy as sp

# ---------------------------
# Optional Numba acceleration
# ---------------------------
try:
    from numba import njit
except Exception:
    def njit(func=None, **kwargs):
        if func is None:
            return lambda f: f
        return func

# =================================
# p-adic valuations via Legendre/Kummer
# =================================

@njit
def legendre_vp_factorial(n, p):
    """v_p(n!) by Legendre's formula."""
    v = 0
    pp = p
    while pp <= n:
        v += n // pp
        pp *= p
    return v

@njit
def vp_binomial(n, k, p):
    """v_p(C(n,k)) using Legendre's formula."""
    return legendre_vp_factorial(n, p) - legendre_vp_factorial(k, p) - legendre_vp_factorial(n - k, p)

# =================================
# CRT helpers
# =================================

def egcd(a, b):
    if a == 0:
        return b, 0, 1
    d, x1, y1 = egcd(b % a, a)
    x = y1 - (b // a) * x1
    y = x1
    return d, x, y

def modinv(a, m):
    d, x, _ = egcd(a, m)
    if d != 1:
        raise ValueError("no inverse")
    return x % m

def crt_combine(remainders, moduli):
    if len(remainders) != len(moduli):
        raise ValueError("length mismatch")
    M = 1
    for m in moduli:
        M *= m
    x = 0
    for ai, mi in zip(remainders, moduli):
        Mi = M // mi
        yi = modinv(Mi, mi)
        x = (x + ai * yi * Mi) % M
    return x

# =================================
# nCr mod p^k (with p-parts handled)
# =================================

def nCr_mod_prime_power(n, r, p, k):
    """
    Compute C(n,r) mod p^k using Legendre valuations plus factorials with p-parts removed.
    """
    pk = p**k
    v = legendre_vp_factorial(n, p) - legendre_vp_factorial(r, p) - legendre_vp_factorial(n - r, p)
    if v >= k:
        return 0

    cache = {}

    def fact_no_p(m):
        if m == 0:
            return 1
        if m in cache:
            return cache[m]

        # Product of units modulo p^k: -1 for odd p or 4, +1 for 2^k (k>=3)
        block_prod = -1 if (p > 2 or pk == 4) else 1

        res = pow(block_prod, m // pk, pk)
        tail = m % pk
        for a in range(1, tail + 1):
            if a % p != 0:
                res = (res * a) % pk
        res = (res * fact_no_p(m // p)) % pk
        cache[m] = res
        return res

    num  = fact_no_p(n)
    den1 = fact_no_p(r)
    den2 = fact_no_p(n - r)
    inv_den = (modinv(den1, pk) * modinv(den2, pk)) % pk
    ans = (num * inv_den) % pk
    ans = (ans * pow(p, v, pk)) % pk
    return ans

# =======================================
# Partition data manager with partial factors (≤ 2n)
# =======================================

class PartitionDataManager:
    """
    Store partition values p(n) in chunked pickle files, including partial factorization
    up to 2n only, which is all that is needed for obstructions and, when residual==1,
    for CRT confirmation as well.
    """

    def __init__(self, data_dir="partition_data", default_chunk_size=1000):
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.default_chunk_size = default_chunk_size

    def compute_and_save_partitions(self, max_n, chunk_size=None):
        if chunk_size is None:
            chunk_size = self.default_chunk_size
        print(f"Computing partition values up to n={max_n}...")

        for chunk_start in range(0, max_n + 1, chunk_size):
            chunk_end = min(chunk_start + chunk_size - 1, max_n)
            chunk_file = self.data_dir / f"partitions_{chunk_start}_{chunk_end}.pkl"
            if chunk_file.exists():
                print(f"Chunk {chunk_start}-{chunk_end} exists, skipping...")
                continue

            print(f"Computing chunk {chunk_start}-{chunk_end}...")
            t0 = time.time()
            rows = []
            for n in range(chunk_start, chunk_end + 1):
                pn = sp.partition(n)

                # Partial factorization up to 2n
                raw = sp.factorint(pn, limit=2 * n) if n > 0 else {}

                # Keep only true primes ≤ 2n
                small_primes = {}
                small_prod = 1
                for base, exp in raw.items():
                    if sp.isprime(base) and base <= 2 * n:
                        b = int(base)
                        e = int(exp)
                        small_primes[b] = e
                        small_prod *= pow(b, e)

                # Residual > 1 certifies a factor > 2n
                residual = pn // small_prod
                residual_gt1 = residual > 1

                # Max small prime among primes ≤ 2n only
                max_small_prime = max(small_primes) if small_primes else None

                rows.append({
                    "n": n,
                    "p_n": pn,                                 # big int, keep as object dtype
                    "p_n_digits": len(str(pn)),
                    "small_factors": small_primes,             # dict of primes ≤ 2n only
                    "residual_gt1": bool(residual_gt1),        # True ⇒ impossible by prime bound
                    "max_small_prime": max_small_prime         # None or ≤ 2n
                })

            df = pd.DataFrame(rows)

            # Optional: keep nullable integer dtype for max_small_prime; safe since ≤ 2n
            if "max_small_prime" in df.columns:
                try:
                    df["max_small_prime"] = pd.array(df["max_small_prime"], dtype="Int64")
                except Exception:
                    pass

            df.to_pickle(chunk_file)
            print(f"  Saved {len(df)} rows in {time.time()-t0:.2f}s")

    def _sorted_chunks(self):
        files = list(self.data_dir.glob("partitions_*.pkl"))
        def keyfunc(p):
            name = p.stem  # partitions_START_END
            _, s, e = name.split("_")
            return (int(s), int(e))
        return sorted(files, key=keyfunc)

    def detect_chunk_size(self):
        files = self._sorted_chunks()
        if not files:
            return self.default_chunk_size
        first = files[0].stem.split("_")
        if len(first) == 3:
            s, e = int(first[1]), int(first[2])
            return e - s + 1
        return self.default_chunk_size

    # -------- schema normalization helpers --------

    def _ensure_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Normalize schema for legacy chunk files:
        - ensure p_n_digits exists
        - ensure small_factors, residual_gt1, max_small_prime exist (compute via partial factorization up to 2n)
        """
        # Ensure p_n and n exist
        if "n" not in df.columns or "p_n" not in df.columns:
            raise ValueError("Chunk file missing required columns 'n' or 'p_n'")

        # p_n_digits
        if "p_n_digits" not in df.columns:
            df["p_n_digits"] = df["p_n"].apply(lambda x: len(str(int(x))))

        need_small = ("small_factors" not in df.columns) or ("residual_gt1" not in df.columns) or ("max_small_prime" not in df.columns)
        if need_small:
            # initialize columns if missing
            if "small_factors" not in df.columns:
                df["small_factors"] = pd.NA
            if "residual_gt1" not in df.columns:
                df["residual_gt1"] = pd.NA
            if "max_small_prime" not in df.columns:
                df["max_small_prime"] = pd.NA

            mask = df["small_factors"].isna() | df["residual_gt1"].isna() | df["max_small_prime"].isna()

            def compute_row(row):
                n = int(row["n"])
                pn = int(row["p_n"])
                if n == 0:
                    return pd.Series({"small_factors": {}, "residual_gt1": False, "max_small_prime": None})
                raw = sp.factorint(pn, limit=2 * n)
                small_primes = {}
                small_prod = 1
                for base, exp in raw.items():
                    if sp.isprime(base) and base <= 2 * n:
                        b = int(base)
                        e = int(exp)
                        small_primes[b] = e
                        small_prod *= pow(b, e)
                residual = pn // small_prod
                residual_gt1 = residual > 1
                max_small_prime = max(small_primes) if small_primes else None
                return pd.Series({
                    "small_factors": small_primes,
                    "residual_gt1": bool(residual_gt1),
                    "max_small_prime": max_small_prime
                })

            if mask.any():
                computed = df.loc[mask].apply(compute_row, axis=1)
                for col in ["small_factors", "residual_gt1", "max_small_prime"]:
                    df.loc[mask, col] = computed[col].values

        # dtype for max_small_prime
        if "max_small_prime" in df.columns:
            try:
                df["max_small_prime"] = pd.array(df["max_small_prime"], dtype="Int64")
            except Exception:
                pass
        return df

    # -------- load methods using normalization --------

    def load_partition_value(self, n):
        chunk_size = self.detect_chunk_size()
        chunk_start = (n // chunk_size) * chunk_size
        chunk_end = chunk_start + chunk_size - 1
        chunk_file = self.data_dir / f"partitions_{chunk_start}_{chunk_end}.pkl"
        if not chunk_file.exists():
            raise ValueError(f"Chunk {chunk_start}-{chunk_end} not found; compute first.")
        df = pd.read_pickle(chunk_file)
        df = self._ensure_columns(df)
        row = df[df["n"] == n]
        if row.empty:
            raise ValueError(f"n={n} not found in chunk {chunk_start}-{chunk_end}.")
        return row.iloc[0].to_dict()

    def load_range(self, start_n, end_n):
        frames = []
        for pkl in self._sorted_chunks():
            name = pkl.stem.split("_")
            s, e = int(name[1]), int(name[2])
            if e < start_n or s > end_n:
                continue
            df = pd.read_pickle(pkl)
            df = self._ensure_columns(df)
            df = df[(df["n"] >= start_n) & (df["n"] <= end_n)]
            if not df.empty:
                frames.append(df)
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    def get_statistics(self):
        stats = {
            "total_values": 0,
            "largest_n": 0,
            "largest_p_n_digits": 0,
            "chunks_stored": 0,
        }
        for pkl in self._sorted_chunks():
            df = pd.read_pickle(pkl)
            df = self._ensure_columns(df)
            stats["total_values"] += len(df)
            stats["largest_n"] = max(stats["largest_n"], int(df["n"].max()))
            stats["largest_p_n_digits"] = max(stats["largest_p_n_digits"], int(df["p_n_digits"].max()))
            stats["chunks_stored"] += 1
        return stats

# =======================================
# Obstruction checks (no full factorization needed)
# =======================================

def quick_divisibility_check(n, row):
    """
    Obstruction screen using only primes ≤ 2n stored in small_factors and residual flag.
    Returns (can_be_divisible, reason).
    """
    pn = int(row["p_n"])
    if n == 0 or pn == 1:
        return True, None  # trivial

    # Prime-bound obstruction: if residual_gt1 then some prime factor of p(n) exceeds 2n
    if bool(row["residual_gt1"]):
        return False, "prime_bound: residual>1 after removing primes<=2n"

    # Otherwise all prime factors of p(n) are ≤ 2n and are fully captured in small_factors
    sf = row["small_factors"]
    if isinstance(sf, dict):
        factors = {int(k): int(v) for k, v in sf.items()}
    else:
        factors = {int(k): int(v) for k, v in dict(sf).items()}

    # Valuation obstruction
    for p, k in factors.items():
        vp = int(vp_binomial(2 * n, n, int(p)))
        if vp < k:
            return False, f"valuation: v_{p}(C(2n,n))={vp} < {k}"

    return True, None

# =======================================
# Optional: CRT confirmation for candidates
# =======================================

def confirm_by_crt(n, row):
    """
    Confirm C(2n,n) mod p(n) using only small_factors when residual==1.
    Returns final remainder modulo p(n) or None if prime-bound obstruction applies.
    """
    pn = int(row["p_n"])
    if n == 0 or pn == 1:
        return 0

    if bool(row["residual_gt1"]):
        return None  # impossible by prime bound

    factors = {int(k): int(v) for k, v in row["small_factors"].items()}
    moduli, remainders = [], []
    for p, k in factors.items():
        mk = pow(int(p), int(k))
        r = nCr_mod_prime_power(2 * n, n, int(p), int(k))
        moduli.append(mk)
        remainders.append(r)
    rem = crt_combine(remainders, moduli)
    return rem % pn

# =======================================
# Range analysis workflow
# =======================================

def analyze_range(start_n, end_n, data_dir="partition_data"):
    pdm = PartitionDataManager(data_dir=data_dir)
    df = pdm.load_range(start_n, end_n)
    if df.empty:
        print("No data found for requested range; compute first.")
        return pd.DataFrame(), []

    results = []
    candidates = []
    for _, row in df.iterrows():
        n = int(row["n"])
        ok, reason = quick_divisibility_check(n, row)
        item = {
            "n": n,
            "p_n": int(row["p_n"]),
            "p_n_digits": int(row["p_n_digits"]),
            "residual_gt1": bool(row["residual_gt1"]),
            "can_divide": bool(ok),
            "failure_reason": None if ok else reason
        }
        results.append(item)
        if ok:
            candidates.append(n)
    res_df = pd.DataFrame(results)
    return res_df, candidates

# =======================================
# Use to run multiple times
# =======================================

if __name__ == "__main__":
    pdm = PartitionDataManager(data_dir="partition_data", default_chunk_size=1000)

    # Compute/store partitions and partial factors (≤ 2n) up to some N
    N = 20000
    pdm.compute_and_save_partitions(max_n=N, chunk_size=1000)

    print("\nStats:", pdm.get_statistics())

    # Analyze a range using only obstructions (no full factorization)
    res_df, cand = analyze_range(0, 20000, data_dir="partition_data")

    print("\nObstruction analysis summary:")
    print("  total analyzed:", len(res_df))
    print("  candidates (pass both obstructions):", cand)

    # Optional: confirm candidates by CRT where residual==1
    confirmed = []
    for n in cand:
        row = pdm.load_partition_value(n)
        rem = confirm_by_crt(n, row)
        confirmed.append((n, rem))
    print("\nCRT confirmations (remainder mod p(n)):")
    print(confirmed)
