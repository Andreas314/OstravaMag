#!/usr/bin/env python3
"""
Input file format  (everything after '#' is a comment, units: Angstrom, mu_B, eV)
---------------------------------------------------------------------------------
    lattice
        a1x a1y a1z
        a2x a2y a2z
        a3x a3y a3z

    atoms cartesian            # or: atoms fractional
        Cr   0.0  0.0  2.1378   0.0 0.0 -3.2576     # symbol x y z  mx my mz
        Cr   0.0  0.0  9.0496   0.0 0.0 -3.2576
        Cr   0.0  0.0  4.7741   0.0 0.0  3.2576
        Cr   0.0  0.0 11.6860   0.0 0.0  3.2576

    exchange
        -0.01127            # J of 1st shell (nearest neighbours)
        -0.00777            # J of 2nd shell
         0.00155            # ...
    # --- or, giving the shell distance explicitly (Angstrom, J in eV): ---
    #   2.6363  -0.01127
    #   2.8744  -0.00777
"""

import argparse
import pickle
import sys

import numpy as np
from ase import Atoms


# --------------------------------------------------------------------------
# parsing
# --------------------------------------------------------------------------
def parse_input(path):
    with open(path) as fh:
        raw = fh.readlines()

    lines = []
    for n, line in enumerate(raw, 1):
        line = line.split("#")[0].strip()
        if line:
            lines.append((n, line))

    cell, symbols, coords, moments, exchange = [], [], [], [], []
    frac = False
    section = None
    for n, line in lines:
        tok = line.split()
        key = tok[0].lower()
        if key in ("lattice", "atoms", "exchange"):
            section = key
            if key == "atoms":
                mode = tok[1].lower() if len(tok) > 1 else "cartesian"
                if mode in ("fractional", "frac", "direct", "crystal"):
                    frac = True
                elif mode in ("cartesian", "cart", "angstrom"):
                    frac = False
                else:
                    sys.exit(f"line {n}: unknown coordinate mode '{tok[1]}'")
            continue
        try:
            if section == "lattice":
                if len(tok) != 3:
                    raise ValueError("need 3 numbers")
                cell.append([float(x) for x in tok])
            elif section == "atoms":
                if len(tok) != 7:
                    raise ValueError("need: symbol x y z mx my mz")
                symbols.append(tok[0])
                coords.append([float(x) for x in tok[1:4]])
                moments.append([float(x) for x in tok[4:7]])
            elif section == "exchange":
                if len(tok) == 1:
                    exchange.append((None, float(tok[0]) / 1000))
                elif len(tok) == 2:
                    exchange.append((float(tok[0]), float(tok[1]) / 1000))
                else:
                    raise ValueError("need 'J' or 'distance J'")
            else:
                raise ValueError("data before any section keyword")
        except ValueError as e:
            sys.exit(f"{path}, line {n}: {e}:  '{line}'")

    if len(cell) != 3:
        sys.exit("Input error: 'lattice' section must contain exactly 3 vectors.")
    if not symbols:
        sys.exit("Input error: 'atoms' section is empty.")
    if not exchange:
        sys.exit("Input error: 'exchange' section is empty.")
    has_dist = [d is not None for d, _ in exchange]
    if any(has_dist) and not all(has_dist):
        sys.exit("Input error: in 'exchange' either give only J on every line "
                 "or 'distance J' on every line - not a mix.")

    cell = np.array(cell, float)
    coords = np.array(coords, float)
    pos = coords @ cell if frac else coords
    return cell, symbols, pos, np.array(moments, float), exchange


# --------------------------------------------------------------------------
# neighbour search
# --------------------------------------------------------------------------
def find_pairs(cell, pos, rcut):
    """All pairs (R, i, j) of atoms with 0 < |r_j + R.cell - r_i| <= rcut."""
    inv = np.linalg.inv(cell)
    spacing = 1.0 / np.linalg.norm(inv, axis=0)       # distance between lattice planes
    span = np.max(np.linalg.norm(pos[:, None] - pos[None, :], axis=2))
    nR = np.ceil((rcut + span) / spacing).astype(int) + 1

    grid = np.array(np.meshgrid(*[np.arange(-n, n + 1) for n in nR],
                                indexing="ij")).reshape(3, -1).T
    shifts = grid @ cell

    Rs, Is, Js, Ds = [], [], [], []
    n = len(pos)
    for i in range(n):
        for j in range(n):
            vec = pos[j] + shifts - pos[i]
            d = np.linalg.norm(vec, axis=1)
            ok = (d > 1e-6) & (d <= rcut)
            Rs.append(grid[ok])
            Is.append(np.full(ok.sum(), i))
            Js.append(np.full(ok.sum(), j))
            Ds.append(d[ok])
    return (np.concatenate(Rs), np.concatenate(Is),
            np.concatenate(Js), np.concatenate(Ds))


def cluster_shells(dist, tol):
    """Return (shell_id per pair, list of shell mean distances), sorted by distance."""
    order = np.argsort(dist)
    ids = np.empty(len(dist), int)
    shell, start = 0, 0
    sd = dist[order]
    for k in range(1, len(sd) + 1):
        if k == len(sd) or sd[k] - sd[k - 1] > tol:
            ids[order[start:k]] = shell
            shell += 1
            start = k
    means = [dist[ids == s].mean() for s in range(shell)]
    return ids, means


def neighbours_up_to(cell, pos, tol, n_shells=None, dmax=None):
    """Grow the search radius until the requested shells are complete."""
    rcut = 1.5 * np.max(np.linalg.norm(cell, axis=1))
    for _ in range(30):
        R, I, J, D = find_pairs(cell, pos, rcut)
        if len(D):
            ids, means = cluster_shells(D, tol)
            done = True
            if n_shells is not None:
                if len(means) < n_shells + 1:       # +1 -> shell N surely complete
                    done = False
            if dmax is not None and rcut < dmax + 2 * tol:
                done = False
            if done:
                return R, I, J, D, ids, means
        elif dmax is not None and rcut >= dmax + 2 * tol:
            return R, I, J, D, np.array([], int), []
        rcut *= 1.5
    sys.exit("Could not find the requested neighbour shells (search radius grew too large).")


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="text file with lattice, atoms+moments, exchange")
    ap.add_argument("-o", "--output", default="TB2J.pickle",
                    help="output pickle file (default: TB2J.pickle)")
    ap.add_argument("--tol", type=float, default=1e-3,
                    help="distance tolerance (Angstrom) defining a shell (default 1e-3)")
    args = ap.parse_args()

    cell, symbols, pos, moments, exchange = parse_input(args.input)
    natoms = len(symbols)

    # --- spin bookkeeping: every atom is a spin ----------------------------
    index_spin = list(range(natoms))
    spinat = moments                       # shape (natoms, 3)
    atoms = Atoms(symbols=symbols, positions=pos, cell=cell, pbc=False)

    # --- neighbour shells --------------------------------------------------
    by_distance = exchange[0][0] is not None
    if by_distance:
        dmax = max(d for d, _ in exchange)
        R, I, J, D, ids, means = neighbours_up_to(cell, pos, args.tol, dmax=dmax)
        shell_of_entry = []
        for dist, _ in exchange:
            hits = [s for s, m in enumerate(means) if abs(m - dist) <= args.tol + 1e-9]
            if not hits:
                sys.exit(f"No neighbour shell at distance {dist} A (tol {args.tol}). "
                         f"Shell distances found: {np.round(means[:15], 4).tolist()} ...")
            shell_of_entry.append(hits[0])
    else:
        R, I, J, D, ids, means = neighbours_up_to(cell, pos, args.tol,
                                                  n_shells=len(exchange))
        shell_of_entry = list(range(len(exchange)))

    # --- fill exchange_Jdict ----------------------------------------------
    Jdict = {}
    print(f"{'shell':>5} {'dist [A]':>10} {'#pairs':>7} {'J [meV]':>14}")
    for (_, Jval), s in zip(exchange, shell_of_entry):
        sel = np.where(ids == s)[0]
        for k in sel:
            key = (tuple(int(x) for x in R[k]), int(I[k]), int(J[k]))
            Jdict[key] = np.float64(Jval)
        print(f"{s + 1:>5} {means[s]:>10.4f} {len(sel):>7d} {Jval * 1000:>14.6g}")

    Jdict = dict(sorted(Jdict.items()))

    out = {
        "spinat": spinat,
        "index_spin": index_spin,
        "atoms": atoms,            # atoms.cell holds the lattice vectors
        "exchange_Jdict": Jdict,
    }
    with open(args.output, "wb") as fh:
        pickle.dump(out, fh)
    print(f"\n{natoms} atom(s), {len(Jdict)} exchange entries "
          f"written to {args.output}")


if __name__ == "__main__":
    main()
