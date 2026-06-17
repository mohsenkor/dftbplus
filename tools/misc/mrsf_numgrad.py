#!/usr/bin/env python3
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2025  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
"""Finite-difference gradients and geometry optimisation for MRSF/SF-TD-DFTB.

Analytic excited-state gradients for (MR)SF-TD-DFTB are not yet available; this
helper provides numerical (central finite-difference) gradients of the total
excited-state energy

    E_tot(state) = E_ref + omega(state)

(the SCC energy of the high-spin reference plus the MRSF/SF excitation energy of
the requested state), and a simple optimiser built on them. It drives the normal
``dftb+`` binary, so the energies are exactly those of the running code.

The DFTB+ input must contain an ``ExcitedState/SpinFlip`` block and request the
tagged output (``Options { WriteAutotestTag = Yes }``). The geometry is taken
from a GenFormat file referenced by the input.

Example:
    mrsf_numgrad.py --state 1 --optimise geo.gen
"""

import argparse
import os
import shutil
import subprocess
import sys

BOHR__AA = 0.529177249
AA__BOHR = 1.0 / BOHR__AA


def read_gen(fname):
    """Read a GenFormat geometry. Returns (natom, types, species, coords[AA])."""
    with open(fname) as fp:
        lines = [ln for ln in fp.readlines()]
    natom, geotype = lines[0].split()[0:2]
    natom = int(natom)
    species = lines[1].split()
    types = []
    coords = []
    for ln in lines[2:2 + natom]:
        w = ln.split()
        types.append(int(w[1]))
        coords.append([float(w[2]), float(w[3]), float(w[4])])
    return natom, types, species, coords, geotype


def write_gen(fname, natom, types, species, coords, geotype):
    """Write a GenFormat geometry (coordinates in Angstrom)."""
    with open(fname, "w") as fp:
        fp.write(f"{natom}  {geotype}\n")
        fp.write(" " + " ".join(species) + "\n")
        for ii in range(natom):
            c = coords[ii]
            fp.write(f"{ii + 1} {types[ii]}  {c[0]:.10f}  {c[1]:.10f}  {c[2]:.10f}\n")


def parse_tag_value(fname, label):
    """Return the list of floats following a tag label in an autotest.tag file."""
    with open(fname) as fp:
        lines = fp.readlines()
    for idx, ln in enumerate(lines):
        if ln.split(":")[0].strip() == label:
            shape = ln.strip().split(":")
            # number of entries from the trailing shape field, e.g. ...:real:1:3
            count = 1
            if shape[-1].strip():
                count = 1
                for dim in shape[-1].split(","):
                    count *= int(dim)
            vals = []
            j = idx + 1
            while len(vals) < count and j < len(lines):
                vals += [float(x) for x in lines[j].split()]
                j += 1
            return vals
    raise KeyError(f"Label '{label}' not found in {fname}")


def energy(coords, ctx):
    """Total excited-state energy E_ref + omega(state) for the given coordinates."""
    write_gen(ctx["genfile"], ctx["natom"], ctx["types"], ctx["species"], coords,
              ctx["geotype"])
    with open(os.path.join(ctx["workdir"], "dftb.log"), "w") as log:
        subprocess.run([ctx["dftb"]], cwd=ctx["workdir"], stdout=log, stderr=log,
                       check=True)
    tag = os.path.join(ctx["workdir"], "autotest.tag")
    eref = parse_tag_value(tag, "mermin_energy")[0]
    omega = parse_tag_value(tag, "exc_energies_sqr")
    return eref + omega[ctx["state"] - 1]


def num_gradient(coords, ctx, step):
    """Central finite-difference gradient (Hartree/Angstrom)."""
    grad = [[0.0, 0.0, 0.0] for _ in range(ctx["natom"])]
    for ia in range(ctx["natom"]):
        for k in range(3):
            cp = [row[:] for row in coords]
            cm = [row[:] for row in coords]
            cp[ia][k] += step
            cm[ia][k] -= step
            ep = energy(cp, ctx)
            em = energy(cm, ctx)
            grad[ia][k] = (ep - em) / (2.0 * step)
    return grad


def gnorm(grad):
    return max(abs(g) for row in grad for g in row)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genfile", help="GenFormat geometry file referenced by dftb_in.hsd")
    ap.add_argument("--state", type=int, default=1,
                    help="index of the excited state (default: 1)")
    ap.add_argument("--step", type=float, default=1.0e-3,
                    help="finite-difference step in Angstrom (default: 1e-3)")
    ap.add_argument("--dftb", default="dftb+", help="path to the dftb+ binary")
    ap.add_argument("--workdir", default=".", help="working directory of the calculation")
    ap.add_argument("--optimise", action="store_true",
                    help="run a steepest-descent optimisation instead of a single gradient")
    ap.add_argument("--maxsteps", type=int, default=50)
    ap.add_argument("--gtol", type=float, default=1.0e-3,
                    help="max-component gradient convergence (Hartree/Angstrom)")
    args = ap.parse_args()

    natom, types, species, coords, geotype = read_gen(args.genfile)
    ctx = dict(genfile=args.genfile, natom=natom, types=types, species=species,
               geotype=geotype, dftb=shutil.which(args.dftb) or args.dftb,
               workdir=args.workdir, state=args.state)

    if not args.optimise:
        e0 = energy(coords, ctx)
        grad = num_gradient(coords, ctx, args.step)
        print(f"# state {args.state}  E_tot = {e0:.8f} H")
        print("# atom        d/dx           d/dy           d/dz   (H/Angstrom)")
        for ia in range(natom):
            g = grad[ia]
            print(f"{ia + 1:4d} {g[0]:15.8f}{g[1]:15.8f}{g[2]:15.8f}")
        print(f"# max|grad| = {gnorm(grad):.3e} H/Angstrom")
        return

    # robust steepest descent: capped backtracking line search
    maxdisp = 0.05  # max per-atom displacement per step (Angstrom)
    e_prev = energy(coords, ctx)
    for it in range(1, args.maxsteps + 1):
        grad = num_gradient(coords, ctx, args.step)
        gn = gnorm(grad)
        print(f"step {it:3d}  E = {e_prev:.8f} H   max|grad| = {gn:.3e} H/Ang")
        if gn < args.gtol:
            print("Converged.")
            break
        # scale the step so the largest displacement is at most maxdisp
        scale = min(1.0, maxdisp / gn)
        accepted = False
        for _ in range(8):
            trial = [[coords[ia][k] - scale * grad[ia][k] for k in range(3)]
                     for ia in range(natom)]
            e_trial = energy(trial, ctx)
            if e_trial < e_prev:
                coords = trial
                e_prev = e_trial
                accepted = True
                break
            scale *= 0.5  # backtrack
        if not accepted:
            print("Line search could not lower the energy; stopping.")
            break
    write_gen("opt.gen", natom, types, species, coords, geotype)
    print(f"Final E_tot = {e_prev:.8f} H ; geometry written to opt.gen")


if __name__ == "__main__":
    main()
