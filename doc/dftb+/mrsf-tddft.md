# Mixed-Reference Spin-Flip TD-DFTB: theory and verification

This note documents the theory behind the collinear spin-flip (SF) and
mixed-reference spin-flip (MRSF) TD-DFTB implementation in
`src/dftbp/timedep/linrespsf.F90`, and gives a concrete recipe for **checking
that it is correct**.

The reference for the *ab initio* method is

> S. Lee, M. Filatov, S. Lee, C. H. Choi, *J. Chem. Phys.* **149**, 104101
> (2018), and its supplementary material (referred to as **SI** below).

Throughout, equation labels `SI Eq. Sx.y` and `Eq. (3.11)` refer to that paper.

---

## Part I — Theory

### 1. Reference and orbital spaces

The reference is a collinear, spin-polarized high-spin state with two unpaired
electrons,

    n_alpha = n_beta + 2,    M_S = +1.

DFTB+ solves the alpha and beta channels separately (UHF-like), giving
spin-resolved MOs `c^sigma` and eigenvalues `eps^sigma_p`. Orbital spaces:

| space               | symbol | indices                       | occ. in M_S=+1 |
|---------------------|--------|-------------------------------|----------------|
| closed (doubly occ) | C      | i,j,... = 1 .. n_beta         | alpha + beta   |
| open (SOMO)         | O      | O1 = n_beta+1, O2 = n_beta+2  | alpha only     |
| virtual             | V      | a,b,... = n_alpha+1 .. N      | empty          |

**ROHF (shared) reference.** For genuine MRSF the alpha/beta spin-orbitals must
share spatial parts. A shared set is built by diagonalising the spin-averaged
converged Hamiltonian

    H = 1/2 (H^alpha + H^beta),   H^sigma = S c^sigma eps^sigma (c^sigma)^T S,

and re-expressing the spin-resolved Fock matrices in that basis,

    F^sigma = O^sigma eps^sigma (O^sigma)^T,   O^sigma = C^T S c^sigma.

`F^alpha`, `F^beta` are in general non-diagonal and enter the response in full.
In `Reference = Unrestricted` mode the native alpha/beta orbitals are used, so
`F^sigma` is diagonal (= eps^sigma); this is cheaper but <S^2> is then only
approximate.

The mixed-reference RDM is `rho^MR = 1/2 (rho^{M_S=+1} + rho^{M_S=-1})` (SI Eq.
S1.1), giving occupation 1/2 in each SOMO; idempotency is restored by the
mixed-spin orbitals s1, s2 (SI Eqs. S2.1-S2.2). In this TDA implementation the
fractional-occupation RDM is not propagated explicitly; the response is built
directly in the determinant basis below.

### 2. Spin-flip manifold and the exchange-only kernel

A spin-flip excitation moves one electron alpha-occ `i` -> beta-vir `a`
(lowers M_S by 1). The transition density `phi_i^alpha phi_a^beta` is
off-diagonal in spin, so the Hartree/Coulomb coupling vanishes identically and
only exchange survives. The TDA spin-flip response matrix is

    A^SF_{ia,jb} = delta_ij delta_ab (eps^beta_a - eps^alpha_i)
                 + delta_ij F^beta_ab|offdiag - delta_ab F^alpha_ij|offdiag
                 - sum_AB q^{alpha,ij}_A gamma^LR_AB q^{beta,ab}_B
                 + 2 sum_A W_A q^{ia}_A q^{jb}_A .

Building blocks:

* **Orbital-energy / Fock diagonal.** In the ROHF basis the full
  `delta_ij F^beta_ab - delta_ab F^alpha_ij` enters; in Unrestricted mode it
  reduces to `delta_ij delta_ab (eps^beta_a - eps^alpha_i)`.

* **Long-range exchange.** The exact-exchange integral `(ij|ab)` in the DFTB
  monopole approximation factorises into same-spin transition charges,

      (ij|ab)_LR ~= sum_AB q^{alpha,ij}_A gamma^LR_AB q^{beta,ab}_B,
      q^{sigma,pq}_A = 1/2 sum_{mu in A} ( c^sigma_{mu p}(S c^sigma_q)_mu
                                         + c^sigma_{mu q}(S c^sigma_p)_mu ),

  scaled by the LC fraction `c_exchange`. No cross-spin charges are needed to
  build A^SF -- the key DFTB simplification.

* **On-site spin (transverse W) term.** `W_A` is the atomic spin constant
  (d^2E/dm_A^2), `q^{ia}_A` the cross-spin spin-flip transition charge. This is
  the local (transverse magnetisation) part of the spin-flip kernel; its
  prefactor 2 is fixed by reference recovery (Section 8).

A^SF is real symmetric; negative eigenvalues (states below the high-spin
reference) are kept.

### 3. Configuration types

Following the manuscript Fig. 3, the SF configurations split by hole/particle
location:

| type | name | hole -> particle | character          |
|------|------|------------------|--------------------|
| I    | O->O | SOMO -> SOMO     | open-shell (G/D/L/R)|
| II   | C->O | closed -> SOMO   | single-SOMO        |
| III  | O->V | SOMO -> virtual  | single-SOMO        |
| IV   | C->V | closed -> virtual| four-open-shell    |

MRSF spin-purifies I, II, III completely and IV only partially.

### 4. Type I -- the SOMO pair (G, D, L, R)

Working out occupations from the M_S=+1 reference (both SOMOs alpha):

* `O1->O1`: O1(down) O2(up)  -> one electron in each SOMO = open  = **L**
* `O2->O2`: O1(up) O2(down)  -> open                            = **R**
* `O2->O1`: O1(up,down), O2 empty -> both in O1 = closed        = **G**
* `O1->O2`: O2(up,down) -> both in O2 = closed                  = **D**

(SI Fig. S2: G = ground-like, D = double, L = left, R = right.) The same
configurations arise from the M_S=-1 reference via beta->alpha flips, related by
the **sign settings** (SI p. 9-10):

    G_{+1} = -G_{-1},  D_{+1} = -D_{-1},  L_{+1} = R_{-1},  R_{+1} = L_{-1}.

Symmetrising the orbital Hessian with these relations (SI p. 10 B-blocks:
`B_GL = -B_GR`, `B_DL = -B_DR`, `B_LL = B_RR`) and applying the singlet/triplet
transform U' (SI Eq. S8.1) gives the pure spin combinations (SI Eqs. S8.7-S8.8):

    triplet:  (L + R)/sqrt(2)                       [one configuration]
    singlet:  G,  D,  (L - R)/sqrt(2)               [three configurations]

**Realization in the code.** The two references are spin partners, so
`A_{X,+1} = A_{X,-1}` and the symmetrisation is the TDA congruence

    A^MRSF = T^T A^SF T,

with the retained open column of T carrying +1/sqrt(2) on `O1->O1` (L) and
+/-1/sqrt(2) on `O2->O2` (R): + for triplet (L+R), - for singlet (L-R). The
singlet keeps G, D as their own columns; the triplet drops them. The symmetrised
open diagonal is the reference-averaged (coupled) element

    A^MRSF_OO = 1/2 (A^SF_LL + A^SF_RR) +/- A^SF_LR,

not an isolated Slater-Condon matrix element.

### 5. Spin-pairing coupling (Eq. 3.11 / Fig. 5) -- types II & III

Single-SOMO configurations are spin-contaminated (<S^2> ~ 1) on a single
reference. MRSF pairs each M_S=+1 configuration |+> with its M_S=-1 partner |->
(SI Eq. S8.1) and couples them by the spin-pairing matrix element (Eq. 3.11):

    C = c_H <Psi^{M_S=+1}| H |Psi^{M_S=-1}>,   c_H = exact-exchange fraction.

**Derivation of C by Slater-Condon.** For a type-II config "closed i -> SOMO u",
spectator SOMO v, the occupied spin-orbital sets are

    Phi_{+1}: {core, i_beta,  u_alpha, u_beta, v_alpha}
    Phi_{-1}: {core, i_alpha, u_alpha, u_beta, v_beta }

which differ in exactly two spin-orbitals {i_beta, v_alpha} <-> {i_alpha,
v_beta}. The two-orbital Slater-Condon rule gives

    <Phi_{+1}|H|Phi_{-1}> = (i_beta i_alpha | v_alpha v_beta)   [= 0, spin]
                          - (i_beta v_beta | v_alpha i_alpha)   [= K_iv]
                          = -K_iv,

with K_iv = (iv|iv) >= 0 the spatial exchange integral. Hence

    C = -c_H K,   K = sum_AB q^{vi}_A gamma^LR_AB q^{vi}_B >= 0,

the long-range exchange between the spectator SOMO v and the partner orbital
(i for type II, the virtual a for type III). In the code
`kExch = cExchange * q.gamma.q` and the off-diagonal is entered as -kExch (= +C).

**2x2 augmented blocks** (SI Eqs. S8.9-S8.10), in the {|+>, |->} basis:

    A = [[E, C], [C, E]],   S^2 = [[1, 1], [1, 1]],   C = -c_H K.

These commute, so the eigenvectors are fixed:

    triplet (symmetric)  (|+>+|->)/sqrt2:  E + C = E - c_H K,  <S^2> = 2
    singlet (antisym.)   (|+>-|->)/sqrt2:  E - C = E + c_H K,  <S^2> = 0

This coupling lifts the residual multiplet degeneracy of single-reference
SF-TDDFT. In practice all single-SOMO configs are augmented together with the
type-I reduced block, the whole matrix is diagonalised, and states are selected
by <S^2> for the requested multiplicity.

### 6. Expectation value of S^2 (SI Eq. S7.11)

From `S^2 = S_z(S_z+1) + S_- S_+`, with M_S = 0 the first term vanishes and

    <S^2> = 1 - X_G^2 - X_D^2 + 2 X_L X_R
          + 2 ( sum_II X^+ X^- + sum_III X^+ X^- + sum_IV X^+ X^- ).

For the plain SF / SOMO-pair-only case this reduces in the implementation to

    <S^2> = (X_{O1->O1} + X_{O2->O2})^2
          + sum_ia X_ia^2 ( [a in V] + [i in C] ),

classifying: C->V -> <S^2> = 2 (triplet); single C->O or O->V -> 1 (mixture);
adapted SOMO pair -> 0 (singlet) or 2 (triplet). For the spin-complete case
<S^2> is the expectation value of the augmented S^2 matrix (the [[1,1],[1,1]]
blocks), giving 2 / 0 for the symmetric / antisymmetric partner combinations.
Contributions of the missing type-IV configurations are not included.

### 7. Type IV (C->V) -- not yet spin-complete

A closed->virtual SF determinant has four open shells
(i_beta, a_beta, O1_alpha, O2_alpha). The two singly-excited determinants
available at the singles/TDA level span only the pure triplet plus a
singlet/quintet mixture; the pure C->V singlet requires the doubly-spin-flipped
determinant (a double excitation), absent here and in single-reference
SF-TDDFT. Only one of five missing type-IV configurations is recovered
(manuscript). These are left un-augmented, carry <S^2> ~ 2, are kept on the
triplet side, and excluded from singlet spectra. Their contribution to the
low-lying states is expected to be small (manuscript Sec. III A).

### 8. Reference recovery (fixes the W prefactor)

The three M_S components of the reference triplet are exactly degenerate (no
SOC), so the M_S=0 component -- the spin-adapted SOMO-pair triplet
(L+R)/sqrt(2) -- must return at zero excitation energy. The deviation is a
parameter-free diagnostic, printed on every run.

The SOMO spin-splitting eps^beta_O - eps^alpha_O has a local part
-sum_A W_A m_A^2 (m_A = atomic spin population) and a non-local long-range
exchange part. Exchange-only response cancels only the latter (for CH2 the
recovery state was stuck at 2.4 eV). The transverse self-coupling of the W term
contributes 1/2 c_W sum_A W_A m_A^2, so exact cancellation requires

    c_W = 2.

Numerical sweep: CH2 recovery = 1.26 / 0.09 / -1.12 eV for c_W = 1 / 2 / 3,
crossing zero near c_W = 2, dropping the error from 2.4 eV to ~0.1 eV.

### 9. Singlet/triplet separation (overall structure)

The full U' transform (SI Eq. S8.1) splits the combined M_S = +/-1 response into
two decoupled sets (SI Eqs. S8.5-S8.6):

    A_T X_T = omega_T X_T,    A_S X_S = omega_S X_S,

with the triplet blocks getting A + C and the singlet blocks A - C for the
pairing-coupled types (SI Eqs. S8.9-S8.10). The implementation realises this by
(i) the congruence for type I, (ii) partner augmentation + pairing for types
II/III, (iii) <S^2>-selection of the requested multiplicity.

### Summary of DFTB-specific approximations (to keep in mind)

1. Monopole factorisation of the exchange `(ij|ab) -> q.gamma^LR.q` and of the
   pairing `K -> q.gamma^LR.q`.
2. On-site W transverse term with c_W = 2, replacing the ab initio non-collinear
   XC kernel; justified by reference recovery, absent from the original paper.
3. TDA only (no B-matrix / de-excitations).
4. Type IV left partial (no double-excitation singlet), per the paper.
5. ROHF shared-orbital construction for exact <S^2>; Unrestricted mode uses
   native alpha/beta orbitals.

---

## Part II — How to check it

There are no published MRSF-TD-**DFTB** numbers, and the tight-binding energies
cannot be matched to *ab initio* MRSF directly. Verification therefore rests on
**internal-consistency checks any correct implementation must satisfy**, plus a
few qualitative comparisons to *ab initio* MRSF trends. The checks below are
ordered from cheapest/most diagnostic to most involved.

### Check 1 — Reference recovery (built in, parameter-free)

The M_S=0 component of the reference triplet must appear at ~0 eV. Every run
prints this as a diagnostic (`recoveryEnergy`).

* **Pass:** |recovery energy| <~ 0.1 eV with `c_W = 2`.
* **Diagnostic value:** turning the W term off, or using c_W != 2, drives it
  away from zero (CH2: 2.4 eV with no W; 1.26 / 0.09 / -1.12 eV for c_W=1/2/3).
  Reproducing this sweep confirms the W prefactor derivation (Section 8).

Run: `test/app/dftb+/hybrid/cluster/CH2-MRSF-SpinComplete`.

### Check 2 — Spin purity <S^2>

The single most important correctness signal. Compare plain SF to MRSF on the
same reference:

* **Plain SF** (`MixedReference = No`): states show <S^2> in {0, 1, 2}
  (contaminated).
* **MRSF singlet** (`Multiplicity = Singlet`, `SpinComplete = Yes`):
  selected states have <S^2> -> 0 (to ~1e-3 where the exchange K is sizable).
* **MRSF triplet** (`Multiplicity = Triplet`): selected states have <S^2> -> 2.

A clean separation of <S^2> into 0 vs 2 (with type-IV C->V states correctly
sitting at ~2 and excluded from the singlet list) is the defining MRSF feature.

Run: `Benzene-SF-LC` vs `Benzene-MRSF-Singlet` / `Benzene-MRSF-Triplet`.

### Check 3 — Spin-pairing sign / the 2-state model (Section 5)

This isolates the Eq.-(3.11) coupling and its sign, which decide the multiplet
ordering. For a single-SOMO pair the spin-complete eigenvalues must obey

    E(triplet) = E - c_H K,    E(singlet) = E + c_H K,    splitting = 2 c_H K > 0,

i.e. the triplet partner lies **below** its singlet partner by `2 c_H K`, with
K the spectator/partner long-range exchange. To verify:

1. Run the spin-complete singlet and triplet for the same system and locate a
   pair of states dominated by the same single-SOMO configuration.
2. Their energy difference should equal `2 c_H K`, with `c_H = c_exchange`.
3. Setting `c_exchange -> 0` (or a hypothetical pairing-off build) must collapse
   the splitting to zero (the MRSF(0) limit).

A dedicated regression that asserts this on a small system (e.g. a two-orbital
diradical) is the most direct test of the matrix elements.

### Check 4 — Multiplet degeneracy lifting (qualitative, vs the paper)

The manuscript's own validation is the Be atom: SF-TDDFT splits the
`1P`/`3P` (x,y) vs z components spuriously; the pairing coupling reduces that
splitting (3P split 0.233 eV, 1P split 0.223 eV in the paper). In DFTB the
absolute numbers differ, but the **trend must hold**: with the pairing coupling
on, the (x,y)/z splitting of a degenerate multiplet is smaller than in plain SF.
Pick a system with symmetry-required degeneracy and confirm MRSF reduces the
artificial splitting relative to SF.

### Check 5 — Singlet-triplet gap of a diradical (sign + magnitude)

For a ground-state triplet diradical (e.g. CH2), MRSF must give
`E(S) - E(T) > 0` (triplet ground state, Hund) with a magnitude in the right
ballpark for the chosen Slater-Koster set. For a ground-state singlet diradical
the sign flips. The sign is a hard, qualitative check; the magnitude is a
softer comparison to *ab initio* / experiment.

### Check 6 — Invariances

* **Translation invariance:** shifting all coordinates leaves every excitation
  energy and oscillator strength unchanged (the monopole dipole construction
  satisfies sum_A q^{pq}_A = delta_pq). Verify f and omega are invariant.
* **Rotation invariance:** total oscillator strengths invariant under rigid
  rotation.
* **Symmetry-forbidden transitions:** in a high-symmetry molecule the
  symmetry-forbidden states must have f ~ 0 (benzene D6h: f = 0; breaking the
  symmetry restores f ~ 1e-3..1e-1). This doubles as a check that f is not
  accidentally nonzero.

### Check 7 — ROHF vs Unrestricted consistency

`Reference = ROHF` and `Reference = Unrestricted` should give close excitation
energies when the UHF reference has little spin contamination. Large
disagreement flags either a broken shared-orbital construction or a strongly
spin-contaminated UHF reference. ROHF is required for exact <S^2>.

### Check 8 — Finite-difference gradients

`tools/misc/mrsf_numgrad.py` parses `autotest.tag` (mermin_energy +
exc_energies_sqr) and forms E_tot = E_ref + omega(state). Central differences of
E_tot vs displaced single-points verify the total excited-state energy surface
is smooth and that the excitation energy responds correctly to nuclear motion
(d omega/dR dominates). A successful geometry optimisation of an excited state
(CH2) on this surface is an end-to-end check.

### Check 9 — Numerical sanity (cheap asserts)

* A^SF is symmetric (`max|A - A^T| < 1e-10`).
* The spin-adaptation T has orthonormal columns (`T^T T = I`).
* The augmented S^2 matrix has the [[1,1],[1,1]] partner structure and
  eigenvalues {0, 2} in each 2x2 block.
* Energies are real (Hermitian A) and the lowest MRSF triplet equals the
  recovery energy (~0).

### What a full validation campaign would add

To move from "internally consistent" to "quantitatively validated" one would:

* Implement the same molecule in an *ab initio* MRSF code (e.g. GAMESS) and
  compare **trends** (state ordering, S-T gaps, conical-intersection topology),
  not absolute energies.
* Add the missing type-IV singlet (double excitation) and confirm the residual
  C->V spin contamination is indeed small, as assumed.
* Replace the monopole exchange with a higher-multipole kernel and check the
  pairing splitting converges.

### Regression tests in the tree

    test/app/dftb+/hybrid/cluster/Benzene-SF-LC          # plain SF, contaminated <S^2>
    test/app/dftb+/hybrid/cluster/Benzene-MRSF-Singlet   # MRSF singlet, <S^2> -> 0
    test/app/dftb+/hybrid/cluster/Benzene-MRSF-Triplet   # MRSF triplet, <S^2> -> 2
    test/app/dftb+/hybrid/cluster/CH2-MRSF-SpinComplete  # spin-complete, recovery
    test/app/dftb+/hybrid/cluster/H2-dissoc-MRSF         # multireference / dissociation

These assert the excitation energies and (where available) <S^2> against stored
`_autotest.tag` references and are the first thing to run after any change.
