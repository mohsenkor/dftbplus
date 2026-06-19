!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Collinear spin-flip (SF-TDDFT) and mixed-reference spin-flip (MRSF-TDDFT) excitations within the
!! Tamm-Dancoff approximation for range-separated (LC-)DFTB.
!!
!! On a spin-polarised, high-spin (M_S = +1) reference, the spin-flip manifold consists of
!! alpha-occupied -> beta-virtual single excitations. Because the spin-flip transition density is
!! off-diagonal in spin, the Hartree (gamma) coupling vanishes and the off-diagonal coupling is
!! carried entirely by the long-range exact (Hartree-Fock) exchange of LC-DFTB ("exchange-only"
!! structure of spin-flip / mixed-reference TDDFT). The TDA spin-flip A-matrix is
!!
!!   A_{ia,jb} = delta_ij F^beta_ab - delta_ab F^alpha_ij
!!                 - sum_{AB} q^{ij}_A gamma^LR_{AB} q^{ab}_B
!!                 + 2 sum_A W_A q^{ia}_A q^{jb}_A ,
!!
!! where the last term is the on-site spin (W) transverse coupling - the local part of the
!! spin-flip kernel. Together with the long-range exchange it cancels the SOMO spin-splitting, so the
!! M_S=0 component of the reference triplet returns at (near) zero excitation energy. This
!! "reference recovery" is printed as a diagnostic and is the method's built-in spin-consistency
!! check (e.g. for CH2 it drops from 2.4 eV without the W term to ~0.1 eV with it).
!!
!! Two references are supported:
!!   * ROHF (default): a shared (restricted-open-shell) molecular orbital set is constructed from the
!!     spin-averaged converged Hamiltonian; alpha and beta share spatial orbitals and the
!!     spin-resolved Fock matrices F^alpha/F^beta enter the orbital-Hamiltonian part in full. This is
!!     the basis required for genuine MRSF and gives an exact spin-flip <S^2>.
!!   * Unrestricted: the native collinear (UHF-like) alpha/beta orbitals are used; F^alpha/F^beta are
!!     diagonal (the canonical eigenvalues) and the SOMO pair is identified by index.
!!
!! With MixedReference = Yes the mixed-reference spin-adaptation (MRSF-TDDFT) is applied. This is NOT
!! a generic Slater-Condon / CIS reduction of the single-reference spin-flip matrix: MRSF couples the
!! response configurations originating from BOTH the M_S = +1 and M_S = -1 components of the mixed
!! reference, with the specific sign conventions and block-coupling matrix elements of Lee et al.,
!! J. Chem. Phys. 149, 104101 (2018) and its supplementary material (referred to below as SI). The
!! configurations are classified into four types (manuscript Fig. 3): I = SOMO -> SOMO (O -> O),
!! II = closed -> SOMO (C -> O), III = SOMO -> virtual (O -> V), IV = closed -> virtual (C -> V).
!!
!! Type I (the SOMO pair O1, O2). The two SOMO orbitals carry four flip determinants, labelled (SI
!! Fig. S2) G (both electrons in O1), D (both in O2), L and R (one electron in each SOMO). In the
!! present alpha-occ -> beta-vir bookkeeping O1->O1 and O2->O2 are the open (L, R) determinants and
!! O2->O1, O1->O2 are the closed (G, D) determinants. Symmetrizing over the M_S = +1 / -1 references
!! with the sign settings of SI p.10 (B_GL = -B_GR, B_DL = -B_DR, B_LL = B_RR) splits them into
!!   triplet:  (L + R)/sqrt(2)                       [one configuration, <S^2> = 2],
!!   singlet:  G, D, (L - R)/sqrt(2)                 [three configurations, <S^2> = 0].
!! In TDA, and because the references are spin partners (A_{X,+1} = A_{X,-1}), this symmetrization is
!! realised here by the congruence A_MRSF = T^T A_SF T with the corresponding +-1/sqrt(2) columns.
!!
!! Types II, III (single-SOMO, C -> O and O -> V). These are spin contaminated (<S^2> ~ 1) unless the
!! spin-complete treatment is requested (SpinComplete = Yes, ROHF reference): each M_S = +1
!! configuration |+> is paired with its M_S = -1 partner |-> (SI Eq. S8.1, the U' transform). The
!! pairing introduces the spin-pairing coupling of manuscript Eq. (3.11) / Fig. 5,
!!   C_pq = c_H <Psi^{M_S=+1} | H | Psi^{M_S=-1}>,   c_H = exact-exchange fraction,
!! a cross-reference exchange matrix element. By Slater-Condon the two partner determinants differ in
!! two spin-orbitals (the hole and the spectator SOMO), so C_pq = -c_H K_{spectator,hole} < 0, with
!! K the long-range (LC) exchange integral. In the {|+>, |->} basis (SI Eq. S8.9/S8.10)
!!   A = [[E, C], [C, E]],   S^2 = [[1, 1], [1, 1]],   C = -c_H K,
!! whose eigenvectors (|+> +- |->)/sqrt(2) are the pure triplet (E + C = E - c_H K, S^2 = 2, the
!! symmetric combination) and singlet (E - C = E + c_H K, S^2 = 0). States are selected by <S^2> for
!! the requested multiplicity. This is the coupling that lifts the residual multiplet degeneracy of
!! single-reference SF-TDDFT (confirmed for systems such as CH2, where most states reach
!! <S^2> = 0 / 2 to ~1e-3); residual contamination remains where the exchange K is small.
!!
!! Type IV (closed -> virtual, C -> V) is NOT yet spin-complete - in agreement with the manuscript,
!! only one of the five missing configurations (the triplet partner) is recoverable at the
!! singles/TDA level; the pure C -> V singlet requires the doubly-spin-flipped determinant (a double
!! excitation) and is absent here (as in single-reference SF-TDDFT). These configurations carry
!! <S^2> ~ 2, are kept on the triplet side, and are excluded from singlet spectra by the multiplicity
!! selection. Their contribution to the low-lying states is expected to be small (manuscript Sec.
!! III A) and they are not augmented.
!!
!! <S^2> is evaluated exactly in the shared-orbital (ROHF) basis (SI Eq. S7.11). For the plain
!! spin-flip and the SOMO-pair-only MRSF,
!!   <S^2> = (X_{O1->O1} + X_{O2->O2})^2 + sum_ia X_ia^2 ([a virtual] + [i closed]);
!! for the spin-complete MRSF, from the augmented S^2 matrix above.
!!
!! Note: currently restricted to the serial (non-MPI) build and integer occupations.
module dftbp_timedep_linrespsf
  use dftbp_common_accuracy, only : dp, elecTolMax
  use dftbp_common_constants, only : Hartree__eV, cExchange
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_file, only : TFileDescr, openFile, closeFile
  use dftbp_common_globalenv, only : stdOut
  use dftbp_dftb_hybridxc, only : THybridXcFunc
  use dftbp_io_message, only : error
  use dftbp_io_taggedoutput, only : TTaggedWriter, tagLabels
  use dftbp_math_blasroutines, only : symm, gemm
  use dftbp_math_eigensolver, only : heev, hegv
  use dftbp_timedep_linresptypes, only : TLinResp
  use dftbp_timedep_transcharges, only : transq
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_densedescr, only : TDenseDescr
  implicit none

  private
  public :: LinRespSF_calcExcitations

  !> Output file for spin-flip excitation energies
  character(*), parameter :: sfExcitationsOut = "SF.DAT"

  !> Prefactor of the on-site spin (W) transverse coupling in the spin-flip kernel. The value 2
  !! cancels the local-spin part of the SOMO splitting, restoring the reference-recovery condition
  !! (the M_S=0 component of the reference triplet at zero excitation energy).
  real(dp), parameter :: sfWfactor = 2.0_dp

contains

  !> Calculate collinear spin-flip / mixed-reference spin-flip excitation energies (TDA, LC-DFTB).
  subroutine LinRespSF_calcExcitations(this, env, denseDesc, grndEigVecs, grndEigVal, SSqrReal,&
      & filling, species0, coords0, orb, hybridXc, fdTagged, taggedWriter, excEnergy,&
      & allExcEnergies)

    !> Linear response settings
    type(TLinResp), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Indexing array for dense H and S
    type(TDenseDescr), intent(in) :: denseDesc

    !> Ground state eigenvectors (nOrb, nOrb, nSpin)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Ground state eigenvalues (nOrb, nSpin)
    real(dp), intent(in) :: grndEigVal(:,:)

    !> Square overlap matrix (nOrb, nOrb)
    real(dp), intent(in) :: SSqrReal(:,:)

    !> Ground state occupations (nOrb, nSpin)
    real(dp), intent(in) :: filling(:,:)

    !> Chemical species of the atoms
    integer, intent(in) :: species0(:)

    !> Central-cell atomic coordinates (3, nAtom), for transition dipoles
    real(dp), intent(in) :: coords0(:,:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Range-separated data (provides long-range exchange gamma)
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> File id for tagged output
    type(TFileDescr), intent(in) :: fdTagged

    !> Tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> Energy of state of interest (0 if all states requested)
    real(dp), intent(out) :: excEnergy

    !> Energies of all solved states
    real(dp), intent(inout), allocatable :: allExcEnergies(:)

    integer :: nOrb, nSpin, nAtom, nOccA, nOccB, nVirA, nVirB, nSF, nMat
    integer :: ii, jj, aa, bb, iT, nState, iState
    real(dp), allocatable :: ovrXev(:,:,:), shVecs(:,:,:), lrGamma(:,:)
    real(dp), allocatable :: qOO(:,:,:), qVV(:,:,:)
    real(dp), allocatable :: faOcc(:,:), fbVir(:,:), faFull(:,:), fbFull(:,:)
    real(dp), allocatable :: aMat(:,:), aMrsf(:,:), tMat(:,:), eval(:)
    real(dp), allocatable :: s2(:), xExp(:), wAtom(:), oscStr(:), xExpMat(:,:)
    integer, allocatable :: getIA(:,:), labIA(:,:), domIA(:,:)
    integer :: it1, it2
    real(dp) :: recoveryEnergy

  #:if WITH_SCALAPACK
    call error("Spin-flip linear response is not yet implemented for MPI/ScaLAPACK builds.")
  #:endif

    nSpin = size(grndEigVal, dim=2)
    if (nSpin /= 2) then
      call error("Spin-flip excitations require a collinear spin-polarised reference (nSpin == 2).")
    end if
    if (.not. allocated(hybridXc)) then
      call error("Spin-flip excitations require a range-separated (LC-)DFTB Hamiltonian to provide&
          & the long-range exchange kernel.")
    end if

    nOrb = orb%nOrb
    nAtom = this%nAtom

    ! Occupied/virtual counts per spin channel (integer occupations assumed)
    call countOccVir(filling, nOccA, nVirA, 1)
    call countOccVir(filling, nOccB, nVirB, 2)

    ! Spin-flip manifold: alpha-occupied -> beta-virtual
    nSF = nOccA * nVirB
    if (nSF < 1) then
      call error("No spin-flip (alpha-occupied -> beta-virtual) transitions available.")
    end if

    ! Long-range exchange gamma in AO (atom) basis
    allocate(lrGamma(nAtom, nAtom))
    call hybridXc%getCamGammaCluster(lrGamma)

    ! Enumerate transitions
    allocate(getIA(nSF, 2))
    iT = 0
    do ii = 1, nOccA
      do aa = nOccB + 1, nOrb
        iT = iT + 1
        getIA(iT, :) = [ii, aa]
      end do
    end do

    allocate(faOcc(nOccA, nOccA), fbVir(nVirB, nVirB))
    allocate(qOO(nAtom, nOccA, nOccA), qVV(nAtom, nVirB, nVirB))

    if (this%tRohfRef) then
      ! Restricted-open-shell (shared MO) reference: spin-resolved Fock in a single MO set
      allocate(shVecs(nOrb, nOrb, 1), ovrXev(nOrb, nOrb, 1))
      allocate(faFull(nOrb, nOrb), fbFull(nOrb, nOrb))
      call buildRohfReference(grndEigVecs, grndEigVal, SSqrReal, shVecs, faFull, fbFull)
      call symm(ovrXev(:,:,1), "L", SSqrReal, shVecs(:,:,1))
      ! occupied/virtual Fock blocks
      faOcc(:,:) = faFull(1:nOccA, 1:nOccA)
      fbVir(:,:) = fbFull(nOccB+1:nOrb, nOccB+1:nOrb)
      ! shared-orbital transition charges (alpha occ-occ and beta vir-vir use the same MO set)
      do ii = 1, nOccA
        do jj = 1, ii
          qOO(:, ii, jj) = transq(ii, jj, env, denseDesc, .true., ovrXev, shVecs)
          qOO(:, jj, ii) = qOO(:, ii, jj)
        end do
      end do
      do aa = 1, nVirB
        do bb = 1, aa
          qVV(:, aa, bb) = transq(nOccB+aa, nOccB+bb, env, denseDesc, .true., ovrXev, shVecs)
          qVV(:, bb, aa) = qVV(:, aa, bb)
        end do
      end do
    else
      ! Native collinear (UHF-like) reference: diagonal Fock = canonical eigenvalues
      allocate(ovrXev(nOrb, nOrb, nSpin))
      do ii = 1, nSpin
        call symm(ovrXev(:,:,ii), "L", SSqrReal, grndEigVecs(:,:,ii))
      end do
      faOcc(:,:) = 0.0_dp
      do ii = 1, nOccA
        faOcc(ii, ii) = grndEigVal(ii, 1)
      end do
      fbVir(:,:) = 0.0_dp
      do aa = 1, nVirB
        fbVir(aa, aa) = grndEigVal(nOccB+aa, 2)
      end do
      do ii = 1, nOccA
        do jj = 1, ii
          qOO(:, ii, jj) = transq(ii, jj, env, denseDesc, .true., ovrXev, grndEigVecs)
          qOO(:, jj, ii) = qOO(:, ii, jj)
        end do
      end do
      do aa = 1, nVirB
        do bb = 1, aa
          qVV(:, aa, bb) = transq(nOccB+aa, nOccB+bb, env, denseDesc, .false., ovrXev, grndEigVecs)
          qVV(:, bb, aa) = qVV(:, aa, bb)
        end do
      end do
    end if

    ! Assemble the spin-flip TDA A-matrix
    call buildSpinFlipA(getIA, nOccB, nVirB, faOcc, fbVir, qOO, qVV, lrGamma, aMat)

    ! Add the on-site spin (W) transverse coupling: the local part of the spin-flip kernel that
    ! cancels the local-spin contribution to the SOMO splitting (restores reference recovery).
    if (allocated(this%spinW)) then
      allocate(wAtom(nAtom))
      do ii = 1, nAtom
        wAtom(ii) = this%spinW(species0(ii))
      end do
      if (this%tRohfRef) then
        call addSpinFlipW(aMat, getIA, denseDesc, shVecs(:,:,1), ovrXev(:,:,1), shVecs(:,:,1),&
            & ovrXev(:,:,1), wAtom, sfWfactor)
      else
        call addSpinFlipW(aMat, getIA, denseDesc, grndEigVecs(:,:,1), ovrXev(:,:,1),&
            & grndEigVecs(:,:,2), ovrXev(:,:,2), wAtom, sfWfactor)
      end if
    end if

    ! Reference-recovery diagnostic: energy of the M_S=0 component of the reference triplet,
    ! i.e. the spin-adapted SOMO-pair triplet (O1->O1 + O2->O2)/sqrt(2). Exact spin symmetry
    ! requires this to be ~0; the deviation measures the residual spin-flip kernel error.
    recoveryEnergy = 0.0_dp
    if (nOccA == nOccB + 2) then
      it1 = (nOccB + 1 - 1) * nVirB + (nOccB + 1 - nOccB)
      it2 = (nOccB + 2 - 1) * nVirB + (nOccB + 2 - nOccB)
      recoveryEnergy = 0.5_dp * (aMat(it1, it1) + aMat(it2, it2) + 2.0_dp * aMat(it1, it2))
    end if

    if (this%tMixedRef) then
      ! Mixed-reference (MRSF) spin-adaptation of the SOMO pair: A_MRSF = T^T A_SF T
      call mrsfReduce(aMat, getIA, nOccA, nOccB, nVirB, this%sfMultiplicity, aMrsf, labIA, tMat)
      if (this%tRohfRef .and. this%tSpinComplete) then
        ! Spin-complete treatment: augment the single-SOMO blocks with their M_S=-1 partner
        ! configurations so that closed->SOMO and SOMO->virtual states become spin pure.
        call mrsfSpinComplete(aMrsf, labIA, nOccA, nOccB, this%sfMultiplicity, this%nExc, env,&
            & denseDesc, ovrXev, shVecs, lrGamma, eval, s2, domIA, nState)
        ! Oscillator strengths not yet available in the augmented (spin-complete) basis
        allocate(oscStr(nState))
        oscStr(:) = 0.0_dp
      else
        ! Unrestricted reference: SOMO-pair adaptation only, within-manifold <S^2>
        nMat = size(aMrsf, dim=1)
        allocate(eval(nMat))
        call heev(aMrsf, eval, "U", "V")
        nState = min(this%nExc, nMat)
        allocate(s2(nState), domIA(nState, 2), xExp(nSF), xExpMat(nSF, nState))
        do iState = 1, nState
          xExp(:) = matmul(tMat, aMrsf(:, iState))
          xExpMat(:, iState) = xExp
          s2(iState) = sfSpinSquare(xExp, getIA, nOccA, nOccB)
          domIA(iState, :) = labIA(maxloc(aMrsf(:, iState)**2, dim=1), :)
        end do
        allocate(oscStr(nState))
        call sfOscillator(xExpMat, eval, qOO, qVV, coords0, nOccA, nOccB, nVirB, oscStr)
      end if
      call writeSFResults(eval, domIA, nState, s2, oscStr, fdTagged, taggedWriter,&
          & this%sfMultiplicity, recoveryEnergy)
    else
      ! Plain (spin-contaminated) spin-flip TDDFT
      allocate(eval(nSF))
      call heev(aMat, eval, "U", "V")
      nState = min(this%nExc, nSF)
      allocate(s2(nState), domIA(nState, 2), oscStr(nState))
      do iState = 1, nState
        s2(iState) = sfSpinSquare(aMat(:, iState), getIA, nOccA, nOccB)
        domIA(iState, :) = getIA(maxloc(aMat(:, iState)**2, dim=1), :)
      end do
      call sfOscillator(aMat(:, 1:nState), eval, qOO, qVV, coords0, nOccA, nOccB, nVirB, oscStr)
      call writeSFResults(eval, domIA, nState, s2, oscStr, fdTagged, taggedWriter, 0,&
          & recoveryEnergy)
    end if

    allocate(allExcEnergies(nState))
    allExcEnergies(:) = eval(1:nState)
    if (this%nStat > 0 .and. this%nStat <= nState) then
      excEnergy = eval(this%nStat)
    else
      excEnergy = 0.0_dp
    end if

  end subroutine LinRespSF_calcExcitations


  !> Build the symmetric spin-flip TDA A-matrix in the expanded (alpha-occ -> beta-vir) basis.
  subroutine buildSpinFlipA(getIA, nOccB, nVirB, faOcc, fbVir, qOO, qVV, lrGamma, aMat)

    !> Transition index map [i, a] (orbital numbers)
    integer, intent(in) :: getIA(:,:)

    !> Number of beta-occupied orbitals and beta-virtuals
    integer, intent(in) :: nOccB, nVirB

    !> Alpha occ-occ and beta vir-vir Fock blocks (in the working MO basis)
    real(dp), intent(in) :: faOcc(:,:), fbVir(:,:)

    !> Same-spin occ-occ and vir-vir transition charges
    real(dp), intent(in) :: qOO(:,:,:), qVV(:,:,:)

    !> Long-range exchange gamma (nAtom, nAtom)
    real(dp), intent(in) :: lrGamma(:,:)

    !> Resulting A-matrix
    real(dp), allocatable, intent(out) :: aMat(:,:)

    integer :: nSF, iT, jT, ii, jj, aa, bb, ap, bp
    real(dp), allocatable :: gqVV(:,:,:)
    integer :: nAtom

    nSF = size(getIA, dim=1)
    nAtom = size(lrGamma, dim=1)

    ! Pre-contract gamma with the vir-vir charges
    allocate(gqVV(nAtom, nVirB, nVirB))
    do aa = 1, nVirB
      do bb = 1, nVirB
        gqVV(:, aa, bb) = matmul(lrGamma, qVV(:, aa, bb))
      end do
    end do

    allocate(aMat(nSF, nSF))
    aMat(:,:) = 0.0_dp
    do iT = 1, nSF
      ii = getIA(iT, 1)
      aa = getIA(iT, 2)
      ap = aa - nOccB
      do jT = 1, iT
        jj = getIA(jT, 1)
        bb = getIA(jT, 2)
        bp = bb - nOccB
        ! exchange-only coupling: - q^{ij} . gamma^LR . q^{ab}
        aMat(iT, jT) = -cExchange * dot_product(qOO(:, ii, jj), gqVV(:, ap, bp))
        ! orbital-Hamiltonian part: delta_ij F^beta_ab - delta_ab F^alpha_ij
        if (ii == jj) aMat(iT, jT) = aMat(iT, jT) + fbVir(ap, bp)
        if (aa == bb) aMat(iT, jT) = aMat(iT, jT) - faOcc(ii, jj)
        aMat(jT, iT) = aMat(iT, jT)
      end do
    end do

  end subroutine buildSpinFlipA


  !> Add the on-site spin (W) transverse coupling to the spin-flip A-matrix:
  !!   A_{ia,jb} += cW * sum_A W_A q^{ia}_A q^{jb}_A,
  !! where q^{ia}_A is the (cross-spin) spin-flip transition charge on atom A. This is the local
  !! transverse part of the spin-flip exchange-correlation kernel (the analogue of the magnetisation
  !! coupling used for triplet excitations); together with the long-range exchange it cancels the
  !! SOMO spin-splitting, so the M_S=0 component of the reference triplet returns at zero energy.
  subroutine addSpinFlipW(aMat, getIA, denseDesc, occVec, occOvr, virVec, virOvr, wAtom, cW)

    !> Spin-flip A-matrix (modified in place)
    real(dp), intent(inout) :: aMat(:,:)

    !> Transition index map [i, a]
    integer, intent(in) :: getIA(:,:)

    !> Dense matrix descriptor (atom -> orbital ranges)
    type(TDenseDescr), intent(in) :: denseDesc

    !> Occupied-channel MO coefficients and S times them
    real(dp), intent(in) :: occVec(:,:), occOvr(:,:)

    !> Virtual-channel MO coefficients and S times them
    real(dp), intent(in) :: virVec(:,:), virOvr(:,:)

    !> Per-atom spin constant W
    real(dp), intent(in) :: wAtom(:)

    !> Prefactor
    real(dp), intent(in) :: cW

    integer :: nSF, nAtom, iT, jT, ii, aa, kk, m1, m2
    real(dp), allocatable :: qIA(:,:)
    real(dp) :: s

    nSF = size(getIA, dim=1)
    nAtom = size(wAtom)
    allocate(qIA(nAtom, nSF))
    do iT = 1, nSF
      ii = getIA(iT, 1)
      aa = getIA(iT, 2)
      do kk = 1, nAtom
        m1 = denseDesc%iAtomStart(kk)
        m2 = denseDesc%iAtomStart(kk + 1) - 1
        qIA(kk, iT) = 0.5_dp * sum(occVec(m1:m2, ii) * virOvr(m1:m2, aa)&
            & + virVec(m1:m2, aa) * occOvr(m1:m2, ii))
      end do
    end do

    do iT = 1, nSF
      do jT = 1, iT
        s = cW * sum(wAtom(:) * qIA(:, iT) * qIA(:, jT))
        aMat(iT, jT) = aMat(iT, jT) + s
        if (iT /= jT) aMat(jT, iT) = aMat(jT, iT) + s
      end do
    end do

  end subroutine addSpinFlipW


  !> Oscillator strengths for spin-flip absorption from the lowest state (state 1) to each state.
  !!
  !! The transition dipole between two spin-flip states m, n (with expanded amplitudes X) is
  !!   d_{mn} = sum_{i,a,b} X^m_{ia} X^n_{ib} mu^beta_{ab} - sum_{i,j,a} X^m_{ia} X^n_{ja} mu^alpha_{ij},
  !! with the DFTB monopole MO dipole mu_{pq} = sum_A q^{pq}_A R_A. The oscillator strength is
  !! f_n = (2/3) (omega_n - omega_0) |d_{0n}|^2. The construction is translation invariant
  !! (sum_A q^{pq}_A = delta_pq).
  subroutine sfOscillator(xExp, eval, qOO, qVV, coords0, nOccA, nOccB, nVirB, osc)

    !> Expanded spin-flip amplitudes (nSF, nState), columns are states (ascending energy)
    real(dp), intent(in) :: xExp(:,:)

    !> Excitation energies (ascending)
    real(dp), intent(in) :: eval(:)

    !> Alpha occ-occ and beta vir-vir transition charges
    real(dp), intent(in) :: qOO(:,:,:), qVV(:,:,:)

    !> Central-cell atomic coordinates (3, nAtom)
    real(dp), intent(in) :: coords0(:,:)

    !> Orbital-space sizes
    integer, intent(in) :: nOccA, nOccB, nVirB

    !> Oscillator strengths on exit
    real(dp), intent(out) :: osc(:)

    integer :: nState, ii, jj, ap, bp, xyz, nn
    real(dp), allocatable :: dipOO(:,:,:), dipVV(:,:,:), x0(:,:), xn(:,:)
    real(dp) :: dvec(3), dw

    nState = size(xExp, dim=2)

    ! monopole MO dipole integrals mu_{pq}[xyz] = sum_A q^{pq}_A R_A[xyz]
    allocate(dipOO(3, nOccA, nOccA), dipVV(3, nVirB, nVirB))
    do ii = 1, nOccA
      do jj = 1, nOccA
        do xyz = 1, 3
          dipOO(xyz, ii, jj) = sum(coords0(xyz, :) * qOO(:, ii, jj))
        end do
      end do
    end do
    do ap = 1, nVirB
      do bp = 1, nVirB
        do xyz = 1, 3
          dipVV(xyz, ap, bp) = sum(coords0(xyz, :) * qVV(:, ap, bp))
        end do
      end do
    end do

    allocate(x0(nOccA, nVirB), xn(nOccA, nVirB))
    call reshapeAmp(xExp(:, 1), nOccA, nVirB, x0)

    osc(1) = 0.0_dp
    do nn = 2, nState
      call reshapeAmp(xExp(:, nn), nOccA, nVirB, xn)
      dvec(:) = 0.0_dp
      do xyz = 1, 3
        do ii = 1, nOccA
          do ap = 1, nVirB
            do bp = 1, nVirB
              dvec(xyz) = dvec(xyz) + x0(ii, ap) * xn(ii, bp) * dipVV(xyz, ap, bp)
            end do
          end do
        end do
        do ap = 1, nVirB
          do ii = 1, nOccA
            do jj = 1, nOccA
              dvec(xyz) = dvec(xyz) - x0(ii, ap) * xn(jj, ap) * dipOO(xyz, ii, jj)
            end do
          end do
        end do
      end do
      dw = eval(nn) - eval(1)
      osc(nn) = (2.0_dp / 3.0_dp) * dw * dot_product(dvec, dvec)
    end do

  end subroutine sfOscillator


  !> Reshape an expanded spin-flip amplitude vector to (alpha-occ, beta-vir) form.
  subroutine reshapeAmp(xvec, nOccA, nVirB, xmat)

    !> Expanded amplitude vector (nSF), ordering iT = (i-1)*nVirB + (a - nOccB)
    real(dp), intent(in) :: xvec(:)

    !> Orbital-space sizes
    integer, intent(in) :: nOccA, nVirB

    !> Reshaped amplitudes (nOccA, nVirB)
    real(dp), intent(out) :: xmat(:,:)

    integer :: ii, ap

    do ii = 1, nOccA
      do ap = 1, nVirB
        xmat(ii, ap) = xvec((ii - 1) * nVirB + ap)
      end do
    end do

  end subroutine reshapeAmp


  !> Construct a restricted-open-shell (shared) MO set and the spin-resolved Fock matrices in it.
  !!
  !! The shared orbitals diagonalise the spin-averaged converged Hamiltonian H = (H^a + H^b)/2,
  !! reconstructed from the (S-orthonormal) collinear orbitals as H^s = S C^s E^s (C^s)^T S. The
  !! spin-resolved Fock matrices are then F^s = O^s diag(E^s) (O^s)^T with O^s = C^T S C^s.
  subroutine buildRohfReference(eigVecs, eigVal, SSqr, shVecs, faFull, fbFull)

    !> Collinear eigenvectors (nOrb, nOrb, 2)
    real(dp), intent(in) :: eigVecs(:,:,:)

    !> Collinear eigenvalues (nOrb, 2)
    real(dp), intent(in) :: eigVal(:,:)

    !> Overlap matrix (nOrb, nOrb)
    real(dp), intent(in) :: SSqr(:,:)

    !> Shared MO coefficients on exit (nOrb, nOrb, 1)
    real(dp), intent(out) :: shVecs(:,:,:)

    !> Alpha and beta Fock matrices in the shared MO basis (nOrb, nOrb)
    real(dp), intent(out) :: faFull(:,:), fbFull(:,:)

    integer :: nOrb, spin, pp
    real(dp), allocatable :: pMat(:,:), hCharge(:,:), sCopy(:,:), tmp(:,:)
    real(dp), allocatable :: wScaled(:,:), oMat(:,:), sC(:,:), eR(:)

    nOrb = size(SSqr, dim=1)
    allocate(pMat(nOrb,nOrb), hCharge(nOrb,nOrb), sCopy(nOrb,nOrb), tmp(nOrb,nOrb))
    allocate(wScaled(nOrb,nOrb), oMat(nOrb,nOrb), sC(nOrb,nOrb), eR(nOrb))

    ! P = 0.5 (C^a E^a C^a^T + C^b E^b C^b^T)
    pMat(:,:) = 0.0_dp
    do spin = 1, 2
      do pp = 1, nOrb
        wScaled(:, pp) = eigVecs(:, pp, spin) * eigVal(pp, spin)
      end do
      call gemm(pMat, wScaled, eigVecs(:,:,spin), alpha=0.5_dp, beta=1.0_dp, transB="T")
    end do

    ! H_charge = S P S
    call gemm(tmp, SSqr, pMat)
    call gemm(hCharge, tmp, SSqr)

    ! Solve H_charge C = S C E  ->  shared MOs (S-orthonormal, ascending energy)
    sCopy(:,:) = SSqr
    call hegv(hCharge, sCopy, eR, "U", "V")
    shVecs(:,:,1) = hCharge

    ! F^s = O^s diag(E^s) (O^s)^T,   O^s = C^T S C^s
    do spin = 1, 2
      call gemm(sC, SSqr, eigVecs(:,:,spin))
      call gemm(oMat, shVecs(:,:,1), sC, transA="T")
      do pp = 1, nOrb
        wScaled(:, pp) = oMat(:, pp) * eigVal(pp, spin)
      end do
      if (spin == 1) then
        call gemm(faFull, wScaled, oMat, transB="T")
      else
        call gemm(fbFull, wScaled, oMat, transB="T")
      end if
    end do

  end subroutine buildRohfReference


  !> Spin-square <S^2> of a spin-flip state from its expanded amplitude vector.
  !!
  !! Exact in a shared-orbital (ROHF) basis:
  !!   <S^2> = (X_{O1->O1} + X_{O2->O2})^2 + sum_ia X_ia^2 ( [a virtual] + [i closed] ).
  function sfSpinSquare(xExp, getIA, nOccA, nOccB) result(s2)

    !> Expanded spin-flip amplitudes (normalised), one per transition
    real(dp), intent(in) :: xExp(:)

    !> Transition index map [i, a] (orbital numbers)
    integer, intent(in) :: getIA(:,:)

    !> Number of alpha-occupied and beta-occupied orbitals
    integer, intent(in) :: nOccA, nOccB

    real(dp) :: s2

    integer :: iT, ii, aa, sia
    real(dp) :: openSum

    s2 = 0.0_dp
    openSum = 0.0_dp
    do iT = 1, size(xExp)
      ii = getIA(iT, 1)
      aa = getIA(iT, 2)
      sia = 0
      if (aa > nOccA) sia = sia + 1   ! a is a true virtual (not a SOMO)
      if (ii <= nOccB) sia = sia + 1  ! i is a closed (doubly occupied) orbital
      s2 = s2 + xExp(iT)**2 * real(sia, dp)
      ! i == a (same orbital index) occurs only for the SOMO-diagonal flips O1->O1, O2->O2
      if (ii == aa) openSum = openSum + xExp(iT)
    end do
    s2 = s2 + openSum**2

  end function sfSpinSquare


  !> Reduce the spin-flip A-matrix to the mixed-reference (MRSF) spin-adapted A-matrix.
  !!
  !! Performs the type I (SOMO -> SOMO) symmetrization of MRSF (SI Sec. 5 / Eqs S8.7-S8.8). The two
  !! SOMOs O1, O2 carry four flip determinants (SI Fig. S2): in the alpha-occ -> beta-vir bookkeeping
  !! used here the open determinants L, R are O1->O1 (ijlr1) and O2->O2 (ijlr2), and the closed
  !! determinants G, D are O2->O1 (ijg) and O1->O2 (ijd). Applying the sign settings of SI p.10
  !! (B_GL = -B_GR, B_DL = -B_DR, B_LL = B_RR) the symmetrized configurations are
  !!   triplet:  (L + R)/sqrt(2) = (ijlr1 + ijlr2)/sqrt(2)            [keep; remove G, D],
  !!   singlet:  G, D, (L - R)/sqrt(2) = ijg, ijd, (ijlr1 - ijlr2)/sqrt(2)  [keep; one open combo].
  !! Because the M_S = +1 and M_S = -1 references are spin partners (A_{X,+1} = A_{X,-1}), this
  !! symmetrization is realised in TDA by the congruence A_MRSF = T^T A_SF T, where the column of T
  !! for the retained open combination carries +-1/sqrt(2) on (ijlr1, ijlr2). The single-reference
  !! SF matrix A_SF therefore supplies the diagonal B-blocks; the cross-reference spin-pairing
  !! coupling C (manuscript Eq. 3.11) enters only for types II/III/IV and is added in
  !! mrsfSpinComplete.
  subroutine mrsfReduce(aSF, getIA, nOccA, nOccB, nVirB, mult, aMrsf, labIA, tMat)

    !> Spin-flip A-matrix in the expanded basis
    real(dp), intent(in) :: aSF(:,:)

    !> Transition index map for the expanded basis
    integer, intent(in) :: getIA(:,:)

    !> Number of alpha-occupied / beta-occupied orbitals and beta-virtuals
    integer, intent(in) :: nOccA, nOccB, nVirB

    !> Target multiplicity (1 = singlet, 3 = triplet)
    integer, intent(in) :: mult

    !> Reduced MRSF A-matrix
    real(dp), allocatable, intent(out) :: aMrsf(:,:)

    !> Representative transition labels for the reduced configurations
    integer, allocatable, intent(out) :: labIA(:,:)

    !> Spin-adaptation transformation (expanded x compressed)
    real(dp), allocatable, intent(out) :: tMat(:,:)

    integer :: nSF, o1, o2, ijlr1, ijlr2, ijg, ijd, nRem, nC, ee, col
    logical, allocatable :: active(:)
    real(dp) :: isq2, signLr2

    if (nOccA /= nOccB + 2) then
      call error("MRSF requires a high-spin triplet reference with exactly two unpaired electrons&
          & (SpinPolarisation = Colinear { UnpairedElectrons = 2 }).")
    end if
    if (mult /= 1 .and. mult /= 3) then
      call error("MRSF Multiplicity must be 1 (singlet) or 3 (triplet).")
    end if

    nSF = size(aSF, dim=1)
    isq2 = 1.0_dp / sqrt(2.0_dp)
    o1 = nOccB + 1
    o2 = nOccB + 2

    ! SI Fig. S2 determinants in alpha-occ -> beta-vir bookkeeping:
    !   ijlr1 = O1->O1, ijlr2 = O2->O2  : open determinants L, R (one electron in each SOMO)
    !   ijg   = O2->O1, ijd   = O1->O2  : closed determinants G, D (both electrons in one SOMO)
    ijlr1 = (o1 - 1) * nVirB + (o1 - nOccB)
    ijlr2 = (o2 - 1) * nVirB + (o2 - nOccB)
    ijg = (o2 - 1) * nVirB + (o1 - nOccB)
    ijd = (o1 - 1) * nVirB + (o2 - nOccB)

    allocate(active(nSF))
    active(:) = .true.
    if (mult == 1) then
      ! singlet type I = { G (ijg), D (ijd), (L - R)/sqrt(2) }: drop ijlr2, keep ijlr1 as (L - R)
      active(ijlr2) = .false.
      nRem = 1
      signLr2 = -isq2
    else
      ! triplet type I = { (L + R)/sqrt(2) }: drop the closed G, D and keep ijlr1 as (L + R)
      active(ijlr2) = .false.
      active(ijg) = .false.
      active(ijd) = .false.
      nRem = 3
      signLr2 = isq2
    end if
    nC = nSF - nRem

    allocate(tMat(nSF, nC))
    tMat(:,:) = 0.0_dp
    allocate(labIA(nC, 2))
    col = 0
    do ee = 1, nSF
      if (.not. active(ee)) cycle
      col = col + 1
      if (ee == ijlr1) then
        tMat(ijlr1, col) = isq2
        tMat(ijlr2, col) = signLr2
        labIA(col, :) = [o1, o2]
      else
        tMat(ee, col) = 1.0_dp
        labIA(col, :) = getIA(ee, :)
      end if
    end do

    allocate(aMrsf(nC, nC))
    aMrsf(:,:) = matmul(transpose(tMat), matmul(aSF, tMat))
    aMrsf(:,:) = 0.5_dp * (aMrsf + transpose(aMrsf))

  end subroutine mrsfReduce


  !> Spin-complete MRSF: augment the type II (closed->SOMO) and type III (SOMO->virtual) blocks of
  !! the reduced MRSF matrix with their M_S = -1 partner configurations (SI Eq. S8.1, the U'
  !! transform), so that these states become spin pure.
  !!
  !! Each single-SOMO spin-flip configuration |+> (an alpha->beta excitation of the M_S=+1 reference)
  !! has a partner |-> (the mirror beta->alpha excitation of the M_S=-1 reference). The two are
  !! coupled by the spin-pairing matrix element of manuscript Eq. (3.11) / Fig. 5,
  !!   C = c_H <Psi^{M_S=+1} | H | Psi^{M_S=-1}>,   c_H = exact-exchange fraction (= cExchange here).
  !! The partner determinants differ in two spin-orbitals - the hole (closed orbital for type II, or
  !! the particle virtual for type III) and the spectator SOMO - so by Slater-Condon C reduces to a
  !! single cross-reference exchange integral, C = -c_H K with K = q^{spec,part} . gamma^LR .
  !! q^{spec,part} >= 0 the long-range exchange between the spectator SOMO and the partner orbital.
  !! In the basis {|+>, |->} the response and spin-square blocks are (SI Eq. S8.9/S8.10)
  !!   A = [[E, C], [C, E]],   S^2 = [[1, 1], [1, 1]],   C = -c_H K,
  !! These commute, so the eigenvectors (|+> +- |->)/sqrt(2) are the pure triplet (E + C = E - c_H K,
  !! S^2=2, the symmetric combination) and singlet (E - C = E + c_H K, S^2=0). The full augmented
  !! matrices are diagonalised together; states are selected by <S^2> for the requested multiplicity.
  !!
  !! Type IV (closed->virtual) is left un-augmented: only its triplet partner is recoverable at the
  !! singles/TDA level (manuscript Sec. III A), so these configurations keep <S^2> ~ 2 and are
  !! excluded from the singlet spectrum by the multiplicity selection.
  subroutine mrsfSpinComplete(aMrsf, labIA, nOccA, nOccB, mult, nExc, env, denseDesc, ovrXev,&
      & shVecs, lrGamma, evalOut, s2Out, domOut, nStateOut)

    !> Reduced (SOMO-pair adapted) MRSF A-matrix
    real(dp), intent(in) :: aMrsf(:,:)

    !> Representative transition labels [i, a] of the reduced configurations
    integer, intent(in) :: labIA(:,:)

    !> Number of alpha-occupied and beta-occupied orbitals
    integer, intent(in) :: nOccA, nOccB

    !> Target multiplicity (1 = singlet, 3 = triplet)
    integer, intent(in) :: mult

    !> Number of excited states requested
    integer, intent(in) :: nExc

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Overlap times shared MOs (nOrb, nOrb, 1)
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Shared MO coefficients (nOrb, nOrb, 1)
    real(dp), intent(in) :: shVecs(:,:,:)

    !> Long-range exchange gamma
    real(dp), intent(in) :: lrGamma(:,:)

    !> Excitation energies of the selected states
    real(dp), allocatable, intent(out) :: evalOut(:)

    !> Spin-square of the selected states
    real(dp), allocatable, intent(out) :: s2Out(:)

    !> Dominant transition label of the selected states
    integer, allocatable, intent(out) :: domOut(:,:)

    !> Number of selected states
    integer, intent(out) :: nStateOut

    integer :: nC, nP, nAug, o1, o2, r, kk, ll, pk, ii, aa, nAtom, jj, nSel
    integer, allocatable :: somoCfg(:), specOrb(:), partOrb(:)
    real(dp), allocatable :: aAug(:,:), s2Aug(:,:), qpq(:), eval(:), s2all(:), energy(:), av(:)
    real(dp) :: kExch, targetS2, lambda
    logical, allocatable :: isTarget(:)

    o1 = nOccB + 1
    o2 = nOccB + 2
    nC = size(aMrsf, dim=1)
    nAtom = size(lrGamma, dim=1)
    targetS2 = merge(0.0_dp, 2.0_dp, mult == 1)

    ! Identify the single-SOMO reduced configurations and their spectator / partner orbitals
    allocate(somoCfg(nC), specOrb(nC), partOrb(nC))
    nP = 0
    do r = 1, nC
      ii = labIA(r, 1)
      aa = labIA(r, 2)
      if (ii <= nOccB .and. (aa == o1 .or. aa == o2)) then
        ! closed -> SOMO: spectator is the other SOMO, partner orbital is the hole i
        nP = nP + 1
        somoCfg(nP) = r
        specOrb(nP) = merge(o2, o1, aa == o1)
        partOrb(nP) = ii
      else if ((ii == o1 .or. ii == o2) .and. aa > nOccA) then
        ! SOMO -> virtual: spectator is the other SOMO, partner orbital is the particle a
        nP = nP + 1
        somoCfg(nP) = r
        specOrb(nP) = merge(o2, o1, ii == o1)
        partOrb(nP) = aa
      end if
    end do

    nAug = nC + nP
    allocate(aAug(nAug, nAug), s2Aug(nAug, nAug))
    aAug(:,:) = 0.0_dp
    s2Aug(:,:) = 0.0_dp

    ! Reduced (M_S=+1) block
    aAug(1:nC, 1:nC) = aMrsf
    do r = 1, nC
      ii = labIA(r, 1)
      aa = labIA(r, 2)
      if (ii == o1 .and. aa == o2) then
        ! SOMO-pair spin-adapted configuration (already pure)
        s2Aug(r, r) = targetS2
      else
        s2Aug(r, r) = merge(1.0_dp, 0.0_dp, aa > nOccA) + merge(1.0_dp, 0.0_dp, ii <= nOccB)
      end if
    end do

    ! Partner (M_S=-1) block: mirror the single-SOMO sub-block of the reduced matrix
    do kk = 1, nP
      do ll = 1, nP
        aAug(nC + kk, nC + ll) = aMrsf(somoCfg(kk), somoCfg(ll))
      end do
      s2Aug(nC + kk, nC + kk) = 1.0_dp
    end do

    ! Cross (M_S=+1 <-> M_S=-1) coupling: the spin-pairing matrix element C of manuscript Eq. (3.11),
    !   C = c_H <Psi^{+1}|H|Psi^{-1}> = -c_H K,   c_H = cExchange,
    ! with K = q^{spec,part} . gamma^LR . q^{spec,part} the long-range exchange between the spectator
    ! SOMO and the partner (hole/particle) orbital. The accompanying S^2 cross term is 1 (SI S7.11).
    allocate(qpq(nAtom))
    do kk = 1, nP
      r = somoCfg(kk)
      pk = nC + kk
      qpq(:) = transq(specOrb(kk), partOrb(kk), env, denseDesc, .true., ovrXev, shVecs)
      kExch = cExchange * dot_product(qpq, matmul(lrGamma, qpq))
      aAug(r, pk) = -kExch
      aAug(pk, r) = -kExch
      s2Aug(r, pk) = 1.0_dp
      s2Aug(pk, r) = 1.0_dp
    end do

    ! Diagonalise A (with a tiny S^2 tie-breaker to fix the eigenvectors in degenerate subspaces),
    ! then evaluate the exact energy and <S^2> as expectation values.
    lambda = 1.0e-6_dp
    aAug(:,:) = aAug + lambda * s2Aug
    allocate(eval(nAug))
    call heev(aAug, eval, "U", "V")

    allocate(energy(nAug), s2all(nAug), av(nAug))
    do jj = 1, nAug
      av(:) = matmul(s2Aug, aAug(:, jj))
      s2all(jj) = dot_product(aAug(:, jj), av)
      ! E = (eval including tie-breaker) - lambda * <S^2>
      energy(jj) = eval(jj) - lambda * s2all(jj)
    end do

    ! Select states of the requested multiplicity (lowest in energy)
    allocate(isTarget(nAug))
    isTarget(:) = abs(s2all - targetS2) < 1.0_dp
    nSel = min(nExc, count(isTarget))
    nStateOut = nSel
    allocate(evalOut(nSel), s2Out(nSel), domOut(nSel, 2))
    nSel = 0
    do jj = 1, nAug
      if (.not. isTarget(jj)) cycle
      nSel = nSel + 1
      if (nSel > nStateOut) exit
      evalOut(nSel) = energy(jj)
      s2Out(nSel) = s2all(jj)
      ! dominant reduced (M_S=+1) configuration
      r = maxloc(aAug(1:nC, jj)**2, dim=1)
      domOut(nSel, :) = labIA(r, :)
    end do

  end subroutine mrsfSpinComplete


  !> Count occupied and virtual orbitals in a spin channel (integer occupations).
  subroutine countOccVir(filling, nOcc, nVir, iSpin)

    !> Occupations (nOrb, nSpin)
    real(dp), intent(in) :: filling(:,:)

    !> Number of occupied orbitals
    integer, intent(out) :: nOcc

    !> Number of virtual orbitals
    integer, intent(out) :: nVir

    !> Spin channel
    integer, intent(in) :: iSpin

    integer :: ii, nOrb

    nOrb = size(filling, dim=1)
    nOcc = 0
    do ii = 1, nOrb
      if (filling(ii, iSpin) > 1.0_dp - elecTolMax) then
        nOcc = nOcc + 1
      else if (filling(ii, iSpin) > elecTolMax) then
        call error("Spin-flip excitations require integer occupations (fractional occupation&
            & detected).")
      end if
    end do
    nVir = nOrb - nOcc

  end subroutine countOccVir


  !> Write spin-flip excitation energies, <S^2>, oscillator strengths and tagged output.
  subroutine writeSFResults(eval, domIA, nState, s2, osc, fdTagged, taggedWriter, mult,&
      & recoveryEnergy)

    !> Excitation energies (ascending)
    real(dp), intent(in) :: eval(:)

    !> Dominant transition label [i, a] per state
    integer, intent(in) :: domIA(:,:)

    !> Number of states to report
    integer, intent(in) :: nState

    !> Spin-square of each reported state
    real(dp), intent(in) :: s2(:)

    !> Oscillator strength of each reported state (from the lowest state)
    real(dp), intent(in) :: osc(:)

    !> File id for tagged output
    type(TFileDescr), intent(in) :: fdTagged

    !> Tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> Target multiplicity (0 = plain SF, 1 = MRSF singlet, 3 = MRSF triplet)
    integer, intent(in) :: mult

    !> Reference-recovery energy (M_S=0 component of the reference triplet; should be ~0)
    real(dp), intent(in) :: recoveryEnergy

    type(TFileDescr) :: fdSF
    integer :: iState
    character(40) :: methodStr

    select case (mult)
    case (1)
      methodStr = "MRSF singlet (TDA, LC-DFTB)"
    case (3)
      methodStr = "MRSF triplet (TDA, LC-DFTB)"
    case default
      methodStr = "spin-flip (TDA, LC-DFTB)"
    end select

    write(stdOut, "(/,A)") " "//trim(methodStr)//" excitations:"
    write(stdOut, "(2X,A,F12.6,A)") "Reference recovery (should be ~0): ",&
        & recoveryEnergy * Hartree__eV, " eV"
    write(stdOut, "(2X,A6,2X,A14,2X,A12,2X,A8,2X,A10,2X,A12)") "State", "energy (eV)",&
        & "omega (au)", "<S^2>", "f_osc", "dominant i->a"

    call openFile(fdSF, sfExcitationsOut, mode="w")
    write(fdSF%unit, "(A)") "# "//trim(methodStr)//" excitations"
    write(fdSF%unit, "(A,F16.8,A)") "# reference recovery (should be ~0): ",&
        & recoveryEnergy * Hartree__eV, " eV"
    write(fdSF%unit, "(A)") "# state    energy(eV)         omega(au)        <S^2>      f_osc&
        &      dominant"

    do iState = 1, nState
      write(stdOut, "(2X,I6,2X,F14.6,2X,F12.6,2X,F8.4,2X,F10.5,2X,I5,A,I5)") iState,&
          & eval(iState) * Hartree__eV, eval(iState), s2(iState), osc(iState),&
          & domIA(iState, 1), " ->", domIA(iState, 2)
      write(fdSF%unit, "(I6,2X,F16.8,2X,F16.8,2X,F10.5,2X,F12.6,2X,I5,A,I5)") iState,&
          & eval(iState) * Hartree__eV, eval(iState), s2(iState), osc(iState),&
          & domIA(iState, 1), " ->", domIA(iState, 2)
    end do
    call closeFile(fdSF)

    if (fdTagged%isConnected()) then
      call taggedWriter%write(fdTagged%unit, tagLabels%excEgy, eval(1:nState))
      call taggedWriter%write(fdTagged%unit, "exc_spinsquared", s2(1:nState))
      call taggedWriter%write(fdTagged%unit, tagLabels%excOsc, osc(1:nState))
    end if

  end subroutine writeSFResults

end module dftbp_timedep_linrespsf
