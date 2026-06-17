!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Collinear spin-flip linear-response excitations within the Tamm-Dancoff approximation
!! (SF-TDDFT) for range-separated (LC-)DFTB.
!!
!! On a spin-polarised, high-spin (M_S = +1) reference, the spin-flip manifold consists of
!! alpha-occupied -> beta-virtual single excitations. Because the spin-flip transition density
!! is off-diagonal in spin, the Hartree (gamma) coupling vanishes and the off-diagonal coupling is
!! carried entirely by the long-range exact (Hartree-Fock) exchange of LC-DFTB, in line with the
!! "exchange-only" structure of spin-flip / mixed-reference TDDFT.
!!
!! TDA A-matrix:
!!   A_{ia,jb} = delta_ij delta_ab (eps^beta_a - eps^alpha_i)
!!                 - sum_{AB} q^{alpha,ij}_A gamma^LR_{AB} q^{beta,ab}_B
!!
!! This is the first ("SF-TDDFT") stage of a mixed-reference (MRSF) implementation; the
!! mixed-reference 1/sqrt(2) spin-adaptation of the singly-occupied orbitals is layered on top of
!! this manifold in a subsequent step.
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
  use dftbp_math_blasroutines, only : symm
  use dftbp_math_eigensolver, only : heev
  use dftbp_timedep_linresptypes, only : TLinResp
  use dftbp_timedep_transcharges, only : transq
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_densedescr, only : TDenseDescr
  implicit none

  private
  public :: LinRespSF_calcExcitations

  !> Output file for spin-flip excitation energies
  character(*), parameter :: sfExcitationsOut = "SF.DAT"

contains

  !> Calculate collinear spin-flip excitation energies (TDA, LC-DFTB exchange kernel).
  subroutine LinRespSF_calcExcitations(this, env, denseDesc, grndEigVecs, grndEigVal, SSqrReal,&
      & filling, orb, hybridXc, fdTagged, taggedWriter, excEnergy, allExcEnergies)

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

    integer :: nOrb, nSpin, nAtom, nOccA, nOccB, nVirA, nVirB, nSF
    integer :: ii, jj, aa, bb, ap, bp, iT, jT, nState
    real(dp), allocatable :: ovrXev(:,:,:), lrGamma(:,:)
    real(dp), allocatable :: qOO(:,:,:), qVV(:,:,:), gqVV(:,:,:)
    real(dp), allocatable :: aMat(:,:), eval(:), wia(:)
    integer, allocatable :: getIA(:,:)

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
    if (this%nExc > nSF) then
      this%nExc = nSF
    end if

    ! Overlap times eigenvectors: ovrXev(:,:,s) = S c_s
    allocate(ovrXev(nOrb, nOrb, nSpin))
    do ii = 1, nSpin
      call symm(ovrXev(:,:,ii), "L", SSqrReal, grndEigVecs(:,:,ii))
    end do

    ! Long-range exchange gamma in AO (atom) basis
    allocate(lrGamma(nAtom, nAtom))
    call hybridXc%getCamGammaCluster(lrGamma)

    ! Enumerate transitions and single-particle energy differences
    allocate(getIA(nSF, 2))
    allocate(wia(nSF))
    iT = 0
    do ii = 1, nOccA
      do aa = nOccB + 1, nOrb
        iT = iT + 1
        getIA(iT, :) = [ii, aa]
        wia(iT) = grndEigVal(aa, 2) - grndEigVal(ii, 1)
      end do
    end do

    ! Same-spin transition charges:
    !   qOO(:,i,j) alpha occ-occ,  qVV(:,a',b') beta vir-vir (a' = a - nOccB)
    allocate(qOO(nAtom, nOccA, nOccA))
    do ii = 1, nOccA
      do jj = 1, ii
        qOO(:, ii, jj) = transq(ii, jj, env, denseDesc, .true., ovrXev, grndEigVecs)
        qOO(:, jj, ii) = qOO(:, ii, jj)
      end do
    end do
    allocate(qVV(nAtom, nVirB, nVirB))
    do aa = 1, nVirB
      do bb = 1, aa
        qVV(:, aa, bb) = transq(nOccB + aa, nOccB + bb, env, denseDesc, .false., ovrXev, grndEigVecs)
        qVV(:, bb, aa) = qVV(:, aa, bb)
      end do
    end do

    ! Pre-contract gamma with the beta vir-vir charges: gqVV(:,a',b') = gamma^LR . qVV(:,a',b')
    allocate(gqVV(nAtom, nVirB, nVirB))
    do aa = 1, nVirB
      do bb = 1, nVirB
        gqVV(:, aa, bb) = matmul(lrGamma, qVV(:, aa, bb))
      end do
    end do

    ! Build the symmetric TDA spin-flip A-matrix
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
        ! exchange-only off-diagonal coupling: - q^{ij} . gamma^LR . q^{ab}
        aMat(iT, jT) = -cExchange * dot_product(qOO(:, ii, jj), gqVV(:, ap, bp))
        aMat(jT, iT) = aMat(iT, jT)
      end do
      ! single-particle (orbital energy difference) on the diagonal
      aMat(iT, iT) = aMat(iT, iT) + wia(iT)
    end do

    ! Diagonalise (eigenvalues ascending; negative roots are physical for spin-flip)
    allocate(eval(nSF))
    call heev(aMat, eval, "U", "V")

    nState = this%nExc
    allocate(allExcEnergies(nState))
    allExcEnergies(:) = eval(1:nState)

    if (this%nStat > 0) then
      excEnergy = eval(this%nStat)
    else
      excEnergy = 0.0_dp
    end if

    call writeSFResults(eval, aMat, getIA, nState, fdTagged, taggedWriter)

  end subroutine LinRespSF_calcExcitations


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


  !> Write spin-flip excitation energies, dominant transition character and tagged output.
  subroutine writeSFResults(eval, eigVec, getIA, nState, fdTagged, taggedWriter)

    !> Excitation energies (all roots, ascending)
    real(dp), intent(in) :: eval(:)

    !> Excitation eigenvectors (columns)
    real(dp), intent(in) :: eigVec(:,:)

    !> Transition index map
    integer, intent(in) :: getIA(:,:)

    !> Number of states to report
    integer, intent(in) :: nState

    !> File id for tagged output
    type(TFileDescr), intent(in) :: fdTagged

    !> Tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    type(TFileDescr) :: fdSF
    integer :: iState, iT, iMax
    real(dp) :: wMax

    write(stdOut, "(/,A)") " Spin-flip excitations (TDA, LC-DFTB exchange kernel):"
    write(stdOut, "(2X,A6,2X,A14,2X,A14,2X,A20)") "State", "energy (au)", "energy (eV)",&
        & "dominant i(a)->a(b)"

    call openFile(fdSF, sfExcitationsOut, mode="w")
    write(fdSF%unit, "(A)") "# Spin-flip excitations (TDA, LC-DFTB)"
    write(fdSF%unit, "(A)") "# state   energy(eV)        omega(au)     dominant transition"

    do iState = 1, nState
      ! dominant single-particle character
      iMax = 1
      wMax = 0.0_dp
      do iT = 1, size(eigVec, dim=1)
        if (eigVec(iT, iState)**2 > wMax) then
          wMax = eigVec(iT, iState)**2
          iMax = iT
        end if
      end do
      write(stdOut, "(2X,I6,2X,F14.6,2X,F14.6,2X,I6,A,I6)") iState, eval(iState),&
          & eval(iState) * Hartree__eV, getIA(iMax, 1), " ->", getIA(iMax, 2)
      write(fdSF%unit, "(I6,2X,F16.8,2X,F16.8,2X,I6,A,I6,A,F6.3,A)") iState,&
          & eval(iState) * Hartree__eV, eval(iState), getIA(iMax, 1), " ->", getIA(iMax, 2),&
          & " (w=", wMax, ")"
    end do
    call closeFile(fdSF)

    if (fdTagged%isConnected()) then
      call taggedWriter%write(fdTagged%unit, tagLabels%excEgy, eval(1:nState))
    end if

  end subroutine writeSFResults

end module dftbp_timedep_linrespsf
