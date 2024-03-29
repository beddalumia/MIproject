program ed_hm_square
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                                       :: iloop,Nb,Lk,Nx,Nso,ik,iorb
  logical                                       :: converged
  real(8)                                       :: wband,wmixing
  real(8),dimension(5)                          :: ts,Dband
  real(8),dimension(:),allocatable              :: dens
  !Bath:
  real(8),allocatable                           :: Bath(:),Bath_(:)
  !The local hybridization function:
  complex(8),allocatable                        :: Hloc(:,:,:,:)
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Greal
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Weiss,Weiss_
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gkmats
  complex(8),allocatable,dimension(:)           :: Gtest
  !MIT markers
  real(8),allocatable                           :: entropy(:),luttinger(:)
  !
  character(len=16)                             :: finput,foutput
  complex(8),allocatable                        :: Hk(:,:,:)
  real(8),allocatable                           :: Wt(:)
  !
  integer                                       :: comm,rank
  logical                                       :: master
  logical                                       :: mixG0,symOrbs

  call init_MPI(comm,.true.)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  call parse_cmd_variable(finput,"FINPUT",default='inputHM.conf')
  call parse_input_variable(wmixing,"wmixing",finput,default=0.5d0,comment="Mixing parameter")
  call parse_input_variable(ts,"TS",finput,default=[0.25d0,0d0,0d0,0d0,0d0],comment="Hopping parameter (4ts=D)")
  call parse_input_variable(Dband,"Dband",finput,default=[0d0,0d0,0d0,0d0,0d0],comment="Crystal field splittig (bands shift)")
  call parse_input_variable(Nx,"Nx",finput,default=100,comment="Number of kx point for 2d BZ integration")
  call parse_input_variable(mixG0,"mixG0",finput,default=.false.,comment="T mixes the Weiss field, F mixes the bath (default behavior)")
  call parse_input_variable(symOrbs,"symOrbs",finput,default=.false.,comment="T imposes same bath for all orbitals, reading from the first one")
  !
  call ed_read_input(trim(finput),comm)
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  if(Nspin/=1.OR.Norb>5)stop "Wrong setup from input file: Nspin/=1 OR Norb>5"
  Nso=Nspin*Norb

  !Allocate Fields:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats),Weiss_(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats),Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  allocate(dens(Norb))
  allocate(Gtest(Lmats))

  !Build Hk
  call TB_set_bk(bkx=[pi2,0d0],bky=[0d0,pi2])
  Lk = Nx*Nx
  allocate(Hk(Nso,Nso,Lk),Wt(Lk))
  call TB_build_model(Hk(:,:,:),hk_model,Nso,[Nx,Nx])
  Wt = 1d0/Lk
  Hloc   = zero
  Hloc(1,1,:,:) = sum(Hk,dim=3)/Lk
  where(abs(dreal(Hloc))<1.d-6) Hloc=0d0
  
  if(master)call TB_write_hk(Hk(:,:,:),"Hk2d.dat",Nlat=1,&
                             Nspin=1,&
                             Norb=Norb,&
                             Nkvec=[Nx,Nx])

  !Setup solver
  Nb=ed_get_bath_dimension()
  allocate(bath(Nb))
  allocate(Bath_(Nb))
  call ed_init_solver(comm,bath)



  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(comm,bath,Hloc) 
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_realaxis(Sreal)
     
     !Compute the local entanglement entropy
     if(master)then
        allocate(entropy(Norb))
        call ed_local_ee(entropy)
        call print_local_ee(entropy)
        deallocate(entropy)
     endif

     !Compute the local gfs on the imaginary axis:
     call dmft_gloc_matsubara(Hk,Gmats,Smats)
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)

     !Get the Weiss field/Delta function to be fitted
     call dmft_self_consistency(Gmats,Smats,Weiss,Hloc,SCtype=cg_scheme)
     call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=1)

     if(mixG0)then
        if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_
        Weiss_=Weiss
     endif
     
     !Perform the SELF-CONSISTENCY by fitting the new bath
     if(symOrbs)then
        call ed_chi2_fitgf(comm,Weiss,bath,ispin=1,iorb=1)
        call ed_orb_equality_bath(bath,save=.true.)
     else
        call ed_chi2_fitgf(comm,Weiss,bath,ispin=1)
     endif

     !MIXING:
     if(.not.mixG0)then
        if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
        Bath_=Bath
     endif

     !Check convergence (if required change chemical potential)     
     Gtest=zero
     do iorb=1,Norb
        Gtest=Gtest+Weiss(1,1,iorb,iorb,:)/Norb
     enddo
     converged = check_convergence(Gtest,dmft_error,nsuccess,nloop,reset=.false.)
     call ed_get_dens(dens)
     if(nread/=0d0)call ed_search_variable(xmu,sum(dens),converged)

     call end_loop
  enddo

  !Get kinetic energy:
  call dmft_kinetic_energy(Hk,Smats)

  !Compute the local gfs on the real axis:
  call dmft_gloc_realaxis(Hk,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=1)

  !Compute the Luttinger invariants:
  if(master)then
     allocate(luttinger(Norb))
     call luttinger_integral(luttinger,Greal,Sreal)
     call print_luttinger(luttinger)
     deallocate(luttinger)
  endif

  call finalize_MPI()

contains

  !-------------------------------------------------------------------------------------------
  !PURPOSE:  Hk model for the 2d square lattice
  !-------------------------------------------------------------------------------------------
  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:) :: kpoint
    integer              :: N,ih
    real(8)              :: kx,ky
    complex(8)           :: hk(N,N)
    kx=kpoint(1)
    ky=kpoint(2)
    Hk = zero
    do ih=1,N
       Hk(ih,ih) = -one*2d0*ts(ih)*(cos(kx)+cos(ky)) + Dband(ih)
    enddo
  end function hk_model

  !+---------------------------------------------------------------------------+
  !PURPOSE : Compute the Luttinger integral IL = 1/π * |Im(∫dw(Gloc*dSloc/dw)|
  !          > Cfr. J. Phys. Cond. Mat. 28, 02560 and PRB B 102, 081110(R)
  ! NB) At ph-symmetry IL should be zero, but only if the T->0 limit is taken 
  !     before the mu->0 limit, and here we apparently invert the order of the 
  !     limits because of the numerical discretization on the imaginary axis. 
  !     > Look at [53] (report of a private communication) in the Phys. Rev. B
  !+---------------------------------------------------------------------------+
  subroutine luttinger_integral(IL,Gloc,Sloc)
    real(8),allocatable,dimension(:),intent(out)            :: IL
    complex(8),allocatable,dimension(:,:,:,:,:),intent(in)  :: Gloc,Sloc
    real(8),dimension(Lreal)                                :: dSreal,dSimag,integrand
    integer                                                 :: iorb,ispin
    !
    if(.not.allocated(IL)) allocate(IL(Norb))
    !
    do iorb=1,Norb
        do ispin=1,Nspin !Nspin is constrained to be equal to 1 above.
            dSreal(:) = deriv(dreal(Sloc(ispin,ispin,iorb,iorb,:)),1d0) !No need to include 1/dw if
            dSimag(:) = deriv(dimag(Sloc(ispin,ispin,iorb,iorb,:)),1d0) !we are integrating in dw...
            integrand = dimag(Gloc(ispin,ispin,iorb,iorb,:)*cmplx(dSreal,dSimag))
            IL(iorb) = sum(integrand)       !Naive Riemann-sum (no need of trapezoidal-sums with typical Lreal)
            IL(iorb) = 1/pi * abs(IL(iorb)) !The theoretical sign is determined by sign(mu)...
        enddo
    enddo
    !
  end subroutine luttinger_integral

  !+---------------------------------------------------------------------------+
  !PURPOSE : print to file (and stdout) the Luttinger invariants IL(Norb)
  !+---------------------------------------------------------------------------+
  subroutine print_luttinger(IL)
    real(8),allocatable,dimension(:),intent(in)  :: IL
    integer                                      :: iorb
    integer                                      :: unit
    !
    if(ed_verbose>0)then
       write(LOGfile,*) " "
       write(LOGfile,*) "Luttinger invariants:"
    endif
    !
    do iorb=1,Norb
       if(ed_verbose>0) write(LOGfile,*) iorb, IL(iorb)
       unit = free_unit()
       foutput = "luttinger_l"//str(iorb)//".dat"
       open(unit,file=foutput,action="write",position="rewind",status='unknown')
       write(unit,*) IL(iorb)
       close(unit)
    enddo
    !
    if(ed_verbose>0)then
       write(LOGfile,*) "          iorb        IL"
       write(LOGfile,*) " "
    endif
    !
  end subroutine print_luttinger

  !+---------------------------------------------------------------------------+
  !PURPOSE : build the local ee according to Eq.4 in Mod.Phys.Lett.B.2013.27:05
  !+---------------------------------------------------------------------------+
  subroutine ed_local_ee(EE)
    real(8),allocatable,dimension(:),intent(inout) :: EE
    real(8),allocatable,dimension(:)               :: pp
    real(8),allocatable,dimension(:)               :: dens_up,dens_dw,mag
    real(8),allocatable,dimension(:)               :: dens,docc
    integer                                        :: i
    !
    if(Norb>1)then
       write(LOGfile,*) "WARNING: for Norb>1 ed_local_ee traces down to single orbitals."
    endif
       !
    allocate(dens(Norb),dens_up(Norb),dens_dw(Norb),docc(Norb),mag(Norb))
    allocate(pp(4**Norb))
    !
    call ed_get_mag(mag)
    call ed_get_dens(dens)
    call ed_get_docc(docc)
    dens_up = 0.5d0*(dens + mag)
    dens_dw = 0.5d0*(dens - mag)
    !
    do i=1,Norb
       pp(1) = abs(1-dens_up(i)-dens_dw(i)+docc(i))
       pp(2) = abs(dens_up(i)-docc(i))
       pp(3) = abs(dens_dw(i)-docc(i))
       pp(4) = abs(docc(i))
       EE(i) = -sum(pp*log(pp)/log(2d0))
    enddo
    !
  end subroutine ed_local_ee

  !+---------------------------------------------------------------------------+
  !PURPOSE : print to file (and stdout) the entanglement entropies EE(Norb)
  !+---------------------------------------------------------------------------+
  subroutine print_local_ee(EE)
    real(8),allocatable,dimension(:),intent(in)  :: EE
    integer                                      :: iorb
    integer                                      :: unit
    !
    if(ed_verbose>0)then
       write(LOGfile,*) " "
       write(LOGfile,*) "Entanglement entropies:"
    endif
    !
    do iorb=1,Norb
       if(ed_verbose>0) write(LOGfile,*) iorb, EE(iorb)
       unit = free_unit()
       foutput = "eentropy_l"//str(iorb)//".dat"
       open(unit,file=foutput,action="write",position="rewind",status='unknown')
       write(unit,*) EE(iorb)
       close(unit)
    enddo
    !
    if(ed_verbose>0)then
       write(LOGfile,*) "          iorb        EE"
       write(LOGfile,*) " "
    endif
    !
  end subroutine print_local_ee

end program ed_hm_square

