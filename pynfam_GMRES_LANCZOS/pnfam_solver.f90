!------------------------------------------------------------------------------
!> This module contains the heart of the pnfam solver. It calls the setup routines
!> runs the iterative solver, and writes the output to file.
!>
!> @authors M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!> @authors D.M. Roh, LBNL Berkeley, 2026-
!               Implemented GMRES and Lanczos
!------------------------------------------------------------------------------

   !---------------------------------------------------------------------------
   ! The FAM equations
   ! (Ep+En-omega)X + dH^{20} = -F^{20}
   ! (Ep+En+omega)Y + dH^{02} = -F^{02}
   ! are equivalent to the linear response equation
   ! (|A  B | - omega|I  0|)|X| = -|F^{20}|
   ! (|B* A*|        |0 -I|)|Y|    |F^{02}|
   ! The strength function disctribution can equivalently be found by
   ! an eigendecomposition
   ! |A  B | |X  Y*| = |I  0| |X Y*| |omega  0|
   ! |B* A*| |Y  X*|   |0 -I| |Y X*| |0 -omega|
   !---------------------------------------------------------------------------

module pnfam_solver
   use pnfam_setup
   use type_extfield
   use type_bigblockmatrix

   implicit none
   integer, parameter, private :: dp = kind(1d0)

   logical, private :: use_diagonal_blocks, bminus

   ! DMR
   ! Optional mapping to undo the time-reversal permutation between stored X (m12)
   ! and stored Y (m21). When enabled, matvec treats the incoming/outgoing m21%elem
   ! vector as being in the same ordering as m12%elem, and converts internally.
   integer, allocatable, private :: perm_xy(:), invperm_yx(:)
   integer, allocatable, private :: sgn_xy(:), sgn_yx(:)
   ! DMR

   !---------------------------------------------------------------------------
   ! 'Green_sign' and 'Greens_power' determine the sign and the power of the
   ! Greens function, respectively.
   ! Greens function is defined as
   !               Greens_sign * (Ep+En-omega)^{Greens_power}
   !
   ! IFAM solver updates X and Y by
   ! X = -(Ep+En-omega)^{-1}(dH^{20}+F^{20})
   ! Y = -(Ep+En+omega)^{-1}(dH^{02}+F^{02})
   ! So, need Greens_sign = -1 and Greens_power = -1
   !
   ! GMRES solver is assigned a preconditioner 
   ! diag((Ep+En-omega)^{-1},(Ep+En+omega)^{-1})
   ! This results in a matvec defined as
   ! X <- X + (Ep+En-omega)^{-1}dH^{20}
   ! Y <- Y + (Ep+En+omega)^{-1}dH^{02}
   ! So, need Greens_sign = 1 and Greens_power = -1
   !
   ! LANCZOS solver performs matvecs with the matrix of the linear response eq,
   ! which requires access to multiplications with (EP+En-omega) and (Ep+En+omega).
   ! So, need Greens_sign = 1 and Greens_power = 1 
   !---------------------------------------------------------------------------
#ifdef IFAM
   integer, parameter :: Greens_sign = -1
#else
   integer, parameter :: Greens_sign = 1
#endif

#if defined(LANCZOS)
   integer, parameter :: Greens_power = 1
#else
   integer, parameter :: Greens_power = -1
#endif

   ! External field
   type(bigblockmatrix), private :: Fsp, Fqp
   type(bigblockmatrix), dimension(:), allocatable, private :: Gsp, Gqp
   ! HFB solution
   type(bigblockmatrix), private :: Wn, Wp
   ! FAM matrices
   type(bigblockmatrix), private :: dHqp_re, dHqp_im, dRqp_re, dRqp_im
   type(bigblockmatrix), private :: dHsp_re, dHsp_im, dRsp_re, dRsp_im
   type(bigblockmatrix), private :: Greens_re, Greens_im, T

   ! DMR
   ! For GMRES
   type(bigblockmatrix), private :: temp_re, temp_im
   type(bigblockmatrix), private :: Init_re, Init_im
   type(bigblockmatrix), private :: Out_re, Out_im

contains
   !---------------------------------------------------------------------------
   ! Subroutine for solving the FAM equations for one value of energy
   !---------------------------------------------------------------------------
   subroutine pnfam_solve
      use pnfam_broyden
      use pnfam_txtoutput
      use pnfam_storage, only : fam_io
      implicit none
      integer :: iexit=1, ierr=0

      ! Initialize FAM inputs (as blockmatrix types, where applicable)
      ! - Logging procedures
      ! - HFB solution (E,U,V)
      ! - External field (F in SP basis)
      ! - Routine for calculating dH
      !---------------------------------------------------------------
      call setup_pnfam
      call openlog(output_txtfile, print_stdout)
      call print_main_header
      call setup_extfield
      !call print_main_header !Testing purposes
      ! Prepare ifam
      !---------------------------------------------------------------
      ! 1. Allocate Block Matices:
      !    - F_QP, G_QP, dH_QP(H20,H02,H11,H11t), dH_SP(h,d),
      !      dR_QP(X,Y,P,Q), dR_SP(rho,kap), Greens(X,Y,P,Q), T(1ff,ff)
      ! 2. Transform F to QP basis
      ! 3. Define block structures (mostly based on F)
      ! 4. Populate statmatrix T and Greens function
      !---------------------------------------------------------------
      call init_pnfam_solver

      ! Initialize fam amplitudes
      !---------------------------------------------------------------
      call fam_io(1,iexit)
      if (iexit > 0) then
         si = 1
         call set_val_bbm(dRqp_re,0d0)
         call set_val_bbm(dRqp_im,0d0)
         call set_broinout_from_dRqp(dRqp_re, dRqp_im,'i')
      else
         call set_dRqp_from_broin(dRqp_re, dRqp_im)
      end if

      ! Run the iterative solver
      !---------------------------------------------------------------
      call log_header
      
      ! Compile directives for calling ifam or gmres or lanczos
      write(st,'(1x,72("-"))'); call writelog(st)
#ifdef IFAM
      write(st,'(30x,a)') 'CALLING IFAM'; call writelog(st)
      call ifam
#elif defined(GMRES)
      write(st,'(30x,a)') 'CALLING GMRES'; call writelog(st)      
      call gmres
#elif defined(LANCZOS)
      write(st,'(30x,a)') 'CALLING LANCZOS'; call writelog(st)
      call lanczos
#endif
      write(st,'(1x,72("-"))'); call writelog(st)

      ! Write the result to footer and store result
      !---------------------------------------------------------------
      call log_result
      call fam_io(2,iexit)
      call closelog

      !call transition_density

   end subroutine pnfam_solve

   !---------------------------------------------------------------------------
   ! Subroutine for solving the FAM equations for one value of energy. Before
   ! calling this the first time, the module data structures must be initialized
   ! with init_pnfam_solver.
   !---------------------------------------------------------------------------
   subroutine ifam
      use pnfam_constants, only : iu, pi, get_timer
      use pnfam_broyden
      implicit none

      ! Finite temperature
      real(dp)    :: ft_real_factor
      complex(dp) :: ft_cmplx_factor

      integer :: iter, ixterm
      real(dp) :: re_s, im_s
      real(dp) :: time1, time2, timei, timei0
      time1=0; time2=0; timei=0; timei0=0

      ! Iteration header
      write(st,'(1x,72("-"))'); call writelog(st)
      write(st,'(1x,5a)') '    i','               si',&
        '               Re(S)','               Im(S)','      Time'; call writelog(st)
      write(st,'(1x,72("-"))'); call writelog(st)

      ! Main iteration loop
      do iter=0,max_iter

         ! Record the iteration (or show where we are starting from)
         write(st,'(1x, i4,a1, 2x, f15.10, 2x, f18.10, 2x, f18.10, 2x, f8.3)') &
            iter, qrpa_bbroyden, si, re_str(1), im_str(1), timei; call writelog(st)
         if (allocated(g)) then
            do ixterm=1, size(g)
               write(st,'(1x,2x,a10,2x,"Re(X)=",f16.8,2x,"Im(X)=",f16.8)') &
                   trim(g(ixterm)%label), re_str(1+ixterm), im_str(1+ixterm); call writelog(st)
            end do
         end if

         ! Exit condition
         if (si < convergence_epsilon) then
            iter_conv = iter; exit
         end if
         if (iter==max_iter) then
            iter_conv = -1; exit
         end if

         ! Start the iteration
         time1 = get_timer()
         if (iter==0) timei0 = time1

         if (abs(quench_residual_int) < 1e-10) then
            call set_val_bbm(dHqp_im,0d0)
            call set_val_bbm(dHqp_re,0d0)
         else

            ! Compute the perturbed densities in SP space
            if (bminus) then
               call triprod_bbm('n',Wp,'n',dRqp_re,'t',Wn,dRsp_re)
               call triprod_bbm('n',Wp,'n',dRqp_im,'t',Wn,dRsp_im)
            else
               call triprod_bbm('n',Wn,'n',dRqp_re,'t',Wp,dRsp_re)
               call triprod_bbm('n',Wn,'n',dRqp_im,'t',Wp,dRsp_im)
            end if

            call calc_hamiltonian(&
                                  dRsp_re%m11, dRsp_im%m11, dRsp_re%m12, dRsp_im%m12,&
                                  dRsp_re%m22, dRsp_im%m22, dRsp_re%m21, dRsp_im%m21,&
                                  dHsp_re%m11, dHsp_im%m11, dHsp_re%m12, dHsp_im%m12,&
                                  dHsp_re%m22, dHsp_im%m22, dHsp_re%m21, dHsp_im%m21 &
                                  )

            ! Transform the Hamiltonian into quasiparticle basis
            if (bminus) then
               call triprod_bbm('t',Wp,'n',dHsp_re,'n',Wn,dHqp_re)
               call triprod_bbm('t',Wp,'n',dHsp_im,'n',Wn,dHqp_im)
            else
               call triprod_bbm('t',Wn,'n',dHsp_re,'n',Wp,dHqp_re)
               call triprod_bbm('t',Wn,'n',dHsp_im,'n',Wp,dHqp_im)
            end if

            ! Quench the residual interaction
            call scalar_mult_bbm(dHqp_re,quench_residual_int)
            call scalar_mult_bbm(dHqp_im,quench_residual_int)

         end if ! End computation of the residual interaction

         ! Add the external field to perturbed hamiltonian (assuming real Fqp)
         call add_bbm(Fqp, dHqp_re)

         ! Solve the new X/Y, P/Q
         call complex_emult_bbm(Greens_re, Greens_im, dHqp_re, dHqp_im,&
                                         dRqp_re, dRqp_im)
         if (use_diagonal_blocks) then
            call emult_bbm(T, dRqp_re, dRqp_re)
            call emult_bbm(T, dRqp_im, dRqp_im)
         end if

         ! Use Broyden mixing to stabilize convergence
         ! T > 0 / odd-A have twice the number of vectors to handle P and Q
         call qrpa_broyden(iter, dRqp_re, dRqp_im)

         ! Strength function
         call contract_bbm(Fqp,dRqp_re,re_s)
         call contract_bbm(Fqp,dRqp_im,im_s)
         re_str(1) = -re_s/pi
         im_str(1) = -im_s/pi

         ! Cross-terms
         if (allocated(g)) then
            do ixterm=1, size(g)
               call contract_bbm(Gqp(ixterm),dRqp_re,re_s)
               call contract_bbm(Gqp(ixterm),dRqp_im,im_s)
               re_str(1+ixterm) = -re_s/pi
               im_str(1+ixterm) = -im_s/pi
            end do
         end if

         ! Record the iteration time
         time2 = get_timer()
         timei = (time2-time1)

      end do

      write(st,'(1x,72("-"))'); call writelog(st)
      if (iter_conv >= 0 ) then
         write(st,'("  *  FAM iteration converged after   ",i4," steps.  si = ",ES12.5)') iter,si; call writelog(st)
      else
         write(st,'("  *  FAM iteration interrupted after ",i4," steps.  si = ",ES12.5)') iter,si; call writelog(st)
      end if
      write(st,'(1x,72("-"))'); call writelog(st)
      write(st,'(1x,a,es12.3,a)') 'Total CPU time = ',(time2-timei0)/60.0_dp,' minutes' ; call writelog(st)
      call writelog("")

   end subroutine ifam


   !---------------------------------------------------------------------------
   ! GMRES subroutine for solving the FAM equations for one value of energy. 
   ! Before calling this the first time, the module data structures must be 
   ! initialized with init_pnfam_solver.
   !
   ! pnFAM equations for the linear response can be desribed as a linear eq
   !                     HZ = b
   ! where Z contains X and Y, and b is -Fqp.
   ! This GMRES solves a preconditioned version, 
   !                    PHZ = Pb
   ! where the preconditioner P is inv(diag((Ep+En-omega);(Ep+En+omega))),
   ! with an option to restart, i.e., reset the Krylov subspace after few its.
   !---------------------------------------------------------------------------
   subroutine gmres
      use pnfam_constants, only : iu, pi, get_timer
      use pnfam_broyden
      implicit none

      complex(dp), parameter :: cone=(1.0d0, 0.d0), czero=(0.d0, 0.d0)  ! declaration of complex one and complex zero
      real(dp), external :: dznrm2                                      ! LAPACK function for euclidean norm

      complex(dp), allocatable :: basis(:, :), Hmat(:, :)               ! Krylov subspace & Hessenberg matrix

      complex(dp), allocatable :: sn(:), cn(:)                          ! arrays holding Givens rotations
      complex(dp), allocatable :: betavec(:), errorvec(:)               ! arrays holding beta & residual norm values
      complex(dp), allocatable :: mres_x(:), Z0(:)                      ! gmres solution vector & initial vector

      real(dp) :: betanrm0, betanrm, error                              ! initial residual norm & beta & residual norm
      
      integer :: restart, ncycle                                        ! number of restarts & cycle counter
      integer :: counter                                                ! counter for inner iteration (gmres)
      
      ! Finite temperature
      real(dp)    :: ft_real_factor
      complex(dp) :: ft_cmplx_factor

      integer :: iter, ixterm, nxy, i, j
      real(dp) :: re_s, im_s
      real(dp) :: time1, time2, timei, timei0
      time1=0; time2=0; timei=0; timei0=0      

      ! Iteration header
      write(st,'(1x,72("-"))'); call writelog(st)
      write(st,'(1x,3a)') '    i', '              Error', '               Time'; call writelog(st)
      write(st,'(1x,72("-"))'); call writelog(st)

      ! Allocate memory
      nxy = size(f%mat%elem)
      allocate(basis(2*nxy, max_iter+1))
      allocate(Hmat(max_iter+1, max_iter))
      allocate(mres_x(2*nxy))
      allocate(Z0(2*nxy))

      allocate(betavec(max_iter+1))
      allocate(cn(max_iter))
      allocate(sn(max_iter))
      allocate(errorvec(max_iter+1))

      ! Obtain preconditioned rhs
      ! Set Init_ as Fqp
      Init_re = Fqp
      Init_im = Fqp
      Out_re = Fqp
      Out_im = Fqp
      do i = 1, nxy
         Init_im%m12%elem(i) = 0.0
         Init_im%m21%elem(i) = 0.0
      end do
      ! Apply P * Init_, which is set to Out_
      call complex_emult_bbm(Greens_re, Greens_im, Init_re, Init_im, Out_re, Out_im)

      ! Set the number of restarts
      restart = 3
      do ncycle=1,restart

         write(st,'(" ****** At cycle number ", i4, " *****.")') ncycle; call writelog(st)

         ! Start the cycle
         time1 = get_timer()
         if (ncycle==1) timei0 = time1

         basis = czero
         Hmat = czero
         cn = czero
         sn = czero
         betavec = czero
         errorvec = 0.d0

         if (ncycle == 1) then
            ! For fist ever initialization, Z0 is zero and r0 is just b=-Fqp
            Z0 = czero
            do i = 1, nxy
               basis(i, 1) = dcmplx(-Out_re%m12%elem(i), -Out_im%m12%elem(i))
               basis(nxy+i, 1) = dcmplx(-Out_re%m21%elem(i), -Out_im%m21%elem(i))
            end do
            betanrm0 = dznrm2(2*nxy, basis(1, 1), 1) ! use to compute relative residual norm
         else
            ! For restart, need b-H*Z0
            ! Keep Z0 for final solution at the end of cycle
            Z0 = mres_x
            ! Do matvec H*Z0 (don't use the subroutine matvec)
            if (abs(quench_residual_int) < 1e-10) then
               call set_val_bbm(dHqp_im,0d0)
               call set_val_bbm(dHqp_re,0d0)
            else

               ! Compute the perturbed densities in SP space
               if (bminus) then
                  call triprod_bbm('n',Wp,'n',dRqp_re,'t',Wn,dRsp_re)
                  call triprod_bbm('n',Wp,'n',dRqp_im,'t',Wn,dRsp_im)
               else
                  call triprod_bbm('n',Wn,'n',dRqp_re,'t',Wp,dRsp_re)
                  call triprod_bbm('n',Wn,'n',dRqp_im,'t',Wp,dRsp_im)
               end if

               call calc_hamiltonian(& 
                                     dRsp_re%m11, dRsp_im%m11, dRsp_re%m12, dRsp_im%m12,&
                                     dRsp_re%m22, dRsp_im%m22, dRsp_re%m21, dRsp_im%m21,&
                                     dHsp_re%m11, dHsp_im%m11, dHsp_re%m12, dHsp_im%m12,&
                                     dHsp_re%m22, dHsp_im%m22, dHsp_re%m21, dHsp_im%m21 &
                                    )

               ! Transform the Hamiltonian into quasiparticle basis
               if (bminus) then
                  call triprod_bbm('t',Wp,'n',dHsp_re,'n',Wn,dHqp_re)
                  call triprod_bbm('t',Wp,'n',dHsp_im,'n',Wn,dHqp_im)
               else
                  call triprod_bbm('t',Wn,'n',dHsp_re,'n',Wp,dHqp_re)
                  call triprod_bbm('t',Wn,'n',dHsp_im,'n',Wp,dHqp_im)
               end if

               ! Quench the residual interaction
               call scalar_mult_bbm(dHqp_re,quench_residual_int)
               call scalar_mult_bbm(dHqp_im,quench_residual_int)

            end if ! End computation of the residual interaction

            ! Apply preconditioner to dH : P*dH
            temp_re = dRqp_re ! arbitrary
            temp_im = dRqp_im ! arbitrary
            call complex_emult_bbm(Greens_re, Greens_im, dHqp_re, dHqp_im,&
                                            temp_re, temp_im)
            ! Combine Z0 + P*dH
            call add_bbm(temp_re, dRqp_re)
            call add_bbm(temp_im, dRqp_im)
            ! End of matvec

            ! Get b - H*Z0 = -Fqp - H*Z0 = -(Fqp + H*Z0)
            call add_bbm(Out_re, dRqp_re)
            call add_bbm(Out_im, dRqp_im)
            do i = 1, nxy
               basis(i, 1) = dcmplx(-dRqp_re%m12%elem(i), -dRqp_im%m12%elem(i))
               basis(nxy+i, 1) = dcmplx(-dRqp_re%m21%elem(i), -dRqp_im%m21%elem(i))
            end do
         end if

         ! Normalize the first basis vector
         betanrm = dznrm2(2*nxy, basis(1, 1), 1)
         basis(:, 1) = basis(:, 1) / betanrm
         ! Initialize residual
         betavec(1) = betanrm
         errorvec(1) = betanrm

         ! Record the iteration time
         time2 = get_timer()
         timei = (time2-time1)

         ! Record the initial result
         write(st,'(1x, i4, 2x, f20.10, 2x, f15.3)') 0, betanrm, timei; call writelog(st)
         ! End of initialization step

         ! Main gmres loop
         do iter=1,max_iter

            ! Start the iteration
            time1 = get_timer()

            ! Set Z (X and Y)
            do i=1,nxy
               dRqp_re%m12%elem(i) = dble(basis(i, iter))
               dRqp_re%m21%elem(i) = dble(basis(nxy+i, iter))
               dRqp_im%m12%elem(i) = dimag(basis(i, iter))
               dRqp_im%m21%elem(i) = dimag(basis(nxy+i, iter))
            end do 

            ! Do matvec H*Z

            if (abs(quench_residual_int) < 1e-10) then
               call set_val_bbm(dHqp_im,0d0)
               call set_val_bbm(dHqp_re,0d0)
            else

               ! Compute the perturbed densities in SP space
               if (bminus) then
                  call triprod_bbm('n',Wp,'n',dRqp_re,'t',Wn,dRsp_re)
                  call triprod_bbm('n',Wp,'n',dRqp_im,'t',Wn,dRsp_im)
               else
                  call triprod_bbm('n',Wn,'n',dRqp_re,'t',Wp,dRsp_re)
                  call triprod_bbm('n',Wn,'n',dRqp_im,'t',Wp,dRsp_im)
               end if

               call calc_hamiltonian(&
                                     dRsp_re%m11, dRsp_im%m11, dRsp_re%m12, dRsp_im%m12,&
                                     dRsp_re%m22, dRsp_im%m22, dRsp_re%m21, dRsp_im%m21,&
                                     dHsp_re%m11, dHsp_im%m11, dHsp_re%m12, dHsp_im%m12,&
                                     dHsp_re%m22, dHsp_im%m22, dHsp_re%m21, dHsp_im%m21 &
                                     )

               ! Transform the Hamiltonian into quasiparticle basis
               if (bminus) then
                  call triprod_bbm('t',Wp,'n',dHsp_re,'n',Wn,dHqp_re)
                  call triprod_bbm('t',Wp,'n',dHsp_im,'n',Wn,dHqp_im)
               else
                  call triprod_bbm('t',Wn,'n',dHsp_re,'n',Wp,dHqp_re)
                  call triprod_bbm('t',Wn,'n',dHsp_im,'n',Wp,dHqp_im)
               end if

               ! Quench the residual interaction
               call scalar_mult_bbm(dHqp_re,quench_residual_int)
               call scalar_mult_bbm(dHqp_im,quench_residual_int)

            end if ! End computation of the residual interaction
        
            ! Apply preconditioner to dH : P*dH 
            temp_re = dRqp_re ! arbitrary
            temp_im = dRqp_im ! arbitrary
            call complex_emult_bbm(Greens_re, Greens_im, dHqp_re, dHqp_im,&
                                            temp_re, temp_im)
            ! Combine Z + P*dH
            call add_bbm(temp_re, dRqp_re)
            call add_bbm(temp_im, dRqp_im)
            ! End of matvec

            ! Store into Krylov subspace
            do i = 1,nxy
               basis(i, iter+1) = dcmplx(dRqp_re%m12%elem(i), dRqp_im%m12%elem(i))
               basis(nxy+i, iter+1) = dcmplx(dRqp_re%m21%elem(i), dRqp_im%m21%elem(i))
            end do
            ! Orthogonalize w.r.t. previous Krylov basis vectors
            call zgemv('C', 2*nxy, iter, cone, basis, 2*nxy, basis(:, iter+1), 1, czero, Hmat(1:iter, iter), 1)
            call zgemv('N', 2*nxy, iter, -cone, basis, 2*nxy, Hmat(1:iter, iter), 1, cone, basis(:, iter+1), 1)

            ! Normalize current basis vector
            betanrm = dznrm2(2*nxy, basis(1, iter+1), 1)
            basis(:, iter+1) = basis(:, iter+1) / betanrm
            ! Set subdiag of Hessenberg
            Hmat(iter+1, iter) = betanrm

            ! Apply Givens rotations
            do i = 1,iter-1
               call zrot(1, Hmat(i, iter), 1, Hmat(i+1, iter), 1, cn(i), sn(i))   
            end do
            call zrotg(Hmat(iter, iter), Hmat(iter+1, iter), cn(iter), sn(iter))
            Hmat(iter+1, iter) = czero
            ! Residual vector
            call zrot(1, betavec(iter), 1, betavec(iter+1), 1, cn(iter), sn(iter))
            error = abs(betavec(iter+1)) / betanrm0
            errorvec(iter+1) = error

            ! Record the iteration time
            time2 = get_timer()
            timei = (time2-time1)

            ! Record the iteration (or show where we are starting from)
            write(st,'(1x, i4, 2x, f20.10, 2x, f15.3)') iter, error, timei; call writelog(st)

            counter = iter ! Need separate counter for correct number of iterations

            ! Exit condition
            if (error < convergence_epsilon) then
               iter_conv = iter; exit
            end if

         end do ! end inner iteration

         ! Solve triangular problem and get GMRES result mres_x
         call ztrsv('U', 'N', 'N', counter, Hmat, max_iter+1, betavec, 1)
         call zgemv('N', 2*nxy, counter, cone, basis, 2*nxy, betavec, 1, czero, mres_x, 1)
         ! Combine Z0 + mres_x
         mres_x = Z0 + mres_x

         do i=1,nxy
            dRqp_re%m12%elem(i) = dble(mres_x(i))
            dRqp_re%m21%elem(i) = dble(mres_x(nxy+i))
            dRqp_im%m12%elem(i) = imag(mres_x(i))
            dRqp_im%m21%elem(i) = imag(mres_x(nxy+i))
         end do

         ! Strength function
         call contract_bbm(Fqp,dRqp_re,re_s)
         call contract_bbm(Fqp,dRqp_im,im_s)
         re_str(1) = -re_s/pi
         im_str(1) = -im_s/pi

         ! Cross-terms
         if (allocated(g)) then
            do ixterm=1, size(g)
               call contract_bbm(Gqp(ixterm),dRqp_re,re_s)
               call contract_bbm(Gqp(ixterm),dRqp_im,im_s)
               re_str(1+ixterm) = -re_s/pi
               im_str(1+ixterm) = -im_s/pi
            end do
         end if

         time2 = get_timer()

         ! Exit restart if converged
         if (iter_conv > 0) then 
            exit
         end if
         if (ncycle==restart) then
            iter_conv = -1; exit
         end if

      end do !end cycle

      write(st,'(1x,72("-"))'); call writelog(st)
      if (iter_conv >= 0) then
         write(st,'("  *  GMRES iteration converged after   ", i4, " steps.")') (ncycle-1)*max_iter+iter; call writelog(st)
      else
         write(st,'("  *  GMRES iteration did not converge within step   ", i4, " steps.")') (ncycle-1)*max_iter+iter; call writelog(st)
      end if
      write(st,'(1x,72("-"))'); call writelog(st)
      write(st,'(1x,a,es12.3,a)') 'Total CPU time = ',(time2-timei0)/60.0_dp,' minutes' ; call writelog(st)
      call writelog("")

      ! Deallocate
      deallocate(basis)
      deallocate(Hmat)
      deallocate(mres_x)
      deallocate(Z0)

      deallocate(betavec)
      deallocate(cn)
      deallocate(sn)
      deallocate(errorvec)

   end subroutine gmres


   !---------------------------------------------------------------------------
   ! LANCZOS subroutine for computing the response function by solving the
   ! (generalized) eigenvalue problem.
   ! Before calling this the first time, the module data structures must be
   ! initialized with init_pnfam_solver.
   !
   !                    | A  B| = |X Y||D  ||X Y|^-1
   !                    |-B -A|   |Y X|| -D||Y X|
   !
   ! where X^TX-Y^TY=I and X^Y-Y^X=0, 
   ! and eigenvalues in D are in non-decreasing order.
   !
   ! Solves the reduced eigenvalue problem
   !                         MK = (X-Y)D^2(X-Y)^-1
   ! where M:=A+B and K:=A-B.
   !
   ! We have access to the matvec with the matrix [A B; B A].
   ! To obtain the matvecs with the matrix M and K, we set the vector as [Z;0]
   ! so that            |A B||Z| = |AZ|
   !                    |B A||0|   |BZ|
   ! then do AZ+BZ for MZ or AZ-BZ for KZ.
   ! After m iterations, returns Theta, containing the Ritz values,
   ! and RF0 and RF1, where
   ! |RF0| = 1/2 |(KZ*Theta^{-1/4})(F20+F02)+(Z*Theta^{1/4})(F20-F02)|
   ! |RF1|       |(KZ*Theta^{-1/4})(F20+F02)-(Z*Theta^{1/4})(F20-F02)| 
   ! where Z is the Ritz vectors.
   !---------------------------------------------------------------------------
   subroutine lanczos
      use pnfam_constants, only : iu, pi, get_timer
      use pnfam_broyden
      implicit none

      ! Parameters to keep from original code
      integer :: iter, ixterm, nxy, i, j
      real(dp) :: re_s, im_s
      real(dp) :: time1, time2, timei, timei0

      integer :: lu1, lu2, lu3, lu4                                     ! Write integer files for Theta_MK, Theta_KM, Tau2_MK, Tau2_KM

      complex(dp), parameter :: cone=(1.0d0, 0.d0), czero=(0.d0, 0.d0)  ! declaration of complex one and complex zero
      real(dp), external :: dznrm2                                      ! LAPACK function for euclidean norm

      ! Krylov subspace parameters
      complex(dp), allocatable :: basis_MK(:, :)                        ! Krylov subspaces
      complex(dp), allocatable :: Hmat_MK(:, :)                         ! Hessenberg matrices
      complex(dp), allocatable :: copyHmat_MK(:, :)                     ! Copy of Hessenberg matrices; for reorthogonalization and zgeev
      real(dp) :: betanrm                                               ! For norm of the basis vectors
      complex(dp), allocatable :: ZZ(:), AZ(:), BZ(:)                   ! For storing intermediate matvecs

      ! Eigen parameters
      complex(dp), allocatable :: Theta_MK(:)                           ! Ritz values
      integer, allocatable :: idx(:)                                    ! For sorting eigenvalues
      integer :: temp_idx, k                                            ! For sorting eigenvalues
      complex(dp) :: temp                                               ! For sorting eigenvalues
      complex(dp), allocatable :: Temp_reigvecs(:, :)                   ! For sorting eigenvalues
      complex(dp), allocatable :: ritzvecs_MK(:, :)                     ! Ritz vectors
      complex(dp), allocatable :: reigvecs_MK(:, :)                     ! Right eigenvectors
      integer :: info_MK                                                ! Info values (=0 indicates convergence)
      real(dp), allocatable :: DE(:), SDE(:)                            ! Diagonal & subdiagonal elements of Hmat
      real(dp), allocatable :: eigvecs(:, :)                            ! Eigenvectors for dstev
      real(dp), allocatable :: workvec(:)                               ! Work vector for dstev

      ! Projected transition amplitudes
      complex(dp), allocatable :: F20(:), F02(:)                        ! Weak external fields
      complex(dp), allocatable :: VTF20(:), UTF02(:)                    ! Holds intermediate values
      complex(dp), allocatable :: Z(:, :)                               ! Holds Ritz vectors
      complex(dp), allocatable :: KZ(:, :)                              ! Holds K * Ritz vectors
      complex(dp), allocatable :: RF0(:), RF1(:)                        ! Projected transition amplitudes

      nxy = size(f%mat%elem)                                            ! Size of A and B

      ! Iteration header
      write(st,'(1x,72("-"))'); call writelog(st)
      write(st,'(1x,3a)') '    i', '              Error', '               Time'; call writelog(st)
      write(st,'(1x,72("-"))'); call writelog(st)

      ! Allocate memory
      allocate(basis_MK(nxy, max_iter+1))
      allocate(Hmat_MK(max_iter+1, max_iter))
      allocate(copyHmat_MK(max_iter+1, max_iter))
      allocate(ZZ(nxy))
      allocate(AZ(nxy))
      allocate(BZ(nxy))
      allocate(reigvecs_MK(max_iter, max_iter))
      allocate(Theta_MK(max_iter))
      allocate(idx(max_iter))
      allocate(Temp_reigvecs(max_iter, max_iter))
      allocate(ritzvecs_MK(nxy, max_iter))
      allocate(F20(nxy))
      allocate(F02(nxy))
      allocate(VTF20(max_iter))
      allocate(UTF02(max_iter))
      allocate(DE(max_iter))
      allocate(SDE(max_iter-1))
      allocate(eigvecs(max_iter, max_iter))
      allocate(workvec(2*max_iter-2))
      allocate(Z(nxy, max_iter))
      allocate(KZ(nxy, max_iter))
      allocate(RF0(max_iter))
      allocate(RF1(max_iter))

      time1=0; time2=0; timei=0; timei0=0
      Theta_MK = czero
      basis_MK = czero
      Hmat_MK = czero
      copyHmat_MK = czero
      reigvecs_MK = czero
      Temp_reigvecs = czero
      ritzvecs_MK = czero

      ! Start Lanczos
      time1 = get_timer()
      timei0 = time1

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!! Initialize !!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      do i = 1, nxy
         F20(i) = Fqp%m12%elem(i)
         if (allocated(perm_xy)) then
            F02(i) = dble(sgn_xy(i))*Fqp%m21%elem(perm_xy(i))
         else
            F02(i) = Fqp%m21%elem(i)
         end if
         basis_MK(i, 1) = F20(i) + F02(i)
      end do
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!! End of Initialization !!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!! Normalize the first basis vector !!!!
      !    For MK problem, normalize w.r.t. K  !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 
      ! a) MK problem
      !
      ! Set [Z;0]
      do i=1,nxy
         dRqp_re%m12%elem(i) = dble(basis_MK(i, 1))
         dRqp_re%m21%elem(i) = 0.d0
         dRqp_im%m12%elem(i) = dimag(basis_MK(i, 1))
         dRqp_im%m21%elem(i) = 0.d0
      end do
      ! Do matvec and get [AZ;BZ]
      call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
      do i = 1, nxy
         AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
         BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
      end do
      ! KZ = AZ - BZ
      ZZ = AZ - BZ
      betanrm = sqrt(real(dot_product(basis_MK(:, 1), ZZ)))
      basis_MK(:, 1) = basis_MK(:, 1) / betanrm
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End of first basis vector normalization!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Record the iteration time
      time2 = get_timer()
      timei = (time2-time1)

      ! Record the initial result
      write(st,'(1x, i4, 2x, f15.3)') iter, timei; call writelog(st)
      ! End of initialization step

      ! Main lanczos loop
      do iter=1,max_iter

         ! Start the iteration
         time1 = get_timer()

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!! Generate Krylov subspace for MK !!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! a) Matvec with MK
         !
         ! K * basis_MK(:, iter) == ZZ / betanrm
         do i=1,nxy
            dRqp_re%m12%elem(i) = dble(ZZ(i)/betanrm)
            dRqp_re%m21%elem(i) = 0.d0
            dRqp_im%m12%elem(i) = dimag(ZZ(i)/betanrm)
            dRqp_im%m21%elem(i) = 0.d0
         end do
         ! Do matvec and get [AKZ;BKZ]
         call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
         do i = 1, nxy
            AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
            BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
         end do
         ! MKZ = AKZ+BKZ
         ZZ = AZ + BZ
         ! Store into basis
         do i = 1,nxy
            basis_MK(i, iter+1) = ZZ(i)
         end do
         !
         ! b) K-inner product of MKZ
         !
         ! Set [MKZ;0]
         do i=1,nxy
            dRqp_re%m12%elem(i) = dble(ZZ(i))
            dRqp_re%m21%elem(i) = 0.d0
            dRqp_im%m12%elem(i) = dimag(ZZ(i))
            dRqp_im%m21%elem(i) = 0.d0
         end do
         ! Do matvec and get [AMKZ;BMKZ]
         call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
         do i = 1, nxy
            AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
            BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
         end do
         ! KMKZ = AMKZ-BMKZ
         ZZ = AZ - BZ
         !
         ! c) Orthogonalize w.r.t. previous Krylov basis vectors
         !
         call zgemv('C', nxy, iter, cone, basis_MK, nxy, ZZ, 1, czero, Hmat_MK(1:iter, iter), 1)
         call zgemv('N', nxy, iter, -cone, basis_MK, nxy, Hmat_MK(1:iter, iter), 1, cone, basis_MK(:, iter+1), 1)
         !
         ! d) Reorthogonalize
         !
         !    d)a) K-inner product of f:=basis_MK(:, iter+1)
         !
         ! Set [f;0]
         do i=1,nxy
            dRqp_re%m12%elem(i) = dble(basis_MK(i, iter+1))
            dRqp_re%m21%elem(i) = 0.d0
            dRqp_im%m12%elem(i) = dimag(basis_MK(i, iter+1))
            dRqp_im%m21%elem(i) = 0.d0
         end do
         ! Do matvec and get [Af;Bf]
         call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
         do i = 1, nxy
            AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
            BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
         end do
         ! Kf = Af-Bf
         ZZ = AZ - BZ
         !
         !    d)b) Orthogonalize w.r.t. previous Krylov basis vectors
         !
         call zgemv('C', nxy, iter, cone, basis_MK, nxy, ZZ, 1, czero, copyHmat_MK(1:iter, iter), 1)
         call zgemv('N', nxy, iter, -cone, basis_MK, nxy, copyHmat_MK(1:iter, iter), 1, cone, basis_MK(:, iter+1), 1)
         Hmat_MK(1:iter, iter) = Hmat_MK(1:iter, iter) + copyHmat_MK(1:iter, iter)
         ! Ensure Hmat is tridiagonal
         if (iter > 2) then
            do i = 1, iter-2
               Hmat_MK(i, iter) = czero
            enddo
         end if
         !
         ! e) Normalize current basis vector w.r.t. K
         !
         ! Set [Z;0]
         do i=1,nxy
            dRqp_re%m12%elem(i) = dble(basis_MK(i, iter+1))
            dRqp_re%m21%elem(i) = 0.d0
            dRqp_im%m12%elem(i) = dimag(basis_MK(i, iter+1))
            dRqp_im%m21%elem(i) = 0.d0
         end do
         ! Do matvec and get [AZ;BZ]
         call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
         do i = 1, nxy
            AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
            BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
         end do
         ! KZ = AZ-BZ
         ZZ = AZ - BZ
         betanrm = sqrt(real(dot_product(basis_MK(:, iter+1), ZZ)))
         basis_MK(:, iter+1) = basis_MK(:, iter+1) / betanrm
         Hmat_MK(iter+1, iter) = betanrm
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! End Krylov subspace generation for MK !!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         ! Store tridiagonal elements for dstev
         DE(iter) = Hmat_MK(iter, iter)
         if (iter < max_iter) then
            SDE(iter) = betanrm
         endif

         time2 = get_timer()
         timei = (time2-time1)
         write(st,'(1x, i4, 2x, f15.3)') iter, timei; call writelog(st)

      end do ! end inner iteration

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!! Solve eigenvalue problems !!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! MK eigenvalue problem
      !
      call dstev('V', max_iter, DE, SDE, eigvecs, max_iter, workvec, info_MK) 
      do j = 1, max_iter
         Theta_MK(j) = DE(j)
         reigvecs_MK(:, j) = eigvecs(:, j)
      enddo
      ! Initialize indices
      do i = 1, max_iter
         idx(i) = i
      end do
      ! Simple bubble sort
      do i = 1, max_iter-1
         k = i
         do j = i+1, max_iter
            if (real(Theta_MK(j)) < real(Theta_MK(k))) k = j
         end do
         if (k /= i) then
            ! swap values
            temp = Theta_MK(i); Theta_MK(i) = Theta_MK(k); Theta_MK(k) = temp
            ! swap indices
            temp_idx = idx(i); idx(i) = idx(k); idx(k) = temp_idx
         end if
      end do
      Temp_reigvecs = reigvecs_MK
      do i = 1, max_iter
         reigvecs_MK(:, i) = Temp_reigvecs(:, idx(i))
      end do
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!! End Solve eigenvalue problems !!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Compute Ritz vectors
      call zgemm('N', 'N', nxy, max_iter, max_iter, cone, basis_MK, nxy, reigvecs_MK, max_iter, czero, ritzvecs_MK, nxy)
      ! Want KZ\Omega^{-1/2} and Z\Omega^{1/2}
      ! where Z is ritzvecs_MK
      do k = 1, max_iter
         ! Set [Z;0]
         do i=1,nxy
            dRqp_re%m12%elem(i) = dble(ritzvecs_MK(i, k))
            dRqp_re%m21%elem(i) = 0.d0
            dRqp_im%m12%elem(i) = dimag(ritzvecs_MK(i, k))
            dRqp_im%m21%elem(i) = 0.d0
         end do
         ! Do matvec and get [AZ;BZ]
         call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
         do i = 1, nxy
            AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
            BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
         end do
         ! KZ = AZ - BZ
         ZZ = AZ - BZ
         KZ(:, k) = ZZ * (real(Theta_MK(k)) ** (-0.25))
         Z(:, k) = ritzvecs_MK(:, k) * (real(Theta_MK(k)) ** (0.25))
      end do
      ! VTF20 = KZ\Omega^{-1/2}(F20+F02)
      ! UTF02 = Z\Omega^{1/2}(F20-F02)
      call zgemv('C', nxy, max_iter, cone, KZ, nxy, F20+F02, 1, czero, VTF20, 1)
      call zgemv('C', nxy, max_iter, cone, Z, nxy, F20-F02, 1, czero, UTF02, 1)
      RF0 = (VTF20 + UTF02) / 2;
      RF1 = (VTF20 - UTF02) / 2;

      ! Keeping below because.. original and for writing into file
      ! Strength function
      call contract_bbm(Fqp,dRqp_re,re_s)
      call contract_bbm(Fqp,dRqp_im,im_s)
      re_str(1) = -re_s/pi
      im_str(1) = -im_s/pi

      ! Cross-terms
      if (allocated(g)) then
         do ixterm=1, size(g)
            call contract_bbm(Gqp(ixterm),dRqp_re,re_s)
            call contract_bbm(Gqp(ixterm),dRqp_im,im_s)
            re_str(1+ixterm) = -re_s/pi
            im_str(1+ixterm) = -im_s/pi
         end do
      end if

      time2 = get_timer()

      lu1 = 99
      lu2 = 100
      lu3 = 101
      open(unit=lu1, file='Theta.txt',   status='unknown', action='write')
      open(unit=lu2, file='RF0.txt',   status='unknown', action='write')
      open(unit=lu3, file='RF1.txt',   status='unknown', action='write')
      do i = 1, max_iter
         write(lu1,'(100(1X,ES20.12))') real(Theta_MK(i))
         write(lu2,'(100(1X,ES20.12))') real(RF0(i))
         write(lu3,'(100(1X,ES20.12))') real(RF1(i))
      end do
      close(lu1)
      close(lu2)
      close(lu3)

      write(st,'(1x,72("-"))'); call writelog(st)
      write(st,'("  *  LANCZOS ran for ", i0, " iterations.")') iter-1; call writelog(st)
      write(st,'(1x,72("-"))'); call writelog(st)
      write(st,'(1x,a,es12.3,a)') 'Total CPU time = ',(time2-timei0)/60.0_dp,' minutes' ; call writelog(st)
      write(st,'("  *  The information below is not relevant to the LANCZOS procedure; use Theta, RF0, and RF1 values instead to compute the strength function distribution.")'); call writelog(st) 
      call writelog("")

      ! Deallocate memory
      deallocate(basis_MK)
      deallocate(Hmat_MK)
      deallocate(copyHmat_MK)
      deallocate(ZZ)
      deallocate(AZ)
      deallocate(BZ)
      deallocate(reigvecs_MK)
      deallocate(Theta_MK)
      deallocate(idx)
      deallocate(Temp_reigvecs)
      deallocate(ritzvecs_MK)
      deallocate(F20)
      deallocate(F02)
      deallocate(VTF20)
      deallocate(UTF02)
      deallocate(DE)
      deallocate(SDE)
      deallocate(eigvecs)
      deallocate(workvec)
      deallocate(Z)
      deallocate(KZ)
      deallocate(RF0)
      deallocate(RF1)

   end subroutine lanczos


   !---------------------------------------------------------------------------
   ! Allocates and initializes the data structures used by pnfam_solve.
   !---------------------------------------------------------------------------
   subroutine init_pnfam_solver
      use pnfam_constants, only : IT_NEUTRON, IT_PROTON
      use pnfam_broyden, only : qrpa_broin, qrpa_broout, nbroy, broyden_history_size

      implicit none

      integer :: ixterm, i
      logical :: a11,a12,a21,a22

      !---------------------------------------------------------------------------
      ! General initializations
      !---------------------------------------------------------------------------
      ! The complex strength function (array [strength, xterms])
      if (allocated(re_str)) deallocate(re_str, im_str)
      allocate(re_str(1:nxterms+1), im_str(1:nxterms+1))
      re_str(:) = 0; im_str(:) = 0

      ! While the various FAM matrices may have different non-trivial block structures
      ! they still share the same number of elements in the non-trivial blocks
      nxy = size(f%mat%elem)
      bminus = f%beta_minus

      if (ft_active .or. hfb_blo_active) then
         use_diagonal_blocks = .true.
         nbroy = nxy*8
      else
         use_diagonal_blocks = .false.
         nbroy = nxy*4
      end if

      ! Broyden arrays
      if (allocated(qrpa_broin)) deallocate(qrpa_broin,qrpa_broout)
      allocate(qrpa_broin(nbroy),qrpa_broout(nbroy))
      qrpa_broin = 0; qrpa_broout=0
      si = 1
      ! Don't do broyden if no residual interaction, should just need 1 iteration
      if (abs(quench_residual_int) < 1e-10) then
         broyden_history_size = -1
      end if

      !---------------------------------------------------------------------------
      ! Allocate and initialize big-blockmatrix structure
      !---------------------------------------------------------------------------
      ! Which blocks are active
      a11 = .false.; a12 = .true.
      a21 = .true. ; a22 = .false.
      if (use_diagonal_blocks) then
         a11 = .true.; a22 = .true.
      end if

      ! Single particle external field (F_ij*(a^\dagger a) --> m11)
      ! Represents.. : [ f_pn, 0 ; 0, 0 ], for 2bc w/ pairing: [f11/2, f12; -f12*/2, -f11*]
      ! Stores...... : [ f_pn, _ ; _, _ ], for 2bc w/ pairing: [f11/2, f12;  f12/2,   f11 ]
      Fsp%m11 = f%mat
      ! THIS DOES NOT WORK B/C IT SPOILS THE 1 BLOCK PER ROW+COL STRUCTURE IN FQP
      if (allocated(f%mat12%elem)) then
         Fsp%m12 = f%mat12
         Fsp%m21 = f%mat12
         Fsp%m22 = f%mat
         Fsp%m11%elem = f%mat%elem*0.5_dp
         Fsp%m22%elem = f%mat%elem*0.5_dp
         call set_sign_bbm(Fsp, 1d0,1d0,-1d0,-1d0)
      end if
      if (allocated(g)) then
         if (allocated(Gsp)) deallocate(Gsp)
         allocate(Gsp(size(g)), Gqp(size(g)))
         do ixterm=1, size(g)
            Gsp(ixterm)%m11 = g(ixterm)%mat
         end do
      end if

      ! Bogoliubov transform (real matrices)
      ! Represents.. : [ U, V* ; V, U* ]
      ! Stores...... : [ U, V  ; V, U  ]
      Wn%m11 = Un; Wn%m12 = Vn
      Wn%m21 = Vn; Wn%m22 = Un
      Wp%m11 = Up; Wp%m12 = Vp
      Wp%m21 = Vp; Wp%m22 = Up

      ! External fields in QP basis
      ! Represents.. : [ F11, F20 ; -F02, -F11t ]
      ! Stores...... : [ F11, F20 ;  F02,  F11t ]
      call allocate_bbm(Fqp, a11,a12,a21,a22, nxy)
      call set_sign_bbm(Fqp, 1d0,1d0,-1d0,-1d0)
      if (allocated(g)) then
         do ixterm=1, size(g)
            call allocate_bbm(Gqp(ixterm), a11,a12,a21,a22, nxy)
            call set_sign_bbm(Gqp(ixterm), 1d0,1d0,-1d0,-1d0)
         end do
      end if

      ! Perturbed hamiltonian in QP basis
      ! Represents.. : [ H11, H20 ; -H02, -H11t ]
      ! Stores...... : [ H11, H20 ;  H02,  H11t ]
      call allocate_bbm(dHqp_re, a11,a12,a21,a22, nxy)
      call allocate_bbm(dHqp_im, a11,a12,a21,a22, nxy); dHqp_im%imag=.true.
      call set_sign_bbm(dHqp_re, 1d0,1d0,-1d0,-1d0)
      call set_sign_bbm(dHqp_im, 1d0,1d0,-1d0,-1d0)

      ! Perturbed generalized density in QP basis
      ! Represents.. : [ P, X ; -Y, -Q ]
      ! Stores...... : [ P, X ;  Y,  Q ]
      call allocate_bbm(dRqp_re, a11,a12,a21,a22, nxy)
      call allocate_bbm(dRqp_im, a11,a12,a21,a22, nxy); dRqp_im%imag=.true.
      call set_sign_bbm(dRqp_re, 1d0,1d0,-1d0,-1d0)
      call set_sign_bbm(dRqp_im, 1d0,1d0,-1d0,-1d0)

      ! FAM greens function: g = -1/(E1 +/- E2 +/- (re(eqrpa)+im(eqrpa)*i))
      ! Represents.. : diag( gP, gX, gY, gQ )
      ! Stores...... : [ gP, gX ; gY, gQ ]
      call allocate_bbm(Greens_re, a11,a12,a21,a22, nxy)
      call allocate_bbm(Greens_im, a11,a12,a21,a22, nxy); Greens_im%imag=.true.

      ! Statistical occupation matrix: 1ff=(1-f1-f2) and ff=(f2-f1)
      ! Represents.. : diag( T1ff Tff Tff T1ff )
      ! Stores...... : [ T1ff, Tff ; Tff, T1ff ]
      if (use_diagonal_blocks) call allocate_bbm(T, a11,a12,a21,a22, nxy)

      ! Perturbed hamiltonian in SP basis
      ! Represents.. : [ h_pn, Delta(+)_pn ; -Delta(-)*_pn, -(h^T)_pn ]
      ! Stores...... : [ h_pn, Delta(+)_pn ;  Delta(-) _pn,   h_np    ]
      call allocate_bbm(dHsp_re, .true.,.true.,.true.,.true., nxy)
      call allocate_bbm(dHsp_im, .true.,.true.,.true.,.true., nxy); dHsp_im%imag=.true.
      call  set_sign_bbm(dHsp_re,+1d0,+1d0,-1d0,-1d0)
      call set_trans_bbm(dHsp_re,'n', 'n', 'n', 't')
      call  set_sign_bbm(dHsp_im,+1d0,+1d0,+1d0,-1d0)
      call set_trans_bbm(dHsp_im,'n', 'n', 'n', 't')

      ! Perturbed density in SP basis (Note kap is antisymmetric so np<-->pn gives a minus sign)
      ! Represents.. : [ rho_pn,  kap(+)_pn ; -kap(-)*_pn, -(rho^T)_pn ]
      ! Stores...... : [ rho_pn,  kap(+)_np ;  kap(-) _np,   rho_np    ]
      !            ?=? [ rho_pn, rhotilde(+)_pn ; rhotilde(-)_pn, rho_np ]
      call allocate_bbm(dRsp_re, .true.,.true.,.true.,.true., nxy)
      call allocate_bbm(dRsp_im, .true.,.true.,.true.,.true., nxy); dRsp_im%imag=.true.
      call  set_sign_bbm(dRsp_re,+1d0,-1d0,+1d0,-1d0)
      call set_trans_bbm(dRsp_re,'n', 't', 't', 't')
      call  set_sign_bbm(dRsp_im,+1d0,-1d0,-1d0,-1d0)
      call set_trans_bbm(dRsp_im,'n', 't', 't', 't')

      !---------------------------------------------------------------------------
      ! Transform external field to QP basis
      !---------------------------------------------------------------------------
      if (bminus) then
         call triprod_bbm('t',Wp,'n',Fsp,'n',Wn,Fqp)
      else
         call triprod_bbm('t',Wn,'n',Fsp,'n',Wp,Fqp)
      end if
      if (allocated(g)) then
         do ixterm=1, size(g)
            if (bminus) then
               call triprod_bbm('t',Wp,'n',Gsp(ixterm),'n',Wn,Gqp(ixterm))
            else
               call triprod_bbm('t',Wn,'n',Gsp(ixterm),'n',Wp,Gqp(ixterm))
            end if
         end do
      end if

      ! DMR
      if (debug_tr_permutation) then
         call debug_write_tr_permutation(Fqp%m12, Fqp%m21)
      end if
      if (undo_tr_permutation) then
         call build_tr_permutation(Fqp%m12, Fqp%m21)
      else
         call clear_tr_permutation()
      end if
      ! DMR

      !---------------------------------------------------------------------------
      ! Initialize individual block structures based on external field
      ! (For arrays whose block structure is not determined by routine triprod)
      !---------------------------------------------------------------------------
      ! dRqp has same structure as Fqp
      call copy_block_structure_bbm(Fqp, dRqp_re)
      call copy_block_structure_bbm(Fqp, dRqp_im)

      ! dHqp has same structure as Fqp
      call copy_block_structure_bbm(Fqp, dHqp_re)
      call copy_block_structure_bbm(Fqp, dHqp_im)

      ! Greens has same structure as dRqp
      call copy_block_structure_bbm(dRqp_re, Greens_re)
      call copy_block_structure_bbm(dRqp_im, Greens_im)

      ! Statmatrix has same structure as Fqp
      if (use_diagonal_blocks) call copy_block_structure_bbm(Fqp, T)

      ! dh_pn has same structure as Fsp
      call copy_block_structure(Fsp%m11, dHsp_re%m11)
      call copy_block_structure(Fsp%m11, dHsp_im%m11)
      ! dDelta(+)_pn has same structure as X
      call copy_block_structure(dRqp_re%m12, dHsp_re%m12)
      call copy_block_structure(dRqp_im%m12, dHsp_im%m12)
      ! dDelta(-)_pn has same structure as Y
      call copy_block_structure(dRqp_re%m21, dHsp_re%m21)
      call copy_block_structure(dRqp_im%m21, dHsp_im%m21)
      ! dh_np has same structure as Fsp
      call copy_block_structure(Fsp%m11, dHsp_re%m22)
      call copy_block_structure(Fsp%m11, dHsp_im%m22)

      !---------------------------------------------------------------------------
      ! Calculate 2QP quantities
      !---------------------------------------------------------------------------
      ! Greens function: -1/(E1 +/- E2 +/- (re(eqrpa)+im(eqrpa)*i))
      ! Calculate the 2QP matrix elements M_{ij} = a*( b*f1_i + c*f2_j + e )**d
      
      ! DMR
      ! Changed the first value to take on Greens_sign value, which is set by preprocessor conditionals
      !
      if (bminus) then
         call matrix_2qp(Greens_sign,+1,+1,Greens_power,-real_eqrpa,-imag_eqrpa,.true., Ep,En,Greens_re%m12)
         call matrix_2qp(Greens_sign,+1,+1,Greens_power,-real_eqrpa,-imag_eqrpa,.false.,Ep,En,Greens_im%m12)
         call matrix_2qp(Greens_sign,+1,+1,Greens_power,+real_eqrpa,+imag_eqrpa,.true., Ep,En,Greens_re%m21)
         call matrix_2qp(Greens_sign,+1,+1,Greens_power,+real_eqrpa,+imag_eqrpa,.false.,Ep,En,Greens_im%m21)
         if (use_diagonal_blocks) then
         call matrix_2qp(Greens_sign,+1,-1,Greens_power,-real_eqrpa,-imag_eqrpa,.true., Ep,En,Greens_re%m11)
         call matrix_2qp(Greens_sign,+1,-1,Greens_power,-real_eqrpa,-imag_eqrpa,.false.,Ep,En,Greens_im%m11)
         call matrix_2qp(Greens_sign,+1,-1,Greens_power,+real_eqrpa,+imag_eqrpa,.true., Ep,En,Greens_re%m22)
         call matrix_2qp(Greens_sign,+1,-1,Greens_power,+real_eqrpa,+imag_eqrpa,.false.,Ep,En,Greens_im%m22)
         end if
      else
         call matrix_2qp(-1,+1,+1,-1,-real_eqrpa,-imag_eqrpa,.true., En,Ep,Greens_re%m12)
         call matrix_2qp(-1,+1,+1,-1,-real_eqrpa,-imag_eqrpa,.false.,En,Ep,Greens_im%m12)
         call matrix_2qp(-1,+1,+1,-1,+real_eqrpa,+imag_eqrpa,.true., En,Ep,Greens_re%m21)
         call matrix_2qp(-1,+1,+1,-1,+real_eqrpa,+imag_eqrpa,.false.,En,Ep,Greens_im%m21)
         if (use_diagonal_blocks) then
         call matrix_2qp(-1,+1,-1,-1,-real_eqrpa,-imag_eqrpa,.true., En,Ep,Greens_re%m11)
         call matrix_2qp(-1,+1,-1,-1,-real_eqrpa,-imag_eqrpa,.false.,En,Ep,Greens_im%m11)
         call matrix_2qp(-1,+1,-1,-1,+real_eqrpa,+imag_eqrpa,.true., En,Ep,Greens_re%m22)
         call matrix_2qp(-1,+1,-1,-1,+real_eqrpa,+imag_eqrpa,.false.,En,Ep,Greens_im%m22)
         end if
      end if

      ! Statistical occupation matrix: (1-f1-f2) and (f2-f1)
      if (use_diagonal_blocks) then
         call set_val_bbm(T,0d0)
         if (bminus) then
            call matrix_2qp(+1,-1,-1,+1, 1.0_dp,0.0_dp, .true., qp_fp,qp_fn, T%m12)
            call matrix_2qp(+1,-1,-1,+1, 1.0_dp,0.0_dp, .true., qp_fp,qp_fn, T%m21)
            call matrix_2qp(+1,-1,+1,+1, 0.0_dp,0.0_dp, .true., qp_fp,qp_fn, T%m11)
            call matrix_2qp(+1,-1,+1,+1, 0.0_dp,0.0_dp, .true., qp_fp,qp_fn, T%m22)
         else
            call matrix_2qp(+1,-1,-1,+1, 1.0_dp,0.0_dp, .true., qp_fn,qp_fp, T%m12)
            call matrix_2qp(+1,-1,-1,+1, 1.0_dp,0.0_dp, .true., qp_fn,qp_fp, T%m21)
            call matrix_2qp(+1,-1,+1,+1, 0.0_dp,0.0_dp, .true., qp_fn,qp_fp, T%m11)
            call matrix_2qp(+1,-1,+1,+1, 0.0_dp,0.0_dp, .true., qp_fn,qp_fp, T%m22)
         end if
      end if

   end subroutine init_pnfam_solver


   subroutine deallocate_pnfam_solver
      implicit none
      integer :: ixterm

      ! External fields in QP basis
      call deallocate_bbm(Fqp)
      if (allocated(g)) then
         do ixterm=1, size(g)
            call deallocate_bbm(Gqp(ixterm))
         end do
      end if

      ! Perturbed hamiltonian in QP basis
      call deallocate_bbm(dHqp_re)
      call deallocate_bbm(dHqp_im)

      ! Perturbed generalized density in QP basis
      call deallocate_bbm(dRqp_re)
      call deallocate_bbm(dRqp_im)

      ! FAM greens function: g = -1/(E1 +/- E2 +/- (re(eqrpa)+im(eqrpa)*i))
      call deallocate_bbm(Greens_re)
      call deallocate_bbm(Greens_im)

      ! Statistical occupation matrix: 1ff=(1-f1-f2) and ff=(f2-f1)
      if (use_diagonal_blocks) call deallocate_bbm(T)

      ! Perturbed hamiltonian in SP basis
      call deallocate_bbm(dHsp_re)
      call deallocate_bbm(dHsp_im)

      ! Perturbed density in SP basis (Note kap is antisymmetric so np<-->pn gives a minus sign)
      call deallocate_bbm(dRsp_re)
      call deallocate_bbm(dRsp_im)

   end subroutine deallocate_pnfam_solver

   !----------------------------------------------------------------------------
   ! Calculate the 2QP matrix elements M_{ij} = a*( b*f1_i + c*f2_j + e )**d
   ! Input arrays are assumed size(f) = dqp (e.g. Ep, En, qp_fp, qp_fn)
   ! - Greens(X): a= -1, b= +1, c= +1, d= -1, e= -(w+hw*i), f1= Ep/En,       f2= En/Ep
   ! - Greens(Y): a= -1, b= +1, c= +1, d= -1, e= +(w+hw*i), f1= Ep/En,       f2= En/Ep
   ! - Greens(P): a= -1, b= +1, c= -1, d= -1, e= -(w+hw*i), f1= Ep/En,       f2= En/Ep
   ! - Greens(Q): a= -1, b= +1, c= -1, d= -1, e= +(w+hw*i), f1= Ep/En,       f2= En/Ep
   ! - T(X/Y)   : a= +1, b= -1, c= -1, d= +1, e= +1,        f1= qp_fp/qp_fn, f2= qp_fn/qp_fp
   ! - T(P/Q)   : a= +1, b= +1, c= -1, d= +1, e=  0,        f1= qp_fp/qp_fn, f2= qp_fn/qp_fp
   !----------------------------------------------------------------------------
   subroutine matrix_2qp(a,b,c,d,ree,ime,re,f1,f2,mat)

      implicit none

      integer, intent(in) :: a,b,c,d
      real(dp), intent(in) :: f1(:), f2(:), ree, ime
      logical, intent(in) :: re
      type(blockmatrix), intent(inout) :: mat

      logical  :: l1, l2
      integer  :: ibr, ibc, ipt, i1, i2, i1a, i1b, i2a, i2b
      complex(dp) :: aux

      mat%elem = 0

      ipt = 1
      do ibr=1, nb
         ibc = mat%ir2c(ibr)
         if (ibc == 0) cycle
            i1a = isstart(ibr) ; i1b = i1a + db(ibr) - 1
            i2a = isstart(ibc) ; i2b = i2a + db(ibc) - 1
            do i2 = i2a,i2b ! true column index
            do i1 = i1a,i1b ! true row index
               aux = a*(b*f1(i1) + c*f2(i2) + cmplx(ree,ime,dp))**d
               if (re) then
                  mat%elem(ipt) = real(aux,dp)
               else
                  mat%elem(ipt) = aimag(aux)
               end if
               ipt = ipt + 1
            end do
            end do
      end do

   end subroutine matrix_2qp


   !----------------------------------------------------------------------------
   ! Does matrix-vector multiplication:
   ! Y_ = H * X_
   ! where H is
   ! (|A  B | - omega|I  0|)
   ! (|B* A*|        |0 -I|)
   !----------------------------------------------------------------------------
   subroutine matvec(X_re, X_im, Y_re, Y_im)
      implicit none

      type(bigblockmatrix), intent(in) :: X_re, X_im
      type(bigblockmatrix), intent(inout) :: Y_re, Y_im
      type(bigblockmatrix) :: Xin_re, Xin_im
      type(bigblockmatrix) :: Yout_re, Yout_im

      Xin_re = X_re
      Xin_im = X_im
      if (undo_tr_permutation) then
         call apply_Pt_to_m21(Xin_re)
         call apply_Pt_to_m21(Xin_im)
      end if
      
      if (abs(quench_residual_int) < 1e-10) then
         call set_val_bbm(dHqp_im,0d0)
         call set_val_bbm(dHqp_re,0d0)
      else
         ! Compute the perturbed densities in SP space
         if (bminus) then
            call triprod_bbm('n',Wp,'n',Xin_re,'t',Wn,dRsp_re)
            call triprod_bbm('n',Wp,'n',Xin_im,'t',Wn,dRsp_im)
         else
            call triprod_bbm('n',Wn,'n',Xin_re,'t',Wp,dRsp_re)
            call triprod_bbm('n',Wn,'n',Xin_im,'t',Wp,dRsp_im)
         end if
         call calc_hamiltonian(&
                                     dRsp_re%m11, dRsp_im%m11, dRsp_re%m12, dRsp_im%m12,&
                                     dRsp_re%m22, dRsp_im%m22, dRsp_re%m21, dRsp_im%m21,&
                                     dHsp_re%m11, dHsp_im%m11, dHsp_re%m12, dHsp_im%m12,&
                                     dHsp_re%m22, dHsp_im%m22, dHsp_re%m21, dHsp_im%m21 &
                              )
         ! Transform the Hamiltonian into quasiparticle basis
         if (bminus) then
            call triprod_bbm('t',Wp,'n',dHsp_re,'n',Wn,dHqp_re)
            call triprod_bbm('t',Wp,'n',dHsp_im,'n',Wn,dHqp_im)
         else
            call triprod_bbm('t',Wn,'n',dHsp_re,'n',Wp,dHqp_re)
            call triprod_bbm('t',Wn,'n',dHsp_im,'n',Wp,dHqp_im)
         end if
         ! Quench the residual interaction
         call scalar_mult_bbm(dHqp_re,quench_residual_int)
         call scalar_mult_bbm(dHqp_im,quench_residual_int)
      end if ! End computation of the residual interaction
      Yout_re = Xin_re
      Yout_im = Xin_im
      call complex_emult_bbm(Greens_re, Greens_im, Xin_re, Xin_im, Yout_re, Yout_im)
      call add_bbm(dHqp_re, Yout_re)
      call add_bbm(dHqp_im, Yout_im)

      if (undo_tr_permutation) then
         call apply_P_to_m21(Yout_re)
         call apply_P_to_m21(Yout_im)
      end if

      Y_re = Yout_re
      Y_im = Yout_im
   end subroutine matvec


   subroutine clear_tr_permutation()
      implicit none
      if (allocated(perm_xy)) deallocate(perm_xy)
      if (allocated(invperm_yx)) deallocate(invperm_yx)
      if (allocated(sgn_xy)) deallocate(sgn_xy)
      if (allocated(sgn_yx)) deallocate(sgn_yx)
   end subroutine clear_tr_permutation

   subroutine build_tr_permutation(m12, m21)
      use type_blockmatrix, only : blockmatrix, nb, db, isstart
      use hfb_solution, only : nl, ns
      implicit none

      type(blockmatrix), intent(in) :: m12, m21

      integer :: dqp, half
      integer :: ibr, ibc, nd_row, nd_col, row_local, col_local
      integer :: base12, idx12, row_state, col_state
      integer :: row_tr, col_tr, idx21
      integer :: i, ib
      integer, allocatable :: state_block(:)

      call clear_tr_permutation()

      dqp = sum(db)
      if (mod(dqp, 2) /= 0) return
      half = dqp/2

      allocate(state_block(dqp))
      state_block(:) = 0
      do ib=1, nb
         do i=isstart(ib), isstart(ib)+db(ib)-1
            state_block(i) = ib
         end do
      end do

      allocate(perm_xy(size(m12%elem)))
      allocate(sgn_xy(size(m12%elem)))
      allocate(invperm_yx(size(m21%elem)))
      allocate(sgn_yx(size(m21%elem)))
      perm_xy(:) = 0
      sgn_xy(:) = 0
      invperm_yx(:) = 0
      sgn_yx(:) = 0

      do ibr=1, nb
         ibc = m12%ir2c(ibr)
         if (ibc == 0) cycle
         base12 = m12%ir2m(ibr)
         nd_row = db(ibr)
         nd_col = db(ibc)

         do col_local=1, nd_col
            do row_local=1, nd_row
               idx12 = base12 + (col_local-1)*nd_row + row_local - 1
               row_state = isstart(ibr) + row_local - 1
               col_state = isstart(ibc) + col_local - 1

               row_tr = merge(row_state-half, row_state+half, row_state > half)
               col_tr = merge(col_state-half, col_state+half, col_state > half)

               idx21 = elem_index_local(m21, state_block, row_tr, col_tr)
               perm_xy(idx12) = idx21
               sgn_xy(idx12)  = - state_phase_local(row_state) * state_phase_local(col_state)
            end do
         end do
      end do

      do i=1, size(perm_xy)
         if (perm_xy(i) <= 0) cycle
         invperm_yx(perm_xy(i)) = i
         sgn_yx(perm_xy(i)) = sgn_xy(i)
      end do

      deallocate(state_block)

   contains

      integer function state_phase_local(i_state) result(phi)
         implicit none
         integer, intent(in) :: i_state
         if (i_state < 1 .or. i_state > size(nl)) then
            phi = 1
         else
            phi = merge(-1, 1, (2*nl(i_state) + ns(i_state)) < 0)
         end if
      end function state_phase_local

      integer function elem_index_local(mat, state_block_local, irow, icol) result(idx)
         use type_blockmatrix, only : blockmatrix, db, isstart
         implicit none
         type(blockmatrix), intent(in) :: mat
         integer, intent(in) :: state_block_local(:)
         integer, intent(in) :: irow, icol
         integer :: br, bc, base, rloc, cloc, ndr

         idx = 0
         if (irow < 1 .or. irow > size(state_block_local)) return
         if (icol < 1 .or. icol > size(state_block_local)) return

         br = state_block_local(irow)
         bc = state_block_local(icol)
         if (br == 0 .or. bc == 0) return
         if (mat%ir2c(br) /= bc) return

         base = mat%ir2m(br)
         ndr  = db(br)
         rloc = irow - isstart(br) + 1
         cloc = icol - isstart(bc) + 1
         idx  = base + (cloc-1)*ndr + rloc - 1
      end function elem_index_local

   end subroutine build_tr_permutation

   subroutine apply_P_to_m21(x)
      implicit none
      type(bigblockmatrix), intent(inout) :: x
      real(dp), allocatable :: tmp(:)
      integer :: i

      if (.not. allocated(x%m21%elem)) return
      if (.not. allocated(perm_xy)) return
      if (use_diagonal_blocks) then
         call writelog("WARNING: undo_tr_permutation only implemented for off-diagonal blocks (no FT/EFA diagonal blocks).")
         return
      end if

      allocate(tmp(size(x%m21%elem)))
      do i=1, size(perm_xy)
         if (perm_xy(i) > 0) then
            tmp(i) = sgn_xy(i) * x%m21%elem(perm_xy(i))
         else
            tmp(i) = 0.0_dp
         end if
      end do
      x%m21%elem = tmp
      deallocate(tmp)
   end subroutine apply_P_to_m21

   subroutine apply_Pt_to_m21(x)
      implicit none
      type(bigblockmatrix), intent(inout) :: x
      real(dp), allocatable :: tmp(:)
      integer :: j

      if (.not. allocated(x%m21%elem)) return
      if (.not. allocated(invperm_yx)) return
      if (use_diagonal_blocks) then
         call writelog("WARNING: undo_tr_permutation only implemented for off-diagonal blocks (no FT/EFA diagonal blocks).")
         return
      end if

      allocate(tmp(size(x%m21%elem)))
      do j=1, size(invperm_yx)
         if (invperm_yx(j) > 0) then
            tmp(j) = sgn_yx(j) * x%m21%elem(invperm_yx(j))
         else
            tmp(j) = 0.0_dp
         end if
      end do
      x%m21%elem = tmp
      deallocate(tmp)
   end subroutine apply_Pt_to_m21


   !----------------------------------------------------------------------------
   ! Debug helper: infer the time-reversal permutation/sign matrix P that maps the
   ! ordering of the "X" sector (m12) to the ordering of the "Y" sector (m21).
   !
   ! Writes a file with one non-zero entry per row of P:
   !   i12  i21  sgn  row  col  rowTR  colTR
   ! where i12 is an index into m12%elem and i21 is the corresponding index into
   ! m21%elem, and sgn is the ±1 time-reversal phase for the two-q.p. pair.
   !----------------------------------------------------------------------------
   subroutine debug_write_tr_permutation(m12, m21)
      use type_blockmatrix, only : blockmatrix, nb, db, isstart
      use hfb_solution, only : nl, ns
      implicit none

      type(blockmatrix), intent(in) :: m12, m21

      integer :: dqp, half
      integer :: ibr, ibc, nd_row, nd_col, row_local, col_local
      integer :: base12, idx12, row_state, col_state
      integer :: row_tr, col_tr, idx21
      integer :: i, ib, lu
      integer :: missing, duplicate_targets, nfilled
      integer, allocatable :: state_block(:)
      integer, allocatable :: perm(:), sgn(:), seen(:)
      character(len=256) :: fn

      dqp = sum(db)
      if (mod(dqp, 2) /= 0) then
         call writelog("WARNING: debug_write_tr_permutation: sum(db) is odd; skipping.")
         return
      end if
      half = dqp/2

      allocate(state_block(dqp))
      state_block(:) = 0
      do ib=1, nb
         do i=isstart(ib), isstart(ib)+db(ib)-1
            state_block(i) = ib
         end do
      end do

      allocate(perm(size(m12%elem)))
      allocate(sgn(size(m12%elem)))
      perm(:) = 0
      sgn(:) = 0

      nfilled = 0
      do ibr=1, nb
         ibc = m12%ir2c(ibr)
         if (ibc == 0) cycle
         base12 = m12%ir2m(ibr)
         nd_row = db(ibr)
         nd_col = db(ibc)

         do col_local=1, nd_col
            do row_local=1, nd_row
               idx12 = base12 + (col_local-1)*nd_row + row_local - 1
               row_state = isstart(ibr) + row_local - 1
               col_state = isstart(ibc) + col_local - 1

               row_tr = tr_state(row_state)
               col_tr = tr_state(col_state)
               idx21  = elem_index(m21, row_tr, col_tr)

               perm(idx12) = idx21
               ! Extra minus comes from the convention used throughout pnFAM bigblockmatrices:
               ! the physical (02) sector corresponds to the stored m21 quadrant with sign -1.
               sgn(idx12)  = - state_phase(row_state) * state_phase(col_state)
               nfilled = nfilled + 1
            end do
         end do
      end do

      missing = count(perm == 0)
      allocate(seen(size(m21%elem)))
      seen(:) = 0
      do i=1, size(perm)
         if (perm(i) > 0 .and. perm(i) <= size(m21%elem)) then
            seen(perm(i)) = seen(perm(i)) + 1
         end if
      end do
      duplicate_targets = count(seen > 1)

      write(fn,'("tr_permutation_",a,"_K",i0,".dat")') trim(operator_name), operator_k
      open(newunit=lu, file=trim(fn), status='replace', action='write')
      write(lu,'(a)') "# Time-reversal permutation/sign map between m12 and m21 element orderings"
      write(lu,'(a,1x,a)') "# operator_name=", trim(operator_name)
      write(lu,'(a,1x,i0)') "# operator_k=", operator_k
      write(lu,'(a,1x,i0)') "# sum(db)=", dqp
      write(lu,'(a,1x,i0)') "# nxy(m12)=", size(m12%elem)
      write(lu,'(a,1x,i0)') "# nxy(m21)=", size(m21%elem)
      write(lu,'(a,1x,i0)') "# filled=", nfilled
      write(lu,'(a,1x,i0)') "# missing=", missing
      write(lu,'(a,1x,i0)') "# duplicate_targets=", duplicate_targets
      write(lu,'(a)') "# qp_tr assumes contiguous TR partner (±sum(db)/2)"
      write(lu,'(a)') "# Columns: i12  i21  sgn  row  col  rowTR  colTR"

      do ibr=1, nb
         ibc = m12%ir2c(ibr)
         if (ibc == 0) cycle
         base12 = m12%ir2m(ibr)
         nd_row = db(ibr)
         nd_col = db(ibc)

         do col_local=1, nd_col
            do row_local=1, nd_row
               idx12 = base12 + (col_local-1)*nd_row + row_local - 1
               row_state = isstart(ibr) + row_local - 1
               col_state = isstart(ibc) + col_local - 1
               row_tr = tr_state(row_state)
               col_tr = tr_state(col_state)
               write(lu,'(2(1x,i0),1x,i0,4(1x,i0))') idx12, perm(idx12), sgn(idx12), &
                  row_state, col_state, row_tr, col_tr
            end do
         end do
      end do
      close(lu)

      write(st,'("  *  Wrote TR permutation file: ",a)') trim(fn); call writelog(st)
      if (missing /= 0 .or. duplicate_targets /= 0) then
         write(st,'("  *  TR map diagnostics: missing=",i0,", duplicate_targets=",i0)') missing, duplicate_targets
         call writelog(st)
      end if

      deallocate(state_block, perm, sgn, seen)

   contains

      integer function tr_state(i_state) result(i_tr)
         implicit none
         integer, intent(in) :: i_state

         if (i_state < 1 .or. i_state > dqp) then
            i_tr = 0
            return
         end if

         if (i_state > half) then
            i_tr = i_state - half
         else
            i_tr = i_state + half
         end if
      end function tr_state

      integer function state_phase(i_state) result(phi)
         implicit none
         integer, intent(in) :: i_state
         if (i_state < 1 .or. i_state > size(nl)) then
            phi = 1
         else
            phi = merge(-1, 1, (2*nl(i_state) + ns(i_state)) < 0)
         end if
      end function state_phase

      integer function elem_index(mat, irow, icol) result(idx)
         implicit none
         type(blockmatrix), intent(in) :: mat
         integer, intent(in) :: irow, icol
         integer :: br, bc, base, rloc, cloc, ndr

         idx = 0
         if (irow < 1 .or. irow > dqp) return
         if (icol < 1 .or. icol > dqp) return

         br = state_block(irow)
         bc = state_block(icol)
         if (br == 0 .or. bc == 0) return
         if (mat%ir2c(br) /= bc) return

         base = mat%ir2m(br)
         ndr  = db(br)
         rloc = irow - isstart(br) + 1
         cloc = icol - isstart(bc) + 1
         idx  = base + (cloc-1)*ndr + rloc - 1
      end function elem_index

   end subroutine debug_write_tr_permutation

   !----------------------------------------------
   ! Get the external field objects
   ! (depends on namelist and basis wfs)
   !
   ! init_external_field:
   ! sets: label, beta_minus, parity, K, rank(J)
   ! calls: init_fam_mapping(depends on blockmatrix, K, parity, Kpi-of-blocks)
   ! calls: ext_field_operator (calcs operator, depends on basis, block structure)
   !
   ! setup_cross_terms just calls init_ext_field for the various cross terms
   !----------------------------------------------
   subroutine setup_extfield()
      use pnfam_constants, only : get_timer
      use pnfam_extfield, only : init_external_field, setup_crossterms, set_use_2bc
      use pnfam_storage, only : fam_io
      use pnfam_extfield_2bc, only : effective_2bc_extfield, debug
      use type_extfield_2bc, only : init_extfield_2bc_type
      use type_gamdel_2bc, only : calc_totgamdel
      implicit none
      integer :: iexit, i
      real(dp), allocatable, dimension(:) :: elem_tmp
      real(dp) :: time1, time2

      iexit=1

      ! Calculate 1 body external field
      !11/5/23: since the rest of the code for the two body current mode is designed to handle the GT part, I won't mess with that and 
      !instead just carve out an exception for the P operator with two body current. also only run the two body part if operator = GT.
      !11/29: Add in two body current mode to setup_crossterms. If it's 0 then nothing should change, since set_use_2bc does nothing if it's 0.
      if ((operator_name == 'P' .or. operator_name == 'PS0' .or. operator_name == 'RS0' .or. operator_name == 'RS1' .or. operator_name == 'RS2' .or. operator_name == 'RS0I' .or. operator_name == 'RS1I') .and. two_body_current_mode /= 0) then
         call init_external_field(beta_type=beta_type, label=operator_name, k=operator_k, op=f, use_2bc = two_body_current_mode)
      else
         call init_external_field(beta_type=beta_type, label=operator_name, k=operator_k, op=f)
      end if
      nxterms = 0
      if (compute_crossterms) then
         call setup_crossterms(op=f, crossterms=g, n=nxterms, use_2bc = two_body_current_mode)
      end if

      ! For production runs (no debug) we only calculate/read the 2BC Yukawa part.
      ! The 1BC and 2BC contact part are calculated using the 1BC code and added in.
      ! use_2bc code is:
      ! * 1st digit: 1=1BC+2BC, 2=2BC only
      ! * 2nd digit: 1=Full FAM, 2=SNM+LDA, 3=ASNM+LDA, 4=DME(exc,dir=0), 5=DME(exc)+FAM(dir)
      ! * 3rd digit: 1=Gamma, 2=Delta, 3=Gamma+Delta
      if (two_body_current_mode /= 0 .and. operator_name == 'GT') then
         ! Set the 2BC part of the external field type
         call set_use_2bc(f, two_body_current_mode)
         if (f%use_2bc(4) == 0) return !not calculating 2bc for the GT operator.
         ! Store the 1-body part (if doing 1BC+2BC)
         if (allocated(elem_tmp)) deallocate(elem_tmp)
         allocate(elem_tmp(size(f%mat%elem)))
         elem_tmp = 0
         if (f%use_2bc(1) == 1) then
            elem_tmp = -f%mat%elem ! Park Sign Convention: -sigma tau + 2BC
         end if

         ! Calculate 2-body contact part using (sigma tau rho)  (if not debug).
         ! Also calculates LDA and DME contributions (which have the same contact part as FAM)
         if (abs(debug)<=1) then
            call init_external_field(beta_type=beta_type, label=operator_name, k=operator_k, op=f,&
                    use_2bc=two_body_current_mode)
            do i=1,size(f%mat%elem)
               elem_tmp(i) = elem_tmp(i) + f%mat%elem(i)
            end do
 
            ! If mode needs no FAM calc, we are done: F = -sigma tau + sigma tau f(rho)
            if (f%use_2bc(2) /= 1 .and. f%use_2bc(2) /= 5) then
               f%mat%elem = elem_tmp
               return
            end if
         end if

         ! Finally, read in the Yukawa part or calculate and store it
         f%mat%elem = 0
         call fam_io(-1,iexit)
         if (iexit > 0) then 
            write(st,'(A)') "  * Calculating 2BC matrix elements..."
            call writelog(st)

            time1 = get_timer()
            call effective_2bc_extfield(f)
            time2 = get_timer()

            write(st,'(A, ES12.3)') "  * Two-body current wall time (min): ", (time2-time1)/60.0_dp
            call writelog(st)
            call writelog("")

            ! Use only the direct part of FAM with DME
            if (f%use_2bc(2) == 5) then
                call calc_totgamdel(f%gam, f%mat%elem, 'd')
            end if

            call fam_io(-2,iexit) ! deallocates f%gam/del
         end if

         ! Add 2-body + 1-body and/or contact part
         do i=1,size(f%mat%elem)
            f%mat%elem(i) = f%mat%elem(i) + elem_tmp(i)
         end do

         ! Pairing only - NB I have not checked the contact parts contribution to pairing is correct
         if (f%use_2bc(3) == 2) then
            f%mat12%elem = f%mat%elem
            call deallocate_blockmatrix(f%mat)
         end if

      end if

   end subroutine setup_extfield


   ! NB: Fsp, dRqp, Fqp, imag_eqrpa are all globals in this module
   subroutine transition_density
      use hfb_solution
      use pnfam_constants, only : pi, IT_NEUTRON, IT_PROTON
      use type_blockmatrix
      use type_extfield_2bc, only : init_extfield_2bc_type, caux, cd
      use pnfam_extfield, only : tbc_nmlda_da1, dme_exc
      implicit none
      type(bigblockmatrix) :: dRsp_QRPA, dRqp_QRPA, dR_tmp
      real(dp) :: prefactor, N, Rij, Fij, Sij, amplitude, integral_td1, integral_td2
 
      integer :: ipt, ibx1, ibx2, nd1, nd2, i2, ix2, i1, ix1
      integer :: xl1, xs1, xl2, xs2

      real(dp) :: td(nghl), wf_1(nghl), wf_2(nghl)
      real(dp) :: rhop(nghl), rhon(nghl), rho_fac(nghl)
      character(len=1200) :: denfile, stde

      ! Redundant, but init the couplings
      call init_extfield_2bc_type(.false.)

      ! Transition ampltiude from QP basis solution:
      !   * -1/pi Im(S) = 1/(hw*pi) |<i|F|0>|^2
      !   * |<i|F|0>|^2 = -hw Im(S)
      call contract_bbm(Fqp, dRqp_im, amplitude)
      amplitude = -imag_eqrpa*amplitude ! +/- |<i| F |0>|^2
      amplitude = sqrt(abs(amplitude))  !     |<i| F |0>|

      ! Assume QRPA matrices are real, then:
      !   * Im(X_FAM) = -1/hw (X_QRPA) <i|F|0>
      !   * X_QRPA = (-hw/<i|F|0>) Im(X_FAM)
      dRqp_QRPA = dRqp_im    ! Make a copy
      dRqp_QRPA%imag=.false. ! Declare real
      prefactor = -imag_eqrpa/amplitude
      call scalar_mult_bbm(dRqp_QRPA, prefactor)

      ! Check normalization to see if we are at an eigenvalue:
      !   * N = \sum_{uv} |X_{uv}|^2 - |Y_{uv}|^2
      dR_tmp = dRqp_QRPA
      if (allocated(dR_tmp%m21%elem)) dR_tmp%m21%elem = -dR_tmp%m21%elem
      if (allocated(dR_tmp%m22%elem)) dR_tmp%m22%elem = -dR_tmp%m22%elem
      call contract_bbm(dRqp_QRPA, dR_tmp, N)

      ! Convert X/Y_QRPA to single particle basis
      dRsp_QRPA = dRsp_im     ! Make a copy
      dRsp_QRPA%imag = .false.! Declare real
      if (bminus) then
         call triprod_bbm('n',Wp,'n',dRqp_QRPA,'t',Wn, dRsp_QRPA)
      else
         call triprod_bbm('n',Wn,'n',dRqp_QRPA,'t',Wp, dRsp_QRPA)
      end if

      ! Calculate density of the transition amplitude
      !   * Since Fsp = (f11, 0; 0, 0), only need to contract with drho_QRPA
      !   * TD = \sum_{ij} rho_{ij} f11_{ij} phi_{i} phi_{j} ??? Also what about spin and lambda?
      ipt = 0; td  = 0
      do ibx1 = 1, nb
         ibx2 = dRsp_QRPA%m11%ir2c(ibx1)
         if (ibx2 == 0) cycle
         nd1 = db(ibx1) ; nd2 = db(ibx2)
         if (nd1 == 0 .or. nd2 == 0) cycle

         do i2=1, nd2
            ! State |2>
            ix2 = i2 + isstart(ibx2) - 1
            wf_2 = wf(:,ix2)
            xl2 = nl(ix2);  xs2 = ns(ix2)

            do i1=1, nd1
               ! State <1|
               ix1 = i1 + isstart(ibx1) - 1
               wf_1 = wf(:,ix1)
               xl1 = nl(ix1);  xs1 = ns(ix1)

               ! Matrix element
               ipt = ipt + 1

               ! Technically non-zero in td, but won't contribute to amplitude
               if (xl1 /= xl2) cycle

               ! Spin matrix element
               select case(f%K)
                  case (0)
                     if (xs1 == xs2) then
                        Sij = xs1
                     end if
                  case (1,-1)
                     if (xs1 == xs2 + 2*f%K) then
                        Sij = -f%K*sqrt(2.0_dp)
                     end if
               end select

               Rij = dRsp_QRPA%m11%elem(ipt)
               td = td + Rij*wf_1*wf_2*Sij

            end do ! i1
         end do ! i2
      end do ! ibx

      ! Get densities and NM density function
      call hfb_density_coord(IT_NEUTRON,rhon)
      call hfb_density_coord(IT_PROTON, rhop)

      rho_fac = (caux*2*cd)*(rhon + rhop)
      if (f%use_2bc(2) == 2) then
         ! Symmetric nuclear matter calculation at P=0
         rho_fac = rho_fac + tbc_nmlda_da1(rhon, rhop, 0.0_dp, .true.)
      else if (f%use_2bc(2) == 3) then
         ! Asymmetric nuclear matter calculation at P=0
         rho_fac = rho_fac + tbc_nmlda_da1(rhon, rhop, 0.0_dp, .false.)
      else if (f%use_2bc(2) == 4) then
         ! DME exchange term
         rho_fac = rho_fac + dme_exc()
      else
         rho_fac = 0
      end if

      integral_td1 = sum(td)
      integral_td2 = sum(td*rho_fac)

      ! Strip integration weights
      td = td*wdcori

      write(denfile,'(a,".den")') trim(fam_output_filename)
      call openlog(denfile, .false.)
      call writelog(" QRPA Transition Amplitude Densities")

      write(stde, '(1x,a20,es20.10e3)') "Re(EQRPA)         = ", real_eqrpa ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "Im(EQRPA)         = ", imag_eqrpa ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "QRPA Norm         = ", N ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "Amplitude         = ", amplitude ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "Int. trans        = ", integral_td1 ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "Int. trans*f(rho) = ", integral_td2 ; call writelog(stde)

      write(stde, '(6a20)') "R", "Z", "rhon", "rhop", "trans", "frho_2bc" ; call writelog(stde, .true.)
      do i1 = 1, nghl
         write(stde, '(6ES20.10E3)') 1/y(i1), z(i1), rhon(i1), rhop(i1), td(i1), rho_fac(i1)
         call writelog(stde, .true.)
      end do
      call closelog

   end subroutine transition_density

end module pnfam_solver
