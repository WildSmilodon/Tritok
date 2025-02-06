C^L           
C----------------------------------------------------------------------C
c----------------------------------------------------------------------c
c                                                                      c
c     SOLVER PROGRAM BLOCK                                             c
c     Version: 1.0                                                     c
c                                                                      c
c----------------------------------------------------------------------c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  SUBROUTINES USED :                                              **
c **                                                                  **
c **                   - SolvEQN                                      **
c **                   - FormPRE                                      **
c **                   - SolvEQNfm                                    **
c **                   - FormPREfm                                    **
c **                   - SolvSLE                                      **
c **                   - FormPRM                                      **
c **                   - FormPRMjr                                    **
C **                   - slvlsqr2                                     **
c **                                                                  **
c **********************************************************************
C^L           
C----------------------------------------------------------------------C
c
      SUBROUTINE SolvEQN (slvtg,pretg,prep,neq,nnz,
     &                    iq,jq,dq,q,ia,ja,da,a,b,x,nit,ierr)
c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Solve system of equations using various solvers                  **
c ** ----            --     -                                         **
c **                                                                  **
c ** slvtg .. solver type                                             **
c **          0 direct                                                **
c **          1 cgs                                                   **
c **          2 rbicgstab                                             **
c **          3 rgmres                                                **
c **          4 lsqr                                                  **
c ** pretg .. preconditioner Q type                                   **
c **          0 none                                                  **
c **          1 diag                                                  **
c **          2 ilu                                                   **
c ** prep ... preconditioner Q place of calculation                   **
c **          1 inside                                                **
c **          2 outside                                               **
c ** neq,nnz ... number of equations, number of nonzeros              **
c ** iq,jq,dq,q ... preconditioner matrix Q                           **
c ** ia,ja,da,a ... matrix A                                          **
c ** b  ..... right hand side vector                                  **
c ** x  ..... vector of unknowns                                      **
c ** nit .... number of iterations                                    **
c ** ierr ... error status (0=no error)                               **
c **                                                                  **
c **********************************************************************
c Arguments
      INTEGER   slvtg,pretg,prep,neq,nnz,nit,ierr,
     &          iq(neq+1),jq(nnz),dq(neq),ia(neq+1),ja(nnz),da(neq)
      REAL*8    q(nnz),a(nnz),b(neq),x(neq)
c Internal
      LOGICAL   lexist,info
      INTEGER   slvt,pret,maxit,stopt
      REAL*8    eps,dnrm2,pdnrm2
      EXTERNAL  matvec,prenon,predia,preilu,pdsum,pdnrm2,progress
c
      DATA      maxit,stopt,eps,info/500, 5, 1.D-6, .FALSE./
c
c Read initial values
c
      slvt=slvtg
      pret=pretg
      INQUIRE (FILE='bem.slv',EXIST=lexist)
      IF (lexist) THEN
        OPEN (99,FILE='bem.slv',STATUS='OLD')
        READ (99,*,IOSTAT=ierr) slvt,pret,maxit,stopt,eps
        IF (ierr.EQ.0) GOTO 10
        REWIND(99)
        READ (99,*) slvt
        READ (99,*) pret
        READ (99,*) maxit
        READ (99,*) stopt
        READ (99,*) eps
        READ (99,*) info
 10     CLOSE (99)
        IF (info) THEN
          WRITE (*,'(A,I2,A,I2,A,I5,A,I2,A,E8.2)') 
     &      ' SolvEQN: Solver =',slvt,' Precon =',pret,
     &      ' Maxit =',maxit,' Stop =',stopt,' Eps =',REAL(eps)
          info= .FALSE.
        END IF
      END IF
c
c Check for b=0
c
      IF (DNRM2(neq,b,1).EQ.0.D0) THEN
        CALL dinit(neq,0.0D0,x,1)
        nit=0
        ierr=0
        RETURN
      END IF
c
c Calculate preconditioner
c
      IF (prep.EQ.1.AND.pret.NE.0) THEN
        CALL dcopy(nnz,a,1,q,1)
        IF (pret.EQ.1) THEN
c diagonal
          CALL frmdia(neq,nnz,1,iq,jq,dq,q)
        ELSE IF (pret.EQ.2) THEN
c ilu
          CALL frmilu(neq,nnz,iq,jq,dq,q)
        ELSE
c error
          WRITE (*,*) 'Error: SolvEQN: wrong preconditioner type'
          STOP
        END IF
      END IF
c
c Run solver
c
      IF (slvt.EQ.0) THEN
c direct
        CALL solveysmp(neq,nnz,ia,ja,a,b,x,ierr)
        nit=1
      ELSE IF (slvt.EQ.1) THEN
c CGS
        IF (pret.EQ.0) THEN
          CALL pimcgs(neq,nnz,pret,maxit,stopt,eps,nit,ierr,
     &                ia,ja,da,q,ia,ja,a,b,x,
     &                matvec,prenon,prenon,pdsum,pdnrm2,progress)
        ELSE IF (pret.EQ.1) THEN
          CALL pimcgs(neq,nnz,pret,maxit,stopt,eps,nit,ierr,
     &                ia,ja,da,q,ia,ja,a,b,x,
     &                matvec,predia,predia,pdsum,pdnrm2,progress)
        ELSE IF (pret.EQ.2) THEN
          CALL pimcgs(neq,nnz,pret,maxit,stopt,eps,nit,ierr,
     &                ia,ja,da,q,ia,ja,a,b,x,
     &                matvec,preilu,preilu,pdsum,pdnrm2,progress)
        END IF
      ELSE IF (slvt.EQ.2) THEN
c RBi-CGSTAB
        IF (pret.EQ.0) THEN
          CALL pimrbicgstab(neq,nnz,pret,maxit,stopt,eps,nit,ierr,
     &                      ia,ja,da,q,ia,ja,a,b,x,
     &                      matvec,prenon,prenon,pdsum,pdnrm2,progress)
        ELSE IF (pret.EQ.1) THEN
          CALL pimrbicgstab(neq,nnz,pret,maxit,stopt,eps,nit,ierr,
     &                      ia,ja,da,q,ia,ja,a,b,x,
     &                      matvec,predia,predia,pdsum,pdnrm2,progress)
        ELSE IF (pret.EQ.2) THEN
          CALL pimrbicgstab(neq,nnz,pret,maxit,stopt,eps,nit,ierr,
     &                      ia,ja,da,q,ia,ja,a,b,x,
     &                      matvec,preilu,preilu,pdsum,pdnrm2,progress)
        END IF
      ELSE IF (slvt.EQ.3) THEN
c RGMRES
        IF (pret.EQ.0) THEN
          CALL pimrgmres(neq,nnz,pret,maxit,stopt,eps,nit,ierr,
     &                      ia,ja,da,q,ia,ja,a,b,x,
     &                      matvec,prenon,prenon,pdsum,pdnrm2,progress)
        ELSE IF (pret.EQ.1) THEN
          CALL pimrgmres(neq,nnz,pret,maxit,stopt,eps,nit,ierr,
     &                      ia,ja,da,q,ia,ja,a,b,x,
     &                      matvec,predia,predia,pdsum,pdnrm2,progress)
        ELSE IF (pret.EQ.2) THEN
          CALL pimrgmres(neq,nnz,pret,maxit,stopt,eps,nit,ierr,
     &                      ia,ja,da,q,ia,ja,a,b,x,
     &                      matvec,preilu,preilu,pdsum,pdnrm2,progress)
        END IF
      ELSE IF (slvt.EQ.4) THEN
c LSQR
          CALL slvlsqr(neq,neq,nnz,pret,maxit,stopt,eps,nit,ierr,
     &                      ia,ja,da,q,ia,ja,a,b,x)
      ELSE
c error
        WRITE(*,*) 'SolvEQN: wrong solver type'
        STOP
      END IF
c
      IF (ierr.EQ.-1) THEN
        WRITE(*,*) 'Warning: SolvEQN: solver reached maximum number of iterations.'
      ELSE IF (ierr.LE.-2) THEN
        WRITE(*,*) 'Warning: SolvEQN: solver finished with error',ierr,'.'
      END IF
c
      RETURN
      END
C^L           
C----------------------------------------------------------------------C
c
      SUBROUTINE FormPRE (pretg,neq,nnz,
     &                    iq,jq,dq,q,ia,ja,da,a,ierr)
c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Form Preconditioning Matrix Q                                    **
c ** ---- ---                                                         **
c **                                                                  **
c ** pretg .. preconditioner Q type                                   **
c **          0 none                                                  **
c **          1 diag                                                  **
c **          2 ilu                                                   **
c ** neq,nnz ... number of equations, number of nonzeros              **
c ** iq,jq,dq,q ... preconditioner matrix Q                           **
c ** ia,ja,da,a ... matrix A                                          **
c ** ierr ... error status (0=no error)                               **
c **                                                                  **
c **********************************************************************
c Arguments
      INTEGER   pretg,neq,nnz,ierr,
     &          iq(neq+1),jq(nnz),dq(neq),ia(neq+1),ja(nnz),da(neq)
      REAL*8    q(nnz),a(nnz)
c Internal
      LOGICAL   lexist,info
      INTEGER   pret,dummy
      EXTERNAL  prenon,predia,preilu
      DATA      info/.FALSE./
c
c Read initial values
c
      pret=pretg
      INQUIRE (FILE='bem.slv',EXIST=lexist)
      IF (lexist) THEN
        OPEN (99,FILE='bem.slv',STATUS='OLD')
        READ (99,*,IOSTAT=ierr) dummy,pret
        IF (ierr.EQ.0) GOTO 10
        REWIND(99)
        READ (99,*) dummy
        READ (99,*) pret
        READ (99,*) dummy
        READ (99,*) dummy
        READ (99,*) dummy
        READ (99,*) info
 10     CLOSE (99)
        IF (info) THEN
          WRITE (*,'(A,I2)') ' FormPRE: Precon =',pret
          info= .FALSE.
        END IF
      END IF
c
c Calculate preconditioner
c
      IF (pret.NE.0) THEN
c copy A to Q
        CALL dcopy(nnz,a,1,q,1)
        IF (pret.EQ.1) THEN
c diagonal
          CALL frmdia(neq,nnz,1,iq,jq,dq,q)
        ELSE IF (pret.EQ.2) THEN
c ilu
          CALL frmilu(neq,nnz,iq,jq,dq,q)
        ELSE
c error
          WRITE (*,*) 'Error: FormPRE: wrong preconditioner type'
          STOP
        END IF
      END IF
c
      RETURN
      END
C^L           
C----------------------------------------------------------------------C
c
      SUBROUTINE SolvEQNfm(slvt,pret,prep,maxit,stopt,eps,neq,
     &                     q,a,b,x,nit,cpu,ierr)
c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Solve system of equations using various solvers (full matrix)    **
c ** ----            --     -                         -    -          **
c **                                                                  **
c ** slvt ... solver type                                             **
c **          0 direct                                                **
c **          1 cgs                                                   **
c **          2 rbicgstab                                             **
c **          3 rgmres                                                **
c ** pret ... preconditioner Q type                                   **
c **          0 none                                                  **
c **          1 diag                                                  **
c **          2 lu                                                    **
c ** prep ... preconditioner Q place of calculation                   **
c **          1 inside                                                **
c **          2 outside                                               **
c ** neq  ... number of equations                                     **
c ** maxit .. maximum number of iterations                            **
c ** stopt .. stopping criterum type                                  **
c ** eps  ... stopping criterium value                                **
c ** q  ..... preconditioner vector Q                                 **
c ** a  ..... matrix A                                                **
c ** b  ..... right hand side vector                                  **
c ** x  ..... vector of unknowns                                      **
c ** nit .... number of iterations                                    **
c ** cpu .... elapsed CPU time                                        **
c ** ierr ... error status (0=no error)                               **
c **                                                                  **
c **********************************************************************
c Arguments
      INTEGER   slvt,pret,prep,maxit,stopt,neq,nit,ierr
      REAL      eps,cpu
      REAL*8    q(neq),a(neq,neq),b(neq),x(neq)
c Internal
      REAL      cpu0,cptime
      REAL*8    dnrm2,pdnrm2
      EXTERNAL  matvecfm,prenon,prediafm,pdsum,pdnrm2,progress
c Dummy
      INTEGER   nnz,ia,ja,da
c
c Measure CPU time
c
      cpu0=cptime(0.)
c
c Check for b=0
c
      IF (DNRM2(neq,b,1).EQ.0.D0) THEN
        CALL dinit(neq,0.0D0,x,1)
        nit=0
        ierr=0
        RETURN
      END IF
c
c Calculate preconditioner
c
      IF (pret.NE.0 .AND. prep.EQ.1) THEN
        IF (pret.EQ.1) THEN
c diagonal
          CALL frmdiafm(neq,1,a,q)
        ELSE IF (pret.EQ.2) THEN
c lu
          CALL dgefa(a,neq,neq,q,ierr)
        ELSE
c error
          CALL WarnErrSolv(2,'SolvEQNfm','Wrong preconditioner type!',0)
        END IF
      END IF
c
c Run solver
c
      IF (slvt.EQ.0) THEN
c Direct solver
        x=b
        CALL dgesl(a,neq,neq,q,x,0)
      ELSE IF (slvt.EQ.1) THEN
c CGS
        IF (pret.EQ.0) THEN
          CALL pimcgs(neq,nnz,pret,maxit,stopt,DBLE(eps),nit,ierr,
     &                ia,ja,da,q,ia,ja,a,b,x,
     &                matvecfm,prenon,prenon,pdsum,pdnrm2,progress)
        ELSE IF (pret.EQ.1) THEN
          CALL pimcgs(neq,nnz,pret,maxit,stopt,DBLE(eps),nit,ierr,
     &                ia,ja,da,q,ia,ja,a,b,x,
     &                matvecfm,prediafm,prediafm,pdsum,pdnrm2,progress)
        ELSE
          CALL WarnErrSolv(2,'SolvEQNfm','Wrong preconditioner type!',0)
        END IF
      ELSE IF (slvt.EQ.2) THEN
c RBi-CGSTAB
        IF (pret.EQ.0) THEN
          CALL pimrbicgstab(neq,nnz,pret,maxit,stopt,DBLE(eps),nit,ierr,
     &                      ia,ja,da,q,ia,ja,a,b,x,
     &                      matvecfm,prenon,prenon,pdsum,pdnrm2,progress)
        ELSE IF (pret.EQ.1) THEN
          CALL pimrbicgstab(neq,nnz,pret,maxit,stopt,DBLE(eps),nit,ierr,
     &                      ia,ja,da,q,ia,ja,a,b,x,
     &                      matvecfm,prediafm,prediafm,pdsum,pdnrm2,progress)
        ELSE
          CALL WarnErrSolv(2,'SolvEQNfm','Wrong preconditioner type!',0)
        END IF
      ELSE IF (slvt.EQ.3) THEN
c RGMRES
        IF (pret.EQ.0) THEN
          CALL pimrgmres(neq,nnz,pret,maxit,stopt,DBLE(eps),nit,ierr,
     &                      ia,ja,da,q,ia,ja,a,b,x,
     &                      matvecfm,prenon,prenon,pdsum,pdnrm2,progress)
        ELSE IF (pret.EQ.1) THEN
          CALL pimrgmres(neq,nnz,pret,maxit,stopt,DBLE(eps),nit,ierr,
     &                      ia,ja,da,q,ia,ja,a,b,x,
     &                      matvecfm,prediafm,prediafm,pdsum,pdnrm2,progress)
        ELSE
          CALL WarnErrSolv(2,'SolvEQNfm','Wrong preconditioner type!',0)
        END IF
      ELSE
c error
        CALL WarnErrSolv(2,'SolvEQNfm','Wrong solver type!',0)
      END IF
c
      IF (ierr.EQ.-1) THEN
        CALL WarnErrSolv(1,'SolvEQNfm','Solver reached maximum number of iterations!',0)
c     ELSE IF (ierr.LE.-2) THEN
c       WRITE(*,*) 'Warning: SolvSLE: solver finished with error',ierr,'.'
      END IF
c
      cpu=cpu+cptime(cpu0)
c
      RETURN
      END
C^L           
C----------------------------------------------------------------------C
c
      SUBROUTINE FormPREfm(pret,neq,q,a,cpu,ierr)
c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Form Preconditioning Matrix Q (full matrix A)                    **
c ** ---- ---                       -    -                            **
c **                                                                  **
c ** pret ... preconditioner Q type                                   **
c **          0 none                                                  **
c **          1 diag                                                  **
c **          2 lu                                                    **
c ** neq  ... number of equations                                     **
c ** q    ... preconditioner matrix Q                                 **
c ** a    ... matrix A                                                **
c ** cpu .... elapsed CPU time                                        **
c ** ierr ... error status (0=no error)                               **
c **                                                                  **
c **********************************************************************
c Arguments
      INTEGER pret,neq,ierr
      REAL    cpu
      REAL*8  q(neq),a(neq,neq)
c Internal
      REAL    cpu0,cptime
c
c Measure CPU time
c
      cpu0=cptime(0.)
c
c Calculate preconditioner
c
      IF (pret.NE.0) THEN
        IF (pret.EQ.1) THEN
c diagonal
          CALL frmdiafm(neq,1,a,q)
        ELSE IF (pret.EQ.2) THEN
c lu
          CALL dgefa(a,neq,neq,q,ierr)
        ELSE
c error
          CALL WarnErrSolv(2,'FormPREfm','Wrong preconditioner type!',0)
        END IF
      END IF
c
      cpu=cpu+cptime(cpu0)
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c
      SUBROUTINE SolvSLE(slvt,pret,prep,maxit,stopt,eps,neq,nnz,
     &                   q,a,b,x,nit,cpu,ierr)
c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Solve system of linear equations using various solvers           **
c ** ----  -         -      -                                         **
c **                                                                  **
c ** slvt ... solver type                                             **
c **          0 direct                                                **
c **          1 cgs                                                   **
c **          2 rbicgstab                                             **
c **          3 rgmres                                                **
c ** pret ... preconditioner Q type                                   **
c **          0 none                                                  **
c **          1 diag                                                  **
c **          2 ilu                                                   **
c ** prep ... preconditioner Q place of calculation                   **
c **          1 inside                                                **
c **          2 outside                                               **
c ** neq  ... number of equations                                     **
c ** maxit .. maximum number of iterations                            **
c ** stopt .. stopping criterum type                                  **
c ** eps  ... stopping criterium value                                **
c ** q  ..... preconditioner vector Q                                 **
c ** a  ..... matrix A                                                **
c ** b  ..... right hand side vector                                  **
c ** x  ..... vector of unknowns                                      **
c ** nit .... number of iterations                                    **
c ** cpu .... elapsed CPU time                                        **
c ** ierr ... error status (0=no error)                               **
c **                                                                  **
c **********************************************************************
c Type definition
      TYPE MATRIX
        INTEGER :: neq,nnz
        INTEGER, POINTER :: i(:),j(:),d(:),ij(:,:)
        REAL(8), POINTER :: v(:)
      END TYPE MATRIX
c Arguments
      INTEGER  slvt,pret,prep,maxit,stopt,neq,nnz,nit,ierr
      REAL     eps,cpu
      TYPE     (MATRIX) q,a
      REAL(8)  b(neq),x(neq)
c Internal
      REAL     cpu0,cptime
      REAL(8)  dnrm2,pdnrm2
      EXTERNAL matvec,prenon,predia,preilu,pdsum,pdnrm2,progress
c
c Measure CPU time
c
      cpu0=cptime(0.)
c
c Check for b=0
c
      IF (DNRM2(neq,b,1).EQ.0.0D0) THEN
        x=0.0D0
        nit=0
        ierr=0
        RETURN
      END IF
c
c Calculate preconditioner
c
      IF (prep.EQ.1.AND.pret.NE.0) THEN
        q%v=a%v
        IF (pret.EQ.1) THEN
c diagonal
          CALL frmdia(neq,nnz,1,q%i,q%j,q%d,q%v)
        ELSE IF (pret.EQ.2) THEN
c ilu
          CALL frmilu(neq,nnz,q%i,q%j,q%d,q%v)
        ELSE
c error
          CALL WarnErrSolv(2,'SolvSLE','Wrong preconditioner type!',0)
        END IF
      END IF
c
c Run solver
c
      IF (slvt.EQ.0) THEN
c Direct solver
        CALL solveysmp(neq,nnz,a%i,a%j,a%v,b,x,ierr)
        nit=1
      ELSE IF (slvt.EQ.1) THEN
c CGS
        IF (pret.EQ.0) THEN
          CALL pimcgs(neq,nnz,pret,maxit,stopt,DBLE(eps),nit,ierr,
     &                q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &                matvec,prenon,prenon,pdsum,pdnrm2,progress)
        ELSE IF (pret.EQ.1) THEN
          CALL pimcgs(neq,nnz,pret,maxit,stopt,DBLE(eps),nit,ierr,
     &                q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &                matvec,predia,predia,pdsum,pdnrm2,progress)
        ELSE IF (pret.EQ.2) THEN
          CALL pimcgs(neq,nnz,pret,maxit,stopt,DBLE(eps),nit,ierr,
     &                q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &                matvec,preilu,preilu,pdsum,pdnrm2,progress)
        ELSE
          CALL WarnErrSolv(3,'SolvSLE','Wrong preconditioner type, using no preconditioner. Pret=',pret)
          CALL pimcgs(neq,nnz,pret,maxit,stopt,DBLE(eps),nit,ierr,
     &                q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &                matvec,prenon,prenon,pdsum,pdnrm2,progress)
        END IF
      ELSE IF (slvt.EQ.2) THEN
c RBi-CGSTAB
        IF (pret.EQ.0) THEN
          CALL pimrbicgstab(neq,nnz,pret,maxit,stopt,DBLE(eps),nit,ierr,
     &                      q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &                      matvec,prenon,prenon,pdsum,pdnrm2,progress)
        ELSE IF (pret.EQ.1) THEN
          CALL pimrbicgstab(neq,nnz,pret,maxit,stopt,DBLE(eps),nit,ierr,
     &                      q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &                      matvec,predia,predia,pdsum,pdnrm2,progress)
        ELSE IF (pret.EQ.2) THEN
          CALL pimrbicgstab(neq,nnz,pret,maxit,stopt,DBLE(eps),nit,ierr,
     &                      q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &                      matvec,preilu,preilu,pdsum,pdnrm2,progress)
        ELSE
          CALL WarnErrSolv(3,'SolvSLE','Wrong preconditioner type, using no preconditioner. Pret=',pret)
          CALL pimrbicgstab(neq,nnz,pret,maxit,stopt,DBLE(eps),nit,ierr,
     &                      q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &                      matvec,prenon,prenon,pdsum,pdnrm2,progress)
        END IF
      ELSE IF (slvt.EQ.3) THEN
c RGMRES
        IF (pret.EQ.0) THEN
          CALL pimrgmres(neq,nnz,pret,maxit,stopt,DBLE(eps),nit,ierr,
     &                      q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &                      matvec,prenon,prenon,pdsum,pdnrm2,progress)
        ELSE IF (pret.EQ.1) THEN
          CALL pimrgmres(neq,nnz,pret,maxit,stopt,DBLE(eps),nit,ierr,
     &                      q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &                      matvec,predia,predia,pdsum,pdnrm2,progress)
        ELSE IF (pret.EQ.2) THEN
          CALL pimrgmres(neq,nnz,pret,maxit,stopt,DBLE(eps),nit,ierr,
     &                      q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &                      matvec,preilu,preilu,pdsum,pdnrm2,progress)
        ELSE
          CALL WarnErrSolv(3,'SolvSLE','Wrong preconditioner type, using no preconditioner. Pret=',pret)
          CALL pimrgmres(neq,nnz,pret,maxit,stopt,DBLE(eps),nit,ierr,
     &                      q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &                      matvec,prenon,prenon,pdsum,pdnrm2,progress)
        END IF
      ELSE IF (slvt.EQ.4) THEN
c LSQR
          CALL slvlsqr(neq,neq,nnz,pret,maxit,stopt,DBLE(eps),nit,ierr,
     &                      q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x)
      ELSE
c Error
        CALL WarnErrSolv(4,'SolvSLE','Wrong solver type!',slvt)
      END IF
c
      IF (ierr.EQ.-1) THEN
        CALL WarnErrSolv(1,'SolvEQN','Solver reached maximum number of iterations!',0)
c     ELSE IF (ierr.LE.-2) THEN
c       WRITE(*,*) 'Warning: SolvSLE: solver finished with error',ierr,'.'
      END IF
c
      cpu=cpu+cptime(cpu0)
c
      RETURN
      END
C^L

C----------------------------------------------------------------------C
c
      SUBROUTINE FormPRMjr(pret,q,a)
c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Form Preconditioning Matrix Q                                    **
c **  ---- --              -                                           **
c **                                                                  **
c **  pret ... preconditioner Q type                                   **
c **          0 none                                                  **
c **          1 diag                                                  **
c **          2 ilu                                                   **
c **  q    ... preconditioner matrix Q                                 **
c **  a    ... matrix A                                                **
c **                                                                  **
c **********************************************************************
c Type definition
      TYPE MATRIX
        INTEGER :: neq,nnz
        INTEGER, POINTER :: i(:),j(:),d(:),ij(:,:)
        REAL(8), POINTER :: v(:)
      END TYPE MATRIX
c Arguments
      INTEGER pret
      TYPE    (MATRIX) q,a
c
c Calculate preconditioner
c       Pred tem mora biti narejeno q=a !!!!!!!!!!!!!!!
c
      IF (pret.EQ.1) THEN
c diagonal
        CALL frmdia(q%neq,q%nnz,1,q%i,q%j,q%d,q%v)
      ELSE IF (pret.EQ.2) THEN
c ilu
        CALL frmilu(q%neq,q%nnz,q%i,q%j,q%d,q%v)
      ELSE
c error
        CALL WarnErrSolv(1,'FormPRM','Wrong preconditioner type!',0)
      END IF
c
      END


C----------------------------------------------------------------------C
c
      SUBROUTINE FormPRM(pret,neq,nnz,q,a,cpu,ierr)
c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Form Preconditioning Matrix Q                                    **
c ** ---- --              -                                           **
c **                                                                  **
c ** pret ... preconditioner Q type                                   **
c **          0 none                                                  **
c **          1 diag                                                  **
c **          2 ilu                                                   **
c ** neq  ... number of equations                                     **
c ** q    ... preconditioner matrix Q                                 **
c ** a    ... matrix A                                                **
c ** cpu .... elapsed CPU time                                        **
c ** ierr ... error status (0=no error)                               **
c **                                                                  **
c **********************************************************************
c Type definition
      TYPE MATRIX
        INTEGER :: neq,nnz
        INTEGER, POINTER :: i(:),j(:),d(:),ij(:,:)
        REAL(8), POINTER :: v(:)
      END TYPE MATRIX
c Arguments
      INTEGER pret,neq,nnz,ierr
      REAL    cpu
      TYPE    (MATRIX) q,a
c Internal
      REAL    cpu0,cptime
c
c Measure CPU time
c
      cpu0=cptime(0.)
c
c Calculate preconditioner
c
      IF (pret.NE.0) THEN
        q%v=a%v
        IF (pret.EQ.1) THEN
c diagonal
          CALL frmdia(neq,nnz,1,q%i,q%j,q%d,q%v)
        ELSE IF (pret.EQ.2) THEN
c ilu
          CALL frmilu(neq,nnz,q%i,q%j,q%d,q%v)
        ELSE
c error
          CALL WarnErrSolv(1,'FormPRM','Wrong preconditioner type!',0)
          ierr=1
        END IF
      END IF
c
      cpu=cpu+cptime(cpu0)
c
      RETURN
      END
C^L           
C----------------------------------------------------------------------C
c----------------------------------------------------------------------c
c                                                                      c
c     GAUSS SOLVER PROGRAM BLOCK                                             c
c                                                                      c
c----------------------------------------------------------------------c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  SUBROUTINES USED :                                              **
c **                                                                  **
c **                   - SolveGAUSS                                   **
c **                   - dgefa                                        **
c **                   - dgesl                                        **
c **                                                                  **
c **********************************************************************
C^L           
C----------------------------------------------------------------------C
c
      SUBROUTINE SolveGAUSS (neq,a,b,x,ierr)

      INTEGER neq,ierr
      REAL*8  a(neq,neq),b(neq),x(neq)
c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Solve system of equations using direct Gaussian ellimination     **
c ** with partial pivoting                                            **
c **                                                                  **
c ** neq  = number of equations                                       **
c ** a    = matrix  A                                                 **
c ** b    = right hand side vector                                    **
c ** x    = solution vector                                           **
c ** ierr = error status                                              **
c **         =0     no errors detected                                **
c **         >0     singular matrix                                   **
c **********************************************************************
      REAL*8 pivot
      ALLOCATABLE pivot(:)
c
c Allocate memory for pivot vector
c
      ALLOCATE(pivot(neq))
c
c Solve system Ax=b
c
      CALL dgefa(a,neq,neq,pivot,ierr)
      IF (ierr.GT.0) STOP 'Error: SolveGAUSS: singular matrix.'
      CALL dcopy(neq,b,1,x,1)
      CALL dgesl(a,neq,neq,pivot,x,0)
c
c Free memory
c
      DEALLOCATE(pivot)
c
      RETURN
      END
C^L           
C----------------------------------------------------------------------C
c
      subroutine dgefa(a,lda,n,ipvt,info)
c
C----------------------------------------------------------------------C
c***begin prologue  dgefa
c***date written   780814   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  d2a1
c***keywords  double precision,factor,linear algebra,linpack,matrix
c***author  moler, c. b., (u. of new mexico)
c***purpose  factors a double precision matrix by gaussian elimination.
c***description
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    double precision(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c***references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
c                 *linpack users  guide*, siam, 1979.
c***routines called  daxpy,dscal,idamax
c***end prologue  dgefa
      integer lda,n,info
      double precision ipvt(*),a(lda,*)
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c     gaussian elimination with partial pivoting
c
c***first executable statement  dgefa
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
C^L           
C----------------------------------------------------------------------C
c
      subroutine dgesl(a,lda,n,ipvt,b,job)
c
C----------------------------------------------------------------------C
c***begin prologue  dgesl
c***date written   780814   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  d2a1
c***keywords  double precision,linear algebra,linpack,matrix,solve
c***author  moler, c. b., (u. of new mexico)
c***purpose  solves the double precision system  a*x=b or  trans(a)*x=b
c            using the factors computed by dgeco or dgefa.
c***description
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    double precision(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c***references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
c                 *linpack users  guide*, siam, 1979.
c***routines called  daxpy,ddot
c***end prologue  dgesl
      integer lda,n,job
      double precision ipvt(*),a(lda,*),b(*)
c
      double precision ddot,t
      integer k,kb,l,nm1
c***first executable statement  dgesl
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
C^L
C----------------------------------------------------------------------C
c                                                                      c      
      SUBROUTINE pimCGS
     &   (neq,nnz,pret,maxit,stopt,eps,nit,ierr,iq,jq,dq,q,ia,ja,a,b,x,
     &    matvec,preconl,preconr,pdsum,pdnrm,progress)
c                                                                      c      
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Solve system of equations using CGS iterative method             **
c ** - --                            ---                              **
c **********************************************************************
c Parameters
      INTEGER  neq,nnz,pret,maxit,stopt,nit,ierr,
     &         iq(neq+1),jq(nnz),dq(neq),
     &         ia(neq+1),ja(nnz),da(neq)
      REAL*8   q(nnz),a(nnz),b(neq),x(neq),eps,pdnrm
      EXTERNAL matvec,preconl,preconr,pdsum,pdnrm,progress
c Internal variables
      INTEGER  nwrk,ipar(13),pres
      REAL*8   wrk,dpar(6)
      EXTERNAL dinit,dvprod,dcgs
c
      ALLOCATABLE wrk(:)
c
c Set initial values
c
c if pret>0 use left preconditioning
      pres=0
      IF (pret.GT.0) pres=1
c
c Allocate memory for matrix and vectors
c
      nwrk = (9)*neq
      ALLOCATE(wrk(nwrk))
      CALL pimdsetpar(ipar,dpar,neq,neq,neq,neq,-1,-1,-1,
     &                pres,stopt,maxit,eps)
c
c CGS
c
      CALL dcgs(neq,nnz,nwrk,iq,jq,dq,q,ia,ja,a,b,x,wrk,ipar,dpar,
     &         matvec,preconl,preconr,pdsum,pdnrm,progress)
c
c Set output values
c
      nit=ipar(11)
      ierr=ipar(12)
c
c Free memory
c
      DEALLOCATE(wrk)
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE dcgs
     &   (neq,nnz,nwrk,iq,jq,dq,q,ia,ja,a,b,x,wrk,ipar,dpar,
     &    matvec,preconl,preconr,pdsum,pdnrm,progress)
c                                                                      c
C----------------------------------------------------------------------C
*     ..
*     .. Parameters ..
      REAL*8  ZERO
      PARAMETER (ZERO=0.0D0)
      REAL*8  ONE
      PARAMETER (ONE=1.0D0)
      INTEGER IPARSIZ
      PARAMETER (IPARSIZ=13)
      INTEGER DPARSIZ
      PARAMETER (DPARSIZ=6)
*     ..
*     .. Matrices A,Q ..
      INTEGER neq,nnz,nwrk,iq(neq+1),jq(nnz),dq(neq),ia(neq+1),ja(nnz)
      REAL*8  q(nnz),a(nnz)
*     ..
*     .. Array Arguments ..
      INTEGER ipar(IPARSIZ)
      REAL*8  b(neq),x(neq),dpar(DPARSIZ),wrk(nwrk)
*     ..
*     .. Function Arguments ..
      REAL*8   PDNRM
      EXTERNAL PDNRM
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL MATVEC,PDSUM,PRECONL,PRECONR,PROGRESS
*     ..
*     .. Local Scalars ..
      INTEGER BASIS,BLKSZ,CNVRTX,IP,IR,IRTILDE,IS,IT,ITNO,IU,IW,
     +        IXOLD,IZ,LDA,LOCLEN,MAXIT,N,NPROCS,PRECONT,PROCID,
     +        STATUS,STEPERR,STOPT
      REAL*8  ALPHA,BETA,EPSILON,EXITNORM,RHO,RHO0,RHSSTOP,XI
*     ..
*     .. Local Arrays ..
      REAL*8  DOTS(1)
*     ..
*     .. External Functions ..
      REAL*8   DDOT,DSETRHSSTOP
      EXTERNAL DDOT,DSETRHSSTOP
*     ..
*     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,PIMDGETPAR,STOPCRIT
*     ..
      CALL PIMDGETPAR(IPAR,DPAR,LDA,N,BLKSZ,LOCLEN,BASIS,NPROCS,
     +                PROCID,PRECONT,STOPT,MAXIT,ITNO,STATUS,
     +                STEPERR,EPSILON,EXITNORM)

*  Check consistency of preconditioning and stop types
      IF (((PRECONT.EQ.0).OR. (PRECONT.EQ.2)) .AND. (STOPT.EQ.6)) THEN
          ITNO = 0
          STATUS = -4
          STEPERR = 0
          GO TO 9999
      END IF

*  Does not need conversion Y=Q2X for residual
      CNVRTX = 0

*  Set indices for mapping local vectors into wrk
      IR = 1
      IRTILDE = IR + LOCLEN
      IP = IRTILDE + LOCLEN
      IS = IP + LOCLEN
      IT = IS + LOCLEN
      IU = IT + LOCLEN
      IW = IU + LOCLEN
      IZ = IW + LOCLEN
      IXOLD = IZ + LOCLEN

*  Set rhs of stopping criteria
      RHSSTOP = DSETRHSSTOP(NEQ,NNZ,IQ,JQ,DQ,Q,B,WRK(IR),EPSILON,IPAR,
     +                      PRECONL,PDNRM)

*  1. r=Q1(b-AQ2x)
      IF (STOPT.NE.6) THEN
          IF (PRECONT.EQ.0) THEN
*     r=b-Ax
              CALL DCOPY(LOCLEN,B,1,WRK(IR),1)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK(IW),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IR),1)

          ELSE IF (PRECONT.EQ.1) THEN
*     r=Q1(b-Ax)
              CALL DCOPY(LOCLEN,B,1,WRK(IZ),1)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK(IW),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IZ),1)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IR),IPAR)

          ELSE IF (PRECONT.EQ.2) THEN
*     r=b-AQ2x
              CALL DCOPY(LOCLEN,B,1,WRK(IR),1)
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK(IW),IPAR)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IZ),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IR),1)

          ELSE IF (PRECONT.EQ.3) THEN
*     r=Q1(b-AQ2x)
              CALL DCOPY(LOCLEN,B,1,WRK(IP),1)
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK(IW),IPAR)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IZ),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IP),1)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IP),WRK(IR),IPAR)
          END IF

      ELSE
*     r has been set to Qb in the call to dsetrhsstop
          IF (PRECONT.EQ.1) THEN
*     r=Q1(b-Ax)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK(IW),IPAR)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IZ),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IR),1)

          ELSE IF (PRECONT.EQ.3) THEN
*     r=Q1(b-AQ2x)
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK(IZ),IPAR)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IZ),WRK(IW),IPAR)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IZ),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IR),1)
          END IF

      END IF

*  2. p=s=rtilde=r
      CALL DCOPY(LOCLEN,WRK(IR),1,WRK(IRTILDE),1)
      CALL DCOPY(LOCLEN,WRK(IR),1,WRK(IP),1)
      CALL DCOPY(LOCLEN,WRK(IR),1,WRK(IS),1)

*  3. rho=dot(rtilde,r)
      DOTS(1) = DDOT(LOCLEN,WRK(IRTILDE),1,WRK(IR),1)
      CALL PDSUM(1,DOTS,IPAR)
      RHO = DOTS(1)

*  Loop
      STATUS = 0
      EXITNORM = -ONE
      STEPERR = -1
      DO 20 ITNO = 1,MAXIT

*  4. w=Q1AQ2p
          IF (PRECONT.EQ.0) THEN
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IP),WRK(IW),IPAR)

          ELSE IF (PRECONT.EQ.1) THEN
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IP),WRK(IZ),IPAR)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IW),IPAR)

          ELSE IF (PRECONT.EQ.2) THEN
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IP),WRK(IZ),IPAR)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IZ),WRK(IW),IPAR)

          ELSE IF (PRECONT.EQ.3) THEN
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IP),WRK(IW),IPAR)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IZ),IPAR)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IW),IPAR)
          END IF

*  5. xi=dot(rtilde,w)
          DOTS(1) = DDOT(LOCLEN,WRK(IRTILDE),1,WRK(IW),1)
          CALL PDSUM(1,DOTS,IPAR)
          XI = DOTS(1)

*  6. alpha=rho/xi
          IF (XI.EQ.ZERO) THEN
              STATUS = -3
              STEPERR = 6
              GO TO 9999

          END IF

          ALPHA = RHO/XI

*  7. t=s-alpha*w
          CALL DCOPY(LOCLEN,WRK(IS),1,WRK(IT),1)
          CALL DAXPY(LOCLEN,-ALPHA,WRK(IW),1,WRK(IT),1)

*  8. w=s+t
          CALL DCOPY(LOCLEN,WRK(IS),1,WRK(IW),1)
          CALL DAXPY(LOCLEN,ONE,WRK(IT),1,WRK(IW),1)

*  9. x=x+alpha*w
          CALL DCOPY(LOCLEN,X,1,WRK(IXOLD),1)
          CALL DAXPY(LOCLEN,ALPHA,WRK(IW),1,X,1)

* 10. r=r-alpha*Q1AQ2w
          IF (PRECONT.EQ.0) THEN
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IU),IPAR)

          ELSE IF (PRECONT.EQ.1) THEN
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IZ),IPAR)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IU),IPAR)

          ELSE IF (PRECONT.EQ.2) THEN
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IZ),IPAR)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IZ),WRK(IU),IPAR)

          ELSE IF (PRECONT.EQ.3) THEN
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IU),IPAR)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IU),WRK(IZ),IPAR)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IU),IPAR)
          END IF

          CALL DAXPY(LOCLEN,-ALPHA,WRK(IU),1,WRK(IR),1)

* 11. check stopping criterion
          CALL STOPCRIT(NEQ,NNZ,IQ,JQ,DQ,Q,IA,JA,A,B,
     +                  WRK(IR),WRK(IZ),X,WRK(IXOLD),WRK(IU),RHSSTOP,
     +                  CNVRTX,EXITNORM,STATUS,IPAR,MATVEC,MATVEC,
     +                  PRECONR,PDSUM,PDNRM)

*  Call monitoring routine
          CALL PROGRESS(LOCLEN,ITNO,EXITNORM,X,WRK(IR),WRK(IZ))

          IF (STATUS.EQ.-5) THEN
              STEPERR = 11
              GO TO 9999
          ELSE IF (STATUS.EQ.0) THEN
              GO TO 9999
          END IF

* 12. rho=dot(rtilde0,r)
          RHO0 = RHO
          DOTS(1) = DDOT(LOCLEN,WRK(IRTILDE),1,WRK(IR),1)
          CALL PDSUM(1,DOTS,IPAR)
          RHO = DOTS(1)

* 13. beta=rho/rho0
          IF (RHO0.EQ.ZERO) THEN
              STATUS = -3
              STEPERR = 13
              GO TO 9999

          END IF

          BETA = RHO/RHO0

* 14. s=r+beta*t
          CALL DCOPY(LOCLEN,WRK(IR),1,WRK(IS),1)
          CALL DAXPY(LOCLEN,BETA,WRK(IT),1,WRK(IS),1)

* 15. w=t+beta*p
          CALL DCOPY(LOCLEN,WRK(IT),1,WRK(IW),1)
          CALL DAXPY(LOCLEN,BETA,WRK(IP),1,WRK(IW),1)

* 16. p=s+beta*w
          CALL DCOPY(LOCLEN,WRK(IS),1,WRK(IP),1)
          CALL DAXPY(LOCLEN,BETA,WRK(IW),1,WRK(IP),1)

   20 CONTINUE

      IF (ITNO.GT.MAXIT) THEN
          STATUS = -1
          ITNO = MAXIT
      END IF

 9999 CONTINUE

      IF ((PRECONT.EQ.2) .OR. (PRECONT.EQ.3)) THEN
          CALL DCOPY(LOCLEN,X,1,WRK(IZ),1)
          CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),X,IPAR)
      END IF

*  Set output parameters
      IPAR(11) = ITNO
      IPAR(12) = STATUS
      IPAR(13) = STEPERR
      DPAR(2) = EXITNORM

      RETURN

      END
C^L
C----------------------------------------------------------------------C
c----------------------------------------------------------------------c
c                                                                      c
c     SOLVER PROGRAM BLOCK                                             c
c                                                                      c
c----------------------------------------------------------------------c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  SUBROUTINES USED :                                              **
c **                                                                  **
c **                   - matvec                                       **
c **                   - matvecfm                                     **
c **                   - prenon                                       **
c **                   - predia                                       **
c **                   - prediafm                                     **
c **                   - preilu                                       **
c **                   - pdsum                                        **
c **                   - pdnrm2                                       **
c **                   - progress                                     **
c **                   - printv                                       **
c **                   - dinit                                        **
c **                   - dvprod                                       **
c **                   - dsetrhsstop                                  **
c **                   - stopcrit                                     **
c **                   - pimdsetpar                                   **
c **                   - pimdgetpar                                   **
c **                                                                  **
c **********************************************************************
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE matvec(neq,nnz,ia,ja,a,u,v,ipar)

      INTEGER neq,nnz,ia(neq+1),ja(nnz),ipar(13)
      REAL*8  a(nnz),u(neq),v(neq)
c                                                                      c
C----------------------------------------------------------------------C
      INTEGER i,ij,n
	REAL*8  sum
c
      n = ipar(2)
      DO i = 1,n
	  sum = 0.0D0
        DO ij = ia(i),ia(i+1) - 1
          sum = sum + a(ij)*u(ja(ij))
        END DO
        v(i) = sum
      END DO
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE matvecfm(neq,nnz,ia,ja,a,u,v,ipar)

      INTEGER neq,nnz,ia(neq+1),ja(nnz),ipar(13)
      REAL*8  a(neq,neq),u(neq),v(neq)
c                                                                      c
C----------------------------------------------------------------------C
      INTEGER i,j,n
c
      n = ipar(2)
      v = 0.0D0
      DO j = 1,n
        DO i = 1,n
          v(i) = v(i) + a(i,j)*u(j)
        END DO
      END DO
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE prenon(neq,nnz,iq,jq,dq,q,u,v,ipar)

      INTEGER neq,nnz,iq(neq+1),jq(nnz),dq(neq),ipar(13)
      REAL*8  q(neq),u(neq),v(neq)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Precondition system with none preconditioner                    **
c **  ---                      ---                                    **
c **********************************************************************
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE predia(neq,nnz,iq,jq,dq,q,u,v,ipar)

      INTEGER neq,nnz,iq(neq+1),jq(nnz),dq(neq),ipar(13)
      REAL*8  q(nnz),u(neq),v(neq)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Precondition system with diagonal preconditioner                **
c **  ---                      ---                                    **
c **********************************************************************
      INTEGER i
c
      DO i=1,neq
        v(i)=q(dq(i))*u(i)
      END DO
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE prediafm(neq,nnz,iq,jq,dq,q,u,v,ipar)

      INTEGER neq,nnz,iq(neq+1),jq(nnz),dq(neq),ipar(13)
      REAL*8  q(neq),u(neq),v(neq)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Precondition system with diagonal preconditioner                **
c **  ---                      ---                                    **
c **********************************************************************
      v=q*u
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE preilu(neq,nnz,iq,jq,dq,q,b,x,ipar)

      INTEGER neq,nnz,iq(neq+1),jq(nnz),dq(neq),ipar(13)
      REAL*8  q(nnz),b(neq),x(neq)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Precondition system with ILU preconditioner                     **
c **  ---                      ---                                    **
c **                                                                  **
c **  Solve Qx=b  (Q=LU; Compressed Row Storage)                      **
c **                                                                  **
c **  neq ... number of equations                                     **
c **  nnz ... number of non zero coefficients in Q                    **
c **  iQ  ... beginning of the rows in Q                              **
c **  jQ  ... column indices of elements in Q                         **
c **  dQ  ... positions of diagonal elements in Q                     **
c **  Q   ... coefficients of LU decomposed matrix Q                  **
c **  b   ... right hand side vector                                  **
c **  x   ... solution vector                                         **
c **                                                                  **
c **********************************************************************
      INTEGER i,j
      REAL*8  sum
c
c  Lx=b; RHS factorization
c
      x(1)=b(1)
      DO i=2,neq
        sum=b(i)
        DO j=iq(i),dq(i)-1
          sum=sum-q(j)*x(jq(j))
        END DO
        x(i)=sum
      END DO
c
c  Ux'=x; Back substitution
c
      x(neq)=x(neq)/q(dq(neq))
      DO i=neq-1,1,-1
        sum=x(i)
        DO j=dq(i)+1,iq(i+1)-1
          sum=sum-q(j)*x(jq(j))
        END DO
        x(i)=sum/q(dq(i))
      END DO
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pdsum(isize,x,ipar)
c                                                                      c
C----------------------------------------------------------------------C
c
*     .. Scalar Arguments ..
      INTEGER ISIZE
*     ..
*     .. Array Arguments ..
      INTEGER ipar(13)
      REAL*8  X(isize)
*     ..
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      REAL*8  FUNCTION pdnrm2(loclen,u,ipar)
c                                                                      c
C----------------------------------------------------------------------C
c
*     .. Scalar Arguments ..
      INTEGER LOCLEN
*     ..
*     .. Array Arguments ..
      INTEGER ipar(13)
      REAL*8  U(loclen)
*     ..
*     .. External Functions ..
      REAL*8  DNRM2
      EXTERNAL DNRM2
*     ..
      PDNRM2 = DNRM2(LOCLEN,U,1)
*     ..
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE progress(loclen,itno,normres,x,res,trueres)
c                                                                      c
C----------------------------------------------------------------------C
c
*     .. Scalar Arguments ..
      INTEGER ITNO,LOCLEN
      REAL*8  NORMRES
*     ..
*     .. Array Arguments ..
      REAL*8  RES(loclen),TRUERES(loclen),X(loclen)
*     ..
*     .. External Subroutines ..
      EXTERNAL PRINTV
*     ..
*     WRITE (6,FMT=9000) ITNO,NORMRES
*     WRITE (6,FMT=9010) 'X:'
*     CALL PRINTV(LOCLEN,X)
*     WRITE (6,FMT=9010) 'RES:'
*     CALL PRINTV(LOCLEN,RES)
*     WRITE (6,FMT=9010) 'TRUE RES:'
*     CALL PRINTV(LOCLEN,TRUERES)
      RETURN
 9000 FORMAT (I5,1X,D16.10)
 9010 FORMAT (A)
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE printv(n,u)
c                                                                      c
C----------------------------------------------------------------------C
c
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER N
*     ..
*     .. Array Arguments ..
      REAL*8  U(n)
*     ..
*     .. Local Scalars ..
      INTEGER I
*     ..
      WRITE (6,FMT=9000) (U(I),I=1,N)
      RETURN

 9000 FORMAT (5(D14.8,1X))
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE dinit(n,alpha,dx,incx)
c                                                                      c
C----------------------------------------------------------------------C
c
*
*     Initialises a vector x with a scalar alpha.
*     Modified from dcopy, BLAS Level 1.
*     Rudnei Dias da Cunha, 14/6/93.
*

*     copies a vector, x, to a vector, y.
*     uses unrolled loops for increments equal to one.
*     jack dongarra, linpack, 3/11/78.
*
*
*     .. Scalar Arguments ..
      REAL*8  ALPHA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      REAL*8  DX(n)
*     ..
*     .. Local Scalars ..
      INTEGER I,IX,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
          DX(IX) = ALPHA
          IX = IX + INCX
   10 CONTINUE
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,7)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DX(I) = ALPHA
   30 CONTINUE
      IF (N.LT.7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
          DX(I) = ALPHA
          DX(I+1) = ALPHA
          DX(I+2) = ALPHA
          DX(I+3) = ALPHA
          DX(I+4) = ALPHA
          DX(I+5) = ALPHA
          DX(I+6) = ALPHA
   50 CONTINUE
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE dvprod(n,dx,incx,dy,incy)
c                                                                      c
C----------------------------------------------------------------------C
*
*     Modified from daxpy level 1 BLAS
*     element-wise vector multiplication, y<-x*y
*     Rudnei Dias da Cunha, 16/6/93
*
*     constant times a vector plus a vector.
*     uses unrolled loops for increments equal to one.
*     jack dongarra, linpack, 3/11/78.
*
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      REAL*8  DX(n),DY(n)
*     ..
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DY(IY) = DY(IY)*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,4)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DY(I) = DY(I)*DX(I)
   30 CONTINUE
      IF (N.LT.4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
          DY(I) = DY(I)*DX(I)
          DY(I+1) = DY(I+1)*DX(I+1)
          DY(I+2) = DY(I+2)*DX(I+2)
          DY(I+3) = DY(I+3)*DX(I+3)
   50 CONTINUE
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      REAL*8 FUNCTION dsetrhsstop
     &          (neq,nnz,iq,jq,dq,q,b,r,epsilon,ipar,preconl,pdnrm)
c                                                                      c
C----------------------------------------------------------------------C
c
*     .. Scalar Arguments ..
      REAL*8  EPSILON
*     ..
*     .. Array Arguments ..
      INTEGER neq,nnz,iq(neq+1),jq(nnz),dq(neq)
      REAL*8  q(neq)
*     ..
*     .. Array Arguments ..
      INTEGER ipar(13)
      REAL*8  b(neq),r(neq)
*     ..
*     .. Function Arguments ..
      REAL*8   PDNRM
      EXTERNAL PDNRM
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL PRECONL
*     ..
*     .. Local Scalars ..
      INTEGER LOCLEN,STOPT
*     ..
      LOCLEN = ipar(4)
      STOPT = ipar(9)
*     ..
      IF ((STOPT.EQ.1) .OR. (STOPT.EQ.4) .OR. (STOPT.EQ.7)) THEN

*  ||r||<epsilon or ||Q1r||<epsilon ||x(k)-x(k-1)||<epsilon
          DSETRHSSTOP = EPSILON

      ELSE IF ((STOPT.EQ.2) .OR. (STOPT.EQ.3) .OR. (STOPT.EQ.5)) THEN

*  ||r||<epsilon||b|| or sqrt(r(Q1r))<epsilon||b|| or ||Q1r||<epsilon||b||
          DSETRHSSTOP = EPSILON*PDNRM(LOCLEN,B,ipar)

      ELSE IF (STOPT.EQ.6) THEN
*  ||Q1r||<epsilon||Q1b||
          CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,B,R,ipar)
          DSETRHSSTOP = EPSILON*PDNRM(LOCLEN,R,ipar)
      END IF

      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE stopcrit
     &   (neq,nnz,iq,jq,dq,q,ia,ja,a,b,r,rtrue,x,xold,wrk,
     &    rhsstop,cnvrtx,exitnorm,status,ipar,
     &    matvec,tmatvec,preconr,pdsum,pdnrm)
c                                                                      c
C----------------------------------------------------------------------C
*     ..
*     .. Arguments ..
      INTEGER  neq,nnz,iq(neq+1),jq(nnz),dq(neq),ia(neq+1),ja(nnz),
     &         cnvrtx,status,ipar(13)
      REAL*8   q(nnz),a(nnz),
     &         b(neq),r(neq),rtrue(neq),x(neq),xold(neq),wrk(neq),
     &         rhsstop,exitnorm,pdnrm
      EXTERNAL MATVEC,TMATVEC,PRECONR,PDSUM,PDNRM
*     ..
*     .. Local ..
      INTEGER  LOCLEN,PRECONT,STOPT
      REAL*8   DOTS(1),DDOT
      EXTERNAL DAXPY,DCOPY,DDOT
      INTRINSIC SQRT
*     ..
*     .. Parameters ..
      REAL*8   ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
*     ..
      LOCLEN = ipar(4)
      PRECONT = ipar(8)
      STOPT = ipar(9)

      IF ((STOPT.EQ.1) .OR. (STOPT.EQ.2) .OR. (STOPT.EQ.3)) THEN

*  Compute true residual if needed
          CALL DCOPY(LOCLEN,B,1,RTRUE,1)

          IF ((PRECONT.EQ.2) .OR. (PRECONT.EQ.3)) THEN
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK,ipar)
              IF (CNVRTX.EQ.1) THEN
*    r=b-AATQ2x
                  CALL TMATVEC(NEQ,NNZ,IA,JA,A,WRK,XOLD,ipar)
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,XOLD,WRK,ipar)
                  CALL DAXPY(LOCLEN,-ONE,WRK,1,RTRUE,1)
              ELSE
*    r=b-AQ2x
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK,XOLD,ipar)
                  CALL DAXPY(LOCLEN,-ONE,XOLD,1,RTRUE,1)
              END IF
          ELSE IF (CNVRTX.EQ.1) THEN
*    r=b-AATx
              CALL TMATVEC(NEQ,NNZ,IA,JA,A,X,XOLD,ipar)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,XOLD,WRK,ipar)
              CALL DAXPY(LOCLEN,-ONE,WRK,1,RTRUE,1)
          ELSE
*    r=b-Ax
              CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK,ipar)
              CALL DAXPY(LOCLEN,-ONE,WRK,1,RTRUE,1)
          END IF
      END IF

      IF ((STOPT.EQ.1) .OR. (STOPT.EQ.2)) THEN

*  ||r||<epsilon or ||r||<epsilon||b||
          EXITNORM = PDNRM(LOCLEN,RTRUE,ipar)
          IF (EXITNORM.LT.RHSSTOP) THEN
              STATUS = 0
          ELSE
              STATUS = -99
          END IF

      ELSE IF (STOPT.EQ.3) THEN

*  sqrt(rT(Q1r))<epsilon||b||
          DOTS(1) = DDOT(LOCLEN,RTRUE,1,R,1)
          CALL PDSUM(1,DOTS(1),ipar)
          IF (DOTS(1).LT.ZERO) THEN
              STATUS = -5
              RETURN
          END IF
          EXITNORM = SQRT(DOTS(1))
          IF (EXITNORM.LT.RHSSTOP) THEN
              STATUS = 0
          ELSE
              STATUS = -99
          END IF

      ELSE IF ((STOPT.EQ.4) .OR. (STOPT.EQ.5) .OR. (STOPT.EQ.6)) THEN

*  ||Q1r||<epsilon or ||Q1r||<epsilon||b|| or ||Q1r||<epsilon||Q1b||
          EXITNORM = PDNRM(LOCLEN,R,ipar)
          IF (EXITNORM.LT.RHSSTOP) THEN
              STATUS = 0
          ELSE
              STATUS = -99
          END IF

      ELSE IF (STOPT.EQ.7) THEN

*  ||x-x0||<epsilon
          CALL DCOPY(LOCLEN,X,1,WRK,1)
          CALL DAXPY(LOCLEN,-ONE,XOLD,1,WRK,1)
          EXITNORM = PDNRM(LOCLEN,WRK,ipar)
          IF (EXITNORM.LT.RHSSTOP) THEN
              STATUS = 0
          ELSE
              STATUS = -99
          END IF
      END IF

      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pimdsetpar
     &   (IPAR,DPAR,LDA,N,BLKSZ,LOCLEN,BASIS,NPROCS,PROCID,
     &    PRECONT,STOPT,MAXIT,EPSILON)
c                                                                      c
C----------------------------------------------------------------------C
C     ..
C     .. Parameters ..
      REAL*8  ONE
      PARAMETER (ONE=1.0D0)
C     ..
C     .. Arguments ..
      INTEGER IPAR(13),BASIS,BLKSZ,LDA,LOCLEN,MAXIT,N,NPROCS,PRECONT,
     &        PROCID,STOPT
      REAL*8  DPAR(6),EPSILON
C     ..
      ipar(1) = LDA
      ipar(2) = N
      ipar(3) = BLKSZ
      ipar(4) = LOCLEN
      ipar(5) = BASIS
      ipar(6) = NPROCS
      ipar(7) = PROCID
      ipar(8) = PRECONT
      ipar(9) = STOPT
      ipar(10) = MAXIT
      ipar(11) = -1
      ipar(12) = -1
      ipar(13) = -1

      dpar(1) = EPSILON
      dpar(2) = -ONE

      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pimdgetpar
     &   (IPAR,DPAR,LDA,N,BLKSZ,LOCLEN,BASIS,NPROCS,PROCID,
     &    PRECONT,STOPT,MAXIT,ITNO,STATUS,STEPERR,EPSILON,EXITNORM)
c                                                                      c
C----------------------------------------------------------------------C
C     ..
C     .. Scalar Arguments ..
      INTEGER IPAR(13),BASIS,BLKSZ,ITNO,LDA,LOCLEN,MAXIT,N,NPROCS,
     &        PRECONT,PROCID,STATUS,STEPERR,STOPT
      REAL*8  DPAR(6),EPSILON,EXITNORM
C     ..
      LDA = ipar(1)
      N = ipar(2)
      BLKSZ = ipar(3)
      LOCLEN = ipar(4)
      BASIS = ipar(5)
      NPROCS = ipar(6)
      PROCID = ipar(7)
      PRECONT = ipar(8)
      STOPT = ipar(9)
      MAXIT = ipar(10)
      ITNO = ipar(11)
      STATUS = ipar(12)
      STEPERR = ipar(13)

      EPSILON = dpar(1)
      EXITNORM = dpar(2)

      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c      
      SUBROUTINE pimRBiCGSTAB
     &   (neq,nnz,pret,maxit,stopt,eps,nit,ierr,iq,jq,dq,q,ia,ja,a,b,x,
     &    matvec,preconl,preconr,pdsum,pdnrm,progress)
c                                                                      c      
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Solve system of equations using RBi-CGSTAB iterative method      **
c ** - --                            --- -----                        **
c **********************************************************************
c Parameters
      INTEGER  neq,nnz,pret,maxit,stopt,nit,ierr,
     &         iq(neq+1),jq(nnz),dq(neq),
     &         ia(neq+1),ja(nnz),da(neq)
      REAL*8   q(nnz),a(nnz),b(neq),x(neq),eps,pdnrm
      EXTERNAL matvec,preconl,preconr,pdsum,pdnrm,progress
c Internal variables
      INTEGER  nwrk,ipar(13),basis,pres
      REAL*8   wrk,dpar(6)
      EXTERNAL dinit,dvprod,rbicgstab
c
      ALLOCATABLE wrk(:)
c
c Set initial values
c
c if pret>0 use left preconditioning
      pres=0
      IF (pret.GT.0) pres=1
c basis dimension for RBi-CGSTAB vectors
      basis=2
c
c Allocate memory for matrix and vectors
c
      nwrk = (6+2*basis)*neq
      ALLOCATE(wrk(nwrk))
      CALL pimdsetpar(ipar,dpar,neq,neq,neq,neq,basis,-1,-1,
     &                pres,stopt,maxit,eps)
c
c RBi-CGSTAB
c
      CALL rbicgstab(neq,nnz,nwrk,iq,jq,dq,q,ia,ja,a,b,x,wrk,ipar,dpar,
     &               matvec,preconl,preconr,pdsum,pdnrm,progress)
c
c Set output values
c
      nit=ipar(11)
      ierr=ipar(12)
c
c Free memory
c
      DEALLOCATE(wrk)
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE rbicgstab
     &   (neq,nnz,nwrk,iq,jq,dq,q,ia,ja,a,b,x,wrk,ipar,dpar,
     &    matvec,preconl,preconr,pdsum,pdnrm,progress)
c                                                                      c
C----------------------------------------------------------------------C
c
*     .. Parameters ..
      REAL*8  ZERO
      PARAMETER (ZERO=0.0D0)
      REAL*8  ONE
      PARAMETER (ONE=1.0D0)
      INTEGER IBDIM
      PARAMETER (IBDIM=8)
*     ..
*     .. Matrices Q,A ..
      INTEGER neq,nnz,nwrk,iq(neq+1),jq(nnz),dq(neq),ia(neq+1),ja(nnz)
      REAL*8  q(nnz),a(nnz)
*     ..
*     .. Array Arguments ..
      INTEGER ipar(13)
      REAL*8  b(neq),x(neq),wrk(nwrk),dpar(6)
*     ..
*     .. Function Arguments ..
      REAL*8   PDNRM
      EXTERNAL PDNRM
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL MATVEC,PDSUM,PRECONL,PRECONR,PROGRESS
*     ..
*     .. Local Scalars ..
      INTEGER BASIS,BLKSZ,CNVRTX,I,I0,I1,I2,I3,I4,IR,IRTILDE,ITNO,IU,
     +        IW,IXOLD,IZ,J,LDA,LOCLEN,MAXIT,N,NPROCS,PRECONT,PROCID,
     +        STATUS,STEPERR,STOPT
      REAL*8  ALPHA,BETA,EPSILON,EXITNORM,KSI,OMEGA,RHO0,RHO1,RHSSTOP,S
*     ..
*     .. Local Arrays ..
      REAL*8  DOTS(IBDIM),GAMMA(IBDIM),GAMMA1(IBDIM),
     +        GAMMA2(IBDIM),SIGMA(IBDIM),TAU(IBDIM,IBDIM)
*     ..
*     .. External Functions ..
      REAL*8   DDOT,DSETRHSSTOP
      EXTERNAL DDOT,DSETRHSSTOP
*     ..
*     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,DINIT,PIMDGETPAR,STOPCRIT
*     ..
      CALL PIMDGETPAR(IPAR,DPAR,LDA,N,BLKSZ,LOCLEN,BASIS,NPROCS,
     +                PROCID,PRECONT,STOPT,MAXIT,ITNO,STATUS,
     +                STEPERR,EPSILON,EXITNORM)

*  Check consistency of preconditioning and stop types
      IF (((PRECONT.EQ.0).OR. (PRECONT.EQ.2)) .AND. (STOPT.EQ.6)) THEN
          ITNO = 0
          STATUS = -4
          STEPERR = 0
          GO TO 9999
      END IF

*  Does not need conversion Y=Q2X for residual
      CNVRTX = 0

*  Set indices for mapping local vectors into wrk
      IRTILDE = 1
      IW = IRTILDE + LOCLEN
      IZ = IW + LOCLEN
      IXOLD = IZ + LOCLEN
      IR = IXOLD + LOCLEN
      IU = IR + (BASIS+1)*LOCLEN

*  Set rhs of stopping criteria
      RHSSTOP = DSETRHSSTOP(NEQ,NNZ,IQ,JQ,DQ,Q,B,WRK(IR),EPSILON,IPAR,
     +                      PRECONL,PDNRM)

*  1. r=Q1(b-AQ2x)
      IF (PRECONT.EQ.0) THEN
*     r=b-Ax
          CALL DCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK(IW),ipar)
          CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IR),1)

      ELSE IF (PRECONT.EQ.1) THEN
*     r=Q1(b-Ax)
          CALL DCOPY(LOCLEN,B,1,WRK(IZ),1)
          CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK(IW),ipar)
          CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IZ),1)
          CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IR),ipar)

      ELSE IF (PRECONT.EQ.2) THEN
*     r=b-AQ2x
          CALL DCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK(IW),ipar)
          CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IZ),ipar)
          CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IR),1)

      ELSE IF (PRECONT.EQ.3) THEN
*     r=Q1(b-AQ2x)
          CALL DCOPY(LOCLEN,B,1,WRK(IW),1)
          CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK(IR),ipar)
          CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IR),WRK(IZ),ipar)
          CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IW),1)
          CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IR),ipar)
      END IF

*  2. rtilde=r
      CALL DCOPY(LOCLEN,WRK(IR),1,WRK(IRTILDE),1)

*  3. u0=0
      CALL DINIT(LOCLEN,ZERO,WRK(IU),1)

*  4. rho0=1, alpha=0, omega=1
      RHO0 = ONE
      ALPHA = ZERO
      OMEGA = ONE

*  Loop
      STATUS = 0
      STEPERR = -1
      EXITNORM = -ONE
      DO 120 ITNO = 1,MAXIT

*  5. rho0=-omega*rho0
          RHO0 = -OMEGA*RHO0

*  BiCG loop
          I1 = 0
          I2 = LOCLEN
          DO 30 J = 0,BASIS - 1

*  6. rho1=r(j)^{T}rtilde
              DOTS(1) = DDOT(LOCLEN,WRK(IR+I1),1,WRK(IRTILDE),1)
              CALL PDSUM(1,DOTS,ipar)
              RHO1 = DOTS(1)

*  7. beta=alpha*rho1/rho0
              IF (RHO0.EQ.ZERO) THEN
                  STATUS = -3
                  STEPERR = 7
                  GO TO 9999

              END IF

              BETA = ALPHA*RHO1/RHO0

*  8. rho0=rho1
              RHO0 = RHO1

*  9. u(i)=r(i)-beta*u(i), i=0:j
              I3 = 0
              DO 10 I = 0,J
                  CALL DCOPY(LOCLEN,WRK(IU+I3),1,WRK(IZ),1)
                  CALL DCOPY(LOCLEN,WRK(IR+I3),1,WRK(IU+I3),1)
                  CALL DAXPY(LOCLEN,-BETA,WRK(IZ),1,WRK(IU+I3),1)
                  I3 = I3 + LOCLEN
   10         CONTINUE

* 10. u(j+1)=Q1AQ2u(j)
              IF (PRECONT.EQ.0) THEN
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IU+I1),WRK(IU+I2),ipar)

              ELSE IF (PRECONT.EQ.1) THEN
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IU+I1),WRK(IW),ipar)
                  CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IU+I2),ipar)

              ELSE IF (PRECONT.EQ.2) THEN
                  CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IU+I1),WRK(IW),ipar)
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IU+I2),ipar)

              ELSE IF (PRECONT.EQ.3) THEN
                  CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IU+I1),WRK(IZ),ipar)
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IZ),WRK(IW),ipar)
                  CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IU+I2),ipar)
              END IF

* 11. ksi=u(j+1)^{T}rtilde
              DOTS(1) = DDOT(LOCLEN,WRK(IU+I2),1,WRK(IRTILDE),1)
              CALL PDSUM(1,DOTS,ipar)
              KSI = DOTS(1)

* 12. alpha=rho0/ksi
              IF (KSI.EQ.ZERO) THEN
                  STATUS = -3
                  STEPERR = 12
                  GO TO 9999
              END IF

              ALPHA = RHO0/KSI

* 13. r(i)=r(i)-alpha*u(i+1), i=0:j
              I3 = 0
              I4 = LOCLEN
              DO 20 I = 0,J
                  CALL DAXPY(LOCLEN,-ALPHA,WRK(IU+I4),1,WRK(IR+I3),1)
                  I3 = I3 + LOCLEN
                  I4 = I4 + LOCLEN
   20         CONTINUE

* 14. r(j+1)=Q1AQ2r(j)
              IF (PRECONT.EQ.0) THEN
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IR+I1),WRK(IR+I2),ipar)

              ELSE IF (PRECONT.EQ.1) THEN
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IR+I1),WRK(IW),ipar)
                  CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IR+I2),ipar)

              ELSE IF (PRECONT.EQ.2) THEN
                  CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IR+I1),WRK(IW),ipar)
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IR+I2),ipar)

              ELSE IF (PRECONT.EQ.3) THEN
                  CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IR+I1),WRK(IZ),ipar)
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IZ),WRK(IW),ipar)
                  CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IR+I2),ipar)
              END IF

* 15. x0=x0+alpha*u0
              CALL DCOPY(LOCLEN,X,1,WRK(IXOLD),1)
              CALL DAXPY(LOCLEN,ALPHA,WRK(IU),1,X,1)
              I1 = I1 + LOCLEN
              I2 = I2 + LOCLEN
   30     CONTINUE

* 16. check stopping criterion
          CALL STOPCRIT(NEQ,NNZ,IQ,JQ,DQ,Q,IA,JA,A,B,
     +                  WRK(IR),WRK(IZ),X,WRK(IXOLD),WRK(IW),
     +                  RHSSTOP,CNVRTX,EXITNORM,STATUS,IPAR,
     +                  MATVEC,MATVEC,PRECONR,PDSUM,PDNRM)

*  Call monitoring routine
          CALL PROGRESS(LOCLEN,ITNO,EXITNORM,X,WRK(IR),WRK(IZ))

          IF (STATUS.EQ.-5) THEN
              STEPERR = 16
              GO TO 9999
          ELSE IF (STATUS.EQ.0) THEN
              GO TO 9999
          END IF
*  MR loop

* 17. sigma(1)=r(1)^{T}r(1), gamma'(1)=r(0)^{T}r(1)/sigma(1)
          DOTS(1) = DDOT(LOCLEN,WRK(IR+LOCLEN),1,WRK(IR+LOCLEN),1)
          DOTS(2) = DDOT(LOCLEN,WRK(IR),1,WRK(IR+LOCLEN),1)
          CALL PDSUM(2,DOTS,ipar)
          SIGMA(1) = DOTS(1)

          IF (SIGMA(1).EQ.ZERO) THEN
              STATUS = -3
              STEPERR = 17
              GO TO 9999
          END IF

          GAMMA1(1) = DOTS(2)/SIGMA(1)

          I0 = LOCLEN + LOCLEN
          DO 60 J = 2,BASIS

* 18. tau(i,j)=r(j)^{T}r(i)/sigma(i), r(j)=r(j)-tau(i,j)r(i)
              I1 = LOCLEN
              DO 40 I = 1,J - 1
                  DOTS(I) = DDOT(LOCLEN,WRK(IR+I0),1,WRK(IR+I1),1)
                  I1 = I1 + LOCLEN
   40         CONTINUE
              CALL PDSUM(J-1,DOTS,ipar)
              I1 = LOCLEN
              DO 50 I = 1,J - 1
                  TAU(I,J) = DOTS(I)/SIGMA(I)
                  CALL DAXPY(LOCLEN,-TAU(I,J),WRK(IR+I1),1,WRK(IR+I0),1)
   50         CONTINUE

* 19. sigma(j)=r(j)^{T}r(j), gamma'(j)=r(0)^{T}r(j)/sigma(j)
              DOTS(1) = DDOT(LOCLEN,WRK(IR+I0),1,WRK(IR+I0),1)
              DOTS(2) = DDOT(LOCLEN,WRK(IR),1,WRK(IR+I0),1)
              CALL PDSUM(2,DOTS,ipar)
              SIGMA(J) = DOTS(1)

              IF (SIGMA(J).EQ.ZERO) THEN
                  STATUS = -3
                  STEPERR = 19
                  GO TO 9999
              END IF

              GAMMA1(J) = DOTS(2)/SIGMA(J)
              I0 = I0 + LOCLEN
   60     CONTINUE

* 20. gamma_{l}=omega=gamma'_{l}
*     gamma_{j}=gamma'_{j}-\sum_{i=j+1}^{l}{tau_{j,i}gamma_{i}}
          GAMMA(BASIS) = GAMMA1(BASIS)
          OMEGA = GAMMA(BASIS)
          DO 80 J = BASIS - 1,1,-1
              S = ZERO
              DO 70 I = J + 1,BASIS
                  S = S + TAU(J,I)*GAMMA(I)
   70         CONTINUE
              GAMMA(J) = GAMMA1(J) - S
   80     CONTINUE

* 21. gamma''=gamma_{j+1}+\sum_{i=j+1}^{l-1}{tau_{j,i}gamma_{i+1}}
          DO 100 J = 1,BASIS - 1
              S = ZERO
              DO 90 I = J + 1,BASIS - 1
                  S = S + TAU(J,I)*GAMMA(I+1)
   90         CONTINUE
              GAMMA2(J) = GAMMA(J+1) + S
  100     CONTINUE

*  Update

* 22. x(0)=x(0)+gamma(1)r(0)
          CALL DAXPY(LOCLEN,GAMMA(1),WRK(IR),1,X,1)

* 23. r(0)=r(0)-gamma'(l)r(l)
          CALL DAXPY(LOCLEN,-GAMMA1(BASIS),WRK(IR+BASIS*LOCLEN),1,
     +               WRK(IR),1)

* 24. u(0)=u(0)-gamma(l)u(l)
          CALL DAXPY(LOCLEN,-GAMMA(BASIS),WRK(IU+BASIS*LOCLEN),1,
     +               WRK(IU),1)

          I0 = LOCLEN
          DO 110 J = 1,BASIS - 1

* 25. u(0)=u(0)-gamma(j)u(j), j=1:l-1
              CALL DAXPY(LOCLEN,-GAMMA(J),WRK(IU+I0),1,WRK(IU),1)

* 26. x(0)=x(0)+gamma''(j)r(j), j=1:l-1
              CALL DAXPY(LOCLEN,GAMMA2(J),WRK(IR+I0),1,X,1)

* 27. r(0)=r(0)-gamma'(j)r(j), j=1:l-1
              CALL DAXPY(LOCLEN,-GAMMA1(J),WRK(IR+I0),1,WRK(IR),1)
              I0 = I0 + LOCLEN
  110     CONTINUE

  120 CONTINUE

      IF (ITNO.GT.MAXIT) THEN
          STATUS = -1
          ITNO = MAXIT
      END IF

 9999 CONTINUE

      IF ((PRECONT.EQ.2) .OR. (PRECONT.EQ.3)) THEN
          CALL DCOPY(LOCLEN,X,1,WRK(IZ),1)
          CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),X,ipar)
      END IF

*  Set output parameters
      ipar(11) = ITNO
      ipar(12) = STATUS
      ipar(13) = STEPERR
      dpar(2) = EXITNORM

      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c      
      SUBROUTINE pimRGMRES
     &   (neq,nnz,pret,maxit,stopt,eps,nit,ierr,iq,jq,dq,q,ia,ja,a,b,x,
     &    matvec,preconl,preconr,pdsum,pdnrm,progress)
c                                                                      c      
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Solve system of equations using RGMRES iterative method          **
c ** - --                            ------                           **
c **********************************************************************
c Parameters
      INTEGER  neq,nnz,pret,maxit,stopt,nit,ierr,
     &         iq(neq+1),jq(nnz),dq(neq),
     &         ia(neq+1),ja(nnz),da(neq)
      REAL*8   q(nnz),a(nnz),b(neq),x(neq),eps,pdnrm
      EXTERNAL matvec,preconl,preconr,pdsum,pdnrm,progress
c Internal variables
      INTEGER  nwrk,ipar(13),pres,rst
      REAL*8   wrk,dpar(6)
      EXTERNAL dinit,dvprod,drgmres
c
      ALLOCATABLE wrk(:)
c
c Set initial values
c
      rst=10
c if pret>0 use left preconditioning
      pres=0
      IF (pret.GT.0) pres=1
c
c Allocate memory for matrix and vectors
c
      nwrk = (4+rst)*neq
      ALLOCATE(wrk(nwrk))
      CALL pimdsetpar(ipar,dpar,neq,neq,neq,neq,rst,-1,-1,
     &                pres,stopt,maxit,eps)
c
c RGMRES
c
      CALL drgmres(neq,nnz,nwrk,iq,jq,dq,q,ia,ja,a,b,x,wrk,ipar,dpar,
     &             matvec,preconl,preconr,pdsum,pdnrm,progress)
c
c Set output values
c
      nit=ipar(11)
      ierr=ipar(12)
c
c Free memory
c
      DEALLOCATE(wrk)
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE drgmres
     &   (neq,nnz,nwrk,iq,jq,dq,q,ia,ja,a,b,x,wrk,ipar,dpar,
     &    matvec,preconl,preconr,pdsum,pdnrm2,progress)
c                                                                      c
C----------------------------------------------------------------------C
*     ..
*     .. Parameters ..
      REAL*8  ZERO
      PARAMETER (ZERO=0.0D0)
      REAL*8  ONE
      PARAMETER (ONE=1.0D0)
      INTEGER IPARSIZ
      PARAMETER (IPARSIZ=13)
      INTEGER DPARSIZ
      PARAMETER (DPARSIZ=6) ! JURE RAVNIK - Changed "=2" to "=6"
      INTEGER IBDIM
      PARAMETER (IBDIM=50)
      INTEGER LDR
      PARAMETER (LDR=IBDIM+1)
      INTEGER LDG
      PARAMETER (LDG=IBDIM+2)
*     ..
*     .. Matrices A,Q ..
      INTEGER neq,nnz,nwrk,iq(neq+1),jq(nnz),dq(neq),ia(neq+1),ja(nnz)
      REAL*8  q(nnz),a(nnz)
*     ..
*     .. Array Arguments ..
      INTEGER ipar(IPARSIZ)
      REAL*8  b(neq),x(neq),dpar(DPARSIZ),wrk(nwrk)
*     ..
*     .. Function Arguments ..
      REAL*8   PDNRM2
      EXTERNAL PDNRM2
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL MATVEC,PDSUM,PRECONL,PRECONR,PROGRESS
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION BETA,EPSILON,ETA,EXITNORM,KSI,RHSSTOP,TAU1,TAU2
      INTEGER BASISDIM,BLKSZ,I,IRES,ITNO,IV,IW,IZ,J,K0,K1,LDA,LOCLEN,
     +        MAXIT,N,NPROCS,PRECONTYPE,PROCID,STATUS,STEPERR,STOPTYPE
      LOGICAL ENDED
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION G(LDG),R(LDR,LDR),RHO(IBDIM)
*     ..
*     .. External Functions ..
      DOUBLE PRECISION DDOT,DSETRHSSTOP
      EXTERNAL DDOT,DSETRHSSTOP
*     ..
*     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,DECODE,DINIT,DSCAL,DTRSV,ENCODE,GIVENS,
     +         PIMDGETPAR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS
*     ..
      CALL PIMDGETPAR(IPAR,DPAR,LDA,N,BLKSZ,LOCLEN,BASISDIM,NPROCS,
     +                PROCID,PRECONTYPE,STOPTYPE,MAXIT,ITNO,STATUS,
     +                STEPERR,EPSILON,EXITNORM)

*  Check consistency of preconditioning and stop types
      IF (((PRECONTYPE.EQ.0).OR. (PRECONTYPE.EQ.2)) .AND.
     +    (STOPTYPE.EQ.6)) THEN
          ITNO = 0
          STATUS = -4
          STEPERR = 0
          GO TO 9999

      END IF

*  Set indices for mapping local vectors into wrk
      IRES = 1
      IZ = IRES + LOCLEN
      IW = IZ + LOCLEN
      IV = IW + LOCLEN

*  Set rhs of stopping criteria
      RHSSTOP = DSETRHSSTOP(NEQ,NNZ,IQ,JQ,DQ,Q,B,WRK(IRES),EPSILON,IPAR,
     +                      PRECONL,PDNRM2)

*  1. r=Q1(b-AQ2x)
      IF (PRECONTYPE.EQ.0) THEN
*     r=b-Ax
          CALL DCOPY(LOCLEN,B,1,WRK(IRES),1)
          CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK(IW),IPAR)
          CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IRES),1)

      ELSE IF (PRECONTYPE.EQ.1) THEN
*     r=Q1(b-Ax)
          CALL DCOPY(LOCLEN,B,1,WRK(IZ),1)
          CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK(IW),IPAR)
          CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IZ),1)
          CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IRES),IPAR)

      ELSE IF (PRECONTYPE.EQ.2) THEN
*     r=b-AQ2x
          CALL DCOPY(LOCLEN,B,1,WRK(IRES),1)
          CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK(IW),IPAR)
          CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IZ),IPAR)
          CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IRES),1)

      ELSE IF (PRECONTYPE.EQ.3) THEN
*     r=Q1(b-AQ2x)
          CALL DCOPY(LOCLEN,B,1,WRK(IW),1)
          CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK(IRES),IPAR)
          CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IRES),WRK(IZ),IPAR)
          CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IW),1)
          CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IRES),IPAR)
      END IF

*  2. beta=||r||_2
      BETA = PDNRM2(LOCLEN,WRK(IRES),IPAR)

*  Loop
      STATUS = 0
      STEPERR = -1
      EXITNORM = -ONE
      ENDED = .FALSE.
      DO 20 ITNO = 1,MAXIT

*  3. g=(beta,beta,...)
          G(1) = BETA
          G(2) = BETA

*  4. V(1)=r/beta
          IF (BETA.EQ.ZERO) THEN
              STATUS = -3
              STEPERR = 4
              GO TO 9999

          END IF

          CALL DCOPY(LOCLEN,WRK(IRES),1,WRK(IV),1)
          CALL DSCAL(LOCLEN,ONE/BETA,WRK(IV),1)

          K0 = 0
          DO 40 J = 1,BASISDIM

*     z=Q1AQ2V(j)
              IF (PRECONTYPE.EQ.0) THEN
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IV+K0),WRK(IZ),IPAR)

              ELSE IF (PRECONTYPE.EQ.1) THEN
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IV+K0),WRK(IW),IPAR)
                  CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IZ),IPAR)

              ELSE IF (PRECONTYPE.EQ.2) THEN
                  CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IV+K0),WRK(IW),IPAR)
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IZ),IPAR)

              ELSE IF (PRECONTYPE.EQ.3) THEN
                  CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IV+K0),WRK(IZ),IPAR)
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IZ),WRK(IW),IPAR)
                  CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IZ),IPAR)
              END IF


*  5. R(i,j)=dot(V(i),Q1AQ2V(j))
              K1 = 0
              DO 50 I = 1,J
                  R(I,J) = DDOT(LOCLEN,WRK(IV+K1),1,WRK(IZ),1)
                  K1 = K1 + LOCLEN
   50         CONTINUE
              CALL PDSUM(J,R(1,J),IPAR)

*  6. Vhat(j)=Q1AQ2V(j)-sum_{i=1}^{j}{R(i,j)V(i)}
              K1 = 0
              CALL DINIT(LOCLEN,ZERO,WRK(IW),1)
              DO 60 I = 1,J
                  CALL DAXPY(LOCLEN,R(I,J),WRK(IV+K1),1,WRK(IW),1)
                  K1 = K1 + LOCLEN
   60         CONTINUE
              CALL DSCAL(LOCLEN,-ONE,WRK(IW),1)
              CALL DAXPY(LOCLEN,ONE,WRK(IZ),1,WRK(IW),1)

*  From this point, w holds the (j+1)-st column of vhat

*  7. R(j+1,j)=||Vhat(j)||_2
              R(J+1,J) = PDNRM2(LOCLEN,WRK(IW),IPAR)

*  8. V(j+1)=Vhat(j)/R(j+1,j)
              IF (R(J+1,J).EQ.ZERO) THEN
                  STATUS = -2
                  STEPERR = 8
                  GO TO 9999

              END IF

              K0 = K0 + LOCLEN
              CALL DSCAL(LOCLEN,ONE/R(J+1,J),WRK(IW),1)
              CALL DCOPY(LOCLEN,WRK(IW),1,WRK(IV+K0),1)

*  9. Apply previous Givens' rotations to column j of R
              DO 70 I = 1,J - 1
                  CALL DECODE(RHO(I),KSI,ETA)
                  TAU1 = R(I,J)
                  TAU2 = R(I+1,J)
                  R(I,J) = KSI*TAU1 - ETA*TAU2
                  R(I+1,J) = ETA*TAU1 + KSI*TAU2
   70         CONTINUE

* 10. Compute Givens' rotation to zero element R(j+1,j)
              CALL GIVENS(R(J,J),R(J+1,J),KSI,ETA)
              TAU1 = R(J,J)
              TAU2 = R(J+1,J)
              R(J,J) = KSI*TAU1 - ETA*TAU2
              R(J+1,J) = ETA*TAU1 + KSI*TAU2
              CALL ENCODE(RHO(J),KSI,ETA)

*  11. Update g
              G(J) = G(J)*KSI
              G(J+1) = G(J+1)*ETA
              G(J+2) = G(J+1)

*  12. If |g(j+1)|<rhsstop stop
              EXITNORM = ABS(G(J+1))
              IF (EXITNORM.LT.RHSSTOP) THEN
                  BASISDIM = J
                  ENDED = .TRUE.
                  GO TO 80

              END IF

   40     CONTINUE
   80     CONTINUE

*  13. Solve Ry=g (solution to least-squares problem)
          CALL DTRSV('U','N','N',BASISDIM,R,LDR,G,1)

*  14. x=x+Vy (Form approximate solution after a c-cycle)
          K1 = 0
          DO 100 I = 1,BASISDIM
              CALL DAXPY(LOCLEN,G(I),WRK(IV+K1),1,X,1)
              K1 = K1 + LOCLEN
  100     CONTINUE

*  Call monitoring routine
          CALL PROGRESS(LOCLEN,ITNO,EXITNORM,X,WRK(IRES),WRK(IRES))

          IF (ENDED) GO TO 9999

*  15. r=Q1(b-AQ2x)
          IF (PRECONTYPE.EQ.0) THEN
*     r=b-Ax
              CALL DCOPY(LOCLEN,B,1,WRK(IRES),1)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK(IW),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IRES),1)

          ELSE IF (PRECONTYPE.EQ.1) THEN
*     r=Q1(b-Ax)
              CALL DCOPY(LOCLEN,B,1,WRK(IZ),1)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK(IW),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IZ),1)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IRES),IPAR)

          ELSE IF (PRECONTYPE.EQ.2) THEN
*     r=b-AQ2x
              CALL DCOPY(LOCLEN,B,1,WRK(IRES),1)
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK(IW),IPAR)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IZ),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IRES),1)

          ELSE IF (PRECONTYPE.EQ.3) THEN
*     r=Q1(b-AQ2x)
              CALL DCOPY(LOCLEN,B,1,WRK(IW),1)
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK(IRES),IPAR)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IRES),WRK(IZ),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IW),1)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IRES),IPAR)
          END IF

*  16. beta=||r||_2
          BETA = PDNRM2(LOCLEN,WRK(IRES),IPAR)

   20 CONTINUE

      IF (ITNO.GT.MAXIT) THEN
          STATUS = -1
          ITNO = MAXIT
      END IF

 9999 CONTINUE

      IF ((PRECONTYPE.EQ.2) .OR. (PRECONTYPE.EQ.3)) THEN
          CALL DCOPY(LOCLEN,X,1,WRK(IZ),1)
          CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),X,IPAR)
      END IF

*  Set output parameters
      IPAR(11) = ITNO
      IPAR(12) = STATUS
      IPAR(13) = STEPERR
      DPAR(2) = EXITNORM

      RETURN

      END

      SUBROUTINE ENCODE(RHO,C,S)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      DOUBLE PRECISION C,RHO,S
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS,SIGN
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D0)
*     ..
      IF (C.EQ.ZERO) THEN
          RHO = ONE

      ELSE IF (ABS(S).LT.ABS(C)) THEN
          RHO = SIGN(ONE,C)*S/TWO

      ELSE
          RHO = TWO*SIGN(ONE,S)/C
      END IF

      RETURN

      END

      SUBROUTINE DECODE(RHO,C,S)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      DOUBLE PRECISION C,RHO,S
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D0)
*     ..
      IF (RHO.EQ.ONE) THEN
          C = ZERO
          S = ONE

      ELSE IF (ABS(RHO).LT.ONE) THEN
          S = TWO*RHO
          C = SQRT(ONE-S**2)

      ELSE
          C = TWO/RHO
          S = SQRT(ONE-C**2)
      END IF

      RETURN

      END

      SUBROUTINE GIVENS(A,B,C,S)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,C,S
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TAU
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
*     ..
      IF (B.EQ.ZERO) THEN
          C = ONE
          S = ZERO

      ELSE IF (ABS(B).GT.ABS(A)) THEN
          TAU = -A/B
          S = ONE/SQRT(ONE+TAU**2)
          C = S*TAU

      ELSE
          TAU = -B/A
          C = ONE/SQRT(ONE+TAU**2)
          S = C*TAU
      END IF

      RETURN

      END
C^L
C----------------------------------------------------------------------C
c----------------------------------------------------------------------c
c                                                                      c
c     PRECONDITIONER PROGRAM BLOCK                                     c
c                                                                      c
c----------------------------------------------------------------------c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  SUBROUTINES USED :                                              **
c **                                                                  **
c **                   - frmdia                                       **
c **                   - frmdiafm                                     **
c **                   - frmilu                                       **
c **                                                                  **
c **********************************************************************
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE frmdia(neq,nnz,pres,iq,jq,dq,q)

      INTEGER neq,nnz,pres,iq(neq+1),jq(nnz),dq(neq)
      REAL*8  q(nnz)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Form Preconditioner with Diagonal matrix                        **
c **  - --                     ---                                    **
c **********************************************************************
      INTEGER i,j
c
      IF (pres.EQ.1.OR.pres.EQ.2) THEN
        DO i = 1,neq
	    j=dq(i)
          q(j) = 1.D0/q(j)
        END DO
      ELSE IF (pres.EQ.3) THEN
        DO i = 1,neq
	    j=dq(i)
          q(j) = 1.D0/SQRT(q(j))
        END DO
      END IF
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE frmdiafm(neq,pres,a,q)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Form Preconditioner with Diagonal matrix (full matrix)          **
c **  - --                     ---              -    -                **
c **                                                                  **
c **  neq  ... number of equations                                    **
c **  pres ... preconditioner side (1=left,2=right,3=both)            **
c **  a    ... matrix A                                               **
c **  q  ..... preconditioner vector Q                                **
c **                                                                  **
c **********************************************************************
c Arguments
      INTEGER neq,pres
      REAL*8  a(neq,neq),q(neq)
c Internal
      INTEGER i
c
      IF (pres.EQ.1.OR.pres.EQ.2) THEN
        DO i = 1,neq
          q(i) = 1.D0/a(i,i)
        END DO
      ELSE IF (pres.EQ.3) THEN
        DO i = 1,neq
           q(i) = 1.D0/SQRT(a(i,i))
        END DO
      END IF
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE frmilu(neq,nnz,iq,jq,dq,q)

      INTEGER neq,nnz,iq(neq+1),jq(nnz),dq(neq)
      REAL*8  q(nnz)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Form Preconditioner with Incomplete LU decomposition            **
c **  - --                     -          --                          **
c **********************************************************************
      INTEGER jj,i,j,k,l,m,jbeg,jend
      ALLOCATABLE jj(:)
c
c allocate memory for index vector
c
      ALLOCATE(jj(neq))
c
c loop over the rows
c
      DO i=2,neq
        jbeg=iq(i)+1
        jend=iq(i+1)-1
c
c build index vector for current row
c
        DO j=jbeg,jend
          jj(jq(j))=j
        END DO
c
c loop over the rows abow the current row
c
        DO j=iq(i),dq(i)-1
          m=jq(j)
c
c RHS factor, q_ik = q_ik/q_kk
c
          q(j)=q(j)/q(dq(m))
c
c elimination, q_ij = q_ij - (q_ik/q_kk)*q_kj
c
          DO l=dq(m)+1,iq(m+1)-1
            k=jj(jq(l))
            IF (k.EQ.0) GOTO 130
            q(k)=q(k)-q(j)*q(l)
130       END DO
        END DO
c
c clear index vector
c
        DO j=jbeg,jend
          jj(jq(j))=0
        END DO
      END DO
c
c free memory
c
      DEALLOCATE(jj)
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c      
      SUBROUTINE slvlsqr
     &   (m,n,nnz,pret,maxit,stopt,eps4,nit,ierr,
     &    iq,jq,dq,q,ia,ja,a,b,x)
c                                                                      c      
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Solve system of equations using Least Square iterative method    **
c ** - --                            -     --  -                      **
c **********************************************************************
c Parameters
      INTEGER  m,n,nnz,pret,maxit,stopt,nit,ierr,
     &         iq(m+1),jq(nnz),dq(m),
     &         ia(m+1),ja(nnz),da(m)
      REAL*8   eps,q(nnz),a(nnz),b(m),x(n)
c Internal variables
      LOGICAL  wantse,extra,hotstart,precond
      INTEGER  istop,iout
      REAL*8   x0,wm,wn,u,se,v,w
      REAL*8   zero,damp,conlim,anorm,acond,rnorm,arnorm,xnorm
      REAL*4   eps4
c
      ALLOCATABLE x0(:),wm(:),wn(:),u(:),se(:),v(:),w(:)
c
c Allocate memory for matrix and vectors
c
      ALLOCATE(x0(n),wm(m),wn(n),u(m),se(1),v(n),w(n))
c
c Set initial values
c
      ierr=0
      eps=DBLE(eps4)
      wantse=.FALSE.
      extra=.FALSE.
      hotstart=.TRUE.
      precond=.TRUE.
c
      zero=0.0D0
      damp=zero
      conlim=1.0D8
c
      IF (pret.EQ.0) THEN
        precond=.FALSE.
      END IF
c
      IF (maxit.EQ.0) THEN
        maxit = 10*(m + n + 50)
        damp=zero
      END IF
c
      IF (extra) THEN
        iout=66
      ELSE
        iout=-1
      END IF
c
      IF (hotstart) THEN
        CALL dcopy(n,x,1,x0,1)                       !x0=x
        CALL dinit(m,zero,wm,1)                      !
        CALL aprod(1,m,n,nnz,x0,wm,iA,jA,A)          !wm=Ax
        CALL dcopy(m,b,1,u,1)                        !
        CALL daxpy(m,-1.0D0,wm,1,u,1)                !u=b-wm
      ELSE
        CALL dcopy(m,b,1,u,1)                        !u=b
      END IF
c
      CALL dcopy(nnz,A,1,Q,1)                        !Q=A
      IF (precond) THEN
        CALL aprod(3,m,n,nnz,wn,wm,iA,jA,A)          !M^-1=1/diag(sqrt(AT*A))
        CALL aprod(4,m,n,nnz,wn,wm,iQ,jQ,Q)          !Q=Q*M^-1
      END IF
c
c LSQR
c
      CALL lsqr(
     &  m,n,nnz,damp,wantse,
     &  iQ,jQ,Q,
     &  u,x,se,v,w,
     &  eps,eps,conlim,maxit,
     &  istop,nit,anorm,acond,rnorm,arnorm,xnorm)
c
      IF (precond) THEN
        CALL dvprod(n,wn,1,x,1)                      !x=M^-1*x
      END IF
c
      IF (hotstart) THEN
        CALL daxpy(n,1.0D0,x0,1,x,1)                 !x=x0+dx
      END IF
c
      IF (istop.GT.3) THEN
        WRITE (*,*) 'Warning: LSQR: stoping condition ISTOP= ',istop
        ierr=-1
      END IF
c      print *,istop,nit
c      print *,anorm,acond,rnorm,arnorm,xnorm
c
c Free memory
c
      DEALLOCATE(x0,wm,wn,u,se,v,w)
c
      RETURN
      END
c
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE LSQR  ( 
     &  m, n, nnz, damp, wantse,
     &  ia, ja, a,
     &  u, x, se, v, w,
     &  atol, btol, conlim, itnlim, 
     &  istop, itn, anorm, acond, rnorm, arnorm, xnorm)

      LOGICAL            wantse
      INTEGER            m, n, nnz, itnlim, nout, istop, itn
      INTEGER            ia(m+1), ja(nnz)
      DOUBLE PRECISION   a(nnz), u(m), x(n), se(*), v(n), w(n),
     &                   atol, btol, conlim, damp,
     &                   anorm, acond, rnorm, arnorm, xnorm
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c
c     LSQR  finds a solution x to the following problems:
c
c     1. Unsymmetric equations --    solve  A*x = b
c
c     2. Linear least squares  --    solve  A*x = b
c                                    in the least-squares sense
c
c     3. Damped least squares  --    solve  (   A    )*x = ( b )
c                                           ( damp*I )     ( 0 )
c                                    in the least-squares sense
c
c     where A is a matrix with m rows and n columns, b is an
c     m-vector, and damp is a scalar.  (All quantities are REAL.)
c     The matrix A is intended to be large and sparse.  It is accessed
c     by means of SUBROUTINE CALLs of the form
c
c                CALL aprod ( mode, m, n, nnz, x, y, ia, ja, a )
c
c     which must perform the following functions:
c
c                If mode = 1, compute  y = y + A*x.
c                If mode = 2, compute  x = x + A(transpose)*y.
c
c     The vectors x and y are input parameters in both cases.
c     If  mode = 1,  y should be altered without changing x.
c     If  mode = 2,  x should be altered without changing y.
c     The parameters, ia, ja, a  may be used for workspace
c     as described below.
c
c     The rhs vector b is input via u, and subsequently overwritten.
c
c
c     Note:  LSQR uses an iterative method to approximate the solution.
c     The number of iterations required to reach a certain accuracy
c     depends strongly on the scaling of the problem.  Poor scaling of
c     the rows or columns of A should therefore be avoided where
c     possible.
c
c     For example, in problem 1 the solution is unaltered by
c     row-scaling.  If a row of A is very small or large compared to
c     the other rows of A, the corresponding row of ( A  b ) should be
c     scaled up or down.
c
c     In problems 1 and 2, the solution x is easily recovered
c     following column-scaling.  Unless better information is known,
c     the nonzero columns of A should be scaled so that they all have
c     the same Euclidean norm (e.g., 1.0).
c
c     In problem 3, there is no freedom to re-scale if damp is
c     nonzero.  However, the value of damp should be assigned only
c     after attention has been paid to the scaling of A.
c
c     The parameter damp is intended to help regularize
c     ill-conditioned systems, by preventing the true solution from
c     being very large.  Another aid to regularization is provided by
c     the parameter acond, which may be used to terminate iterations
c     before the computed solution becomes very large.
c
c     Note that x is not an input parameter.
c     If some initial estimate x0 is known and if damp = 0,
c     one could proceed as follows:
c
c       1. Compute a residual vector     r0 = b - A*x0.
c       2. Use LSQR to solve the system  A*dx = r0.
c       3. Add the correction dx to obtain a final solution x = x0 + dx.
c
c     This requires that x0 be available before and after the CALL
c     to LSQR.  To judge the benefits, suppose LSQR takes k1 iterations
c     to solve A*x = b and k2 iterations to solve A*dx = r0.
c     If x0 is "good", norm(r0) will be smaller than norm(b).
c     If the same stopping tolerances atol and btol are used for each
c     system, k1 and k2 will be similar, but the final solution x0 + dx
c     should be more accurate.  The only way to reduce the total work
c     is to use a larger stopping tolerance for the second system.
c     If some value btol is suitable for A*x = b, the larger value
c     btol*norm(b)/norm(r0)  should be suitable for A*dx = r0.
c
c     Preconditioning is another way to reduce the number of iterations.
c     If it is possible to solve a related system M*x = b efficiently,
c     where M approximates A in some helpful way
c     (e.g. M - A has low rank or its elements are small relative to
c     those of A), LSQR may converge more rapidly on the system
c           A*M(inverse)*z = b,
c     after which x can be recovered by solving M*x = z.
c
c     NOTE: If A is symmetric, LSQR should not be used!
c     Alternatives are the symmetric conjugate-gradient method (cg)
c     and/or SYMMLQ.
c     SYMMLQ is an implementation of symmetric cg that applies to
c     any symmetric A and will converge more rapidly than LSQR.
c     If A is positive definite, there are other implementations of
c     symmetric cg that require slightly less work per iteration
c     than SYMMLQ (but will take the same number of iterations).
c
c
c     Notation
c     --------
c
c     The following quantities are used in discussing the SUBROUTINE
c     parameters:
c
c     Abar   =  (   A    ),          bbar  =  ( b )
c               ( damp*I )                    ( 0 )
c
c     r      =  b  -  A*x,           rbar  =  bbar  -  Abar*x
c
c     rnorm  =  sqrt( norm(r)**2  +  damp**2 * norm(x)**2 )
c            =  norm( rbar )
c
c     relpr  =  the relative precision of floating-point arithmetic
c               on the machine being used.  On most machines,
c               relpr is about 1.0e-7 and 1.0D-16 in single and double
c               precision respectively.
c
c     LSQR  minimizes the function rnorm with respect to x.
c
c
c     Parameters
c     ----------
c
c     m       input      m, the number of rows in A.
c
c     n       input      n, the number of columns in A.
c
c     damp    input      The damping parameter for problem 3 above.
c                        (damp should be 0.0 for problems 1 and 2.)
c                        If the system A*x = b is incompatible, values
c                        of damp in the range 0 to sqrt(relpr)*norm(A)
c                        will probably have a negligible effect.
c                        Larger values of damp will tend to decrease
c                        the norm of x and reduce the number of
c                        iterations required by LSQR.
c
c                        The work per iteration and the storage needed
c                        by LSQR are the same for all values of damp.
c
c     wantse  input      A logical variable to say if the array se(*)
c                        of standard error estimates should be computed.
c                        If m .GT. n  or  damp .GT. 0,  the system is
c                        overdetermined and the standard errors may be
c                        useful.  (See the first LSQR reference.)
c                        Otherwise (m .LE. n  and  damp = 0) they DO not
c                        mean much.  Some time and storage can be saved
c                        by setting wantse = .FALSE. and using any
c                        convenient array for se(*), which won't be
c                        touched.
c
c     nnz     input      The length of the workspace array rw.
c     ia      workspace  An INTEGER array of length m+1.
c     ja      workspace  An INTEGER array of length nnz.
c     a       workspace  A REAL array of length nnz.
c
c             Note:  LSQR  does not explicitly use the previous four
c             PARAMETERs, but passes them to SUBROUTINE aprod for
c             possible use as workspace.  
c
c     u(m)    input      The rhs vector b.  Beware that u is
c                        over-written by LSQR.
c
c     v(n)    workspace
c
c     w(n)    workspace
c
c     x(n)    output     Returns the computed solution x.
c
c     se(*)   output     If wantse is true, the DIMENSION of se must be
c             (maybe)    n or more.  se(*) THEN RETURNs standard error
c                        estimates for the components of x.
c                        For each i, se(i) is set to the value
c                           rnorm * sqrt( sigma(i,i) / t ),
c                        where sigma(i,i) is an estimate of the i-th
c                        diagonal of the inverse of Abar(transpose)*Abar
c                        and  t = 1      if  m .LE. n,
c                             t = m - n  if  m .GT. n  and  damp = 0,
c                             t = m      if  damp .NE. 0.
c
c                        If wantse is false, se(*) will not be touched.
c                        The actual parameter can be any suitable array
c                        of any length.
c
c     atol    input      An estimate of the relative error in the data
c                        defining the matrix A.  For example,
c                        if A is accurate to about 6 digits, set
c                        atol = 1.0e-6 .
c
c     btol    input      An estimate of the relative error in the data
c                        defining the rhs vector b.  For example,
c                        if b is accurate to about 6 digits, set
c                        btol = 1.0e-6 .
c
c     conlim  input      An upper limit on cond(Abar), the apparent
c                        condition number of the matrix Abar.
c                        Iterations will be terminated if a computed
c                        estimate of cond(Abar) exceeds conlim.
c                        This is intended to prevent certain small or
c                        zero singular values of A or Abar from
c                        coming into effect and causing unwanted growth
c                        in the computed solution.
c
c                        conlim and damp may be used separately or
c                        together to regularize ill-conditioned systems.
c
c                        Normally, conlim should be in the range
c                        1000 to 1/relpr.
c                        Suggested value:
c                        conlim = 1/(100*relpr)  for compatible systems,
c                        conlim = 1/(10*sqrt(relpr)) for least squares.
c
c             Note:  If the user is not concerned about the parameters
c             atol, btol and conlim, any or all of them may be set
c             to zero.  The effect will be the same as the values
c             relpr, relpr and 1/relpr respectively.
c
c     itnlim  input      An upper limit on the number of iterations.
c                        Suggested value:
c                        itnlim = n/2   for well-conditioned systems
c                                       with clustered singular values,
c                        itnlim = 4*n   otherwise.
c
c     nout    input      File number for printed output.  If positive,
c                        a summary will be printed on file nout.
c
c     istop   output     An integer giving the reason for termination:
c
c                0       x = 0  is the exact solution.
c                        No iterations were performed.
c
c                1       The equations A*x = b are probably
c                        compatible.  Norm(A*x - b) is sufficiently
c                        small, given the values of atol and btol.
c
c                2       damp is zero.  The system A*x = b is probably
c                        not compatible.  A least-squares solution has
c                        been obtained that is sufficiently accurate,
c                        given the value of atol.
c
c                3       damp is nonzero.  A damped least-squares
c                        solution has been obtained that is sufficiently
c                        accurate, given the value of atol.
c
c                4       An estimate of cond(Abar) has exceeded
c                        conlim.  The system A*x = b appears to be
c                        ill-conditioned.  Otherwise, there could be an
c                        error in SUBROUTINE aprod.
c
c                5       The iteration limit itnlim was reached.
c
c     itn     output     The number of iterations performed.
c
c     anorm   output     An estimate of the Frobenius norm of  Abar.
c                        This is the square-root of the sum of squares
c                        of the elements of Abar.
c                        If damp is small and if the columns of A
c                        have all been scaled to have length 1.0,
c                        anorm should increase to roughly sqrt(n).
c                        A radiCALLy different value for anorm may
c                        indicate an error in SUBROUTINE aprod (there
c                        may be an inconsistency between modes 1 and 2).
c
c     acond   output     An estimate of cond(Abar), the condition
c                        number of Abar.  A very high value of acond
c                        may again indicate an error in aprod.
c
c     rnorm   output     An estimate of the final value of norm(rbar),
c                        the function being minimized (see notation
c                        above).  This will be small if A*x = b has
c                        a solution.
c
c     arnorm  output     An estimate of the final value of
c                        norm( Abar(transpose)*rbar ), the norm of
c                        the residual for the usual normal equations.
c                        This should be small in all cases.  (arnorm
c                        will often be smaller than the true value
c                        computed from the output vector x.)
c
c     xnorm   output     An estimate of the norm of the final
c                        solution vector x.
c
c
c     Subroutines and functions used
c     ------------------------------
c
c     USER               aprod
c     LSQR               d2norm
c     BLAS               dcopy, dnrm2, dscal (see Lawson et al. below)
c
c
c     Precision
c     ---------
c
c     The number of iterations required by LSQR will usually decrease
c     if the computation is performed in higher precision.
c     At least 15-digit arithmetic should normally be used.
c     To convert LSQR and D2NORM between single and double precision,
c     change
c                        DOUBLE PRECISION
c                        dcopy, dnrm2, dscal
c     to the appropriate FORTRAN and BLAS equivalents.
c     Also change 'd+' or 'e+' in the parameter statement.
c
c
c     References
c     ----------
c
c     C.C. Paige and M.A. Saunders,  LSQR: An algorithm for sparse
c          linear equations and sparse least squares,
c          ACM Transactions on Mathematical Software 8, 1 (March 1982),
c          pp. 43-71.
c
c     C.C. Paige and M.A. Saunders,  Algorithm 583, LSQR: Sparse
c          linear equations and least-squares problems,
c          ACM Transactions on Mathematical Software 8, 2 (June 1982),
c          pp. 195-209.
c
c     C.L. Lawson, R.J. Hanson, D.R. Kincaid and F.T. Krogh,
c          Basic linear algebra subprograms for Fortran usage,
c          ACM Transactions on Mathematical Software 5, 3 (Sept 1979),
c          pp. 308-323 and 324-325.
c     ------------------------------------------------------------------
c
c
c     LSQR development:
c     22 Feb 1982: LSQR sent to ACM TOMS to become Algorithm 583.
c     15 Sep 1985: Final F66 version.  LSQR sent to "misc" in netlib.
c     13 Oct 1987: Bug (Robert Davies, DSIR).  Have to delete
c                     IF ( (one + DABS(t)) .LE. one ) GO TO 200
c                  from loop 200.  The test was an attempt to reduce
c                  underflows, but caused w(i) not to be updated.
c     17 Mar 1989: First F77 version.
c     04 May 1989: Bug (David Gay, AT&T).  When the second beta is zero,
c                  rnorm = 0 and
c                  test2 = arnorm / (anorm * rnorm) overflows.
c                  Fixed by testing for rnorm = 0.
c     05 May 1989: Sent to "misc" in netlib.
c     14 Mar 1990: Bug (John Tomlin via IBM OSL testing).
c                  Setting rhbar2 = rhobar**2 + dampsq can give zero
c                  if rhobar underflows and damp = 0.
c                  Fixed by testing for damp = 0 specially.
c     15 Mar 1990: Converted to lower case.
c     21 Mar 1990: d2norm introduced to avoid overflow in numerous
c                  items like  c = sqrt( a**2 + b**2 ).
c     04 Sep 1991: wantse added as an argument to LSQR, to make
c                  standard errors optional.  This saves storage and
c                  time when se(*) is not wanted.
c     13 Feb 1992: istop now RETURNs a value in [1,5], not [1,7].
c                  1, 2 or 3 means that x solves one of the problems
c                  Ax = b,  min norm(Ax - b)  or  damped least squares.
c                  4 means the limit on cond(A) was reached.
c                  5 means the limit on iterations was reached.
c     07 Dec 1994: Keep track of dxmax = max_k norm( phi_k * d_k ).
c                  So far, this is just printed at the end.
c                  A large value (relative to norm(x)) indicates
c                  significant cancellation in forming
c                  x  =  D*f  =  sum( phi_k * d_k ).
c                  A large column of D need NOT be serious if the
c                  corresponding phi_k is small.
c     27 Dec 1994: Include estimate of alfa_opt in iteration log.
c                  alfa_opt is the optimal scale factor for the
c                  residual in the "augmented system", as described by
c                  A. Bjorck (1992),
c                  Pivoting and stability in the augmented system method,
c                  in D. F. Griffiths and G. A. Watson (eds.),
c                  "Numerical Analysis 1991",
c                  Proceedings of the 14th Dundee Conference,
c                  Pitman Research Notes in Mathematics 260,
c                  Longman Scientific and Technical, Harlow, Essex, 1992.
c
c
c     Michael A. Saunders                  mike@sol-michael.stanford.edu
c     Dept of Operations Research          na.Msaunders@na-net.ornl.gov
c     Stanford University
c     Stanford, CA 94305-4022              (415) 723-1875
c **                                                                  **
c **********************************************************************

c     EXTERNAL           d2norm, dnrm2, dcopy, dscal
      DOUBLE PRECISION   d2norm, dnrm2

c     Local variables

      LOGICAL            damped, extra
      INTEGER            i, maxdx, nconv, nstop
      DOUBLE PRECISION   alpha, beta, bnorm,
     &                   cs, cs1, cs2, ctol,
     &                   delta, dknorm, dnorm, dxk, dxmax,
     &                   gamma, gambar, phi, phibar, psi,
     &                   res2, rho, rhobar, rhbar1,
     &                   rhs, rtol, sn, sn1, sn2,
     &                   t, tau, temp, test1, test2, test3,
     &                   theta, t1, t2, t3, xnorm1, z, zbar, alfopt

      DOUBLE PRECISION   zero,           one
      PARAMETER        ( zero = 0.0d+0,  one = 1.0d+0 )

      CHARACTER*14       enter, exit
      CHARACTER*53       msg(0:5)

      DATA               enter /' Enter LSQR.  '/
      DATA               exit  /' Exit  LSQR.  '/
      DATA               msg
     &/ 'The exact solution is  x = 0',
     &'A solution to Ax = b was found, given atol, btol',
     &'A least-squares solution was found, given atol',
     &'A damped least-squares solution was found, given atol',
     &'Cond(Abar) seems to be too large, given conlim',
     &'The iteration limit was reached' /
c-----------------------------------------------------------------------

c     Initialize.
      nout=-1
      IF (nout .GT. 0) THEN
        WRITE(nout, 1000) enter, m, n, damp, wantse,
     &  atol, conlim, btol, itnlim
      END IF

      damped =   damp .GT. zero
      extra  =  .true.     ! true for extra printing below.
      itn    =   0
      istop  =   0
      nstop  =   0
      maxdx  =   0
      ctol   =   zero
      IF (conlim .GT. zero) ctol = one / conlim
      anorm  =   zero
      acond  =   zero
      dnorm  =   zero
      dxmax  =   zero
      res2   =   zero
      psi    =   zero
      xnorm  =   zero
      xnorm1 =   zero
      cs2    = - one
      sn2    =   zero
      z      =   zero

c     ------------------------------------------------------------------
c     Set up the first vectors u and v for the bidiagonalization.
c     These satisfy  beta*u = b,  alpha*v = A(transpose)*u.
c     ------------------------------------------------------------------
      DO 10  i = 1, n
        v(i)  =  zero
        x(i)  =  zero
   10 CONTINUE

      IF ( wantse ) THEN
        DO 20  i = 1, n
          se(i) =  zero
   20   CONTINUE
      END IF

      alpha  =   zero
      beta   =   dnrm2 ( m, u, 1 )  ! to paraleliziraj u(zac:kon), racuna normo

      IF (beta .GT. zero) THEN
        CALL dscal ( m, (one / beta), u, 1 ) ! u=u*(one/beta)
        CALL aprod ( 2, m, n, nnz, v, u, ia, ja, a ) ! paraleliziraj
        alpha  =   dnrm2 ( n, v, 1 )
      END IF

      IF (alpha .GT. zero) THEN
        CALL dscal ( n, (one / alpha), v, 1 )
        CALL dcopy ( n, v, 1, w, 1 )
      END IF

      arnorm =   alpha * beta
      IF (arnorm .EQ. zero) GOTO 800

      rhobar =   alpha
      phibar =   beta
      bnorm  =   beta
      rnorm  =   beta

      IF (nout   .GT.  0  ) THEN
        IF ( damped ) THEN
          WRITE(nout, 1300)
        ELSE
          WRITE(nout, 1200)
        END IF
        test1  = one
        test2  = alpha / beta

        IF ( extra ) THEN
          WRITE(nout, 1400)
        END IF
        WRITE(nout, 1500) itn, x(1), rnorm, test1, test2
        WRITE(nout, 1600)
      END IF

c     ==================================================================
c     Main iteration loop.
c     ==================================================================
  100 itn    = itn + 1

c     ------------------------------------------------------------------
c     Perform the next step of the bidiagonalization to obtain the
c     next  beta, u, alpha, v.  These satisfy the relations
c                beta*u  =  A*v  -  alpha*u,
c               alpha*v  =  A(transpose)*u  -  beta*v.
c     ------------------------------------------------------------------
      CALL dscal ( m, (- alpha), u, 1 )
      CALL aprod ( 1, m, n, nnz, v, u, ia, ja, a )
      beta   =   dnrm2 ( m, u, 1 )

c     Accumulate  anorm = || Bk ||
c                       =  sqrt( sum of  alpha**2 + beta**2 + damp**2 ).

      temp   =   d2norm( alpha, beta )
      temp   =   d2norm( temp , damp )
      anorm  =   d2norm( anorm, temp )

      IF (beta .GT. zero) THEN
        CALL dscal ( m, (one / beta), u, 1 )
        CALL dscal ( n, (- beta), v, 1 )
        CALL aprod ( 2, m, n, nnz, v, u, ia, ja, a )
        alpha  =   dnrm2 ( n, v, 1 )
        IF (alpha .GT. zero) THEN
          CALL dscal ( n, (one / alpha), v, 1 )
        END IF
      END IF

c     ------------------------------------------------------------------
c     Use a plane rotation to eliminate the damping parameter.
c     This alters the diagonal (rhobar) of the lower-bidiagonal matrix.
c     ------------------------------------------------------------------
      rhbar1 = rhobar
      IF ( damped ) THEN
        rhbar1 = d2norm( rhobar, damp )
        cs1    = rhobar / rhbar1
        sn1    = damp   / rhbar1
        psi    = sn1 * phibar
        phibar = cs1 * phibar
      END IF

c     ------------------------------------------------------------------
c     Use a plane rotation to eliminate the subdiagonal element (beta)
c     of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
c     ------------------------------------------------------------------
      rho    =   d2norm( rhbar1, beta )
      cs     =   rhbar1 / rho
      sn     =   beta   / rho
      theta  =   sn * alpha
      rhobar = - cs * alpha
      phi    =   cs * phibar
      phibar =   sn * phibar
      tau    =   sn * phi

c     ------------------------------------------------------------------
c     Update  x, w  and (perhaps) the standard error estimates.
c     ------------------------------------------------------------------
      t1     =   phi   / rho
      t2     = - theta / rho
      t3     =   one   / rho
      dknorm =   zero

      IF ( wantse ) THEN
        DO 200  i =  1, n
          t      =  w(i)
          x(i)   =  t1*t  +  x(i)
          w(i)   =  t2*t  +  v(i)
          t      = (t3*t)**2
          se(i)  =  t     +  se(i)
          dknorm =  t     +  dknorm
  200   CONTINUE
      ELSE
        DO 220  i =  1, n
          t      =  w(i)
          x(i)   =  t1*t  +  x(i)
          w(i)   =  t2*t  +  v(i)
          dknorm = (t3*t)**2  +  dknorm
  220   CONTINUE
      END IF

c     ------------------------------------------------------------------
c     Monitor the norm of d_k, the update to x.
c     dknorm = norm( d_k )
c     dnorm  = norm( D_k ),        where   D_k = (d_1, d_2, ..., d_k )
c     dxk    = norm( phi_k d_k ),  where new x = x_k + phi_k d_k.
c     ------------------------------------------------------------------
      dknorm = DSQRT( dknorm )
      dnorm  = d2norm( dnorm, dknorm )
      dxk    = DABS( phi * dknorm )
      IF (dxmax .LT. dxk ) THEN
        dxmax   =  dxk
        maxdx   =  itn
      END IF

c     ------------------------------------------------------------------
c     Use a plane rotation on the right to eliminate the
c     super-diagonal element (theta) of the upper-bidiagonal matrix.
c     Then use the result to estimate  norm(x).
c     ------------------------------------------------------------------
      delta  =   sn2 * rho
      gambar = - cs2 * rho
      rhs    =   phi    - delta * z
      zbar   =   rhs    / gambar
      xnorm  =   d2norm( xnorm1, zbar  )
      gamma  =   d2norm( gambar, theta )
      cs2    =   gambar / gamma
      sn2    =   theta  / gamma
      z      =   rhs    / gamma
      xnorm1 =   d2norm( xnorm1, z     )

c     ------------------------------------------------------------------
c     Test for convergence.
c     First, estimate the norm and condition of the matrix  Abar,
c     and the norms of  rbar  and  Abar(transpose)*rbar.
c     ------------------------------------------------------------------
      acond  =   anorm * dnorm
      res2   =   d2norm( res2 , psi    )
      rnorm  =   d2norm( res2 , phibar )
      arnorm =   alpha * DABS( tau )

c     Now use these norms to estimate certain other quantities,
c     some of which will be small near a solution.

      alfopt =   DSQRT( rnorm / (dnorm * xnorm) )
      test1  =   rnorm /  bnorm
      test2  =   zero
      IF (rnorm .GT. zero) test2 = arnorm / (anorm * rnorm)
      test3  =   one   /  acond
      t1     =   test1 / (one  +  anorm * xnorm / bnorm)
      rtol   =   btol  +  atol *  anorm * xnorm / bnorm

c     The following tests guard against extremely small values of
c     atol, btol  or  ctol.  (The user may have set any or all of
c     the parameters  atol, btol, conlim  to zero.)
c     The effect is equivalent to the normal tests using
c     atol = relpr,  btol = relpr,  conlim = 1/relpr.

      t3     =   one + test3
      t2     =   one + test2
      t1     =   one + t1
      IF (itn .GE. itnlim) istop = 5
      IF (t3  .LE. one   ) istop = 4
      IF (t2  .LE. one   ) istop = 2
      IF (t1  .LE. one   ) istop = 1

c     Allow for tolerances set by the user.

      IF (test3 .LE. ctol) istop = 4
      IF (test2 .LE. atol) istop = 2
      IF (test1 .LE. rtol) istop = 1

c     ------------------------------------------------------------------
c     See if it is time to print something.
c     ------------------------------------------------------------------
      IF (nout  .LE.  0       ) GOTO 600
      IF (n     .LE. 40       ) GOTO 400
      IF (itn   .LE. 10       ) GOTO 400
      IF (itn   .GE. itnlim-10) GOTO 400
      IF (MOD(itn,10) .EQ. 0  ) GOTO 400
      IF (test3 .LE.  2.0*ctol) GOTO 400
      IF (test2 .LE. 10.0*atol) GOTO 400
      IF (test1 .LE. 10.0*rtol) GOTO 400
      IF (istop .NE.  0       ) GOTO 400
      GOTO 600

c     Print a line for this iteration.
c     "extra" is for experimental purposes.

  400 IF ( extra ) THEN
        WRITE(nout, 1500) itn, x(1), rnorm, test1, test2, anorm, acond
     &  , phi, dknorm, dxk, alfopt
      ELSE
        WRITE(nout, 1500) itn, x(1), rnorm, test1, test2, anorm, acond
      END IF
      IF (MOD(itn,10) .EQ. 0) WRITE(nout, 1600)

c     ------------------------------------------------------------------
c     Stop if appropriate.
c     The convergence criteria are required to be met on  nconv
c     consecutive iterations, where  nconv  is set below.
c     Suggested value:  nconv = 1, 2  or  3.
c     ------------------------------------------------------------------
  600 IF (istop .EQ. 0) THEN
        nstop  = 0
      ELSE
        nconv  = 1
        nstop  = nstop + 1
        IF (nstop .LT. nconv  .AND.  itn .LT. itnlim) istop = 0
      END IF
      IF (istop .EQ. 0) GOTO 100

c     ==================================================================
c     End of iteration loop.
c     ==================================================================

c     Finish off the standard error estimates.

      IF ( wantse ) THEN
        t    =   one
        IF (m .GT. n)  t = m - n
        IF ( damped )  t = m
        t    =   rnorm / DSQRT( t )

        DO 700  i = 1, n
          se(i)  = t * DSQRT( se(i) )
  700   CONTINUE
      END IF

c     Decide if istop = 2 or 3.
c     Print the stopping condition.

  800 IF (damped  .AND.  istop .EQ. 2) istop = 3
      IF (nout .GT. 0) THEN
        WRITE(nout, 2000) EXIT, istop, itn,
     &  EXIT, anorm, acond,
     &  EXIT, bnorm, xnorm,
     &  EXIT, rnorm, arnorm
        WRITE(nout, 2100) EXIT, dxmax, maxdx,
     &  EXIT, dxmax/(xnorm + 1.0D-20)
        WRITE(nout, 3000) EXIT, msg(istop)
      END IF

  900 RETURN

c     ------------------------------------------------------------------
 1000 format(// 1p, a, '     Least-squares solution of  Ax = b'
     &/ ' The matrix  A  has', i7, ' rows   and', i7, ' columns'
     &/ ' damp   =', e22.14, 3x,        'wantse =', l10
     &/ ' atol   =', e10.2, 15x,        'conlim =', e10.2
     &/ ' btol   =', e10.2, 15x,        'itnlim =', i10)
 1200 format(// '   Itn       x(1)           Function',
     &'     Compatible   LS        Norm A    Cond A')
 1300 format(// '   Itn       x(1)           Function',
     &'     Compatible   LS     Norm Abar Cond Abar')
 1400 format(80x, '    phi    dknorm   dxk  alfa_opt')
 1500 format(1p, i6, 2e17.9, 4e10.2, e9.1, 3e8.1)
 1600 format(1x)
 2000 format(/ 1p, a, 5x, 'istop  =', i2,   15x, 'itn    =', i8
     &/     a, 5x, 'anorm  =', e12.5, 5x, 'acond  =', e12.5
     &/     a, 5x, 'bnorm  =', e12.5, 5x, 'xnorm  =', e12.5
     &/     a, 5x, 'rnorm  =', e12.5, 5x, 'arnorm =', e12.5)
 2100 format(  1p, a, 5x, 'max dx =', e8.1 , ' occurred at itn ', i8,
     &/     a, 5x, '       =', e8.1 , '*xnorm' )
 3000 format( a, 5x, a )

c     End of LSQR
      END
c
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE aprod(MODE,M,N,NNZ,X,Y,IA,JA,A )

      INTEGER mode,m,n,nnz,ia(m+1),ja(nnz)
      REAL*8  x(n),y(m),a(nnz)
c                                                                      c
C----------------------------------------------------------------------C
c                If mode = 1, compute  y = y + A*x.
c                If mode = 2, compute  x = x + AT*y.
c                If mode = 3, compute  x = 1/diag(sqrt((AT*A)).
c                If mode = 4, compute  A_ij = A_ij*x_j.
C----------------------------------------------------------------------C
      IF (mode.EQ.1) THEN
        CALL matvec1(m,n,nnz,ia,ja,a,x,y)
      ELSE IF (mode.EQ.2) THEN
        CALL tmatvec2(m,n,nnz,ia,ja,a,y,x)
      ELSE IF (mode.EQ.3) THEN
        CALL tmatvec3(m,n,nnz,ia,ja,a,x)
      ELSE IF (mode.EQ.4) THEN
        CALL matvec4(m,n,nnz,ia,ja,a,x)
      END IF
c
      RETURN
      END
c
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE matvec1(m,n,nnz,ia,ja,a,u,v)

      INTEGER m,n,nnz,ia(m+1),ja(nnz)
      REAL*8  a(nnz),u(n),v(m)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Matrix vector product: v = v + A*u                              **
c **                                                                  **
c **********************************************************************
      INTEGER i,ii,j
      REAL*8  sum
c
      DO i = 1,m
        sum=0.0D0
        DO j = ia(i),ia(i+1) - 1
          ii=ja(j)
          sum=sum+a(j)*u(ii)
        END DO
        v(i)=v(i)+sum
      END DO
c
      RETURN
      END
c
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE tmatvec2(m,n,nnz,ia,ja,a,u,v)

      INTEGER m,n,nnz,ia(m+1),ja(nnz)
      REAL*8  a(nnz),u(m),v(n)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Transpose matrix vector product: v = v + AT*u                   **
c **                                                                  **
c **********************************************************************
      INTEGER i,ii,j
c
      DO j = 1,m
        DO i = ia(j),ia(j+1) - 1
          ii=ja(i)
          v(ii) = v(ii)+a(i)*u(j)
        END DO
      END DO
c
      RETURN
      END
c
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE tmatvec3(m,n,nnz,ia,ja,a,v)

      INTEGER m,n,nnz,ia(m+1),ja(nnz)
      REAL*8  a(nnz),v(n)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Create diagonal preconditioner: v = 1 / diag(sqrt(AT * A))      **
c **                                                                  **
c **********************************************************************
      INTEGER i,ii,j
c
      CALL dinit(n,0.D0,v,1)
c
      DO j = 1,m
        DO i = ia(j),ia(j+1) - 1
          ii=ja(i)
          v(ii) = v(ii)+a(i)**2
        END DO
      END DO
c
      DO i=1,n
        v(i)=1.0D0/sqrt(v(i))
      END DO
c
      RETURN
      END
c
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE matvec4(m,n,nnz,ia,ja,a,v)

      INTEGER m,n,nnz,ia(m+1),ja(nnz)
      REAL*8  a(nnz),v(n)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Apply diagonal preconditioner: A_ij = A_ij*v_j                  **
c **                                                                  **
c **********************************************************************
      INTEGER i,ij,j
c
      DO i = 1,m
        DO ij = ia(i),ia(i+1) - 1
          j=ja(ij)
          a(ij)=a(ij)*v(j)
        END DO
      END DO
c
      RETURN
      END
c
C----------------------------------------------------------------------C
c                                                                      c
      FUNCTION          d2norm( a, b )
      DOUBLE PRECISION  d2norm, a, b
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c     d2norm  returns  sqrt( a**2 + b**2 )  with precautions
c     to avoid overflow.
c
c     21 Mar 1990: First version.
c **                                                                  **
c **********************************************************************
      INTRINSIC         abs, sqrt
      DOUBLE PRECISION  scale
      DOUBLE PRECISION  zero
      PARAMETER       ( zero = 0.0d+0 )

      scale  = abs( a ) + abs( b )
      IF (scale .EQ. zero) THEN
        d2norm = zero
      ELSE
        d2norm = scale * sqrt( (a/scale)**2   +  (b/scale)**2 )
      END IF
c
      RETURN
      END
C^L           
C----------------------------------------------------------------------C
c----------------------------------------------------------------------c
c                                                                      c
c     DIRECT SOLVER PROGRAM BLOCK                                             c
c                                                                      c
c----------------------------------------------------------------------c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  SUBROUTINES USED :                                              **
c **                                                                  **
c **                   - SolveYSMP                                    **
c **                   - tdrv                                         **
c **                   - trk                                          **
c **                                                                  **
c **********************************************************************
c$OPTION RANGE OFF
C^L           
C----------------------------------------------------------------------C
c
      SUBROUTINE SolveYSMP (neq,nnz,ia,ja,a,b,x,ierr)

      INTEGER neq,nnz,ia(neq+1),ja(nnz),ierr
      REAL*8  a(nnz),b(neq),x(neq)
c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Solve system of equations using Yale's YSMP solver (tdrv.f)      **
c **                                                                  **
c ** neq  = number of equations                                       **
c ** ia   = row index vector of A                                     **
c ** ja   = column index vector of A                                  **
c ** a    = coefficients vector of A                                  **
c ** b    = right hand side vector                                    **
c ** x    = solution vector                                           **
c ** ierr = error status                                              **
c **          0     no errors detected                                **
c **          N+K   null row in  A  --  ROW = K                       **
c **         2N+K   duplicate entry in A  --  ROW = K                 **
c **         8N+K   zero pivot  --  ROW = K                           **
c **        10N+1   insufficient storage in tdrv                      **
c **        12N+K   insufficient storage in trk  --  ROW = K          **
c **                                                                  **
c **********************************************************************
      INTEGER iwrk
      INTEGER nwrk,ewrk
      REAL*8  rwrk
      POINTER (piwrk,iwrk)
      ALLOCATABLE rwrk(:)
c
c Allocate memory for work vector
c
      nwrk = (6*neq+2)+(9*nnz)
      ALLOCATE(rwrk(nwrk))
      piwrk = LOC(rwrk(1))
c
c Solve system Ax=b
c
      CALL tdrv(neq,nnz,ia,ja,a,b,x,nwrk,iwrk,rwrk,ewrk,ierr)
c     CALL tdrv(neq,nnz,ia,ja,a,b,x,nwrk,rwrk,rwrk,ewrk,ierr)
      IF (ierr.NE.0) THEN
        WRITE(*,*) ' Error: solveysmp: ERR =',ierr,' FREE =',ewrk
        STOP
      END IF
c
c Free memory
c
      DEALLOCATE(rwrk)
c
      RETURN
      END
C^L           
C----------------------------------------------------------------------C
c
        SUBROUTINE  TDRV
     *     (N,NNZ,IA,JA,A,B,Z,NSP,ISP,RSP,ESP,FLAG)              
c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Driver for subroutine for solving sparse nonsymmetric systems    **
c ** of linear equations (track nonzeroes dynamically)                **
c **                                                                  **
c **********************************************************************
C                                                                       
        INTEGER N,NNZ,NSP,LRATIO,IJU,IU,IM,JU
        INTEGER IA(N+1), JA(NNZ), ISP(NSP),
     *          U, ROW, TMP, Q, ESP, MXX, FLAG
C       REAL    A(NNZ),  B(N),  Z(N),  RSP(NSP)
        REAL*8  A(NNZ),  B(N),  Z(N),  RSP(NSP)
C                                                                       
C  SET LRATIO EQUAL TO THE RATIO BETWEEN THE LENGTH OF FLOATING POINT   
C  AND INTEGER ARRAY DATA;  E. G., LRATIO = 1 FOR (REAL, INTEGER),      
C  LRATIO = 2 FOR (DOUBLE PRECISION, INTEGER)                           
C                                                                       
        DATA LRATIO/2/                                                  
C                                                                       
C  ******  INITIALIZE AND DIVIDE UP TEMPORARY STORAGE  *****************
        IJU = 1                                                         
        IU  = IJU +  N                                                  
        Q   = IU  + N+1                                                 
        IM  = Q   + N+1                                                 
        U   = (IM + N - 2 + LRATIO) / LRATIO   +   1                    
        JU  = LRATIO * (U - 1) + 1                                      
        ROW = NSP + 1 -  N                                              
        TMP = ROW -  N                                                  
        MXX = TMP - U                                                   
        ESP = MXX                                                       
C                                                                       
C  ******  CALL ZERO-TRACKING SUBROUTINE  ******************************
        IF (MXX.LT.0)  GO TO 110                                        
        CALL  TRK                                                       
     *     (N, NNZ, IA, JA, A, Z, B,                              
     *      ISP(IJU), ISP(JU), ISP(IU), RSP(U), MXX,                    
     *      ISP(Q), ISP(IM), RSP(ROW), RSP(TMP), FLAG, ESP, LRATIO)  
        IF (FLAG.NE.0)  GO TO 100                                       
        RETURN                                                          
C                                                                       
C ** ERROR:  ERROR DETECTED IN TRK                                      
 100    RETURN                                                          
C ** ERROR:  INSUFFICIENT STORAGE                                       
 110    FLAG = 10*N + 1                                                 
        RETURN                                                          
        END                                                             
C^L           
C----------------------------------------------------------------------C
c
        SUBROUTINE  TRK
     *     (N,NNZ,IA,JA,A,Z,B,IJU,JU,IU,U,MAX,                    
     *      Q,IM,ROW,TMP,FLAG,ESP,LRATIO)                         
c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Subroutine for solving sparse nonsymmetric systems               **
c ** of linear equations (track nonzeroes dynamically)                **
c **                                                                  **
c **********************************************************************
C                                                                       
        INTEGER N,NNZ,MAX,LRATIO,JUMIN,JUMAX,K,LUK,JMIN,JMAX,M,LMAX,
     &          I,J,LUI,JTMP,MU
        INTEGER IA(N+1), JA(NNZ), IJU(N), JU(MAX), IU(N+1), Q(N+1), 
     *          IM(N), FLAG, ESP, VJ, QM                           
C       REAL    A(NNZ), Z(N), B(N), U(MAX), ROW(N), TMP(N), SUM, AKI, DK    
        REAL*8  A(NNZ), Z(N), B(N), U(MAX), ROW(N), TMP(N), SUM, AKI, DK  
C                                                                       
C  ******  INITIALIZE  *************************************************
        JUMIN = 1                                                       
        JUMAX = 0                                                       
        IU(1) = MAX                                                     
C                                                                       
C  ******  FOR EACH ROW  ***********************************************
        DO 20 K=1,N                                                     
C  ******  INITIALIZE Q AND ROW TO THE KTH ROW OF REORDERED A  *********
          LUK = 0                                                       
          Q(N+1) = N+1                                                  
          JMIN = IA(K)                                               
          JMAX = IA(K+1) - 1                                         
          IF (JMIN.GT.JMAX)  GO TO 101                                  
          DO 2 J=JMIN,JMAX                                              
            VJ = JA(J)
            QM = N+1                                                    
   1        M = QM                                                      
            QM = Q(M)                                                   
            IF (QM.LT.VJ)  GO TO 1                                      
            IF (QM.EQ.VJ)  GO TO 102                                    
              LUK = LUK+1                                               
              Q(M) = VJ                                                 
              Q(VJ) = QM                                                
              ROW(VJ) = A(J)                                            
   2        CONTINUE                                                    
C                                                                       
C  ******  LINK THROUGH Q  *********************************************
          LMAX = 0                                                      
          IJU(K) = JUMAX                                                
          I = N+1                                                       
   3      I = Q(I)                                                      
          LUK = LUK-1                                                   
          IF (I.GE.K)  GO TO 8                                          
            QM = I                                                      
            LUI = 0                                                     
            JMIN = IJU(I)                                               
            JMAX = IM(I)                                                
            IF (JMIN.GT.JMAX)  GO TO 7                                  
C  ******  AND FIND NONZERO STRUCTURE OF KTH ROW OF L AND U  ***********
              DO 5 J=JMIN,JMAX                                          
                VJ = JU(J)                                              
                IF (VJ.GT.K)  LUI = LUI+1                               
   4            M = QM                                                  
                QM = Q(M)                                               
                IF (QM.LT.VJ)  GO TO 4                                  
                IF (QM.EQ.VJ)  GO TO 5                                  
                  LUK = LUK+1                                           
                  Q(M) = VJ                                             
                  Q(VJ) = QM                                            
                  ROW(VJ) = 0                                           
                  QM = VJ                                               
   5            CONTINUE                                                
C  ******  ADJUST IJU AND IM  ******************************************
              JTMP = JMAX - LUI                                         
              IF (LUI.LE.LMAX)  GO TO 6                                 
                LMAX = LUI                                              
                IJU(K) = JTMP+1                                         
   6          IF (JTMP.LT.JMIN)  GO TO 7                                
                IF (JU(JTMP).EQ.K)  IM(I) = JTMP                        
   7        GO TO 3                                                     
C                                                                       
C  ******  SEE IF JU STORAGE CAN BE COMPRESSED  ************************
   8      IF (I.NE.K)  GO TO 105                                        
          IF (LUK.EQ.LMAX)  GO TO 14                                    
            I = Q(K)                                                    
            IF (JUMIN.GT.JUMAX)  GO TO 12                               
              DO 9 JMIN=JUMIN,JUMAX                                     
                IF (JU(JMIN)-I)  9, 10, 12                              
   9            CONTINUE                                                
              GO TO 12                                                  
  10          IJU(K) = JMIN                                             
              DO 11 J=JMIN,JUMAX                                        
                IF (JU(J).NE.I)  GO TO 12                               
                I = Q(I)                                                
                IF (I.GT.N)  GO TO 14                                   
  11            CONTINUE                                                
              JUMAX = JMIN - 1                                          
C  ******  STORE POINTERS IN JU  ***************************************
  12        I = K                                                       
            JUMIN = JUMAX +  1                                          
            JUMAX = JUMAX + LUK                                         
            IF (JUMAX.GT.LRATIO*(IU(K)-1)+1)  GO TO 112                 
            DO 13 J=JUMIN,JUMAX                                         
              I = Q(I)                                                  
  13          JU(J) = I                                                 
            IJU(K) = JUMIN                                              
  14      IU(K+1) = IU(K) - LUK                                         
          IF (JUMAX.GT.LRATIO*(IU(K+1)-1)+1)  GO TO 112                 
          IM(K) = IJU(K) + LUK - 1                                      
C                                                                       
C  ******  CALCULATE NUMERICAL VALUES FOR KTH ROW  *********************
          SUM = B(K)                                                 
          I = N+1                                                       
  15      I = Q(I)                                                      
          IF (I.GE.K)  GO TO 18                                         
            AKI = - ROW(I)                                              
            SUM = SUM + AKI * TMP(I)                                    
            JMIN = IU(I+1) + 1                                          
            JMAX = IU(I)                                                
            IF (JMIN.GT.JMAX)  GO TO 17                                 
            MU = IJU(I) - JMIN                                          
            DO 16 J=JMIN,JMAX                                           
  16          ROW(JU(MU+J)) = ROW(JU(MU+J)) + AKI * U(J)                
  17        GO TO 15                                                    
C  ******  STORE VALUES IN TMP AND U  **********************************
  18      IF (ROW(K).EQ.0)  GO TO 108                                   
          DK = 1 / ROW(K)                                               
          TMP(K) = SUM * DK                                             
          JMIN = IU(K+1) + 1                                            
          JMAX = IU(K)                                                  
          IF (JMIN.GT.JMAX)  GO TO 20                                   
          MU = IJU(K) - JMIN                                            
          DO 19 J=JMIN,JMAX                                             
  19        U(J) = ROW(JU(MU+J)) * DK                                   
  20      CONTINUE                                                      
C                                                                       
C  ******  SOLVE  UX = TMP  BY BACK SUBSTITUTION  **********************
        K = N                                                           
        DO 23 I=1,N                                                     
          SUM = TMP(K)                                                  
          JMIN = IU(K+1) + 1                                            
          JMAX = IU(K)                                                  
          IF (JMIN.GT.JMAX)  GO TO 22                                   
          MU = IJU(K) - JMIN                                            
          DO 21 J=JMIN,JMAX                                             
  21        SUM = SUM - U(J) * TMP(JU(MU+J))                            
  22      TMP(K) = SUM                                                  
  23      K = K-1                                                       
        DO 24 K=1,N                                                     
  24      Z(K) = TMP(K)                                             
C                                                                       
        FLAG = 0                                                        
        ESP = IU(N+1) - 1 - (JUMAX - 1)/LRATIO                          
        RETURN                                                          
C                                                                       
C ** ERROR:  NULL ROW IN A                                              
 101    FLAG = N + K                                                 
        RETURN                                                          
C ** ERROR:  DUPLICATE ENTRY IN A                                       
 102    FLAG = 2*N + K                                               
        RETURN                                                          
C ** ERROR:  NULL PIVOT                                                 
 105    FLAG = 5*N + K                                                  
        RETURN                                                          
C ** ERROR:  ZERO PIVOT                                                 
 108    FLAG = 8*N + K                                                  
        RETURN                                                          
C ** ERROR:  INSUFFICIENT STORAGE FOR JU AND U                          
 112    FLAG = 12*N + K                                                 
C
        RETURN                                                          
        END                                                             

        
C______________________________________________________________________C  
      SUBROUTINE WarnErrSolv(WorE,routine,text,ierr)
C                       
C     $: Writes warning or error to error file
C______________________________________________________________________C
      
      CHARACTER*(*) routine,text
      CHARACTER*255 msgtype
      INTEGER WorE,ierr

C
C       Message, warning or error
C
      IF (WorE.EQ.0) THEN
      WRITE(msgtype,'(A7,A4,A,A4,A)') 
     &      "Message"," :: ",trim(routine)," :: ",trim(text)
      ELSE IF (WorE.EQ.5) THEN
       WRITE(msgtype,'(A7,A4,A,A4,A,I8)') 
     &      "Message"," :: ",trim(routine)," :: ",trim(text),ierr
      ELSE IF (WorE.EQ.1) THEN
       WRITE(msgtype,'(A7,A4,A,A4,A)') 
     &      "Warning"," :: ",trim(routine)," :: ",trim(text)
      ELSE IF (WorE.EQ.2) THEN
       WRITE(msgtype,'(A7,A4,A,A4,A)') 
     &      "Error  "," :: ",trim(routine)," :: ",trim(text)
      ELSE IF (WorE.EQ.3) THEN
       WRITE(msgtype,'(A7,A4,A,A4,A,I8)') 
     &      "Warning"," :: ",trim(routine)," :: ",trim(text),ierr
      ELSE IF (WorE.EQ.4) THEN
          WRITE(msgtype,'(A7,A4,A,A4,A,I8)') 
     &      "Error  "," :: ",trim(routine)," :: ",trim(text),ierr
       ELSE
         msgtype="Unkonwn"
        END IF
C
C       write to error file
C   
        WRITE(*,'(A)') trim(msgtype)
C
C       stop if error ocured
C     
      IF (WorE.EQ.2.OR.WorE.EQ.4) THEN
        STOP
      END IF

      END      
      
      
C----------------------------------------------------------------------C
c                                                                      c      
      SUBROUTINE slvlsqr2
     &   (m,n,nnz,maxit,eps4,nit,ierr,
     &    q,ia,ja,a,b,x)
c                                                                      c      
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Solve system of equations using Least Square iterative method    **
c ** - --                            -     --  -                      **
c ** !!! With built-in diagonal preconditioning !!!                   **
c ** a call to prelsqr2 must preceede a call to slvlsqr2              **
c **********************************************************************
c Parameters
      INTEGER  m,n,nnz,maxit,nit,ierr,
     &         ia(m+1),ja(nnz)
      REAL*8   eps,
     &         q(n),a(nnz),b(m),x(n)
c Internal variables
      LOGICAL  wantse,extra
      INTEGER  istop,iout,i
      REAL*8   x0,wm,u,se,v,w
      REAL*8   atol,btol,nrm_b,nrm_r
      REAL*8   zero,damp,conlim,anorm,acond,rnorm,arnorm,xnorm
      REAL*8   dnrm2
      REAL*4   eps4
c
      ALLOCATABLE x0(:),wm(:),u(:),se(:),v(:),w(:)
c
c Allocate memory for matrix and vectors
c
      ALLOCATE(x0(n),wm(m),u(m),se(1),v(n),w(n))
c
c Set initial values
c
      wantse=.FALSE.
      extra=.FALSE.
c
      zero=0.0D0
      damp=zero
      conlim=1.0D8
      eps=DBLE(eps4)
c
      IF (maxit.EQ.0) THEN
        maxit = 10*(m + n + 50)
        damp=zero
      END IF
c
      IF (extra) THEN
        iout=66
      ELSE
        iout=-1
      END IF
c
      CALL dcopy(n,x,1,x0,1)                       !x0=x
      wm=zero
      DO i=1,n
        x(i)=x(i)/q(i)                             !x=Q^-1*x
      END DO
      CALL aprod(1,m,n,nnz,x,wm,iA,jA,A)           !wm=A*x
      CALL dcopy(m,b,1,u,1)                        !
      CALL daxpy(m,-1.0D0,wm,1,u,1)                !u=u-wm=b-Ax
      nrm_r=dnrm2(m,u,1)
      nrm_b=dnrm2(m,b,1)
      atol=eps
      btol=eps*nrm_b/nrm_r
c      print *,btol,nrm_b,nrm_r
c
c LSQR
c
      CALL lsqr(
     &  m,n,nnz,damp,wantse,
     &  iA,jA,A,
     &  u,x,se,v,w,
     &  atol,btol,conlim,maxit,
     &  istop,nit,anorm,acond,rnorm,arnorm,xnorm)
c
c Preconditioninig
c
      CALL dvprod(n,Q,1,x,1)                         !x=Q*x
c
c Hotstart
c
      CALL daxpy(n,1.0D0,x0,1,x,1)                   !x=x0+dx
c
      IF (istop.GT.3) THEN
c       WRITE (*,*) 'Warning: LSQR: stoping condition ISTOP= ',istop
        ierr=istop
      ELSE
        ierr=0
      END IF
c
c Free memory
c
      DEALLOCATE(x0,wm,u,se,v,w)
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c      
      SUBROUTINE prelsqr2
     &   (m,n,nnz,
     &    q,ia,ja,a)
c                                                                      c      
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Solve system of equations using Least Square iterative method    **
c ** - --                            -     --  -                      **
c ** !!! With built-in diagonal preconditioning !!!                   **
c **********************************************************************
c Parameters
      INTEGER  m,n,nnz,
     &         ia(m+1),ja(nnz)
      REAL*8   q(n),a(nnz)
c Internal variables
      REAL*8   wm
c
      ALLOCATABLE wm(:)
c
c Allocate memory for matrix and vectors
c
      ALLOCATE(wm(m))
c
c Compute preconditioner
c
      CALL aprod(3,m,n,nnz,Q,wm,iA,jA,A)            !Q=M^-1=1/diag(sqrt(AT*A))
      CALL aprod(4,m,n,nnz,Q,wm,iA,jA,A)            !A=A*M^-1
c
c Free memory
c
      DEALLOCATE(wm)
c
      RETURN
      END        
