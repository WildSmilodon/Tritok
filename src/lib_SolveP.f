C^L           
C----------------------------------------------------------------------C
c----------------------------------------------------------------------c
c                                                                      c
c     PARALLEL SOLVER PROGRAM BLOCK                                    c
c                                                                      c
c----------------------------------------------------------------------c
C----------------------------------------------------------------------C
c $Date: 15/08/1998 $
c $Source: ~/bem/flow/libpsolve/psolve.f $
c $Revision: 0.3 $
c $State: Exp $
c
c                   Popacal JureR, Marca 2006 in Julija 2009
c
c **********************************************************************
c **                                                                  **
c **  SUBROUTINES USED :                                              **
c **                                                                  **
c **                   - pSolvEQN                                     **
c **                   - pFormPRE                                     **
c **                   - pSolvEQNfm                                   **
c **                   - pFormPREfm                                   **
c **                   - pSolvSLE                                     **
c **                   - pFormPRM                                     **
c **                                                                  **
c **********************************************************************
C^L           
C----------------------------------------------------------------------C
c
      SUBROUTINE pSolvEQN (env,slvt,pret,prep,neq,nnz,peq,pnz,npr,
     &                     iq,jq,dq,q,ia,ja,da,a,b,y,nit,ierr)

      USE inc_types
      TYPE(penv)  env
      INTEGER   slvt,pret,prep,neq,nnz,peq,pnz,npr,nit,ierr,
     &          iq(peq+1),jq(pnz),dq(peq),ia(peq+1),ja(pnz),da(peq)
      REAL*8    q(pnz),a(pnz),b(peq),y(neq)
c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  SOLVe system of EQuatioNs using various solvers                 **
c **  ----            --     -                                        **
c **********************************************************************
      LOGICAL   lexist,info
      INTEGER   mpr,inp,nnp,maxit,stopt
      REAL*8    eps,PARDNRM2,x(peq)
      EXTERNAL  pmatvec,pprenon,ppredia,ppreilu,PARDSUM,PARDNRM2,PPROGRESS
c
      DATA      maxit,stopt,eps,info/500, 5, 1.d-6, .TRUE./
c
C***  READ INITIAL VALUES
c
      INQUIRE (FILE='bem.slv',EXIST=lexist)
      IF (lexist) THEN
        OPEN (99,FILE='bem.slv',STATUS='OLD')
        READ (99,*) slvt,pret,maxit,stopt,eps
        CLOSE (99)
        IF (info) THEN
          CALL pgetmyproc(env,mpr)
          IF (mpr.EQ.1) THEN
            WRITE (*,'(A,I2,A,I2,A,I5,A,I2,A,E11.4)') 
     &        ' SolvEQN: Solver =',slvt,' Precon =',pret,
     &        ' Maxit =',maxit,' Stop =',stopt,' Eps =',REAL(eps)
          END IF
          info= .FALSE.
        END IF
      END IF
c
C***  CALCULATE PRECONDITIONER
c
      IF (prep.EQ.1.AND.pret.NE.0) THEN
        CALL dcopy(pnz,a,1,q,1)
        IF (pret.EQ.1) THEN
C***  diagonal
          CALL pformdia(env,neq,peq,pnz,npr,iq,jq,dq,q)
        ELSE IF (pret.EQ.2) THEN
C***  ilu
          CALL pformilu(env,neq,peq,pnz,npr,iq,jq,dq,q)
        ELSE
C***  error
          CALL WarnErrSolv(2,'pSolvEQN','Wrong preconditioner type!',0)
        END IF
      END IF
c
C***  RUN SOLVER
c
      CALL pgetnodes(env,inp,nnp)
      CALL dcopy(nnp,y(inp+1),1,x,1)
      IF (slvt.EQ.0) THEN
C***  direct
        IF (neq.NE.peq) THEN
          CALL WarnErrSolv(2,'pSolvEQN','Direct solver runs only on one processor!',0)
        END IF
        CALL solveysmp(neq,nnz,ia,ja,a,b,x,ierr)
        nit=1
      ELSE IF (slvt.EQ.1) THEN
C***  CGS
        IF (pret.EQ.0) THEN
          CALL pcgs(env,neq,peq,pnz,npr,pret,maxit,stopt,eps,nit,ierr,
     &              ia,ja,da,q,ia,ja,a,b,x,
     &              pmatvec,pprenon,pprenon,PARDSUM,PARDNRM2,PPROGRESS)
        ELSE IF (pret.EQ.1) THEN
          CALL pcgs(env,neq,peq,pnz,npr,pret,maxit,stopt,eps,nit,ierr,
     &              ia,ja,da,q,ia,ja,a,b,x,
     &              pmatvec,ppredia,ppredia,PARDSUM,PARDNRM2,PPROGRESS)
        ELSE IF (pret.EQ.2) THEN
          CALL pcgs(env,neq,peq,pnz,npr,pret,maxit,stopt,eps,nit,ierr,
     &              ia,ja,da,q,ia,ja,a,b,x,
     &              pmatvec,ppreilu,ppreilu,PARDSUM,PARDNRM2,PPROGRESS)
        END IF
        CALL pgetvec(neq,peq,npr,env,x,y)
      ELSE IF (slvt.EQ.2) THEN
C***  RBi-CGSTAB
        IF (pret.EQ.0) THEN
          CALL prbicgstab(env,neq,peq,pnz,npr,pret,maxit,stopt,eps,nit,ierr,
     &                    ia,ja,da,q,ia,ja,a,b,x,
     &                    pmatvec,pprenon,pprenon,PARDSUM,PARDNRM2,PPROGRESS)
        ELSE IF (pret.EQ.1) THEN
          CALL prbicgstab(env,neq,peq,pnz,npr,pret,maxit,stopt,eps,nit,ierr,
     &                    ia,ja,da,q,ia,ja,a,b,x,
     &                    pmatvec,ppredia,ppredia,PARDSUM,PARDNRM2,PPROGRESS)
        ELSE IF (pret.EQ.2) THEN
          CALL prbicgstab(env,neq,peq,pnz,npr,pret,maxit,stopt,eps,nit,ierr,
     &                    ia,ja,da,q,ia,ja,a,b,x,
     &                    pmatvec,ppreilu,ppreilu,PARDSUM,PARDNRM2,PPROGRESS)
        END IF
        CALL pgetvec(neq,peq,npr,env,x,y)
      ELSE
C***  error
        CALL WarnErrSolv(2,'pSolvEQN','wrong solver type',0)
      END IF
c
      RETURN
      END
C^L           
C----------------------------------------------------------------------C
c
      SUBROUTINE pFormPRE (env,pret,neq,peq,pnz,npr,
     &                     iq,jq,dq,q,ia,ja,da,a,ierr)

      USE inc_types
      TYPE(penv)  env
      INTEGER   pret,neq,peq,pnz,npr,ierr,
     &          iq(peq+1),jq(pnz),dq(peq),ia(peq+1),ja(pnz),da(peq)
      REAL*8    q(pnz),a(pnz)
c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  FORM PREconditioning Matrix Q                                   **
c **  ---- ---                                                        **
c **********************************************************************
      LOGICAL   lexist
      INTEGER   dummy
      EXTERNAL  pprenon,ppredia,ppreilu
c
C***  READ INITIAL VALUES
c
      ierr=0
      INQUIRE (FILE='bem.slv',EXIST=lexist)
      IF (lexist) THEN
        OPEN (99,FILE='bem.slv',STATUS='OLD')
        READ (99,*) dummy,pret
        CLOSE (99)
c       IF (info) THEN
c         WRITE (*,'(AI2)') ' FormPRE: Precon =',pret
c         info= .FALSE.
c       END IF
      END IF
c
C***  CALCULATE PRECONDITIONER
c
      IF (pret.NE.0) THEN
C***  copy A to Q
        CALL dcopy(pnz,a,1,q,1)
        IF (pret.EQ.1) THEN
C***  diagonal
          CALL pformdia(env,neq,peq,pnz,npr,iq,jq,dq,q)
        ELSE IF (pret.EQ.2) THEN
C***  ilu
          CALL pformilu(env,neq,peq,pnz,npr,iq,jq,dq,q)
        ELSE
C***  error
          CALL WarnErrSolv(2,'pFormPRE','Wrong preconditioner type!',0)
        END IF
      END IF
c
      RETURN
      END
C^L           
C----------------------------------------------------------------------C
c
      SUBROUTINE pSolvEQNfm (env,slvt,pret,prep,maxit,stopt,eps,neq,peq,
     &                       q,a,b,y,nit,cpu,ierr)

      USE inc_types
	TYPE(penv)  env
      INTEGER   slvt,pret,prep,maxit,stopt,neq,peq,nit,ierr
      REAL(4)   eps,cpu
      REAL(8)   q(peq),a(peq,neq),b(peq),y(neq)
c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  SOLVe system of EQuatioNs using various solvers (Full Matrix)   **
c **  ----            --     -                                        **
c **********************************************************************
c Internal
      INTEGER   npr
      REAL      cpu0,cptime
      REAL*8    PARDNRM2,x(peq)
      EXTERNAL  pmatvecfm,pprenon,pprediafm,PARDSUM,PARDNRM2,PPROGRESS
c Dummy
      INTEGER   pnz,ia,ja,da
c
C***  Measure CPU time
c
      cpu0=cptime(0.)
c
C***  Read environment
c
      npr=env%nproc
c
C***  CALCULATE PRECONDITIONER
c
      IF (prep.EQ.1.AND.pret.NE.0) THEN
        IF (pret.EQ.1) THEN
C***  diagonal
          CALL pformdiafm(neq,peq,env%zac,a,q)
        ELSE IF (pret.EQ.2) THEN
C***  error
          CALL WarnErrSolv(2,'pSolvEQNfm','Wrong preconditioner type!',0)
        END IF
      END IF
c
C***  RUN SOLVER
c
C***  initial guess
      x(1:env%nmn)=y(env%zac:env%kon)
c
      IF (slvt.EQ.0) THEN
C***  direct
        IF (neq.NE.peq) THEN
          CALL WarnErrSolv(2,'pSolvEQNfm','direct solver runs only on one processor!',0)
        END IF
        x=b
        CALL dgesl(a,neq,neq,q,x,0)
	  y=x
        nit=1
      ELSE IF (slvt.EQ.1) THEN
C***  CGS
        IF (pret.EQ.0) THEN
          CALL pcgs(env,neq,peq,pnz,npr,pret,maxit,stopt,DBLE(eps),
     &              nit,ierr,ia,ja,da,q,ia,ja,a,b,x,
     &              pmatvecfm,pprenon,pprenon,PARDSUM,PARDNRM2,PPROGRESS)
        ELSE IF (pret.EQ.1) THEN
          CALL pcgs(env,neq,peq,pnz,npr,pret,maxit,stopt,DBLE(eps),
     &              nit,ierr,ia,ja,da,q,ia,ja,a,b,x,
     &              pmatvecfm,pprediafm,pprediafm,PARDSUM,PARDNRM2,PPROGRESS)
        END IF
        CALL pgetvec(neq,peq,npr,env,x,y)
      ELSE IF (slvt.EQ.2) THEN
C***  RBi-CGSTAB
        IF (pret.EQ.0) THEN
          CALL prbicgstab(env,neq,peq,pnz,npr,pret,maxit,stopt,DBLE(eps),
     &              nit,ierr,ia,ja,da,q,ia,ja,a,b,x,
     &              pmatvecfm,pprenon,pprenon,PARDSUM,PARDNRM2,PPROGRESS)
        ELSE IF (pret.EQ.1) THEN
          CALL prbicgstab(env,neq,peq,pnz,npr,pret,maxit,stopt,DBLE(eps),
     &              nit,ierr,ia,ja,da,q,ia,ja,a,b,x,
     &              pmatvecfm,pprediafm,pprediafm,PARDSUM,PARDNRM2,PPROGRESS)
        END IF
        CALL pgetvec(neq,peq,npr,env,x,y)
      ELSE
C***  error
          CALL WarnErrSolv(2,'pSolvEQN','Wrong solver type!',0)
      END IF
c
      cpu=cpu+cptime(cpu0)
c
      RETURN
      END
C^L           
C----------------------------------------------------------------------C
c
      SUBROUTINE pFormPREfm (env,pret,neq,peq,q,a,cpu,ierr)

      USE inc_types
	TYPE(penv)  env
      INTEGER   pret,neq,peq,ierr
      REAL(4)   cpu
      REAL(8)   q(peq),a(peq,neq)
c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  FORM PREconditioning Matrix Q                                   **
c **  ---- ---                                                        **
c **********************************************************************
c Internal
      REAL(4)   cpu0,cptime
      EXTERNAL  pprenon,ppredia,ppreilu
c
c Measure CPU time
c
      cpu0=cptime(0.)
c
C***  CALCULATE PRECONDITIONER
c
      ierr=0
      IF (pret.NE.0) THEN
C***  copy A to Q
        IF (pret.EQ.1) THEN
C***  diagonal
          CALL pformdiafm(neq,peq,env%zac,a,q)
        ELSE IF (pret.EQ.2.AND.env%nproc.EQ.1) THEN
c lu
          CALL dgefa(a,neq,neq,q,ierr)
        ELSE
C***  error
          CALL WarnErrSolv(2,'pFormPREfm','Wrong preconditioner type!',0)
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
      SUBROUTINE pSolvSLE (env,slvt,pret,prep,maxit,stopt,eps,
     &                     neq,peq,pnz,npr,q,a,b,y,nit,cpu,ierr)
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
c ** peq  ... number of equations on one processor                    **
c ** pnz  ... number of non-zeroes on one processor                   **
c ** npr  ... number of processors                                    **
c ** env  ... parallel environment variable                           **
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
      USE inc_types
c      INCLUDE 'FemInc.f'
	TYPE(penv)  env
c Arguments
      INTEGER  slvt,pret,prep,maxit,stopt,neq,peq,pnz,npr,nit,ierr
      REAL     eps,cpu
      TYPE     (MATRIX) q,a
      REAL(8)  b(peq),x(peq),y(neq)
c Internal
      REAL     cpu0,cptime
      REAL(8)  PARDNRM2
      EXTERNAL pmatsvec,pmatvec,pprenon,ppredia,ppreilu,PARDSUM,PARDNRM2,PPROGRESS
c
c Measure CPU time
c
      cpu0=cptime(0.)
c
c Calculate preconditioner
c
      IF (prep.EQ.1.AND.pret.NE.0) THEN
        q%v=a%v
        IF (pret.EQ.1) THEN
c diagonal
          CALL pformdia(env,neq,peq,pnz,npr,q%i,q%j,q%d,q%v)
        ELSE IF (pret.EQ.2) THEN
c ilu
          CALL pformilu(env,neq,peq,pnz,npr,q%i,q%j,q%d,q%v)
        ELSE
c error
          CALL WarnErrSolv(2,'pSolvSLE','Wrong preconditioner type!',0)
        END IF
      END IF
c
c Run solver
c
c Initial guess
      x(1:env%nmn)=y(env%zac:env%kon)
c
      IF (slvt.EQ.0) THEN
c Direct solver
        IF (neq.NE.peq) THEN
          CALL WarnErrSolv(2,'pSolvSLE','direct solver runs only on one processor!',0)
        END IF
        CALL solveysmp(peq,pnz,a%i,a%j,a%v,b,x,ierr)
        nit=1
      ELSE IF (slvt.EQ.1) THEN
c CGS
        IF (pret.EQ.0) THEN
          CALL pcgs(env,neq,peq,pnz,npr,pret,maxit,stopt,DBLE(eps),
     &              nit,ierr,q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &              pmatsvec,pprenon,pprenon,PARDSUM,PARDNRM2,PPROGRESS)
        ELSE IF (pret.EQ.1) THEN
          CALL pcgs(env,neq,peq,pnz,npr,pret,maxit,stopt,DBLE(eps),
     &              nit,ierr,q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &              pmatsvec,ppredia,ppredia,PARDSUM,PARDNRM2,PPROGRESS)
        ELSE IF (pret.EQ.2) THEN
          CALL pcgs(env,neq,peq,pnz,npr,pret,maxit,stopt,DBLE(eps),
     &              nit,ierr,q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &              pmatsvec,ppreilu,ppreilu,PARDSUM,PARDNRM2,PPROGRESS)
        ELSE
          CALL WarnErrSolv(1,'SolvSLE','Wrong preconditioner type',pret)
          CALL pcgs(env,neq,peq,pnz,npr,pret,maxit,stopt,DBLE(eps),
     &              nit,ierr,q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &              pmatsvec,pprenon,pprenon,PARDSUM,PARDNRM2,PPROGRESS)
        END IF
        CALL pgetvec(neq,peq,npr,env,x,y)
      ELSE IF (slvt.EQ.2) THEN
c RBi-CGSTAB
        IF (pret.EQ.0) THEN
          CALL prbicgstab(env,neq,peq,pnz,npr,pret,maxit,stopt,DBLE(eps),
     &                    nit,ierr,q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &                    pmatsvec,pprenon,pprenon,PARDSUM,PARDNRM2,PPROGRESS)
        ELSE IF (pret.EQ.1) THEN
          CALL prbicgstab(env,neq,peq,pnz,npr,pret,maxit,stopt,DBLE(eps),
     &                    nit,ierr,q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &                    pmatsvec,ppredia,ppredia,PARDSUM,PARDNRM2,PPROGRESS)
        ELSE IF (pret.EQ.2) THEN
          CALL prbicgstab(env,neq,peq,pnz,npr,pret,maxit,stopt,DBLE(eps),
     &                    nit,ierr,q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &                    pmatsvec,ppreilu,ppreilu,PARDSUM,PARDNRM2,PPROGRESS)
        ELSE
          CALL WarnErrSolv(1,'SolvSLE','Wrong preconditioner type',pret)
          CALL prbicgstab(env,neq,peq,pnz,npr,pret,maxit,stopt,DBLE(eps),
     &                    nit,ierr,q%i,q%j,q%d,q%v,a%i,a%j,a%v,b,x,
     &                    pmatsvec,pprenon,pprenon,PARDSUM,PARDNRM2,PPROGRESS)
        END IF
	  CALL pgetvec(neq,peq,npr,env,x,y)
      ELSE
c Error
        IF (env%mpr.EQ.0) THEN
          CALL WarnErrSolv(2,'pSolvSLE','Wrong solver type!',0)
	  END IF
	  CALL pquitpenv
        STOP
      END IF
c
      IF (env%mpr.EQ.0) THEN
        IF (ierr.EQ.-1) THEN
          CALL WarnErrSolv(1,'pSolvSLE','Solver reached maximum number of iterations!',0)
c        ELSE IF (ierr.LE.-2) THEN
c          CALL WarnErrSolv(3,'pSolvSLE','Solver finished with error',ierr)
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
      SUBROUTINE pFormPRM (env,pret,neq,peq,pnz,npr,q,a,cpu,ierr)
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
c ** peq  ... number of equations on one processor                    **
c ** pnz  ... number of non-zeroes on one processor                   **
c ** q    ... preconditioner matrix Q                                 **
c ** env  ... parallel environment variable                           **
c ** a    ... matrix A                                                **
c ** cpu .... elapsed CPU time                                        **
c ** ierr ... error status (0=no error)                               **
c **                                                                  **
c **********************************************************************
c Type definition
      USE inc_types
c      INCLUDE 'FemInc.f'
	TYPE(penv)  env
c Arguments
      INTEGER pret,neq,peq,pnz,npr,ierr
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
      ierr=0
      IF (pret.NE.0) THEN
        q%v=a%v
        IF (pret.EQ.1) THEN
c diagonal
          CALL pformdia(env,neq,peq,pnz,npr,q%i,q%j,q%d,q%v)
        ELSE IF (pret.EQ.2) THEN
c ilu
          CALL pformilu(env,neq,peq,pnz,npr,q%i,q%j,q%d,q%v)
        ELSE
c error
          IF (env%mpr.EQ.0) THEN
	      CALL WarnErrSolv(2,'pFormPRM','Wrong preconditioner type!',0)
	    END IF
	    CALL pquitpenv
          STOP
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
c     PARALLEL INTERFACE ROUTINES PROGRAM BLOCK                        c
c                                                                      c
c----------------------------------------------------------------------c
C----------------------------------------------------------------------C
c $Date: 15/08/1998 $
c $Source: ~/bem/flow/libpsolve/parallel.f $
c $Revision: 0.3 $
c $State: Exp $
c **********************************************************************
c **                                                                  **
c **  SUBROUTINES USED :                                              **
c **                                                                  **
c **                   - pinitpenv                                    **
c **                   - pgetcomm                                     **
c **                   - pgetnproc                                    **
c **                   - pgetmyproc                                   **
c **                   - pgetmacro                                    **
c **                   - pgetnodes                                    **
c **                   - pgetmatrix                                   **
c **                   - pgetvec                                      **
c **                   - pgetsvec                                     **
c **                   - pgetsvec2                                    **
c **                   - pgetsvec3                                    **
c **                   - pgetarg                                      **
c **                   - pquitpenv                                    **
c **                                                                  **
c **********************************************************************     
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pinitpenv (env)

      USE inc_types
      TYPE(penv)  env
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  INITialize Parallel ENVironment                                 **
c **  ----       -        ---                                         **
c **********************************************************************
c      INCLUDE 'mpif.h'
      INTEGER comm,npr,mpr,ierr
c
      comm = MPI_COMM_WORLD
c
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(comm,npr,ierr)
      CALL MPI_COMM_RANK(comm,mpr,ierr)
c
      env%comm = comm
      env%nproc = npr
      env%myproc = mpr+1
      env%mpr = mpr
c
      RETURN
      END SUBROUTINE
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pgetcomm (env,comm)

      USE inc_types
      TYPE(penv)  env
      INTEGER comm
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  GET COMMunicator number                                         **
c **  --- ----                                                        **
c **********************************************************************
c
      comm = env%comm
c
      RETURN
      END SUBROUTINE
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pgetnproc (env,nproc)

      USE inc_types
      TYPE(penv)  env
      INTEGER nproc
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  GET Number of PROCesses                                         **
c **  --- -         ----                                              **
c **********************************************************************
c
      nproc = env%nproc
c
      RETURN
      END SUBROUTINE
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pgetmyproc (env,myproc)

      USE inc_types
      TYPE(penv)  env
      INTEGER myproc
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  GET MY PROCess number                                           **
c **  --- -- ----                                                     **
c **********************************************************************
c
      myproc = env%myproc
c
      RETURN
      END SUBROUTINE
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pgetmacro (env,ima,nma)

      USE inc_types
      TYPE(penv)  env
      INTEGER ima,nma
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  GET MACRO distribution from parallel environment vector         **
c **  --- -----                                                       **
c **********************************************************************
      INTEGER mpr
c
      mpr = env%myproc
c
      ima = env%imac(mpr)
      nma = env%nmac(mpr)
c
      RETURN
      END SUBROUTINE
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pgetnodes (env,inp,nnp)

      USE inc_types
      TYPE(penv)  env
      INTEGER inp,nnp
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  GET NODES distribution from parallel environment vector         **
c **  --- -----                                                       **
c **********************************************************************
      INTEGER mpr
c
      mpr = env%myproc
c
      inp = env%inod(mpr)
      nnp = env%nnod(mpr)
c
      RETURN
      END SUBROUTINE
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pgetmatrix (env,imx,nmx)

      USE inc_types
      TYPE(penv)  env
      INTEGER imx,nmx
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  GET MATRIX distribution from parallel environment vector        **
c **  --- ------                                                      **
c **********************************************************************
      INTEGER mpr
c
      mpr = env%myproc
c
      imx = env%imat(mpr)
      nmx = env%nmat(mpr)
c
      RETURN
      END SUBROUTINE
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pgetsvec (neq,peq,env,x,y)

      USE inc_types
      TYPE(penv)  env
      INTEGER neq,peq
      REAL*8  x(peq),y(neq)
	REAL*8 poslji,prejmi
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  GET whole Sparce VECtor                                         **
c **  ---       _      ---                                            **
c **********************************************************************
c      INCLUDE 'mpif.h'
      INTEGER mpr,comm,ierr,i,j,stat(MPI_STATUS_SIZE),tag
      INTEGER, POINTER :: inod(:),nnod(:)

	ALLOCATABLE poslji(:),prejmi(:)
c
C***  GET ENVIRONMENT VARIABLES
c
      comm = env%comm
      mpr  = env%myproc
      inod => env%inod
      nnod => env%nnod
  
      tag=1

      y=0.00D00
c
C***  Copy my nodes to output vector y
c
	DO i=1,peq
	  y(inod(mpr)+i)=x(i)
	END DO
c
C***  Loop through all my communication pairs
c
      DO i=1,env%CLstkom
c
C***    Allocate send and recieve vectors
c 
        ALLOCATE (poslji(env%CLste1(i)),prejmi(env%CLste2(i)))
c
C***    Compose send vector
c 
	  DO j=1,env%CLste1(i)
	    poslji(j)=x(env%CLlist(env%CLzac1(i)-1+j)-env%inod(mpr))
	  END DO
c
C***    Send POSLJI vector and recieve PREJMI vector
c 
        CALL MPI_SENDRECV(poslji,env%CLste1(i),MPI_DOUBLE_PRECISION,env%CLp2(i)-1,tag,
     &                    prejmi,env%CLste2(i),MPI_DOUBLE_PRECISION,env%CLp2(i)-1,tag,
     &                    comm,stat,ierr)
c
C***    Distribute recieved vector
c 
	  DO j=1,env%CLste2(i)
	    y(env%CLlist(env%CLzac2(i)-1+j))=prejmi(j)
	  END DO
c
C***    Free memory
c 
        DEALLOCATE (poslji,prejmi)
	END DO

	
	RETURN
      END SUBROUTINE


C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE SizeCommLists(env,             ! parallel environment
     &                         zvr,sto,nnz,neq, ! matrix structure
     &                         CLst,CLstkom)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Size of Communication List for matrix * vector multplication    **
c **  ----- -         -                         ---                   **
c **********************************************************************
c      INCLUDE 'mpif.h'
      USE inc_types
      TYPE(penv) :: env


	INTEGER i,l,o,proc,col,icol
      INTEGER nnz,neq
      INTEGER Sto(nnz),Zvr(neq+1),kateri(neq)


	INTEGER CLst,CLstkom

      CLstkom=env%nproc*(env%nproc-1)/2

      CLst=0
c     gremo po procesorjih
      DO i=1,env%nproc
c       vrstica swithcov, ki pove katere komuniciramo (da se resimo podvojenih)
	  kateri=0
c       po vrsticah znotraj enega procesorja
        DO l=env%inod(i)+1,env%inod(i)+env%nnod(i)
c         po vrsti
          DO icol=zvr(l),zvr(l+1)-1  
	      col=sto(icol)
c           po procesorjih ki oddajajo
            DO proc=1,i-1
              IF (col.GE.env%inod(proc)+1.AND.col.LE.env%inod(proc)+env%nnod(proc)) THEN
c               nasel! procesor PROC poslje clen col procesurju I
                kateri(col)=1
	        END IF
	      END DO
c           po procesorjih ki sprejemajo	     
            DO proc=i+1,env%nproc
              IF (col.GE.env%inod(proc)+1.AND.col.LE.env%inod(proc)+env%nnod(proc)) THEN
c               nasel! procesor PROC poslje clen col procesurju I
                kateri(col)=1
	        END IF
	      END DO
 	    END DO
	  END DO
c       stevilo za komunikacijo na procesor i
        DO o=1,neq
	    CLst=CLst+kateri(o)
	  END DO 	
	END DO

	RETURN 
	END


C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pgetvec (neq,peq,npr,env,x,y)

      USE inc_types
      TYPE(penv)  env
      INTEGER neq,peq,npr
      REAL*8  x(peq),y(neq)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  GET whole VECtor                                                **
c **  ---       ---                                                   **
c **********************************************************************
c      INCLUDE 'mpif.h'
      INTEGER mpr,comm,ierr
      INTEGER, POINTER :: inod(:),nnod(:)
c
C***  GET ENVIRONMENT VARIABLES
c
      comm = env%comm
      mpr  = env%myproc
      inod => env%inod
      nnod => env%nnod
c
C***  GATHER x IN y
c
      CALL MPI_ALLGATHERV(x,peq,MPI_DOUBLE_PRECISION,
     &                    y,nnod,inod,MPI_DOUBLE_PRECISION,comm,ierr)
c
      RETURN
      END SUBROUTINE
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pgetsvec2 (neq,peq,nex,nsipex,env,x,y)

      USE inc_types
      TYPE(penv)  env
      INTEGER neq,peq,nex,nsipex
      REAL*8  x(peq),y(neq)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  GET Sparse VECtor from processes in group                       **
c **  --- -      ---                                                  **
c **********************************************************************
c      INCLUDE 'mpif.h'
      INTEGER  mpr,comm,stat(MPI_STATUS_SIZE),ierr,ib,pr,bi,bo,xb,yb
      INTEGER, POINTER :: nmac(:),inod(:),prop(:),beop(:),bein(:)
c
C***  GET ENVIRONMENT VARIABLES
c
      mpr  =  env%myproc
      comm =  env%comm
      nmac => env%nmac
      inod => env%inod
      prop => env%prop
      beop => env%beop
      bein => env%bein
c
C***  FIND CONNECTIVITY POINTS BETWEEN LOCAL MATRIX AND GLOBAL VECTOR
c
      DO ib = 1,nex*nmac(mpr)
        pr = prop(ib)
        IF (pr.EQ.0) GOTO 10
        bi = bein(ib)
        bo = beop(ib)
        xb = nsipex*(bi-1)+1
        yb = inod(pr)+nsipex*(bo-1)+1
c       IF (mpr.LE.0) THEN
c         WRITE(*,*) 'MYPROC  = ',mpr,ib,pr,bi,bo,xb,yb,nex*nmac(mpr)
c       END IF
        CALL MPI_SENDRECV(x(xb),nsipex,MPI_DOUBLE_PRECISION,pr-1,bi,
     &                    y(yb),nsipex,MPI_DOUBLE_PRECISION,pr-1,bo,
     &                    comm,stat,ierr)
    5 END DO
c     IF (mpr.LE.0) THEN
c       WRITE(*,*) 'MYPROC  = ',mpr
c       WRITE(*,'(8E9.3)') (x(i),i = 1,peq)
c       WRITE(*,'(8E9.3)') (y(i),i = 1,neq)
c     END IF
c
   10 RETURN
      END SUBROUTINE
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pgetsvec4 (neq,peq,nex,nsipex,env,x,y)

      USE inc_types
      TYPE(penv)  env
      INTEGER neq,peq,nex,nsipex
      REAL*8  x(peq),y(neq)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  GET Sparse VECtor from processes in group                       **
c **  --- -      ---                                                  **
c **********************************************************************
c      INCLUDE 'mpif.h'
      INTEGER  mpr,comm,stat(MPI_STATUS_SIZE),ierr,ib,pr,bi,bo,xb,yb
      INTEGER, POINTER :: nmac(:),inod(:),prop(:),beop(:),bein(:)
c
C***  GET ENVIRONMENT VARIABLES
c
      mpr  =  env%myproc
      comm =  env%comm
      nmac => env%nmac
      inod => env%inod
      prop => env%prop
      beop => env%beop
      bein => env%bein
c
C***  FIND CONNECTIVITY POINTS BETWEEN LOCAL MATRIX AND GLOBAL VECTOR
c
      DO ib=1,nex*nmac(mpr)
        pr=prop(ib)
        IF (pr.EQ.0) GOTO 10
        bi=bein(ib)
        xb=nsipex*(bi-1)+1
        CALL MPI_SEND(x(xb),nsipex,MPI_DOUBLE_PRECISION,pr-1,bi,
     &                comm,ierr)
    5 END DO
c
   10 DO ib=1,nex*nmac(mpr)
        pr=prop(ib)
        IF (pr.EQ.0) GOTO 20
        bo=beop(ib)
        yb=inod(pr)+nsipex*(bo-1)+1
        CALL MPI_RECV(y(yb),nsipex,MPI_DOUBLE_PRECISION,pr-1,bo,
     &                comm,stat,ierr)
   15 END DO
c
   20 RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pgetsvec3 (neq,peq,nsipex,env,x,y)

      USE inc_types
      TYPE(penv)  env
      INTEGER neq,peq,nsipex
      REAL*8  x(peq),y(neq)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  GET Sparse VECtor from processes in group                       **
c **  --- -      ---                                                  **
c **********************************************************************
c      INCLUDE 'mpif.h'
      INTEGER  mpr,ncp,ncex,comm,ierr,stat(MPI_STATUS_SIZE),
     &         ip,ie,je,pr,bi,bo,b,xb,yb,ncn
      INTEGER, POINTER :: inp(:),ice(:),jce(:),prop(:),beop(:),bein(:)
      REAL*8,  POINTER :: send(:),recv(:)
c
C***  GET ENVIRONMENT VARIABLES
c
      mpr  =  env%myproc
      comm =  env%comm
      inp  => env%inod
      ncp  =  env%ncp
      ncex =  env%ncex
      ice  => env%ice
      jce  => env%jce
      prop => env%prop
      beop => env%beop
      bein => env%bein
c
C***  FIND CONNECTIVITY POINTS BETWEEN LOCAL MATRIX AND GLOBAL VECTOR
c
      DO ip=1,ncp
        ncn=(ice(ip+1)-ice(ip))*nsipex
        ALLOCATE(send(ncn),recv(ncn))
        b=1
        DO ie=ice(ip),ice(ip+1)-1
          je=jce(ie)
          pr=prop(je)
          bi=bein(je)
          xb=nsipex*(bi-1)+1
          CALL dcopy(nsipex,x(xb),1,send(b),1)
          b=b+nsipex
        END DO
        pr=prop(jce(ice(ip)))-1
        CALL MPI_SENDRECV(send,ncn,MPI_DOUBLE_PRECISION,pr,mpr-1,
     &                    recv,ncn,MPI_DOUBLE_PRECISION,pr,pr,
     &                    comm,stat,ierr)
        b=1
        DO ie=ice(ip),ice(ip+1)-1
          je=jce(ie)
          pr=prop(je)
          bo=beop(je)
          yb=inp(pr)+nsipex*(bo-1)+1
          CALL dcopy(nsipex,recv(b),1,y(yb),1)
          b=b+nsipex
        END DO
        DEALLOCATE(send,recv)
      END DO
      IF (mpr.LE.0) THEN
c        WRITE(*,*) 'MYPROC =',mpr
        CALL PPRINTV(neq,y)
      END IF
c
   20 RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pgetarg (env,narg,fileA,fileB,fileX)

      USE inc_types
      TYPE(penv)  env
      CHARACTER narg*(4)
      CHARACTER fileA*(20),fileB*(20),fileX*(20)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  GET command line ARGuments                                      **
c **  ---                    ---                                      **
c **********************************************************************
c      INCLUDE 'mpif.h'
      INTEGER comm
c
      comm = env%comm
c
c     CALL MPI_BCAST(narg,1,MPI_INTEGER,0,comm,ierr)
c     CALL MPI_BCAST(fileA,20,MPI_CHARACTER,0,comm,ierr)
c     CALL MPI_BCAST(fileB,20,MPI_CHARACTER,0,comm,ierr)
c     CALL MPI_BCAST(fileX,20,MPI_CHARACTER,0,comm,ierr)
c
      RETURN
      END SUBROUTINE
c
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pquitpenv ()
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  QUIT Parallel ENVironment                                       **
c **  ---- -        ---                                               **
c **********************************************************************
      INTEGER ierr
c
      CALL MPI_FINALIZE(ierr)
c
      RETURN
      END SUBROUTINE
C^L
C----------------------------------------------------------------------C
c----------------------------------------------------------------------c
c                                                                      c
c     PARALLEL PRECONDITIONER FORMATION PROGRAM BLOCK                  c
c                                                                      c
c----------------------------------------------------------------------c
C----------------------------------------------------------------------C
c $Date: 15/08/1998 $
c $Source: ~/bem/flow/libpsolve/pprecon.f $
c $Revision: 0.3 $
c $State: Exp $
c **********************************************************************
c **                                                                  **
c **  SUBROUTINES USED :                                              **
c **                                                                  **
c **                   - pformdia                                     **
c **                   - pformdiafm                                   **
c **                   - pformilu                                     **
c **                                                                  **
c **********************************************************************
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pformdia (env,neq,peq,pnz,npr,iq,jq,dq,q)

      USE inc_types
	TYPE(penv)  env
      INTEGER neq,peq,pnz,npr,iq(peq+1),jq(pnz),dq(peq)
      REAL*8  q(pnz)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  FORM preconditioner with DIAgonal matrix                        **
c **  ----                     ---                                    **
c **********************************************************************
      INTEGER i
c
      DO i = 1,peq
        q(dq(i)) = 1.D0/q(dq(i))
      END DO
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pformdiafm (neq,peq,zac,a,q)

      INTEGER neq,peq,zac
      REAL(8) a(peq,neq),q(peq)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  FORM preconditioner with DIAgonal matrix (Full Matrix)          **
c **  ----                     ---                                    **
c **********************************************************************
      INTEGER i
c
      DO i = 1,peq
        q(i) = 1.D0/a(i,zac+i-1)
      END DO
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pformilu (env,neq,peq,pnz,npr,iq,jq,dq,q)

      USE inc_types
	TYPE(penv)  env
      INTEGER neq,peq,pnz,npr,iq(peq+1),jq(pnz),dq(peq)
      REAL*8  q(pnz)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  FORM preconditioner with Incomplete LU decomposed matrix        **
c **  ----                     -          --                          **
c **********************************************************************
      INTEGER jj,i,j,k,l,m,inp,nnp,jbeg,jend
      ALLOCATABLE jj(:)
      ALLOCATE(jj(neq))
c
C***  INITIALIZE INDEX VECTOR
c
      jj=0
c
C***  GET LOCAL NODES
c
      CALL pgetnodes(env,inp,nnp)
c
C***  LOOP OVER THE ROWS
c
      DO i=1,peq
        jbeg=iq(i)+1
        jend=iq(i+1)-1
c
C***  BUILD INDEX VECTOR FOR CURRENT ROW
c
        DO j=jbeg,jend
          jj(jq(j))=j
        END DO
c
C***  LOOP OVER THE ROWS ABOW THE CURRENT ROW
c
        DO j=iq(i),dq(i)-1
          m=jq(j)-inp
c
C***  CHECK IF M IS OUTSIDE OF LOCAL PART
c
          IF (m.GT.0) THEN
c
C***  RHS FACTOR, q_ik = q_ik/q_kk
c
            q(j)=q(j)/q(dq(m))
c
C***  ELIMINATION, q_ij = q_ij - (q_ik/q_kk)*q_kj
c
            DO l=dq(m)+1,iq(m+1)-1
              k=jj(jq(l))
              IF (k.EQ.0) GOTO 10
              q(k)=q(k)-q(j)*q(l)
   10       END DO
          ELSE
c
C***  ZERO IF M IS OUTSIDE OF LOCAL PART
c
            q(j)=0.d0
          END IF
        END DO
c
C***  CLEAR INDEX VECTOR
c
        DO j=jbeg,jend
          jj(jq(j))=0
        END DO
      END DO
c
C***  FREE MEMORY
c
      DEALLOCATE(jj)
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c----------------------------------------------------------------------c
c                                                                      c
c     PARALLEL CGS SOLVER PROGRAM BLOCK                                c
c                                                                      c
c----------------------------------------------------------------------c
C----------------------------------------------------------------------C
c $Date: 15/08/1998 $
c $Source: ~/bem/flow/libpsolve/pcgs.f $
c $Revision: 0.3 $
c $State: Exp $
c **********************************************************************
c **                                                                  **
c **  SUBROUTINES USED :                                              **
c **                                                                  **
c **                   - pcgs                                         **
c **                   - pimdcgs                                      **
c **                                                                  **
c **********************************************************************
C^L
C----------------------------------------------------------------------C
c                                                                      c      
      SUBROUTINE pCGS
     &   (env,neq,peq,pnz,npr,pret,maxit,stopt,eps,nit,ierr,
     &    iq,jq,dq,q,ia,ja,a,b,x,
     &    matvec,preconl,preconr,PARDSUM,pdnrm,PPROGRESS)

      USE inc_types
      TYPE(penv)  env
      INTEGER  neq,peq,pnz,npr,pret,maxit,stopt,nit,ierr,
     &         iq(peq+1),jq(pnz),dq(peq),ia(peq+1),ja(pnz),da(peq)
      REAL*8   q(pnz),a(pnz),b(peq),x(peq),eps,pdnrm
      EXTERNAL matvec,preconl,preconr,PARDSUM,pdnrm,PPROGRESS
c                                                                      c      
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Solve system of equations using CGS iterative method            **
c **                                  ---                             **
c **********************************************************************
      TYPE(pslv) slv
      INTEGER  mpr,nwrk,pres
      REAL*8   wrk
      EXTERNAL PDINIT,PDVPROD,pimdcgs
      ALLOCATABLE wrk(:)
c
C***  SET INITIAL VALUES
c
c     if pret>0 use left preconditioning
      pres=0
      IF (pret.GT.0) pres=1
c
C***  ALLOCATE MEMORY FOR MATRIX AND VECTORS
c
      nwrk = (9)*peq
      ALLOCATE(wrk(nwrk),STAT=ierr)
      IF (ierr.NE.0) STOP 'PCGS: not enough memory'
      CALL pgetmyproc(env,mpr)
      CALL PDSETPAR(slv,neq,peq,peq,-1,-1,npr,mpr,pres,stopt,maxit,
     &                eps,env)
c
C***  CGS
c
	CALL pimdcgs(neq,peq,pnz,nwrk,iq,jq,dq,q,ia,ja,a,b,x,wrk,slv,
     &             matvec,preconl,preconr,PARDSUM,pdnrm,PPROGRESS)
c
C***  SET OUTPUT VALUES
c
      nit=slv%itno
      ierr=slv%status
c
C***  FREE MEMORY
c
      DEALLOCATE(wrk)
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pimdcgs
     &   (neq,peq,pnz,nwrk,iq,jq,dq,q,ia,ja,a,b,x,wrk,slv,
     &    matvec,preconl,preconr,PARDSUM,pdnrm,PPROGRESS)
c                                                                      c
C----------------------------------------------------------------------C
      USE inc_types
	TYPE(pslv) SLV
*     ..
*     .. Parameters ..
      REAL*8  ZERO
      PARAMETER (ZERO=0.0D0)
      REAL*8  ONE
      PARAMETER (ONE=1.0D0)
*     ..
*     .. Matrices A,Q ..
      INTEGER neq,peq,pnz,nwrk,iq(peq+1),jq(pnz),dq(peq),ia(peq+1),ja(pnz)
      REAL*8  q(pnz),a(pnz)
*     ..
*     .. Array Arguments ..
      REAL*8  b(peq),x(peq),wrk(nwrk)
*     ..
*     .. Function Arguments ..
      REAL*8   PDNRM
      EXTERNAL PDNRM
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL MATVEC,PARDSUM,PRECONL,PRECONR,PPROGRESS
*     ..
*     .. Local Scalars ..
      INTEGER BASIS,BLKSZ,CNVRTX,IP,IR,IRTILDE,IS,IT,ITNO,IU,IW,
     +        IXOLD,IZ,LDA,LOCLEN,MAXIT,N,NPROCS,PRECONT,PROCID,
     +        STATUS,STEPERR,STOPT
      REAL*8  ALPHA,BETA,EPS,EXITNORM,RHO,RHO0,RHSSTOP,XI
*     ..
*     .. Local Arrays ..
      REAL*8  DOTS(1)
*     ..
*     .. External Functions ..
      REAL*8   DDOT,PDSETRHSSTOP
      EXTERNAL DDOT,PDSETRHSSTOP
*     ..
*     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,PDGETPAR,PSTOPCRIT
*     ..
      CALL PDGETPAR(SLV,LDA,N,LOCLEN,BLKSZ,BASIS,NPROCS,
     +                PROCID,PRECONT,STOPT,MAXIT,ITNO,STATUS,
     +                STEPERR,EPS,EXITNORM)

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
      RHSSTOP = PDSETRHSSTOP(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,B,WRK(IR),EPS,SLV,
     +                      PRECONL,PDNRM)

*  1. r=Q1(b-AQ2x)
      IF (STOPT.NE.6) THEN
          IF (PRECONT.EQ.0) THEN
*     r=b-Ax
              CALL DCOPY(LOCLEN,B,1,WRK(IR),1)
              CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,X,WRK(IW),SLV)
              CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IR),1)

          ELSE IF (PRECONT.EQ.1) THEN
*     r=Q1(b-Ax)
              CALL DCOPY(LOCLEN,B,1,WRK(IZ),1)
              CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,X,WRK(IW),SLV)
              CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IZ),1)
              CALL PRECONL(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IR),SLV)

          ELSE IF (PRECONT.EQ.2) THEN
*     r=b-AQ2x
              CALL DCOPY(LOCLEN,B,1,WRK(IR),1)
              CALL PRECONR(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,X,WRK(IW),SLV)
              CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IW),WRK(IZ),SLV)
              CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IR),1)

          ELSE IF (PRECONT.EQ.3) THEN
*     r=Q1(b-AQ2x)
              CALL DCOPY(LOCLEN,B,1,WRK(IP),1)
              CALL PRECONR(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,X,WRK(IW),SLV)
              CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IW),WRK(IZ),SLV)
              CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IP),1)
              CALL PRECONL(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IP),WRK(IR),SLV)
          END IF

      ELSE
*     r has been set to Qb in the call to pdsetrhsstop
          IF (PRECONT.EQ.1) THEN
*     r=Q1(b-Ax)
              CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,X,WRK(IW),SLV)
              CALL PRECONL(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IZ),SLV)
              CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IR),1)

          ELSE IF (PRECONT.EQ.3) THEN
*     r=Q1(b-AQ2x)
              CALL PRECONR(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,X,WRK(IZ),SLV)
              CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IZ),WRK(IW),SLV)
              CALL PRECONL(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IZ),SLV)
              CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IR),1)
          END IF

      END IF

*  2. p=s=rtilde=r
      CALL DCOPY(LOCLEN,WRK(IR),1,WRK(IRTILDE),1)
      CALL DCOPY(LOCLEN,WRK(IR),1,WRK(IP),1)
      CALL DCOPY(LOCLEN,WRK(IR),1,WRK(IS),1)

*  3. rho=dot(rtilde,r)
      DOTS(1) = DDOT(LOCLEN,WRK(IRTILDE),1,WRK(IR),1)
      CALL PARDSUM(1,DOTS,SLV)
      RHO = DOTS(1)

*  Loop
      STATUS = 0
      EXITNORM = -ONE
      STEPERR = -1
      DO 20 ITNO = 1,MAXIT

*  4. w=Q1AQ2p
          IF (PRECONT.EQ.0) THEN
              CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IP),WRK(IW),SLV)

          ELSE IF (PRECONT.EQ.1) THEN
              CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IP),WRK(IZ),SLV)
              CALL PRECONL(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IW),SLV)

          ELSE IF (PRECONT.EQ.2) THEN
              CALL PRECONR(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IP),WRK(IZ),SLV)
              CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IZ),WRK(IW),SLV)

          ELSE IF (PRECONT.EQ.3) THEN
              CALL PRECONR(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IP),WRK(IW),SLV)
              CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IW),WRK(IZ),SLV)
              CALL PRECONL(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IW),SLV)
          END IF

*  5. xi=dot(rtilde,w)
          DOTS(1) = DDOT(LOCLEN,WRK(IRTILDE),1,WRK(IW),1)
          CALL PARDSUM(1,DOTS,SLV)
          XI = DOTS(1)

*  6. alpha=rho/xi
c         IF (XI.EQ.ZERO) THEN
          IF (ABS(XI).LE.ABS(EPSILON(RHO)*RHO)) THEN
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
              CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IW),WRK(IU),SLV)

          ELSE IF (PRECONT.EQ.1) THEN
              CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IW),WRK(IZ),SLV)
              CALL PRECONL(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IU),SLV)

          ELSE IF (PRECONT.EQ.2) THEN
              CALL PRECONR(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IZ),SLV)
              CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IZ),WRK(IU),SLV)

          ELSE IF (PRECONT.EQ.3) THEN
              CALL PRECONR(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IU),SLV)
              CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IU),WRK(IZ),SLV)
              CALL PRECONL(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IU),SLV)
          END IF

          CALL DAXPY(LOCLEN,-ALPHA,WRK(IU),1,WRK(IR),1)

* 11. check stopping criterion
          CALL PSTOPCRIT(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,IA,JA,A,B,
     +                  WRK(IR),WRK(IZ),X,WRK(IXOLD),WRK(IU),RHSSTOP,
     +                  CNVRTX,EXITNORM,STATUS,SLV,
     +                  MATVEC,MATVEC,PRECONR,PARDSUM,PDNRM)

*  Call monitoring routine
          CALL PPROGRESS(LOCLEN,ITNO,EXITNORM,X,WRK(IR),WRK(IZ),SLV)

          IF (STATUS.EQ.-5) THEN
              STEPERR = 11
              GO TO 9999
          ELSE IF (STATUS.EQ.0) THEN
              GO TO 9999
          END IF

* 12. rho=dot(rtilde0,r)
          RHO0 = RHO
          DOTS(1) = DDOT(LOCLEN,WRK(IRTILDE),1,WRK(IR),1)
          CALL PARDSUM(1,DOTS,SLV)
          RHO = DOTS(1)

* 13. beta=rho/rho0
c         IF (RHO0.EQ.ZERO) THEN
          IF (ABS(RHO0).LT.ABS(EPSILON(RHO)*RHO)) THEN
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

      IF (ITNO.GE.MAXIT) THEN
          STATUS = -1
          ITNO = MAXIT
      END IF

 9999 CONTINUE

      IF ((PRECONT.EQ.2) .OR. (PRECONT.EQ.3)) THEN
          CALL DCOPY(LOCLEN,X,1,WRK(IZ),1)
          CALL PRECONR(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IZ),X,SLV)
      END IF

*  Set output parameters
      SLV%ITNO      = ITNO
      SLV%STATUS    = STATUS
      SLV%STEPERR   = STEPERR
      SLV%EXITNORM  = EXITNORM

      RETURN
      END
C^L
C----------------------------------------------------------------------C
c----------------------------------------------------------------------c
c                                                                      c
c     PARALLEL BiCGSTAB(L) SOLVER PROGRAM BLOCK                        c
c                                                                      c
c----------------------------------------------------------------------c
C----------------------------------------------------------------------C
c $Date: 15/08/1998 $
c $Source: ~/bem/flow/libpsolve/prbicgstab.f $
c $Revision: 0.3 $
c $State: Exp $
c **********************************************************************
c **                                                                  **
c **  SUBROUTINES USED :                                              **
c **                                                                  **
c **                   - prbicgstab                                   **
c **                   - pimdrbicgstab                                **
c **                                                                  **
c **********************************************************************
C^L
C----------------------------------------------------------------------C
c                                                                      c      
      SUBROUTINE pRBiCGSTAB
     &   (env,neq,peq,pnz,npr,pret,maxit,stopt,eps,nit,ierr,
     &    iq,jq,dq,q,ia,ja,a,b,x,
     &    matvec,preconl,preconr,PARDSUM,pdnrm,PPROGRESS)

      USE inc_types
      TYPE(penv)  env
      INTEGER  neq,peq,pnz,npr,pret,maxit,stopt,nit,ierr,
     &         iq(peq+1),jq(pnz),dq(peq),ia(peq+1),ja(pnz),da(peq)
      REAL*8   q(pnz),a(pnz),b(peq),x(peq),eps,pdnrm
      EXTERNAL matvec,preconl,preconr,PARDSUM,pdnrm,PPROGRESS
c                                                                      c      
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Solve system of equations using RBi-CGSTAB iterative method     **
c **                                  --- -----                       **
c **********************************************************************
      TYPE(pslv) slv
      INTEGER  mpr,nwrk,basis,pres
      REAL*8   wrk
      EXTERNAL PDINIT,PDVPROD,pimdrbicgstab
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
      nwrk = (6+2*basis)*peq
      ALLOCATE(wrk(nwrk),STAT=ierr)
      IF (ierr.NE.0) STOP 'PRBiCGSTAB: not enough memory'
      CALL pgetmyproc(env,mpr)
      CALL PDSETPAR(slv,neq,peq,peq,-1,basis,npr,mpr,pres,stopt,maxit,
     &                eps,env)
c
c RBi-CGSTAB
c
      CALL pimdrbicgstab(neq,peq,pnz,nwrk,iq,jq,dq,q,ia,ja,a,b,x,wrk,slv,
     &                   matvec,preconl,preconr,PARDSUM,pdnrm,PPROGRESS)
c
c Set output values
c
      nit=slv%itno
      ierr=slv%status
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
      SUBROUTINE pimdrbicgstab
     &   (neq,peq,pnz,nwrk,iq,jq,dq,q,ia,ja,a,b,x,wrk,slv,
     &    matvec,preconl,preconr,PARDSUM,pdnrm,PPROGRESS)
c                                                                      c
C----------------------------------------------------------------------C
      USE inc_types
	TYPE(pslv) SLV
*
*     .. Parameters ..
      REAL*8  ZERO
      PARAMETER (ZERO=0.0D0)
      REAL*8  ONE
      PARAMETER (ONE=1.0D0)
      INTEGER IBDIM
      PARAMETER (IBDIM=8)
*     ..
*     .. Matrices Q,A ..
      INTEGER neq,peq,pnz,nwrk,iq(peq+1),jq(pnz),dq(peq),ia(peq+1),ja(pnz)
      REAL*8  q(pnz),a(pnz)
*     ..
*     .. Array Arguments ..
      REAL*8  b(peq),x(peq),wrk(nwrk)
*     ..
*     .. Function Arguments ..
      REAL*8   PDNRM
      EXTERNAL PDNRM
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL MATVEC,PARDSUM,PRECONL,PRECONR,PPROGRESS
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
      REAL*8   DDOT,PDSETRHSSTOP
      EXTERNAL DDOT,PDSETRHSSTOP
*     ..
*     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,PDINIT,PDGETPAR,PSTOPCRIT
*     ..
      CALL PDGETPAR(SLV,LDA,N,LOCLEN,BLKSZ,BASIS,NPROCS,
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
      RHSSTOP = PDSETRHSSTOP(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,B,WRK(IR),EPSILON,SLV,
     +                      PRECONL,PDNRM)

*  1. r=Q1(b-AQ2x)
      IF (PRECONT.EQ.0) THEN
*     r=b-Ax
          CALL DCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,X,WRK(IW),SLV)
          CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IR),1)

      ELSE IF (PRECONT.EQ.1) THEN
*     r=Q1(b-Ax)
          CALL DCOPY(LOCLEN,B,1,WRK(IZ),1)
          CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,X,WRK(IW),SLV)
          CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IZ),1)
          CALL PRECONL(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IR),SLV)

      ELSE IF (PRECONT.EQ.2) THEN
*     r=b-AQ2x
          CALL DCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL PRECONR(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,X,WRK(IW),SLV)
          CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IW),WRK(IZ),SLV)
          CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IR),1)

      ELSE IF (PRECONT.EQ.3) THEN
*     r=Q1(b-AQ2x)
          CALL DCOPY(LOCLEN,B,1,WRK(IW),1)
          CALL PRECONR(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,X,WRK(IR),SLV)
          CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IR),WRK(IZ),SLV)
          CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IW),1)
          CALL PRECONL(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IR),SLV)
      END IF

*  2. rtilde=r
      CALL DCOPY(LOCLEN,WRK(IR),1,WRK(IRTILDE),1)

*  3. u0=0
      CALL PDINIT(LOCLEN,ZERO,WRK(IU),1)

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
              CALL PARDSUM(1,DOTS,SLV)
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
                  CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IU+I1),WRK(IU+I2),SLV)

              ELSE IF (PRECONT.EQ.1) THEN
                  CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IU+I1),WRK(IW),SLV)
                  CALL PRECONL(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IU+I2),SLV)

              ELSE IF (PRECONT.EQ.2) THEN
                  CALL PRECONR(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IU+I1),WRK(IW),SLV)
                  CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IW),WRK(IU+I2),SLV)

              ELSE IF (PRECONT.EQ.3) THEN
                  CALL PRECONR(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IU+I1),WRK(IZ),SLV)
                  CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IZ),WRK(IW),SLV)
                  CALL PRECONL(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IU+I2),SLV)
              END IF

* 11. ksi=u(j+1)^{T}rtilde
              DOTS(1) = DDOT(LOCLEN,WRK(IU+I2),1,WRK(IRTILDE),1)
              CALL PARDSUM(1,DOTS,SLV)
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
                  CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IR+I1),WRK(IR+I2),SLV)

              ELSE IF (PRECONT.EQ.1) THEN
                  CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IR+I1),WRK(IW),SLV)
                  CALL PRECONL(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IR+I2),SLV)

              ELSE IF (PRECONT.EQ.2) THEN
                  CALL PRECONR(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IR+I1),WRK(IW),SLV)
                  CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IW),WRK(IR+I2),SLV)

              ELSE IF (PRECONT.EQ.3) THEN
                  CALL PRECONR(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IR+I1),WRK(IZ),SLV)
                  CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK(IZ),WRK(IW),SLV)
                  CALL PRECONL(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IR+I2),SLV)
              END IF

* 15. x0=x0+alpha*u0
              CALL DCOPY(LOCLEN,X,1,WRK(IXOLD),1)
              CALL DAXPY(LOCLEN,ALPHA,WRK(IU),1,X,1)
              I1 = I1 + LOCLEN
              I2 = I2 + LOCLEN
   30     CONTINUE

* 16. check stopping criterion
          CALL PSTOPCRIT(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,IA,JA,A,B,
     +                  WRK(IR),WRK(IZ),X,WRK(IXOLD),WRK(IW),
     +                  RHSSTOP,CNVRTX,EXITNORM,STATUS,SLV,
     +                  MATVEC,MATVEC,PRECONR,PARDSUM,PDNRM)

*  Call monitoring routine
          CALL PPROGRESS(LOCLEN,ITNO,EXITNORM,X,WRK(IR),WRK(IZ),SLV)

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
          CALL PARDSUM(2,DOTS,SLV)
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
              CALL PARDSUM(J-1,DOTS,SLV)
              I1 = LOCLEN
              DO 50 I = 1,J - 1
                  TAU(I,J) = DOTS(I)/SIGMA(I)
                  CALL DAXPY(LOCLEN,-TAU(I,J),WRK(IR+I1),1,WRK(IR+I0),1)
   50         CONTINUE

* 19. sigma(j)=r(j)^{T}r(j), gamma'(j)=r(0)^{T}r(j)/sigma(j)
              DOTS(1) = DDOT(LOCLEN,WRK(IR+I0),1,WRK(IR+I0),1)
              DOTS(2) = DDOT(LOCLEN,WRK(IR),1,WRK(IR+I0),1)
              CALL PARDSUM(2,DOTS,SLV)
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

      IF (ITNO.GE.MAXIT) THEN
          STATUS = -1
          ITNO = MAXIT
      END IF

 9999 CONTINUE

      IF ((PRECONT.EQ.2) .OR. (PRECONT.EQ.3)) THEN
          CALL DCOPY(LOCLEN,X,1,WRK(IZ),1)
          CALL PRECONR(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,WRK(IZ),X,SLV)
      END IF

*  Set output parameters
      SLV%ITNO      = ITNO
      SLV%STATUS    = STATUS
      SLV%STEPERR   = STEPERR
      SLV%EXITNORM  = EXITNORM

      RETURN
      END
C^L
C----------------------------------------------------------------------C
c----------------------------------------------------------------------c
c                                                                      c
c     READ PARALLEL DATA PROGRAM BLOCK                                 c
c                                                                      c
c----------------------------------------------------------------------c
C----------------------------------------------------------------------C
c $Date: 15/08/1998 $
c $Source: ~/bem/flow/libpsolve/pread.f $
c $Revision: 0.3 $
c $State: Exp $
c **********************************************************************
c **                                                                  **
c **  SUBROUTINES USED :                                              **
c **                                                                  **
c **                   - preadcrs                                     **
c **                   - preadrhs                                     **
c **                   - pureadcrs                                    **
c **                   - pureadrhs                                    **
c **                                                                  **
c **********************************************************************     
C^L
C-----------------------------------------------------------------------C
c                                                                       c
      SUBROUTINE preadcrs (peq,pnz,npr,env,ia,ja,a,lun,ierr)

      USE inc_types
      TYPE(penv)  env
      INTEGER peq,pnz,npr,ia(peq+1),ja(pnz),lun,ierr
      REAL*8  a(pnz)
c                                                                       c
C-----------------------------------------------------------------------C
      INTEGER  mpr,isb,ise,jsb,jse,i,neq,nnz,foo
      REAL*8   dfoo
      INTEGER, POINTER :: nod(:),nnod(:),mat(:),nmat(:)
c
C***  GET PARALLEL ENVIRONMENT
c
      mpr  =  env%myproc
      nod  => env%inod
      nnod => env%nnod
      mat  => env%imat
      nmat => env%nmat
c
      IF (mpr.LT.npr) THEN
        isb = nod(mpr)
        ise = nod(npr)+nnod(npr)-nod(mpr+1)
        jsb = mat(mpr)
        jse = mat(npr)+nmat(npr)-mat(mpr+1)
      ELSE
        isb = nod(mpr)
        ise = 0
        jsb = mat(mpr)
        jse = 0
      END IF
c     
C***  READ THE SIZE OF THE MATRIX
c
      REWIND(lun)
      READ (lun, *, err=1000, end=1010) neq
      IF (neq.LE.0) GOTO 1020
c     
C***  READ THE POINTER ARRAY iA(*)
c     
      READ (lun, *, err=1001, end=1010) 
     &  (foo,i=1,isb),(ia(i),i=1,peq+1),(foo,i=1,ise)
c
      DO i=1,peq+1
        ia(i)=ia(i)-mat(mpr)
      END DO
c     
C***  NUMBER OF NONE-ZEROS
c     
      nnz = ia(peq+1) - 1
      IF (nnz.LE.0) GOTO 1030
c
C***  READ THE COLUMN INDICES ARRAY jA(*)
c
      READ (lun, *, err=1002, end=1010) 
     &  (foo,i=1,jsb),(ja(i),i=1,pnz),(foo,i=1,jse)
c
C***  READ THE MATRIX ELEMENTS A(*)
c
      READ (lun, *, err=1003, end=1010) (dfoo,i=1,jsb),(a(i),i=1,pnz)
c     
C***  NORMAL RETURN
c
      ierr=0
      RETURN
c     
C***  ERROR HANDLING CODE
c     
c     Error in reading I/O unit
 1000 ierr=1
      GOTO 2000
 1001 ierr=2
      GOTO 2000
 1002 ierr=3
      GOTO 2000
 1003 ierr=4
      GOTO 2000
c
c     EOF reached in reading
 1010 ierr=5
      GOTO 2000
c
c     NEQ non-positive or too large
 1020 ierr=6
      GOTO 2000
c
c     NNZ non-positive or too large
 1030 ierr=7
c     
c     The real return statement
c     
 2000 RETURN
      END SUBROUTINE
C^L
C-----------------------------------------------------------------------C
c                                                                       c
      SUBROUTINE preadrhs (peq,npr,env,b,lun,ierr)

      USE inc_types
      TYPE(penv)  env
      INTEGER peq,npr,lun,ierr
      REAL*8  b(peq)
c                                                                       c
C-----------------------------------------------------------------------C
      INTEGER  mpr,sb,i,neq
      REAL*8   dfoo
      INTEGER, POINTER :: nod(:)
c
C***  GET PARALLEL ENVIRONMENT
c
      mpr =  env%myproc
      nod => env%inod
c
      sb = nod(mpr)
c
C***  READ THE SIZE OF THE MATRIX
c
      REWIND(lun)
      READ (lun, *, err=1000, end=1010) neq
      IF (neq.LE.0) GOTO 1020
c
C***  READ THE VECTOR ELEMENTS A(*)
c
      READ (lun, *, err=1000, end=1010) (dfoo,i=1,sb),(b(i),i=1,peq)
c
C***  NORMAL RETURN
c
      ierr=0
      RETURN
c
C***  ERROR HANDLING CODE
c
c     Error in reading I/O unit
 1000 ierr=1
      GOTO 2000
c
c     EOF reached in reading
 1010 ierr=2
      GOTO 2000
c
c     NEQ non-positive or too large
 1020 ierr=3
      GOTO 2000
c
c     The real return statement
c
 2000 RETURN
      END SUBROUTINE
C^L
C-----------------------------------------------------------------------C
c                                                                       c
      SUBROUTINE pureadcrs (peq,pnz,npr,env,ia,ja,a,lun,ierr)

      USE inc_types
      TYPE(penv)  env
      INTEGER peq,pnz,npr,ia(peq+1),ja(pnz),lun,ierr
      REAL*8  a(pnz)
c                                                                       c
C-----------------------------------------------------------------------C
      INTEGER  mpr,isb,ise,i,j,neq,nnz
      INTEGER, POINTER :: nod(:),nnod(:),mat(:),nmat(:)
c
C***  GET PARALLEL ENVIRONMENT
c
      mpr  =  env%myproc
      nod  => env%inod
      nnod => env%nnod
      mat  => env%imat
      nmat => env%nmat
c
      IF (mpr.LT.npr) THEN
        isb = nod(mpr)
        ise = nod(npr)+nnod(npr)-nod(mpr+1)
      ELSE
        isb = nod(mpr)
        ise = 0
      END IF
c     
C***  READ THE SIZE OF THE MATRIX
c
      REWIND(lun)
      READ (lun, err=1000, end=1010) neq,nnz
      IF (neq.LE.0) GOTO 1020
c     
C***  READ THE POINTER ARRAY iA(*)
c     
      DO i=1,isb
        READ (lun, err=1000, end=1010) 
      END DO
      DO i=1,peq+1
        READ (lun, err=1000, end=1010) ia(i)
      END DO
      DO i=1,ise
        READ (lun, err=1000, end=1010) 
      END DO
c
      DO i=1,peq+1
        ia(i)=ia(i)-mat(mpr)
      END DO
c     
C***  NUMBER OF NONE-ZEROS
c     
      nnz = ia(peq+1) - 1
      IF (nnz.LE.0) GOTO 1030
c
C***  READ THE COLUMN INDICES ARRAY jA(*)
c
      DO i=1,isb
        READ (lun, err=1000, end=1010) 
      END DO
      DO i=1,peq
        READ (lun, err=1000, end=1010) (ja(j), j=ia(i),ia(i+1)-1)
      END DO
      DO i=1,ise
        READ (lun, err=1000, end=1010) 
      END DO
c
C***  READ THE MATRIX ELEMENTS A(*)
c
      DO i=1,isb
        READ (lun, err=1000, end=1010) 
      END DO
      DO i=1,peq
        READ (lun, err=1000, end=1010) (a(j), j=ia(i),ia(i+1)-1)
      END DO
c     
C***  NORMAL RETURN
c
      ierr=0
      RETURN
c     
C***  ERROR HANDLING CODE
c     
c     Error in reading I/O unit
 1000 ierr=1
      GOTO 2000
c
c     EOF reached in reading
 1010 ierr=2
      GOTO 2000
c
c     NEQ non-positive or too large
 1020 ierr=3
      GOTO 2000
c
c     NNZ non-positive or too large
 1030 ierr=4
c     
c     The real return statement
c     
 2000 RETURN
      END SUBROUTINE
C^L
C-----------------------------------------------------------------------C
c                                                                       c
      SUBROUTINE pureadrhs (peq,npr,env,b,lun,ierr)

      USE inc_types
      TYPE(penv)  env
      INTEGER peq,npr,lun,ierr
      REAL*8  b(peq)
c                                                                       c
C-----------------------------------------------------------------------C
      INTEGER  mpr,isb,i,neq
      INTEGER, POINTER :: nod(:)
c
C***  GET PARALLEL ENVIRONMENT
c
      mpr =  env%myproc
      nod => env%inod
c
      isb = nod(mpr)
c
C***  READ THE SIZE OF THE MATRIX
c
      REWIND(lun)
      READ (lun, err=1000, end=1010) neq
      IF (neq.LE.0) GOTO 1020
c
C***  READ THE VECTOR ELEMENTS A(*)
c
      DO i=1,isb
        READ (lun, err=1000, end=1010) 
      END DO
      DO i=1,peq
        READ (lun, err=1000, end=1010) b(i)
      END DO
c
C***  NORMAL RETURN
c
      ierr=0
      RETURN
c
C***  ERROR HANDLING CODE
c
c     Error in reading I/O unit
 1000 ierr=1
      GOTO 2000
c
c     EOF reached in reading
 1010 ierr=2
      GOTO 2000
c
c     NEQ non-positive or too large
 1020 ierr=3
      GOTO 2000
c
c     The real return statement
c
 2000 RETURN
      END SUBROUTINE
C^L
C----------------------------------------------------------------------C
c----------------------------------------------------------------------c
c                                                                      c
c     BASIC PARALLEL ROUTINES PROGRAM BLOCK                            c
c                                                                      c
c----------------------------------------------------------------------c
C----------------------------------------------------------------------C
c $Date: 15/08/1998 $
c $Source: ~/bem/flow/libpsolve/plib.f $
c $Revision: 0.3 $
c $State: Exp $
c **********************************************************************
c **                                                                  **
c **  SUBROUTINES USED :                                              **
c **                                                                  **
c **                   - pmatvec                                      **
c **                   - pmatvecfm                                    **
c **                   - pprenon                                      **
c **                   - ppredia                                      **
c **                   - pprediafm                                    **
c **                   - ppreilu                                      **
c **                   - PARDSUM                                        **
c **                   - PARDNRM2                                       **
c **                   - PPROGRESS                                     **
c **                   - PPRINTV                                       **
c **                   - PDINIT                                        **
c **                   - PDVPROD                                       **
c **                   - pdsetrhsstop                                 **
c **                   - pstopcrit                                    **
c **                   - PDSETPAR                                   **
c **                   - PDGETPAR                                   **
c **                                                                  **
c **********************************************************************
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pmatvec2(neq,peq,pnz,ia,ja,a,u,v,slv)

      USE inc_types
	TYPE(pslv) slv
      INTEGER neq,peq,pnz,ia(peq+1),ja(pnz)
      REAL*8  a(pnz),u(peq),v(peq)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Calculate Matrix Vector product v=Au                            **
c **            ---    ---                                            **
c **********************************************************************
      INTEGER i,j,ii,inp,nnp
      REAL*8  sum,w(neq)
c
      CALL pgetnodes(slv%env,inp,nnp)
c
      CALL dcopy(slv%loclen,u,1,w(inp+1),1)
      CALL pgetsvec3(slv%lda,slv%n,3,slv%env,u,w)
c
      DO i = inp+1,inp+nnp
        ii=i-inp
        sum=0.d0
        DO j = ia(ii),ia(ii+1) - 1
          sum = sum + a(j)*w(ja(j))
        END DO
        v(ii)=sum
      END DO
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pmatvec (neq,peq,pnz,ia,ja,a,u,v,slv)

      USE inc_types
	TYPE(pslv) slv
      INTEGER neq,peq,pnz,ia(peq+1),ja(pnz)
      REAL*8  a(pnz),u(peq),v(peq)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Calculate Matrix Vector product v=Au                            **
c **            ---    ---                                            **
c **********************************************************************
      INTEGER i,j,inp,nnp,imtr,nmtr
      REAL*8  sum,w(neq)
c
      CALL pgetnodes(slv%env,inp,nnp)
      CALL pgetmatrix(slv%env,imtr,nmtr)
c
      CALL pgetvec(neq,peq,slv%nprocs,slv%env,u,w)
c
      DO i = 1,nnp
        sum=0.0D00
        DO j = ia(i),ia(i+1)-1
          sum = sum + a(j)*w(ja(j))
        END DO
        v(i)=sum
      END DO
c
      RETURN
      END
C^L           
c ---------------------------------------------------------------------
      SUBROUTINE pmatsvec (neq,peq,pnz,ia,ja,a,u,v,slv)

      USE inc_types
	TYPE(pslv) slv
      INTEGER neq,peq,pnz,ia(peq+1),ja(pnz)
      REAL*8  a(pnz),u(peq),v(peq)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Calculate Matrix Vector product v=Au                            **
c **            ---    ---                                            **
c **********************************************************************
      INTEGER i,j,inp,nnp,imtr,nmtr
      REAL*8  sum,w(neq)
c
      CALL pgetnodes(slv%env,inp,nnp)
      CALL pgetmatrix(slv%env,imtr,nmtr)
c
      CALL pgetsvec(neq,peq,slv%env,u,w)
c

      DO i = 1,nnp
        sum=0.0D00
        DO j = ia(i),ia(i+1)-1
          sum = sum + a(j)*w(ja(j))
        END DO
        v(i)=sum
      END DO
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pmatvecfm (neq,peq,pnz,ia,ja,a,u,v,slv)

      USE inc_types
	TYPE(pslv) slv
      INTEGER neq,peq,pnz,ia(peq+1),ja(pnz)
      REAL*8  a(peq,neq),u(peq),v(peq)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Calculate Matrix Vector product v=Au                            **
c **            ---    ---                                            **
c **********************************************************************
      INTEGER i,j
      REAL*8  w(neq)
c
      CALL pgetvec(neq,peq,slv%nprocs,slv%env,u,w)
c
      v = 0.0D0
      DO j = 1,neq
        DO i = 1,peq
          v(i) = v(i) + a(i,j)*w(j)
        END DO
      END DO
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pprenon (neq,peq,pnz,iq,jq,dq,q,u,v,SLV)

      INTEGER neq,peq,pnz,iq(peq+1),jq(pnz),dq(peq),SLV(13)
      REAL*8  q(peq),u(peq),v(peq)
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
      SUBROUTINE ppredia (neq,peq,pnz,iq,jq,dq,q,u,v,slv)

      USE inc_types
	TYPE(pslv) slv
      INTEGER neq,peq,pnz,iq(peq+1),jq(pnz),dq(peq)
      REAL*8  q(pnz),u(peq),v(peq)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Precondition system with diagonal preconditioner                **
c **  ---                      ---                                    **
c **********************************************************************
      INTEGER i,inp,nnp
c
      CALL pgetnodes(slv%env,inp,nnp)
c
      DO i=1,nnp
        v(i)=q(dq(i))*u(i)
      END DO
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pprediafm (neq,peq,pnz,iq,jq,dq,q,u,v,slv)

      USE inc_types
	TYPE(pslv) slv
      INTEGER neq,peq,pnz,iq(peq+1),jq(pnz),dq(peq)
      REAL*8  q(peq),u(peq),v(peq)
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
      SUBROUTINE ppreilu (neq,peq,pnz,iq,jq,dq,q,b,x,slv)

      USE inc_types
	TYPE(pslv) slv
      INTEGER neq,peq,pnz,iq(peq+1),jq(pnz),dq(peq)
      REAL*8  q(pnz),b(peq),x(peq)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Precondition system with ILU preconditioner                     **
c **  ---                      ---                                    **
c **********************************************************************
      INTEGER i,j,ij,inp,nnp
      REAL*8  sum
c
C***  GET ENVIRONMENT VARIABLES
c
      CALL pgetnodes(slv%env,inp,nnp)
c
C***  Lx'=b; RHS FACTORIZATION
c
      x(1)=b(1)
      DO i=2,nnp
        sum=b(i)
        DO ij=iq(i),dq(i)-1
          j=jq(ij)-inp
          IF (j.GT.0) THEN
            sum=sum-q(ij)*x(j)
          END IF
        END DO
        x(i)=sum
      END DO
c
C***  Ux=x'; BACK SUBSTITUTION
c
      x(nnp)=x(nnp)/q(dq(nnp))
      DO i=nnp-1,1,-1
        sum=x(i)
        DO ij=dq(i)+1,iq(i+1)-1
          j=jq(ij)-inp
          IF (j.LE.nnp) THEN
            sum=sum-q(ij)*x(j)
          END IF
        END DO
        x(i)=sum/q(dq(i))
      END DO
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE PARDSUM (isize,x,slv)

      USE inc_types
	TYPE(pslv) slv
      INTEGER isize
      REAL*8  x(isize)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Parallel Double precision SUM of a vector                       **
c **  -        -                ---                                   **
c **********************************************************************
c      INCLUDE 'mpif.h'
      INTEGER comm,ierr
      REAL*8  wrk(isize)
c
C***  GET ENVIRONMENT VARIABLE
c
      CALL pgetcomm(slv%env,comm)
c
C***  GLOBAL SUM
c
      CALL MPI_ALLREDUCE(x(1),wrk(1),isize,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                   comm,ierr)
c
      CALL DCOPY(isize,wrk,1,x,1)
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      REAL*8  FUNCTION PARDNRM2 (loclen,u,slv)

      USE inc_types
	TYPE(pslv) slv
      INTEGER loclen
      REAL*8  u(loclen)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Parallel Double precision NORM of a vector                      **
c **  -        -                ----                                  **
c **********************************************************************
c      INCLUDE 'mpif.h'
      INTEGER comm,ierr
      REAL*8  sum,wrk(1),ddot
      EXTERNAL ddot
c
C***  GET ENVIRONMENT VARIABLE
c
      CALL pgetcomm(slv%env,comm)
c
C***  LOCAL DOT PRODUCT
c
      sum = ddot(loclen,u,1,u,1)
c
C***  GLOBAL SUM
c
      CALL MPI_ALLREDUCE(sum,wrk,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                   comm,ierr)
c
C***  NORM = SQRT(GLOBAL SUM)
c
      PARDNRM2=SQRT(wrk(1))
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE PPROGRESS (loclen,itno,normres,x,res,trueres,slv)

      USE inc_types
	TYPE(pslv) slv
      INTEGER loclen,itno
      REAL*8  normres,x(loclen),res(loclen),trueres(loclen)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Write information about PPROGRESS of iterative method            **
c **                          --------                                **
c **********************************************************************
*     INTEGER mpr
c
C***  GET ENVIRONMENT VARIABLE
c
*     CALL pgetmyproc(SLV(4),mpr)
c
C***  WRITE PPROGRESS INFORMATION
c
*     IF (mpr.EQ.1) THEN
*       WRITE (6,FMT=9000) ITNO,NORMRES
*       WRITE (6,FMT=9010) 'X:'
*       CALL PPRINTV(LOCLEN,X)
*       WRITE (6,FMT=9010) 'RES:'
*       CALL PPRINTV(LOCLEN,RES)
*       WRITE (6,FMT=9010) 'TRUE RES:'
*       CALL PPRINTV(LOCLEN,TRUERES)
*     END IF
 9000 FORMAT (I5,1X,E16.10)
 9010 FORMAT (A)
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE PPRINTV (n,u)

      INTEGER n
      REAL*8  u(n)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  PRINT Vector                                                    **
c **  ----- -                                                         **
c **********************************************************************
      INTEGER i
c
      WRITE (6,FMT=9010) (u(i),i=1,n)
 9000 FORMAT (5(D14.8,1X))
 9010 FORMAT (8(D9.3))
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE PDINIT (n,alpha,dx,incx)
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
*        code for upequal increments or equal increments
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
      SUBROUTINE PDVPROD (n,dx,incx,dy,incy)
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
*        code for upequal increments or equal increments
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
      REAL*8 FUNCTION PDSETRHSSTOP
     &         (neq,peq,pnz,iq,jq,dq,q,b,r,eps,slv,preconl,pdnrm)

      USE inc_types
	TYPE(pslv) slv
      INTEGER neq,peq,pnz,iq(peq+1),jq(pnz),dq(peq)
      REAL*8  q(peq),b(peq),r(peq),eps,pdnrm
      EXTERNAL preconl,pdnrm
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Double precision SET RHS STOPing criterion                      **
c **  -                --- --- ----                                   **
c **********************************************************************
      INTEGER loclen,stopt
c
      loclen = slv%loclen
      stopt  = slv%stopt
c
      IF ((stopt.EQ.1) .OR. (stopt.EQ.4) .OR. (stopt.EQ.7)) THEN
c
C***  ||r||<eps or ||Q1r||<eps ||x(k)-x(k-1)||<eps
c
        pdsetrhsstop = eps

      ELSE IF ((stopt.EQ.2) .OR. (stopt.EQ.3) .OR. (stopt.EQ.5)) THEN
c
C***  ||r||<eps||b|| or sqrt(r(Q1r))<eps||b|| or ||Q1r||<eps||b||
c
        pdsetrhsstop = eps*pdnrm(loclen,b,slv)

      ELSE IF (stopt.EQ.6) THEN
c
C***  ||Q1r||<eps||Q1b||
c
        CALL preconl(neq,peq,pnz,iq,jq,dq,q,b,r,slv)
        pdsetrhsstop = eps*pdnrm(loclen,r,slv)
      END IF
c
      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pstopcrit
     &   (neq,peq,pnz,iq,jq,dq,q,ia,ja,a,b,r,rtrue,x,xold,wrk,
     &    rhsstop,cnvrtx,exitnorm,status,slv,
     &    matvec,tmatvec,preconr,PARDSUM,pdnrm)

      USE inc_types
	TYPE(pslv) slv
      INTEGER neq,peq,pnz,iq(peq+1),jq(pnz),dq(peq),ia(peq+1),ja(pnz),
     &        cnvrtx,status
      REAL*8  q(pnz),a(pnz),b(peq),r(peq),rtrue(peq),x(peq),xold(peq),
     &        wrk(peq),rhsstop,exitnorm,pdnrm
      EXTERNAL matvec,tmatvec,preconr,PARDSUM,pdnrm
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Check STOPing CRITerion                                         **
c **        ----    ----                                              **
c **********************************************************************
      INTEGER  loclen,precont,stopt
      REAL*8   dots(1),ddot
      EXTERNAL daxpy,dcopy,ddot
      INTRINSIC sqrt
c
      REAL*8   ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
c
      LOCLEN  = SLV%LOCLEN
      PRECONT = SLV%PRECONT
      STOPT   = SLV%STOPT
c
      IF ((STOPT.EQ.1) .OR. (STOPT.EQ.2) .OR. (STOPT.EQ.3)) THEN
c
C***  Compute true residual if needed
c
          CALL DCOPY(LOCLEN,B,1,RTRUE,1)
          IF ((PRECONT.EQ.2) .OR. (PRECONT.EQ.3)) THEN
              CALL PRECONR(NEQ,PEQ,PNZ,IQ,JQ,DQ,Q,X,WRK,SLV)
              IF (CNVRTX.EQ.1) THEN
C***  r=b-AATQ2x
                  CALL TMATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK,XOLD,SLV)
                  CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,XOLD,WRK,SLV)
                  CALL DAXPY(LOCLEN,-ONE,WRK,1,RTRUE,1)
              ELSE
C***  r=b-AQ2x
                  CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,WRK,XOLD,SLV)
                  CALL DAXPY(LOCLEN,-ONE,XOLD,1,RTRUE,1)
              END IF
          ELSE IF (CNVRTX.EQ.1) THEN
C***  r=b-AATx
              CALL TMATVEC(NEQ,PEQ,PNZ,IA,JA,A,X,XOLD,SLV)
              CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,XOLD,WRK,SLV)
              CALL DAXPY(LOCLEN,-ONE,WRK,1,RTRUE,1)
          ELSE
C***  r=b-Ax
              CALL MATVEC(NEQ,PEQ,PNZ,IA,JA,A,X,WRK,SLV)
              CALL DAXPY(LOCLEN,-ONE,WRK,1,RTRUE,1)
          END IF
      END IF

      IF ((STOPT.EQ.1) .OR. (STOPT.EQ.2)) THEN
c
C***  ||r||<epsilon or ||r||<epsilon||b||
c
          EXITNORM = PDNRM(LOCLEN,RTRUE,SLV)
          IF (EXITNORM.LT.RHSSTOP) THEN
              STATUS = 0
          ELSE
              STATUS = -99
          END IF

      ELSE IF (STOPT.EQ.3) THEN
c
C***  sqrt(rT(Q1r))<epsilon||b||
c
          DOTS(1) = DDOT(LOCLEN,RTRUE,1,R,1)
          CALL PARDSUM(1,DOTS(1),SLV)
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
c
C***  ||Q1r||<epsilon or ||Q1r||<epsilon||b|| or ||Q1r||<epsilon||Q1b||
c
          EXITNORM = PDNRM(LOCLEN,R,SLV)
          IF (EXITNORM.LT.RHSSTOP) THEN
              STATUS = 0
          ELSE
              STATUS = -99
          END IF

      ELSE IF (STOPT.EQ.7) THEN
c
C***  ||x-x0||<epsilon
c
          CALL DCOPY(LOCLEN,X,1,WRK,1)
          CALL DAXPY(LOCLEN,-ONE,XOLD,1,WRK,1)
          EXITNORM = PDNRM(LOCLEN,WRK,SLV)
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
      SUBROUTINE PDSETPAR (slv,lda,n,loclen,blksz,basis,
     &    nprocs,procid,precont,stopt,maxit,epsilon,env)

      USE inc_types
	TYPE(pslv) slv
	TYPE(penv) env
      INTEGER lda,n,loclen,blksz,basis,nprocs,procid,precont,stopt,maxit
      REAL(8) epsilon
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  SET PARameters for parallel calculation                         **
c **  --- ---                                                         **
c **********************************************************************
      REAL*8  ONE
      PARAMETER (ONE=1.0D0)

      SLV%LDA     = LDA
      SLV%N       = N
      SLV%LOCLEN  = LOCLEN
      SLV%BLKSZ   = BLKSZ
      SLV%BASIS   = BASIS
      SLV%NPROCS  = NPROCS
      SLV%PROCID  = PROCID
      SLV%PRECONT = PRECONT
      SLV%STOPT   = STOPT
      SLV%MAXIT   = MAXIT
      SLV%ITNO    = -1
      SLV%STATUS  = -1
      SLV%STEPERR = -1

      SLV%EPSILON = EPSILON
      SLV%EXITNORM = -ONE

	SLV%ENV     = ENV

      RETURN
      END
C^L
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE PDGETPAR (slv,lda,n,loclen,blksz,basis,
     &    nprocs,procid,precont,stopt,maxit,itno,status,steperr,
     &    epsilon,exitnorm)

      USE inc_types
	TYPE(pslv) slv
      INTEGER lda,n,loclen,blksz,basis,nprocs,procid,
     &        precont,stopt,maxit,itno,status,steperr
      REAL*8  epsilon,exitnorm
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  GET PARameters for parallel calculation                         **
c **  --- ---                                                         **
c **********************************************************************

      LDA     = SLV%LDA
      N       = SLV%N
      LOCLEN  = SLV%LOCLEN
      BLKSZ   = SLV%BLKSZ
      BASIS   = SLV%BASIS
      NPROCS  = SLV%NPROCS
      PROCID  = SLV%PROCID
      PRECONT = SLV%PRECONT
      STOPT   = SLV%STOPT
      MAXIT   = SLV%MAXIT
      ITNO    = SLV%ITNO
      STATUS  = SLV%STATUS
      STEPERR = SLV%STEPERR

      EPSILON = SLV%EPSILON
      EXITNORM = SLV%EXITNORM

      RETURN
      END
