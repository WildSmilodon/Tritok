C______________________________________________________________________C
      SUBROUTINE CommErr(env,ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)
C     
C     Communicate errors
C______________________________________________________________________C
      USE inc_types        

      REAL(8) ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz ! ocena natancnosti integralov
      REAL(8) ontHtxDz,ontHtyDx,ontHtzDy

      REAL(8) myontH,myontHtx,myontHty,myontHtz,myontDx,myontDy,myontDz ! ocena natancnosti integralov
      REAL(8) myontHtxDz,myontHtyDx,myontHtzDy

      
      INTEGER ierr
     
      TYPE(penv) :: env      

      myontH=ontH
      myontHtx=ontHtx
      myontHty=ontHty
      myontHtz=ontHtz
      myontDx=ontDx
      myontDy=ontDy
      myontDz=ontDz
      myontHtxDz=ontHtxDz
      myontHtyDx=ontHtyDx
      myontHtzDy=ontHtzDy

c             MPI_REDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, ROOT, COMM, IERROR) 
c          MPI_ALLREDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)


      CALL MPI_ALLREDUCE(myontH, ontH, 1, MPI_DOUBLE_PRECISION, MPI_MAX, env%comm, ierr)
      CALL MPI_ALLREDUCE(myontHtx, ontHtx, 1, MPI_DOUBLE_PRECISION, MPI_MAX, env%comm, ierr)
      CALL MPI_ALLREDUCE(myontHty, ontHty, 1, MPI_DOUBLE_PRECISION, MPI_MAX, env%comm, ierr)
      CALL MPI_ALLREDUCE(myontHtz, ontHtz, 1, MPI_DOUBLE_PRECISION, MPI_MAX, env%comm, ierr)
      CALL MPI_ALLREDUCE(myontDx, ontDx, 1, MPI_DOUBLE_PRECISION, MPI_MAX, env%comm, ierr)
      CALL MPI_ALLREDUCE(myontDy, ontDy, 1, MPI_DOUBLE_PRECISION, MPI_MAX, env%comm, ierr)
      CALL MPI_ALLREDUCE(myontDz, ontDz, 1, MPI_DOUBLE_PRECISION, MPI_MAX, env%comm, ierr)
      CALL MPI_ALLREDUCE(myontHtxDz, ontHtxDz, 1, MPI_DOUBLE_PRECISION, MPI_MAX, env%comm, ierr)
      CALL MPI_ALLREDUCE(myontHtyDx, ontHtyDx, 1, MPI_DOUBLE_PRECISION, MPI_MAX, env%comm, ierr)
      CALL MPI_ALLREDUCE(myontHtzDy, ontHtzDy, 1, MPI_DOUBLE_PRECISION, MPI_MAX, env%comm, ierr)

      END


C______________________________________________________________________C
      SUBROUTINE EstRAM(env,mesh,io,inp)
C
C     $: Estimate requred RAM per processor
C
C______________________________________________________________________C
      USE inc_types   
      
      TYPE(meshType) :: mesh
      TYPE(penv) :: env
      TYPE(IOtype)    :: io
      TYPE(inputtype) :: inp
      
      INTEGER i
      REAL(8) dnb,dnn,nicnsp,npoc,nicnpoc,totmem,dnl
            
      dnb=DBLE(mesh%nbnodes)
      dnn=DBLE(mesh%nnodes)
      nicnsp=DBLE(mesh%nicnsp)
      npoc=DBLE(mesh%npoc)
      nicnpoc=DBLE(mesh%nicnpoc)
      dnl=DBLE(mesh%nnodes8)
      

      IF (inp%iBw.EQ.1.OR.inp%iBw.EQ.3) THEN
        totmem=(7.0D0*dnb*dnb+3.0D0*dnb*dnn)*8.0D0/1024.0D0/1024.0D0
      ELSE IF (inp%iBw.EQ.2) THEN
        totmem=(10.0D0*dnb*dnb+3.0D0*dnb*dnl)*8.0D0/1024.0D0/1024.0D0
      ELSE
        totmem=-9.99D00
      END IF
      WRITE (io%l,'(A,F12.4)') "Total amount of memory for kinematics matrices [Mb] =",totmem
      WRITE (io%l,'(A)') "Divided across processors: (no., nnod, mem [Mb])"
      DO i=1,env%nproc
        WRITE (io%l,'(I5,I7,F12.4)') env%myproc,env%nnod(i),totmem/dnb*DBLE(env%nnod(i))
      END DO

      IF (inp%iBw.EQ.1.OR.inp%iBw.EQ.3) THEN
        totmem=(7.0D0*dnb*dnb+3.0D0*dnb*dnn+7.0D0*nicnsp*npoc+
     &          3.0D0*nicnpoc*npoc)*8.0D0/1024.0D0/1024.0D0
      ELSE IF (inp%iBw.EQ.2) THEN
        totmem=(10.0D0*dnb*dnb+3.0D0*dnb*dnl+7.0D0*nicnsp*npoc+
     &          3.0D0*nicnpoc*npoc)*8.0D0/1024.0D0/1024.0D0
      ELSE
        totmem=-9.99D00
      END IF

      WRITE (io%l,'(A,F12.4/)') "Total amount of memory for matrices [Mb] =",totmem
        
      END 


C______________________________________________________________________C
      SUBROUTINE SetEnvGlob(env,mesh)
C     
C     Set Environment Global size
C______________________________________________________________________C
      USE inc_types        

      TYPE(meshType) :: mesh      
      TYPE(penv) :: env      
C
C
      env%nbnodes=mesh%nbnodes
      env%nnodes=mesh%nnodes
      env%nbelem=mesh%nbelem
      env%nicell=mesh%nicell
      env%npob=mesh%npob
      env%npoc=mesh%npoc
      env%npx=mesh%npx
c
      RETURN
      END

C     __________________________________________________________________
C    
      SUBROUTINE DivideNodes(env,mesh)
C     __________________________________________________________________
      USE inc_types        

      TYPE(meshType) :: mesh      
      TYPE(penv) :: env
      
      INTEGER npp  
      INTEGER npr,mpr,i

      npr=env%nproc
      mpr=env%mpr

      ALLOCATE(env%inod(npr),env%nnod(npr))


      npp=INT(REAL(mesh%nbnodes)/REAL(env%nproc)+0.5D00)

      DO i=1,npr-1
        env%inod(i)=(i-1)*npp
        env%nnod(i)=npp
      END DO

      env%inod(npr)=(i-1)*npp
      env%nnod(npr)=mesh%nbnodes-(i-1)*npp

      env%zac=env%inod(mpr+1)+1
      env%kon=env%inod(mpr+1)+env%nnod(mpr+1)
      env%nmn=env%nnod(mpr+1)

      RETURN
      END

C -----------------------------------------------------------------------------
      SUBROUTINE par_init(env)
C
C     $: Init parallel environment
C
C -----------------------------------------------------------------------------
      USE inc_types 
      INTEGER ierr
      TYPE(penv) :: env
      
      ierr=0
      
      env%comm = MPI_COMM_WORLD
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(env%comm,env%nproc,ierr)
      CALL MPI_COMM_RANK(env%comm,env%mpr,ierr)
      env%myproc = env%mpr+1
      
      IF (ierr.NE.0) THEN
        WRITE (*,*) "Could not init parallel environment!!"
      END IF
      
      END

C -----------------------------------------------------------------------------
      SUBROUTINE StopProgram(env,io,inp,ierr)
C
C     $: stops excecution of the programm
C
C -----------------------------------------------------------------------------      
      USE inc_types    
      TYPE(IOtype) :: io
      TYPE(InputType) :: inp
      TYPE(penv) :: env
      INTEGER ierr
       
      IF (env%myproc.EQ.1) THEN
        CALL DatumInUra(inp%endtime)
        WRITE (io%l,'(//A,A)') "START:",inp%StartTime
        WRITE (io%l,'(A,A)') "END  :",inp%EndTime        
        IF (ierr.NE.0) THEN
          WRITE (io%l,'(/A)') "Finished with an error!"
        ELSE
          WRITE (io%l,'(/A)') "Exited normally!"        
        END IF
        CLOSE(io%l)
c       error/warning file      
        CLOSE (io%errwarn)
        CLOSE (io%res) ! tri.dat
        CLOSE (io%rawl) ! tri.walls.dat
        CLOSE (io%mhr) ! tri.mH-rms.dat
        CLOSE (io%dct) ! tri.dc-rms.dat
      END IF
      
      CALL MPI_FINALIZE(env%comm,ierr)
      STOP
      END           
