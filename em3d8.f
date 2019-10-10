!> \file em3d8.f
c!! \brief File with subroutine UEL for defining of user-defined element for ABAQUS.
c!!
c!! This file contains
c!! UEL-subroutine for ABAQUS for implementation of an
c!! user-defined piezoelectric element with 8 nodes and 8 integration points
C!!
C!! TU Bergakademie Freiberg
C!!
c!! Institute of mechanics and fluid dynamics
C!!
c!! Chair of Applied Mechanics - Solid Mechanics
C!!
C!! D-09596 Freiberg Germany
C!!
C!! http://tu-freiberg.de/fakult4/imfd/fkm/index.en.html
c> \date 14.04.2011
c> \author  Alexey Bratskikh bratskik@tu-freiberg.de
C!!
C!!
C!! Defining of a piezoelectric material
C!!
C!! mechanical behavior:
C!!
C!! DVstress = DMs*DVstrain - DMp*DVepotg
C!!
C!! electrical behavior:
C!!
C!! DVeflx   = transpose(DMp)*DVstrain + DMd*DVepotg
C!!
C!!(DVstress - the mechanical stress vector,
C!!
C!! DVstrain - the strain vector,
C!!
c!! DVeflx   - the electrical flux vector,
C!!
c!! DVepotg  - the electrical potential gradient vector -d(phi)/d(x(i)),
C!!
c!! DMs      - the matrix of material's elastic stiffness
c!!             defined at zero electrical potential gradient
c!!             (short circuit condition),
C!!
c!! DMp      - the matrix of material's piezoelectric stress coefficient ,
C!!
c!! DMd      - the matrix of material's dielectric (permittivity)
c!!             property for the fully constrained material)
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
C
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
C     {                      MODULE KMyConst
!>    Numerical constants
!!    @see KmatrKadd
      MODULE KMyConst
c
      IMPLICIT NONE
      save
      DOUBLE PRECISION, PARAMETER ::
     $   d0   = 0.0D0, d1   = 1.0D0
     $ , d1_2 = 0.5D0, d1_4 = 0.25D0, d1_8 = 0.125D0
      INTEGER( KIND = 4), PARAMETER ::
     $   i0   = 0                   , i2   = 2      , i3   = 3
     $ , i4   = 4    , i5   = 5     , i6   = 6      , i7   = 7
     $ , i8   = 8    , i9   = 9     , i10  = 10     , i11  = 11
     $ , i12  = 12   , i13  = 13    , i14  = 14     , i15  = 15
     $ , i16  = 16   , i17  = 17    , i18  = 18     , i19  = 19
     $ , i20  = 20   , i21  = 21    , i22  = 22     , i23  = 23
     $ , i24  = 24   , i25  = 25    , i26  = 26     , i27  = 27
     $ , i28  = 28   , i29  = 29    , i30  = 30     , i31  = 31
     $ , i32  = 32   , i33  = 33    , i34  = 34     , i35  = 35
     $ , i36  = 36   , i37  = 37    , i38  = 38     , i39  = 39
     $ , i40  = 40   , i41  = 41    , i42  = 42     , i43  = 43
     $ , i44  = 44   , i45  = 45    , i46  = 46     , i47  = 47
     $ , i48  = 48   , i49  = 49    , i50  = 50     , i51  = 51
     $ , i52  = 52   , i100 = 100   , i200 = 200    , i350 = 350
     $ , i701 = 701  , i998 = 998   , i999 = 999    , i9999= 9999
      DOUBLE PRECISION, PARAMETER, dimension( i3, i3) :: dmunity3
     $       = RESHAPE([ 1, i0, i0, i0, 1, i0, i0, i0, 1], ( / i3, i3/))
c
      END MODULE KMyConst
C     }                     END MODULE KMyConst
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
C     {                         MODULE KMyMatProp
!>    Material properties of element.
!!
!!    The user inputs the values for this variables in Abaqus *.inp-file.
!!    This variables will be initialized in KinitUsrProp.
!!    @see KinitUsrProp
      MODULE KMyMatProp
c
      USE KMyConst, only : i3, i6
      IMPLICIT NONE
!>    materials elastic stiffness matrix 6x6 for 3D case
      DOUBLE PRECISION, DIMENSION( i6, i6) :: DMs
!>    density
      DOUBLE PRECISION                     :: Drho
!>    matrix for dielectric (permittivity) property 3x3 for 3D case
      DOUBLE PRECISION, DIMENSION( i3, i3) :: DMd
!>    matrix for piezoelectric stress coefficient 6x3 for 3D case
      DOUBLE PRECISION, DIMENSION( i6, i3) :: DMp
c
      END MODULE KMyMatProp
c     }                    END MODULE KMyMatProp
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     {                        MODULE KMyAxes
!>    Local coordinate system for orientation of element properties
**    from module KMyMatProp in global coordinate system
**
**    The user inputs values for vectors of x'- and z'-axes in Abaqus *.inp-file.
**    vector of y'-axis will be derived from this values with cross product.
!!    This variables will be initialized in KinitUsrProp.
!!    @see KinitUsrProp
      MODULE KMyAxes
c
      USE KMyConst, only : i3
      IMPLICIT NONE
!>    vector of x'-axis in global coordinate system for orientation of element properties
      DOUBLE PRECISION, DIMENSION( i3)     :: DVaxX
!>    vector of y'-axis in global coordinate system for orientation of element properties
      DOUBLE PRECISION, DIMENSION( i3)     :: DVaxY
!>    vector of z'-axis in global coordinate system for orientation of element properties
      DOUBLE PRECISION, DIMENSION( i3)     :: DVaxZ
!>    matrix x'y'z'-axes in global coordinate system for orientation of element properties
      DOUBLE PRECISION, DIMENSION( i3, i3) :: DMax
c
      END MODULE KMyAxes
c     }                     END MODULE KMyAxes
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     {                  MODULE KMyShapeFn
!>    Number of integration points, weights of integration points
!!    , N-Tensor for displacements, N-Matrix for electrical potentials
!!    and local partial derivatives of shape functions
      MODULE KMyShapeFn
c
      USE KMyConst, only : i2, i3, i4, i5, i8, i24, d0, d1, d1_8
      IMPLICIT NONE
c
c     LOCAL
!>    (local) number of nodes
      INTEGER*1, private, PARAMETER                       :: Nnodes = i8
c
c     NOT LOCAL
!>    number of integration points
      INTEGER*1, PARAMETER                                ::IGNintPt=i8
!>    weights of Gaussian integration points
      DOUBLE PRECISION, save, DIMENSION( IGNintPt)         :: DVGw
!>    N-Tensor for displacements couples nodal displacements
!!    with displasements at integration point (3)*(3*Nnodes)*(IGNintPt)
      DOUBLE PRECISION, save, DIMENSION( i3, i24, IGNintPt):: DTNdis
!>    (out) N-Matrix for electrical potentials couples
!!    nodal electrical potentials with electrical potentials at integration point
      DOUBLE PRECISION, save, DIMENSION( Nnodes, IGNintPt):: DMNepot
!>    local partial derivatives of shape functions with respect to xi
      DOUBLE PRECISION, save, DIMENSION( Nnodes, IGNintPt):: DMGshp_xi
!>    local partial derivatives of shape functions with respect to eta
      DOUBLE PRECISION, save, DIMENSION( Nnodes, IGNintPt):: DMGshp_eta
!>    local partial derivatives of shape functions with respect to zeta
      DOUBLE PRECISION, save, DIMENSION( Nnodes, IGNintPt):: DMGshp_zeta
c
      contains
!>    Initializing of values of variables in module KMyShapeFn
      subroutine KsetMyShapeVal( )
!>    (local) vector for local coordinates of nodes in xi-direction
!!    with values -1.0 or +1.0
      DOUBLE PRECISION, dimension( Nnodes)           :: DVGnodeXi
!>    (local) vector for local coordinates of nodes in eta-direction
!!    with values -1.0 or +1.0
      DOUBLE PRECISION, dimension( Nnodes)           :: DVGnodeEta
!>    (local) vector for local coordinates of nodes in zeta-direction
!!    with values -1.0 or +1.0
      DOUBLE PRECISION, dimension( Nnodes)           :: DVGnodeZeta
!>    (local) variable for coordinates of 8 Gaussian integration points (2 supporting points^3D)
      DOUBLE PRECISION, PARAMETER          :: Dcr1 = 0.577350269189626D0
!>    (local) variable for weights of 8 Gaussian integration points
      DOUBLE PRECISION, PARAMETER                    :: Dwr1 = d1
!>    (local) xi-coordinates of Gaussian integration points
      DOUBLE PRECISION, DIMENSION( IGNintPt)         :: DVGxi
!>    (local) eta-coordinates of Gaussian integration points
      DOUBLE PRECISION, DIMENSION( IGNintPt)         :: DVGeta
!>    (local) zeta-coordinates of Gaussian integration points
      DOUBLE PRECISION, DIMENSION( IGNintPt)         :: DVGzeta
!>    (local) weighted coordinate
      DOUBLE PRECISION                               :: Dxi0, Deta0
     $                                                , Dzeta0
!>    (local) linear shape functions
      DOUBLE PRECISION, DIMENSION( Nnodes, IGNintPt) :: DMGshp
c
!>    (local) loop increment for number of Gaussian point
      INTEGER*1                                      :: iGp
!>    (local) loop increment for node number
      INTEGER*1                                      :: iNnode
c
c
c     { initializing of vectors for local coordinates of nodes with
c     values -1.0 or 1.0
c         eta
c         ^
c         |
c         .->xi
c        /
c       v
c     zeta
C         8 nodes numbering
c      behind            ahead
C     4      3         8       7
C
C     1      2         5       6
      DVGnodeXi(    1: i4) = [ -d1, d1, d1, -d1 ]
      DVGnodeXi(   i5: i8) = [ -d1, d1, d1, -d1 ]
c
      DVGnodeEta(   1: i4) = [ -d1, -d1, d1, d1 ]
      DVGnodeEta(  i5: i8) = [ -d1, -d1, d1, d1 ]
c
      DVGnodeZeta(  1: i4) = -d1
      DVGnodeZeta( i5: i8) =  d1
c     } end initializing of vectors for local coordinates of nodes
c
c
c     { initializing of Gaussian integration points
c         eta
c         ^
c         |
c         .->xi
c        /
c       v
c     zeta
C         8 integration points numbering
c      behind            ahead
C     4      3         8       7
C
C     1      2         5       6
      DVGxi(   1:i4) =  [ -Dcr1, Dcr1, Dcr1, -Dcr1 ]
      DVGxi(  i5:i8) =  [ -Dcr1, Dcr1, Dcr1, -Dcr1 ]
c
      DVGeta(  1:i4) =  [ -Dcr1, -Dcr1, Dcr1, Dcr1 ]
      DVGeta( i5:i8) =  [ -Dcr1, -Dcr1, Dcr1, Dcr1 ]
c
      DVGzeta( 1:i4) = -Dcr1
      DVGzeta(i5:i8) =  Dcr1
c     } end initializing of Gaussian integration points
c
c     { initializing of weights of Gaussian integration points
      DVGw(    1:i8) =  Dwr1
c     } end initializing of weights of Gaussian integration points
c
c
C     { initializing of shape functions and their local partial derivatives,
c       initializing of N-Tensor for displacements and N-Matrix for electrical potentials
      DO iGp = 1, IGNintPt
         DO iNnode = 1, Nnodes
C           weighted coordinates
            Dxi0   = DVGnodeXi(   iNnode) * DVGxi(   iGp)
            Deta0  = DVGnodeEta(  iNnode) * DVGeta(  iGp)
            Dzeta0 = DVGnodeZeta( iNnode) * DVGzeta( iGp)
C           linear shape-function evaluated at integration point
            DMGshp( iNnode, iGp)      = d1_8 * ( d1 + Dxi0)
     $                                       * ( d1 + Deta0)
     $                                       * ( d1 + Dzeta0)
C           local derivative of linear shape-function with respect to xi at integration point
            DMGshp_xi(   iNnode, iGp) = d1_8 * DVGnodeXi(   iNnode)
     $                                       * ( d1 + Deta0)
     $                                       * ( d1 + Dzeta0)
C           local derivative of linear shape-function with respect to eta at integration point
            DMGshp_eta(  iNnode, iGp) = d1_8 * ( d1 + Dxi0)
     $                                       * DVGnodeEta(  iNnode)
     $                                       * ( d1 + Dzeta0)
C           local derivative of linear shape-function with respect to zeta at integration point
            DMGshp_zeta( iNnode, iGp) = d1_8 * ( d1 + Dxi0)
     $                                       * ( d1 + Deta0)
     $                                       * DVGnodeZeta( iNnode)
c           N-Tensor for displacements
            DTNdis( :,  1 + i3*( iNnode - 1), iGp) =
     1         [ DMGshp( iNnode, iGp), d0             , d0             ]
            DTNdis( :, i2 + i3*( iNnode - 1), iGp) =
     1         [ d0             , DMGshp( iNnode, iGp), d0             ]
            DTNdis( :, i3*iNnode, iGp)              =
     1         [ d0             , d0             , DMGshp( iNnode, iGp)]

c           N-Matrix for electrical potentials
            DMNepot( iNnode, iGp) = DMGshp( iNnode, iGp)
         END DO
      END DO
C     } end initializing of shape functions and their local partial derivatives,
c       initializing of N-Tensor for displacements and N-Matrix for electrical potentials
c
      end subroutine KsetMyShapeVal
c
      END MODULE KMyShapeFn
c     }               END MODULE KMyShapeFn
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
c
c
c
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     {                   MAIN SUBROUTINE UEL
c  Only relevant variables are commented. Refer to ABAQUS manual for others.
!> \brief Subroutine to define an user element (interface to ABAQUS)
!!
!! UEL-subroutine for ABAQUS for implementation of an
!! user-defined piezoelectric element with 8 nodes and 8 integration points
!! @param RHS   (out) element residual force or load vector
!!                    (element contribution to the right-hand-side vector of the overall system of equations)
!!                    to solve the equation AMATRX * U= RHS
!! @param AMATRX (out) element stiffness/Jacobian/mass matrix
!!                    (element contribution to the stiffness matrix
!!                    of the overall system of equations)
!!                    to solve the equation AMATRX * U= RHS
!! @param SVARS (out) solution-dependent state variables associated with this element:
!!    6 stress components {s11 s22 s33 s12 s13 s23},
!!    6 strain components {e11 e22 e33 e12 e13 e23},
!!    3 electrical flux and
!!    3 electrical potential gradients
!! @param COORDS (in) nodal coordinate array: COORDS(J,N) is Jth coord of Nth node
!! @param DU     (in) increments of accumulated DOF array
!!                    (of displacements and of electrical potential)
!! @param JPROPS (in) integer valued Element property array
!! @param JTYPE  (in) integer valued element type
!!                  (1 - input piezoelectric stress coefficient matrix conforming Abaqus PIEZOELECTRIC TYPE=S,
!!                  2 - input piezoelectric strain coefficient matrix conforming Abaqus PIEZOELECTRIC TYPE=E)
!! @param LFLAGS (in) flag variables from Abaqus for solving procedure
!! @param MCRD   (in) maximum coordinates <=3.
!! @param NNODE  (in) number of nodes per element 8
!! @param NDOFEL (in) number of degrees of freedom for element 4*8
!!                   (4 = 3 displacements + 1 electrical potential)
!! @param NPROPS (in) number of real valued element properties
!! @param NJPROP (in) number of integer valued element properties
!! @param MLVARX (in) Dimensioning variable for RHS and DU >= NDOFEL
!! @param NRHS   (in) number of RHS vectors
!! @param NSVARS (in) number of solution-dependent state variables
!! @param PROPS  (in) real valued Element property array
!! @param U      (in) Total accumulated DOF array.  Contains accumulated displacements u(j) (j=1,2,3)
!!                    and electrical potentials phi
!!                    ordered as (u_i_x, u_i_y, u_i_z, phi_i) i = 1, ..., NNODE like
!!                    u_1_x, u_1_y, u_1_z, phi_1, u_2_x, u_2_y, u_2_z, phi_2, ...
!!                    to solve the equation AMATRX * U= RHS
!! @see KctrlUsrParam, KDshape, Kjacobi
      SUBROUTINE UEL( RHS, AMATRX
     2              , SVARS, ENERGY, NDOFEL, NRHS, NSVARS, PROPS, NPROPS
     3              , COORDS, MCRD, NNODE, U, DU, V
     4              , A, JTYPE, TIME, DTIME
     5              , KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG
     6              , PREDEF, NPREDF, LFLAGS, MLVARX, DdLMAG, MDLOAD
     7              , PNEWDT, JPROPS, NJPROP, PERIOD)
C     
c      INCLUDE 'ABA_PARAM.INC'
      USE KMyConst
      USE KMyMatProp, only : DMs, Drho, DMd, DMp
      USE KMyShapeFn, only : IGNintPt, DVGw, DTNdis, DMNepot, DMGshp_xi
     $                     , DMGshp_eta, DMGshp_zeta, KsetMyShapeVal
      IMPLICIT NONE
C
c     ABAQUS defined variables (subroutine parameters)
      INTEGER ::
     2      NDOFEL, NRHS, NSVARS, NPROPS
     3    , MCRD, NNODE
     4    , JTYPE
     5    , KSTEP, KINC, JELEM, NDLOAD
     6    , NPREDF, MLVARX, MDLOAD
     7    , NJPROP
c
      INTEGER ::
     5      JDLTYP( MDLOAD, *)
     6    , LFLAGS( *)
     7    , JPROPS( NJPROP)
c
      DOUBLE PRECISION ::
     4      DTIME
     7    , PNEWDT, PERIOD
c
      DOUBLE PRECISION :: RHS( MLVARX, *)    , AMATRX( NDOFEL, NDOFEL)
     2    , SVARS( NSVARS)      , ENERGY( i8), PROPS( NPROPS)
     3    , COORDS( MCRD, NNODE), U( NDOFEL) , DU( MLVARX, *),V( NDOFEL)
     4    , A( NDOFEL)          , TIME( i2)
     5    , PARAMS( *)          , ADLMAG( MDLOAD, *)
     6    , PREDEF( i2, NPREDF, NNODE)       , DdLMAG( MDLOAD, *)
c
c
c     LOCAL OWN DEFINED variables
C     INudofel = MCRD * NNODE = 3*8
      INTEGER*1, PARAMETER                             :: INudofel = i24! max 127
C     global partial derivatives of linear shape function evaluated at integration point in loop
c     calculated in KDshape
      DOUBLE PRECISION, DIMENSION( NNODE)              :: DVGshp_x
     $                                                  , DVGshp_y
     $                                                  , DVGshp_z
c     determinant of Jacobian matrix for derivates of linear shape functions
c     times weight of integration point in loop calculated in KDshape
      DOUBLE PRECISION                               :: DdetJ_t_Gw
c     vector for displacements of nodes initialized in KgetNodeVal
      DOUBLE PRECISION, dimension( MCRD * NNODE)       :: DVdis
c     vector for electric potentials of nodes initialized in KgetNodeVal
      DOUBLE PRECISION, dimension( NNODE)              :: DVepot
c
c     strain components vector at Gaussian integration point in loop
c     calculated in KmatrB
c     notation {e}={e11 e22 e33 e12 e13 e23}
      DOUBLE PRECISION, dimension( i6)                 :: DVstrain
c     stress components vector at Gaussian integration point in loop
c     calculated in KmatrB
c     notation {s}={s11 s22 s33 s12 s13 s23}
      DOUBLE PRECISION, dimension( i6)                 :: DVstress
c     the electrical flux vector at Gaussian integration point in loop
c     calculated in KmatrB
      DOUBLE PRECISION, dimension( MCRD)               :: DVeflx
c     the electrical potential gradient vector at Gaussian integration point in loop
c     calculated in KmatrB     -d(phi)/d(x(i))
      DOUBLE PRECISION, dimension( MCRD)               :: DVepotg
c
c     B-Matrix for displacements in loop 6*(3*8) initialized in KmatrB
c     notation {e}={e11 e22 e33 e12 e13 e23}
      DOUBLE PRECISION, DIMENSION( i6, MCRD * NNODE)   :: DMBdis
c     B-Matrix for electrical potentials in nodes in loop 3*8 initialized in KmatrB
      DOUBLE PRECISION, DIMENSION( MCRD, NNODE)        :: DMBepot
c
c     mass matrix
      DOUBLE PRECISION, DIMENSION( INudofel, INudofel) :: DMMmass
c     K-matrices actualized in KmatrKadd
c     displacement stiffness matrix
      DOUBLE PRECISION, DIMENSION( INudofel, INudofel) :: DMKdisDis
c     piezoelectric coupling matrix
      DOUBLE PRECISION, DIMENSION( INudofel, NNODE   ) :: DMKdisEpot
c     transposed piezoelectric coupling matrix
      DOUBLE PRECISION, DIMENSION( NNODE   , INudofel) :: DMKepotDis
c     dielectric "stiffness" matrix
      DOUBLE PRECISION, DIMENSION( NNODE   , NNODE   ) :: DMKepotEpot
c     vector of mechanical forces
      DOUBLE PRECISION, DIMENSION( INudofel)           :: DVFmech
c     vector of electrical charge
      DOUBLE PRECISION, DIMENSION( NNODE)              :: DVQel
c
c     length of message (string Cmsg)
      INTEGER*2, PARAMETER                            :: ILenMsg = i200
c     variable for messages
      CHARACTER( len = ILenMsg)                       :: Cmsg
c     error number
      INTEGER*2                                       :: Ierr
c     flag for first call of this routine
      LOGICAL, save                                   :: Lcallf = .true.
c     loop increment for number of Gaussian point
      INTEGER*1                                       :: iGp
c     temp variable
      INTEGER*2                                       :: iGp_times_18
c     LOCAL increment variables for loops
      INTEGER*1                                       :: p, p3, p4!column
     $                                                 , k, k3, k4!row
c     (in) number of output medium in Abaqus Standard
c     (7 - output into message file (.msg),
c      6 - output into printed output file (.dat))
      INTEGER*1, PARAMETER                            :: INoutopt = i6
c
c     called functions
      INTEGER*1                                       :: KctrlUsrParam
     $                                                 , KctrlAbqParam
     $                                                 , KctrlUsrProp
     $                                                 , Kamatrx
      LOGICAL                                         :: KDshape
      DOUBLE PRECISION                                :: Kjacobi
c
c
      Ierr = i0
c     {if UEL called for the first time
      if( Lcallf) then
         write( Cmsg, *) JTYPE
         write( Cmsg, *)'First call of UEL for user element with JTYPE='
     $               , trim( adjustl( Cmsg))
         call KoutMsg(Ierr, trim( Cmsg), len( trim( Cmsg)), INoutopt)
c        control of inputs parameteres from user
         Ierr = KctrlUsrParam( Cmsg  , ILenMsg , MCRD  , NNODE, NPROPS
     2                       , NJPROP, IGNintPt, NSVARS, JTYPE)
         call KoutMsg(Ierr, trim( Cmsg), len( trim( Cmsg)), INoutopt)
c        if user input data are wrong
         if( Ierr .ne. i0) goto 900
c
c        if user input data are OK go ahead
c        output information about parameters by first call of UEL
         call KoutInfo( MCRD  , NNODE , NPROPS, NJPROP, IGNintPt
     2        , NSVARS, NDOFEL, MLVARX, JTYPE , NRHS  , INoutopt)
c
c        control of inputs properties from user
         Ierr = KctrlUsrProp( Cmsg , ILenMsg, NPROPS, NJPROP
     $                      , PROPS, JPROPS)
         call KoutMsg(Ierr, trim( Cmsg), len( trim( Cmsg)), INoutopt)
c        if user input data are wrong
         if( Ierr .ne. i0) goto 900
c
c        initialize user properties: DMs, DMp, DMd, Drho
c        , DVaxX, DVaxY, DVaxZ, DMax
c        and put the information about it out
         call KinitUsrProp( PROPS, JPROPS, JTYPE, INoutopt)
c
c        initialize values of variables in module KMyShapeFn
         call KsetMyShapeVal( )
c
c        set variable for first call Lcallf to .false.
         Lcallf = .false.
      end if
c     } end if UEL called for the first time
c
c     control of Abaqus parameter inputs
      Ierr = KctrlAbqParam( Cmsg, ILenMsg, LFLAGS, NRHS)
c
c     if Abaqus solving parameters are not supported by element.
      if( Ierr .ne. i0) then
         call KoutMsg( Ierr, trim( Cmsg), len( trim( Cmsg)), INoutopt)
         goto 900
      end if
c
c     if Abaqus parameter inputs OK go ahead
c
c     RHS and AMATRX should be defined in depence on
c     flag variable LFLAGS(i3) from Abaqus for solving procedure.
c
c     if (LFLAGS( i3) .eq. i4) define current mass matrix only
c     current mass matrix is always requested
c     as initial mass matrix at the start of the analysis
c     if (LFLAGS( i3) .eq. i6) define the current mass matrix and the residual vector
c     for the initial acceleration calculation
c     { define current mass matrix
      if( ( LFLAGS( i3) .eq. i4) .or. ( LFLAGS( i3) .eq. i6)) then
c        initializing of mass matrix with zeros
         DMMmass = d0
c        Gaussian integration loop
         DO iGp = 1, IGNintPt
         !DMBepot will be used hier as just placeholder for 3x3 Jacobian matrix
            DdetJ_t_Gw = DVGw( iGp) *
     $                         Kjacobi( NNODE , MCRD, COORDS
     $                             , DMGshp_xi(   :, iGp)
     $                             , DMGshp_eta(  :, iGp)
     $                             , DMGshp_zeta( :, iGp), DMBepot)
            DMMmass = DMMmass + Drho * DdetJ_t_Gw * (
     2                MATMUL( TRANSPOSE( DTNdis( :, :, iGp))
     3                                 , DTNdis( :, :, iGp)))
         end do
c        sorting of M- (mass) matrix DMMmass into AMATRX in the wright order
         do p = 1, NNODE  !column
           do k = 1, NNODE !row
c             column
              p3 = i3 * p
              p4 = i4 * p
c             row
              k3 = i3 * k
              k4 = i4 * k
c
              AMATRX( ( k4 - i3):( k4 - 1), ( p4 - i3):( p4 - 1)) =
     $                         DMMmass( ( k3 - i2):k3, ( p3 - i2):p3)
           end do
         end do
c     } end define current mass matrix
      else
c     {if stiffness matrix and/or residual vector should should be defined
c        assignment of nodal values
         call KgetNodeVal( MCRD, NNODE, Inudofel, NDOFEL
     2                   , U   , DVdis, DVepot)
c
c        initializing of variables with zeros
c        K-matrices
         DMKdisDis   = d0
         DMKdisEpot  = d0
         DMKepotDis  = d0
         DMKepotEpot = d0
c        forces vectors
         DVFmech     = d0
         DVQel       = d0
c
c        { Gaussian integration loop
         DO iGp = 1, IGNintPt
c
c        initializing of global partial derivatives of linear shape functions and
c        calculating of determinant of Jacobian matrix for derivates of shape functions
c        times weight of integration point
            if( KDshape( NNODE               , MCRD    , COORDS
     $                 , DMGshp_xi(   :, iGp)
     $                 , DMGshp_eta(  :, iGp)
     $                 , DMGshp_zeta( :, iGp)
     $                 , DVGw( iGp)
     $                 , DVGshp_x            , DVGshp_y, DVGshp_z
     $                 , DdetJ_t_Gw)) then
               Ierr = i11
               write( Cmsg, *) 'Calculating of global derivates of'
     $           , ' linear shape functions failed'
     $           , ' at integration point = ', iGp
     $           , ' (determinant of Jacobian matrix is equal to zero).'
               call KoutMsg( Ierr, trim( Cmsg), len( trim( Cmsg))
     $                     , INoutopt)
               goto 900
            end if
c
c           Determination of values at one Gaussian integration point
c           flag for nonlinear geometry:
c           LFLAGS( i2) .eq. i0: small displacement analysis (geometrically linear)
c           LFLAGS( i2) .eq. 1: large-displacement analysis (nonlinear geometric effects)
            call KmatrB( MCRD    , NNODE   , INudofel
     $                      , ( LFLAGS( i2) .eq. 1)
     $                      , DVGshp_x, DVGshp_y, DVGshp_z
     $                      , DVdis   , DMBdis  , DVstrain, DVstress
     $                      , DVepot  , DMBepot , DVepotg , DVeflx)
c
            iGp_times_18 = iGp * i18
c           save solution-dependent state variables associated with this element
c           SVARS:
c           6 stress components, 6 strain components (notation {s}={s11 s22 s33 s12 s13 s23})
c           3 electric flux, 3 electric potentials
            call KsaveSVARS( DVstress, DVstrain
     2                   , DVeflx  , DVepotg
     3                   , SVARS( ( iGp_times_18 - i17) : iGp_times_18))
c
c           add values to matrices M- (mass), K- (stiffness), F- (mechanical force), Q- (electrical charge)
            call KmatrKadd(MCRD    , NNODE      , INudofel
     1                  , DdetJ_t_Gw
     2                  , DMBdis   , DVstress   , DMBepot   , DVeflx
     3                  , DMKdisDis, DMKepotEpot, DMKdisEpot, DMKepotDis
     4                  , DVFmech  , DVQel)
c
         end do
c        } end Gaussian integration loop
c
c        defining residual/load vector RHS and stiffness matrix AMATRX for element
         Ierr = Kamatrx( Cmsg     , ILenMsg
     1                 , MCRD     , NNODE, INudofel, NDOFEL, MLVARX
     2                 , LFLAGS( i3)
     3                 , DMKdisDis, DMKepotEpot, DMKdisEpot, DMKepotDis
     4                 , DVFmech  , DVQel
     5                 , RHS      , AMATRX)
c
         if( Ierr .ne. i0)
     $     call KoutMsg( Ierr, trim( Cmsg), len( trim( Cmsg)), INoutopt)
c
c     }if stiffness matrix and/or residual vector should should be defined
      end if
c
 900  if( Ierr .ne. i0 ) then
         call KoutMsg( i998, 'STOP', i4, INoutopt)
         STOP
      end if
c
      RETURN
      END SUBROUTINE UEL
c     }                  END MAIN SUBROUTINE UEL
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
c
c
c
c
c
c
c
c     {                  SUBROUTINES AND FUNCTIONS
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     {                   INTEGER*1 FUNCTION KctrlUsrParam
!>    Controling of user input parameters
!!
!!    User input parameters will be checked for reasonability.
!!    @return error number for wrong user input parameters
!!    @see UEL
      INTEGER*1 FUNCTION KctrlUsrParam( Cmsg  , ILenMsg , MCRD  , NNODE
     $                        , NPROPS, NJPROP, IGNintPt, NSVARS, JTYPE)
c
      USE KMyConst, only : i0, i2, i3, i4, i5, i6, i7, i8, i18, i52
     $                   , i9999
      IMPLICIT NONE
c
c     INPUT
!>    (in) length of message (length of string Cmsg)
      INTEGER*2, intent( in) :: ILenMsg
!>    (in) user input parameter
      INTEGER  , intent( in) :: MCRD, NNODE, NPROPS, NJPROP, IGNintPt
     $                                             , NSVARS, JTYPE
c
c     OUTPUT
!>    (out) variable for message
      CHARACTER( len = ILenMsg), intent( out) :: Cmsg
c
c     LOCAL variables
c     variable for calculations
      INTEGER*2                 :: iGp_times_18
c     temporary message
      CHARACTER( len = ILenMsg) :: Cmsg_tmp
c
c
      KctrlUsrParam = i0
c
      if( MCRD .ne. i3) then
         write( Cmsg_tmp, *) MCRD
         write( Cmsg, *) 'Model has to be 3-dimensional, not '
     2                 , trim( adjustl( Cmsg_tmp))
     3                 , ' like in Your input.'
         KctrlUsrParam = 1
c
      else if( NNODE .ne. i8) then
         write( Cmsg_tmp, *) NNODE
         write( Cmsg, *)'In *USER ELEMENT statement in .inp-file'
     1                , ' nodes number (NODES=) should be 8, not '
     2                , trim( adjustl( Cmsg_tmp))
     3                , ' like in Your input.'
         KctrlUsrParam = i2
c
      else if( NPROPS .ne. i52) then
         write( Cmsg_tmp, *) NPROPS
         write( Cmsg, *) 'In *USER ELEMENT statement in .inp-file'
     1                 , ' number of user-defined real properties'
     2                 , ' associated with the element (PROPERTIES=)'
     3                 , ' should be 52, not '
     4                 , trim( adjustl( Cmsg_tmp))
     5                 , ' like in Your input.'
         KctrlUsrParam = i3
c
      else if( NJPROP .ne. 1) then
         write( Cmsg_tmp, *) NJPROP
         write( Cmsg, *) 'In *USER ELEMENT statement in .inp-file'
     1                 , ' number of user-defined integer properties'
     2                 , ' associated with the element (IPROPERTIES=)'
     3                 , ' should be 1, not '
     4                 , trim( adjustl( Cmsg_tmp))
     5                 , ' like in Your input.'
         KctrlUsrParam = i4
c
      else if( IGNintPt .ne. i8) then
         write( Cmsg_tmp, *) IGNintPt
         write( Cmsg, *) ' Number of integration points should be 8.'
         KctrlUsrParam = i5
c
      else if( ( JTYPE .lt. 1) .and. ( JTYPE .gt. i9999)) then
         write( Cmsg_tmp, *) JTYPE
         write( Cmsg, *) 'This type TYPE=U'
     2                 , trim( adjustl( Cmsg_tmp))
     3                 , ' from Your input is not supported.'
     4                 ,' See  *USER ELEMENT and *ELEMENT statements'
     5                 , ' in Your .inp-file!'
         KctrlUsrParam = i7
c
      else
         iGp_times_18 = IGNintPt * i18
         if( NSVARS .ne. iGp_times_18) then
            write( Cmsg_tmp, *) IGNintPt
            write( Cmsg, *) 'In *USER ELEMENT statement in .inp-file'
     1                    , ' number of the solution dependent state'
     2                    , ' variables (VARIABLES=) should be'
     3                    , ' 18*number of integration points'
     4                    , ' ( 18*', trim( adjustl( Cmsg_tmp)), ' ='
            write( Cmsg_tmp, *) iGp_times_18
            write( Cmsg( len( trim( Cmsg)) + 1: ), *) ' '
     1                                       , trim( adjustl( Cmsg_tmp))
     2                                       ,' is not equal to'
            write( Cmsg_tmp, *) NSVARS
            write( Cmsg( len( trim( Cmsg)) + 1: ), *) ' '
     1                                       , trim( adjustl( Cmsg_tmp))
     2                                       ,' from Your input).'
            KctrlUsrParam = i6
         else
c           all parameters are OK!
            write( Cmsg, *) 'User parameters are OK.'
            KctrlUsrParam = i0
         end if
      end if
c
      RETURN
      END FUNCTION KctrlUsrParam
c     }                END INTEGER*1 FUNCTION KctrlUsrParam
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
c
c
c
c
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     {                   INTEGER*1 FUNCTION KctrlAbqParam
!>    Controling of ABAQUS parameters for element solving procedure
!!
!!    ABAQUS parameters for element solving procedure will be checked for
!!    supporting by implemented user-defined element.
!!    @return error number for procedures of ABAQUS which are not supported by element
      INTEGER*1 FUNCTION KctrlAbqParam( Cmsg, ILenMsg, LFLAGS, NRHS)
c
      USE KMyConst, only : i0, i2, i3, i4, i8, i9, i10, i100
      IMPLICIT NONE
c
c     INPUT
!>    (in) length of message (length of string Cmsg)
      INTEGER*2, intent( in)              :: ILenMsg
!>    (in) flag variables from Abaqus for solving procedure
      INTEGER, DIMENSION( *), intent( in) :: LFLAGS
!>    (in) number of load vectors
      INTEGER*1, intent( in)              :: NRHS
c
c     OUTPUT
!>    (out) variable for message
      CHARACTER( len = ILenMsg), intent( out) :: Cmsg
c
c
      if( ( LFLAGS( 1) .ne. 1) .and. ( LFLAGS( 1) .ne. i2)) then
         write( Cmsg, *)'Only static procedure is supported by element.'
         KctrlAbqParam = i8
      else if( ( LFLAGS( i4) .eq. 1) .or. ( LFLAGS( i3) .eq. i100)) then
         write( Cmsg, *)'Perturbation step is not supported by element.'
         KctrlAbqParam = i9
      else if( NRHS .ne. 1) then
         write( Cmsg, *)    'RIKS solution is not supported by element.'
         KctrlAbqParam = i10
      else
c        all parameters are OK!
         write( Cmsg, *) 'Abaqus parameters for element solving are OK.'
         KctrlAbqParam = i0
      end if
c
      RETURN
      END FUNCTION KctrlAbqParam
c     }               END INTEGER*1 FUNCTION KctrlAbqParam
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
c
c
c
c
c

c
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     {                   INTEGER*1 FUNCTION KctrlUsrProp
!>    Controling of user input properties
!!
!!    User input properties will be checked for reasonability.
!!    @return error number for wrong user input properties
!!    @see UEL
      INTEGER*1 FUNCTION KctrlUsrProp( Cmsg , ILenMsg, NPROPS, NJPROP
     $                               , PROPS, JPROPS)
c
      USE KMyConst, only : i0 , i2 , i3 , i13, i14, i15, i47, i48, i49
     $                   , i50, i51, i52, d0
      IMPLICIT NONE
c
c     INPUT
!>    (in) length of message (length of string Cmsg)
      INTEGER*2, intent( in)                            :: ILenMsg
!>    (in) number of float properties inputed from user
      INTEGER,   intent( in)                            :: NPROPS
!>    (in) number of integer properties inputed from user
      INTEGER,   intent( in)                            :: NJPROP
!>    (in) from user inputed float properties
      DOUBLE PRECISION, dimension( NPROPS), intent( in) :: PROPS
!>    (in) from user inputed integer properties
      INTEGER, dimension( NJPROP), intent( in)          :: JPROPS
c
c     OUTPUT
!>    (out) variable for message
      CHARACTER( len = ILenMsg), intent( out)           :: Cmsg
c
c     LOCAL variables
c     3-vector for y'-axis
      DOUBLE PRECISION, DIMENSION( i3)                  :: DVaxY
c
c
      if( ( ( PROPS( i47) .eq. d0) .and. ( PROPS( i48) .eq. d0))
     $                             .and. ( PROPS( i49) .eq. d0))    then
         write( Cmsg, *) 'In *UEL PROPERTY statement in .inp-file'
     $                 , ' vector for x-axis should be not ( 0, 0, 0)'
     $                 , ' like in Your input.'
         KctrlUsrProp = i13
      else if( ( ( PROPS( i50) .eq. d0) .and. ( PROPS( i51) .eq. d0))
     $                             .and. ( PROPS( i52) .eq. d0))    then
         write( Cmsg, *) 'In *UEL PROPERTY statement in .inp-file'
     $                 , ' vector for z-axis should be not ( 0, 0, 0)'
     $                 , ' like in Your input.'
         KctrlUsrProp = i14
      else
         call KcrossProd3( PROPS( i47:i49), PROPS( i50:i52), DVaxY)
         if( ( ( DVaxY( 1) .eq. d0) .and. ( DVaxY( i2) .eq. d0))
     $                              .and. ( DVaxY( i3) .eq. d0))    then
            write( Cmsg, *) 'The vector of y-axis should not yield'
     $         , ' ( 0, 0, 0) by deriving from vectors of x- and z-axes'
     $         , ' from *UEL PROPERTY statement in Your .inp-file.'
         KctrlUsrProp = i15
         else
c           all parameters are OK!
            write( Cmsg, *) 'User properties are OK.'
         KctrlUsrProp = i0
         end if
      end if
c
      RETURN
      END FUNCTION KctrlUsrProp
c     }               END INTEGER*1 FUNCTIONKctrlUsrProp
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
c
c
c
c
c
C     assignment of element properties
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
C     declaration of used variables:
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     material properties
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c DOUBLE PRECISION
C     Ds1111, Ds1122, Ds2222, Ds1133, Ds2233, Ds3333, Ds1112, Ds2212
C     Ds3312, Ds1212, Ds1113, Ds2213, Ds3313, Ds1213, Ds1313, Ds1123
C     Ds2223, Ds3323, Ds1223, Ds1323, Ds2323
C        -- 21 elements of material's elastic stiffness matrix
C           defined at zero electrical potential gradient
c           (short circuit condition)
c
c                                           , Drho
c          -- density
c
c                                                   , Dd11,   Dd12
c     Dd22,   Dd13,   Dd23,   Dd33
c          -- 6 elements of matrix for dielectric (permittivity) property
c             for the fully constrained material
c
c                                  ,  Dp111,  Dp122,  Dp133,  Dp112
c     Dp113,  Dp123,  Dp211,  Dp222,  Dp233,  Dp212,  Dp213,  Dp223
c     Dp311,  Dp322,  Dp333,  Dp312,  Dp313,  Dp323
c          -- 18 elements of piezoelectric stress coefficient matrix
c                (if JPROPS( 1) = 0) or
c             or 18 elements of piezoelectric strain coefficient matrix
c                (if JPROPS( 1) = 1)
c                                                     DVaxX1,	DVaxX2
c     DVaxX3
c          -- three values of vector to decribe local x'-axis
c             for orientation of element properties
c             in global coordinates system
c
c           , DVaxZ1, DVaxZ2, DVaxZ3
c          -- three values of vector to decribe local z'-axis
c             for orientation of element properties
c             in global coordinates system
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     {                     SUBROUTINE KinitUsrProp
!>    Assign element properties and put the information about it out
!!
!!    Values from user inputed float properties PROPS
!!    will be copied into DMs, DMd and Drho from module KMyMatProp.
!!    Depence of first value from integer properties JPROPS
!!    values from user inputed float properties PROPS
!!    will be copied into DMp from module KMyMatProp
!!    or this values will be calculated regarding equation:
!!
!!    Dp(mij) = Ds(ijkl) * Dpstrain(mkl)
!!
!!    Values DVaxX, DVaxY, DVaxZ and DMax from module KMyAxes will be
!!    copied or derived from user inputed float properties PROPS
!!    and will be normalized.
      SUBROUTINE KinitUsrProp( PROPS, JPROPS, JTYPE, INoutopt)
c
      USE KMyConst
      USE KMyMatProp, only : DMs, Drho, DMd, DMp
      USE KMyAxes
      IMPLICIT NONE
c
c     INPUT
!>    (in) from user inputed float properties
      DOUBLE PRECISION, dimension( *), intent( in) :: PROPS
!>    (in) from user inputed integer properties
      INTEGER, dimension( *), intent( in)          :: JPROPS
!>    (in) integer valued element type
      INTEGER, intent( in)                         :: JTYPE
!>    (in) number of output medium in Abaqus Standard
!!    (7 - output into message file (.msg),
!!    6 - output into printed output file (.dat))
      INTEGER*1, intent( in)                       :: INoutopt
c
c     CALLED FUNCTION
      DOUBLE PRECISION                             ::KmagVec
c
c     LOCAL
c     21 elements of material's elastic stiffness matrix
c     defined at zero electrical potential gradient
c     (short circuit condition)
      DOUBLE PRECISION :: Ds1111, Ds1122, Ds2222, Ds1133, Ds2233, Ds3333
     $  , Ds1112, Ds2212, Ds3312, Ds1212, Ds1113, Ds2213, Ds3313, Ds1213
     $  , Ds1313, Ds1123, Ds2223, Ds3323, Ds1223, Ds1323, Ds2323
c     6 elements of matrix for dielectric (permittivity) property
c     for the fully constrained material
      DOUBLE PRECISION :: Dd11,  Dd12,  Dd22,  Dd13,  Dd23,  Dd33
c     18 elements of piezoelectric stress coefficient matrix
      DOUBLE PRECISION ::Dp111, Dp122, Dp133, Dp112, Dp113, Dp123, Dp211
     $   , Dp222, Dp233, Dp212, Dp213, Dp223, Dp311, Dp322, Dp333, Dp312
     $   , Dp313, Dp323
c     temporar variable for calculation of reverse value of vector magnitude
      DOUBLE PRECISION :: DRlenV
c     increments for loops
      INTEGER*1                       :: i, j
c
c
      write( INoutopt, *) 'User defined float properties'
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     { DOUBLE PRECISION material properties
c     21 elements of material's elastic stiffness matrix
      Ds1111 = PROPS( 1  ); Ds1122 = PROPS( i2 ); Ds2222 = PROPS( i3 )
      Ds1133 = PROPS( i4 ); Ds2233 = PROPS( i5 ); Ds3333 = PROPS( i6 )
      Ds1112 = PROPS( i7 ); Ds2212 = PROPS( i8 ); Ds3312 = PROPS( i9 )
      Ds1212 = PROPS( i10); Ds1113 = PROPS( i11); Ds2213 = PROPS( i12)
      Ds3313 = PROPS( i13); Ds1213 = PROPS( i14); Ds1313 = PROPS( i15)
      Ds1123 = PROPS( i16); Ds2223 = PROPS( i17); Ds3323 = PROPS( i18)
      Ds1223 = PROPS( i19); Ds1323 = PROPS( i20); Ds2323 = PROPS( i21)
c
c     material's elastic stiffness matrix
      DMs = transpose( reshape( (
     1                 / Ds1111, Ds1122, Ds1133, Ds1112, Ds1113, Ds1123,
     2                   Ds1122, Ds2222, Ds2233, Ds2212, Ds2213, Ds2223,
     3                   Ds1133, Ds2233, Ds3333, Ds3312, Ds3313, Ds3323,
     4                   Ds1112, Ds2212, Ds3312, Ds1212, Ds1213, Ds1223,
     5                   Ds1113, Ds2213, Ds3313, Ds1213, Ds1313, Ds1323,
     6                   Ds1123, Ds2223, Ds3323, Ds1223, Ds1323, Ds2323
     7                 /), shape( DMs)))
c
      write( INoutopt, *) 'material elastic stiffness matrix:'
      do j = 1, i6
         do i = 1, i6
           write( INoutopt, *) 'DMs(' , i, ', ', j, ') = ', DMs( i, j)
         end do
      end do
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     density
      Drho = PROPS( i22)
      write( INoutopt, *) 'density:'
      write( INoutopt, *) 'Drho = ', Drho
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     6 elements of matrix for dielectric (permittivity) properties
      Dd11 = PROPS( i23); Dd12 = PROPS( i24); Dd22 = PROPS( i25)
      Dd13 = PROPS( i26); Dd23 = PROPS( i27); Dd33 = PROPS( i28)
c
c     matrix for dielectric (permittivity) properties
      DMd = transpose( reshape( (
     1                 / Dd11, Dd12, Dd13,
     2                   Dd12, Dd22, Dd23,
     3                   Dd13, Dd23, Dd33
     4                 /), shape( DMd)))
      write( INoutopt, *) 'matrix for dielectric'
     $                    ,' (permittivity) properties:'
      do j = 1, i3
         do i = 1, i3
           write( INoutopt, *) 'DMd(' , i, ', ' ,j, ') = ', DMd( i, j)
         end do
      end do
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     18 elements of piezoelectric stress (if JPROPS( 1) = 0) or
c                                  strain (if JPROPS( 1) = 1) coefficient matrix
      if( JPROPS( 1) .eq. i0) then
c        input is piezoelectric stress coefficient matrix
         Dp111 = PROPS( i29); Dp122 = PROPS( i30); Dp133 = PROPS( i31)
         Dp112 = PROPS( i32); Dp113 = PROPS( i33); Dp123 = PROPS( i34)
         Dp211 = PROPS( i35); Dp222 = PROPS( i36); Dp233 = PROPS( i37)
         Dp212 = PROPS( i38); Dp213 = PROPS( i39); Dp223 = PROPS( i40)
         Dp311 = PROPS( i41); Dp322 = PROPS( i42); Dp333 = PROPS( i43)
         Dp312 = PROPS( i44); Dp313 = PROPS( i45); Dp323 = PROPS( i46)
      else if( JPROPS( 1) .eq. 1) then
c        input is piezoelectric strain coefficient matrix
c        piezoelectric stress coefficient matrix will be calculated
c        regarding equation:
c        Dp(mij) = Ds(ijkl) * Dpstrain(mkl)
c        Dpstrain(mkl) = PROPS(i29:i46)
c        Dpstrain(mkl) - piezoelectric strain coefficient matrix
c        Dp(mij)       - piezoelectric stress coefficient matrix
c        Ds(ijkl)      - material's elastic stiffness matrix
         Dp111= Ds1111*PROPS(i29)+ Ds1122*PROPS(i30) + Ds1133*PROPS(i31)
     $        + Ds1112*PROPS(i32)+ Ds1113*PROPS(i33) + Ds1123*PROPS(i34)
         Dp122= Ds1122*PROPS(i29)+ Ds2222*PROPS(i30) + Ds2233*PROPS(i31)
     $        + Ds2212*PROPS(i32)+ Ds2213*PROPS(i33) + Ds2223*PROPS(i34)
         Dp133= Ds1133*PROPS(i29)+ Ds2233*PROPS(i30) + Ds3333*PROPS(i31)
     $        + Ds3312*PROPS(i32)+ Ds3313*PROPS(i33) + Ds3323*PROPS(i34)
         Dp112= Ds1112*PROPS(i29)+ Ds2212*PROPS(i30) + Ds3312*PROPS(i31)
     $        + Ds1212*PROPS(i32)+ Ds1213*PROPS(i33) + Ds1223*PROPS(i34)
         Dp113= Ds1113*PROPS(i29)+ Ds2213*PROPS(i30) + Ds3313*PROPS(i31)
     $        + Ds1213*PROPS(i32)+ Ds1313*PROPS(i33) + Ds1323*PROPS(i34)
         Dp123= Ds1123*PROPS(i29)+ Ds2223*PROPS(i30) + Ds3323*PROPS(i31)
     $        + Ds1223*PROPS(i32)+ Ds1323*PROPS(i33) + Ds2323*PROPS(i34)
c
         Dp211= Ds1111*PROPS(i35)+ Ds1122*PROPS(i36) + Ds1133*PROPS(i37)
     $        + Ds1112*PROPS(i38)+ Ds1113*PROPS(i39) + Ds1123*PROPS(i40)
         Dp222= Ds1122*PROPS(i35)+ Ds2222*PROPS(i36) + Ds2233*PROPS(i37)
     $        + Ds2212*PROPS(i38)+ Ds2213*PROPS(i39) + Ds2223*PROPS(i40)
         Dp233= Ds1133*PROPS(i35)+ Ds2233*PROPS(i36) + Ds3333*PROPS(i37)
     $        + Ds3312*PROPS(i38)+ Ds3313*PROPS(i39) + Ds3323*PROPS(i40)
         Dp212= Ds1112*PROPS(i35)+ Ds2212*PROPS(i36) + Ds3312*PROPS(i37)
     $        + Ds1212*PROPS(i38)+ Ds1213*PROPS(i39) + Ds1223*PROPS(i40)
         Dp213= Ds1113*PROPS(i35)+ Ds2213*PROPS(i36) + Ds3313*PROPS(i37)
     $        + Ds1213*PROPS(i38)+ Ds1313*PROPS(i39) + Ds1323*PROPS(i40)
         Dp223= Ds1123*PROPS(i35)+ Ds2223*PROPS(i36) + Ds3323*PROPS(i37)
     $        + Ds1223*PROPS(i38)+ Ds1323*PROPS(i39) + Ds2323*PROPS(i40)
c
         Dp311= Ds1111*PROPS(i41)+ Ds1122*PROPS(i42) + Ds1133*PROPS(i43)
     $        + Ds1112*PROPS(i44)+ Ds1113*PROPS(i45) + Ds1123*PROPS(i46)
         Dp322= Ds1122*PROPS(i41)+ Ds2222*PROPS(i42) + Ds2233*PROPS(i43)
     $        + Ds2212*PROPS(i44)+ Ds2213*PROPS(i45) + Ds2223*PROPS(i46)
         Dp333= Ds1133*PROPS(i41)+ Ds2233*PROPS(i42) + Ds3333*PROPS(i43)
     $        + Ds3312*PROPS(i44)+ Ds3313*PROPS(i45) + Ds3323*PROPS(i46)
         Dp312= Ds1112*PROPS(i41)+ Ds2212*PROPS(i42) + Ds3312*PROPS(i43)
     $        + Ds1212*PROPS(i44)+ Ds1213*PROPS(i45) + Ds1223*PROPS(i46)
         Dp313= Ds1113*PROPS(i41)+ Ds2213*PROPS(i42) + Ds3313*PROPS(i43)
     $        + Ds1213*PROPS(i44)+ Ds1313*PROPS(i45) + Ds1323*PROPS(i46)
         Dp323= Ds1123*PROPS(i41)+ Ds2223*PROPS(i42) + Ds3323*PROPS(i43)
     $        + Ds1223*PROPS(i44)+ Ds1323*PROPS(i45) + Ds2323*PROPS(i46)
      end if
c     matrix for piezoelectric stress coefficients
      DMp =  transpose( reshape( (
     1                  / Dp111, Dp211, Dp311
     2                  , Dp122, Dp222, Dp322
     3                  , Dp133, Dp233, Dp333
     4                  , Dp112, Dp212, Dp312
     5                  , Dp113, Dp213, Dp313
     6                  , Dp123, Dp223, Dp323
     7                  /), ( / size( DMp, i2), size( DMp, 1) /)))
      write( INoutopt, *)'matrix for piezoelectric stress coefficients:'
      do j = 1, i3
         do i = 1, i6
           write( INoutopt, *) 'DMp(' , i, ', ', j, ') = ', DMp( i, j)
         end do
      end do
c
c     } end DOUBLE PRECISION material properties
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     { DOUBLE PRECISION values for vectors of x'y'z'-axes
c       for orientation of element properties in global coordinate system
c     normalized vector for x'-axis
      DRlenV = d1/KmagVec( i3, PROPS( i47:i49))
      DVaxX(  1) = PROPS( i47) * DRlenV
      DVaxX( i2) = PROPS( i48) * DRlenV
      DVaxX( i3) = PROPS( i49) * DRlenV
c     normalized vector for z'-axis
      DRlenV = d1/KmagVec( i3, PROPS( i50:i52))
      DVaxZ(  1) = PROPS( i50) * DRlenV
      DVaxZ( i2) = PROPS( i51) * DRlenV
      DVaxZ( i3) = PROPS( i52) * DRlenV
c     normalized vector for y'-axis
      call KcrossProd3( DVaxX, DVaxZ, DVaxY)
      DVaxY = -DVaxY
      DMax = RESHAPE( [DVaxX, DVaxY, DVaxZ]
     $              , ( / i3, i3/))
      write( INoutopt, *) 'vectors of xyz-axes for orientation of'
     $              , 'element properties in global coordinate system:'
      do j = 1, i3
         do i = 1, i3
           write( INoutopt, *) 'DMax(' , i, ', ', j, ') = ', DMax( i, j)
         end do
      end do
c     } end DOUBLE PRECISION values for vectors of x'y'z'-axes
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
      Return
      END SUBROUTINE KinitUsrProp
c     }                  END SUBROUTINE KinitUsrProp
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
c
c
c
c
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     {                  SUBROUTINE KgetNodeVal
!>    Assignment of nodal values
!!
!!    Values from vector U will be copied into vectors DVdis and DVepot
      SUBROUTINE KgetNodeVal( MCRD, NNODE, Inudofel, NDOFEL
     2                      , U   , DVdis, DVepot)
c
      USE KMyConst, only : i3, i4
      IMPLICIT NONE
c
c     INPUT
!>    (in) number of coordinates
      INTEGER*1, intent( in)                                 :: MCRD
!>    (in) nodes number
      INTEGER*1, intent( in)                                 :: NNODE
!>    (in) INudofel = MCRD * NNODE = 3*8
      INTEGER*1, intent( in)                                 :: INudofel
!>    (in) number of degrees of freedom for element 4*8 =
!!    (3 displacements + 1 electrical potential)*number of nodes
      INTEGER*1, intent( in)                                 :: NDOFEL
!>    (in) Total accumulated DOF array.  Contains accumulated displacements
!!         and electrical potentials ordered as
!!         (u_i_x, u_i_y, u_i_z, phi_i) i = 1, ..., NNODE
!!         (u_1_x, u_1_y, u_1_z, phi_1, u_2_x, u_2_y, u_2_z, phi_2, ...)
      DOUBLE PRECISION, DIMENSION( NDOFEL), intent( in)      :: U
c
c     OUTPUT
!>    (out) vector for displacements of nodes (3*8) initialized here
      DOUBLE PRECISION, dimension( INudofel), intent(out)    :: DVdis
!>    (out) vector for electric potentials of nodes initialized here
      DOUBLE PRECISION, dimension( NNODE), intent(out)       :: DVepot
c
c     LOCAL variables
c     loop increment for number of node
      INTEGER*1                                 :: iNnode
c     loop increment of coordinate
      INTEGER*1                                 :: iNcoord
c     numbering variables
      INTEGER*1                                 :: m, k
c
      DO iNnode = 1, NNODE
        do iNcoord = 1, MCRD
          k              = ( iNnode - 1) * i4 + iNcoord
          m              = ( iNnode - 1) * i3 + iNcoord! k - ( iNnode-1)
          DVdis(  m)     = U(  k)
        end do
        DVepot(  iNnode) = U(  k + 1)
      END DO
c
      RETURN
      END SUBROUTINE KgetNodeVal
c     }               END SUBROUTINE KgetNodeVal
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
c
c
c
c
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
C     {                   FUNCTION KDshape
!>    Initializing of global partial derivatives of shape function
!!    and calculating of determinant of Jacobian matrix for
!!    derivates of shape function times Gaussian weight
!!    @return failure of initializing of global derivates of shape function
!!    @see Kjacobi, Kinverse3det, UEL
      LOGICAL FUNCTION KDshape(  NNODE    , MCRD      , COORDS
     $                         , DVGshp_xi, DVGshp_eta, DVGshp_zeta
     $                         , DGw
     $                         , DVGshp_x , DVGshp_y  , DVGshp_z
     $                         , DdetJ_t_Gw)
c
      USE KMyConst, only: i2, i3
      IMPLICIT NONE
c
c     INPUT
!>    (in) number of nodes
      INTEGER*1, intent( in)                                 :: NNODE
!>    (in) number of coordinates
      INTEGER*1, intent( in)                                 :: MCRD
!>    (in) coordinates of nodes
      DOUBLE PRECISION, dimension( MCRD, NNODE), intent( in) :: COORDS
!>    (in) local derivative of shape-functions with respect to xi
!!    at integration point which will be done evaluation for
      DOUBLE PRECISION, DIMENSION( NNODE), intent( in)    :: DVGshp_xi
!>    (in) local derivative of shape-functions with respect to eta
!!    at integration point which will be done evaluation for
      DOUBLE PRECISION, DIMENSION( NNODE), intent( in)    :: DVGshp_eta
!>    (in) local derivative of shape-functions with respect to zeta
!!    at integration point which will be done evaluation for
      DOUBLE PRECISION, DIMENSION( NNODE), intent( in)    :: DVGshp_zeta
!>    (in) weight of integration point which will be done evaluation for
      DOUBLE PRECISION, intent( in)                       :: DGw
c
c     OUTPUT
!>    (out) global derivative of linear shape function evaluated at
!!    integration point which will be done evaluation for
      DOUBLE PRECISION, DIMENSION( NNODE), intent( out)   :: DVGshp_x
     $                                                     , DVGshp_y
     $                                                     , DVGshp_z
!>    (out) determinant of Jacobian matrix for derivates of linear shape functions
!!    times weight of integration point which will be done evaluation for
      DOUBLE PRECISION, intent( out)            :: DdetJ_t_Gw
c
c     LOCAL
c     determinant of Jacobian matrix for derivates of linear shape functions
      DOUBLE PRECISION                          :: DdetJac
c     loop increment for number of node
      INTEGER*1                                 :: iNnode
c     Jacobian and inverse Jacobian matrices for derivates of shape functions
      DOUBLE PRECISION, DIMENSION ( i3, i3)     :: DMjacobi
     $                                           , DMRjacobi
c     CALLED functions
      DOUBLE PRECISION                          :: Kjacobi
      LOGICAL                                   :: Kinverse3det
c     error
      LOGICAL                                   :: Lerr
c
c
      KDshape = .false.
c
c     { local - global coordinate transformation
C     calculation of Jacobian matrix for derivates of linear shape functions
c     and its determinant
      DdetJac = Kjacobi( NNODE    , MCRD      , COORDS
     $                 , DVGshp_xi, DVGshp_eta, DVGshp_zeta
     $                 , DMjacobi)
c
C     calculate inverse Jacobian matrix for derivates of linear shape functions
      Lerr = Kinverse3det( DMjacobi, DdetJac, DMRjacobi)
      if( Lerr) then
         KDshape = .true.
         goto 900
      end if
c     multiplication of determinant of Jacobian matrix for derivates
c     of shape functions with weight of integration point
      DdetJ_t_Gw = DdetJac * DGw
C     } end local - global coordinate transformation
c
C     { global partial derivatives of linear shape functions
      DO iNnode = 1, NNODE
         DVGshp_x( iNnode)= DOT_PRODUCT( DMRjacobi( 1, :),
     $                                     [ DVGshp_xi(   iNnode)
     $                                     , DVGshp_eta(  iNnode)
     $                                     , DVGshp_zeta( iNnode)])
c
         DVGshp_y( iNnode)= DOT_PRODUCT( DMRjacobi( i2, :),
     $                                     [ DVGshp_xi(   iNnode)
     $                                     , DVGshp_eta(  iNnode)
     $                                     , DVGshp_zeta( iNnode)])
c
         DVGshp_z( iNnode)= DOT_PRODUCT( DMRjacobi( i3, :),
     $                                     [ DVGshp_xi(   iNnode)
     $                                     , DVGshp_eta(  iNnode)
     $                                     , DVGshp_zeta( iNnode)])
      END DO
C     } end of global partial derivatives of linear shape functions
c
 900  RETURN
      END FUNCTION KDshape
C     }                    END FUNCTION KDshape
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
c
c
c
c
c
c
c
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
C     {                   SUBROUTINE KmatrB
!> Determination of values at one Gaussian integration point
!!
!! B-matrices (DMBdis and DMBepot) for coupling of nodal displacements DVdis
!! and electrical potentials DVepot
!! with strains DVstrain and electrical potential gradients DVepotg
!! in integration point respectively
!! will be set up with help of global derivates of shape functions
!! DVGshp_x, DVGshp_y, DVGshp_z evaluated at Gaussian integration point.
!!
!! Strain components DVstrain and electrical potential gradients DVepotg
!! at evaluated Gaussian integration point will be calculated
!! with help of B-matrices DMBdis and DMBepot respectively.
!!
!! Stress components DVstress and electrical flux vector DVeflx
!! at evaluating Gaussian integration point will be calculated
!! with help of strain components DVstrain
!! and electrical potential gradients DVepotg and 
!! with help of material properties DMs, DMp and DMd from module KMyMatProp
!! in accordance with equations of defining material behavior:
!!
!! DVstress = DMs*DVstrain - DMp*DVepotg
!!
!! DVeflx = transpose(DMp)*DVstrain + DMd*DVepotg
      SUBROUTINE KmatrB( MCRD    , NNODE   , INudofel
     $                      , Lnonlin
     $                      , DVGshp_x, DVGshp_y, DVGshp_z
     $                      , DVdis   , DMBdis  , DVstrain, DVstress
     $                      , DVepot  , DMBepot , DVepotg , DVeflx)
c
      USE KMyConst, only: i2, i3, i4, i5, i6, i9, d0, d1, d1_2, dmunity3
      USE KMyMatProp, only : DMs, DMd, DMp
      IMPLICIT NONE
c
c     INPUT
!>    (in) number of coordinates
      INTEGER*1, intent( in)                              :: MCRD
!>    (in) number of nodes
      INTEGER*1, intent( in)                              :: NNODE
!>    (in) MCRD * NNODE = 3*8
      INTEGER*1, intent( in)                              :: INudofel
!>    (in) flag for nonlinear geometry:
!!    .false. - small displacement analysis (geometrically linear),
!!    .true. - large-displacement analysis (nonlinear geometric effects)
      LOGICAL, intent( in)                                :: Lnonlin
!>    (in) global derivative of shape functions evaluated at Gaussian integration point
      DOUBLE PRECISION, DIMENSION( NNODE)   , intent( in) :: DVGshp_x
     $                                                     , DVGshp_y
     $                                                     , DVGshp_z
!>    (in) vector for displacements of nodes (3*8)
      DOUBLE PRECISION, dimension( INudofel), intent( in) :: DVdis
!>    (in) vector for electric potentials of nodes
      DOUBLE PRECISION, dimension( NNODE)   , intent( in) :: DVepot
c
c     OUTPUT
!>    (out) B-Matrix for displacements.
!!    This matrix couples the nodal displacements
!!    with the strains in integration point in the element (6*(3*8))
      DOUBLE PRECISION, DIMENSION( i6, INudofel), intent( out):: DMBdis
!>    (out) B-Matrix for electrical potentials at nodes.
!!    This matrix couples the electrical nodal potentials
!!    with the electrical field in integration point in the element (3*8)
      DOUBLE PRECISION, DIMENSION( MCRD, NNODE), intent( out):: DMBepot
!>    (out) strain components vector at Gaussian integration point which will be done evaluation for,
!!    notation {e}={e11 e22 e33 e12 e13 e23}
      DOUBLE PRECISION, dimension( i6), intent( out)         :: DVstrain
!>    (out) stress components vector at Gaussian integration point which will be done evaluation for,
!!    notation {s}={s11 s22 s33 s12 s13 s23}
      DOUBLE PRECISION, dimension( i6), intent( out)         :: DVstress
!>    (out) the electrical potential gradient vector
!!    at Gaussian integration point which will be done evaluation for, -d(phi)/d(x(i))
      DOUBLE PRECISION, dimension( MCRD), intent( out)       :: DVepotg
!>    (out) the electrical flux vector at Gaussian integration point which will be done evaluation for
      DOUBLE PRECISION, dimension( MCRD), intent( out)       :: DVeflx
c
c     LOCAL
c     deformation gradient            3*3
      DOUBLE PRECISION, DIMENSION( MCRD, MCRD)     :: DMdefgr
c     vector with 9 strain elements which will be transfered into vector DVstrain
      DOUBLE PRECISION, dimension( i9)             :: DV9strain
c     loop increment for number of node
      INTEGER*1                                    :: iNnode
c     increment for loops
      INTEGER*1                                    :: p
c
c
c     { if (Lnonlin .eq. .false.)
      IF( Lnonlin .eq. .false.) THEN
C     { geometrically linear B-matrix for displacements 
c     This matrix couples the node displacements
c     with the strains in integration points in the element,
c     notation {e}={e11 e22 e33 e12 e13 e23}
         DO iNnode = 1, NNODE
          DMBdis( :, 1 + i3*( iNnode - 1)) =
     1     [ DVGshp_x( iNnode), d0               , d0
     2     , DVGshp_y( iNnode), DVGshp_z( iNnode), d0                ]
          DMBdis( :, i2 + i3*( iNnode - 1)) =
     1     [ d0               , DVGshp_y( iNnode), d0
     2     , DVGshp_x( iNnode), d0               , DVGshp_z( iNnode)]
          DMBdis( :, i3*iNnode)              =
     1     [ d0               , d0               , DVGshp_z( iNnode)
     2     , d0               , DVGshp_x( iNnode), DVGshp_y( iNnode)]
         END DO
C     } end geometrically linear displacement B-matrices
      ELSE
C     { geometrically non-linear
C     { deformation gradient
         DMdefgr          = d0
         DMdefgr(  1,  1) = d1
         DMdefgr( i2, i2) = d1
         DMdefgr( i3, i3) = d1
         DO p = 1, i3
            DMdefgr( p,  1) = DMdefgr( p,  1) + DOT_PRODUCT(
     $                      DVGshp_x, DVdis( p:i3 * NNODE - i3 + p:i3))
            DMdefgr( p, i2) = DMdefgr( p, i2) + DOT_PRODUCT(
     $                      DVGshp_y, DVdis( p:i3 * NNODE - i3 + p:i3))
            DMdefgr( p, i3) = DMdefgr( p, i3) + DOT_PRODUCT(
     $                      DVGshp_z, DVdis( p:i3 * NNODE - i3 + p:i3))
         END DO
C     } end deformation gradient
C     { geometrically non-linear displacement B-matrix
c       (total Lagrangian formulation)
         DO iNnode = 1, NNODE
            DMBdis( :, 1 + i3 * ( iNnode - 1)) =
     1            [ DMdefgr( 1, 1 ) * DVGshp_x( iNnode),
     2              DMdefgr( 1, i2) * DVGshp_y( iNnode),
     3              DMdefgr( 1, i3) * DVGshp_z( iNnode),

     4              DMdefgr( 1, i2) * DVGshp_x( iNnode)
     4                      + DMdefgr( 1, 1 ) * DVGshp_y( iNnode),

     5              DMdefgr( 1, i3) * DVGshp_x( iNnode)
     5                      + DMdefgr( 1, 1 ) * DVGshp_z( iNnode),

     6              DMdefgr( 1, i3) * DVGshp_y( iNnode)
     6                      + DMdefgr( 1, i2) * DVGshp_z( iNnode)]

            DMBdis( :, i2 + i3 * ( iNnode - 1)) =
     1            [ DMdefgr( i2, 1) * DVGshp_x( iNnode),
     2              DMdefgr( i2, i2) * DVGshp_y( iNnode),
     3              DMdefgr( i2, i3) * DVGshp_z( iNnode),

     4              DMdefgr( i2, i2) * DVGshp_x( iNnode)
     4                      + DMdefgr( i2, 1 ) * DVGshp_y( iNnode),

     5              DMdefgr( i2, i3) * DVGshp_x( iNnode)
     5                      + DMdefgr( i2, 1 ) * DVGshp_z( iNnode),

     6              DMdefgr( i2, i3) * DVGshp_y( iNnode)
     6                      + DMdefgr( i2, i2) * DVGshp_z( iNnode)]

            DMBdis( :, i3 * iNnode) =
     1            [ DMdefgr( i3, 1 ) * DVGshp_x( iNnode),
     2              DMdefgr( i3, i2) * DVGshp_y( iNnode),
     3              DMdefgr( i3, i3) * DVGshp_z( iNnode),

     4              DMdefgr( i3, i2) * DVGshp_x( iNnode)
     4                      + DMdefgr( i3, 1 ) * DVGshp_y( iNnode),

     5              DMdefgr( i3, i3) * DVGshp_x( iNnode)
     5                      + DMdefgr( i3, 1 ) * DVGshp_z( iNnode),

     6              DMdefgr( i3, i3) * DVGshp_y( iNnode)
     6                      + DMdefgr( i3, i2) * DVGshp_z( iNnode)]
         END DO
C     } end geometrically non-linear displacement B-matrix
C     } end geometrically non-linear
      END IF
C     } end if (Lnonlin .eq. i0)
c
c
C     { strain
c     { IF (Lnonlin .eq. .false.)
      IF( Lnonlin .eq. .false.) THEN
C     { geometrically linear: small strain tensor
         DO p = 1, i6
            DVstrain( p) = DOT_PRODUCT( DMBdis( p, :), DVdis)
         END DO
C     } end geometrically linear: small strain tensor
      ELSE
C     { geometrically non-linear: Green-Lagrange strain
c       (total Lagrangian formulation)
             DV9strain = RESHAPE( d1_2 * ( MATMUL(
     $                    TRANSPOSE( DMdefgr), DMdefgr) - dmunity3),
     $                   (/ i9/))
             DVstrain(  1:i3) = DV9strain( [ 1, i5, i9])
             DVstrain( i4:i6) = i2 * DV9strain( [ i6, i3, i2])
C     } end geometrically non-linear: Green-Lagrange strain
      END IF
c     } end IF (Lnonlin .eq. .false.)
c     } end strain
c
c
c     { electrical behavior
C     { B-matrix for electrical potentials
c     This matrix couples the electrical node potentials
c     with the electrical field in integration points in the element.
      DO iNnode = 1, NNODE
        DMBepot( :, iNnode) = [ DVGshp_x( iNnode),
     2                          DVGshp_y( iNnode),
     3                          DVGshp_z( iNnode)]
      end do
c     } end B-matrix for electrical potential gradients
c
c     { the electrical potential gradient vector, -d(phi)/d(x(i))
      DO p = 1, MCRD
         DVepotg( p) = -d1 * DOT_PRODUCT( DMBepot( p, :), DVepot)
      END DO
c     } end the electrical potential gradient vector
c
c     { electrical flux
      DVeflx = MATMUL( DMd, DVepotg)
     2         + MATMUL( TRANSPOSE( DMp), DVstrain)
c     } end electrical flux
c     } end electrical behavior
c
c
c     { mechanical stress
      DVstress = MATMUL( DMs, DVstrain)
     2           - MATMUL( DMp, DVepotg)
c     } end mechanical stress
c
      Return
      END SUBROUTINE KmatrB
C     }                 END SUBROUTINE KmatrB
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
c
c
c
c
c
c
c
c
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     {                  SUBROUTINE KmatrKadd
!> Add values to K- (stiffness), F- (mechanical force) and Q- (electrical charge) matrices
!!
!! Contributions to K- (stiffness) matrices DMKdisDis, DMKdisEpot, DMKepotDis and DMKepotEpot
!! at evaluated integration point will be calculated
!! with help of material properties DMs, DMp and DMd from module KMyMatProp,
!! with help of weight of this integration point
!! multiplied with determinant of Jacobian matrix DdetJ_t_Gw
!! and with help of B-matrices DMBdis and DMBepot
!! for coupling of nodal displacements and electrical potentials
!! with strains and electrical potential gradients in integration point respectively.
      SUBROUTINE KmatrKadd(MCRD, NNODE, INudofel, DdetJ_t_Gw
     2                  , DMBdis   , DVstress   , DMBepot   , DVeflx
     3                  , DMKdisDis, DMKepotEpot, DMKdisEpot, DMKepotDis
     4                  , DVFmech  , DVQel)
c
      USE KMyConst  , only : i6
      USE KMyMatProp, only : DMs, DMd, DMp
      IMPLICIT NONE
c
c     INPUT
!>    (in) number of coordinates
      INTEGER*1, intent( in)                                 :: MCRD
!>    (in) number of nodes
      INTEGER*1, intent( in)                                 :: NNODE
!>    (in) MCRD * NNODE = 3*8
      INTEGER*1, intent( in)                                 :: INudofel
!>    (in) Gaussian integration weight at integration point which will be done evaluation for
!!    multiplied with determinant of Jacobian matrix
      DOUBLE PRECISION, intent( in)                        :: DdetJ_t_Gw
!>    (in) B-Matrix for displacements (6*(3*8))
!!    with notation {e}={e11 e22 e33 e12 e13 e23}
      DOUBLE PRECISION, DIMENSION( i6, INudofel), intent( in):: DMBdis
!>    (in) B-Matrix for electrical potentials in nodes (3*8)
      DOUBLE PRECISION, DIMENSION( MCRD, NNODE), intent( in) :: DMBepot
!>    (in) stress components vector at Gaussian integration point
!!    with notation {s}={s11 s22 s33 s12 s13 s23}
      DOUBLE PRECISION, dimension( i6), intent( in)          :: DVstress
!>    (in) the electrical flux vector at Gaussian integration point
      DOUBLE PRECISION, dimension( MCRD), intent( in)        :: DVeflx
c
c     INPUT/OUTPUT actualizing
c     K-matrices
!>    (in,out) displacement stiffness matrix
      DOUBLE PRECISION, DIMENSION( INudofel, INudofel), intent( inout)
     $                                                 :: DMKdisDis
!>    (in,out) piezoelectric coupling matrix
      DOUBLE PRECISION, DIMENSION( INudofel, NNODE   ), intent( inout)
     $                                                 :: DMKdisEpot
!>    (in,out) transposed piezoelectric coupling matrix
      DOUBLE PRECISION, DIMENSION( NNODE   , INudofel), intent( inout)
     $                                                 :: DMKepotDis
!>    (in,out) dielectric "stiffness" matrix
      DOUBLE PRECISION, DIMENSION( NNODE   , NNODE   ), intent( inout)
     $                                                 :: DMKepotEpot
!>    (in,out) vector of mechanical forces (3*8)
      DOUBLE PRECISION, DIMENSION( INudofel), intent( inout) :: DVFmech
!>    (in,out) vector of electrical charge
      DOUBLE PRECISION, DIMENSION( NNODE)   , intent( inout) :: DVQel
c
c
c     K-matrices
      DMKdisDis   = DMKdisDis   + DdetJ_t_Gw *(
     2  MATMUL( MATMUL( TRANSPOSE( DMBdis),  DMs)            , DMBdis))
c
      DMKdisEpot  = DMKdisEpot  + DdetJ_t_Gw *(
     2  MATMUL( MATMUL( TRANSPOSE( DMBdis),  DMp)            , DMBepot))
c
      DMKepotDis  = DMKepotDis  + DdetJ_t_Gw *(
     2  MATMUL( MATMUL( TRANSPOSE( DMBepot), TRANSPOSE( DMp)), DMBdis))
c
      DMKepotEpot = DMKepotEpot - DdetJ_t_Gw *(
     2  MATMUL( MATMUL( TRANSPOSE( DMBepot), TRANSPOSE( DMd)), DMBepot))
c
c     forces vectors
c      DVFmech = d0
c
c      DVQel   = d0
c
      RETURN
      END SUBROUTINE KmatrKadd
c
c     }                  END SUBROUTINE KmatrKadd
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
c
c
c
c
c
c
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     {                  SUBROUTINE KsaveSVARS
!>    Assembly of 18 solution-dependent state variables associated with this element
!!    at Gaussian integration point which will be done evaluation for
!!
!!    6 stress components DVstress {s11 s22 s33 s12 s13 s23},
!!    6 strain components DVstrain {e11 e22 e33 e12 e13 e23},
!!    3 electrical flux DVeflx and
!!    3 electrical potential gradients DVepotg will be copied into the sector of
!!    vector field for solution-dependent state variables associated with this element.
      SUBROUTINE KsaveSVARS( DVstress, DVstrain
     2                       , DVeflx  , DVepotg
     3                       , SVARS18)
c
      USE KMyConst, only : i3, i6, i12, i15, i18
      IMPLICIT NONE
c
c     INPUT
!>    (in) stress components vector at Gaussian integration point,
!!    notation {s}={s11 s22 s33 s12 s13 s23}
      DOUBLE PRECISION, dimension( i6), intent( in)  :: DVstress
!>    (in) strain components vector at Gaussian integration point,
!!    notation {e}={e11 e22 e33 e12 e13 e23}
      DOUBLE PRECISION, dimension( i6), intent( in)  :: DVstrain
!>    (in) the electrical flux vector at Gaussian integration point
      DOUBLE PRECISION, dimension( i3), intent( in)  :: DVeflx
!>    (in) the electrical potential gradient vector at Gaussian integration point, -d(phi)/d(x(i))
      DOUBLE PRECISION, dimension( i3), intent( in)  :: DVepotg
c
c     OUTPUT
!>    (out) solution-dependent state variables associated with this element
      DOUBLE PRECISION, dimension( i18), intent(out) :: SVARS18
c
c     LOCAL increment for loop
      INTEGER*1                                      :: p
c
         DO p = 1, i6
            SVARS18( p)       = DVstress( p)
            SVARS18( p + i6)  = DVstrain( p)
         END DO
         DO p = 1, i3
            SVARS18( p + i12) = DVeflx(  p)
            SVARS18( p + i15) = DVepotg( p)
         END DO
c
       END SUBROUTINE KsaveSVARS
c     }                 END SUBROUTINE KsaveSVARS
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
c
c
c
c
c
c
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     {                INTEGER*1 FUNCTION Kamatrx
!>    Defining residual or load vector RHS and stiffness or mass matrix AMATRX for element
!!
!!    In dependence on flag variable Itodefine=LFLAGS(i3) from Abaqus for solving procedure
!!    RHS and AMATRX will be defined.
!!
!!    If stiffness matrix have to be difened then
!!    K- (stiffness) matrices DMKdisDis, DMKdisEpot, DMKepotDis and DMKepotEpot
!!    will be sorted into AMATRX in the wright order.
!!
!!    For defining of residual or load vector
!!    F- (mechanical force) and Q- (electrical charge) matrices DVFmech and DVQel
!!    will be sorted into RHS in the wright order.
!!    @return error number for procedures of ABAQUS which are not supported by element
      INTEGER*1 FUNCTION Kamatrx( Cmsg          , ILenMsg
     1                  , MCRD, NNODE, INudofel , NDOFEL    , MLVARX
     2                  , Itodefine
     3                  , DMKdisDis, DMKepotEpot, DMKdisEpot, DMKepotDis
     4                  , DVFmech  , DVQel
     5                  , RHS      , AMATRX)
c
c
      USE KMyConst, only : i0, i2, i3, i4, i5, i6, i9, i100, d0, d1_4
      IMPLICIT NONE
c
c     INPUT
!>    (in) length of message (string Cmsg)
      INTEGER*2, intent( in)                                :: ILenMsg
!>    (in) number of coordinates
      INTEGER*1, intent( in)                                :: MCRD
!>    (in) number of nodes
      INTEGER*1, intent( in)                                :: NNODE
!>    (in) MCRD * NNODE = 3*8
      INTEGER*1, intent( in)                                :: INudofel
!>    (in) 4*8 = number of DOF of element * NNODE =
!!    (3 displacements + 1 electric displacement)* number of nodes
      INTEGER*1, intent( in)                                :: NDOFEL
!>    (in) dimensioning variable for RHS
      INTEGER  , intent( in)                                :: MLVARX
!>    (in) flag variable LFLAGS( i3) from Abaqus for solving procedure
!!    how to define RHS, AMATRX
      INTEGER*1, intent( in)                                :: Itodefine
c     K-matrices
!>    (in) displacement stiffness matrix
      DOUBLE PRECISION, DIMENSION( INudofel, INudofel), intent( in)
     $                                                    :: DMKdisDis
!>    (in) piezoelectric coupling matrix
      DOUBLE PRECISION, DIMENSION( INudofel, NNODE   ), intent( in)
     $                                                    :: DMKdisEpot
!>    (in) transposed piezoelectric coupling matrix
      DOUBLE PRECISION, DIMENSION( NNODE   , INudofel), intent( in)
     $                                                    :: DMKepotDis
!>    (in) dielectric "stiffness" matrix
      DOUBLE PRECISION, DIMENSION( NNODE   , NNODE   ), intent( in)
     $                                                    :: DMKepotEpot
c     forces
!>    (in) vector of mechanical forces (3*8)
      DOUBLE PRECISION, DIMENSION( INudofel), intent( in) :: DVFmech
!>    (in) vector of electrical charge
      DOUBLE PRECISION, DIMENSION( NNODE)   , intent( in) :: DVQel
c
c     OUTPUT
!>    (out) variable for message
      CHARACTER( len = ILenMsg), intent( out)             :: Cmsg
!>    (out) residual or load vector RHS
!!    ordered as (f_i_x, f_i_y, f_i_z, q_i) i = 1, ..., NNODE :
!!               f_1_x, f_1_y, f_1_z, q_1, f_2_x, f_2_y, f_2_z, q_2, ...
!!               (f - mechanical forces, q - electrical charges)
      DOUBLE PRECISION,DIMENSION(MLVARX, *)    ,intent(inout):: RHS
!>    (out) stiffness or mass matrix (4*8)*(4*8)
      DOUBLE PRECISION,DIMENSION(NDOFEL,NDOFEL),intent(inout):: AMATRX
c
c     LOCAL increment variables for loops
      INTEGER*1                                        :: p, p3, p4!column
     $                                                  , k, k3, k4!row
c
c
      Kamatrx = i0
c
      if( Itodefine .eq. 1) then
c     { normal implicit time incremental procedure
c     { define the residual vector in RHS
c       and the stiffness matrix in AMATRX
c
c          RHS( NDOFEL:MLVARX, 1) = d0
c          RHS( 1:NDOFEL     , 1) = d0
          do p = 1, NNODE  !column
           do k = 1, NNODE !row
c             column
              p3 = i3 * p
              p4 = i4 * p
c             row
              k3 = i3 * k
              k4 = i4 * k
c              RHS( ( k4 - i3):( k4 - 1), 1) = DVFmech( ( k3 - i2):k3)
c              RHS( k4                  , 1) = DVQel( k)
c
              AMATRX( ( k4 - i3):( k4 - 1), ( p4 - i3):( p4 - 1)) =
     $                         DMKdisDis( ( k3 - i2):k3, ( p3 - i2):p3)
c
              AMATRX( ( k4 - i3):( k4 - 1), p4)                   =
     $                         DMKdisEpot( ( k3 - i2):k3, p)
c
              AMATRX( k4, ( p4 - i3):( p4 - 1))                   =
     $                         DMKepotDis( k, ( p3 - i2):p3)
c
              AMATRX( k4, p4)                                     =
     $                         DMKepotEpot( k, p)
           end do
          end do
c
c     } end define the residual vector in RHS and the stiffness matrix in AMATRX
c     } end normal implicit time incremental procedure
      else if ( Itodefine .eq. i2) then
c     { define the current stiffness matrix only
         do p = 1, NNODE  !column
           do k = 1, NNODE !row
c             column
              p3 = i3 * p
              p4 = i4 * p
c             row
              k3 = i3 * k
              k4 = i4 * k
c
              AMATRX( ( k4 - i3):( k4 - 1), ( p4 - i3):( p4 - 1)) =
     $                         DMKdisDis( ( k3 - i2):k3, ( p3 - i2):p3)
c
              AMATRX( ( k4 - i3):( k4 - 1), p4)                   =
     $                         DMKdisEpot( ( k3 - i2):k3, p)
c
              AMATRX( k4, ( p4 - i3):( p4 - 1))                   =
     $                         DMKepotDis( k, ( p3 - i2):p3)
c
              AMATRX( k4, p4)                                     =
     $                         DMKepotEpot( k, p)
           end do
         end do
c     } end define the current stiffness matrix only
c
      else if( Itodefine .eq. i5) then
c     { define the current residual or load vector RHS only
c       Half step residual calculation which is needed only for automatic time
c       incrementation
c          RHS( NDOFEL:MLVARX, 1) = d0
c          RHS( 1:NDOFEL     , 1) = d0
!          do k = 1, NNODE
!              k4 = i4 * k
!              k3 = i3 * k
!              RHS( ( k4 - i3):( k4 - 1), 1) = DVFmech( ( k3 - i2):k3)
!              RHS( k4                  , 1) = DVQel( k)
!          end do
c     } end define the current residual or load vector RHS
c
c     next procedures are NOT SUPPORTED by element at this moment
      else if ( Itodefine .eq. i3) then
c     { define the current damping matrix only 
c         AMATRX = d0
         do k = 1, NDOFEL
            AMATRX( k, k) = d1_4
         end do
c     } end define the current damping matrix
c
      else if( Itodefine .eq. i100) then
c     { define perturbation quantities for output
c       force vector that contains
c       the difference between external perturbation loads and internal
c       perturbation forces
c
c     } end define perturbation quantities for output
         write( Cmsg, *)'Perturbation step is not supported by element.'
         Kamatrx = i9
      else
         write( Cmsg, *)'This procedure is not supported by element.'
         Kamatrx = i100
      end if
c
 900  RETURN
      END FUNCTION Kamatrx
c
c     }                END INTEGER*1 FUNCTION Kamatrx
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
c
c
c
c
c
c
c
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     {                     SUBROUTINE KoutInfo
!>    Ouput of information about parameters by first call of UEL
!!    in one of Abaqus output media (.msg or .dat file dependent on value of INoutopt)
      SUBROUTINE KoutInfo(MCRD  , NNODE , NPROPS, NJPROP     , IGNintPt
     2                  , NSVARS, NDOFEL, MLVARX, JTYPE, NRHS, INoutopt)
c
      IMPLICIT NONE
c
c     INPUT
!>    (in) parameter by first call of UEL
      INTEGER  , intent( in) :: MCRD  , NNODE,  NPROPS, NJPROP, IGNintPt
     2                        , NSVARS, NDOFEL, MLVARX, JTYPE,  NRHS
!>    (in) number of output medium in Abaqus Standard
!!    (7 - output into message file (.msg),
!!     6 - output into printed output file (.dat))
      INTEGER*1, intent( in) :: INoutopt
c
c
      write( INoutopt, *) MCRD  , ' - number of dimensions'
      write( INoutopt, *) NRHS  , ' - number of RHS vectors'
      write( INoutopt, *) NNODE , ' - number of nodes on element'
      write( INoutopt, *) NPROPS, ' - number of real user properties'
      write( INoutopt, *) NJPROP, ' - number of integer user properties'
      write( INoutopt, *)IGNintPt,' - number of integration points'
      write( INoutopt, *) NSVARS
     $               , ' - number of solution-dependent state variables'
      write( INoutopt, *) NDOFEL, ' - degrees of freedom for element'
      write( INoutopt, *) MLVARX
     $                      , ' - dimensioning parameter for RHS and DU'
      write( INoutopt, *) JTYPE , ' - element type JTYPE'
c
      return
      end SUBROUTINE KoutInfo
c     }                  END SUBROUTINE KoutInfo
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
c
c
c
c
c
c
c
c
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     {                     SUBROUTINE KoutMsg
!>    Ouput of message Cmsg in one of Abaqus output media
!!    (.msg or .dat file dependent on value of INoutopt)
      SUBROUTINE KoutMsg( INmsg, Cmsg, ILenMsg, INoutopt)
c
      USE KMyConst, only: i0, i12, i350, i701, i999
      IMPLICIT NONE
c
c     INPUT
!>    (in) number of message
      INTEGER*2, intent( in)                 :: INmsg
!>    (in) length of message
      INTEGER*2 , intent( in)                :: ILenMsg
!>    (in)  message
      CHARACTER( len = ILenMsg), intent( in) :: Cmsg
!>    (in) number of output medium in Abaqus Standard
!!    (7 - output into message file (.msg),
!!     6 - output into printed output file (.dat))
      INTEGER*1, intent( in)                 :: INoutopt
c
c     LOCAL temporary message
      CHARACTER( len = i12)                  :: Cmsg_tmp
c
c
      write( Cmsg_tmp, *) INmsg
c
      if( ( INmsg .eq. i0) .or.( INmsg .eq. i999)) then
c     0 and 999 mean OK
c         write( INoutopt, *) 'OK: '
      else if(( INmsg .gt. i0) .and.( INmsg .le. i350)) then
c        error messages from 1 to 350
         write( INoutopt, *) 'error ', trim( adjustl( Cmsg_tmp)), ':'
      else if(( INmsg .gt. i350) .and.( INmsg .lt. i701)) then
c     warning messages from 351 to 700
         write( INoutopt, *) 'warning ', trim( adjustl( Cmsg_tmp)), ':'
      else
c     information messages from 701 to 998
         write( INoutopt, *) 'info ', trim( adjustl( Cmsg_tmp)), ':'
      end if
      write( INoutopt, *) trim( Cmsg)
c
      RETURN
      END SUBROUTINE KoutMsg
c     }                  END SUBROUTINE KoutMsg
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
c     }              END SUBROUTINES AND FUNCTIONS
c
c
c
c
c
c
c
c
c
C     { UTILITY SUBROUTINES
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     {                      FUNCTION KDET3
!>    Calculating determinant of 3x3 matrix DM
!!    @return determinant of 3x3 matrix DM
!!    @see Kjacobi
      DOUBLE PRECISION FUNCTION Kdet3( DM)
c
      USE KMyConst, only: i2, i3
      IMPLICIT NONE
c
c     INPUT
!>    (in) 3x3 matrix
      DOUBLE PRECISION, DIMENSION( i3, i3), intent( in) :: DM
c
      Kdet3 =  DM( 1, 1 ) * DM( i2, i2) * DM( i3, i3)
     2       - DM( 1, i3) * DM( i2, i2) * DM( i3, 1 )
     3       + DM( 1, i2) * DM( i2, i3) * DM( i3, 1 )
     4       + DM( 1, i3) * DM( i2, 1 ) * DM( i3, i2)
     5       - DM( 1, 1 ) * DM( i2, i3) * DM( i3, i2)
     6       - DM( 1, i2) * DM( i2, 1 ) * DM( i3, i3)
c
      RETURN
      END FUNCTION Kdet3
C     }                      END FUNCTION KDET3
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     {                 FUNCTION KINVERSE3DET
!>    Inverting 3x3 matrix DM with its determinant Ddet
!!    and storing result in 3x3 matrix DMr
!!    @return failure of matrix inverting
!!    @see KDshape
      LOGICAL FUNCTION Kinverse3det( DM, Ddet, DMr)
c
      USE KMyConst, only: i2, i3, d0, d1
      IMPLICIT NONE
c
c     INPUT
!>    (in) 3x3 matrix to invert
      DOUBLE PRECISION, DIMENSION ( i3, i3), intent( in)  :: DM
!>    (in) determinant of 3x3 matrix DM which is to invert
      DOUBLE PRECISION, intent( in)                       :: Ddet
c     OUTPUT
!>    (out) 3x3 inverted matrix
      DOUBLE PRECISION, DIMENSION ( i3, i3), intent( out) :: DMr
c     LOCAL variable 1/determinant
      DOUBLE PRECISION                                    :: Drdet
c
c
      if( Ddet .eq. d0) then 
         Kinverse3det = .true.
         DMr = d0
      else
      Kinverse3det = .false.
      DMr( 1, 1) = DM(i2,i2)*DM(i3,i3) - DM(i2,i3)*DM(i3,i2)
      DMr(i2, 1) = DM(i2,i3)*DM(i3, 1) - DM(i2, 1)*DM(i3,i3)
      DMr(i3, 1) = DM(i2, 1)*DM(i3,i2) - DM(i2,i2)*DM(i3, 1)
      DMr( 1,i2) = DM( 1,i3)*DM(i3,i2) - DM( 1,i2)*DM(i3,i3)
      DMr(i2,i2) = DM( 1, 1)*DM(i3,i3) - DM( 1,i3)*DM(i3, 1)
      DMr(i3,i2) = DM( 1,i2)*DM(i3, 1) - DM( 1, 1)*DM(i3,i2)
      DMr( 1,i3) = DM( 1,i2)*DM(i2,i3) - DM( 1,i3)*DM(i2,i2)
      DMr(i2,i3) = DM( 1,i3)*DM(i2, 1) - DM( 1, 1)*DM(i2,i3)
      DMr(i3,i3) = DM( 1, 1)*DM(i2,i2) - DM( 1,i2)*DM(i2, 1)
c
      Drdet =    d1 / Ddet
      DMr   = Drdet * DMr
c
      end if
c
      RETURN
      END FUNCTION Kinverse3det
C     }              END FUNCTION KINVERSE3DET
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
C     {                   INTEGER FUNCTION Kjacobi
!>    Calculating of Jacobian matrix DMjacobi and its determinant
!!    @return determinant of Jacobian matrix
!!    @see Kdet3, KDshape, UEL
      DOUBLE PRECISION FUNCTION Kjacobi( NNODE, MCRD      , COORDS
     $                             , DVGshp_xi, DVGshp_eta, DVGshp_zeta
     $                             , DMjacobi)
c
      IMPLICIT NONE
c
c     INPUT
!>    (in) number of nodes
      INTEGER*1, intent( in)                                 :: NNODE
!>    (in) number of coordinate axes
      INTEGER*1, intent( in)                                 :: MCRD
!>    (in) coordinates of nodes
      DOUBLE PRECISION, DIMENSION( MCRD, NNODE), intent( in) :: COORDS
!>    (in) partial derivatives with respect to xi
      DOUBLE PRECISION, DIMENSION( NNODE), intent( in)    :: DVGshp_xi
!>    (in) partial derivatives with respect to eta
      DOUBLE PRECISION, DIMENSION( NNODE), intent( in)    :: DVGshp_eta
!>    (in) partial derivatives with respect to zeta
      DOUBLE PRECISION, DIMENSION( NNODE), intent( in)    :: DVGshp_zeta
c     OUTPUT
!>    (out) Jacobian matrix for local partial derivatives of shape functions
      DOUBLE PRECISION, DIMENSION( MCRD, MCRD), intent( out) :: DMjacobi
c
c     CALLED functions
      DOUBLE PRECISION                                       :: Kdet3
c
c
C     calculation of Jacobian matrix for derivates of shape functions
      DMjacobi = MATMUL( TRANSPOSE(
     $                    RESHAPE( [ DVGshp_xi, DVGshp_eta, DVGshp_zeta]
     $                  , ( / NNODE, MCRD/)))
     $                  , TRANSPOSE( COORDS))
c     calculate determinant of Jacobian matrix for derivates of shape functions
      Kjacobi = Kdet3( DMjacobi)
c
      RETURN
      END FUNCTION Kjacobi
C     }                    END FUNCTION Kjacobi
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c     {                 SUBROUTINE KcrossProd3
!>    Calculating the right-handed 3-vector cross product DVc
!!    of two 3-vectors DVa and DVb: DVc = DVa x DVb.
      SUBROUTINE KcrossProd3( DVa, DVb, DVc)
c
      USE KMyConst, only : i2, i3
      IMPLICIT NONE
c
!>    (in) multiplicand 3-vector
      DOUBLE PRECISION, DIMENSION( i3), intent (  in) :: DVa
!>    (in) multiplicand 3-vector
      DOUBLE PRECISION, DIMENSION( i3), intent (  in) :: DVb
!>    (out) 3-vector as result of 3-vector cross product
      DOUBLE PRECISION, DIMENSION( i3), intent ( out) :: DVc
c
c
      DVc(  1) = DVa( i2) * DVb( i3) - DVa( i3) * DVb( i2)
      DVc( i2) = DVa( i3) * DVb(  1) - DVa(  1) * DVb( i3)
      DVc( i3) = DVa(  1) * DVb( i2) - DVa( i2) * DVb(  1)
c
      RETURN
      END SUBROUTINE KcrossProd3
c     }               END SUBROUTINE KcrossProd3
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
c
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
C     {                   INTEGER FUNCTION KmagVec3
!>    Calculating of magnitude of vector DV of size Isize
!!    @return  magnitude of vector DV of size Isize
!!    @see KinitUsrProp
      DOUBLE PRECISION FUNCTION KmagVec( Isize, DV)
c
      USE KMyConst, only : d0
      IMPLICIT NONE
c
c     INPUT
!>    (in) size of vector
      INTEGER*1, intent( in)                           :: Isize!max 127
!>    (in) vector
      DOUBLE PRECISION, DIMENSION( Isize), intent( in) :: DV
c
c     LOCAL increment variable for loop
      INTEGER*1                                        :: i
c
      KmagVec = d0
      do i = 1, Isize
         KmagVec = KmagVec + DV( i)*DV( i)
      end do
      KmagVec = sqrt( KmagVec)
c
      RETURN
      END FUNCTION KmagVec
c     }                END FUNCTION KmagVec
C---.|---1----.----2----.----3----.----4----.----5----.----6----.----7--|-.----8
C     } END UTILITY SUBROUTINES
c
c
c
c
c
c
c
c
c
c
