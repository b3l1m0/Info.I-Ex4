//
// Function for performing the roto-translational fitting of
// a protein structure (target) onto another one (reference).
// Adapted from GROMOS96 Fortran routine lsqstr.f (using eigen.f)
// --- Phil, Nov 2017 ---
//
// Notes:
// - For now, the ref and trg structures will *both* be translated
//   to the origin (maybe later we want to translate trg on ref instead)
// - The function lsqstr and eigen still have Fortran-like arrays
//   starting from 1 (maybe clean up later)
//
// INCLUDES
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "fit.h" // prototype for StructureFit
//
// NAMESPACE
using namespace std;
//
// MACROS
#define NUM_CRD 3 // number of Cartesian dimensions
//
// PROTOTYPES
void str_cen ( int num_atm, double atm_crd[][NUM_CRD] );
void str_rot ( int num_atm, double atm_crd_ref[][NUM_CRD], double atm_crd_trg[][NUM_CRD] );
void lsqstr ( int NR, double W[], double XP[], double X[], double &E, int &IROT, int NDIM );
void eig_sys (double A[], double R[], int N, bool get_vec ); // get eigenvalues (and eigenvectors) of a symmetric matrix
//
// --------------------------------------------------------------------------------
//
void StructureFit ( int num_atm, double atm_crd_ref[][NUM_CRD], double atm_crd_trg[][NUM_CRD] ) {
  str_cen ( num_atm, atm_crd_ref );
  str_cen ( num_atm, atm_crd_trg );
  str_rot ( num_atm, atm_crd_ref, atm_crd_trg );
  return;
}
//
// --------------------------------------------------------------------------------
//
void str_cen ( int num_atm, double atm_crd[][NUM_CRD] ) {
  double x, y, z;
  x = y = z = 0.0;
  for ( int n=0; n<num_atm; n++ ) {
    x += atm_crd[n][0];
    y += atm_crd[n][1];
    z += atm_crd[n][2];
  }
  x /= num_atm;
  y /= num_atm;
  z /= num_atm;
  for ( int n=0; n<num_atm; n++ ) {
    atm_crd[n][0] -= x;
    atm_crd[n][1] -= y;
    atm_crd[n][2] -= z;
  }
}
//
// --------------------------------------------------------------------------------
//
void str_rot ( int num_atm, double atm_crd_ref[][NUM_CRD], double atm_crd_trg[][NUM_CRD] ) {
  double wei[num_atm+1];
  for ( int n=0; n<num_atm+1; n++ ) wei[n] = 1.0;
  double e;
  int irot = 2;
  lsqstr ( num_atm, wei, ((double*)(atm_crd_ref)) - 1, ((double*)(atm_crd_trg)) - 1, e, irot, NUM_CRD );
}
//
// --------------------------------------------------------------------------------
//
void lsqstr ( int NR, double W[], double XP[], double X[], double &E, int &IROT, int NDIM ) {
//
//     lsqstr rotates the atoms with coordinates x about the origin
//     such that the function e = 0.5 * sum over all nr atoms of
//     w*(x-xp)**2 has a minimum. here w denotes the weight factors and
//     xp are the reference coordinates. for a description of the
//     method see a.d. mclachlan, j. mol. biol. 128 (1979) 49.
//     if the subroutine fails, it is returned with a message and
//     irot=0.
//
//     nr = number of atoms
//     w(1..nr) = atomic weight factors
//     xp(1..3*nr) = reference atom coordinates
//     x(1..3*nr) = atom coordinates; delivered with the rotated ones
//     e = delivered with the minimum value of the function e (irot=2)
//     irot = 1 : x is rotated, e is not computed
//          = 2 : x is rotated and e is computed
//     ndim = munber of dimensions (restricted to 3)
//
//     lsqstr uses function eig_sys()
//
  int KDIM, MAXOM, MAXC;
  KDIM = 3;
  MAXOM = 6;
  MAXC = 21;
  //
  double atmp;
  //
  int M1, M2, I, J, M;
  double XJ[KDIM+1], XPJ[KDIM+1], U[KDIM+1][KDIM+1], COM[MAXC+1],
    VH[KDIM+1][KDIM+1], VK[KDIM+1][KDIM+1], TWO, SQT2, DU, SIG, SIGD;
  //
  double OM[MAXOM*MAXOM+1];
  //
  TWO = 2.E0;
  SQT2 = sqrt(TWO);
  //
  // CALCULATE THE MATRIX U AND ITS DETERMINANT
  //
  for ( M2=1; M2<=3; M2++ ) {
    for ( M1=1; M1<=3; M1++ ) {
      U[M1][M2] = 0.E0;
    }
  }
  I = 0;
  for ( J=1; J<=NR; J++ ) {
    for ( M=1; M<=3; M++ ) {
      XJ[M] = X[I+M];
      XPJ[M] = XP[I+M];
    }
    for ( M2=1; M2<=3; M2++ ) {
      for ( M1=1; M1<=3; M1++ ) {
	U[M1][M2] = U[M1][M2]+W[J]*XJ[M1]*XPJ[M2];
      }
    }
    I = I+3;
  }
  //
  DU= U[1][1]*U[2][2]*U[3][3]+U[1][3]*U[2][1]*U[3][2]
    +U[1][2]*U[2][3]*U[3][1]-U[3][1]*U[2][2]*U[1][3]
    -U[3][3]*U[2][1]*U[1][2]-U[3][2]*U[2][3]*U[1][1];
  //
  if ( DU < 1.E-10 && DU > -1.E-10 ) {
    cout << "Error in eig_sys: determinant of U is zero (exiting!)" << endl;
    exit(1);
  }
  if ( DU >= 0.0 ) SIGD = +1.0; else SIGD = -1.0;
  //
  // CONSTRUCT OMEGA, DIAGONALIZE IT AND DETERMINE H AND K
  //
  M = 0;
  for ( M1=1; M1<=6; M1++ ) {
    for ( M2=1; M2<=M1; M2++ ) {
      M=M+1;
      if (M1 > 3 && M2 <= 3) {
	COM[M]=U[M2][M1-3];
      } else
	COM[M]=0.E0;
    }
  }
  //
  eig_sys(COM,OM,6,true);
  //
  atmp = COM[3]-COM[6];
  if (atmp<0) atmp = -atmp;
  if ( DU < 0.0 && atmp < 1.E-5) {
    cout << "Error in eig_sys: determinant of U is negative and Omega has degenerate eigenvalues (exiting!)" << endl;
    exit(1);
  }
  //
  for ( M2=1; M2<=3; M2++ ) {
    for ( M1=1; M1<=3; M1++ ) {
      VH[M1][M2] = SQT2 * OM[(M1-1)+MAXOM*(M2-1)+1];
      VK[M1][M2] = SQT2 * OM[(M1+3-1)+MAXOM*(M2-1)+1];
    }
  }
  SIG = (VH[2][1]*VH[3][2]-VH[3][1]*VH[2][2])*VH[1][3]
    +(VH[3][1]*VH[1][2]-VH[1][1]*VH[3][2])*VH[2][3]
    +(VH[1][1]*VH[2][2]-VH[2][1]*VH[1][2])*VH[3][3];

  if (SIG <= 0.E0) {
    for ( M=1; M<=3; M++ ) {
      VH[M][3]=-VH[M][3];
      VK[M][3]=-VK[M][3];
    }
  }
  //
  // DETERMINE R AND ROTATE X
  //
  for ( M2=1; M2<=3; M2++ ) {
    for ( M1=1; M1<=3; M1++ ) {
      U[M1][M2]=VK[M1][1]*VH[M2][1]+VK[M1][2]*VH[M2][2] + SIGD*VK[M1][3]*VH[M2][3];
    }
  }
  //
  I = 0;
  for ( J=1; J<=NR; J++ ) {
    for ( M=1; M<=3; M++ ) {
      XPJ[M]=X[I+M];
    }
    for ( M1=1; M1<=3; M1++ ) {
      XJ[M1] = 0.E0;
      for ( M2=1; M2<=3; M2++ ) {
	XJ[M1]=XJ[M1]+U[M1][M2]*XPJ[M2];
      }
      X[I+M1]=XJ[M1];
    }
    I=I+3;
  }
  //
  // CALCULATE E, WHEN REQUIRED
  //
  if (IROT < 2) return;
  E = 0.E0;
  I = 0;
  for ( J=1; J<=NR; J++ ) {
    for ( M1=1; M1<=3; M1++ ) {
      I=I+1;
      E = E + W[J]* (X[I]-XP[I])*(X[I]-XP[I]);
    }
  }
  E = E/2.E0;
  //
  return;
}
//
// --------------------------------------------------------------------------------
//
void eig_sys (double A[], double R[], int N, bool get_vec ) {
//
//     subroutine eig_sys (a,r,n,mv)
//
//     eig_sys computes eigenvalues and eigenvectors of the real
//     symmetric n*n matrix a, using the diagonalization method
//     described in "mathematical methods for digital computers", eds.
//     a.ralston and h.s.wilf, wiley, new york, 1962, chapter 7.
//     it has been copied from the ibm scientific subroutine package.
//
//     a(1..n*(n+1)/2) = matrix to be diagonalized, stored in symmetric
//                       storage mode, viz. the i,j-th element (i.ge.j)
//                       is stored at the location k=i*(i-1)/2+j in a;
//                       the eigenvalues are delivered in descending
//                       order on the diagonal, viz. at the locations
//                       k=i*(i+1)/2
//     r(1..n,1..n) = delivered with the corresponding eigenvectors
//                    stored columnwise
//     n = order of matrices a and r
//     mv = 0 : eigenvalues and eigenvectors are computed
//        = 1 : only eigenvalues are computed
//

//     ARGS
//      INTEGER N, MV
//      REAL A(N*(N + 1)/2), R(N * N)

  int I, IA, IQ, IJ, IL, ILR, ILQ, IM, IMR, IMQ, IND, J, JQ, K, L, LL, LQ, LM, MQ, M, MM;
  double RANGE, ANORM, ANRMX, THR, X, Y, SINCS, SINX2, SINX, COSX, COSX2;

  double atmp;

  RANGE = 1.E-12;

  //
  // GENERATE IDENTITY MATRIX
  //
  if ( get_vec ) {
    IQ = -N;
    for ( J = 1; J<=N; J++ ) {
      IQ += N;
      for ( I = 1; I<=N; I++ ) {
	IJ = IQ+I;
	R[IJ] = 0.E0;
	if ( I == J ) {
	  R[IJ] = 1.E0;
	}
      }
    }
  }
  //
  // COMPUTE INITIAL AND FINAL NORMS (ANORM AND ANRMX)
  //
  ANORM = 0.E0;
  for ( I=1; I<=N; I++ ) {
    for ( J=I; J<=N; J++ ) {
      if ( I != J ) {
	IA = I+(J*J-J)/2;
	ANORM = ANORM + A[IA]*A[IA];
      }
    }
  }
  if ( ANORM > 0.0 ) {
    //
    ANORM = 1.414E0 * sqrt(ANORM);
    ANRMX = ANORM*RANGE/(double)N;
    //
    // INITIALIZE INDICATORS AND COMPUTE THRESHOLD, THR
    //
    IND = 0;
    THR = ANORM;
    //
    do {
      THR = THR / (double) N;
      do{
	L = 1;
	do{
	  M = L+1;
	  //
	  // COMPUT SIN AND COS
	  //
	  do {
	    MQ = (M*M-M)/2;
	    LQ = (L*L-L)/2;
	    LM = L+MQ;
	    //
	    atmp = A[LM];
	    if ( atmp < 0.0 ) atmp = -atmp;
	    //
	    if ( atmp-THR >= 0 ) {
	      IND = 1;
	      LL = L+LQ;
	      MM = M+MQ;
	      X = 0.5*(A[LL]-A[MM]);
	      Y = -A[LM]/ sqrt(A[LM]*A[LM]+X*X);
	      if ( X < 0.0 ) {
		Y = -Y;
	      }
	      SINX = Y/sqrt(2.E0*(1.E0+(sqrt(1.E0-Y*Y))));
	      SINX2 = SINX*SINX;
	      COSX = sqrt(1.E0-SINX2);
	      COSX2 = COSX*COSX;
	      SINCS = SINX*COSX;
	      //
	      // ROTATE L AND M COLUMNS
	      //
	      ILQ = N*(L-1);
	      IMQ = N*(M-1);
	      for ( I=1; I<=N; I++ ) {
		IQ = (I*I-I)/2;
		if ( I != L ) {
		  if ( I != M ) {
		    if ( I < M ) {
		      IM = I+MQ;
		    } else {
		      IM = M+IQ;
		    }
		    if ( I < L ) {
		      IL = I+LQ;
		    } else {
		      IL = L+IQ;
		    }
		    X = A[IL]*COSX - A[IM]*SINX;
		    A[IM] = A[IL]*SINX + A[IM]*COSX;
		    A[IL] = X;
		  }
		}
		if ( get_vec ) {
		  ILR = ILQ+I;
		  IMR = IMQ+I;
		  X = R[ILR]*COSX - R[IMR]*SINX;
		  R[IMR]=R[ILR]*SINX + R[IMR]*COSX;
		  R[ILR]=X;
		}
	      }
	      X = 2.E0 * A[LM]*SINCS;
	      Y = A[LL]*COSX2 + A[MM]*SINX2-X;
	      X = A[LL]*SINX2 + A[MM]*COSX2+X;
	      A[LM] = (A[LL]-A[MM])*SINCS + A[LM]*(COSX2-SINX2);
	      A[LL] = Y;
	      A[MM] = X;
	      //
	      // TEST FOR M = LAST COLUMN
	      //
	    }
	  } while ( M++ != N );
	  M--;
	  //
	  // TEST FOR L = SECOND FROM LAST COLUMN
	  //
	} while ( L++ != N-1 );
	L--;
      } while ( IND-- == 1 );
      IND++;
      //
      // COMPARE THRESHOLD WITH FINAL NORM
      //
    } while ( THR > ANRMX );
    //
    // SORT EIGENVALUES AND EIGENVECTORS
    //
  }
  IQ = -N;
  for ( I=1; I<=N; I++ ) {
    IQ = IQ+N;
    LL = I+(I*I-I)/2;
    JQ = N*(I-2);
    for ( J=I; J<=N; J++ ) {
      JQ = JQ+N;
      MM = J+(J*J-J)/2;
      if ( A[LL] < A[MM] ) {
	X = A[LL];
	A[LL] = A[MM];
	A[MM] = X;
	if ( get_vec ) {
	  for ( K=1; K<=N; K++ ) {
	    ILR = IQ+K;
	    IMR = JQ+K;
	    X = R[ILR];
	    R[ILR] = R[IMR];
	    R[IMR] = X;
	  }
	}
      }
    }
  }
  //
  return;
}
