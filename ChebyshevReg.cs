using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ConsoleApplication1
{
    class ChebyshevReg
    {
        static double r_sign(double a, double b)
        {
            double x;
            x = (a >= 0 ? a : -a);
            return (b >= 0 ? x : -x);
        }

        /* Table of constant values */

        static double c_b44 = 1.0;

        static  int cheb_(int m, int n, int mdim, int ndim,
            double[,] a, double[] b, double tol, double relerr, double[] x, ref int
            rank, ref double resmax, ref int iter)
        {
            /* Initialized data */

            double big = 1e32f;

            /* System generated locals */
            //int a_dim1, a_offset,
            int i__1, i__2;
            double r__1;

            /* Local variables */
            double d__;
            int i__, j, k;
            double dd;
            int mm1, kp1, mp1, np1, np2, np3;
            mm1 = 0;
            double val;
            int lev, mode, pcol, prow, np1mk, np1mr;
            pcol = 0;
            double pivot;
            int rankp1;
            double reltmp, tpivot;

            /* THIS SUBROUTINE USES A MODIFICATION OF THE SIMPLEX METHOD */
            /* OF LINEAR PROGRAMMING TO CALCULATE A CHEBYSHEV SOLUTION TO */
            /* AN OVER-DETERMINED SYSTEM OF LINEAR EQUATIONS. */
            /* DESCRIPTION OF PARAMETERS. */
            /* M      NUMBER OF EQUATIONS. */
            /* N      NUMBER OF UNKNOWNS (N MUST NOT EXCEED M). */
            /* MDIM   THE NUMBER OF COLUMNS OF A, AT LEAST M+1. */
            /* NDIM   THE NUMBER OF ROWS OF A, AT LEAST N+3. */
            /* A      TWO DIMENSIONAL double ARRAY OF SIZE (NDIM,MDIM). */
            /*        ON ENTRY,THE TRANSPOSE OF THE MATRIX OF */
            /*        COEFFICIENTS OF THE OVER-DETERMINED SYSTEM MUST */
            /*        BE STORED IN THE FIRST M COLUMNS AND N ROWS OF A. */
            /*        THESE VALUES ARE DESTROYED BY THE SUBROUTINE. */
            /* B      ONE DIMENSIONAL double ARRAY OF SIZE MDIM. ON ENTRY, */
            /*        B MUST CONTAIN THE RIGHT-HAND SIDES OF THE */
            /*        EQUATIONS IN ITS FIRST M LOCATIONS. ON EXIT, B */
            /*        CONTAINS THE RESIDUALS FOR THE EQUATIONS IN ITS */
            /*        FIRST M LOCATIONS (SEE DESCRIPTION). */
            /* TOL    A SMALL POSITIVE TOLERANCE. EMPIRICAL EVIDENCE */
            /*        SUGGESTS TOL=10**(-D+1) WHERE D REPRESENTS THE */
            /*        NUMBER OF DECIMAL DIGITS OF ACCURACY AVAILABLE */
            /*        (SEE DESCRIPTION). */
            /* RELERR A double VARIABLE WHICH ON ENTRY MUST HAVE THE VALUE */
            /*        0.0 IF A CHEBYSHEV SOLUTION IS REQUIRED. IF RELERR */
            /*        IS POSITIVE, THE SUBROUTINE CALCULATES AN */
            /*        APPROXIMATE SOLUTION WITH RELERR AS AN UPPER BOUND */
            /*        ON THE RELATIVE ERROR OF ITS LARGEST RESIDUAL (SEE */
            /*        DESCRIPTION-INEQUALITY (2)). ON EXIT, THE VALUE OF */
            /*        RELERR GIVES A SMALLER UPPER BOUND FOR THIS */
            /*        RELATIVE ERROR. */
            /* X      ONE DIMENSIONAL double ARRAY OF SIZE NDIM. ON EXIT, */
            /*        THIS ARRAY CONTAINS A SOLUTION TO THE PROBLEM IN */
            /*        ITS FIRST N LOCATIONS. */
            /* RANK   AN int WHICH GIVES ON EXIT THE RANK OF THE */
            /*        MATRIX OF COEFFICIENTS. */
            /* RESMAX THE LARGEST RESIDUAL IN MAGNITUDE. */
            /* ITER   THE NUMBER OF SIMPLEX ITERATIONS PERFORMED. */
            /* OCODE  AN EXIT CODE WITH VALUES.. */
            /*              0 - OPTIMAL SOLUTION WHICH IS PROBABLY */
            /*                  NON-UNIQUE (SEE DESCRIPTION). */
            /*              1 - UNIQUE OPTIMAL SOLUTION. */
            /*              2 - CALCULATIONS TERMINATED PREMATURELY */
            /*                  DUE TO ROUNDING ERRORS. */
            /* IF YOUR FORTRAN COMPILER PERMITS A SINGLE COLUMN OF A TWO */
            /* DIMENSIONAL ARRAY TO BE PASSED TO A ONE DIMENSIONAL ARRAY */
            /* THROUGH A SUBROUTINE CALL, CONSIDERABLE SAVINGS IN */
            /* EXECUTION TIME MAY BE ACHIEVED THROUGH THE USE OF THE */
            /* FOLLOWING SUBROUTINE WHICH OPERATES ON COLUMN VECTORS. */
            /*     SUBROUTINE COL (V1,V2,MLT,NOTROW,I1,NP2) */
            /* THIS SUBROUTINE SUBTRACTS FROM THE VECTOR V1 A MULTIPLE OF */
            /* THE VECTOR V2 STARTING AT THE I1*TH ELEMENT UP TO THE */
            /* NP2*TH ELEMENT, EXCEPT FOR THE NOTROW*TH ELEMENT. */
            /*     double V1(NP2),V2(NP2),MLT */
            /*     DO 1 I=I1,NP2 */
            /*       IF(I.EQ.NOTROW) GO TO 1 */
            /*       V1(I)=V1(I)-MLT*V2(I) */
            /*   1   CONTINUE */
            /*     RETURN */
            /*     END */
            /* SEE COMMENTS FOLLOWING STATEMENT NUMBER 340 FOR */
            /* INSTRUCTIONS ON THE IMPLEMENTATION OF THIS MODIFICATION. */
            /* BIG MUST BE SET EQUAL TO ANY VERY LARGE double CONSTANT. */
            /* ITS VALUE HERE IS APPROPRIATE FOR THE IBM/370/145. */
            /* Parameter adjustments */
            //--b;
            //--x;
            //a_dim1 = ndim;
            //a_offset = 1 + a_dim1;
            //a -= a_offset;
            int ocode = 1;
            /* Function Body */
            pcol = 0;
            prow = 0;
            rankp1 = 0;
            /* INITIALIZATION. */
            mp1 = m + 1;
            np1 = n + 1;
            np2 = n + 2;
            np3 = n + 3;
            np1mr = 1;
            rank = n;
            reltmp = relerr;
            relerr = 0.0;

            for (j = 1; j <= m; ++j)
            {
                a[np1, j] = 1.0;
                a[np2, j] = -b[j];
                a[np3, j] = (double)(n + j);
                /* L10: */
            }
            a[np1, mp1] = 0.0;
            iter = 0;
            ocode = 1;
            for (i__ = 1; i__ <= n; ++i__)
            {
                x[i__] = 0.0;
                a[i__, mp1] = (double)i__;
                /* L20: */
            }
            /* LEVEL 1. */
            lev = 1;
            k = 0;
        L30:
            ++k;
            kp1 = k + 1;
            np1mk = np1 - k;
            mode = 0;

            for (j = k; j <= m; ++j)
            {
                b[j] = 1.0;
                /* L40: */
            }
        /* DETERMINE THE VECTOR TO ENTER THE BASIS. */
        L50:
            d__ = -big;

            for (j = k; j <= m; ++j)
            {
                if (b[j] == 0.0)
                {
                    goto L60;
                }
                dd = Math.Abs(a[np2, j]);
                if (dd <= d__)
                {
                    goto L60;
                }
                pcol = j;
                d__ = dd;
            L60:
                ;
            }
            if (k > 1)
            {
                goto L70;
            }
            /* TEST FOR ZERO RIGHT-HAND SIDE. */
            if (d__ > tol)
            {
                goto L70;
            }
            resmax = 0.0;
            mode = 2;
            goto L380;
        /* DETERMINE THE VECTOR TO LEAVE THE BASIS. */
        L70:
            d__ = tol;

            for (i__ = 1; i__ <= np1mk; ++i__)
            {
                dd = Math.Abs(a[i__, pcol]);
                if (dd <= d__)
                {
                    goto L80;
                }
                prow = i__;
                d__ = dd;
            L80:
                ;
            }
            if (d__ > tol)
            {
                goto L330;
            }
            /* CHECK FOR LINEAR DEPENDENCE IN LEVEL 1. */
            b[pcol] = 0.0;
            if (mode == 1)
            {
                goto L50;
            }

            for (j = k; j <= m; ++j)
            {
                if (b[j] == 0.0)
                {
                    goto L100;
                }
                i__2 = np1mk;
                for (i__ = 1; i__ <= i__2; ++i__)
                {
                    if (Math.Abs(a[i__, j]) <= tol)
                    {
                        goto L90;
                    }
                    mode = 1;
                    goto L50;
                L90:
                    ;
                }
            L100:
                ;
            }
            rank = k - 1;
            np1mr = np1 - rank;
            ocode = 0;
            goto L160;
        L110:
            if (pcol == k)
            {
                goto L130;
            }
            /* INTERCHANGE COLUMNS IN LEVEL 1. */

            for (i__ = 1; i__ <= np3; ++i__)
            {
                d__ = a[i__, pcol];
                a[i__, pcol] = a[i__, k];
                a[i__, k] = d__;
                /* L120: */
            }
        L130:
            if (prow == np1mk)
            {
                goto L150;
            }
            /* INTERCHANGE ROWS IN LEVEL 1. */

            for (j = 1; j <= mp1; ++j)
            {
                d__ = a[prow, j];
                a[prow, j] = a[np1mk, j];
                a[np1mk, j] = d__;
                /* L140: */
            }
        L150:
            if (k < n)
            {
                goto L30;
            }
        L160:
            if (rank == m)
            {
                goto L380;
            }
            rankp1 = rank + 1;
            /* LEVEL 2. */
            lev = 2;
            /* DETERMINE THE VECTOR TO ENTER THE BASIS */
            d__ = tol;

            for (j = rankp1; j <= m; ++j)
            {
                dd = Math.Abs(a[np2, j]);
                if (dd <= d__)
                {
                    goto L170;
                }
                pcol = j;
                d__ = dd;
            L170:
                ;
            }
            /* COMPARE CHEBYSHEV ERROR WITH TOL. */
            if (d__ > tol)
            {
                goto L180;
            }
            resmax = 0.0;
            mode = 3;
            goto L380;
        L180:
            if (a[np2, pcol] < -(tol))
            {
                goto L200;
            }
            a[np1, pcol] = 2.0 - a[np1, pcol];

            for (i__ = np1mr; i__ <= np3; ++i__)
            {
                if (i__ == np1)
                {
                    goto L190;
                }
                a[i__, pcol] = -a[i__, pcol];
            L190:
                ;
            }
        /* ARRANGE FOR ALL ENTRIES IN PIVOT COLUMN */
        /* (EXCEPT PIVOT) TO BE NEGATIVE. */
        L200:

            for (i__ = np1mr; i__ <= n; ++i__)
            {
                if (a[i__, pcol] < tol)
                {
                    goto L220;
                }

                for (j = 1; j <= m; ++j)
                {
                    a[np1, j] += a[i__, j] * 2.0;
                    a[i__, j] = -a[i__, j];
                    /* L210: */
                }
                a[i__, mp1] = -a[i__, mp1];
            L220:
                ;
            }
            prow = np1;
            goto L330;
        L230:
            if (rankp1 == m)
            {
                goto L380;
            }
            if (pcol == m)
            {
                goto L250;
            }
            /* INTERCHANGE COLUMNS IN LEVEL 2. */

            for (i__ = np1mr; i__ <= np3; ++i__)
            {
                d__ = a[i__, pcol];
                a[i__, pcol] = a[i__, m];
                a[i__, m] = d__;
                /* L240: */
            }
        L250:
            mm1 = m - 1;
            /* LEVEL 3. */
            lev = 3;
        /* DETERMINE THE VECTOR TO ENTER THE BASIS. */
        L260:
            d__ = -(tol);
            val = a[np2, m] * 2.0;

            for (j = rankp1; j <= mm1; ++j)
            {
                if (a[np2, j] >= d__)
                {
                    goto L270;
                }
                pcol = j;
                d__ = a[np2, j];
                mode = 0;
                goto L280;
            L270:
                dd = val - a[np2, j];
                if (dd >= d__)
                {
                    goto L280;
                }
                mode = 1;
                pcol = j;
                d__ = dd;
            L280:
                ;
            }
            if (d__ >= -(tol))
            {
                goto L380;
            }
            dd = -d__ / a[np2, m];
            if (dd >= reltmp)
            {
                goto L290;
            }
            relerr = dd;
            mode = 4;
            goto L380;
        L290:
            if (mode == 0)
            {
                goto L310;
            }

            for (i__ = np1mr; i__ <= np1; ++i__)
            {
                a[i__, pcol] = a[i__, m] * 2.0 - a[i__, pcol];
                /* L300: */
            }
            a[np2, pcol] = d__;
            a[np3, pcol] = -a[np3, pcol];
        /* DETERMINE THE VECTOR TO LEAVE THE BASIS. */
        L310:
            d__ = big;

            for (i__ = np1mr; i__ <= np1; ++i__)
            {
                if (a[i__, pcol] <= tol)
                {
                    goto L320;
                }
                dd = a[i__, m] / a[i__, pcol];
                if (dd >= d__)
                {
                    goto L320;
                }
                prow = i__;
                d__ = dd;
            L320:
                ;
            }
            if (d__ < big)
            {
                goto L330;
            }
            ocode = 2;
            goto L380;
        /* PIVOT ON A(PROW,PCOL). */
        L330:
            pivot = a[prow, pcol];

            for (j = 1; j <= m; ++j)
            {
                a[prow, j] /= pivot;
                /* L340: */
            }
            /* IF PERMITTED, USE SUBROUTINE COL IN THE DESCRIPTION */
            /* SECTION AND REPLACE THE FOLLOWING EIGHT STATEMENTS DOWN TO */
            /* AND INCLUDING STATEMENT NUMBER 360 BY.. */
            /*     DO 360 J=1,M */
            /*       IF(J.EQ.PCOL) GO TO 360 */
            /*       CALL COL (A(1,J),A(1,PCOL),A(PROW,J),PROW,NP1MR,NP2) */
            /* 360 CONTINUE */

            for (j = 1; j <= m; ++j)
            {
                if (j == pcol)
                {
                    goto L360;
                }
                d__ = a[prow, j];
                i__2 = np2;
                for (i__ = np1mr; i__ <= i__2; ++i__)
                {
                    if (i__ == prow)
                    {
                        goto L350;
                    }
                    a[i__, j] -= d__ * a[i__, pcol];
                L350:
                    ;
                }
            L360:
                ;
            }
            tpivot = -pivot;

            for (i__ = np1mr; i__ <= np2; ++i__)
            {
                a[i__, pcol] /= tpivot;
                /* L370: */
            }
            a[prow, pcol] = 1.0 / pivot;
            d__ = a[prow, mp1];
            a[prow, mp1] = a[np3, pcol];
            a[np3, pcol] = d__;
            ++(iter);
            switch (lev)
            {
                case 1: goto L110;
                case 2: goto L230;
                case 3: goto L260;
            }
        /* PREPARE OUTPUT. */
        L380:

            for (j = 1; j <= m; ++j)
            {
                b[j] = 0.0;
                /* L390: */
            }
            if (mode == 2)
            {
                goto L450;
            }
            i__1 = rank;
            for (j = 1; j <= rank; ++j)
            {
                k = (int)a[np3, j];
                x[k] = a[np2, j];
                /* L400: */
            }
            if (mode == 3 || rank == m)
            {
                goto L450;
            }
            i__1 = np1;
            for (i__ = np1mr; i__ <= i__1; ++i__)
            {
                k = (int)((Math.Abs(a[i__, mp1]) - (double)(n)));
                b[k] = a[np2, m] * r_sign(c_b44, a[i__, mp1]);
                /* L410: */
            }
            if (rankp1 == m)
            {
                goto L430;
            }
            i__1 = mm1;
            for (j = rankp1; j <= i__1; ++j)
            {
                k = (int)(Math.Abs(a[np3, j]) - (double)(n));
                b[k] = (a[np2, m] - a[np2, j]) * r_sign(c_b44, a[np3, j]);
                /* L420: */
            }
        /* TEST FOR NON-UNIQUE SOLUTION. */
        L430:
            i__1 = np1;
            for (i__ = np1mr; i__ <= i__1; ++i__)
            {
                if (Math.Abs(a[i__, m]) > tol)
                {
                    goto L440;
                }
                ocode = 0;
                goto L450;
            L440:
                ;
            }
        L450:
            if (mode != 2 && mode != 3)
            {
                resmax = a[np2, m];
            }
            if (rank == m)
            {
                resmax = 0.0;
            }
            if (mode == 4)
            {
                resmax -= d__;
            }
            return 0;
        } /* cheb_ */
        static public void test()
        {
            /* M      NUMBER OF EQUATIONS. */
            /* N      NUMBER OF UNKNOWNS (N MUST NOT EXCEED M). */
            /* MDIM   THE NUMBER OF COLUMNS OF A, AT LEAST M+1. */
            /* NDIM   THE NUMBER OF ROWS OF A, AT LEAST N+3. */
            /* A      TWO DIMENSIONAL double ARRAY OF SIZE (NDIM,MDIM). */
            /*        ON ENTRY,THE TRANSPOSE OF THE MATRIX OF */
            /*        COEFFICIENTS OF THE OVER-DETERMINED SYSTEM MUST */
            /*        BE STORED IN THE FIRST M COLUMNS AND N ROWS OF A. */
            /*        THESE VALUES ARE DESTROYED BY THE SUBROUTINE. */
            /* B      ONE DIMENSIONAL double ARRAY OF SIZE MDIM. ON ENTRY, */
            /*        B MUST CONTAIN THE RIGHT-HAND SIDES OF THE */
            /*        EQUATIONS IN ITS FIRST M LOCATIONS. ON EXIT, B */
            /*        CONTAINS THE RESIDUALS FOR THE EQUATIONS IN ITS */
            /*        FIRST M LOCATIONS (SEE DESCRIPTION). */
            /* TOL    A SMALL POSITIVE TOLERANCE. EMPIRICAL EVIDENCE */
            /*        SUGGESTS TOL=10**(-D+1) WHERE D REPRESENTS THE */
            /*        NUMBER OF DECIMAL DIGITS OF ACCURACY AVAILABLE */
            /*        (SEE DESCRIPTION). */
            /* RELERR A double VARIABLE WHICH ON ENTRY MUST HAVE THE VALUE */
            /*        0.0 IF A CHEBYSHEV SOLUTION IS REQUIRED. IF RELERR */
            /*        IS POSITIVE, THE SUBROUTINE CALCULATES AN */
            /*        APPROXIMATE SOLUTION WITH RELERR AS AN UPPER BOUND */
            /*        ON THE RELATIVE ERROR OF ITS LARGEST RESIDUAL (SEE */
            /*        DESCRIPTION-INEQUALITY (2)). ON EXIT, THE VALUE OF */
            /*        RELERR GIVES A SMALLER UPPER BOUND FOR THIS */
            /*        RELATIVE ERROR. */
            /* X      ONE DIMENSIONAL double ARRAY OF SIZE NDIM. ON EXIT, */
            /*        THIS ARRAY CONTAINS A SOLUTION TO THE PROBLEM IN */
            /*        ITS FIRST N LOCATIONS. */
            /* RANK   AN int WHICH GIVES ON EXIT THE RANK OF THE */
            /*        MATRIX OF COEFFICIENTS. */
            /* RESMAX THE LARGEST RESIDUAL IN MAGNITUDE. */
            /* ITER   THE NUMBER OF SIMPLEX ITERATIONS PERFORMED. */
            int m = 3;
            int n = 2;
            int mdim = m + 1 + 1;
            int ndim = n + 3 + 1;
            double[,] a = new double[,] { { 0, 0, 0, 0, 0 }, { 0, 1, 1, 1, 0 }, { 0, 1, 3, 4, 0 }, { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0 } };
            double[] b = new double[] { 0, 6, 7, 10, 0 };
            double[] x = new double[ndim];
            int rank = 0;
            double resmax = 0.0;
            int iter = 0;

            cheb_(m, n, mdim, ndim, a, b, 1e-15, 0.0, x, ref rank, ref  resmax, ref  iter);
        }
        static public double [] Solve( double []bx,int start=0)
        {
            int m = bx.Length;
            int n = 2;
            int mdim = m + 1 + 1;
            int ndim = n + 3 + 1;
            double[,] a = new double[ndim, mdim];
            double[] x = new double[ndim];
            double[] b = new double[mdim];
            int rank = 0;
            double resmax = 0.0;
            int iter = 0;           
            
             for (int j = 0; j < m; j++)
                {
                    a[0 + 1, j + 1] = j+start;
                    a[1 + 1, j + 1] = 1;
                }
            for (int i = 0; i < m; i++)
                b[i + 1] = bx[i];
            cheb_(m, n, mdim, ndim, a, b, 1e-15, 0.0, x, ref rank, ref  resmax, ref  iter);
            double[] bb = new double[n];
            for (int i = 0; i < n; i++)
                bb[i] = x[i + 1];
            return bb;
        }
        
    }

}
