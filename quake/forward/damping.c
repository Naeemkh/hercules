/* -*- C -*- */

/* @copyright_notice_start
 *
 * This file is part of the CMU Hercules ground motion simulator developed
 * by the CMU Quake project.
 *
 * Copyright (C) Carnegie Mellon University. All rights reserved.
 *
 * This program is covered by the terms described in the 'LICENSE.txt' file
 * included with this software package.
 *
 * This program comes WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * 'LICENSE.txt' file for more details.
 *
 *  @copyright_notice_end
 */

#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "psolve.h"
#include "quake_util.h"
#include "damping.h"
#include "stiffness.h"

void damping_addforce(mesh_t *myMesh, mysolver_t *mySolver, fmatrix_t (*theK1)[8], fmatrix_t (*theK2)[8]){

    fvector_t localForce[8];
    int       i,j;
    int32_t   eindex;

    fvector_t deltaDisp[8];

    /* loop on the number of elements */
    for (eindex = 0; eindex < myMesh->lenum; eindex++)
    {
        elem_t *elemp;
        e_t    *ep;

        elemp = &myMesh->elemTable[eindex];
        ep = &mySolver->eTable[eindex];

        /* Step 1. calculate the force due to the element stiffness */
        memset(localForce, 0, 8 * sizeof(fvector_t));

        /* compute the diff between disp(tm1) and disp(tm2) */

        /* the_E1_timer -= MPI_Wtime();*/
        for (i = 0; i < 8; i++) {
            fvector_t *tm1Disp, *tm2Disp;
            int32_t    lnid;

            lnid = elemp->lnid[i];

            tm1Disp = mySolver->tm1 + lnid;
            tm2Disp = mySolver->tm2 + lnid;

            deltaDisp[i].f[0] = tm1Disp->f[0] - tm2Disp->f[0];
            deltaDisp[i].f[1] = tm1Disp->f[1] - tm2Disp->f[1];
            deltaDisp[i].f[2] = tm1Disp->f[2] - tm2Disp->f[2];
        }

        if(vector_is_zero( deltaDisp ) != 0) {
            /*the_E3_timer += MPI_Wtime();  */

            for (i = 0; i < 8; i++)
            {
                fvector_t* toForce;
                toForce = &localForce[i];

                for (j = 0; j < 8; j++)
                {
                    fvector_t *myDeltaDisp;

                    /* contribution by ( - b * deltaT * Ke * ( Ut - Ut-1 ) ) */
                    /* But if myDeltaDisp is zero avoids multiplications     */
                    myDeltaDisp = &deltaDisp[j];

                    MultAddMatVec(&theK1[i][j], myDeltaDisp, -ep->c3, toForce);
                    MultAddMatVec(&theK2[i][j], myDeltaDisp, -ep->c4, toForce);

                }
            }
        }
        for (i = 0; i < 8; i++) {
            int32_t lnid;
            fvector_t *nodalForce;

            lnid = elemp->lnid[i];

            nodalForce = mySolver->force + lnid;
            nodalForce->f[0] += localForce[i].f[0];
            nodalForce->f[1] += localForce[i].f[1];
            nodalForce->f[2] += localForce[i].f[2];
        }

    } /* for all the elements */

    return;
}

/**
 * Introduce BKT Model: Compute and add the force due to the element
 *              damping.
 */

void calc_conv(mesh_t *myMesh, mysolver_t *mySolver, double theFreq, double theDeltaT, double theDeltaTSquared, damping_type_t typeOfDamping){

    int32_t eindex;
    int i;
    double cdt;

    if (typeOfDamping == BKT) {
        cdt = 2. * M_PI * theFreq * theDeltaT;
    } else {
    	cdt = theDeltaT;
    }

    for (eindex = 0; eindex < myMesh->lenum; eindex++)
    {
    	elem_t *elemp;
    	edata_t *edata;

    	elemp = &myMesh->elemTable[eindex];
    	edata = (edata_t *)elemp->data;

        double g0, g02, cg0, eg0;
        double g1, g12, cg1, eg1;
        double g2, g22, cg2, eg2;

        // SHEAR RELATED CONVOLUTION

    	if ( (edata->g0_shear != 0) && (edata->g1_shear != 0) ) {

            g0  = cdt * edata->g0_shear;
            g02 = g0 / 2.;
            cg0 = g02 * ( 1. - g0 );
            eg0 = exp( -g0 );

            g1  = cdt * edata->g1_shear;
            g12 = g1 / 2.;
            cg1 = g12 * ( 1. - g1 );
            eg1 = exp( -g1 );

            if (typeOfDamping >= BKT3) {
                g2  = cdt * edata->g2_shear;
                g22 = g2 / 2.;
                cg2 = g22 * ( 1. - g2 );
                eg2 = exp( -g2 );
            }

            for(i = 0; i < 8; i++)
            {
                int32_t     lnid, cindex;
                fvector_t  *tm1Disp, *tm2Disp;
                fvector_t  *f0_tm1, *f1_tm1, *f2_tm1;

                lnid = elemp->lnid[i];

                /* cindex is the index of the node in the convolution vector */
                cindex = eindex * 8 + i;

                tm1Disp = mySolver->tm1 + lnid;
                tm2Disp = mySolver->tm2 + lnid;

                f0_tm1 = mySolver->conv_shear_1 + cindex;
                f1_tm1 = mySolver->conv_shear_2 + cindex;

                f0_tm1->f[0] = cg0 * tm1Disp->f[0] + g02 * tm2Disp->f[0] + eg0 * f0_tm1->f[0];
                f0_tm1->f[1] = cg0 * tm1Disp->f[1] + g02 * tm2Disp->f[1] + eg0 * f0_tm1->f[1];
                f0_tm1->f[2] = cg0 * tm1Disp->f[2] + g02 * tm2Disp->f[2] + eg0 * f0_tm1->f[2];

                f1_tm1->f[0] = cg1 * tm1Disp->f[0] + g12 * tm2Disp->f[0] + eg1 * f1_tm1->f[0];
                f1_tm1->f[1] = cg1 * tm1Disp->f[1] + g12 * tm2Disp->f[1] + eg1 * f1_tm1->f[1];
                f1_tm1->f[2] = cg1 * tm1Disp->f[2] + g12 * tm2Disp->f[2] + eg1 * f1_tm1->f[2];

                if (typeOfDamping >= BKT3) {
                    f2_tm1 = mySolver->conv_shear_3 + cindex;
                    f2_tm1->f[0] = cg2 * tm1Disp->f[0] + g22 * tm2Disp->f[0] + eg2 * f1_tm1->f[0];
                    f2_tm1->f[1] = cg2 * tm1Disp->f[1] + g22 * tm2Disp->f[1] + eg2 * f1_tm1->f[1];
                    f2_tm1->f[2] = cg2 * tm1Disp->f[2] + g22 * tm2Disp->f[2] + eg2 * f1_tm1->f[2];
                }

            } // For local nodes (0:7)

    	} // end if null coefficients

        // DILATATION RELATED CONVOLUTION

        if ( (edata->g0_kappa != 0) && (edata->g1_kappa != 0) ) {

            g0  = cdt * edata->g0_kappa;
            g02 = g0 / 2.;
            cg0 = g02 * ( 1. - g0 );
            eg0 = exp( -g0 );

            g1  = cdt * edata->g1_kappa;
            g12 = g1 / 2.;
            cg1 = g12 * ( 1. - g1 );
            eg1 = exp( -g1 );

            if (typeOfDamping >= BKT3) {
                g2  = cdt * edata->g2_kappa;
                g22 = g2 / 2.;
                cg2 = g22 * ( 1. - g2 );
                eg2 = exp( -g2 );
            }

            for(i = 0; i < 8; i++)
            {
                int32_t     lnid, cindex;
                fvector_t  *tm1Disp, *tm2Disp;
                fvector_t  *f0_tm1, *f1_tm1, *f2_tm1;

                lnid = elemp->lnid[i];

                /* cindex is the index of the node in the convolution vector */
                cindex = eindex * 8 + i;

                tm1Disp = mySolver->tm1 + lnid;
                tm2Disp = mySolver->tm2 + lnid;

                f0_tm1 = mySolver->conv_kappa_1 + cindex;
                f1_tm1 = mySolver->conv_kappa_2 + cindex;

                f0_tm1->f[0] = cg0 * tm1Disp->f[0] + g02 * tm2Disp->f[0] + eg0 * f0_tm1->f[0];
                f0_tm1->f[1] = cg0 * tm1Disp->f[1] + g02 * tm2Disp->f[1] + eg0 * f0_tm1->f[1];
                f0_tm1->f[2] = cg0 * tm1Disp->f[2] + g02 * tm2Disp->f[2] + eg0 * f0_tm1->f[2];

                f1_tm1->f[0] = cg1 * tm1Disp->f[0] + g12 * tm2Disp->f[0] + eg1 * f1_tm1->f[0];
                f1_tm1->f[1] = cg1 * tm1Disp->f[1] + g12 * tm2Disp->f[1] + eg1 * f1_tm1->f[1];
                f1_tm1->f[2] = cg1 * tm1Disp->f[2] + g12 * tm2Disp->f[2] + eg1 * f1_tm1->f[2];

                if (typeOfDamping >= BKT3) {
                    f2_tm1 = mySolver->conv_kappa_3 + cindex;
                    f2_tm1->f[0] = cg2 * tm1Disp->f[0] + g22 * tm2Disp->f[0] + eg2 * f1_tm1->f[0];
                    f2_tm1->f[1] = cg2 * tm1Disp->f[1] + g22 * tm2Disp->f[1] + eg2 * f1_tm1->f[1];
                    f2_tm1->f[2] = cg2 * tm1Disp->f[2] + g22 * tm2Disp->f[2] + eg2 * f1_tm1->f[2];
                }

            } // For local nodes (0:7)

        } // end if null coefficients

    } // For all elements

    return;

}

/**
 * new_damping: Compute and add the force due to the element
 *              damping.
 */
void constant_Q_addforce(mesh_t *myMesh, mysolver_t *mySolver, double theFreq, double theDeltaT, double theDeltaTSquared, damping_type_t typeOfDamping)
{
	/* \todo use mu_and_lamda to compute first,second and third coefficients */

	int i;
	fvector_t localForce[8];
	int32_t   eindex;
	fvector_t damping_vector_shear[8], damping_vector_kappa[8];
    double rmax;

    if (typeOfDamping == BKT) {
        rmax = 2. * M_PI * theFreq * theDeltaT;
    } else {
    	rmax = theDeltaT;
    }

	/* theAddForceETime -= MPI_Wtime(); */

	/* loop on the number of elements */
	for (eindex = 0; eindex < myMesh->lenum; eindex++)
	{
		elem_t *elemp;
		e_t    *ep;
		edata_t *edata;

		double a0_shear, a1_shear, a2_shear,
		       a0_kappa, a1_kappa, a2_kappa,
		       b_shear, b_kappa,
		       csum;

		elemp = &myMesh->elemTable[eindex];
		edata = (edata_t *)elemp->data;
		ep = &mySolver->eTable[eindex];

		// SHEAR CONTRIBUTION

        a0_shear = edata->a0_shear;
        a1_shear = edata->a1_shear;
        a2_shear = edata->a2_shear;
        b_shear  = edata->b_shear;

        csum = a0_shear + a1_shear + b_shear;
        if ( typeOfDamping >= BKT3 ) {
            csum += a2_shear;
        }

		if ( csum != 0 ) {

	        double coef_shear = b_shear / rmax;

            for (i = 0; i < 8; i++) {

                fvector_t *tm1Disp, *tm2Disp, *f0_tm1, *f1_tm1;
                int32_t    lnid, cindex;

                cindex = eindex * 8 + i;
                lnid = elemp->lnid[i];

                tm1Disp = mySolver->tm1 + lnid;
                tm2Disp = mySolver->tm2 + lnid;
                f0_tm1  = mySolver->conv_shear_1 + cindex;
                f1_tm1  = mySolver->conv_shear_2 + cindex;

                damping_vector_shear[i].f[0] = coef_shear * (tm1Disp->f[0] - tm2Disp->f[0])
                                             - (a0_shear * f0_tm1->f[0] + a1_shear * f1_tm1->f[0])
                                             + tm1Disp->f[0];

                damping_vector_shear[i].f[1] = coef_shear * (tm1Disp->f[1] - tm2Disp->f[1])
                                             - (a0_shear * f0_tm1->f[1] + a1_shear * f1_tm1->f[1])
                                             + tm1Disp->f[1];

                damping_vector_shear[i].f[2] = coef_shear * (tm1Disp->f[2] - tm2Disp->f[2])
                                             - (a0_shear * f0_tm1->f[2] + a1_shear * f1_tm1->f[2])
                                             + tm1Disp->f[2];

                if ( typeOfDamping >= BKT3 ) {
                    damping_vector_shear[i].f[0] -= a2_shear * f0_tm1->f[0];
                    damping_vector_shear[i].f[1] -= a2_shear * f0_tm1->f[1];
                    damping_vector_shear[i].f[2] -= a2_shear * f0_tm1->f[2];
                }

            } // end for nodes in the element

		} else {

		    for (i = 0; i < 8; i++) {

                fvector_t *tm1Disp, *tm2Disp;
                int32_t    lnid;

                lnid = elemp->lnid[i];
                tm1Disp = mySolver->tm1 + lnid;
                tm2Disp = mySolver->tm2 + lnid;

                damping_vector_shear[i].f[0] = tm1Disp->f[0];
                damping_vector_shear[i].f[1] = tm1Disp->f[1];
                damping_vector_shear[i].f[2] = tm1Disp->f[2];

            } // end for nodes in the element

		} // end if for coefficients

		// DILATATION CONTRIBUTION

        a0_kappa   = edata->a0_kappa;
        a1_kappa   = edata->a1_kappa;
        a2_kappa   = edata->a2_kappa;
        b_kappa    = edata->b_kappa;

        csum = a0_kappa + a1_kappa + b_kappa;
        if ( typeOfDamping >= BKT3 ) {
            csum += a2_kappa;
        }

        if ( csum != 0 ) {

            double coef_kappa = b_kappa / rmax;

            for (i = 0; i < 8; i++) {

                fvector_t *tm1Disp, *tm2Disp, *f0_tm1, *f1_tm1;
                int32_t    lnid, cindex;

                cindex = eindex * 8 + i;

                lnid = elemp->lnid[i];

                tm1Disp = mySolver->tm1 + lnid;
                tm2Disp = mySolver->tm2 + lnid;

                f0_tm1  = mySolver->conv_kappa_1 + cindex;
                f1_tm1  = mySolver->conv_kappa_2 + cindex;

                damping_vector_kappa[i].f[0] = coef_kappa * (tm1Disp->f[0] - tm2Disp->f[0])
                                             - (a0_kappa * f0_tm1->f[0] + a1_kappa * f1_tm1->f[0])
                                             + tm1Disp->f[0];

                damping_vector_kappa[i].f[1] = coef_kappa * (tm1Disp->f[1] - tm2Disp->f[1])
                                             - (a0_kappa * f0_tm1->f[1] + a1_kappa * f1_tm1->f[1])
                                             + tm1Disp->f[1];

                damping_vector_kappa[i].f[2] = coef_kappa * (tm1Disp->f[2] - tm2Disp->f[2])
                                             - (a0_kappa * f0_tm1->f[2] + a1_kappa * f1_tm1->f[2])
                                             + tm1Disp->f[2];

                if ( typeOfDamping >= BKT3 ) {
                    damping_vector_kappa[i].f[0] -= a2_kappa * f0_tm1->f[0];
                    damping_vector_kappa[i].f[1] -= a2_kappa * f0_tm1->f[1];
                    damping_vector_kappa[i].f[2] -= a2_kappa * f0_tm1->f[2];
                }

            } // end for nodes in the element

        } else {

            for (i = 0; i < 8; i++) {

                fvector_t *tm1Disp, *tm2Disp;
                int32_t    lnid;

                lnid = elemp->lnid[i];
                tm1Disp = mySolver->tm1 + lnid;
                tm2Disp = mySolver->tm2 + lnid;

                damping_vector_kappa[i].f[0] = tm1Disp->f[0];
                damping_vector_kappa[i].f[1] = tm1Disp->f[1];
                damping_vector_kappa[i].f[2] = tm1Disp->f[2];

            } // end for nodes in the element

        } // end if for coefficients

		double kappa = -0.5625 * (ep->c2 + 2. / 3. * ep->c1);
		double mu = -0.5625 * ep->c1;

		double atu[24];
		double firstVec[24];

		for(i = 0; i<24; i++)
			firstVec[i] = 0.;

		memset(localForce, 0, 8 * sizeof(fvector_t));

        if(vector_is_zero( damping_vector_shear ) != 0) {

            aTransposeU( damping_vector_shear, atu );
            firstVector_mu( atu, firstVec, mu);

        }

		if(vector_is_zero( damping_vector_kappa ) != 0) {

			aTransposeU( damping_vector_kappa, atu );
			firstVector_kappa( atu, firstVec, kappa);

		}

		au( localForce, firstVec );

		for (i = 0; i < 8; i++) {
			int32_t lnid;
			fvector_t *nodalForce;

			lnid = elemp->lnid[i];

			nodalForce = mySolver->force + lnid;

			//Naeem: Commented three follow line

			nodalForce->f[0] += localForce[i].f[0];
			nodalForce->f[1] += localForce[i].f[1];
			nodalForce->f[2] += localForce[i].f[2];
		}

	} /* for all the elements */

	return;

}

/* setting damping parameters*/


void set_damping_params(edata_t *edata, double *theQTABLE, int QTable_Size,int  useInfQk,  double theFreq_Vel, double theFreq)
{

 	double   vs, vp, rho;
    double Qs, Qp, Qk, L, vs_vp_Ratio, vksquared, w;
    int index_Qs, index_Qk,i,k,j,l;
    double QTABLE[26][6];

    int array_size=(QTable_Size*6)+1;

    // regenertate the Q_table
    k=0;
    j=0;
	for(l = 1; l < array_size; l++)
	{
		i = l-1;
		QTABLE[k][j] = theQTABLE[i];
        if (l % 6 ==0){
        	k=k+1;
        }
        j=j+1;
        if (j==6){
        	j=0;
        }
	}



//    vksquared = edata->Vp * edata->Vp - 4. / 3. * edata->Vs * edata->Vs;
	vs_vp_Ratio = edata->Vs / edata->Vp;
	vs = edata->Vs;
	L = 4. / 3. * vs_vp_Ratio * vs_vp_Ratio;

	Qs = set_Qs(vs);
	Qp = set_Qp(Qs);

	if (useInfQk == YES) {
	    Qk = 1000;
	} else {
        Qk = (1. - L) / (1. / Qp - L / Qs);
	}


	index_Qs = Search_Quality_Table(Qs, &(QTABLE[0][0]), QTable_Size);
	index_Qk = Search_Quality_Table(Qk, &(QTABLE[0][0]), QTable_Size);
	update_Q_params(edata,Qs,Qk,theFreq);
    control_correction_factor(edata,theFreq_Vel,theFreq);


}


void update_Q_params(edata_t *edata,double Qs, double Qk, double theFreq){

    int i,k,j,l;
    //double myQTABLE[26][6];

    //int array_size=(QTable_Size*6)+1;

    // regenertate the Q_table
//    k=0;
//    j=0;
//	for(l = 1; l < array_size; l++)
//	{
//		i = l-1;
//		myQTABLE[k][j] = QTABLE[i];
//        if (l % 6 ==0){
//        	k=k+1;
//        }
//        j=j+1;
//        if (j==6){
//        	j=0;
//        }
//	}

//	if(index_Qs == -2 || index_Qs >= QTable_Size)
//		        	{
//		        		fprintf(stderr,"Problem with the Quality Factor Table\n Qs : %lf \n Vs : %lf\n",Qs,edata->Vs);
//		        		exit(1);
//		        	}
//		        	else if(index_Qs == -1)
//		        	{
//		        		edata->a0_shear = 0;
//		        		edata->a1_shear = 0;
//		        		edata->g0_shear = 0;
//		        		edata->g1_shear = 0;
//		        		edata->b_shear  = 0;
//		        	}
//		        	else
//		        	{
//		        		edata->a0_shear = myQTABLE[index_Qs][1];
//		        		edata->a1_shear = myQTABLE[index_Qs][2];
//		        		edata->g0_shear = myQTABLE[index_Qs][3];
//		        		edata->g1_shear = myQTABLE[index_Qs][4];
//		        		edata->b_shear  = myQTABLE[index_Qs][5];
//		        	}
//
//	if (index_Qk != 0){
//	if(index_Qk == -2 || index_Qk >= QTable_Size)
//	        	{
//	        		fprintf(stderr,"Problem with the Quality Factor Table\n Qk : %lf \n Vs : %lf\n",Qk,edata->Vs);
//	        		exit(1);
//	        	}
//	        	else if(index_Qk == -1)
//	        	{
//	        		edata->a0_kappa = 0;
//	        		edata->a1_kappa = 0;
//	        		edata->g0_kappa = 0;
//	        		edata->g1_kappa = 0;
//	        		edata->b_kappa  = 0;
//	        	}
//	        	else
//	        	{
//	        		edata->a0_kappa = myQTABLE[index_Qk][1];
//	        		edata->a1_kappa = myQTABLE[index_Qk][2];
//	        		edata->g0_kappa = myQTABLE[index_Qk][3];
//	        		edata->g1_kappa = myQTABLE[index_Qk][4];
//	        		edata->b_kappa  = myQTABLE[index_Qk][5];
//	        	}
//	}


    if ( Qs >= 1000 ) {
        edata->g0_shear = 0;
        edata->g1_shear = 0;
        edata->a0_shear = 0;
        edata->a1_shear = 0;
        edata->b_shear  = 0;
    } else {
        edata->g0_shear =   0.0373 * (2. * M_PI * theFreq);
        edata->g1_shear =   0.3082 * (2. * M_PI * theFreq);
        edata->a0_shear = (-2.656  * pow(Qs, -0.8788) + 1.677 ) / Qs;
        edata->a1_shear = (-0.5623 * pow(Qs, -1.0300) + 1.262 ) / Qs;
        edata->b_shear  = ( 0.1876 * pow(Qs, -0.9196) + 0.6137) / (Qs * (2. * M_PI * theFreq));
    }

//    if ( ( Param.useInfQk == YES ) || ( Qk >= 1000 ) ) {
//        edata->g0_kappa = 0;
//        edata->g1_kappa = 0;
//        edata->a0_kappa = 0;
//        edata->a1_kappa = 0;
//        edata->b_kappa  = 0;
//    } else {
//        edata->g0_kappa =   0.0373 * (2. * M_PI * Param.theFreq);
//        edata->g1_kappa =   0.3082 * (2. * M_PI * Param.theFreq);
//        edata->a0_kappa = (-2.656  * pow(Qk, -0.8788) + 1.677 ) / Qk;
//        edata->a1_kappa = (-0.5623 * pow(Qk, -1.0300) + 1.262 ) / Qk;
//        edata->b_kappa  = ( 0.1876 * pow(Qk, -0.9196) + 0.6137) / (Qk * (2. * M_PI * Param.theFreq));
//    }



}

void control_correction_factor(edata_t *edata, double theFreq_Vel,double theFreq){

	double vksquared = edata->Vp * edata->Vp - 4. / 3. * edata->Vs * edata->Vs;
	double w;

	if(theFreq_Vel != 0.)
	{
		w = theFreq_Vel / theFreq;

		if ( (edata->a0_shear != 0) && (edata->a1_shear != 0) ) {
		    double shear_vel_corr_factor;
		    shear_vel_corr_factor = sqrt(1. - (edata->a0_shear * edata->g0_shear * edata->g0_shear / (edata->g0_shear * edata->g0_shear + w * w) + edata->a1_shear * edata->g1_shear * edata->g1_shear / (edata->g1_shear * edata->g1_shear + w * w)));
            edata->Vs = shear_vel_corr_factor * edata->Vs;
		}

		if ( (edata->a0_kappa != 0) && (edata->a0_kappa != 0) ) {
		    double kappa_vel_corr_factor;
		    kappa_vel_corr_factor = sqrt(1. - (edata->a0_kappa * edata->g0_kappa * edata->g0_kappa / (edata->g0_kappa * edata->g0_kappa + w * w) + edata->a1_kappa * edata->g1_kappa * edata->g1_kappa / (edata->g1_kappa * edata->g1_kappa + w * w)));
            edata->Vp = sqrt(kappa_vel_corr_factor * kappa_vel_corr_factor * vksquared + 4. / 3. * edata->Vs * edata->Vs);
		}
	}


}

