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
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <gsl/gsl_poly.h>

#include "geometrics.h"
//#include "nonlinear.h"
#include "eqlinear.h"
#include "octor.h"
#include "psolve.h"
#include "quake_util.h"
#include "util.h"
#include "damping.h"

#define  QC  qc = 0.577350269189 /* sqrt(3.0)/3.0; */

#define MAX(a, b) ((a)>(b)?(a):(b))

#define  XI  xi[3][8] = { {-1,  1, -1,  1, -1,  1, -1, 1} , \
                          {-1, -1,  1,  1, -1, -1,  1, 1} , \
                          {-1, -1, -1, -1,  1,  1,  1, 1} }


/* -------------------------------------------------------------------------- */
/*                             Global Variables                               */
/* -------------------------------------------------------------------------- */
static int32_t              *myStationsElementIndices;
static elsolver_t           *myEqlinSolver;
static int32_t              *myEqlinStationsMapping;
static int32_t               myNumberOfEqlinStations;
static int32_t               myEqlinElementsCount = 0;
static int32_t              *myEqlinElementsMapping;
static double                theGDTABLE[11][3];




/* -------------------------------------------------------------------------- */
/*                                 Utilities                                  */
/* -------------------------------------------------------------------------- */


int isThisElementEqLinear(mesh_t *myMesh, int32_t eindex) {

	elem_t  *elemp;
	edata_t *edata;



	return YES;
}






/* -------------------------------------------------------------------------- */
/*       Initialization of parameters, structures and memory allocations      */
/* -------------------------------------------------------------------------- */

/*
 * Counts the number of nonlinear elements in my local mesh
 */
void eqlinear_elements_count(int32_t myID, mesh_t *myMesh) {

    int32_t eindex;
    int32_t count = 0;

    for (eindex = 0; eindex < myMesh->lenum; eindex++) {


        if ( isThisElementEqLinear(myMesh, eindex) == YES ) {
            count++;
        }
    }

    if ( count > myMesh-> lenum ) {
        fprintf(stderr,"Thread %d: nl_elements_count: "
                "more elements than expected\n", myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    myEqlinElementsCount = count;

    return;
}


void eqlinear_elements_mapping(int32_t myID, mesh_t *myMesh) {

    int32_t eindex;
    int32_t count = 0;

    XMALLOC_VAR_N(myEqlinElementsMapping, int32_t, myEqlinElementsCount);

    for (eindex = 0; eindex < myMesh->lenum; eindex++) {

        if ( isThisElementEqLinear(myMesh, eindex) == YES ) {
            myEqlinElementsMapping[count] = eindex;
            count++;
        }
    }

    if ( count != myEqlinElementsCount ) {
        fprintf(stderr,"Thread %d: el_elements_mapping: "
                "more elements than the count\n", myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    return;
}




void eqlinear_solver_init(int32_t myID, mesh_t *myMesh, double depth) {

      int32_t eindex, el_eindex;   // el_eindex = equivalent linear element index
//
      eqlinear_elements_count(myID, myMesh);
      eqlinear_elements_mapping(myID, myMesh);


//
//    if ( theGeostaticLoadingT > 0 ) {
//        bottom_elements_count(myID, myMesh, depth);
//        bottom_elements_mapping(myID, myMesh, depth);
//    }

    /* Memory allocation for mother structure */



      myEqlinSolver = (elsolver_t *)malloc(sizeof(elsolver_t));

    if (myEqlinSolver == NULL) {
        fprintf(stderr, "Thread %d: nonlinear_init: out of memory\n", myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    /* Memory allocation for internal structures */

    myEqlinSolver->constants =
        (elconstants_t *)calloc(myEqlinElementsCount, sizeof(elconstants_t));
    myEqlinSolver->stresses =
        (eq_qptensors_t *)calloc(myEqlinElementsCount, sizeof(eq_qptensors_t));
    myEqlinSolver->strains =
        (eq_qptensors_t *)calloc(myEqlinElementsCount, sizeof(eq_qptensors_t));
    myEqlinSolver->maxstrains =
        (eq_qptensors_t *)calloc(myEqlinElementsCount, sizeof(eq_qptensors_t));
//    myEqlinSolver->pstrains1 =
//        (qptensors_t *)calloc(myEqlinElementsCount, sizeof(qptensors_t));
//    myEqlinSolver->pstrains2 =
//        (qptensors_t *)calloc(myEqlinElementsCount, sizeof(qptensors_t));
//    myEqlinSolver->alphastress1 =
//        (qptensors_t *)calloc(myEqlinElementsCount, sizeof(qptensors_t));
//    myEqlinSolver->alphastress2 =
//        (qptensors_t *)calloc(myEqlinElementsCount, sizeof(qptensors_t));
//    myEqlinSolver->ep1 =
//        (qpvectors_t *)calloc(myEqlinElementsCount, sizeof(qpvectors_t));
//    myEqlinSolver->ep2 =
//        (qpvectors_t *)calloc(myEqlinElementsCount, sizeof(qpvectors_t));

    if ( (myEqlinSolver->constants           == NULL) ||
         (myEqlinSolver->stresses            == NULL) ||
         (myEqlinSolver->strains             == NULL) ||
		 (myEqlinSolver->maxstrains          == NULL)
//         (myEqlinSolver->ep1                 == NULL) ||
//         (myEqlinSolver->ep2                 == NULL) ||
//         (myEqlinSolver->alphastress1        == NULL) ||
//         (myEqlinSolver->alphastress2        == NULL) ||
//         (myEqlinSolver->pstrains1           == NULL) ||
//         (myEqlinSolver->pstrains2           == NULL)
		 ) {

        fprintf(stderr, "Thread %d: eqlinear_init: out of memory\n", myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    /* Initialization of element constants
     * Tensors have been initialized to 0 by calloc
     */

    for (el_eindex = 0; el_eindex < myEqlinElementsCount; el_eindex++) {

        elem_t     *elemp;
        edata_t    *edata;
        elconstants_t *ecp;
        double      mu, lambda;
        double      elementVs, elementVp;

        eindex = myEqlinElementsMapping[el_eindex];

        elemp = &myMesh->elemTable[eindex];
        edata = (edata_t *)elemp->data;
        ecp   = myEqlinSolver->constants + el_eindex;

        int32_t lnid0 = elemp->lnid[0];
        double  zo    = myMesh->ticksize * myMesh->nodeTable[lnid0].z;

        /* get element Vs */

        elementVs   = (double)edata->Vs;
        elementVp   = (double)edata->Vp;



        /* Calculate the lame constants and store in element */

        mu_and_lambda(&mu, &lambda, edata, eindex);
        ecp->lambda = lambda;
        ecp->mu     = mu;

        /* Calculate the vertical stress as a homogeneous half-space */
//        if ( theApproxGeoState == YES )
//        	ecp->sigmaZ_st = edata->rho * 9.80 * ( zo + edata->edgesize / 2.0 );

        /* Calculate the yield function constants */
//        switch (theMaterialModel) {
//
//            case LINEAR:
//                ecp->alpha    = 0.;
//                ecp->k        = 0.;
//                ecp->phi      = 0.;
//                ecp->beta     = 0.;
//                ecp->h        = 0.;
//                ecp->Sstrain0 = 0.;
//                break;
//
//            case VONMISES:
//            	ecp->c     = get_cohesion(elementVs);
//            	ecp->phi   = get_phi(elementVs);
//            	ecp->dil_angle = 0.0;
//
//            	ecp->alpha = 0.;
//            	ecp->beta  = 0.;
//            	ecp->gamma = 1.0;
//
//            	ecp->Sstrain0    = interpolate_property_value(elementVs, theGamma0);
//
//            	if (ecp->Sstrain0 == 0)
//            		ecp->k     = ecp->c;
//            	else
//            		ecp->k     = ecp->Sstrain0 * ecp->mu;
//
//            	ecp->h         = 0; /*  no isotropic hardening  in von Mises model */
//            	break;
//            case DRUCKERPRAGER:
//                ecp->c         = get_cohesion(elementVs);
//                ecp->phi       = get_phi(elementVs);
//                ecp->dil_angle = get_dilatancy(elementVs);
//
//                ecp->alpha = get_alpha(elementVs, ecp->phi);
//                ecp->beta  = get_beta(elementVs);
//                ecp->gamma = get_gamma(elementVs,ecp->phi);
//
//                ecp->k         = ecp->gamma * ecp->c;
//                ecp->h         = get_hardmod(elementVs);
//                ecp->Sstrain0  = 0.0;
//
//            	break;
//            case MOHR_COULOMB:
//                ecp->c     = get_cohesion(elementVs);
//                ecp->phi   = get_phi(elementVs);
//                ecp->dil_angle = get_dilatancy(elementVs);
//
//                ecp->alpha = 0.;
//                ecp->beta  = 0.;
//                ecp->gamma = 0.;
//
//                ecp->h         = get_hardmod(elementVs);
//                ecp->Sstrain0  = 0.0;
//            	break;
//
//            default:
//                fprintf(stderr, "Thread %d: nonlinear_solver_init:\n"
//                        "\tUnexpected error with the material model\n", myID);
//                MPI_Abort(MPI_COMM_WORLD, ERROR);
//                exit(1);
//                break;
//        }


//        ecp->strainrate  =
//        		interpolate_property_value(elementVs, theStrainRates  );
//        ecp->sensitivity =
//        		interpolate_property_value(elementVs, theSensitivities );
////        ecp->hardmodulus =
////        		interpolate_property_value(elementVs, theHardeningModulus );
//
////        if ( theApproxGeoState == NO ) {
////            ecp->I1_st        =  -1.;
////            ecp->J2square_st  =   0.;
////        } else {
////            ecp->I1_st        = -S_zz * ( 3.0 - 2.0 * sin ( ecp->phi ) );
////            ecp->J2square_st  =  S_zz * S_zz * sin ( ecp->phi ) * sin ( ecp->phi ) / 3.0;
////        }
//
//
   } /* for all elements */

}



void constract_GD_Table()
{
	// Construct Shear modulus degredation and damping table

		int i,j;
		double local_GDtable[11][3] = {{ 0.0001,	1.000, 0.24},
				                       { 0.0003,	1.000, 0.42},
				                       { 0.0010,	1.000, 0.80},
				                       { 0.0030,	0.981, 1.40},
				                       { 0.0100,	0.941, 2.80},
			                           { 0.0300,	0.847, 5.10},
				                       { 0.1000,	0.656, 9.80},
				                       { 0.3000,	0.438, 15.50},
				                       { 1.0000,	0.238, 21.00},
				                       { 3.0000,	0.144, 25.00},
				                       { 10.0000,	0.110, 28.00},
		};

		for(i = 0; i < 11; i++)
		{
			for(j = 0; j < 2; j++)
			{
				theGDTABLE[i][j] = local_GDtable[i][j];
		//		printf("%f ",theQTABLE[i][j]);
			}
		//	printf("\n");
		}
return;
}


void eqlinear_stats(int32_t myID, int32_t theGroupSize) {

//    int32_t *nonlinElementsCount = NULL;
//    int32_t *nonlinStationsCount = NULL;
//    int32_t *bottomElementsCount = NULL;
//
//    if ( myID == 0 ) {
//        XMALLOC_VAR_N( nonlinElementsCount, int32_t, theGroupSize);
//        XMALLOC_VAR_N( nonlinStationsCount, int32_t, theGroupSize);
//        XMALLOC_VAR_N( bottomElementsCount, int32_t, theGroupSize);
//    }
//
//    MPI_Gather( &myNonlinElementsCount,    1, MPI_INT,
//                nonlinElementsCount,       1, MPI_INT, 0, comm_solver );
//    MPI_Gather( &myNumberOfNonlinStations, 1, MPI_INT,
//                nonlinStationsCount,       1, MPI_INT, 0, comm_solver );
//    MPI_Gather( &myBottomElementsCount,    1, MPI_INT,
//                bottomElementsCount,       1, MPI_INT, 0, comm_solver );
//
//    if ( myID == 0 ) {
//
//        nonlinear_print_stats( nonlinElementsCount, nonlinStationsCount,
//                               bottomElementsCount, theGroupSize);
//
//        xfree_int32_t( &nonlinElementsCount );
//    }

    return;
}

void eqlinear_init( int32_t     myID,
                     const char *parametersin,
                     double      theDeltaT,
                     double      theEndT )
{


//    double  double_message[2];
//    int     int_message[7];
//
//    /* Capturing data from file --- only done by PE0 */
//    if (myID == 0) {
//        if (nonlinear_initparameters(parametersin, theDeltaT, theEndT) != 0) {
//            fprintf(stderr,"Thread 0: nonlinear_local_init: "
//                    "nonlinear_initparameters error\n");
//            MPI_Abort(MPI_COMM_WORLD, ERROR);
//            exit(1);
//        }
//    }
//
//    /* Broadcasting data */
//    double_message[0] = theGeostaticLoadingT;
//    double_message[1] = theGeostaticCushionT;
//
//    int_message[0] = (int)theMaterialModel;
//    int_message[1] = thePropertiesCount;
//    int_message[2] = theGeostaticFinalStep;
//    int_message[3] = (int)thePlasticityModel;
//    int_message[4] = (int)theApproxGeoState;
//    int_message[5] = (int)theNonlinearFlag;
//    int_message[6] = (int)theTensionCutoff;
//
//    MPI_Bcast(double_message, 2, MPI_DOUBLE, 0, comm_solver);
//    MPI_Bcast(int_message,    7, MPI_INT,    0, comm_solver);
//
//    theGeostaticLoadingT  = double_message[0];
//    theGeostaticCushionT  = double_message[1];
//
//    theMaterialModel      = int_message[0];
//    thePropertiesCount    = int_message[1];
//    theGeostaticFinalStep = int_message[2];
//    thePlasticityModel    = int_message[3];
//    theApproxGeoState     = int_message[4];
//    theNonlinearFlag      = int_message[5];
//    theTensionCutoff      = int_message[6];
//
//    /* allocate table of properties for all other PEs */
//
//    if (myID != 0) {
//        theVsLimits         = (double*)malloc(sizeof(double) * thePropertiesCount);
//        theAlphaCohes       = (double*)malloc(sizeof(double) * thePropertiesCount);
//        theKayPhis          = (double*)malloc(sizeof(double) * thePropertiesCount);
//        theStrainRates      = (double*)malloc(sizeof(double) * thePropertiesCount);
//        theSensitivities    = (double*)malloc(sizeof(double) * thePropertiesCount);
//        theHardeningModulus = (double*)malloc(sizeof(double) * thePropertiesCount);
//        theBetaDilatancy    = (double*)malloc(sizeof(double) * thePropertiesCount);
//        theGamma0           = (double*)malloc(sizeof(double) * thePropertiesCount);
//    }
//
//    /* Broadcast table of properties */
//    MPI_Bcast(theVsLimits,         thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
//    MPI_Bcast(theAlphaCohes,       thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
//    MPI_Bcast(theKayPhis,          thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
//    MPI_Bcast(theStrainRates,      thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
//    MPI_Bcast(theSensitivities,    thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
//    MPI_Bcast(theHardeningModulus, thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
//    MPI_Bcast(theBetaDilatancy,    thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
//    MPI_Bcast(theGamma0,           thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
}




/* -------------------------------------------------------------------------- */
/*                   Auxiliary tensor manipulation methods                    */
/* -------------------------------------------------------------------------- */
eq_tensor_t point_strain_eq (fvector_t *u, double lx, double ly, double lz, double h) {

    int i;

    eq_tensor_t strain = init_tensor_eq();

    /* Contribution of each node */
    for (i = 0; i < 8; i++) {

        double dx, dy, dz;

        point_dxi(&dx, &dy, &dz, lx, ly, lz, h, i);

        strain.xx += dx * u[i].f[0];
        strain.yy += dy * u[i].f[1];
        strain.zz += dz * u[i].f[2];

        strain.xy += 0.5 * ( dy * u[i].f[0] + dx * u[i].f[1] );
        strain.yz += 0.5 * ( dz * u[i].f[1] + dy * u[i].f[2] );
        strain.xz += 0.5 * ( dz * u[i].f[0] + dx * u[i].f[2] );

    } /* nodes contribution */

    return strain;
}

/*
 * Resets a tensor to zero in all its components.
 */
eq_tensor_t init_tensor_eq() {

    eq_tensor_t tensor;

    tensor.xx = 0.0;
    tensor.yy = 0.0;
    tensor.zz = 0.0;
    tensor.xy = 0.0;
    tensor.yz = 0.0;
    tensor.xz = 0.0;

    return tensor;
}



/* -------------------------------------------------------------------------- */
/*                              Stability methods                             */
/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */
/*                   Nonlinear core computational methods                     */
/* -------------------------------------------------------------------------- */


GD_t  search_GD_table(double strain){

      GD_t GD;
      int  table_r;

	  double thGDtable[11][3] = {{ 0.0001,	1.000, 0.24},
				                       { 0.0003,	1.000, 0.42},
				                       { 0.0010,	1.000, 0.80},
				                       { 0.0030,	0.981, 1.40},
				                       { 0.0100,	0.941, 2.80},
			                           { 0.0300,	0.847, 5.10},
				                       { 0.1000,	0.656, 9.80},
				                       { 0.3000,	0.438, 15.50},
				                       { 1.0000,	0.238, 21.00},
				                       { 3.0000,	0.144, 25.00},
				                       { 10.0000,	0.110, 28.00},
		};


	  int GDTable_Size = (int)(sizeof(thGDtable)/( 3 * sizeof(double)));

	  if (strain <= thGDtable[0][0]){
		  GD.g = thGDtable[0][1];
		  GD.d = thGDtable[0][2];
	  } else if (strain >= thGDtable[GDTable_Size-1][0]) {
		  GD.g = thGDtable[GDTable_Size-1][1];
		  GD.d = thGDtable[GDTable_Size-1][2];
	  } else {

		  for  (table_r = 0; table_r < GDTable_Size; table_r++) {

			  if (strain >= thGDtable[table_r][0] && strain < thGDtable[table_r+1][0]){
				  GD.g = thGDtable[table_r+1][1] - ((thGDtable[table_r+1][0] - strain)*(thGDtable[table_r+1][1]-thGDtable[table_r][1])/(thGDtable[table_r+1][0]-thGDtable[table_r][0]));
				  GD.d = thGDtable[table_r+1][2] - ((thGDtable[table_r+1][0] - strain)*(thGDtable[table_r+1][2]-thGDtable[table_r][2])/(thGDtable[table_r+1][0]-thGDtable[table_r][0]));
				  return GD;
			  }

		  }


	  }

	  return GD;

}

void compute_eqlinear_state ( mesh_t     *myMesh,
                               mysolver_t *mySolver,
                               int32_t     theNumberOfStations,
                               int32_t     myNumberOfStations,
                               station_t  *myStations,
                               double      theDeltaT,
                               int         step )
{
	/* In general, j-index refers to the quadrature point in a loop (0 to 7 for
	 * eight points), and i-index refers to the tensor component (0 to 5), with
	 * the following order xx[0], yy[1], zz[2], xy[3], yz[4], xz[5]. i-index is
	 * also some times used for the number of nodes (8, 0 to 7).
	 */

	int     i;
	int32_t eindex, el_eindex;


	/* Loop over the number of local elements */
	for (el_eindex = 0; el_eindex < myEqlinElementsCount; el_eindex++) {

		elem_t        *elemp;
		edata_t       *edata;
		elconstants_t *enlcons;

		double         h;          /* Element edge-size in meters   */
// 		double         alpha, k;   /* Drucker-Prager constants      */
		double         mu, lambda; /* Elasticity material constants */
//		double		   hrd;        /* Hardening Modulus  */
//		double         beta;       /* Plastic flow rule constant */
		double         XI, QC;
		fvector_t      u[8];
		eq_qptensors_t   *stresses, *tstrains, *maxstrains; //, *pstrains1, *pstrains2, *alphastress1, *alphastress2;
//		qpvectors_t   *epstr1, *epstr2;

		/* Capture data from the element and mesh */

		eindex = myEqlinElementsMapping[el_eindex];

		elemp = &myMesh->elemTable[eindex];
		edata = (edata_t *)elemp->data;
		h     = edata->edgesize;

		/* Capture data from the nonlinear element structure */

		enlcons = myEqlinSolver->constants + el_eindex;

		mu     = enlcons->mu;
		lambda = enlcons->lambda;
//		alpha  = enlcons->alpha;
//		beta   = enlcons->beta;
//		k      = enlcons->k;
//		hrd    = enlcons->h;

		/* Capture the current state in the element */
		tstrains     = myEqlinSolver->strains      + el_eindex;
		stresses     = myEqlinSolver->stresses     + el_eindex;
		maxstrains   = myEqlinSolver->maxstrains   + el_eindex;
//		pstrains1    = myNonlinSolver->pstrains1    + nl_eindex;   /* Previous plastic tensor  */
//		pstrains2    = myNonlinSolver->pstrains2    + nl_eindex;   /* Current  plastic tensor  */
//		alphastress1 = myNonlinSolver->alphastress1 + nl_eindex;   /* Previous backstress tensor  */
//		alphastress2 = myNonlinSolver->alphastress2 + nl_eindex;   /* Current  backstress tensor  */
//		epstr1       = myNonlinSolver->ep1          + nl_eindex;
//		epstr2       = myNonlinSolver->ep2          + nl_eindex;

		// temp code to control (naeem)
		// printf("This is strain: %f \n", tstrains[1]);

		if ( get_displacements(mySolver, elemp, u) == 0 ) {
			/* If all displacements are zero go for next element */
			continue;
		}

		/* Loop over the quadrature points */
		for (i = 0; i < 8; i++) {

			eq_tensor_t  sigma0;

			/* Quadrature point local coordinates */
			double lx = xi[0][i] * qc ;
			double ly = xi[1][i] * qc ;
			double lz = xi[2][i] * qc ;

			/* Calculate total strains */
			tstrains->qp[i] = point_strain_eq(u, lx, ly, lz, h);

			//For each element you need to define another strain tensor, to only keep the maximum value of strain.

            maxstrains->qp[i].xx = MAX(fabs(maxstrains->qp[i].xx),fabs(tstrains->qp[i].xx));
            maxstrains->qp[i].yy = MAX(fabs(maxstrains->qp[i].yy),fabs(tstrains->qp[i].yy));
            maxstrains->qp[i].zz = MAX(fabs(maxstrains->qp[i].zz),fabs(tstrains->qp[i].zz));
            maxstrains->qp[i].xy = MAX(fabs(maxstrains->qp[i].xy),fabs(tstrains->qp[i].xy));
            maxstrains->qp[i].yz = MAX(fabs(maxstrains->qp[i].yz),fabs(tstrains->qp[i].yz));
            maxstrains->qp[i].xz = MAX(fabs(maxstrains->qp[i].xz),fabs(tstrains->qp[i].xz));


//            printf("This is xy max strain: %.20f (node %i) \n", maxstrains->qp[i].xy, i);


			//printf("This is strain: %f \n", tstrains[1])
			/* strain and backstress predictor  */
//	        pstrains1->qp[i]    = copy_tensor ( pstrains2->qp[i] );     /* The strain predictor assumes that the current plastic strains equal those from the previous step   */
//	        alphastress1->qp[i] = copy_tensor ( alphastress2->qp[i] );

			/* Calculate stresses */
//			if ( ( theMaterialModel == LINEAR ) || ( step <= theGeostaticFinalStep ) ){
//				stresses->qp[i]  = point_stress ( tstrains->qp[i], mu, lambda );
//				continue;
//			} else {

//				if ( theApproxGeoState == YES )
//					sigma0 = ApproxGravity_tensor(enlcons->sigmaZ_st, enlcons->phi, h, lz, edata->rho);
//				else
//					sigma0 = zero_tensor();

//				material_update ( *enlcons,           tstrains->qp[i],      pstrains1->qp[i], alphastress1->qp[i], epstr1->qv[i], sigma0, theDeltaT,
//						           &pstrains2->qp[i], &alphastress2->qp[i], &stresses->qp[i], &epstr2->qv[i],      &enlcons->fs[i]);




//			}
		} /* for all quadrature points */
	} /* for all nonlinear elements */
}



void material_update_eq ( mesh_t     *myMesh,
                               mysolver_t *mySolver,
                               int32_t     theNumberOfStations,
                               int32_t     myNumberOfStations,
                               station_t  *myStations,
                               double      theDeltaT,
                               int         eq_it,
							   double      theBBase,
							   double      theThresholdVpVs,
							   double      *theQTABLE,
							   int         QTable_Size,
							   double      theFreq_Vel,
							   double      theFreq
							   )
{
	/* In general, j-index refers to the quadrature point in a loop (0 to 7 for
	 * eight points), and i-index refers to the tensor component (0 to 5), with
	 * the following order xx[0], yy[1], zz[2], xy[3], yz[4], xz[5]. i-index is
	 * also some times used for the number of nodes (8, 0 to 7).
	 */

	int     i;
	int32_t eindex, el_eindex;


	/* Loop over the number of local elements */
	for (el_eindex = 0; el_eindex < myEqlinElementsCount; el_eindex++) {

		elem_t        *elemp;
		edata_t       *edata;
		elconstants_t *enlcons;
		e_t           *ep;    /* pointer to the element constant table */

		double         h;          /* Element edge-size in meters   */
		double         mu, lambda; /* Elasticity material constants */
		double         XI, QC;
		fvector_t      u[8];
		eq_qptensors_t   *maxstrains;
		double         zeta, a, b, updated_mu, updated_lambda;
		double         requested_Q, new_Q;

		/* Capture data from the element and mesh */
		eindex = myEqlinElementsMapping[el_eindex];

		elemp = &myMesh->elemTable[eindex];
		edata = (edata_t *)elemp->data;
		h     = edata->edgesize;
		ep    = &mySolver->eTable[eindex];

		/* Capture data from the eqlinear element structure */
		enlcons = myEqlinSolver->constants + el_eindex;

		mu     = enlcons->mu;
		lambda = enlcons->lambda;


		/* Capture the maximum strain of the element */

		maxstrains   = myEqlinSolver->maxstrains   + el_eindex;
//		pstrains1    = myNonlinSolver->pstrains1    + nl_eindex;   /* Previous plastic tensor  */
//		pstrains2    = myNonlinSolver->pstrains2    + nl_eindex;   /* Current  plastic tensor  */
//		alphastress1 = myNonlinSolver->alphastress1 + nl_eindex;   /* Previous backstress tensor  */
//		alphastress2 = myNonlinSolver->alphastress2 + nl_eindex;   /* Current  backstress tensor  */
//		epstr1       = myNonlinSolver->ep1          + nl_eindex;
//		epstr2       = myNonlinSolver->ep2          + nl_eindex;
//
//		 temp code to control (naeem)
//		 printf("This is strain: %f \n", tstrains[1]);
//
//		if ( get_displacements(mySolver, elemp, u) == 0 ) {
//			/* If all displacements are zero go for next element */
//			continue;
//		}
//         define a temprory matrix for strain

		double  strain_mat[8][6]={{0,  0,  0,  0,  0,  0},
								  {0,  0,  0,  0,  0,  0},
								  {0,  0,  0,  0,  0,  0},
								  {0,  0,  0,  0,  0,  0},
								  {0,  0,  0,  0,  0,  0},
								  {0,  0,  0,  0,  0,  0},
								  {0,  0,  0,  0,  0,  0},
								  {0,  0,  0,  0,  0,  0}};

		/* Loop over the quadrature points */
		for (i = 0; i < 8; i++) {


			strain_mat[i][0]=maxstrains->qp[i].xx;
			strain_mat[i][1]=maxstrains->qp[i].yy;
			strain_mat[i][2]=maxstrains->qp[i].zz;
			strain_mat[i][3]=maxstrains->qp[i].xy;
			strain_mat[i][4]=maxstrains->qp[i].yz;
			strain_mat[i][5]=maxstrains->qp[i].xz;

		} /* for all quadrature points */

		//printf("Element %i - XX strain is: %.20f \n", el_eindex, strain_mat[0][0]);

		// Now update the material based on strain level

		double myteststrain = 0;

		for (i=0; i < 8; i++){
		myteststrain = strain_mat[i][3] + myteststrain;
		}

		// mean of all xy strains.
		myteststrain = 100* (myteststrain/8);

//		printf("Here is the max strain of element %i : %.10f \n", el_eindex, myteststrain);

		 GD_t GD = search_GD_table(myteststrain);

		//printf("Results : strain = %f , G = %f, D = %f \n", myteststrain, GD.g, GD.d);


		 updated_mu = mu * GD.g;

         /* control the poisson ratio and lambda value */


         if ( edata->Vp > (edata->Vs * theThresholdVpVs) ) {
             updated_lambda = edata->rho * edata->Vs * edata->Vs * theThresholdVpVs
                    * theThresholdVpVs - 2 * updated_mu;
         } else {
        	 updated_lambda = edata->rho * edata->Vp * edata->Vp - 2 * updated_mu;
         }

         /* Adjust Vs, Vp to fix Poisson ratio problem, formula provided by Jacobo */
         if ( updated_lambda < 0 ) {
             if ( edata->Vs < 500 )
                 edata->Vp = 2.45 * edata->Vs;
             else if ( edata->Vs < 1200 )
                 edata->Vp = 2 * edata->Vs;
             else
                 edata->Vp = 1.87 * edata->Vs;

             updated_lambda = edata->rho * edata->Vp * edata->Vp;
         }

         if ( updated_lambda < 0) {
//             fprintf(stderr, "\nThread %d: %d element produces negative lambda = %.6f; Vp = %f; Vs = %f; Rho = %f",
//                     Global.myID, eindex, lambda, edata->Vp, edata->Vs, edata->rho);
//             MPI_Abort(MPI_COMM_WORLD, ERROR);


             fprintf(stderr, "\nThread : %d element produces negative lambda = %.6f; Vp = %f; Vs = %f; Rho = %f",
                      eindex, lambda, edata->Vp, edata->Vs, edata->rho);
             MPI_Abort(MPI_COMM_WORLD, ERROR);
         }

         /* update mu and lambda */

 		enlcons->mu = updated_mu;
 		enlcons->lambda = updated_lambda;

 		double old_c1 = ep->c1;



         /* update c1 - c4 */

         /* coefficients for term (deltaT_squared * Ke * Ut) */
         ep->c1 = theDeltaT *  theDeltaT *  edata->edgesize * updated_mu / 9;
         ep->c2 = theDeltaT *  theDeltaT *  edata->edgesize * updated_lambda / 9;

         zeta = 10 / edata->Vs;


         //printf("Element number: %i, old c: %f, new c: %f \n", el_eindex, old_c1,ep->c1);



         /* Damping update */

   	     if (eq_it == 0){
   	     enlcons->Qs_value = set_Qs(edata->Vs);
   	     enlcons->Qp_value = set_Qp(edata->Vs);
   	     }


 	    double anelastic_damping = 1/(2*enlcons->Qs_value);
 	    double total_requested_damping = GD.d/100 + anelastic_damping;



		 b = zeta * theBBase;
	     ep->c3 = b *theDeltaT *  theDeltaT * edata->edgesize * updated_mu / 9;
	     ep->c4 = b * theDeltaT *  theDeltaT * edata->edgesize * updated_lambda / 9;


//			/* Calculate total strains */
//			tstrains->qp[i] = point_strain(u, lx, ly, lz, h);
//			For each element you need to define another strain tensor, to only keep the maximum value of strain.
//            maxstrains->qp[i].xx = MAX(fabs(maxstrains->qp[i].xx),fabs(tstrains->qp[i].xx));
//            maxstrains->qp[i].yy = MAX(fabs(maxstrains->qp[i].yy),fabs(tstrains->qp[i].yy));
//            maxstrains->qp[i].zz = MAX(fabs(maxstrains->qp[i].zz),fabs(tstrains->qp[i].zz));
//            maxstrains->qp[i].xy = MAX(fabs(maxstrains->qp[i].xy),fabs(tstrains->qp[i].xy));
//            maxstrains->qp[i].yz = MAX(fabs(maxstrains->qp[i].yz),fabs(tstrains->qp[i].yz));
//            maxstrains->qp[i].xz = MAX(fabs(maxstrains->qp[i].xz),fabs(tstrains->qp[i].xz));
//            printf("This is xy max strain: %.20f (node %i) \n", maxstrains->qp[i].xy, i);
//
//			printf("This is strain: %f \n", tstrains[1])
//			/* strain and backstress predictor  */
//	        pstrains1->qp[i]    = copy_tensor ( pstrains2->qp[i] );     /* The strain predictor assumes that the current plastic strains equal those from the previous step   */
//	        alphastress1->qp[i] = copy_tensor ( alphastress2->qp[i] );
//			/* Calculate stresses */
//			if ( ( theMaterialModel == LINEAR ) || ( step <= theGeostaticFinalStep ) ){
//				stresses->qp[i]  = point_stress ( tstrains->qp[i], mu, lambda );
//				continue;
//			} else {
//				if ( theApproxGeoState == YES )
//					sigma0 = ApproxGravity_tensor(enlcons->sigmaZ_st, enlcons->phi, h, lz, edata->rho);
//				else
//					sigma0 = zero_tensor();
//				material_update ( *enlcons,           tstrains->qp[i],      pstrains1->qp[i], alphastress1->qp[i], epstr1->qv[i], sigma0, theDeltaT,
//						           &pstrains2->qp[i], &alphastress2->qp[i], &stresses->qp[i], &epstr2->qv[i],      &enlcons->fs[i]);



	     /* convert damping to Q */





	    new_Q = 1 / (2*total_requested_damping);


	    /* regenerate the Q table */

	    int i,k,j,l;
	    double myQTABLE[26][6];

	    int array_size=(QTable_Size*6)+1;

	    // regenertate the Q_table
	    k=0;
	    j=0;
		for(l = 1; l < array_size; l++)
		{
			i = l-1;
			myQTABLE[k][j] = theQTABLE[i];
	        if (l % 6 ==0){
	        	k=k+1;
	        }
	        j=j+1;
	        if (j==6){
	        	j=0;
	        }
		}

	    double index_Qs = Search_Quality_Table(new_Q,&(myQTABLE[0][0]), QTable_Size);

		update_Q_params(edata,index_Qs,0,QTable_Size,&(myQTABLE[0][0]),new_Q,0);
	    control_correction_factor(edata,theFreq_Vel,theFreq);



        //printf("strain_level= %f, table_damping: %f, Q_orig: %f, Q_new: %f. \n",myteststrain,GD.d,enlcons->Qs_value ,new_Q);
//			}
//		}
	} /* for all nonlinear elements */
}


/* -------------------------------------------------------------------------- */
/*                        Nonlinear Finalize and Stats                        */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/*                        Nonlinear Output to Stations                        */
/* -------------------------------------------------------------------------- */


void eqlinear_stations_init(mesh_t    *myMesh,
                             station_t *myStations,
                             int32_t    myNumberOfStations)
{

    if ( myNumberOfStations == 0 ) {
        return;
    }

    int32_t     eindex, el_eindex;
    int32_t     iStation=0;
    vector3D_t  point;
    octant_t   *octant;
    int32_t     lnid0;

    myNumberOfEqlinStations = 0;
    for (iStation = 0; iStation < myNumberOfStations; iStation++) {

        for ( el_eindex = 0; el_eindex < myEqlinElementsCount; el_eindex++ ) {

            /* capture the stations coordinates */
            point = myStations[iStation].coords;

            /* search the octant */
            if ( search_point(point, &octant) != 1 ) {
                fprintf(stderr,
                        "eqlinear_stations_init: "
                        "No octant with station coords\n");
                MPI_Abort(MPI_COMM_WORLD, ERROR);
                exit(1);
            }

            eindex = myEqlinElementsMapping[el_eindex];

            lnid0 = myMesh->elemTable[eindex].lnid[0];

            if ( (myMesh->nodeTable[lnid0].x == octant->lx) &&
                 (myMesh->nodeTable[lnid0].y == octant->ly) &&
                 (myMesh->nodeTable[lnid0].z == octant->lz) ) {

                /* I have a match for the element's origin */

                /* Now, perform level sanity check */
                if (myMesh->elemTable[eindex].level != octant->level) {
                    fprintf(stderr,
                            "eqlinear_stations_init: First pass: "
                            "Wrong level of octant\n");
                    MPI_Abort(MPI_COMM_WORLD, ERROR);
                    exit(1);
                }

                myNumberOfEqlinStations++;

                break;
            }
        }
    }

    XMALLOC_VAR_N( myStationsElementIndices, int32_t, myNumberOfEqlinStations);
    XMALLOC_VAR_N( myEqlinStationsMapping, int32_t, myNumberOfEqlinStations);
// // XMALLOC_VAR_N( myNonlinStations, nlstation_t, myNumberOfNonlinStations);
//
    int32_t elStationsCount = 0;
    for (iStation = 0; iStation < myNumberOfStations; iStation++) {

        for ( el_eindex = 0; el_eindex < myEqlinElementsCount; el_eindex++ ) {

            /* capture the stations coordinates */
            point = myStations[iStation].coords;

            /* search the octant */
            if ( search_point(point, &octant) != 1 ) {
                fprintf(stderr,
                        "eqlinear_stations_init: "
                        "No octant with station coords\n");
                MPI_Abort(MPI_COMM_WORLD, ERROR);
                exit(1);
            }

            eindex = myEqlinElementsMapping[el_eindex];

            lnid0 = myMesh->elemTable[eindex].lnid[0];

            if ( (myMesh->nodeTable[lnid0].x == octant->lx) &&
                 (myMesh->nodeTable[lnid0].y == octant->ly) &&
                 (myMesh->nodeTable[lnid0].z == octant->lz) ) {

                /* I have a match for the element's origin */

                /* Now, perform level sanity check */
                if (myMesh->elemTable[eindex].level != octant->level) {
                    fprintf(stderr,
                            "nonlinear_stations_init: Second pass: "
                            "Wrong level of octant\n");
                    MPI_Abort(MPI_COMM_WORLD, ERROR);
                    exit(1);
                }

                if ( elStationsCount >= myNumberOfEqlinStations ) {
                    fprintf(stderr,
                            "eqlinear_stations_init: Second pass: "
                            "More stations than initially counted\n");
                    MPI_Abort(MPI_COMM_WORLD, ERROR);
                    exit(1);
                }

                /* Store the element index and mapping to stations */
                myStationsElementIndices[elStationsCount] = el_eindex;
                myEqlinStationsMapping[elStationsCount] = iStation;

                elStationsCount++;

                break;
            }

        } /* for all my elements */

    } /* for all my stations */

/*    for ( iStation = 0; iStation < myNumberOfNonlinStations; iStation++ ) {

        tensor_t *stress, *strain, *pstrain1, *pstrain2;
        double   *ep1;

        strain   = &(myNonlinStations[iStation].strain);
        stress   = &(myNonlinStations[iStation].stress);
        pstrain1 = &(myNonlinStations[iStation].pstrain1);
        pstrain2 = &(myNonlinStations[iStation].pstrain2);
        ep1      = &(myNonlinStations[iStation].ep );
        *ep1     = 0.;

        init_tensorptr(strain);
        init_tensorptr(stress);
        init_tensorptr(pstrain1);
        init_tensorptr(pstrain2);

    }*/

}

