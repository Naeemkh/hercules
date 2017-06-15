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

//	elem_t  *elemp;
//	edata_t *edata;

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
       // printf("element_num : %i, mu: %f \n",el_eindex,mu);
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


noyesflag_t isThisElementsAtTheBottom_eq( mesh_t  *myMesh,
                                       int32_t  eindex,
                                       double   depth )
{
    elem_t  *elemp;
    int32_t  nindex;
    double   x_m,y_m,z_m;

    /* Capture the element's last node at the bottom */
    elemp  = &myMesh->elemTable[eindex];
    nindex = elemp->lnid[7];

    x_m = (myMesh->ticksize)*(double)myMesh->nodeTable[nindex].x;
    y_m = (myMesh->ticksize)*(double)myMesh->nodeTable[nindex].y;
    z_m = (myMesh->ticksize)*(double)myMesh->nodeTable[nindex].z;


  //  if ( z_m == depth && x_m < 1024 + 32 && x_m > 1024 - 32 && y_m < 1024 + 32 && y_m > 1024 - 32  ) {
  // important: This shows all elements in the middle column not bottom
  //  if ( z_m == depth && x_m == 992 && y_m == 992  ) {
    if ( x_m == 1000 && y_m == 1000  ) {
    return YES;
    }

    return NO;
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


int get_displacements_eq(mysolver_t *solver, elem_t *elemp, fvector_t *u) {

    int i;
    int res = 0;

    /* Capture displacements for each node */
    for (i = 0; i < 8; i++) {

        int32_t    lnid;
        fvector_t *dis;

        lnid = elemp->lnid[i];
        dis  = solver->tm1 + lnid;

        res += vector_is_all_zero( dis );

        u[i].f[0] = dis->f[0];
        u[i].f[1] = dis->f[1];
        u[i].f[2] = dis->f[2];

    }

    return res;
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



		if ( get_displacements_eq(mySolver, elemp, u) == 0 ) {
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



		    int32_t  nindex;
		    double   x_m,y_m,z_m;

		    /* Capture the element's last node at the bottom */
		    elemp  = &myMesh->elemTable[eindex];
		    nindex = elemp->lnid[7];

//		    x_m = (myMesh->ticksize)*(double)myMesh->nodeTable[nindex].x;
//		    y_m = (myMesh->ticksize)*(double)myMesh->nodeTable[nindex].y;
		    z_m = (myMesh->ticksize)*(double)myMesh->nodeTable[nindex].z;



// the following part of code is written to extract the strain.Uncomment it. Alos double check the isThisElementsAtTheBottom_eq function
//			if ( isThisElementsAtTheBottom_eq(myMesh, el_eindex, 512) == YES  && (step==3228 || step==4807)) {
//				printf("STST el_eindex = %i,node =%i,depth = %f, timestep = %i ,xz = %.20f \n",el_eindex,i,z_m,step,tstrains->qp[i].xz);
//			        }




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



void material_update_eq (      mesh_t     *myMesh,
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
	int32_t eindex, el_eindex, nindex;
//    int btn=0;


	/* Loop over the number of local elements */
	for (el_eindex = 0; el_eindex < myEqlinElementsCount; el_eindex++) {

		elem_t        *elemp;
		edata_t       *edata;
		elconstants_t *enlcons;
		e_t           *ep;    /* pointer to the element constant table */

		double         h;          /* Element edge-size in meters   */
		double         mu, lambda, original_mu; /* Elasticity material constants */
		double         XI, QC;
		fvector_t      u[8];
		eq_qptensors_t   *maxstrains;
		double         zeta, a, b, updated_mu, updated_lambda;
		double         requested_Q, new_Q;
        double         z_m;
        noyesflag_t    bottom_element=NO;

		/* Capture data from the element and mesh */
		eindex = myEqlinElementsMapping[el_eindex];

		elemp = &myMesh->elemTable[eindex];
		edata = (edata_t *)elemp->data;
		h     = edata->edgesize;
		ep    = &mySolver->eTable[eindex];

	    /* Capture the element's last node at the bottom */
	    elemp  = &myMesh->elemTable[eindex];
	    nindex = elemp->lnid[7];


		z_m = (myMesh->ticksize)*(double)myMesh->nodeTable[nindex].z;

		if  (z_m == 512){
			bottom_element=YES;
		}

		if (bottom_element==YES){
//			printf("this is bottom element number: %i\n",btn);
//			btn=btn+1;
			continue;
		}

		/* Capture data from the eqlinear element structure */
		enlcons = myEqlinSolver->constants + el_eindex;



		/* Keep the original mu in other variable and don't change it.
		 *
		 */
		if (eq_it == 0){
			enlcons->original_mu = enlcons -> mu;
		}

//		printf("Element number: %i, old Vs: %f \n", el_eindex, edata->Vs);

		original_mu     = enlcons->original_mu;
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
		myteststrain = strain_mat[i][4] + myteststrain;
		}

		// mean of all xy strains.
		myteststrain = 100* (myteststrain/8);

//		printf("Here is the max strain of element %i : %.10f \n", el_eindex, myteststrain);

		 GD_t GD = search_GD_table(myteststrain);

		//printf("Results : strain = %f , G = %f, D = %f \n", myteststrain, GD.g, GD.d);


		 updated_mu = original_mu * GD.g;

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

         /* update mu and lambda - Don't update the DRM element. */

 		enlcons->mu = updated_mu;
 		enlcons->lambda = updated_lambda;
 		edata -> Vs =  sqrt(updated_mu/edata ->rho);
 		edata -> Vp =  sqrt((updated_lambda+2*updated_mu)/edata ->rho);

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

//   	     printf("Element number: %i, new Vs: %f \n", el_eindex, edata->Vs);
//   	     printf("Element number: %i, Mu: %f, and updated Mu: %f \n", el_eindex, original_mu, enlcons->mu);

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

		update_Q_params(edata,index_Qs,0,QTable_Size,&(myQTABLE[0][0]),new_Q,0, theFreq);
	    control_correction_factor(edata,theFreq_Vel,theFreq);

	} /* for all nonlinear elements */

}


void    compute_addforce_bottom(int32_t timestep, mesh_t *myMesh, mysolver_t *mySolver)
{

    // dt = 0.01
    //double Force_1[1000]={1.35111255205411E-24,	3.31034338608617E-24,	6.14774662326503E-24,	1.02516565323711E-23,	1.61797605160242E-23,	2.47319038914087E-23,	3.70537343261079E-23,	5.47840515413847E-23,	8.02640110544514E-23,	1.16833748165512E-22,	1.69252387779081E-22,	2.44291965431406E-22,	3.51576141017056E-22,	5.04763006932892E-22,	7.23210915764468E-22,	1.03432141811737E-21,	1.47683009736597E-21,	2.10542256860622E-21,	2.99720053284144E-21,	4.26072715409923E-21,	6.0486635768945E-21,	8.5753984748259E-21,	1.21416102767648E-20,	1.71684420215414E-20,	2.42449864966829E-20,	3.41941763880942E-20,	4.81640893060695E-20,	6.77542991500259E-20,	9.5190488838295E-20,	1.33565430987247E-19,	1.87171110598859E-19,	2.61955839032891E-19,	3.6615254337669E-19,	5.11141156651131E-19,	7.12630723826186E-19,	9.92277545799426E-19,	1.37989718196362E-18,	1.91648413371098E-18,	2.65832849646046E-18,	3.682620415017E-18,	5.09506962403629E-18,	7.04025043229777E-18,	9.7156277975747E-18,	1.33905483033423E-17,	1.8431918256102E-17,	2.53388759593255E-17,	3.47895442634112E-17,	4.77039857230719E-17,	6.53288590160379E-17,	8.93511006442267E-17,	1.22050373325902E-16,	1.66503219226438E-16,	2.26856077338103E-16,	3.08689890021523E-16,	4.19506351456677E-16,	5.69375484170172E-16,	7.71796875740479E-16,	1.04484348379901E-15,	1.41267861281199E-15,	1.90756501986615E-15,	2.5725221072992E-15,	3.46483537275146E-15,	4.66068501609413E-15,	6.26124209656813E-15,	8.40068634756923E-15,	1.12567366838138E-14,	1.50644625260262E-14,	2.01343727579731E-14,	2.68760740092946E-14,	3.582916960792E-14,	4.77035585984831E-14,	6.34319206789198E-14,	8.42379757928853E-14,	1.11725134551318E-13,	1.47991468436758E-13,	1.95778603548699E-13,	2.58664275034035E-13,	3.41310992698064E-13,	4.49786713951593E-13,	5.9197778909583E-13,	7.78119974738096E-13,	1.02148030201538E-12,	1.33923140895941E-12,	1.75357106332519E-12,	2.29315358820295E-12,	2.99491747048686E-12,	3.90641546424848E-12,	5.08878108795893E-12,	6.62049990205187E-12,	8.60219699899902E-12,	1.11627057797399E-11,	1.44667498290573E-11,	1.87246526024766E-11,	2.42045924457432E-11,	3.12480477592831E-11,	4.02892344681907E-11,	5.18795321427313E-11,	6.67181343657044E-11,	8.56904532351597E-11,	1.09916169287681E-10,	1.40809261192854E-10,	1.80152891853103E-10,	2.30192690150975E-10,	2.9375277585319E-10,	3.74379859459252E-10,	4.76521945362924E-10,	6.05749618883788E-10,	7.69029657097426E-10,	9.75062831243446E-10,	1.23470033752398E-09,	1.56145638937665E-09,	1.97213822946126E-09,	2.48761929344654E-09,	3.13378662113627E-09,	3.94270002891016E-09,	4.95400822453848E-09,	6.21667618744348E-09,	7.79108901822214E-09,	9.75161038754459E-09,	1.21896890424433E-08,	1.52176249696887E-08,	1.897312824741E-08,	2.36248288857415E-08,	2.93789256951777E-08,	3.64871971498006E-08,	4.52566381572313E-08,	5.60610345494614E-08,	6.93548430389403E-08,	8.5689809555945E-08,	1.05734834668221E-07,	1.302996827284E-07,	1.60363233236273E-07,	1.97107090587624E-07,	2.41955504081723E-07,	2.96622706192584E-07,	3.63168956354562E-07,	4.44066782829773E-07,	5.42279149869267E-07,	6.61351544918343E-07,	8.05520284939971E-07,	9.79839686263406E-07,	1.19033113327213E-06,	1.44415752263147E-06,	1.74982705686299E-06,	2.11743091946632E-06,	2.55891998933804E-06,	3.08842645108543E-06,	3.72263693635386E-06,	4.48122469602727E-06,	5.38734926077733E-06,	6.46823310459243E-06,	7.75582598918887E-06,	9.2875689429452E-06,	1.11072712221821E-05,	1.32661151206558E-05,	1.58238051397037E-05,	0.00001884987981024,	2.24252063711497E-05,	0.00002664368055741,	3.16141559344243E-05,	3.74626295292922E-05,	4.43347129490715E-05,	5.23984207315135E-05,	6.18473103326552E-05,	7.29040109014669E-05,	8.58241808033325E-05,	0.000100900936705338,	0.000118469799895386,	0.000138914208336198,	0.000162671645709822,	0.000190240441336642,	0.000222187297295956,	0.000259155601265509,	0.000301874585459455,	0.000351169393493627,	0.00040797211795069,	0.000473333871753659,	0.000548437956073796,	0.000634614186278792,	0.00073335443524235,	0.000846329450052448,	0.000975406993632119,	0.00112267135687785,	0.0012904442794777,	0.00148130730844311,	0.00169812561242513,	0.00194407325694054,	0.00222265993056453,	0.0025377590948236,	0.00289363751082551,	0.00329498607349048,	0.00374695185951989,	0.00425517126789866,	0.00482580410175243,	0.00546556840778054,	0.00618177585431291,	0.0069823673913907,	0.00787594889629625,	0.00887182646586039,	0.00998004097292257,	0.0112114014588398,	0.0125775168873351,	0.0140908257377163,	0.0157646228681269,	0.0176130830326376,	0.0196512803903425,	0.0218952033009686,	0.024361763660678,	0.027068799994662,	0.0300350734907684,	0.0332802561317892,	0.0368249100642666,	0.0406904573298339,	0.0448991390823394,	0.0494739634214452,	0.054438640992151,	0.0598175075308782,	0.0656354325833704,	0.0719177136787083,	0.0786899553180234,	0.0859779322267872,	0.0938074364264016,	0.10220410780464,	0.111193248005461,	0.120799617616832,	0.131047216810128,	0.141959049775888,	0.153556873507335,	0.165860931703912,	0.178889674800692,	0.192659467374028,	0.207184284427028,	0.222475398317956,	0.23854105835747,	0.255386165363778,	0.273011943724686,	0.291415613768501,	0.310590067487947,	0.330523550888371,	0.351199356439436,	0.372595529293766,	0.394684591092216,	0.417433285299343,	0.440802348099881,	0.464746308933637,	0.489213324748381,	0.514145052004557,	0.539476560369099,	0.565136291885586,	0.591046069202778,	0.617121156181683,	0.643270373882431,	0.669396274556464,	0.695395375838056,	0.721158456843762,	0.746570917352214,	0.771513200653001,	0.795861280027415,	0.819487208160658,	0.842259728091535,	0.86404494358855,	0.884707046109044,	0.90410909475876,	0.922113844932701,	0.938584620593488,	0.95338622444111,	0.966385879557766,	0.977454195483903,	0.986466151106528,	0.993302086228303,	0.997848693245159,	1,	0.999658334608501,	0.996735262877173,	0.99115248885963,	0.982842709129457,	0.971750411490379,	0.957832609098869,	0.941059501341448,	0.921415053287494,	0.898897486126068,	0.873519671687591,	0.845309424942457,	0.814309689251451,	0.780578610108199,	0.744189494151644,	0.70523065132514,	0.663805119205529,	0.620030269706927,	0.57403729956559,	0.525970607219161,	0.475987059890558,	0.42425515585854,	0.370954088028045,	0.316272715989177,	0.260408454759576,	0.203566089327616,	0.145956524940774,	0.0877954838034294,	0.0293021594514236,	-0.0293021594514245,	-0.0877954838034304,	-0.145956524940774,	-0.203566089327617,	-0.260408454759578,	-0.316272715989178,	-0.370954088028045,	-0.424255155858542,	-0.47598705989056,	-0.525970607219161,	-0.57403729956559,	-0.620030269706927,	-0.663805119205531,	-0.705230651325142,	-0.744189494151644,	-0.780578610108199,	-0.814309689251451,	-0.845309424942458,	-0.873519671687591,	-0.898897486126069,	-0.921415053287494,	-0.941059501341448,	-0.957832609098869,	-0.971750411490379,	-0.982842709129457,	-0.99115248885963,	-0.996735262877174,	-0.999658334608503,	-1,	-0.997848693245161,	-0.993302086228305,	-0.986466151106529,	-0.977454195483903,	-0.966385879557766,	-0.95338622444111,	-0.938584620593488,	-0.922113844932701,	-0.904109094758762,	-0.884707046109044,	-0.86404494358855,	-0.842259728091535,	-0.81948720816066,	-0.795861280027415,	-0.771513200653001,	-0.746570917352214,	-0.721158456843762,	-0.695395375838056,	-0.669396274556466,	-0.643270373882431,	-0.617121156181683,	-0.591046069202778,	-0.565136291885586,	-0.539476560369099,	-0.514145052004557,	-0.489213324748381,	-0.464746308933637,	-0.440802348099881,	-0.417433285299343,	-0.394684591092216,	-0.372595529293766,	-0.351199356439436,	-0.330523550888371,	-0.310590067487949,	-0.291415613768501,	-0.273011943724686,	-0.255386165363778,	-0.23854105835747,	-0.222475398317956,	-0.207184284427028,	-0.192659467374028,	-0.178889674800692,	-0.165860931703912,	-0.153556873507335,	-0.141959049775888,	-0.131047216810128,	-0.120799617616832,	-0.111193248005461,	-0.10220410780464,	-0.0938074364264019,	-0.0859779322267875,	-0.0786899553180237,	-0.0719177136787086,	-0.0656354325833708,	-0.0598175075308785,	-0.0544386409921513,	-0.0494739634214455,	-0.0448991390823399,	-0.0406904573298344,	-0.0368249100642671,	-0.0332802561317895,	-0.0300350734907687,	-0.0270687999946623,	-0.0243617636606782,	-0.021895203300969,	-0.0196512803903429,	-0.0176130830326379,	-0.0157646228681272,	-0.0140908257377166,	-0.0125775168873355,	-0.0112114014588401,	-0.00998004097292293,	-0.00887182646586075,	-0.00787594889629661,	-0.00698236739139106,	-0.00618177585431327,	-0.0054655684077809,	-0.00482580410175279,	-0.00425517126789902,	-0.00374695185952026,	-0.00329498607349084,	-0.00289363751082587,	-0.00253775909482396,	-0.00222265993056489,	-0.00194407325694091,	-0.00169812561242549,	-0.00148130730844348,	-0.00129044427947807,	-0.00112267135687822,	-0.000975406993632486,	-0.000846329450052815,	-0.000733354435242717,	-0.000634614186279159,	-0.000548437956074162,	-0.000473333871754026,	-0.000407972117951058,	-0.000351169393493994,	-0.000301874585459822,	-0.000259155601265877,	-0.000222187297296324,	-0.00019024044133701,	-0.00016267164571019,	-0.000138914208336566,	-0.000118469799895754,	-0.000100900936705706,	-8.58241808037005E-05,	-7.29040109018349E-05,	-6.18473103330234E-05,	-5.23984207318817E-05,	-4.43347129494397E-05,	-3.74626295296604E-05,	-3.16141559347925E-05,	-2.66436805577782E-05,	-2.24252063715181E-05,	-1.88498798106082E-05,	-1.58238051400721E-05,	-1.32661151210241E-05,	-1.11072712225503E-05,	-9.28756894331343E-06,	-7.75582598955712E-06,	-6.46823310496068E-06,	-5.38734926114557E-06,	-4.48122469639552E-06,	-3.72263693672212E-06,	-3.08842645145369E-06,	-2.55891998970628E-06,	-2.11743091983458E-06,	-1.74982705723125E-06,	-1.44415752299972E-06,	-1.19033113364039E-06,	-9.79839686631661E-07,	-8.05520285308228E-07,	-6.61351545286602E-07,	-5.42279150237526E-07,	-4.44066783198031E-07,	-3.6316895672282E-07,	-2.96622706560844E-07,	-2.41955504449981E-07,	-1.97107090955884E-07,	-1.60363233604531E-07,	-1.30299683096659E-07,	-1.0573483503648E-07,	-8.56898099242044E-08,	-6.93548434071995E-08,	-5.60610349177208E-08,	-4.52566385254907E-08,	-3.648719751806E-08,	-2.93789260634371E-08,	-2.36248292540009E-08,	-1.89731286156694E-08,	-1.52176253379481E-08,	-1.21896894107027E-08,	-9.75161075580404E-09,	-7.79108938648158E-09,	-6.21667655570292E-09,	-4.95400859279792E-09,	-3.9427003971696E-09,	-3.13378698939575E-09,	-2.487619661706E-09,	-1.97213859772073E-09,	-1.56145675763611E-09,	-1.23470070578345E-09,	-9.75063199502915E-10,	-7.69030025356895E-10,	-6.05749987143255E-10,	-4.76522313622393E-10,	-3.74380227718723E-10,	-2.93753144112661E-10,	-2.30193058410446E-10,	-1.80153260112575E-10,	-1.40809629452325E-10,	-1.09916537547152E-10,	-8.56908214946309E-11,	-6.67185026251758E-11,	-5.18799004022028E-11,	-4.0289602727662E-11,	-3.12484160187545E-11,	-2.42049607052147E-11,	-1.87250208619481E-11,	-1.44671180885289E-11,	-1.11630740392114E-11,	-8.60256525847056E-12,	-6.62086816152344E-12,	-5.0891493474305E-12,	-3.90678372372004E-12,	-2.99528572995842E-12,	-2.29352184767451E-12,	-1.75393932279676E-12,	-1.33959966843097E-12,	-1.02184856148694E-12,	-7.78488234209659E-13,	-5.92346048567396E-13,	-4.50154973423158E-13,	-3.4167925216963E-13,	-2.590325345056E-13,	-1.96146863020266E-13,	-1.48359727908323E-13,	-1.12093394022884E-13,	-8.46062352644517E-14,	-6.38001801504864E-14,	-4.80718180700495E-14,	-3.61974290794863E-14,	-2.72443334808611E-14,	-2.05026322295396E-14,	-1.54327219975927E-14,	-1.16249961553804E-14,	-8.76894581913584E-15,	-6.62950156813477E-15,	-5.02894448766078E-15,	-3.83309484431811E-15,	-2.94078157886584E-15,	-2.2758244914328E-15,	-1.78093808437863E-15,	-1.41310295536566E-15,	-1.14005634730712E-15,	-9.37634955736819E-16,	-7.87765823023325E-16,	-6.76949361588171E-16,	-5.95115548904751E-16,	-5.34762690793086E-16,	-4.90309844892551E-16,	-4.57610572210875E-16,	-4.33588330582687E-16,	-4.15963457289721E-16,	-4.03049015830059E-16,	-3.93598347525974E-16,	-3.8669138982275E-16,	-3.81650019869991E-16,	-3.77975099364223E-16,	-3.75299721998947E-16,	-3.73354541190686E-16,	-3.71942091981665E-16,	-3.70917800063109E-16,	-3.7017595570036E-16,	-3.69639368748613E-16,	-3.69251749112448E-16,	-3.68972102290475E-16,	-3.687706127233E-16,	-3.68625624110026E-16,	-3.68521427405681E-16,	-3.68446642677248E-16,	-3.68393036997636E-16,	-3.68354662055488E-16,	-3.68327225865799E-16,	-3.68307635655955E-16,	-3.68293665743037E-16,	-3.68283716553145E-16,	-3.6827664000867E-16,	-3.68271613176926E-16,	-3.68268046965124E-16,	-3.68265520230226E-16,	-3.68263732293804E-16,	-3.68262468767182E-16,	-3.68261576989217E-16,	-3.68260948396746E-16,	-3.68260505888067E-16,	-3.68260194777565E-16,	-3.68259976329656E-16,	-3.68259823142789E-16,	-3.68259715858614E-16,	-3.68259640819036E-16,	-3.68259588400396E-16,	-3.6825955183066E-16,	-3.682595263507E-16,	-3.68259508620383E-16,	-3.68259496298553E-16,	-3.68259487746409E-16,	-3.68259481818306E-16,	-3.68259477714395E-16,	-3.68259474876992E-16,	-3.68259472917762E-16,	-3.68259471566649E-16,	-3.68259470636101E-16,	-3.6825946999603E-16,	-3.68259469556328E-16,	-3.6825946925466E-16,	-3.68259469047959E-16,	-3.68259468906509E-16,	-3.68259468809838E-16,	-3.68259468743854E-16,	-3.68259468698874E-16,	-3.68259468668252E-16,	-3.68259468647431E-16,	-3.68259468633291E-16,	-3.68259468623704E-16,	-3.68259468617208E-16,	-3.68259468612815E-16,	-3.68259468609847E-16,	-3.68259468607844E-16,	-3.68259468606495E-16,	-3.68259468605587E-16,	-3.68259468604977E-16,	-3.68259468604566E-16,	-3.68259468604293E-16,	-3.68259468604109E-16,	-3.68259468603986E-16,	-3.68259468603903E-16,	-3.68259468603848E-16,	-3.68259468603812E-16,	-3.68259468603788E-16,	-3.68259468603772E-16,	-3.6825946860376E-16,	-3.68259468603754E-16,	-3.68259468603749E-16,	-3.68259468603746E-16,	-3.68259468603744E-16,	-3.68259468603743E-16,	-3.68259468603741E-16,	-3.68259468603741E-16,	-3.68259468603741E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16};

    // dt = 0.001
    double Force_1[8000] ={0.0106459760945398,	0.0107703980759154,	0.0108961256240057,	0.0110231705082228,	0.0111515445802402,	0.0112812597742825,	0.0114123281074118,	0.0115447616798117,	0.0116785726750686,	0.0118137733604501,	0.0119503760871806,	0.0120883932907139,	0.0122278374910027,	0.0123687212927654,	0.0125110573857494,	0.0126548585449916,	0.0128001376310751,	0.0129469075903832,	0.0130951814553497,	0.0132449723447058,	0.0133962934637236,	0.0135491581044561,	0.0137035796459731,	0.013859571554594,	0.0140171473841166,	0.014176320776042,	0.0143371054597956,	0.0144995152529445,	0.01466356406141,	0.0148292658796772,	0.0149966347909995,	0.0151656849675989,	0.0153364306708628,	0.0155088862515354,	0.0156830661499052,	0.0158589848959883,	0.016036657109706,	0.0162160975010589,	0.0163973208702957,	0.0165803421080771,	0.0167651761956351,	0.0169518382049273,	0.017140343298786,	0.017330706731062,	0.017522943846764,	0.0177170700821918,	0.0179131009650649,	0.018111052114645,	0.0183109392418537,	0.0185127781493843,	0.0187165847318076,	0.0189223749756729,	0.0191301649596021,	0.019339970854379,	0.0195518089230317,	0.0197656955209094,	0.0199816470957533,	0.020199680187761,	0.0204198114296444,	0.0206420575466819,	0.0208664353567636,	0.0210929617704301,	0.0213216537909051,	0.0215525285141205,	0.0217856031287356,	0.0220208949161491,	0.0222584212505038,	0.0224981995986849,	0.022740247520311,	0.0229845826677178,	0.0232312227859345,	0.0234801857126529,	0.0237314893781895,	0.0239851518054395,	0.0242411911098232,	0.0244996254992253,	0.0247604732739257,	0.0250237528265234,	0.0252894826418511,	0.025557681296883,	0.0258283674606343,	0.0261015598940514,	0.0263772774498956,	0.0266555390726168,	0.0269363637982198,	0.0272197707541214,	0.0275057791589995,	0.0277944083226329,	0.0280856776457326,	0.0283796066197641,	0.0286762148267608,	0.0289755219391284,	0.0292775477194399,	0.0295823120202214,	0.0298898347837291,	0.0302001360417161,	0.0305132359151899,	0.0308291546141611,	0.0311479124373812,	0.0314695297720713,	0.0317940270936411,	0.0321214249653976,	0.0324517440382437,	0.0327850050503673,	0.0331212288269192,	0.0334604362796816,	0.0338026484067259,	0.0341478862920598,	0.0344961711052644,	0.0348475241011197,	0.0352019666192207,	0.0355595200835812,	0.0359202060022281,	0.0362840459667837,	0.0366510616520372,	0.0370212748155051,	0.0373947072969807,	0.0377713810180715,	0.0381513179817254,	0.0385345402717461,	0.038921070052296,	0.0393109295673877,	0.0397041411403639,	0.0401007271733653,	0.0405007101467864,	0.0409041126187198,	0.0413109572243878,	0.041721266675562,	0.0421350637599708,	0.0425523713406941,	0.0429732123555462,	0.0433976098164458,	0.0438255868087729,	0.0442571664907142,	0.0446923720925942,	0.0451312269161953,	0.0455737543340627,	0.0460199777887987,	0.0464699207923417,	0.0469236069252337,	0.0473810598358733,	0.0478423032397565,	0.0483073609187031,	0.0487762567200698,	0.0492490145559507,	0.0497256584023625,	0.0502062122984173,	0.0506907003454812,	0.0511791467063188,	0.0516715756042239,	0.0521680113221362,	0.0526684782017442,	0.0531730006425731,	0.0536816031010594,	0.0541943100896106,	0.0547111461756506,	0.0552321359806511,	0.0557573041791475,	0.0562866754977413,	0.0568202747140872,	0.0573581266558656,	0.0579002561997403,	0.0584466882703014,	0.0589974478389933,	0.0595525599230276,	0.0601120495842811,	0.0606759419281785,	0.06124426210256,	0.0618170352965339,	0.0623942867393133,	0.0629760416990385,	0.0635623254815825,	0.0641531634293427,	0.0647485809200157,	0.0653486033653575,	0.0659532562099271,	0.0665625649298156,	0.0671765550313585,	0.0677952520498328,	0.0684186815481384,	0.0690468691154626,	0.0696798403659305,	0.0703176209372374,	0.0709602364892665,	0.0716077127026904,	0.0722600752775557,	0.0729173499318527,	0.0735795624000676,	0.0742467384317195,	0.0749189037898807,	0.0755960842496807,	0.0762783055967938,	0.0769655936259106,	0.0776579741391929,	0.0783554729447123,	0.079058115854872,	0.0797659286848128,	0.0804789372508016,	0.0811971673686045,	0.0819206448518422,	0.0826493955103299,	0.0833834451483994,	0.0841228195632057,	0.0848675445430166,	0.0856176458654848,	0.0863731492959047,	0.0871340805854516,	0.0879004654694042,	0.0886723296653511,	0.0894496988713795,	0.0902325987642482,	0.0910210549975432,	0.0918150931998166,	0.0926147389727094,	0.0934200178890564,	0.0942309554909757,	0.0950475772879402,	0.0958699087548335,	0.0966979753299877,	0.0975318024132063,	0.0983714153637684,	0.0992168394984177,	0.100068100089334,	0.100925222362089,	0.101788231493584,	0.10265715260997,	0.103532010784558,	0.104412831035701,	0.105299638324671,	0.106192457553513,	0.107091313562882,	0.107996231129868,	0.108907234965802,	0.109824349714042,	0.110747599947751,	0.11167701016765,	0.112612604799761,	0.11355440819313,	0.114502444617535,	0.115456738261177,	0.11641731322836,	0.117384193537144,	0.118357403116993,	0.119336965806403,	0.120322905350511,	0.121315245398693,	0.122314009502143,	0.123319221111437,	0.124330903574083,	0.125349080132054,	0.126373773919304,	0.127405007959271,	0.128442805162366,	0.129487188323438,	0.130538180119239,	0.131595803105858,	0.132660079716149,	0.133731032257144,	0.134808682907448,	0.135893053714618,	0.136984166592534,	0.138082043318747,	0.139186705531819,	0.140298174728646,	0.141416472261765,	0.142541619336652,	0.143673637008998,	0.144812546181981,	0.145958367603512,	0.147111121863481,	0.148270829390977,	0.149437510451505,	0.150611185144177,	0.151791873398904,	0.152979594973565,	0.154174369451168,	0.15537621623699,	0.156585154555721,	0.157801203448573,	0.159024381770399,	0.16025470818678,	0.161492201171114,	0.162736879001689,	0.163988759758736,	0.165247861321485,	0.166514201365197,	0.16778779735819,	0.169068666558854,	0.170356826012652,	0.171652292549115,	0.172955082778819,	0.174265213090358,	0.175582699647306,	0.17690755838516,	0.17823980500829,	0.179579454986861,	0.180926523553754,	0.182281025701482,	0.183642976179087,	0.185012389489032,	0.186389279884086,	0.187773661364201,	0.189165547673372,	0.190564952296502,	0.191971888456246,	0.193386369109854,	0.194808406946009,	0.196238014381647,	0.197675203558783,	0.199119986341319,	0.200572374311848,	0.20203237876846,	0.203500010721526,	0.204975280890488,	0.206458199700642,	0.207948777279908,	0.2094470234556,	0.210952947751192,	0.212466559383078,	0.21398786725732,	0.215516879966406,	0.217053605785989,	0.218598052671631,	0.220150228255543,	0.221710139843316,	0.223277794410652,	0.224853198600096,	0.226436358717758,	0.228027280730037,	0.229625970260342,	0.231232432585814,	0.232846672634036,	0.23446869497976,	0.236098503841611,	0.237736103078809,	0.239381496187881,	0.241034686299375,	0.24269567617457,	0.244364468202197,	0.24604106439515,	0.247725466387201,	0.24941767542972,	0.251117692388392,	0.252825517739941,	0.254541151568848,	0.256264593564079,	0.257995843015816,	0.259734898812183,	0.261481759435986,	0.263236422961449,	0.264998887050959,	0.266769148951809,	0.268547205492955,	0.270333053081771,	0.272126687700811,	0.273928104904578,	0.2757372998163,	0.277554267124706,	0.279379001080818,	0.281211495494744,	0.283051743732481,	0.284899738712723,	0.286755472903679,	0.288618938319903,	0.290490126519125,	0.292369028599096,	0.294255635194445,	0.296149936473539,	0.298051922135359,	0.299961581406386,	0.301878903037494,	0.303803875300858,	0.305736485986873,	0.307676722401084,	0.309624571361131,	0.311580019193702,	0.313543051731505,	0.315513654310248,	0.317491811765639,	0.319477508430393,	0.32147072813126,	0.323471454186065,	0.325479669400765,	0.327495356066519,	0.329518495956778,	0.331549070324388,	0.333587059898713,	0.335632444882774,	0.337685204950403,	0.339745319243422,	0.341812766368835,	0.34388752439604,	0.345969570854058,	0.348058882728791,	0.350155436460288,	0.352259207940039,	0.354370172508286,	0.356488304951362,	0.358613579499043,	0.360745969821926,	0.362885449028831,	0.365031989664227,	0.367185563705672,	0.369346142561291,	0.371513697067266,	0.373688197485356,	0.375869613500445,	0.378057914218104,	0.380253068162195,	0.382455043272485,	0.384663806902298,	0.38687932581619,	0.389101566187649,	0.391330493596829,	0.393566073028303,	0.395808268868856,	0.398057044905296,	0.400312364322302,	0.402574189700299,	0.404842483013363,	0.407117205627157,	0.409398318296901,	0.411685781165365,	0.413979553760905,	0.416279594995524,	0.418585863162968,	0.420898315936852,	0.423216910368828,	0.425541602886775,	0.427872349293033,	0.430209104762668,	0.432551823841773,	0.4349004604458,	0.437254967857938,	0.439615298727516,	0.441981405068451,	0.44435323825773,	0.446730749033928,	0.449113887495768,	0.451502603100715,	0.453896844663613,	0.456296560355357,	0.458701697701609,	0.46111220358155,	0.463528024226673,	0.465949105219619,	0.468375391493054,	0.470806827328584,	0.473243356355714,	0.475684921550854,	0.478131465236356,	0.480582929079608,	0.483039254092161,	0.485500380628904,	0.487966248387284,	0.490436796406568,	0.492911963067153,	0.495391686089917,	0.49787590253562,	0.500364548804352,	0.502857560635019,	0.505354873104885,	0.507856420629159,	0.510362136960627,	0.512871955189331,	0.515385807742303,	0.517903626383339,	0.520425342212827,	0.522950885667626,	0.525480186520988,	0.528013173882535,	0.530549776198286,	0.533089921250735,	0.535633536158974,	0.538180547378878,	0.540730880703333,	0.543284461262514,	0.545841213524231,	0.548401061294306,	0.550963927717018,	0.553529735275597,	0.55609840579277,	0.558669860431363,	0.561244019694956,	0.56382080342859,	0.566400130819536,	0.568981920398113,	0.571566090038558,	0.574152556959965,	0.576741237727262,	0.579332048252261,	0.581924903794751,	0.584519718963657,	0.587116407718251,	0.589714883369423,	0.592315058581005,	0.59491684537116,	0.597520155113819,	0.60012489854019,	0.60273098574031,	0.605338326164667,	0.607946828625876,	0.610556401300415,	0.613166951730422,	0.615778386825547,	0.618390612864871,	0.621003535498877,	0.623617059751488,	0.62623109002216,	0.628845530088041,	0.631460283106186,	0.634075251615833,	0.636690337540745,	0.639305442191612,	0.641920466268506,	0.644535309863409,	0.647149872462799,	0.649764052950295,	0.652377749609365,	0.6549908601261,	0.657603281592048,	0.660214910507107,	0.662825642782486,	0.665435373743728,	0.668043998133791,	0.670651410116198,	0.673257503278247,	0.675862170634285,	0.678465304629041,	0.681066797141034,	0.683666539486027,	0.686264422420562,	0.688860336145547,	0.69145417030991,	0.694045814014318,	0.696635155814957,	0.69922208372738,	0.701806485230413,	0.704388247270131,	0.70696725626389,	0.709543398104433,	0.712116558164052,	0.714686621298813,	0.717253471852856,	0.719816993662743,	0.722377070061882,	0.724933583885009,	0.727486417472738,	0.730035452676168,	0.732580570861562,	0.735121652915081,	0.737658579247592,	0.740191229799529,	0.742719484045822,	0.745243221000892,	0.747762319223707,	0.750276656822898,	0.752786111461943,	0.755290560364414,	0.757789880319281,	0.760283947686286,	0.762772638401376,	0.765255827982197,	0.767733391533656,	0.770205203753538,	0.772671138938192,	0.775131070988274,	0.777584873414553,	0.780032419343781,	0.78247358152462,	0.784908232333636,	0.787336243781348,	0.789757487518341,	0.792171834841441,	0.794579156699944,	0.796979323701918,	0.799372206120548,	0.801757673900555,	0.804135596664665,	0.806505843720145,	0.808868284065391,	0.811222786396578,	0.81356921911437,	0.815907450330685,	0.818237347875516,	0.820558779303818,	0.822871611902442,	0.825175712697132,	0.82747094845958,	0.829757185714529,	0.832034290746941,	0.834302129609216,	0.836560568128466,	0.838809471913848,	0.841048706363943,	0.843278136674201,	0.845497627844426,	0.847707044686328,	0.84990625183112,	0.852095113737164,	0.854273494697682,	0.856441258848505,	0.858598270175883,	0.860744392524339,	0.862879489604578,	0.865003425001445,	0.867116062181933,	0.869217264503233,	0.871306895220848,	0.873384817496739,	0.875450894407526,	0.87750498895274,	0.879546964063111,	0.881576682608913,	0.88359400740835,	0.885598801235981,	0.887590926831206,	0.889570246906775,	0.891536624157359,	0.893489921268152,	0.89543000092352,	0.897356725815696,	0.899269958653504,	0.901169562171137,	0.903055399136969,	0.904927332362404,	0.90678522471077,	0.908628939106248,	0.910458338542838,	0.912273286093363,	0.914073644918511,	0.915859278275909,	0.91763004952924,	0.919385822157381,	0.921126459763591,	0.922851826084721,	0.924561785000459,	0.926256200542613,	0.927934936904412,	0.929597858449853,	0.93124482972307,	0.932875715457728,	0.934490380586457,	0.936088690250303,	0.937670509808209,	0.939235704846532,	0.940784141188565,	0.942315684904107,	0.943830202319043,	0.945327560024946,	0.946807624888717,	0.948270264062224,	0.949715344991985,	0.951142735428857,	0.952552303437748,	0.953943917407351,	0.955317446059891,	0.956672758460895,	0.958009724028974,	0.959328212545622,	0.96062809416503,	0.961909239423918,	0.96317151925137,	0.964414804978696,	0.96563896834929,	0.966843881528514,	0.96802941711358,	0.969195448143445,	0.970341848108722,	0.971468490961585,	0.972575251125687,	0.973662003506093,	0.974728623499198,	0.975774987002669,	0.976800970425375,	0.977806450697333,	0.97879130527964,	0.979755412174421,	0.980698649934768,	0.981620897674678,	0.98252203507899,	0.983401942413323,	0.984260500534005,	0.985097590897999,	0.985913095572824,	0.986706897246468,	0.987478879237298,	0.988228925503954,	0.98895692065524,	0.989662749960007,	0.990346299357012,	0.991007455464783,	0.991646105591456,	0.99226213774461,	0.992855440641074,	0.993425903716733,	0.993973417136307,	0.994497871803116,	0.994999159368826,	0.995477172243177,	0.995931803603691,	0.996362947405351,	0.996770498390272,	0.997154352097338,	0.997514404871815,	0.997850553874948,	0.998162697093517,	0.998450733349386,	0.998714562309003,	0.998954084492888,	0.999169201285082,	0.999359814942568,	0.99952582860466,	0.999667146302361,	0.999783672967683,	0.99987531444294,	0.999941977489997,	0.999983569799489,	1,	0.999991177667203,	0.999957013332964,	0.999897418494402,	0.99981230562291,	0.999701588173132,	0.999565180591903,	0.999402998327135,	0.999214957836667,	0.999000976597065,	0.998760973112375,	0.998494866922832,	0.998202578613514,	0.997884029822952,	0.997539143251689,	0.997167842670781,	0.996770052930259,	0.996345699967519,	0.995894710815678,	0.995417013611859,	0.994912537605427,	0.994381213166168,	0.993822971792409,	0.993237746119078,	0.992625469925707,	0.991986078144372,	0.991319506867572,	0.990625693356047,	0.989904576046532,	0.989156094559442,	0.988380189706505,	0.98757680349831,	0.986745879151806,	0.985887361097723,	0.985001194987927,	0.984087327702706,	0.983145707357984,	0.982176283312468,	0.981179006174716,	0.98015382781014,	0.979100701347927,	0.978019581187895,	0.976910423007262,	0.97577318376735,	0.974607821720206,	0.973414296415143,	0.972192568705208,	0.970942600753568,	0.969664356039814,	0.968357799366183,	0.967022896863707,	0.965659615998264,	0.964267925576562,	0.962847795752025,	0.961399198030601,	0.959922105276486,	0.958416491717755,	0.95688233295191,	0.955319605951339,	0.953728289068686,	0.952108362042132,	0.950459806000586,	0.948782603468782,	0.947076738372291,	0.945342196042433,	0.9435789632211,	0.941787028065486,	0.93996638015272,	0.938117010484406,	0.936238911491065,	0.934332077036484,	0.932396502421965,	0.930432184390478,	0.928439121130714,	0.926417312281039,	0.924366758933352,	0.922287463636841,	0.920179430401634,	0.918042664702354,	0.915877173481572,	0.913682965153154,	0.911460049605507,	0.909208438204722,	0.906928143797612,	0.904619180714645,	0.902281564772772,	0.89991531327815,	0.897520445028759,	0.89509698031691,	0.89264494093165,	0.890164350161052,	0.887655232794407,	0.885117615124299,	0.882551524948574,	0.879956991572201,	0.877334045809019,	0.874682719983379,	0.872003047931669,	0.869295065003736,	0.866558808064187,	0.863794315493587,	0.861001627189536,	0.858180784567643,	0.855331830562377,	0.852454809627817,	0.84954976773827,	0.846616752388799,	0.843655812595612,	0.840666998896358,	0.837650363350291,	0.834605959538332,	0.831533842563008,	0.828434069048276,	0.825306697139236,	0.822151786501721,	0.818969398321778,	0.815759595305026,	0.812522441675903,	0.809258003176791,	0.805966347067026,	0.802647542121796,	0.799301658630912,	0.795928768397466,	0.792528944736379,	0.789102262472815,	0.785648797940495,	0.782168628979877,	0.778661834936232,	0.775128496657593,	0.771568696492589,	0.767982518288159,	0.764370047387153,	0.76073137062581,	0.757066576331118,	0.753375754318056,	0.749658995886725,	0.745916393819348,	0.742148042377163,	0.738354037297192,	0.734534475788896,	0.730689456530709,	0.726819079666455,	0.722923446801649,	0.719002660999678,	0.715056826777867,	0.711086050103429,	0.70709043838929,	0.703070100489808,	0.699025146696365,	0.694955688732855,	0.690861839751038,	0.686743714325798,	0.682601428450267,	0.678435099530847,	0.674244846382104,	0.670030789221564,	0.665793049664373,	0.661531750717861,	0.65724701677598,	0.652938973613631,	0.648607748380882,	0.644253469597068,	0.639876267144773,	0.635476272263715,	0.6310536175445,	0.626608436922279,	0.622140865670282,	0.617651040393254,	0.613139099020764,	0.60860518080042,	0.604049426290961,	0.599471977355251,	0.594872977153151,	0.590252570134298,	0.585610902030762,	0.580948119849605,	0.576264371865328,	0.571559807612214,	0.566834577876565,	0.562088834688834,	0.557322731315646,	0.552536422251729,	0.547730063211724,	0.542903811121907,	0.538057824111798,	0.533192261505674,	0.528307283813981,	0.523403052724639,	0.518479731094255,	0.513537482939235,	0.508576473426789,	0.503596868865855,	0.498598836697907,	0.493582545487678,	0.488548164913788,	0.483495865759268,	0.478425819901998,	0.473338200305048,	0.468233181006925,	0.463110937111732,	0.457971644779231,	0.45281548121482,	0.447642624659413,	0.442453254379241,	0.437247550655558,	0.432025694774259,	0.426787869015418,	0.421534256642734,	0.416265041892892,	0.41098040996485,	0.405680547009025,	0.400365640116419,	0.39503587730764,	0.389691447521864,	0.384332540605705,	0.378959347302007,	0.373572059238563,	0.368170868916756,	0.362755969700123,	0.357327555802841,	0.351885822278148,	0.346430965006682,	0.340963180684751,	0.335482666812536,	0.329989621682221,	0.324484244366047,	0.318966734704315,	0.313437293293302,	0.307896121473129,	0.30234342131555,	0.296779395611687,	0.291204247859696,	0.285618182252375,	0.280021403664707,	0.274414117641352,	0.268796530384066,	0.263168848739079,	0.257531280184402,	0.25188403281709,	0.246227315340444,	0.24056133705116,	0.234886307826432,	0.229202438110994,	0.223509938904122,	0.217809021746583,	0.212099898707534,	0.206382782371382,	0.200657885824586,	0.194925422642432,	0.18918560687575,	0.183438653037594,	0.177684776089887,	0.171924191430019,	0.166157114877408,	0.160383762660029,	0.154604351400898,	0.148819098104532,	0.143028220143366,	0.137231935244143,	0.131430461474268,	0.125624017228142,	0.119812821213451,	0.113997092437445,	0.108177050193176,	0.10235291404572,	0.0965249038183706,	0.0906932395788124,	0.0848581416252698,	0.0790198304726378,	0.0731785268385933,	0.0673344516296877,	0.0614878259274234,	0.0556388709743148,	0.0497878081599346,	0.0439348590069479,	0.0380802451571339,	0.0322241883573973,	0.0263669104457702,	0.0205086333374064,	0.0146495790105684,	0.00878996949260909,	0.00293002684594913,	-0.00293002684594911,	-0.00878996949260907,	-0.0146495790105683,	-0.0205086333374063,	-0.0263669104457702,	-0.0322241883573973,	-0.0380802451571339,	-0.0439348590069479,	-0.0497878081599346,	-0.0556388709743147,	-0.0614878259274234,	-0.0673344516296877,	-0.0731785268385933,	-0.0790198304726378,	-0.0848581416252697,	-0.0906932395788124,	-0.0965249038183705,	-0.10235291404572,	-0.108177050193176,	-0.113997092437445,	-0.119812821213451,	-0.125624017228142,	-0.131430461474268,	-0.137231935244143,	-0.143028220143366,	-0.148819098104532,	-0.154604351400898,	-0.160383762660029,	-0.166157114877408,	-0.171924191430019,	-0.177684776089887,	-0.183438653037594,	-0.18918560687575,	-0.194925422642432,	-0.200657885824586,	-0.206382782371382,	-0.212099898707534,	-0.217809021746583,	-0.223509938904122,	-0.229202438110994,	-0.234886307826432,	-0.24056133705116,	-0.246227315340444,	-0.25188403281709,	-0.257531280184402,	-0.263168848739079,	-0.268796530384066,	-0.274414117641352,	-0.280021403664707,	-0.285618182252375,	-0.291204247859696,	-0.296779395611687,	-0.30234342131555,	-0.307896121473129,	-0.313437293293302,	-0.318966734704315,	-0.324484244366047,	-0.329989621682221,	-0.335482666812536,	-0.340963180684751,	-0.346430965006681,	-0.351885822278148,	-0.357327555802841,	-0.362755969700123,	-0.368170868916756,	-0.373572059238563,	-0.378959347302007,	-0.384332540605705,	-0.389691447521864,	-0.39503587730764,	-0.400365640116418,	-0.405680547009025,	-0.41098040996485,	-0.416265041892892,	-0.421534256642734,	-0.426787869015418,	-0.432025694774259,	-0.437247550655558,	-0.442453254379241,	-0.447642624659413,	-0.45281548121482,	-0.457971644779231,	-0.463110937111732,	-0.468233181006925,	-0.473338200305048,	-0.478425819901998,	-0.483495865759268,	-0.488548164913788,	-0.493582545487678,	-0.498598836697907,	-0.503596868865855,	-0.508576473426789,	-0.513537482939235,	-0.518479731094255,	-0.523403052724639,	-0.528307283813981,	-0.533192261505674,	-0.538057824111798,	-0.542903811121907,	-0.547730063211724,	-0.552536422251729,	-0.557322731315646,	-0.562088834688834,	-0.566834577876565,	-0.571559807612214,	-0.576264371865327,	-0.580948119849605,	-0.585610902030762,	-0.590252570134298,	-0.594872977153151,	-0.599471977355251,	-0.604049426290961,	-0.60860518080042,	-0.613139099020764,	-0.617651040393253,	-0.622140865670282,	-0.626608436922278,	-0.6310536175445,	-0.635476272263715,	-0.639876267144773,	-0.644253469597067,	-0.648607748380882,	-0.652938973613631,	-0.657247016775979,	-0.661531750717861,	-0.665793049664373,	-0.670030789221564,	-0.674244846382104,	-0.678435099530846,	-0.682601428450267,	-0.686743714325798,	-0.690861839751038,	-0.694955688732854,	-0.699025146696365,	-0.703070100489807,	-0.70709043838929,	-0.711086050103429,	-0.715056826777867,	-0.719002660999677,	-0.722923446801648,	-0.726819079666455,	-0.730689456530709,	-0.734534475788896,	-0.738354037297192,	-0.742148042377162,	-0.745916393819348,	-0.749658995886725,	-0.753375754318056,	-0.757066576331117,	-0.76073137062581,	-0.764370047387153,	-0.767982518288158,	-0.771568696492588,	-0.775128496657593,	-0.778661834936232,	-0.782168628979876,	-0.785648797940494,	-0.789102262472815,	-0.792528944736378,	-0.795928768397466,	-0.799301658630911,	-0.802647542121796,	-0.805966347067026,	-0.80925800317679,	-0.812522441675903,	-0.815759595305026,	-0.818969398321778,	-0.82215178650172,	-0.825306697139235,	-0.828434069048276,	-0.831533842563007,	-0.834605959538332,	-0.83765036335029,	-0.840666998896357,	-0.843655812595612,	-0.846616752388798,	-0.84954976773827,	-0.852454809627816,	-0.855331830562377,	-0.858180784567642,	-0.861001627189535,	-0.863794315493586,	-0.866558808064187,	-0.869295065003735,	-0.872003047931669,	-0.874682719983378,	-0.877334045809019,	-0.879956991572201,	-0.882551524948574,	-0.885117615124298,	-0.887655232794406,	-0.890164350161051,	-0.892644940931649,	-0.89509698031691,	-0.897520445028759,	-0.89991531327815,	-0.902281564772771,	-0.904619180714645,	-0.906928143797612,	-0.909208438204722,	-0.911460049605506,	-0.913682965153153,	-0.915877173481571,	-0.918042664702353,	-0.920179430401633,	-0.922287463636841,	-0.924366758933352,	-0.926417312281038,	-0.928439121130713,	-0.930432184390478,	-0.932396502421965,	-0.934332077036483,	-0.936238911491064,	-0.938117010484405,	-0.93996638015272,	-0.941787028065485,	-0.9435789632211,	-0.945342196042432,	-0.94707673837229,	-0.948782603468781,	-0.950459806000585,	-0.952108362042131,	-0.953728289068685,	-0.955319605951339,	-0.95688233295191,	-0.958416491717755,	-0.959922105276486,	-0.9613991980306,	-0.962847795752024,	-0.964267925576561,	-0.965659615998264,	-0.967022896863706,	-0.968357799366183,	-0.969664356039813,	-0.970942600753568,	-0.972192568705208,	-0.973414296415142,	-0.974607821720205,	-0.97577318376735,	-0.976910423007261,	-0.978019581187894,	-0.979100701347927,	-0.980153827810139,	-0.981179006174716,	-0.982176283312467,	-0.983145707357983,	-0.984087327702705,	-0.985001194987926,	-0.985887361097722,	-0.986745879151805,	-0.987576803498309,	-0.988380189706504,	-0.989156094559442,	-0.989904576046531,	-0.990625693356047,	-0.991319506867572,	-0.991986078144371,	-0.992625469925706,	-0.993237746119077,	-0.993822971792408,	-0.994381213166167,	-0.994912537605426,	-0.995417013611858,	-0.995894710815678,	-0.996345699967518,	-0.996770052930258,	-0.997167842670781,	-0.997539143251688,	-0.997884029822952,	-0.998202578613513,	-0.998494866922831,	-0.998760973112374,	-0.999000976597064,	-0.999214957836666,	-0.999402998327134,	-0.999565180591902,	-0.999701588173131,	-0.999812305622909,	-0.999897418494402,	-0.999957013332964,	-0.999991177667203,	-0.999999999999999,	-0.999983569799488,	-0.999941977489996,	-0.999875314442939,	-0.999783672967683,	-0.99966714630236,	-0.999525828604659,	-0.999359814942567,	-0.999169201285081,	-0.998954084492887,	-0.998714562309002,	-0.998450733349385,	-0.998162697093516,	-0.997850553874947,	-0.997514404871815,	-0.997154352097337,	-0.996770498390271,	-0.99636294740535,	-0.99593180360369,	-0.995477172243177,	-0.994999159368825,	-0.994497871803115,	-0.993973417136306,	-0.993425903716732,	-0.992855440641073,	-0.992262137744609,	-0.991646105591455,	-0.991007455464782,	-0.990346299357011,	-0.989662749960006,	-0.98895692065524,	-0.988228925503953,	-0.987478879237297,	-0.986706897246467,	-0.985913095572823,	-0.985097590897998,	-0.984260500534004,	-0.983401942413322,	-0.982522035078989,	-0.981620897674677,	-0.980698649934767,	-0.97975541217442,	-0.978791305279639,	-0.977806450697332,	-0.976800970425375,	-0.975774987002668,	-0.974728623499197,	-0.973662003506092,	-0.972575251125687,	-0.971468490961584,	-0.970341848108722,	-0.969195448143445,	-0.968029417113579,	-0.966843881528513,	-0.965638968349289,	-0.964414804978695,	-0.963171519251369,	-0.961909239423917,	-0.960628094165029,	-0.959328212545621,	-0.958009724028973,	-0.956672758460894,	-0.95531744605989,	-0.95394391740735,	-0.952552303437747,	-0.951142735428856,	-0.949715344991984,	-0.948270264062223,	-0.946807624888716,	-0.945327560024946,	-0.943830202319042,	-0.942315684904107,	-0.940784141188564,	-0.939235704846531,	-0.937670509808209,	-0.936088690250302,	-0.934490380586457,	-0.932875715457728,	-0.931244829723069,	-0.929597858449852,	-0.927934936904411,	-0.926256200542612,	-0.924561785000459,	-0.92285182608472,	-0.92112645976359,	-0.91938582215738,	-0.917630049529239,	-0.915859278275909,	-0.91407364491851,	-0.912273286093362,	-0.910458338542837,	-0.908628939106248,	-0.906785224710769,	-0.904927332362403,	-0.903055399136968,	-0.901169562171136,	-0.899269958653503,	-0.897356725815695,	-0.895430000923519,	-0.893489921268151,	-0.891536624157358,	-0.889570246906774,	-0.887590926831205,	-0.88559880123598,	-0.883594007408349,	-0.881576682608912,	-0.87954696406311,	-0.877504988952739,	-0.875450894407525,	-0.873384817496738,	-0.871306895220847,	-0.869217264503232,	-0.867116062181932,	-0.865003425001444,	-0.862879489604577,	-0.860744392524338,	-0.858598270175882,	-0.856441258848504,	-0.854273494697681,	-0.852095113737163,	-0.849906251831119,	-0.847707044686327,	-0.845497627844425,	-0.8432781366742,	-0.841048706363942,	-0.838809471913847,	-0.836560568128465,	-0.834302129609215,	-0.83203429074694,	-0.829757185714528,	-0.827470948459579,	-0.825175712697131,	-0.82287161190244,	-0.820558779303817,	-0.818237347875515,	-0.815907450330684,	-0.813569219114369,	-0.811222786396577,	-0.808868284065389,	-0.806505843720144,	-0.804135596664664,	-0.801757673900554,	-0.799372206120547,	-0.796979323701917,	-0.794579156699943,	-0.792171834841439,	-0.78975748751834,	-0.787336243781347,	-0.784908232333635,	-0.782473581524619,	-0.78003241934378,	-0.777584873414552,	-0.775131070988273,	-0.772671138938191,	-0.770205203753537,	-0.767733391533655,	-0.765255827982196,	-0.762772638401375,	-0.760283947686285,	-0.75778988031928,	-0.755290560364413,	-0.752786111461942,	-0.750276656822896,	-0.747762319223706,	-0.745243221000891,	-0.74271948404582,	-0.740191229799528,	-0.737658579247591,	-0.73512165291508,	-0.73258057086156,	-0.730035452676167,	-0.727486417472737,	-0.724933583885008,	-0.722377070061881,	-0.719816993662742,	-0.717253471852855,	-0.714686621298812,	-0.71211655816405,	-0.709543398104432,	-0.706967256263889,	-0.70438824727013,	-0.701806485230412,	-0.699222083727379,	-0.696635155814956,	-0.694045814014317,	-0.691454170309909,	-0.688860336145546,	-0.686264422420561,	-0.683666539486026,	-0.681066797141033,	-0.67846530462904,	-0.675862170634283,	-0.673257503278246,	-0.670651410116197,	-0.668043998133789,	-0.665435373743726,	-0.662825642782485,	-0.660214910507105,	-0.657603281592047,	-0.654990860126099,	-0.652377749609364,	-0.649764052950294,	-0.647149872462798,	-0.644535309863408,	-0.641920466268504,	-0.639305442191611,	-0.636690337540744,	-0.634075251615831,	-0.631460283106184,	-0.62884553008804,	-0.626231090022159,	-0.623617059751487,	-0.621003535498876,	-0.61839061286487,	-0.615778386825546,	-0.613166951730421,	-0.610556401300414,	-0.607946828625875,	-0.605338326164666,	-0.602730985740309,	-0.600124898540189,	-0.597520155113818,	-0.594916845371158,	-0.592315058581004,	-0.589714883369422,	-0.58711640771825,	-0.584519718963656,	-0.581924903794749,	-0.57933204825226,	-0.576741237727261,	-0.574152556959964,	-0.571566090038557,	-0.568981920398111,	-0.566400130819535,	-0.563820803428589,	-0.561244019694954,	-0.558669860431362,	-0.556098405792769,	-0.553529735275596,	-0.550963927717017,	-0.548401061294305,	-0.54584121352423,	-0.543284461262513,	-0.540730880703331,	-0.538180547378877,	-0.535633536158973,	-0.533089921250733,	-0.530549776198285,	-0.528013173882533,	-0.525480186520986,	-0.522950885667625,	-0.520425342212826,	-0.517903626383337,	-0.515385807742302,	-0.51287195518933,	-0.510362136960626,	-0.507856420629158,	-0.505354873104884,	-0.502857560635018,	-0.500364548804351,	-0.497875902535619,	-0.495391686089915,	-0.492911963067151,	-0.490436796406567,	-0.487966248387282,	-0.485500380628902,	-0.483039254092159,	-0.480582929079607,	-0.478131465236355,	-0.475684921550853,	-0.473243356355713,	-0.470806827328582,	-0.468375391493053,	-0.465949105219618,	-0.463528024226672,	-0.461112203581549,	-0.458701697701608,	-0.456296560355356,	-0.453896844663611,	-0.451502603100714,	-0.449113887495767,	-0.446730749033927,	-0.444353238257729,	-0.44198140506845,	-0.439615298727515,	-0.437254967857937,	-0.434900460445799,	-0.432551823841771,	-0.430209104762667,	-0.427872349293032,	-0.425541602886773,	-0.423216910368827,	-0.420898315936851,	-0.418585863162967,	-0.416279594995523,	-0.413979553760904,	-0.411685781165364,	-0.4093983182969,	-0.407117205627156,	-0.404842483013362,	-0.402574189700297,	-0.400312364322301,	-0.398057044905295,	-0.395808268868855,	-0.393566073028302,	-0.391330493596828,	-0.389101566187648,	-0.386879325816189,	-0.384663806902297,	-0.382455043272484,	-0.380253068162194,	-0.378057914218103,	-0.375869613500444,	-0.373688197485355,	-0.371513697067265,	-0.36934614256129,	-0.367185563705671,	-0.365031989664226,	-0.36288544902883,	-0.360745969821925,	-0.358613579499042,	-0.356488304951361,	-0.354370172508285,	-0.352259207940038,	-0.350155436460287,	-0.34805888272879,	-0.345969570854057,	-0.343887524396038,	-0.341812766368834,	-0.339745319243421,	-0.337685204950401,	-0.335632444882772,	-0.333587059898712,	-0.331549070324387,	-0.329518495956777,	-0.327495356066518,	-0.325479669400764,	-0.323471454186064,	-0.321470728131258,	-0.319477508430391,	-0.317491811765638,	-0.315513654310247,	-0.313543051731504,	-0.311580019193701,	-0.30962457136113,	-0.307676722401083,	-0.305736485986872,	-0.303803875300857,	-0.301878903037493,	-0.299961581406385,	-0.298051922135358,	-0.296149936473537,	-0.294255635194444,	-0.292369028599095,	-0.290490126519124,	-0.288618938319902,	-0.286755472903678,	-0.284899738712721,	-0.28305174373248,	-0.281211495494743,	-0.279379001080817,	-0.277554267124704,	-0.275737299816299,	-0.273928104904577,	-0.27212668770081,	-0.270333053081769,	-0.268547205492953,	-0.266769148951807,	-0.264998887050958,	-0.263236422961448,	-0.261481759435985,	-0.259734898812182,	-0.257995843015815,	-0.256264593564078,	-0.254541151568847,	-0.25282551773994,	-0.251117692388391,	-0.249417675429719,	-0.2477254663872,	-0.246041064395149,	-0.244364468202196,	-0.242695676174569,	-0.241034686299374,	-0.23938149618788,	-0.237736103078808,	-0.23609850384161,	-0.234468694979759,	-0.232846672634035,	-0.231232432585813,	-0.229625970260341,	-0.228027280730036,	-0.226436358717757,	-0.224853198600095,	-0.223277794410651,	-0.221710139843315,	-0.220150228255542,	-0.21859805267163,	-0.217053605785988,	-0.215516879966405,	-0.213987867257319,	-0.212466559383077,	-0.210952947751191,	-0.209447023455598,	-0.207948777279906,	-0.206458199700641,	-0.204975280890487,	-0.203500010721524,	-0.202032378768459,	-0.200572374311847,	-0.199119986341317,	-0.197675203558782,	-0.196238014381646,	-0.194808406946008,	-0.193386369109853,	-0.191971888456245,	-0.190564952296501,	-0.189165547673371,	-0.1877736613642,	-0.186389279884085,	-0.18501238948903,	-0.183642976179086,	-0.182281025701481,	-0.180926523553753,	-0.17957945498686,	-0.178239805008289,	-0.176907558385159,	-0.175582699647304,	-0.174265213090357,	-0.172955082778818,	-0.171652292549114,	-0.170356826012651,	-0.169068666558853,	-0.167787797358189,	-0.166514201365195,	-0.165247861321484,	-0.163988759758735,	-0.162736879001687,	-0.161492201171113,	-0.160254708186779,	-0.159024381770398,	-0.157801203448572,	-0.15658515455572,	-0.155376216236989,	-0.154174369451167,	-0.152979594973564,	-0.151791873398903,	-0.150611185144176,	-0.149437510451504,	-0.148270829390976,	-0.14711112186348,	-0.145958367603511,	-0.14481254618198,	-0.143673637008997,	-0.142541619336651,	-0.141416472261764,	-0.140298174728645,	-0.139186705531818,	-0.138082043318745,	-0.136984166592532,	-0.135893053714617,	-0.134808682907446,	-0.133731032257143,	-0.132660079716148,	-0.131595803105856,	-0.130538180119238,	-0.129487188323437,	-0.128442805162365,	-0.12740500795927,	-0.126373773919303,	-0.125349080132053,	-0.124330903574082,	-0.123319221111436,	-0.122314009502142,	-0.121315245398692,	-0.12032290535051,	-0.119336965806402,	-0.118357403116992,	-0.117384193537142,	-0.116417313228359,	-0.115456738261176,	-0.114502444617534,	-0.113554408193129,	-0.11261260479976,	-0.111677010167649,	-0.11074759994775,	-0.109824349714041,	-0.108907234965801,	-0.107996231129867,	-0.107091313562881,	-0.106192457553512,	-0.10529963832467,	-0.104412831035699,	-0.103532010784556,	-0.102657152609969,	-0.101788231493583,	-0.100925222362088,	-0.100068100089333,	-0.0992168394984166,	-0.0983714153637673,	-0.0975318024132052,	-0.0966979753299867,	-0.0958699087548324,	-0.0950475772879391,	-0.0942309554909746,	-0.0934200178890553,	-0.0926147389727083,	-0.0918150931998155,	-0.0910210549975421,	-0.0902325987642472,	-0.0894496988713785,	-0.08867232966535,	-0.0879004654694031,	-0.0871340805854505,	-0.0863731492959036,	-0.0856176458654837,	-0.0848675445430155,	-0.0841228195632047,	-0.0833834451483983,	-0.0826493955103288,	-0.0819206448518412,	-0.0811971673686034,	-0.0804789372508005,	-0.0797659286848117,	-0.0790581158548709,	-0.0783554729447112,	-0.0776579741391918,	-0.0769655936259095,	-0.0762783055967927,	-0.0755960842496797,	-0.0749189037898797,	-0.0742467384317184,	-0.0735795624000665,	-0.0729173499318516,	-0.0722600752775546,	-0.0716077127026893,	-0.0709602364892655,	-0.0703176209372364,	-0.0696798403659295,	-0.0690468691154616,	-0.0684186815481373,	-0.0677952520498318,	-0.0671765550313574,	-0.0665625649298145,	-0.0659532562099261,	-0.0653486033653564,	-0.0647485809200147,	-0.0641531634293416,	-0.0635623254815814,	-0.0629760416990374,	-0.0623942867393123,	-0.0618170352965328,	-0.061244262102559,	-0.0606759419281775,	-0.0601120495842801,	-0.0595525599230266,	-0.0589974478389922,	-0.0584466882703003,	-0.0579002561997392,	-0.0573581266558646,	-0.0568202747140862,	-0.0562866754977403,	-0.0557573041791465,	-0.0552321359806501,	-0.0547111461756496,	-0.0541943100896095,	-0.0536816031010583,	-0.053173000642572,	-0.0526684782017432,	-0.0521680113221352,	-0.0516715756042228,	-0.0511791467063178,	-0.0506907003454802,	-0.0502062122984163,	-0.0497256584023615,	-0.0492490145559497,	-0.0487762567200688,	-0.048307360918702,	-0.0478423032397555,	-0.0473810598358723,	-0.0469236069252326,	-0.0464699207923407,	-0.0460199777887977,	-0.0455737543340617,	-0.0451312269161942,	-0.0446923720925932,	-0.0442571664907131,	-0.0438255868087719,	-0.0433976098164448,	-0.0429732123555452,	-0.0425523713406931,	-0.0421350637599698,	-0.041721266675561,	-0.0413109572243868,	-0.0409041126187188,	-0.0405007101467854,	-0.0401007271733643,	-0.0397041411403629,	-0.0393109295673867,	-0.038921070052295,	-0.0385345402717451,	-0.0381513179817244,	-0.0377713810180704,	-0.0373947072969797,	-0.0370212748155041,	-0.0366510616520362,	-0.0362840459667827,	-0.0359202060022271,	-0.0355595200835802,	-0.0352019666192197,	-0.0348475241011187,	-0.0344961711052634,	-0.0341478862920588,	-0.0338026484067248,	-0.0334604362796805,	-0.0331212288269181,	-0.0327850050503663,	-0.0324517440382427,	-0.0321214249653966,	-0.0317940270936401,	-0.0314695297720703,	-0.0311479124373802,	-0.0308291546141601,	-0.0305132359151889,	-0.0302001360417151,	-0.0298898347837281,	-0.0295823120202204,	-0.0292775477194389,	-0.0289755219391274,	-0.0286762148267598,	-0.0283796066197631,	-0.0280856776457316,	-0.0277944083226319,	-0.0275057791589985,	-0.0272197707541204,	-0.0269363637982188,	-0.0266555390726158,	-0.0263772774498946,	-0.0261015598940504,	-0.0258283674606333,	-0.025557681296882,	-0.02528948264185,	-0.0250237528265224,	-0.0247604732739247,	-0.0244996254992243,	-0.0242411911098222,	-0.0239851518054385,	-0.0237314893781885,	-0.0234801857126519,	-0.0232312227859334,	-0.0229845826677168,	-0.02274024752031,	-0.0224981995986839,	-0.0222584212505028,	-0.0220208949161481,	-0.0217856031287346,	-0.0215525285141195,	-0.0213216537909041,	-0.0210929617704291,	-0.0208664353567625,	-0.0206420575466809,	-0.0204198114296433,	-0.0201996801877599,	-0.0199816470957523,	-0.0197656955209084,	-0.0195518089230307,	-0.019339970854378,	-0.0191301649596011,	-0.0189223749756719,	-0.0187165847318066,	-0.0185127781493833,	-0.0183109392418527,	-0.018111052114644,	-0.0179131009650639,	-0.0177170700821908,	-0.017522943846763,	-0.017330706731061,	-0.017140343298785,	-0.0169518382049263,	-0.0167651761956341,	-0.0165803421080761,	-0.0163973208702947,	-0.016216097501058,	-0.016036657109705,	-0.0158589848959873,	-0.0156830661499042,	-0.0155088862515344,	-0.0153364306708618,	-0.0151656849675979,	-0.0149966347909985,	-0.0148292658796763,	-0.014663564061409,	-0.0144995152529435,	-0.0143371054597946,	-0.014176320776041,	-0.0140171473841156,	-0.013859571554593,	-0.0137035796459721,	-0.0135491581044551,	-0.0133962934637227,	-0.0132449723447048,	-0.0130951814553487,	-0.0129469075903822,	-0.0128001376310741,	-0.0126548585449906,	-0.0125110573857484,	-0.0123687212927644,	-0.0122278374910017,	-0.0120883932907129,	-0.0119503760871796,	-0.0118137733604491,	-0.0116785726750676,	-0.0115447616798107,	-0.0114123281074108,	-0.0112812597742815,	-0.0111515445802392,	-0.0110231705082218,	-0.0108961256240047,	-0.0107703980759144,	-0.0106459760945388,	-0.0105228479924357,	-0.0104010021638379,	-0.0102804270843569,	-0.0101611113106834,	-0.0100430434802856,	-0.00992621231110553,	-0.00981060660125272,	-0.00969621522869606,	-0.00958302715095331,	-0.00947103140477862,	-0.00936021710584801,	-0.00925057344844294,	-0.00914208970513188,	-0.00903475522645012,	-0.00892855944057768,	-0.00882349185301557,	-0.00871954204626023,	-0.00861669967947642,	-0.00851495448816841,	-0.00841429628384972,	-0.00831471495371128,	-0.00821620046028813,	-0.00811874284112478,	-0.00802233220843912,	-0.00792695874878505,	-0.00783261272271391,	-0.00773928446443455,	-0.00764696438147236,	-0.00755564295432705,	-0.00746531073612947,	-0.00737595835229724,	-0.00728757650018951,	-0.0072001559487607,	-0.00711368753821333,	-0.00702816217965,	-0.00694357085472451,	-0.00685990461529223,	-0.0067771545830597,	-0.0066953119492335,	-0.00661436797416853,	-0.00653431398701559,	-0.0064551413853684,	-0.00637684163491005,	-0.006299406269059,	-0.00622282688861451,	-0.00614709516140176,	-0.00607220282191644,	-0.00599814167096905,	-0.0059249035753289,	-0.00585248046736776,	-0.00578086434470324,	-0.00571004726984205,	-0.00564002136982297,	-0.00557077883585971,	-0.00550231192298367,	-0.00543461294968656,	-0.00536767429756302,	-0.00530148841095317,	-0.00523604779658524,	-0.00517134502321815,	-0.00510737272128423,	-0.00504412358253202,	-0.0049815903596692,	-0.00491976586600569,	-0.00485864297509694,	-0.00479821462038741,	-0.00473847379485438,	-0.00467941355065192,	-0.00462102699875526,	-0.00456330730860547,	-0.00450624770775444,	-0.00444984148151032,	-0.00439408197258336,	-0.00433896258073211,	-0.00428447676241019,	-0.00423061803041345,	-0.00417737995352772,	-0.00412475615617702,	-0.00407274031807241,	-0.00402132617386134,	-0.00397050751277765,	-0.0039202781782922,	-0.00387063206776413,	-0.0038215631320928,	-0.00377306537537041,	-0.00372513285453537,	-0.00367775967902635,	-0.00363094001043713,	-0.00358466806217222,	-0.00353893809910324,	-0.00349374443722614,	-0.00344908144331925,	-0.00340494353460221,	-0.00336132517839569,	-0.00331822089178209,	-0.0032756252412671,	-0.00323353284244214,	-0.00319193835964786,	-0.00315083650563847,	-0.00311022204124709,	-0.00307008977505213,	-0.00303043456304464,	-0.00299125130829662,	-0.00295253496063049,	-0.0029142805162895,	-0.00287648301760926,	-0.00283913755269034,	-0.00280223925507193,	-0.00276578330340667,	-0.00272976492113652,	-0.00269417937616984,	-0.00265902198055954,	-0.00262428809018245,	-0.0025899731044198,	-0.00255607246583894,	-0.0025225816598762,	-0.00248949621452101,	-0.00245681170000115,	-0.00242452372846931,	-0.00239262795369083,	-0.00236112007073272,	-0.00232999581565392,	-0.0022992509651968,	-0.00226888133647996,	-0.00223888278669236,	-0.00220925121278862,	-0.00217998255118574,	-0.00215107277746107,	-0.00212251790605159,	-0.00209431398995457,	-0.00206645712042947,	-0.00203894342670129,	-0.00201176907566516,	-0.00198493027159237,	-0.00195842325583771,	-0.00193224430654817,	-0.00190638973837306,	-0.00188085590217545,	-0.00185563918474507,	-0.0018307360085125,	-0.00180614283126484,	-0.00178185614586274,	-0.00175787247995882,	-0.00173418839571759,	-0.00171080048953662,	-0.0016877053917693,	-0.00166489976644891,	-0.00164238031101417,	-0.00162014375603616,	-0.00159818686494679,	-0.00157650643376855,	-0.00155509929084582,	-0.00153396229657758,	-0.00151309234315155,	-0.00149248635427982,	-0.00147214128493586,	-0.00145205412109308,	-0.00143222187946472,	-0.00141264160724534,	-0.0013933103818536,	-0.00137422531067669,	-0.00135538353081601,	-0.00133678220883452,	-0.00131841854050538,	-0.00130028975056213,	-0.00128239309245038,	-0.00126472584808083,	-0.00124728532758389,	-0.00123006886906569,	-0.00121307383836556,	-0.001196297628815,	-0.0011797376609981,	-0.00116339138251339,	-0.00114725626773724,	-0.00113132981758859,	-0.0011156095592953,	-0.00110009304616182,	-0.00108477785733839,	-0.00106966159759168,	-0.00105474189707691,	-0.00104001641111138,	-0.00102548281994949,	-0.00101113882855921,	-0.000996982166399964,	-0.000983010587201993,	-0.000969221868747178,	-0.000955613812651267,	-0.000942184244147566,	-0.000928931011872069,	-0.000915851987650011,	-0.000902945066283871,	-0.000890208165342797,	-0.000877639224953454,	-0.000865236207592308,	-0.000852997097879327,	-0.000840919902373096,	-0.000829002649367355,	-0.000817243388688945,	-0.000805640191497166,	-0.000794191150084532,	-0.00078289437767894,	-0.000771748008247229,	-0.00076075019630013,	-0.000749899116698618,	-0.000739192964461637,	-0.000728629954575219,	-0.000718208321802972,	-0.000707926320497948,	-0.000697782224415878,	-0.000687774326529774,	-0.000677900938845889,	-0.00066816039222104,	-0.000658551036181277,	-0.000649071238741901,	-0.000639719386228828,	-0.000630493883101285,	-0.000621393151775847,	-0.000612415632451795,	-0.000603559782937805,	-0.000594824078479947,	-0.000586207011591006,	-0.000577707091881107,	-0.000569322845889642,	-0.000561052816918502,	-0.000552895564866591,	-0.000544849666065639,	-0.000536913713117291,	-0.000529086314731472,	-0.00052136609556603,	-0.000513751696067636,	-0.000506241772313954,	-0.000498834995857058,	-0.000491530053568103,	-0.000484325647483238,	-0.000477220494650759,	-0.000470213326979492,	-0.000463302891088402,	-0.000456487948157427,	-0.000449767273779525,	-0.000443139657813925,	-0.000436603904240592,	-0.00043015883101588,	-0.00042380326992938,	-0.000417536066461955,	-0.000411356079644951,	-0.000405262181920587,	-0.000399253259003505,	-0.000393328209743488,	-0.00038748594598933,	-0.000381725392453856,	-0.000376045486580087,	-0.000370445178408539,	-0.000364923430445654,	-0.000359479217533354,	-0.000354111526719721,	-0.000348819357130778,	-0.000343601719843384,	-0.000338457637759227,	-0.000333386145479908,	-0.000328386289183118,	-0.000323457126499887,	-0.00031859772639292,	-0.000313807169035989,	-0.000309084545694398,	-0.000304428958606493,	-0.000299839520866232,	-0.000295315356306795,	-0.000290855599385224,	-0.000286459395068104,	-0.000282125898718263,	-0.000277854275982491,	-0.000273643702680267,	-0.000269493364693501,	-0.000265402457857263,	-0.000261370187851514,	-0.000257395770093819,	-0.00025347842963304,	-0.000249617401044004,	-0.000245811928323141,	-0.000242061264785073,	-0.000238364672960175,	-0.000234721424493065,	-0.000231130800042051,	-0.000227592089179507,	-0.000224104590293179,	-0.000220667610488419,	-0.000217280465491327,	-0.000213942479552814,	-0.000210652985353562,	-0.000207411323909887,	-0.000204216844480491,	-0.000201068904474101,	-0.000197966869357987,	-0.000194910112567358,	-0.000191898015415616,	-0.00018892996700548,	-0.000186005364140961,	-0.000183123611240183,	-0.000180284120249058,	-0.000177486310555779,	-0.000174729608906162,	-0.000172013449319802,	-0.000169337273007048,	-0.000166700528286794,	-0.000164102670505074,	-0.000161543161954455,	-0.000159021471794232,	-0.0001565370759714,	-0.000154089457142417,	-0.000151678104595744,	-0.000149302514175145,	-0.000146962188203764,	-0.000144656635408953,	-0.000142385370847855,	-0.000140147915833738,	-0.00013794379786306,	-0.000135772550543278,	-0.000133633713521389,	-0.00013152683241318,	-0.000129451458733216,	-0.000127407149825525,	-0.000125393468794996,	-0.000123409984439478,	-0.000121456271182573,	-0.00011953190900712,	-0.000117636483389361,	-0.00011576958523379,	-0.000113930810808671,	-0.000112119761682222,	-0.000110336044659469,	-0.000108579271719748,	-0.000106849059954864,	-0.00010514503150789,	-0.000103466813512615,	-0.000101814038033616,	-0.000100186342006972,	-9.85833671815882E-05,	-9.70047600611532E-05,	-9.54501718467006E-05,	-9.39192583797826E-05,	-9.24116800862466E-05,	-9.09271019206098E-05,	-8.94651933110269E-05,	-8.80256281048454E-05,	-8.66080845147441E-05,	-8.52122450654486E-05,	-8.38377965410207E-05,	-8.24844299327136E-05,	-8.11518403873915E-05,	-7.98397271565059E-05,	-7.85477935456249E-05,	-7.72757468645101E-05,	-7.60232983777369E-05,	-7.47901632558526E-05,	-7.35760605270676E-05,	-7.23807130294755E-05,	-7.12038473637968E-05,	-7.00451938466412E-05,	-6.89044864642845E-05,	-6.77814628269553E-05,	-6.66758641236261E-05,	-6.55874350773051E-05,	-6.45159239008233E-05,	-6.34610822531134E-05,	-6.24226651959747E-05,	-6.14004311513202E-05,	-6.03941418589015E-05,	-0.000059403562334507,	-5.84284608286282E-05,	-0.000057468608785591,	-5.65237808031466E-05,	-0.000055593754592518,	-5.46783109388983E-05,	-5.37772336623954E-05,	-0.00005289030957942,	-5.20173284645116E-05,	-5.11580830125996E-05,	-0.000050312368801694,	-4.94799842560025E-05,	-4.86607306094695E-05,	-4.78544118697325E-05,	-4.70608347824935E-05,	-4.62798087962986E-05,	-0.000045511146027725,	-4.47546612269687E-05,	-4.40101717438305E-05,	-4.32774974940962E-05,	-4.25564609263065E-05,	-4.18468869889135E-05,	-4.11486030978194E-05,	-4.04614391042941E-05,	-3.97852272632673E-05,	-0.000039119802201992,	-3.84650008890754E-05,	-3.78206626038733E-05,	-3.71866289062449E-05,	-3.65627436066639E-05,	-3.59488527366823E-05,	-3.53448045197434E-05,	-3.47504493423413E-05,	-3.41656397255213E-05,	-3.35902302967199E-05,	-3.30240777619402E-05,	-3.24670408782583E-05,	-3.19189804266585E-05,	-0.000031379759185194,	-3.08492419024682E-05,	-3.03272952714353E-05,	-2.98137879035152E-05,	-2.93085903030209E-05,	-2.88115748418934E-05,	-2.83226157347427E-05,	-2.78415890141905E-05,	-2.73683725065118E-05,	-2.69028458075724E-05,	-2.64448902590588E-05,	-0.000025994388924998,	-2.55512265685631E-05,	-2.51152896291632E-05,	-0.000024686466199813,	-0.00002426464600478,	-2.38497203775067E-05,	-0.000023441582238803,	-2.30401260753087E-05,	-0.000022645247918221,	-2.22568453222843E-05,	-2.18748173450411E-05,	-2.14990645263393E-05,	-2.11294888680941E-05,	-2.07659938143015E-05,	-2.04084842313007E-05,	-2.00568663882829E-05,	-1.97110479380427E-05,	-1.93709378979717E-05,	-1.90364466312886E-05,	-1.87074858285063E-05,	-1.83839684891308E-05,	-1.80658089035912E-05,	-1.77529226353969E-05,	-1.74452265035204E-05,	-1.71426385650029E-05,	-1.68450780977798E-05,	-1.65524655837242E-05,	-0.000016264722691906,	-1.59817722620634E-05,	-1.57035382882852E-05,	-1.54299459029017E-05,	-1.51609213605804E-05,	-1.48963920226262E-05,	-1.46362863414824E-05,	-0.00001438053384543,	-1.41290651234849E-05,	-1.38818118104883E-05,	-1.36387065723899E-05,	-1.33996830917206E-05,	-1.31646760532539E-05,	-1.29336211298513E-05,	-1.27064549684932E-05,	-1.24831151764895E-05,	-1.22635403078701E-05,	-1.20476698499529E-05,	-1.18354442100864E-05,	-1.16268047025659E-05,	-1.14216935357199E-05,	-0.000011220053799167,	-1.10218294512383E-05,	-1.08269653065664E-05,	-1.06354070238366E-05,	-1.04471010936996E-05,	-1.02619948268443E-05,	-1.00800363422269E-05,	-9.90117455545663E-06,	-9.72535916733479E-06,	-9.55254065254605E-06,	-9.38267024849966E-06,	-9.21569994431914E-06,	-9.05158246997867E-06,	-8.89027128558417E-06,	-8.73172057079763E-06,	-8.57588521440278E-06,	-8.42272080401045E-06,	-8.27218361590201E-06,	-8.12423060500914E-06,	-7.97881939502822E-06,	-7.83590826866797E-06,	-7.69545615802838E-06,	-7.55742263510963E-06,	-7.4217679024492E-06,	-7.28845278388577E-06,	-7.15743871544827E-06,	-7.02868773636856E-06,	-6.90216248021624E-06,	-6.77782616615412E-06,	-6.65564259031278E-06,	-6.53557611728282E-06,	-6.4175916717234E-06,	-6.30165473008549E-06,	-6.18773131244848E-06,	-6.07578797446885E-06,	-5.96579179943928E-06,	-5.85771039045709E-06,	-5.75151186270035E-06,	-5.64716483581068E-06,	-0.000005544638426381,	-5.44390224054728E-06,	-5.34492636668272E-06,	-5.24768136819324E-06,	-5.15213827641291E-06,	-5.05826858359812E-06,	-4.96604423601917E-06,	-4.87543762714814E-06,	-4.78642159094177E-06,	-4.69896939521808E-06,	-4.61305473512573E-06,	-4.52865172670468E-06,	-4.44573490053723E-06,	-4.36427919548811E-06,	-4.28425995253257E-06,	-4.20565290867134E-06,	-4.12843419093127E-06,	-4.05258031045068E-06,	-3.9780681566482E-06,	-3.90487499147408E-06,	-3.83297844374292E-06,	-3.76235650354677E-06,	-3.69298751674746E-06,	-3.62485017954731E-06,	-3.55792353313704E-06,	-3.49218695841993E-06,	-3.42762017081129E-06,	-3.36420321511219E-06,	-3.30191646045647E-06,	-3.24074059533014E-06,	-3.18065662266216E-06,	-3.1216458549857E-06,	-3.06368990966892E-06,	-3.00677070421434E-06,	-2.95087045162604E-06,	-2.89597165584351E-06,	-2.84205710724158E-06,	-2.78910987819536E-06,	-2.73711331870935E-06,	-2.68605105210992E-06,	-2.63590697080034E-06,	-2.58666523207738E-06,	-2.53831025400884E-06,	-2.49082671137107E-06,	-2.44419953164573E-06,	-2.39841389107497E-06,	-2.3534552107743E-06,	-2.3093091529023E-06,	-2.26596161688641E-06,	-2.22339873570415E-06,	-2.18160687221886E-06,	-2.14057261556937E-06,	-2.10028277761277E-06,	-2.0607243894196E-06,	-2.02188469782074E-06,	-1.98375116200533E-06,	-1.94631145016892E-06,	-1.90955343621136E-06,	-1.87346519648345E-06,	-1.8380350065821E-06,	-1.80325133819283E-06,	-1.76910285597948E-06,	-1.73557841452002E-06,	-1.70266705528819E-06,	-1.67035800368009E-06,	-1.63864066608524E-06,	-1.60750462700139E-06,	-1.57693964619266E-06,	-1.54693565589012E-06,	-1.5174827580345E-06,	-1.48857122156029E-06,	-1.46019147972069E-06,	-1.43233412745288E-06,	-1.40498991878294E-06,	-1.37814976427002E-06,	-1.35180472848908E-06,	-1.32594602755171E-06,	-1.30056502666452E-06,	-1.27565323772453E-06,	-1.25120231695106E-06,	-1.22720406255361E-06,	-1.20365041243519E-06,	-1.18053344193068E-06,	-1.15784536157959E-06,	-1.1355785149329E-06,	-1.11372537639334E-06,	-1.09227854908874E-06,	-1.07123076277795E-06,	-1.05057487178883E-06,	-1.03030385298792E-06,	-1.01041080378128E-06,	-9.90888940146074E-07,	-9.71731594692495E-07,	-9.52932214755521E-07,	-9.34484360516143E-07,	-9.163817031516E-07,	-8.9861802301423E-07,	-8.81187207838504E-07,	-8.64083250975851E-07,	-8.4730024965687E-07,	-8.30832403280523E-07,	-8.14674011729937E-07,	-7.9881947371441E-07,	-7.83263285137254E-07,	-7.68000037489097E-07,	-7.53024416266268E-07,	-7.38331199413905E-07,	-7.23915255793425E-07,	-7.09771543674002E-07,	-6.95895109247693E-07,	-6.82281085167875E-07,	-6.68924689110654E-07,	-6.55821222358893E-07,	-6.42966068408547E-07,	-6.3035469159696E-07,	-6.17982635752813E-07,	-6.0584552286739E-07,	-5.93939051786867E-07,	-5.8225899692529E-07,	-5.70801206997952E-07,	-5.59561603774865E-07,	-5.48536180854021E-07,	-5.37721002454154E-07,	-5.27112202226713E-07,	-5.16705982086753E-07,	-5.0649861106248E-07,	-4.9648642416314E-07,	-4.86665821265017E-07,	-4.77033266015231E-07,	-4.675852847531E-07,	-4.5831846544878E-07,	-4.49229456658941E-07,	-4.40314966499206E-07,	-4.31571761633117E-07,	-4.22996666277368E-07,	-4.14586561223063E-07,	-4.06338382872753E-07,	-3.98249122293025E-07,	-3.90315824282393E-07,	-3.82535586454263E-07,	-3.74905558334756E-07,	-3.67422940475141E-07,	-3.60084983578672E-07,	-3.52888987641609E-07,	-3.45832301108195E-07,	-3.38912320039386E-07,	-3.32126487295126E-07,	-3.25472291729947E-07,	-3.18947267401703E-07,	-3.12548992793228E-07,	-3.06275090046727E-07,	-3.00123224210693E-07,	-2.94091102499171E-07,	-2.88176473563165E-07,	-2.82377126774014E-07,	-2.76690891518541E-07,	-2.71115636505796E-07,	-2.65649269085222E-07,	-2.60289734576048E-07,	-2.55035015607762E-07,	-2.49883131471461E-07,	-2.44832137481933E-07,	-2.39880124350294E-07,	-2.35025217567015E-07,	-2.30265576795177E-07,	-2.25599395273803E-07,	-2.21024899231089E-07,	-2.16540347307407E-07,	-2.12144029987897E-07,	-2.07834269044527E-07,	-2.03609416987444E-07,	-1.99467856525492E-07,	-1.9540800003574E-07,	-1.91428289041888E-07,	-1.87527193701395E-07,	-1.83703212301207E-07,	-1.79954870761947E-07,	-1.76280722150421E-07,	-1.72679346200318E-07,	-1.69149348840982E-07,	-1.656893617341E-07,	-1.62298041818218E-07,	-1.58974070860924E-07,	-1.55716155018602E-07,	-1.52523024403626E-07,	-1.49393432658868E-07,	-1.4632615653942E-07,	-1.43319995501404E-07,	-1.40373771297749E-07,	-1.37486327580842E-07,	-1.34656529511929E-07,	-1.31883263377159E-07,	-1.2916543621017E-07,	-1.26501975421105E-07,	-1.23891828431956E-07,	-1.21333962318139E-07,	-1.18827363456191E-07,	-1.16371037177498E-07,	-1.13964007427948E-07,	-1.11605316433424E-07,	-1.09294024371026E-07,	-1.07029209045947E-07,	-1.04809965573901E-07,	-1.02635406069008E-07,	-1.00504659337064E-07,	-9.84168705740874E-08,	-9.63712010700699E-08,	-9.43668279178435E-08,	-9.24029437269769E-08,	-9.04787563426215E-08,	-8.85934885692253E-08,	-8.67463778990348E-08,	-8.49366762453055E-08,	-8.31636496801449E-08,	-8.14265781769104E-08,	-7.97247553570875E-08,	-7.80574882415749E-08,	-7.64240970063027E-08,	-7.48239147421128E-08,	-7.32562872188304E-08,	-7.1720572653458E-08,	-7.02161414824223E-08,	-6.87423761378079E-08,	-6.72986708275107E-08,	-6.58844313192459E-08,	-6.44990747283474E-08,	-6.31420293092935E-08,	-6.18127342508985E-08,	-6.05106394751082E-08,	-5.92352054393387E-08,	-5.79859029422999E-08,	-5.67622129332449E-08,	-5.5563626324588E-08,	-5.43896438078343E-08,	-5.32397756727664E-08,	-5.21135416298314E-08,	-5.1010470635677E-08,	-4.99301007217808E-08,	-4.88719788261225E-08,	-4.78356606278477E-08,	-4.68207103848708E-08,	-4.58267007743693E-08,	-4.48532127361199E-08,	-4.38998353186276E-08,	-4.29661655280011E-08,	-4.20518081795278E-08,	-4.11563757519016E-08,	-4.02794882440599E-08,	-3.9420773034584E-08,	-3.85798647436195E-08,	-3.77564050972744E-08,	-3.69500427944513E-08,	-3.61604333760732E-08,	-3.53872390966614E-08,	-3.4630128798225E-08,	-3.38887777864228E-08,	-3.31628677089583E-08,	-3.24520864361697E-08,	-3.17561279437762E-08,	-3.10746921977451E-08,	-3.04074850412416E-08,	-2.97542180836262E-08,	-2.91146085914639E-08,	-2.84883793815109E-08,	-2.78752587156441E-08,	-2.72749801977002E-08,	-2.66872826721906E-08,	-2.61119101248605E-08,	-2.55486115850593E-08,	-2.49971410298913E-08,	-2.44572572901156E-08,	-2.39287239577644E-08,	-2.3411309295451E-08,	-2.29047861473354E-08,	-2.24089318517219E-08,	-2.19235281552575E-08,	-2.1448361128704E-08,	-2.0983221084257E-08,	-2.05279024943829E-08,	-2.00822039121496E-08,	-1.96459278930218E-08,	-1.92188809180982E-08,	-1.88008733187629E-08,	-1.83917192027269E-08,	-1.79912363814356E-08,	-1.75992462988176E-08,	-1.72155739613512E-08,	-1.68400478694256E-08,	-1.64724999499735E-08,	-1.6112765490353E-08,	-1.5760683073456E-08,	-1.54160945140224E-08,	-1.50788447961372E-08,	-1.47487820118912E-08,	-1.4425757301183E-08,	-1.41096247926427E-08,	-1.38002415456579E-08,	-1.3497467493481E-08,	-1.32011653873993E-08,	-1.29112007419491E-08,	-1.26274417811548E-08,	-1.23497593857738E-08,	-1.20780270415315E-08,	-1.18121207883256E-08,	-1.1551919170385E-08,	-1.12973031873643E-08,	-1.10481562463581E-08,	-1.0804364114818E-08,	-1.05658148743566E-08,	-1.03323988754212E-08,	-1.01040086928236E-08,	-9.88053908210867E-09,	-9.66188693674715E-09,	-9.44795124613808E-09,	-9.23863305440579E-09,	-9.03383541997721E-09,	-8.8334633759252E-09,	-8.63742389106421E-09,	-8.44562583178435E-09,	-8.2579799246106E-09,	-8.07439871947383E-09,	-7.89479655368075E-09,	-7.71908951657005E-09,	-7.54719541484199E-09,	-7.37903373854947E-09,	-7.21452562773814E-09,	-7.05359383972391E-09,	-6.89616271699606E-09,	-6.74215815573447E-09,	-6.59150757492976E-09,	-6.44413988609522E-09,	-6.29998546355964E-09,	-6.15897611533035E-09,	-6.021045054516E-09,	-5.88612687129871E-09,	-5.75415750544547E-09,	-5.62507421934887E-09,	-5.49881557158727E-09,	-5.37532139099493E-09,	-5.25453275123255E-09,	-5.13639194584898E-09,	-5.02084246382506E-09,	-4.90782896559046E-09,	-4.79729725950498E-09,	-4.68919427879538E-09,	-4.58346805893955E-09,	-4.4800677154895E-09,	-4.37894342232498E-09,	-4.28004639032986E-09,	-4.18332884648317E-09,	-4.08874401335724E-09,	-3.99624608901512E-09,	-3.90579022729995E-09,	-3.81733251850898E-09,	-3.73082997044482E-09,	-3.64624048983698E-09,	-3.56352286412681E-09,	-3.48263674360881E-09,	-3.40354262392176E-09,	-3.32620182888305E-09,	-3.25057649365966E-09,	-3.17662954826959E-09,	-3.10432470140729E-09,	-3.0336264245872E-09,	-2.96449993659922E-09,	-2.89691118827023E-09,	-2.83082684752593E-09,	-2.76621428474724E-09,	-2.70304155841571E-09,	-2.64127740104238E-09,	-2.5808912053748E-09,	-2.52185301087689E-09,	-2.46413349047635E-09,	-2.4077039375747E-09,	-2.3525362533147E-09,	-2.29860293410055E-09,	-2.24587705936573E-09,	-2.19433227958391E-09,	-2.14394280451833E-09,	-2.09468339170487E-09,	-2.04652933516461E-09,	-1.99945645434121E-09,	-1.95344108325902E-09,	-1.9084600598975E-09,	-1.86449071577788E-09,	-1.82151086575802E-09,	-1.77949879803136E-09,	-1.73843326432602E-09,	-1.69829347030035E-09,	-1.65905906613097E-09,	-1.62071013728955E-09,	-1.58322719550497E-09,	-1.54659116990684E-09,	-1.51078339834726E-09,	-1.47578561889715E-09,	-1.44157996151378E-09,	-1.40814893987624E-09,	-1.37547544338559E-09,	-1.34354272932637E-09,	-1.31233441518652E-09,	-1.28183447113241E-09,	-1.25202721263614E-09,	-1.22289729325207E-09,	-1.19442969753958E-09,	-1.16660973412937E-09,	-1.13942302893036E-09,	-1.1128555184745E-09,	-1.08689344339678E-09,	-1.06152334204775E-09,	-1.03673204423604E-09,	-1.01250666509821E-09,	-9.88834599093532E-10,	-9.65703514121194E-10,	-9.43101345757522E-10,	-9.2101629161087E-10,	-8.99436805791845E-10,	-8.78351593496601E-10,	-8.57749605700968E-10,	-8.37620033963221E-10,	-8.1795230533335E-10,	-7.98736077366717E-10,	-7.79961233240034E-10,	-7.6161787696763E-10,	-7.43696328716037E-10,	-7.26187120214912E-10,	-7.09080990262412E-10,	-6.92368880323122E-10,	-6.76041930216709E-10,	-6.60091473895487E-10,	-6.44509035309127E-10,	-6.29286324354769E-10,	-6.14415232910835E-10,	-5.99887830952876E-10,	-5.85696362749798E-10,	-5.71833243138876E-10,	-5.58291053877976E-10,	-5.45062540073422E-10,	-5.32140606682015E-10,	-5.19518315085702E-10,	-5.07188879737437E-10,	-4.95145664876804E-10,	-4.83382181314011E-10,	-4.71892083280854E-10,	-4.60669165347328E-10,	-4.49707359402549E-10,	-4.39000731698686E-10,	-4.28543479956644E-10,	-4.18329930532236E-10,	-4.08354535641636E-10,	-3.98611870644896E-10,	-3.89096631386374E-10,	-3.79803631590893E-10,	-3.70727800314528E-10,	-3.61864179448888E-10,	-3.53207921277825E-10,	-3.44754286085496E-10,	-3.36498639814733E-10,	-3.28436451774702E-10,	-3.20563292396841E-10,	-3.12874831038097E-10,	-3.05366833830495E-10,	-2.98035161576097E-10,	-2.90875767686419E-10,	-2.83884696165394E-10,	-2.77058079635006E-10,	-2.70392137402695E-10,	-2.63883173569707E-10,	-2.57527575179526E-10,	-2.51321810405571E-10,	-2.45262426777366E-10,	-2.39346049444369E-10,	-2.33569379476703E-10,	-2.27929192202026E-10,	-2.22422335577785E-10,	-2.17045728598145E-10,	-2.11796359734853E-10,	-2.0667128541136E-10,	-2.01667628509496E-10,	-1.96782576908037E-10,	-1.92013382052496E-10,	-1.87357357555498E-10,	-1.828118778271E-10,	-1.7837437673444E-10,	-1.74042346290099E-10,	-1.69813335368594E-10,	-1.65684948450396E-10,	-1.61654844392926E-10,	-1.57720735227943E-10,	-1.53880384984794E-10,	-1.50131608538975E-10,	-1.46472270485479E-10,	-1.42900284036419E-10,	-1.39413609942404E-10,	-1.36010255437194E-10,	-1.32688273205123E-10,	-1.29445760370831E-10,	-1.26280857510836E-10,	-1.2319174768647E-10,	-1.20176655497762E-10,	-1.172338461578E-10,	-1.14361624587156E-10,	-1.11558334527951E-10,	-1.08822357677145E-10,	-1.06152112838642E-10,	-1.03546055093831E-10,	-1.01002674990151E-10,	-9.85204977473175E-11,	-9.60980824808284E-11,	-9.37340214423924E-11,	-9.14269392769118E-11,	-8.91754922956777E-11,	-8.69783677654313E-11,	-8.4834283212955E-11,	-8.27419857448647E-11,	-8.07002513822813E-11,	-7.87078844100639E-11,	-7.67637167402966E-11,	-7.48666072897264E-11,	-7.30154413708534E-11,	-7.12091300963845E-11,	-6.94466097967652E-11,	-6.77268414505104E-11,	-6.6048810127062E-11,	-6.44115244419046E-11,	-6.28140160236788E-11,	-6.12553389930343E-11,	-5.97345694529728E-11,	-5.82508049904339E-11,	-5.68031641888835E-11,	-5.53907861516684E-11,	-5.40128300359071E-11,	-5.26684745966892E-11,	-5.13569177413635E-11,	-5.00773760936959E-11,	-4.88290845676874E-11,	-4.76112959508417E-11,	-4.64232804966803E-11,	-4.52643255263054E-11,	-4.41337350388152E-11,	-4.3030829330381E-11,	-4.19549446217979E-11,	-4.09054326943278E-11,	-3.98816605336531E-11,	-3.88830099817671E-11,	-3.79088773966289E-11,	-3.69586733194136E-11,	-3.60318221491946E-11,	-3.51277618248948E-11,	-3.42459435143509E-11,	-3.33858313103338E-11,	-3.25469019333764E-11,	-3.17286444412582E-11,	-3.09305599450038E-11,	-3.01521613312514E-11,	-2.93929729908539E-11,	-2.86525305535761E-11,	-2.7930380628754E-11,	-2.72260805517878E-11,	-2.65391981363386E-11,	-2.58693114321061E-11,	-2.52160084880642E-11,	-2.45788871210343E-11,	-2.39575546894804E-11,	-2.3351627872411E-11,	-2.27607324532749E-11,	-2.21845031087426E-11,	-2.16225832022641E-11,	-2.10746245823005E-11,	-2.05402873851231E-11,	-2.00192398420833E-11,	-1.95111580912505E-11,	-1.9015725993326E-11,	-1.85326349517348E-11,	-1.80615837368049E-11,	-1.76022783139438E-11,	-1.71544316757233E-11,	-1.67177636777853E-11,	-1.62920008784858E-11,	-1.58768763821928E-11,	-1.54721296861578E-11,	-1.50775065308808E-11,	-1.4692758753892E-11,	-1.43176441468741E-11,	-1.39519263160508E-11,	-1.35953745457684E-11,	-1.32477636652001E-11,	-1.29088739181033E-11,	-1.25784908355613E-11,	-1.22564051116433E-11,	-1.19424124819175E-11,	-1.16363136047531E-11,	-1.13379139453501E-11,	-1.10470236624336E-11,	-1.07634574975555E-11,	-1.0487034666944E-11,	-1.02175787558429E-11,	-9.95491761528674E-12,	-9.69888326125573E-12,	-9.44931177615755E-12,	-9.20604321258361E-12,	-8.96892149928863E-12,	-8.73779434934359E-12,	-8.51251317041294E-12,	-8.2929329771084E-12,	-8.07891230537241E-12,	-7.87031312884548E-12,	-7.66700077717261E-12,	-7.46884385620506E-12,	-7.27571417005454E-12,	-7.08748664495796E-12,	-6.9040392549118E-12,	-6.72525294903587E-12,	-6.55101158062748E-12,	-6.38120183786751E-12,	-6.21571317614101E-12,	-6.05443775193556E-12,	-5.89727035828162E-12,	-5.74410836169985E-12,	-5.59485164062095E-12,	-5.44940252524469E-12,	-5.30766573880523E-12,	-5.16954834021067E-12,	-5.03495966802557E-12,	-4.90381128576571E-12,	-4.77601692847512E-12,	-4.65149245055617E-12,	-4.53015577482395E-12,	-4.41192684275694E-12,	-4.29672756591673E-12,	-4.18448177850975E-12,	-4.07511519106507E-12,	-3.96855534520252E-12,	-3.86473156946616E-12,	-3.76357493619859E-12,	-3.66501821943228E-12,	-3.56899585377435E-12,	-3.47544389426211E-12,	-3.38429997716693E-12,	-3.2955032817246E-12,	-3.20899449277079E-12,	-3.12471576426077E-12,	-3.04261068365297E-12,	-2.96262423713641E-12,	-2.88470277568251E-12,	-2.80879398190215E-12,	-2.73484683768947E-12,	-2.66281159263412E-12,	-2.59263973318403E-12,	-2.52428395254161E-12,	-2.45769812127599E-12,	-2.39283725863499E-12,	-2.32965750454036E-12,	-2.26811609225051E-12,	-2.20817132167516E-12,	-2.14978253332681E-12,	-2.09291008289404E-12,	-2.03751531642233E-12,	-1.98356054608816E-12,	-1.9310090265525E-12,	-1.87982493188024E-12,	-1.8299733330123E-12,	-1.7814201757775E-12,	-1.73413225943161E-12,	-1.68807721571113E-12,	-1.64322348838993E-12,	-1.59954031332676E-12,	-1.55699769899227E-12,	-1.51556640746423E-12,	-1.47521793587998E-12,	-1.43592449833533E-12,	-1.39765900821945E-12,	-1.36039506097552E-12,	-1.32410691727707E-12,	-1.28876948661019E-12,	-1.25435831125221E-12,	-1.22084955063728E-12,	-1.18821996609992E-12,	-1.15644690598742E-12,	-1.12550829113264E-12,	-1.09538260067842E-12,	-1.06604885824541E-12,	-1.03748661843527E-12,	-1.00967595366117E-12,	-9.82597441297826E-13,	-9.56232151143683E-13,	-9.30561633187586E-13,	-9.05567905672917E-13,	-8.81233443452019E-13,	-8.57541166624059E-13,	-8.34474429449567E-13,	-8.12017009535078E-13,	-7.90153097281449E-13,	-7.68867285589558E-13,	-7.48144559817267E-13,	-7.27970287981641E-13,	-7.08330211200576E-13,	-6.89210434368126E-13,	-6.7059741705793E-13,	-6.52477964649302E-13,	-6.34839219670647E-13,	-6.17668653354999E-13,	-6.00954057402619E-13,	-5.8468353594566E-13,	-5.68845497710086E-13,	-5.5342864837009E-13,	-5.38421983090396E-13,	-5.23814779251937E-13,	-5.09596589356492E-13,	-4.95757234105984E-13,	-4.82286795652241E-13,	-4.69175611013094E-13,	-4.56414265650836E-13,	-4.43993587209098E-13,	-4.31904639404339E-13,	-4.20138716068217E-13,	-4.08687335337179E-13,	-3.97542233985744E-13,	-3.8669536189997E-13,	-3.76138876687747E-13,	-3.65865138422569E-13,	-3.55866704517587E-13,	-3.46136324726752E-13,	-3.36666936269999E-13,	-3.27451659079429E-13,	-3.18483791163581E-13,	-3.09756804086903E-13,	-3.01264338561633E-13,	-2.93000200149355E-13,	-2.84958355069562E-13,	-2.77132926112606E-13,	-2.69518188654518E-13,	-2.62108566771191E-13,	-2.54898629449511E-13,	-2.47883086893075E-13,	-2.41056786920176E-13,	-2.34414711451805E-13,	-2.27951973087475E-13,	-2.21663811766702E-13,	-2.15545591514068E-13,	-2.09592797265802E-13,	-2.03801031775883E-13,	-1.98166012599729E-13,	-1.92683569153551E-13,	-1.87349639847521E-13,	-1.82160269290944E-13,	-1.77111605567662E-13,	-1.72199897579963E-13,	-1.67421492459307E-13,	-1.62772833042222E-13,	-1.5825045540978E-13,	-1.53850986489061E-13,	-1.49571141715096E-13,	-1.4540772275179E-13,	-1.41357615270365E-13,	-1.37417786783909E-13,	-1.33585284536634E-13,	-1.29857233446508E-13,	-1.26230834099912E-13,	-1.22703360797065E-13,	-1.19272159646934E-13,	-1.15934646710416E-13,	-1.12688306190592E-13,	-1.0953068866888E-13,	-1.06459409385953E-13,	-1.03472146566306E-13,	-1.00566639785389E-13,	-9.77406883782476E-14,	-9.49921498886327E-14,	-9.23189385575784E-14,	-8.97190238504606E-14,	-8.71904290215788E-14,	-8.47312297153221E-14,	-8.23395526030103E-14,	-8.00135740545153E-14,	-7.77515188437964E-14,	-7.55516588874998E-14,	-7.34123120157961E-14,	-7.13318407746483E-14,	-6.93086512587232E-14,	-6.73411919741785E-14,	-6.54279527305763E-14,	-6.35674635611921E-14,	-6.17582936710066E-14,	-5.99990504116858E-14,	-5.82883782828694E-14,	-5.66249579591077E-14,	-5.50075053418009E-14,	-5.34347706355112E-14,	-5.1905537448035E-14,	-5.04186219136343E-14,	-4.89728718388458E-14,	-4.75671658702963E-14,	-4.62004126839692E-14,	-4.48715501953811E-14,	-4.35795447901393E-14,	-4.23233905743649E-14,	-4.11021086444793E-14,	-3.99147463758633E-14,	-3.87603767299115E-14,	-3.76380975790159E-14,	-3.65470310490225E-14,	-3.54863228787199E-14,	-3.44551417959258E-14,	-3.34526789097511E-14,	-3.24781471186287E-14,	-3.15307805337089E-14,	-3.06098339172273E-14,	-2.97145821354667E-14,	-2.88443196259398E-14,	-2.79983598784308E-14,	-2.71760349295426E-14,	-2.6376694870405E-14,	-2.55997073672078E-14,	-2.48444571942311E-14,	-2.41103457790545E-14,	-2.33967907596323E-14,	-2.27032255529323E-14,	-2.20290989348416E-14,	-2.13738746310517E-14,	-2.07370309186397E-14,	-2.0118060238074E-14,	-1.95164688153752E-14,	-1.89317762941715E-14,	-1.83635153773963E-14,	-1.7811231478378E-14,	-1.72744823810823E-14,	-1.67528379092704E-14,	-1.62458796043447E-14,	-1.57532004116571E-14,	-1.52744043750629E-14,	-1.48091063395074E-14,	-1.43569316614378E-14,	-1.39175159268387E-14,	-1.34905046766945E-14,	-1.30755531396858E-14,	-1.26723259719354E-14,	-1.22804970036175E-14,	-1.1899748992257E-14,	-1.15297733825424E-14,	-1.11702700724847E-14,	-1.08209471857583E-14,	-1.04815208500629E-14,	-1.01517149813504E-14,	-9.8312610737647E-15,	-9.51989799514562E-15,	-9.2173717879528E-15,	-8.9234354754684E-15,	-8.63784887314132E-15,	-8.36037840493906E-15,	-8.0907969245769E-15,	-7.82888354149719E-15,	-7.57442345147507E-15,	-7.32720777172975E-15,	-7.08703338042386E-15,	-6.85370276043634E-15,	-6.62702384729706E-15,	-6.40680988117452E-15,	-6.19287926281049E-15,	-5.98505541329846E-15,	-5.78316663760511E-15,	-5.58704599173684E-15,	-5.3965311534557E-15,	-5.21146429645185E-15,	-5.0316919678816E-15,	-4.85706496918288E-15,	-4.6874382400819E-15,	-4.52267074570727E-15,	-4.36262536672976E-15,	-4.20716879244812E-15,	-4.05617141674352E-15,	-3.90950723682683E-15,	-3.76705375470545E-15,	-3.62869188129775E-15,	-3.49430584312531E-15,	-3.36378309151511E-15,	-3.23701421424512E-15,	-3.11389284956902E-15,	-2.99431560255686E-15,	-2.87818196369069E-15,	-2.7653942296553E-15,	-2.65585742626603E-15,	-2.54947923347704E-15,	-2.44616991241485E-15,	-2.34584223438351E-15,	-2.2484114117891E-15,	-2.15379503093245E-15,	-2.06191298662074E-15,	-1.97268741854935E-15,	-1.88604264940716E-15,	-1.80190512465922E-15,	-1.72020335396234E-15,	-1.64086785416998E-15,	-1.5638310938842E-15,	-1.48902743951335E-15,	-1.41639310279539E-15,	-1.34586608974768E-15,	-1.2773861510053E-15,	-1.21089473351053E-15,	-1.14633493351772E-15,	-1.08365145087813E-15,	-1.02279054457053E-15,	-9.63699989444377E-16,	-9.06329034142866E-16,	-8.50628360174419E-16,	-7.96550042101724E-16,	-7.44047508818393E-16,	-6.93075505884036E-16,	-6.43590058889325E-16,	-5.95548437823406E-16,	-5.4890912241668E-16,	-5.03631768432751E-16,	-4.59677174883994E-16,	-4.17007252145876E-16,	-3.75584990945834E-16,	-3.35374432203123E-16,	-2.96340637696714E-16,	-2.58449661538896E-16,	-2.21668522432823E-16,	-1.85965176692864E-16,	-1.51308492007126E-16,	-1.17668221922095E-16,	-8.50149810298787E-17,	-5.33202208390093E-17,	-2.25562063103369E-17,	7.30400696002473E-18,	3.62865949282281E-17,	6.44169866777655E-17,	9.17198853986977E-17,	1.18219288787218E-16,	1.4393850888128E-16,	1.68900191347628E-16,	1.93126334235121E-16,	2.16638306208877E-16,	2.39456864279374E-16,	2.61602171040235E-16,	2.83093811428094E-16,	3.03950809017563E-16,	3.24191641863982E-16,	3.4383425790627E-16,	3.62896089941911E-16,	3.81394070185727E-16,	3.99344644423839E-16,	4.16763785773856E-16,	4.33667008062081E-16,	4.50069378828196E-16,	4.65985531967641E-16,	4.81429680021595E-16,	4.96415626124222E-16,	5.10956775616589E-16,	5.25066147336371E-16,	5.38756384592283E-16,	5.52039765831851E-16,	5.64928215010979E-16,	5.77433311673496E-16,	5.8956630074866E-16,	6.01338102074374E-16,	6.12759319653672E-16,	6.23840250651824E-16,	6.34590894141189E-16,	6.45020959600791E-16,	6.55139875177372E-16,	6.64956795714509E-16,	6.74480610556189E-16,	6.83719951131094E-16,	6.92683198323623E-16,	7.01378489637587E-16,	7.09813726158277E-16,	7.17996579318512E-16,	7.25934497474078E-16,	7.33634712293853E-16,	7.41104244969735E-16,	7.48349912251402E-16,	7.55378332310724E-16,	7.62195930440598E-16,	7.68808944592764E-16,	7.75223430759118E-16,	7.81445268200834E-16,	7.87480164529553E-16,	7.93333660644753E-16,	7.99011135531296E-16,	8.04517810921049E-16,	8.09858755822377E-16,	8.15038890921177E-16,	8.20062992857055E-16,	8.24935698378107E-16,	8.29661508377722E-16,	8.34244791816683E-16,	8.38689789533774E-16,	8.43000617948023E-16,	8.47181272655603E-16,	8.5123563192434E-16,	8.55167460088703E-16,	8.58980410848054E-16,	8.62678030470885E-16,	8.66263760907657E-16,	8.69740942814834E-16,	8.73112818492577E-16,	8.76382534738546E-16,	8.79553145620155E-16,	8.82627615167572E-16,	8.85608819989702E-16,	8.8849955181532E-16,	8.91302519961448E-16,	8.94020353731045E-16,	8.96655604741986E-16,	8.99210749189278E-16,	9.01688190042394E-16,	9.04090259179557E-16,	9.06419219460742E-16,	9.08677266741156E-16,	9.10866531826847E-16,	9.12989082374102E-16,	9.15046924734206E-16,	9.1704200574513E-16,	9.18976214471629E-16,	9.20851383895225E-16,	9.22669292555499E-16,	9.24431666144055E-16,	9.26140179052522E-16,	9.27796455875886E-16,	9.29402072872411E-16,	9.30958559381413E-16,	9.32467399200053E-16,	9.33930031920332E-16,	9.35347854227419E-16,	9.36722221160409E-16,	9.38054447336583E-16,	9.39345808140212E-16,	9.40597540876911E-16,	9.41810845894536E-16,	9.42986887671562E-16,	9.4412679587389E-16,	9.4523166638097E-16,	9.46302562282124E-16,	9.4734051484392E-16,	9.48346524449423E-16,	9.49321561510136E-16,	9.50266567351396E-16,	9.51182455072008E-16,	9.52070110378834E-16,	9.52930392397073E-16,	9.53764134456918E-16,	9.54572144857269E-16,	9.55355207607167E-16,	9.56114083145589E-16,	9.5684950904021E-16,	9.57562200665767E-16,	9.58252851862578E-16,	9.58922135575817E-16,	9.59570704476075E-16,	9.60199191561764E-16,	9.60808210743883E-16,	9.61398357413639E-16,	9.61970208993447E-16,	9.62524325471764E-16,	9.63061249922232E-16,	9.63581509007588E-16,	9.64085613468778E-16,	9.64574058599694E-16,	9.65047324707971E-16,	9.65505877562224E-16,	9.65950168826135E-16,	9.66380636479762E-16,	9.66797705228445E-16,	9.67201786899662E-16,	9.67593280828199E-16,	9.67972574229954E-16,	9.68340042564719E-16,	9.68696049888257E-16,	9.69040949193986E-16,	9.69375082744563E-16,	9.69698782393678E-16,	9.70012369898332E-16,	9.70316157221877E-16,	9.7061044682809E-16,	9.7089553196654E-16,	9.71171696949508E-16,	9.71439217420697E-16,	9.71698360615978E-16,	9.719493856164E-16,	9.72192543593703E-16,	9.72428078048529E-16,	9.72656225041571E-16,	9.72877213417851E-16,	9.73091265024324E-16,	9.73298594921024E-16,	9.73499411585919E-16,	9.73693917113673E-16,	9.73882307408487E-16,	9.74064772371199E-16,	9.74241496080802E-16,	9.74412656970555E-16,	9.74578427998839E-16,	9.74738976814913E-16,	9.7489446591972E-16,	9.75045052821898E-16,	9.75190890189113E-16,	9.75332125994884E-16,	9.75468903661009E-16,	9.75601362195725E-16,	9.75729636327744E-16,	9.75853856636259E-16,	9.7597414967707E-16,	9.7609063810492E-16,	9.76203440792163E-16,	9.76312672943872E-16,	9.76418446209498E-16,	9.76520868791166E-16,	9.76620045548724E-16,	9.76716078101631E-16,	9.76809064927789E-16,	9.76899101459388E-16,	9.7698628017588E-16,	9.77070690694148E-16,	9.77152419855956E-16,	9.77231551812769E-16,	9.77308168108012E-16,	9.77382347756853E-16,	9.77454167323568E-16,	9.77523700996583E-16,	9.77591020661242E-16,	9.77656195970377E-16,	9.77719294412748E-16,	9.77780381379409E-16,	9.77839520228068E-16,	9.77896772345501E-16,	9.77952197208073E-16,	9.78005852440426E-16,	9.78057793872389E-16,	9.78108075594163E-16,	9.7815675000983E-16,	9.78203867889241E-16,	9.78249478418324E-16,	9.7829362924787E-16,	9.78336366540831E-16,	9.78377735018183E-16,	9.78417778003394E-16,	9.78456537465537E-16,	9.78494054061095E-16,	9.78530367174485E-16,	9.78565514957356E-16,	9.78599534366686E-16,	9.78632461201713E-16,	9.78664330139744E-16,	9.78695174770873E-16,	9.78725027631633E-16,	9.78753920237624E-16,	9.78781883115143E-16,	9.78808945831845E-16,	9.78835137026472E-16,	9.78860484437663E-16,	9.78885014931893E-16,	9.78908754530548E-16,	9.78931728436171E-16,	9.78953961057909E-16,	9.78975476036174E-16,	9.78996296266551E-16,	9.79016443922972E-16,	9.79035940480172E-16,	9.79054806735468E-16,	9.79073062829854E-16,	9.79090728268459E-16,	9.79107821940372E-16,	9.7912436213785E-16,	9.79140366574949E-16,	9.79155852405567E-16,	9.7917083624094E-16,	9.79185334166594E-16,	9.79199361758781E-16,	9.79212934100394E-16,	9.79226065796409E-16,	9.79238770988835E-16,	9.79251063371212E-16,	9.79262956202653E-16,	9.7927446232146E-16,	9.79285594158311E-16,	9.79296363749048E-16,	9.79306782747057E-16,	9.79316862435277E-16,	9.7932661373783E-16,	9.7933604723129E-16,	9.79345173155606E-16,	9.79354001424685E-16,	9.79362541636645E-16,	9.79370803083751E-16,	9.79378794762042E-16,	9.79386525380659E-16,	9.79394003370885E-16,	9.79401236894898E-16,	9.79408233854258E-16,	9.79415001898127E-16,	9.79421548431231E-16,	9.79427880621581E-16,	9.79434005407942E-16,	9.79439929507082E-16,	9.79445659420785E-16,	9.79451201442648E-16,	9.79456561664669E-16,	9.79461745983625E-16,	9.79466760107253E-16,	9.79471609560237E-16,	9.79476299690007E-16,	9.79480835672359E-16,	9.79485222516896E-16,	9.79489465072303E-16,	9.79493568031449E-16,	9.7949753593634E-16,	9.79501373182905E-16,	9.79505084025643E-16,	9.79508672582113E-16,	9.79512142837295E-16,	9.79515498647804E-16,	9.79518743745976E-16,	9.79521881743826E-16,	9.79524916136878E-16,	9.79527850307884E-16,	9.79530687530408E-16,	9.79533430972319E-16,	9.79536083699159E-16,	9.79538648677404E-16,	9.79541128777638E-16,	9.79543526777606E-16,	9.79545845365188E-16,	9.79548087141266E-16,	9.79550254622512E-16,	9.79552350244079E-16,	9.79554376362212E-16,	9.79556335256772E-16,	9.79558229133688E-16,	9.7956006012732E-16,	9.7956183030276E-16,	9.79563541658049E-16,	9.79565196126328E-16,	9.79566795577928E-16,	9.79568341822379E-16,	9.79569836610366E-16,	9.79571281635623E-16,	9.79572678536761E-16,	9.79574028899042E-16,	9.79575334256094E-16,	9.79576596091576E-16,	9.79577815840785E-16,	9.79578994892216E-16,	9.79580134589069E-16,	9.79581236230707E-16,	9.79582301074076E-16,	9.79583330335069E-16,	9.79584325189852E-16,	9.79585286776148E-16,	9.79586216194478E-16,	9.79587114509365E-16,	9.79587982750495E-16,	9.79588821913848E-16,	9.79589632962784E-16,	9.79590416829103E-16,	9.79591174414063E-16,	9.79591906589375E-16,	9.79592614198154E-16,	9.79593298055853E-16,	9.79593958951155E-16,	9.79594597646844E-16,	9.79595214880644E-16,	9.79595811366037E-16,	9.79596387793046E-16,	9.79596944828999E-16,	9.79597483119267E-16,	9.7959800328798E-16,	9.79598505938714E-16,	9.79598991655165E-16,	9.79599461001792E-16,	9.79599914524447E-16,	9.7960035275098E-16,	9.79600776191826E-16,	9.79601185340572E-16,	9.7960158067451E-16,	9.79601962655166E-16,	9.79602331728816E-16,	9.79602688326981E-16,	9.79603032866917E-16,	9.79603365752074E-16,	9.7960368737255E-16,	9.79603998105528E-16,	9.79604298315699E-16,	9.79604588355672E-16,	9.79604868566367E-16,	9.796051392774E-16,	9.79605400807455E-16,	9.79605653464639E-16,	9.79605897546831E-16,	9.7960613334202E-16,	9.79606361128625E-16,	9.79606581175812E-16,	9.79606793743799E-16,	9.7960699908415E-16,	9.79607197440057E-16,	9.79607389046617E-16,	9.79607574131099E-16,	9.79607752913203E-16,	9.79607925605306E-16,	9.79608092412705E-16,	9.79608253533852E-16,	9.79608409160577E-16,	9.79608559478307E-16,	9.7960870466628E-16,	9.79608844897746E-16,	9.79608980340166E-16,	9.79609111155404E-16,	9.79609237499912E-16,	9.79609359524908E-16,	9.79609477376548E-16,	9.79609591196099E-16,	9.79609701120093E-16,	9.79609807280489E-16,	9.79609909804825E-16,	9.79610008816362E-16,	9.79610104434225E-16,	9.79610196773545E-16,	9.79610285945586E-16,	9.79610372057878E-16,	9.79610455214337E-16,	9.79610535515387E-16,	9.79610613058078E-16,	9.79610687936195E-16,	9.79610760240366E-16,	9.79610830058173E-16,	9.79610897474245E-16,	9.79610962570364E-16,	9.79611025425555E-16,	9.79611086116179E-16,	9.79611144716022E-16,	9.79611201296382E-16,	9.79611255926152E-16,	9.79611308671897E-16,	9.79611359597935E-16,	9.79611408766412E-16,	9.79611456237373E-16,	9.79611502068833E-16,	9.79611546316845E-16,	9.79611589035566E-16,	9.79611630277319E-16,	9.79611670092657E-16,	9.79611708530418E-16,	9.79611745637787E-16,	9.7961178146035E-16,	9.79611816042146E-16,	9.79611849425721E-16,	9.79611881652177E-16,	9.79611912761219E-16,	9.79611942791207E-16,	9.79611971779194E-16,	9.79611999760976E-16,	9.79612026771131E-16,	9.79612052843063E-16,	9.79612078009036E-16,	9.79612102300219E-16,	9.79612125746717E-16,	9.79612148377612E-16,	9.79612170220993E-16,	9.79612191303993E-16,	9.79612211652818E-16,	9.79612231292783E-16,	9.79612250248337E-16,	9.79612268543093E-16,	9.7961228619986E-16,	9.79612303240667E-16,	9.79612319686791E-16,	9.79612335558778E-16,	9.79612350876476E-16,	9.79612365659051E-16,	9.79612379925014E-16,	9.79612393692242E-16,	9.79612406977998E-16,	9.79612419798957E-16,	9.7961243217122E-16,	9.79612444110335E-16,	9.79612455631319E-16,	9.79612466748672E-16,	9.79612477476397E-16,	9.79612487828016E-16,	9.79612497816587E-16,	9.79612507454717E-16,	9.79612516754582E-16,	9.79612525727939E-16,	9.79612534386138E-16,	9.79612542740142E-16,	9.79612550800534E-16,	9.79612558577533E-16,	9.79612566081006E-16,	9.79612573320481E-16,	9.79612580305156E-16,	9.79612587043911E-16,	9.79612593545322E-16,	9.79612599817667E-16,	9.79612605868939E-16,	9.79612611706856E-16,	9.79612617338865E-16,	9.79612622772161E-16,	9.79612628013687E-16,	9.79612633070146E-16,	9.7961263794801E-16,	9.79612642653526E-16,	9.79612647192725E-16,	9.7961265157143E-16,	9.7961265579526E-16,	9.79612659869642E-16,	9.79612663799812E-16,	9.79612667590826E-16,	9.79612671247564E-16,	9.79612674774736E-16,	9.79612678176889E-16,	9.79612681458413E-16,	9.79612684623544E-16,	9.79612687676371E-16,	9.7961269062084E-16,	9.79612693460762E-16,	9.79612696199813E-16,	9.79612698841541E-16,	9.79612701389373E-16,	9.79612703846613E-16,	9.79612706216453E-16,	9.79612708501971E-16,	9.7961271070614E-16,	9.79612712831828E-16,	9.79612714881803E-16,	9.79612716858735E-16,	9.79612718765204E-16,	9.79612720603697E-16,	9.79612722376614E-16,	9.79612724086273E-16,	9.7961272573491E-16,	9.79612727324681E-16,	9.79612728857671E-16,	9.79612730335886E-16,	9.79612731761266E-16,	9.79612733135683E-16,	9.7961273446094E-16,	9.7961273573878E-16,	9.79612736970883E-16,	9.79612738158872E-16,	9.79612739304311E-16,	9.79612740408709E-16,	9.79612741473525E-16,	9.79612742500162E-16,	9.79612743489977E-16,	9.79612744444279E-16,	9.79612745364329E-16,	9.79612746251346E-16,	9.79612747106505E-16,	9.7961274793094E-16,	9.79612748725744E-16,	9.79612749491972E-16,	9.79612750230641E-16,	9.79612750942735E-16,	9.796127516292E-16,	9.79612752290951E-16,	9.79612752928868E-16,	9.79612753543803E-16,	9.79612754136575E-16,	9.79612754707976E-16,	9.79612755258769E-16,	9.79612755789692E-16,	9.79612756301453E-16,	9.79612756794739E-16,	9.7961275727021E-16,	9.79612757728504E-16,	9.79612758170235E-16,	9.79612758595997E-16,	9.79612759006361E-16,	9.7961275940188E-16,	9.79612759783086E-16,	9.79612760150491E-16,	9.79612760504591E-16,	9.79612760845864E-16,	9.79612761174769E-16,	9.79612761491751E-16,	9.79612761797238E-16,	9.79612762091643E-16,	9.79612762375365E-16,	9.79612762648788E-16,	9.79612762912282E-16,	9.79612763166204E-16,	9.796127634109E-16,	9.79612763646701E-16,	9.79612763873927E-16,	9.79612764092888E-16,	9.79612764303882E-16,	9.79612764507195E-16,	9.79612764703106E-16,	9.7961276489188E-16,	9.79612765073776E-16,	9.79612765249043E-16,	9.79612765417918E-16,	9.79612765580634E-16,	9.79612765737413E-16,	9.7961276588847E-16,	9.79612766034012E-16,	9.79612766174238E-16,	9.7961276630934E-16,	9.79612766439505E-16,	9.7961276656491E-16,	9.7961276668573E-16,	9.79612766802129E-16,	9.79612766914268E-16,	9.79612767022302E-16,	9.79612767126379E-16,	9.79612767226643E-16,	9.79612767323233E-16,	9.79612767416282E-16,	9.79612767505918E-16,	9.79612767592266E-16,	9.79612767675446E-16,	9.79612767755572E-16,	9.79612767832755E-16,	9.79612767907104E-16,	9.7961276797872E-16,	9.79612768047704E-16,	9.79612768114151E-16,	9.79612768178154E-16,	9.79612768239802E-16,	9.7961276829918E-16,	9.79612768356373E-16,	9.79612768411458E-16,	9.79612768464514E-16,	9.79612768515615E-16,	9.79612768564831E-16,	9.79612768612232E-16,	9.79612768657884E-16,	9.79612768701852E-16,	9.79612768744196E-16,	9.79612768784976E-16,	9.7961276882425E-16,	9.79612768862072E-16,	9.79612768898496E-16,	9.79612768933573E-16,	9.79612768967352E-16,	9.79612768999881E-16,	9.79612769031206E-16,	9.79612769061371E-16,	9.79612769090418E-16,	9.79612769118389E-16,	9.79612769145324E-16,	9.79612769171259E-16,	9.79612769196232E-16,	9.79612769220279E-16,	9.79612769243434E-16,	9.79612769265728E-16,	9.79612769287194E-16,	9.79612769307863E-16,	9.79612769327764E-16,	9.79612769346925E-16,	9.79612769365373E-16,	9.79612769383135E-16,	9.79612769400235E-16,	9.79612769416699E-16,	9.7961276943255E-16,	9.7961276944781E-16,	9.79612769462502E-16,	9.79612769476646E-16,	9.79612769490263E-16,	9.79612769503371E-16,	9.79612769515991E-16,	9.79612769528139E-16,	9.79612769539834E-16,	9.79612769551092E-16,	9.79612769561929E-16,	9.79612769572361E-16,	9.79612769582403E-16,	9.7961276959207E-16,	9.79612769601374E-16,	9.79612769610331E-16,	9.79612769618952E-16,	9.79612769627251E-16,	9.79612769635238E-16,	9.79612769642926E-16,	9.79612769650326E-16,	9.79612769657449E-16,	9.79612769664304E-16,	9.79612769670902E-16,	9.79612769677252E-16,	9.79612769683364E-16,	9.79612769689246E-16,	9.79612769694907E-16,	9.79612769700356E-16,	9.79612769705599E-16,	9.79612769710645E-16,	9.79612769715502E-16,	9.79612769720176E-16,	9.79612769724673E-16,	9.79612769729001E-16,	9.79612769733167E-16,	9.79612769737175E-16,	9.79612769741032E-16,	9.79612769744744E-16,	9.79612769748315E-16,	9.79612769751752E-16,	9.79612769755059E-16,	9.79612769758241E-16,	9.79612769761304E-16,	9.7961276976425E-16,	9.79612769767085E-16,	9.79612769769813E-16,	9.79612769772437E-16,	9.79612769774963E-16,	9.79612769777392E-16,	9.7961276977973E-16,	9.79612769781979E-16,	9.79612769784143E-16,	9.79612769786225E-16,	9.79612769788228E-16,	9.79612769790155E-16,	9.79612769792009E-16,	9.79612769793793E-16,	9.79612769795509E-16,	9.7961276979716E-16,	9.79612769798748E-16,	9.79612769800276E-16,	9.79612769801745E-16,	9.79612769803159E-16,	9.79612769804519E-16,	9.79612769805827E-16,	9.79612769807086E-16,	9.79612769808296E-16,	9.79612769809461E-16,	9.79612769810581E-16,	9.79612769811658E-16,	9.79612769812695E-16,	9.79612769813692E-16,	9.79612769814651E-16,	9.79612769815573E-16,	9.7961276981646E-16,	9.79612769817313E-16,	9.79612769818133E-16,	9.79612769818922E-16,	9.79612769819681E-16,	9.79612769820412E-16,	9.79612769821114E-16,	9.79612769821789E-16,	9.79612769822438E-16,	9.79612769823063E-16,	9.79612769823663E-16,	9.79612769824241E-16,	9.79612769824797E-16,	9.79612769825331E-16,	9.79612769825844E-16,	9.79612769826338E-16,	9.79612769826814E-16,	9.79612769827271E-16,	9.7961276982771E-16,	9.79612769828132E-16,	9.79612769828538E-16,	9.79612769828929E-16,	9.79612769829305E-16,	9.79612769829666E-16,	9.79612769830013E-16,	9.79612769830347E-16,	9.79612769830668E-16,	9.79612769830977E-16,	9.79612769831273E-16,	9.79612769831559E-16,	9.79612769831833E-16,	9.79612769832097E-16,	9.7961276983235E-16,	9.79612769832594E-16,	9.79612769832828E-16,	9.79612769833054E-16,	9.79612769833271E-16,	9.79612769833479E-16,	9.79612769833679E-16,	9.79612769833872E-16,	9.79612769834057E-16,	9.79612769834235E-16,	9.79612769834406E-16,	9.7961276983457E-16,	9.79612769834728E-16,	9.7961276983488E-16,	9.79612769835026E-16,	9.79612769835166E-16,	9.79612769835301E-16,	9.79612769835431E-16,	9.79612769835556E-16,	9.79612769835675E-16,	9.79612769835791E-16,	9.79612769835901E-16,	9.79612769836008E-16,	9.7961276983611E-16,	9.79612769836208E-16,	9.79612769836303E-16,	9.79612769836394E-16,	9.79612769836481E-16,	9.79612769836565E-16,	9.79612769836646E-16,	9.79612769836723E-16,	9.79612769836797E-16,	9.79612769836869E-16,	9.79612769836938E-16,	9.79612769837004E-16,	9.79612769837067E-16,	9.79612769837128E-16,	9.79612769837187E-16,	9.79612769837243E-16,	9.79612769837297E-16,	9.7961276983735E-16,	9.79612769837399E-16,	9.79612769837448E-16,	9.79612769837494E-16,	9.79612769837538E-16,	9.79612769837581E-16,	9.79612769837622E-16,	9.79612769837661E-16,	9.79612769837699E-16,	9.79612769837735E-16,	9.7961276983777E-16,	9.79612769837804E-16,	9.79612769837836E-16,	9.79612769837867E-16,	9.79612769837897E-16,	9.79612769837925E-16,	9.79612769837953E-16,	9.79612769837979E-16,	9.79612769838004E-16,	9.79612769838029E-16,	9.79612769838052E-16,	9.79612769838075E-16,	9.79612769838096E-16,	9.79612769838117E-16,	9.79612769838137E-16,	9.79612769838156E-16,	9.79612769838174E-16,	9.79612769838192E-16,	9.79612769838209E-16,	9.79612769838225E-16,	9.79612769838241E-16,	9.79612769838256E-16,	9.79612769838271E-16,	9.79612769838285E-16,	9.79612769838298E-16,	9.79612769838311E-16,	9.79612769838323E-16,	9.79612769838335E-16,	9.79612769838346E-16,	9.79612769838357E-16,	9.79612769838368E-16,	9.79612769838378E-16,	9.79612769838387E-16,	9.79612769838397E-16,	9.79612769838405E-16,	9.79612769838414E-16,	9.79612769838422E-16,	9.7961276983843E-16,	9.79612769838438E-16,	9.79612769838445E-16,	9.79612769838452E-16,	9.79612769838459E-16,	9.79612769838465E-16,	9.79612769838471E-16,	9.79612769838477E-16,	9.79612769838483E-16,	9.79612769838489E-16,	9.79612769838494E-16,	9.79612769838499E-16,	9.79612769838504E-16,	9.79612769838508E-16,	9.79612769838513E-16,	9.79612769838517E-16,	9.79612769838521E-16,	9.79612769838525E-16,	9.79612769838529E-16,	9.79612769838533E-16,	9.79612769838536E-16,	9.7961276983854E-16,	9.79612769838543E-16,	9.79612769838546E-16,	9.79612769838549E-16,	9.79612769838552E-16,	9.79612769838555E-16,	9.79612769838557E-16,	9.7961276983856E-16,	9.79612769838562E-16,	9.79612769838564E-16,	9.79612769838567E-16,	9.79612769838569E-16,	9.79612769838571E-16,	9.79612769838573E-16,	9.79612769838575E-16,	9.79612769838577E-16,	9.79612769838578E-16,	9.7961276983858E-16,	9.79612769838582E-16,	9.79612769838583E-16,	9.79612769838585E-16,	9.79612769838586E-16,	9.79612769838588E-16,	9.79612769838589E-16,	9.7961276983859E-16,	9.79612769838591E-16,	9.79612769838592E-16,	9.79612769838594E-16,	9.79612769838595E-16,	9.79612769838596E-16,	9.79612769838597E-16,	9.79612769838598E-16,	9.79612769838599E-16,	9.79612769838599E-16,	9.796127698386E-16,	9.79612769838601E-16,	9.79612769838602E-16,	9.79612769838603E-16,	9.79612769838603E-16,	9.79612769838604E-16,	9.79612769838605E-16,	9.79612769838605E-16,	9.79612769838606E-16,	9.79612769838606E-16,	9.79612769838607E-16,	9.79612769838608E-16,	9.79612769838608E-16,	9.79612769838608E-16,	9.79612769838609E-16,	9.79612769838609E-16,	9.7961276983861E-16,	9.7961276983861E-16,	9.79612769838611E-16,	9.79612769838611E-16,	9.79612769838611E-16,	9.79612769838612E-16,	9.79612769838612E-16,	9.79612769838612E-16,	9.79612769838613E-16,	9.79612769838613E-16,	9.79612769838613E-16,	9.79612769838614E-16,	9.79612769838614E-16,	9.79612769838614E-16,	9.79612769838614E-16,	9.79612769838615E-16,	9.79612769838615E-16,	9.79612769838615E-16,	9.79612769838615E-16,	9.79612769838615E-16,	9.79612769838616E-16,	9.79612769838616E-16,	9.79612769838616E-16,	9.79612769838616E-16,	9.79612769838616E-16,	9.79612769838616E-16,	9.79612769838616E-16,	9.79612769838617E-16,	9.79612769838617E-16,	9.79612769838617E-16,	9.79612769838617E-16,	9.79612769838617E-16,	9.79612769838617E-16,	9.79612769838618E-16,	9.79612769838618E-16,	9.79612769838618E-16,	9.79612769838618E-16,	9.79612769838618E-16,	9.79612769838618E-16,	9.79612769838618E-16,	9.79612769838618E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16,	9.79612769838619E-16};
    
    
	double force_1 = Force_1[timestep];
	double force_2 = -1*Force_1[timestep+5];

	double force_coefficient =    2211.84;

	double force_1x = force_1 * force_coefficient;
	double force_1z = force_1 * force_coefficient/2;

	double force_2x = force_2 * force_coefficient;
	double force_2z = force_2 * force_coefficient/2;



	int32_t nindex;
	int32_t k1=0,k2=0;

	double f_l_depth = 512;   //first layer depth
	double el_size   = 2;     //element size
	double s_l_depth = f_l_depth - el_size;    //second layer depth
	double d_width_x   = 8192;    //domain width x
	double d_width_y   = 8192;     //domain width y


	for ( nindex = 0; nindex < myMesh->nharbored; nindex++ ) {

		double z_m = (myMesh->ticksize)*(double)myMesh->nodeTable[nindex].z;
		double x_m = (myMesh->ticksize)*(double)myMesh->nodeTable[nindex].x;
		double y_m = (myMesh->ticksize)*(double)myMesh->nodeTable[nindex].y;

		if ( z_m == f_l_depth ){

			if (((x_m == 0 && y_m == 0) || (x_m == 0 && y_m == d_width_y) || (x_m == d_width_x && y_m == 0) || (x_m == d_width_x && y_m == d_width_y))) {
				fvector_t *nodalForce;
				nodalForce = mySolver->force + nindex;
				nodalForce->f[0] += force_1x/4;


			} else if (((x_m == 0 && (y_m != 0 && y_m != d_width_y )) || (x_m == d_width_x && (y_m != 0 && y_m != d_width_y )) || (y_m == 0 && (x_m != 0 && x_m != d_width_x )) || (y_m == d_width_y && (x_m != 0 && x_m != d_width_x )))) {

				fvector_t *nodalForce;
				nodalForce = mySolver->force + nindex;
				nodalForce->f[0] += force_1x/2;
				k1=k1+1;

			}else {
				fvector_t *nodalForce;
				nodalForce = mySolver->force + nindex;
				nodalForce->f[0] += force_1x;
				k2=k2+1;

			}
		}

		if ( z_m == s_l_depth ) {

			if (((x_m == 0 && y_m == 0) || (x_m == 0 && y_m == d_width_y) || (x_m == d_width_x && y_m == 0) || (x_m == d_width_x && y_m == d_width_y))) {
				fvector_t *nodalForce;
				nodalForce = mySolver->force + nindex;
				nodalForce->f[0] += force_2x/4;


			} else if (((x_m == 0 && (y_m != 0 && y_m != d_width_y )) || (x_m == d_width_x && (y_m != 0 && y_m != d_width_y )) || (y_m == 0 && (x_m != 0 && x_m != d_width_x )) || (y_m == d_width_y && (x_m != 0 && x_m != d_width_x )))) {

				fvector_t *nodalForce;
				nodalForce = mySolver->force + nindex;
				nodalForce->f[0] += force_2x/2;
				k1=k1+1;

			} else {
				fvector_t *nodalForce;
				nodalForce = mySolver->force + nindex;
				nodalForce->f[0] += force_2x;
				k2=k2+1;

			}
		}


		// Add vertical force to keep the balance

		if ( z_m == f_l_depth && x_m == 0 && y_m!=0 && y_m != d_width_y){

			fvector_t *nodalForce;
			nodalForce = mySolver->force + nindex;
			nodalForce->f[2] -= force_1z/2;
		}

		if ( z_m == s_l_depth && x_m == 0 && y_m!=0 && y_m != d_width_y ){

			fvector_t *nodalForce;
			nodalForce = mySolver->force + nindex;
			nodalForce->f[2] -= force_2z/2;
		}

		if ( z_m == f_l_depth && x_m == d_width_x && y_m!=0 && y_m != d_width_y ){

			fvector_t *nodalForce;
			nodalForce = mySolver->force + nindex;
			nodalForce->f[2] += force_1z/2;
		}

		if ( z_m == s_l_depth && x_m == d_width_x && y_m!=0 && y_m != d_width_y ){

			fvector_t *nodalForce;
			nodalForce = mySolver->force + nindex;
			nodalForce->f[2] += force_2z/2;
		}

		// Add vertical force at the corners
		if ( z_m == f_l_depth && x_m == 0 && (y_m==0 || y_m == d_width_y)){

					fvector_t *nodalForce;
					nodalForce = mySolver->force + nindex;
					nodalForce->f[2] -= force_1z/4;
				}

				if ( z_m == s_l_depth && x_m == 0 && (y_m==0 || y_m == d_width_y)){

					fvector_t *nodalForce;
					nodalForce = mySolver->force + nindex;
					nodalForce->f[2] -= force_2z/4;
				}

				if ( z_m == f_l_depth && x_m == d_width_x && (y_m==0 || y_m == d_width_y)){

					fvector_t *nodalForce;
					nodalForce = mySolver->force + nindex;
					nodalForce->f[2] += force_1z/4;
				}

				if ( z_m == s_l_depth && x_m == d_width_x && (y_m==0 || y_m == d_width_y)){

					fvector_t *nodalForce;
					nodalForce = mySolver->force + nindex;
					nodalForce->f[2] += force_2z/4;
				}



	}

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

