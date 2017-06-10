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
    double Force_1[1000]={1.35111255205411E-24,	3.31034338608617E-24,	6.14774662326503E-24,	1.02516565323711E-23,	1.61797605160242E-23,	2.47319038914087E-23,	3.70537343261079E-23,	5.47840515413847E-23,	8.02640110544514E-23,	1.16833748165512E-22,	1.69252387779081E-22,	2.44291965431406E-22,	3.51576141017056E-22,	5.04763006932892E-22,	7.23210915764468E-22,	1.03432141811737E-21,	1.47683009736597E-21,	2.10542256860622E-21,	2.99720053284144E-21,	4.26072715409923E-21,	6.0486635768945E-21,	8.5753984748259E-21,	1.21416102767648E-20,	1.71684420215414E-20,	2.42449864966829E-20,	3.41941763880942E-20,	4.81640893060695E-20,	6.77542991500259E-20,	9.5190488838295E-20,	1.33565430987247E-19,	1.87171110598859E-19,	2.61955839032891E-19,	3.6615254337669E-19,	5.11141156651131E-19,	7.12630723826186E-19,	9.92277545799426E-19,	1.37989718196362E-18,	1.91648413371098E-18,	2.65832849646046E-18,	3.682620415017E-18,	5.09506962403629E-18,	7.04025043229777E-18,	9.7156277975747E-18,	1.33905483033423E-17,	1.8431918256102E-17,	2.53388759593255E-17,	3.47895442634112E-17,	4.77039857230719E-17,	6.53288590160379E-17,	8.93511006442267E-17,	1.22050373325902E-16,	1.66503219226438E-16,	2.26856077338103E-16,	3.08689890021523E-16,	4.19506351456677E-16,	5.69375484170172E-16,	7.71796875740479E-16,	1.04484348379901E-15,	1.41267861281199E-15,	1.90756501986615E-15,	2.5725221072992E-15,	3.46483537275146E-15,	4.66068501609413E-15,	6.26124209656813E-15,	8.40068634756923E-15,	1.12567366838138E-14,	1.50644625260262E-14,	2.01343727579731E-14,	2.68760740092946E-14,	3.582916960792E-14,	4.77035585984831E-14,	6.34319206789198E-14,	8.42379757928853E-14,	1.11725134551318E-13,	1.47991468436758E-13,	1.95778603548699E-13,	2.58664275034035E-13,	3.41310992698064E-13,	4.49786713951593E-13,	5.9197778909583E-13,	7.78119974738096E-13,	1.02148030201538E-12,	1.33923140895941E-12,	1.75357106332519E-12,	2.29315358820295E-12,	2.99491747048686E-12,	3.90641546424848E-12,	5.08878108795893E-12,	6.62049990205187E-12,	8.60219699899902E-12,	1.11627057797399E-11,	1.44667498290573E-11,	1.87246526024766E-11,	2.42045924457432E-11,	3.12480477592831E-11,	4.02892344681907E-11,	5.18795321427313E-11,	6.67181343657044E-11,	8.56904532351597E-11,	1.09916169287681E-10,	1.40809261192854E-10,	1.80152891853103E-10,	2.30192690150975E-10,	2.9375277585319E-10,	3.74379859459252E-10,	4.76521945362924E-10,	6.05749618883788E-10,	7.69029657097426E-10,	9.75062831243446E-10,	1.23470033752398E-09,	1.56145638937665E-09,	1.97213822946126E-09,	2.48761929344654E-09,	3.13378662113627E-09,	3.94270002891016E-09,	4.95400822453848E-09,	6.21667618744348E-09,	7.79108901822214E-09,	9.75161038754459E-09,	1.21896890424433E-08,	1.52176249696887E-08,	1.897312824741E-08,	2.36248288857415E-08,	2.93789256951777E-08,	3.64871971498006E-08,	4.52566381572313E-08,	5.60610345494614E-08,	6.93548430389403E-08,	8.5689809555945E-08,	1.05734834668221E-07,	1.302996827284E-07,	1.60363233236273E-07,	1.97107090587624E-07,	2.41955504081723E-07,	2.96622706192584E-07,	3.63168956354562E-07,	4.44066782829773E-07,	5.42279149869267E-07,	6.61351544918343E-07,	8.05520284939971E-07,	9.79839686263406E-07,	1.19033113327213E-06,	1.44415752263147E-06,	1.74982705686299E-06,	2.11743091946632E-06,	2.55891998933804E-06,	3.08842645108543E-06,	3.72263693635386E-06,	4.48122469602727E-06,	5.38734926077733E-06,	6.46823310459243E-06,	7.75582598918887E-06,	9.2875689429452E-06,	1.11072712221821E-05,	1.32661151206558E-05,	1.58238051397037E-05,	0.00001884987981024,	2.24252063711497E-05,	0.00002664368055741,	3.16141559344243E-05,	3.74626295292922E-05,	4.43347129490715E-05,	5.23984207315135E-05,	6.18473103326552E-05,	7.29040109014669E-05,	8.58241808033325E-05,	0.000100900936705338,	0.000118469799895386,	0.000138914208336198,	0.000162671645709822,	0.000190240441336642,	0.000222187297295956,	0.000259155601265509,	0.000301874585459455,	0.000351169393493627,	0.00040797211795069,	0.000473333871753659,	0.000548437956073796,	0.000634614186278792,	0.00073335443524235,	0.000846329450052448,	0.000975406993632119,	0.00112267135687785,	0.0012904442794777,	0.00148130730844311,	0.00169812561242513,	0.00194407325694054,	0.00222265993056453,	0.0025377590948236,	0.00289363751082551,	0.00329498607349048,	0.00374695185951989,	0.00425517126789866,	0.00482580410175243,	0.00546556840778054,	0.00618177585431291,	0.0069823673913907,	0.00787594889629625,	0.00887182646586039,	0.00998004097292257,	0.0112114014588398,	0.0125775168873351,	0.0140908257377163,	0.0157646228681269,	0.0176130830326376,	0.0196512803903425,	0.0218952033009686,	0.024361763660678,	0.027068799994662,	0.0300350734907684,	0.0332802561317892,	0.0368249100642666,	0.0406904573298339,	0.0448991390823394,	0.0494739634214452,	0.054438640992151,	0.0598175075308782,	0.0656354325833704,	0.0719177136787083,	0.0786899553180234,	0.0859779322267872,	0.0938074364264016,	0.10220410780464,	0.111193248005461,	0.120799617616832,	0.131047216810128,	0.141959049775888,	0.153556873507335,	0.165860931703912,	0.178889674800692,	0.192659467374028,	0.207184284427028,	0.222475398317956,	0.23854105835747,	0.255386165363778,	0.273011943724686,	0.291415613768501,	0.310590067487947,	0.330523550888371,	0.351199356439436,	0.372595529293766,	0.394684591092216,	0.417433285299343,	0.440802348099881,	0.464746308933637,	0.489213324748381,	0.514145052004557,	0.539476560369099,	0.565136291885586,	0.591046069202778,	0.617121156181683,	0.643270373882431,	0.669396274556464,	0.695395375838056,	0.721158456843762,	0.746570917352214,	0.771513200653001,	0.795861280027415,	0.819487208160658,	0.842259728091535,	0.86404494358855,	0.884707046109044,	0.90410909475876,	0.922113844932701,	0.938584620593488,	0.95338622444111,	0.966385879557766,	0.977454195483903,	0.986466151106528,	0.993302086228303,	0.997848693245159,	1,	0.999658334608501,	0.996735262877173,	0.99115248885963,	0.982842709129457,	0.971750411490379,	0.957832609098869,	0.941059501341448,	0.921415053287494,	0.898897486126068,	0.873519671687591,	0.845309424942457,	0.814309689251451,	0.780578610108199,	0.744189494151644,	0.70523065132514,	0.663805119205529,	0.620030269706927,	0.57403729956559,	0.525970607219161,	0.475987059890558,	0.42425515585854,	0.370954088028045,	0.316272715989177,	0.260408454759576,	0.203566089327616,	0.145956524940774,	0.0877954838034294,	0.0293021594514236,	-0.0293021594514245,	-0.0877954838034304,	-0.145956524940774,	-0.203566089327617,	-0.260408454759578,	-0.316272715989178,	-0.370954088028045,	-0.424255155858542,	-0.47598705989056,	-0.525970607219161,	-0.57403729956559,	-0.620030269706927,	-0.663805119205531,	-0.705230651325142,	-0.744189494151644,	-0.780578610108199,	-0.814309689251451,	-0.845309424942458,	-0.873519671687591,	-0.898897486126069,	-0.921415053287494,	-0.941059501341448,	-0.957832609098869,	-0.971750411490379,	-0.982842709129457,	-0.99115248885963,	-0.996735262877174,	-0.999658334608503,	-1,	-0.997848693245161,	-0.993302086228305,	-0.986466151106529,	-0.977454195483903,	-0.966385879557766,	-0.95338622444111,	-0.938584620593488,	-0.922113844932701,	-0.904109094758762,	-0.884707046109044,	-0.86404494358855,	-0.842259728091535,	-0.81948720816066,	-0.795861280027415,	-0.771513200653001,	-0.746570917352214,	-0.721158456843762,	-0.695395375838056,	-0.669396274556466,	-0.643270373882431,	-0.617121156181683,	-0.591046069202778,	-0.565136291885586,	-0.539476560369099,	-0.514145052004557,	-0.489213324748381,	-0.464746308933637,	-0.440802348099881,	-0.417433285299343,	-0.394684591092216,	-0.372595529293766,	-0.351199356439436,	-0.330523550888371,	-0.310590067487949,	-0.291415613768501,	-0.273011943724686,	-0.255386165363778,	-0.23854105835747,	-0.222475398317956,	-0.207184284427028,	-0.192659467374028,	-0.178889674800692,	-0.165860931703912,	-0.153556873507335,	-0.141959049775888,	-0.131047216810128,	-0.120799617616832,	-0.111193248005461,	-0.10220410780464,	-0.0938074364264019,	-0.0859779322267875,	-0.0786899553180237,	-0.0719177136787086,	-0.0656354325833708,	-0.0598175075308785,	-0.0544386409921513,	-0.0494739634214455,	-0.0448991390823399,	-0.0406904573298344,	-0.0368249100642671,	-0.0332802561317895,	-0.0300350734907687,	-0.0270687999946623,	-0.0243617636606782,	-0.021895203300969,	-0.0196512803903429,	-0.0176130830326379,	-0.0157646228681272,	-0.0140908257377166,	-0.0125775168873355,	-0.0112114014588401,	-0.00998004097292293,	-0.00887182646586075,	-0.00787594889629661,	-0.00698236739139106,	-0.00618177585431327,	-0.0054655684077809,	-0.00482580410175279,	-0.00425517126789902,	-0.00374695185952026,	-0.00329498607349084,	-0.00289363751082587,	-0.00253775909482396,	-0.00222265993056489,	-0.00194407325694091,	-0.00169812561242549,	-0.00148130730844348,	-0.00129044427947807,	-0.00112267135687822,	-0.000975406993632486,	-0.000846329450052815,	-0.000733354435242717,	-0.000634614186279159,	-0.000548437956074162,	-0.000473333871754026,	-0.000407972117951058,	-0.000351169393493994,	-0.000301874585459822,	-0.000259155601265877,	-0.000222187297296324,	-0.00019024044133701,	-0.00016267164571019,	-0.000138914208336566,	-0.000118469799895754,	-0.000100900936705706,	-8.58241808037005E-05,	-7.29040109018349E-05,	-6.18473103330234E-05,	-5.23984207318817E-05,	-4.43347129494397E-05,	-3.74626295296604E-05,	-3.16141559347925E-05,	-2.66436805577782E-05,	-2.24252063715181E-05,	-1.88498798106082E-05,	-1.58238051400721E-05,	-1.32661151210241E-05,	-1.11072712225503E-05,	-9.28756894331343E-06,	-7.75582598955712E-06,	-6.46823310496068E-06,	-5.38734926114557E-06,	-4.48122469639552E-06,	-3.72263693672212E-06,	-3.08842645145369E-06,	-2.55891998970628E-06,	-2.11743091983458E-06,	-1.74982705723125E-06,	-1.44415752299972E-06,	-1.19033113364039E-06,	-9.79839686631661E-07,	-8.05520285308228E-07,	-6.61351545286602E-07,	-5.42279150237526E-07,	-4.44066783198031E-07,	-3.6316895672282E-07,	-2.96622706560844E-07,	-2.41955504449981E-07,	-1.97107090955884E-07,	-1.60363233604531E-07,	-1.30299683096659E-07,	-1.0573483503648E-07,	-8.56898099242044E-08,	-6.93548434071995E-08,	-5.60610349177208E-08,	-4.52566385254907E-08,	-3.648719751806E-08,	-2.93789260634371E-08,	-2.36248292540009E-08,	-1.89731286156694E-08,	-1.52176253379481E-08,	-1.21896894107027E-08,	-9.75161075580404E-09,	-7.79108938648158E-09,	-6.21667655570292E-09,	-4.95400859279792E-09,	-3.9427003971696E-09,	-3.13378698939575E-09,	-2.487619661706E-09,	-1.97213859772073E-09,	-1.56145675763611E-09,	-1.23470070578345E-09,	-9.75063199502915E-10,	-7.69030025356895E-10,	-6.05749987143255E-10,	-4.76522313622393E-10,	-3.74380227718723E-10,	-2.93753144112661E-10,	-2.30193058410446E-10,	-1.80153260112575E-10,	-1.40809629452325E-10,	-1.09916537547152E-10,	-8.56908214946309E-11,	-6.67185026251758E-11,	-5.18799004022028E-11,	-4.0289602727662E-11,	-3.12484160187545E-11,	-2.42049607052147E-11,	-1.87250208619481E-11,	-1.44671180885289E-11,	-1.11630740392114E-11,	-8.60256525847056E-12,	-6.62086816152344E-12,	-5.0891493474305E-12,	-3.90678372372004E-12,	-2.99528572995842E-12,	-2.29352184767451E-12,	-1.75393932279676E-12,	-1.33959966843097E-12,	-1.02184856148694E-12,	-7.78488234209659E-13,	-5.92346048567396E-13,	-4.50154973423158E-13,	-3.4167925216963E-13,	-2.590325345056E-13,	-1.96146863020266E-13,	-1.48359727908323E-13,	-1.12093394022884E-13,	-8.46062352644517E-14,	-6.38001801504864E-14,	-4.80718180700495E-14,	-3.61974290794863E-14,	-2.72443334808611E-14,	-2.05026322295396E-14,	-1.54327219975927E-14,	-1.16249961553804E-14,	-8.76894581913584E-15,	-6.62950156813477E-15,	-5.02894448766078E-15,	-3.83309484431811E-15,	-2.94078157886584E-15,	-2.2758244914328E-15,	-1.78093808437863E-15,	-1.41310295536566E-15,	-1.14005634730712E-15,	-9.37634955736819E-16,	-7.87765823023325E-16,	-6.76949361588171E-16,	-5.95115548904751E-16,	-5.34762690793086E-16,	-4.90309844892551E-16,	-4.57610572210875E-16,	-4.33588330582687E-16,	-4.15963457289721E-16,	-4.03049015830059E-16,	-3.93598347525974E-16,	-3.8669138982275E-16,	-3.81650019869991E-16,	-3.77975099364223E-16,	-3.75299721998947E-16,	-3.73354541190686E-16,	-3.71942091981665E-16,	-3.70917800063109E-16,	-3.7017595570036E-16,	-3.69639368748613E-16,	-3.69251749112448E-16,	-3.68972102290475E-16,	-3.687706127233E-16,	-3.68625624110026E-16,	-3.68521427405681E-16,	-3.68446642677248E-16,	-3.68393036997636E-16,	-3.68354662055488E-16,	-3.68327225865799E-16,	-3.68307635655955E-16,	-3.68293665743037E-16,	-3.68283716553145E-16,	-3.6827664000867E-16,	-3.68271613176926E-16,	-3.68268046965124E-16,	-3.68265520230226E-16,	-3.68263732293804E-16,	-3.68262468767182E-16,	-3.68261576989217E-16,	-3.68260948396746E-16,	-3.68260505888067E-16,	-3.68260194777565E-16,	-3.68259976329656E-16,	-3.68259823142789E-16,	-3.68259715858614E-16,	-3.68259640819036E-16,	-3.68259588400396E-16,	-3.6825955183066E-16,	-3.682595263507E-16,	-3.68259508620383E-16,	-3.68259496298553E-16,	-3.68259487746409E-16,	-3.68259481818306E-16,	-3.68259477714395E-16,	-3.68259474876992E-16,	-3.68259472917762E-16,	-3.68259471566649E-16,	-3.68259470636101E-16,	-3.6825946999603E-16,	-3.68259469556328E-16,	-3.6825946925466E-16,	-3.68259469047959E-16,	-3.68259468906509E-16,	-3.68259468809838E-16,	-3.68259468743854E-16,	-3.68259468698874E-16,	-3.68259468668252E-16,	-3.68259468647431E-16,	-3.68259468633291E-16,	-3.68259468623704E-16,	-3.68259468617208E-16,	-3.68259468612815E-16,	-3.68259468609847E-16,	-3.68259468607844E-16,	-3.68259468606495E-16,	-3.68259468605587E-16,	-3.68259468604977E-16,	-3.68259468604566E-16,	-3.68259468604293E-16,	-3.68259468604109E-16,	-3.68259468603986E-16,	-3.68259468603903E-16,	-3.68259468603848E-16,	-3.68259468603812E-16,	-3.68259468603788E-16,	-3.68259468603772E-16,	-3.6825946860376E-16,	-3.68259468603754E-16,	-3.68259468603749E-16,	-3.68259468603746E-16,	-3.68259468603744E-16,	-3.68259468603743E-16,	-3.68259468603741E-16,	-3.68259468603741E-16,	-3.68259468603741E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16};


	double force_1 = Force_1[timestep];
	double force_2 = -1*Force_1[timestep+5];

	double force_coefficient =    3538944;

	double force_1x = force_1 * force_coefficient;
	double force_1z = force_1 * force_coefficient/2;

	double force_2x = force_2 * force_coefficient;
	double force_2z = force_2 * force_coefficient/2;



	int32_t nindex;
	int32_t k1=0,k2=0;

	double f_l_depth = 512;   //first layer depth
	double el_size   = 32;     //element size
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

		if ( z_m == f_l_depth && x_m == 0 ){

			fvector_t *nodalForce;
			nodalForce = mySolver->force + nindex;
			nodalForce->f[2] += force_1z;
		}

		if ( z_m == s_l_depth && x_m == 0 ){

			fvector_t *nodalForce;
			nodalForce = mySolver->force + nindex;
			nodalForce->f[2] += force_2z;
		}

		if ( z_m == f_l_depth && x_m == d_width_x ){

			fvector_t *nodalForce;
			nodalForce = mySolver->force + nindex;
			nodalForce->f[2] -= force_1z;
		}

		if ( z_m == s_l_depth && x_m == d_width_x ){

			fvector_t *nodalForce;
			nodalForce = mySolver->force + nindex;
			nodalForce->f[2] -= force_2z;
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

