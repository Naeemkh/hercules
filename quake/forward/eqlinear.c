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

//double fc =0.8,zp=0.04,Vs=640.0,Ts=3.0;
	//double force_coef = 1;
	//double ri_dt = 0.005;
    //double Force_1[1000]={9.33888995979801E-17,	2.28810934846276E-16,	4.24932246600079E-16,	7.08594499517493E-16,	1.11834504686759E-15,	1.70946919697417E-15,	2.56115411662058E-15,	3.78667364254051E-15,	5.54784844408368E-15,	8.07554867320022E-15,	1.16987250432901E-14,	1.68854606506188E-14,	2.43009428670989E-14,	3.48892190392015E-14,	4.998833849764E-14,	7.14922964202729E-14,	1.02078496329936E-13,	1.45526807942062E-13,	2.0716650083E-13,	2.94501460891339E-13,	4.18083626434948E-13,	5.92731542579966E-13,	8.39228102329984E-13,	1.18668271252894E-12,	1.67581346665072E-12,	2.36350147194507E-12,	3.32910185283552E-12,	4.68317715724979E-12,	6.57956658850295E-12,	9.23204258983848E-12,	1.29372671645931E-11,	1.81063875939534E-11,	2.53084637981968E-11,	3.53300767477262E-11,	4.9257035630866E-11,	6.85862239656563E-11,	9.53784932173251E-11,	1.32467383322103E-10,	1.83743665675347E-10,	2.54542723085975E-10,	3.52171212413388E-10,	4.86622109880422E-10,	6.71544193368363E-10,	9.25554698727023E-10,	1.27401418986177E-09,	1.75142310630858E-09,	2.40465329948698E-09,	3.29729949317873E-09,	4.51553073518854E-09,	6.17594807652895E-09,	8.43612180428634E-09,	1.15087025129314E-08,	1.56802920656097E-08,	2.13366451982877E-08,	2.89962790126855E-08,	3.93552334658423E-08,	5.33466000511819E-08,	7.22195816001876E-08,	9.76443457175644E-08,	1.31850894173148E-07,	1.77812728056521E-07,	2.39489420964581E-07,	3.22146548312426E-07,	4.32777053714789E-07,	5.80655440343985E-07,	7.78065639585209E-07,	1.04125564979893E-06,	1.3916878450311E-06,	1.85767423552244E-06,	2.47651220329943E-06,	3.29726997032715E-06,	4.38441435732694E-06,	5.82252888680423E-06,	7.72244130018708E-06,	1.02291702983487E-05,	1.35322170772861E-05,	1.78788746903525E-05,	2.35914158152902E-05,	3.10892576683341E-05,	4.09175047823038E-05,	5.37836526538972E-05,	7.06047184753032E-05,	9.25676749872742E-05,	0.000121206831897037,	0.000158502776016588,	0.000207008695560052,	0.000270011436888855,	0.000351736548799721,	0.000457608953229825,	0.000594583856570812,	0.00077156622349562,	0.000999941748184442,	0.00129424798788318,	0.00167302142984977,	0.00215986506112165,	0.00278479188644134,	0.00358591326170559,	0.00461155744735749,	0.00592292412761424,	0.00759740562116449,	0.00973273613365008,	0.0124521678848865,	0.0159109187432354,	0.0203041918669725,	0.0258771358858235,	0.0329371968634853,	0.0418694136572474,	0.0531553298985741,	0.067396342895547,	0.0853424873296573,	0.107927865633714,	0.136314194420362,	0.171944245563025,	0.216607331252939,	0.27251942599827,	0.3424210484801,	0.429696658076093,	0.538520072939514,	0.674031309987082,	0.842551306613683,	1.05184223790488,	1.31142262446098,	1.63294817258245,	2.03067134405068,	2.52199506699422,	3.12813882942783,	3.87493870805877,	4.79380675085155,	5.92287963650692,	7.30839177226742,	9.00631407018702,	11.0843066812912,	13.6240421014166,	16.7239644421287,	20.5025614520314,	25.1022382632273,	30.6938960291939,	37.4823348389637,	45.7126187847559,	55.6775620950508,	67.7265191145266,	82.2756879317699,	99.820167964287,	120.94804617037,	146.356825153512,	176.872549663045,	213.472036299025,	257.308665040779,	309.742250989405,	372.373580904929,	447.084272189429,	536.082692372735,	641.956765336372,	767.734586877227,	916.953877139732,	1093.74141125632,	1302.90369248379,	1550.03026437387,	1841.61120012818,	2185.17045818741,	2589.41695306468,	3064.41535903982,	3621.77884096221,	4274.88609019313,	5039.12523350939,	5932.16737712634,	6974.27274507297,	8188.63256876907,	9601.75008019803,	11243.8641514629,	13149.4193051887,	15357.5859890965,	17912.835159472,	20865.5713469575,	24272.8284782795,	28199.0327927517,	32716.8372156129,	37908.0315238208,	43864.5325555901,	50689.4585639512,	58498.2915876252,	67420.1313998521,	77599.0441873972,	89195.5085974989,	102387.961159588,	117374.442330825,	134374.34351973,	153630.25440062,	175409.908634207,	200008.224748259,	227749.437399662,	258989.312530015,	294117.438037155,	333559.579513128,	377780.088345791,	427284.347050108,	482621.234092925,	544385.587711997,	613220.64532027,	689820.432048408,	774932.068835006,	869357.967252602,	973957.874990949,	1089650.73264493,	1217416.29921591,	1358296.50058047,	1513396.45216295,	1683885.10422606,	1870995.45563104,	2076024.27968191,	2300331.30382927,	2545337.78364211,	2812524.41063812,	3103428.4933713,	3419640.35169029,	3762798.86537748,	4134586.1205343,	4536721.10016256,	4970952.36947232,	5439049.71158178,	5942794.67551553,	6483970.00579288,	7064347.93145671,	7685677.30213747,	8349669.56967544,	9057983.62591603,	9812209.52050936,	10613851.096827,	11464307.5993744,	12364854.3222238,	13316622.3848928,	14320577.7395962,	15377499.5317371,	16487957.9536683,	17652291.7499443,	18870585.5502503,	20142647.2236788,	21467985.4647669,	22845787.8374042,	24274899.5170938,	25753802.9847851,	27280598.936294,	28852988.6798906,	30468258.3006638,	32123264.873493,	33814425.0066081,	35537705.994555,	37288619.8527121,	39062220.4951317,	40853104.303296,	42655414.3152779,	44462848.2427536,	46268670.4973428,	48065728.3779264,	49846472.5370408,	51602981.807385,	53326992.4291354,	55009931.6754949,	56642955.8280647,	58216992.4056869,	59722786.5008406,	61150951.0270571,	62492020.6297255,	63736508.9617483,	64874968.9754219,	65898055.8333695,	66796591.9950328,	67561633.9918474,	68184540.3644832,	68657040.2001003,	68971301.6771054,	69120000,	69096384.0881396,	68894341.3700702,	68508460.0299776,	67934088.0550281,	67167388.442215,	66205389.9409138,	65046032.7327209,	63688208.4832316,	62131794.2410338,	60377679.7070463,	58427787.4520226,	56285085.7210603,	53953593.5306787,	51438377.8357616,	48745542.6195937,	45882209.8394862,	42856492.2421428,	39677458.1459736,	36355088.3709884,	32900225.5796354,	29324516.3729423,	25640346.5644985,	21860770.1291719,	17999432.3929819,	14070488.0943248,	10088515.0039063,	6068423.84049304,	2025365.2612824,	-2025365.26128246,	-6068423.84049311,	-10088515.0039063,	-14070488.0943249,	-17999432.392982,	-21860770.129172,	-25640346.5644985,	-29324516.3729424,	-32900225.5796355,	-36355088.3709884,	-39677458.1459736,	-42856492.2421428,	-45882209.8394863,	-48745542.6195938,	-51438377.8357616,	-53953593.5306787,	-56285085.7210603,	-58427787.4520227,	-60377679.7070463,	-62131794.2410339,	-63688208.4832316,	-65046032.7327209,	-66205389.9409138,	-67167388.442215,	-67934088.0550281,	-68508460.0299776,	-68894341.3700703,	-69096384.0881397,	-69120000,	-68971301.6771055,	-68657040.2001004,	-68184540.3644833,	-67561633.9918474,	-66796591.9950328,	-65898055.8333695,	-64874968.9754219,	-63736508.9617483,	-62492020.6297256,	-61150951.0270571,	-59722786.5008406,	-58216992.4056869,	-56642955.8280648,	-55009931.6754949,	-53326992.4291354,	-51602981.807385,	-49846472.5370408,	-48065728.3779264,	-46268670.4973429,	-44462848.2427536,	-42655414.3152779,	-40853104.303296,	-39062220.4951317,	-37288619.8527121,	-35537705.994555,	-33814425.0066081,	-32123264.873493,	-30468258.3006638,	-28852988.6798906,	-27280598.936294,	-25753802.9847851,	-24274899.5170938,	-22845787.8374042,	-21467985.464767,	-20142647.2236788,	-18870585.5502503,	-17652291.7499443,	-16487957.9536683,	-15377499.5317371,	-14320577.7395962,	-13316622.3848928,	-12364854.3222238,	-11464307.5993744,	-10613851.096827,	-9812209.52050937,	-9057983.62591605,	-8349669.56967545,	-7685677.30213748,	-7064347.93145673,	-6483970.0057929,	-5942794.67551555,	-5439049.7115818,	-4970952.36947234,	-4536721.10016259,	-4134586.12053432,	-3762798.8653775,	-3419640.35169031,	-3103428.49337133,	-2812524.41063815,	-2545337.78364214,	-2300331.30382929,	-2076024.27968193,	-1870995.45563106,	-1683885.10422608,	-1513396.45216298,	-1358296.5005805,	-1217416.29921593,	-1089650.73264495,	-973957.874990973,	-869357.967252627,	-774932.06883503,	-689820.432048433,	-613220.645320295,	-544385.587712022,	-482621.23409295,	-427284.347050133,	-377780.088345816,	-333559.579513153,	-294117.43803718,	-258989.31253004,	-227749.437399687,	-200008.224748284,	-175409.908634232,	-153630.254400645,	-134374.343519756,	-117374.44233085,	-102387.961159613,	-89195.5085975243,	-77599.0441874225,	-67420.1313998774,	-58498.2915876506,	-50689.4585639766,	-43864.5325556155,	-37908.0315238461,	-32716.8372156383,	-28199.0327927771,	-24272.8284783049,	-20865.5713469829,	-17912.8351594974,	-15357.5859891219,	-13149.4193052141,	-11243.8641514883,	-9601.75008022345,	-8188.6325687945,	-6974.27274509841,	-5932.16737715178,	-5039.12523353483,	-4274.88609021858,	-3621.77884098766,	-3064.41535906527,	-2589.41695309013,	-2185.17045821286,	-1841.61120015363,	-1550.03026439933,	-1302.90369250924,	-1093.74141128178,	-916.953877165185,	-767.73458690268,	-641.956765361824,	-536.082692398188,	-447.084272214882,	-372.373580930382,	-309.742251014858,	-257.308665066233,	-213.472036324479,	-176.872549688498,	-146.356825178966,	-120.948046195824,	-99.8201679897409,	-82.2756879572237,	-67.7265191399804,	-55.6775621205047,	-45.7126188102099,	-37.4823348644178,	-30.6938960546479,	-25.1022382886813,	-20.5025614774855,	-16.7239644675827,	-13.6240421268707,	-11.0843067067452,	-9.0063140956411,	-7.3083917977215,	-5.92287966196101,	-4.79380677630563,	-3.87493873351286,	-3.12813885488192,	-2.52199509244831,	-2.03067136950477,	-1.63294819803654,	-1.31142264991507,	-1.05184226335897,	-0.842551332067774,	-0.674031335441175,	-0.538520098393607,	-0.429696683530186,	-0.342421073934192,	-0.272519451452363,	-0.216607356707034,	-0.171944271017119,	-0.136314219874457,	-0.107927891087808,	-0.0853425127837518,	-0.0673963683496415,	-0.0531553553526686,	-0.0418694391113418,	-0.0329372223175798,	-0.0258771613399181,	-0.0203042173210671,	-0.01591094419733,	-0.0124521933389812,	-0.00973276158774472,	-0.00759743107525915,	-0.00592294958170889,	-0.00461158290145215,	-0.00358593871580026,	-0.002784817340536,	-0.00215989051521631,	-0.00167304688394444,	-0.00129427344197785,	-0.000999967202279115,	-0.000771591677590291,	-0.000594609310665485,	-0.0004576344073245,	-0.000351762002894396,	-0.000270036890983529,	-0.000207034149654726,	-0.000158528230111262,	-0.000121232285991712,	-9.25931290819487E-05,	-7.06301725699776E-05,	-5.38091067485716E-05,	-4.09429588769784E-05,	-3.11147117630087E-05,	-2.36168699099648E-05,	-1.79043287850271E-05,	-1.35576711719608E-05,	-1.02546243930233E-05,	-7.74789539486173E-06,	-5.8479829814789E-06,	-4.40986845200162E-06,	-3.32272406500182E-06,	-2.50196629797409E-06,	-1.88312833019712E-06,	-1.41714193970578E-06,	-1.06670974447361E-06,	-8.03519734259891E-07,	-6.06109535018669E-07,	-4.58231148389475E-07,	-3.47600642987113E-07,	-2.64943515639268E-07,	-2.03266822731207E-07,	-1.57304988847835E-07,	-1.23098440392251E-07,	-9.76736762748742E-08,	-7.88006947258683E-08,	-6.48093281405289E-08,	-5.44503736873722E-08,	-4.67907398729744E-08,	-4.11343867402964E-08,	-3.69627971876181E-08,	-3.38902164789731E-08,	-3.16300427512157E-08,	-2.99696254098753E-08,	-2.87513941678655E-08,	-2.78587479741737E-08,	-2.72055177809953E-08,	-2.67281088645485E-08,	-2.63796493734138E-08,	-2.61256388680551E-08,	-2.59407167845672E-08,	-2.58062658871002E-08,	-2.57086373977727E-08,	-2.56378383403621E-08,	-2.55865620580089E-08,	-2.55494731679041E-08,	-2.55226808986524E-08,	-2.55033517103176E-08,	-2.54894247514345E-08,	-2.5479403138485E-08,	-2.54722010622807E-08,	-2.54670319418514E-08,	-2.54633267172766E-08,	-2.54606742412753E-08,	-2.5458777851844E-08,	-2.54574237765396E-08,	-2.54564581761587E-08,	-2.54557704881534E-08,	-2.54552813573993E-08,	-2.54549339027891E-08,	-2.54546874062294E-08,	-2.54545127583132E-08,	-2.54543891761477E-08,	-2.54543018411876E-08,	-2.54542402014947E-08,	-2.54541967531831E-08,	-2.54541661669832E-08,	-2.54541446630253E-08,	-2.54541295639058E-08,	-2.54541189756296E-08,	-2.54541115601474E-08,	-2.54541063734118E-08,	-2.54541027502354E-08,	-2.54541002225352E-08,	-2.54540984613604E-08,	-2.54540972358409E-08,	-2.5454096384156E-08,	-2.54540957930318E-08,	-2.54540953832813E-08,	-2.5454095099619E-08,	-2.54540949034977E-08,	-2.54540947680757E-08,	-2.54540946746868E-08,	-2.54540946103673E-08,	-2.54540945661256E-08,	-2.54540945357334E-08,	-2.54540945148821E-08,	-2.54540945005949E-08,	-2.54540944908179E-08,	-2.5454094484136E-08,	-2.54540944795752E-08,	-2.54540944764662E-08,	-2.54540944743496E-08,	-2.54540944729104E-08,	-2.54540944719331E-08,	-2.54540944712704E-08,	-2.54540944708214E-08,	-2.54540944705178E-08,	-2.54540944703126E-08,	-2.54540944701742E-08,	-2.54540944700809E-08,	-2.54540944700182E-08,	-2.5454094469976E-08,	-2.54540944699476E-08,	-2.54540944699287E-08,	-2.5454094469916E-08,	-2.54540944699075E-08,	-2.54540944699018E-08,	-2.5454094469898E-08,	-2.54540944698955E-08,	-2.54540944698938E-08,	-2.54540944698927E-08,	-2.54540944698919E-08,	-2.54540944698915E-08,	-2.54540944698911E-08,	-2.54540944698909E-08,	-2.54540944698908E-08,	-2.54540944698907E-08,	-2.54540944698906E-08,	-2.54540944698906E-08,	-2.54540944698906E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08,	-2.54540944698905E-08};

    // dt = 0.01
    double Force_1[1000]={1.35111255205411E-24,	3.31034338608617E-24,	6.14774662326503E-24,	1.02516565323711E-23,	1.61797605160242E-23,	2.47319038914087E-23,	3.70537343261079E-23,	5.47840515413847E-23,	8.02640110544514E-23,	1.16833748165512E-22,	1.69252387779081E-22,	2.44291965431406E-22,	3.51576141017056E-22,	5.04763006932892E-22,	7.23210915764468E-22,	1.03432141811737E-21,	1.47683009736597E-21,	2.10542256860622E-21,	2.99720053284144E-21,	4.26072715409923E-21,	6.0486635768945E-21,	8.5753984748259E-21,	1.21416102767648E-20,	1.71684420215414E-20,	2.42449864966829E-20,	3.41941763880942E-20,	4.81640893060695E-20,	6.77542991500259E-20,	9.5190488838295E-20,	1.33565430987247E-19,	1.87171110598859E-19,	2.61955839032891E-19,	3.6615254337669E-19,	5.11141156651131E-19,	7.12630723826186E-19,	9.92277545799426E-19,	1.37989718196362E-18,	1.91648413371098E-18,	2.65832849646046E-18,	3.682620415017E-18,	5.09506962403629E-18,	7.04025043229777E-18,	9.7156277975747E-18,	1.33905483033423E-17,	1.8431918256102E-17,	2.53388759593255E-17,	3.47895442634112E-17,	4.77039857230719E-17,	6.53288590160379E-17,	8.93511006442267E-17,	1.22050373325902E-16,	1.66503219226438E-16,	2.26856077338103E-16,	3.08689890021523E-16,	4.19506351456677E-16,	5.69375484170172E-16,	7.71796875740479E-16,	1.04484348379901E-15,	1.41267861281199E-15,	1.90756501986615E-15,	2.5725221072992E-15,	3.46483537275146E-15,	4.66068501609413E-15,	6.26124209656813E-15,	8.40068634756923E-15,	1.12567366838138E-14,	1.50644625260262E-14,	2.01343727579731E-14,	2.68760740092946E-14,	3.582916960792E-14,	4.77035585984831E-14,	6.34319206789198E-14,	8.42379757928853E-14,	1.11725134551318E-13,	1.47991468436758E-13,	1.95778603548699E-13,	2.58664275034035E-13,	3.41310992698064E-13,	4.49786713951593E-13,	5.9197778909583E-13,	7.78119974738096E-13,	1.02148030201538E-12,	1.33923140895941E-12,	1.75357106332519E-12,	2.29315358820295E-12,	2.99491747048686E-12,	3.90641546424848E-12,	5.08878108795893E-12,	6.62049990205187E-12,	8.60219699899902E-12,	1.11627057797399E-11,	1.44667498290573E-11,	1.87246526024766E-11,	2.42045924457432E-11,	3.12480477592831E-11,	4.02892344681907E-11,	5.18795321427313E-11,	6.67181343657044E-11,	8.56904532351597E-11,	1.09916169287681E-10,	1.40809261192854E-10,	1.80152891853103E-10,	2.30192690150975E-10,	2.9375277585319E-10,	3.74379859459252E-10,	4.76521945362924E-10,	6.05749618883788E-10,	7.69029657097426E-10,	9.75062831243446E-10,	1.23470033752398E-09,	1.56145638937665E-09,	1.97213822946126E-09,	2.48761929344654E-09,	3.13378662113627E-09,	3.94270002891016E-09,	4.95400822453848E-09,	6.21667618744348E-09,	7.79108901822214E-09,	9.75161038754459E-09,	1.21896890424433E-08,	1.52176249696887E-08,	1.897312824741E-08,	2.36248288857415E-08,	2.93789256951777E-08,	3.64871971498006E-08,	4.52566381572313E-08,	5.60610345494614E-08,	6.93548430389403E-08,	8.5689809555945E-08,	1.05734834668221E-07,	1.302996827284E-07,	1.60363233236273E-07,	1.97107090587624E-07,	2.41955504081723E-07,	2.96622706192584E-07,	3.63168956354562E-07,	4.44066782829773E-07,	5.42279149869267E-07,	6.61351544918343E-07,	8.05520284939971E-07,	9.79839686263406E-07,	1.19033113327213E-06,	1.44415752263147E-06,	1.74982705686299E-06,	2.11743091946632E-06,	2.55891998933804E-06,	3.08842645108543E-06,	3.72263693635386E-06,	4.48122469602727E-06,	5.38734926077733E-06,	6.46823310459243E-06,	7.75582598918887E-06,	9.2875689429452E-06,	1.11072712221821E-05,	1.32661151206558E-05,	1.58238051397037E-05,	0.00001884987981024,	2.24252063711497E-05,	0.00002664368055741,	3.16141559344243E-05,	3.74626295292922E-05,	4.43347129490715E-05,	5.23984207315135E-05,	6.18473103326552E-05,	7.29040109014669E-05,	8.58241808033325E-05,	0.000100900936705338,	0.000118469799895386,	0.000138914208336198,	0.000162671645709822,	0.000190240441336642,	0.000222187297295956,	0.000259155601265509,	0.000301874585459455,	0.000351169393493627,	0.00040797211795069,	0.000473333871753659,	0.000548437956073796,	0.000634614186278792,	0.00073335443524235,	0.000846329450052448,	0.000975406993632119,	0.00112267135687785,	0.0012904442794777,	0.00148130730844311,	0.00169812561242513,	0.00194407325694054,	0.00222265993056453,	0.0025377590948236,	0.00289363751082551,	0.00329498607349048,	0.00374695185951989,	0.00425517126789866,	0.00482580410175243,	0.00546556840778054,	0.00618177585431291,	0.0069823673913907,	0.00787594889629625,	0.00887182646586039,	0.00998004097292257,	0.0112114014588398,	0.0125775168873351,	0.0140908257377163,	0.0157646228681269,	0.0176130830326376,	0.0196512803903425,	0.0218952033009686,	0.024361763660678,	0.027068799994662,	0.0300350734907684,	0.0332802561317892,	0.0368249100642666,	0.0406904573298339,	0.0448991390823394,	0.0494739634214452,	0.054438640992151,	0.0598175075308782,	0.0656354325833704,	0.0719177136787083,	0.0786899553180234,	0.0859779322267872,	0.0938074364264016,	0.10220410780464,	0.111193248005461,	0.120799617616832,	0.131047216810128,	0.141959049775888,	0.153556873507335,	0.165860931703912,	0.178889674800692,	0.192659467374028,	0.207184284427028,	0.222475398317956,	0.23854105835747,	0.255386165363778,	0.273011943724686,	0.291415613768501,	0.310590067487947,	0.330523550888371,	0.351199356439436,	0.372595529293766,	0.394684591092216,	0.417433285299343,	0.440802348099881,	0.464746308933637,	0.489213324748381,	0.514145052004557,	0.539476560369099,	0.565136291885586,	0.591046069202778,	0.617121156181683,	0.643270373882431,	0.669396274556464,	0.695395375838056,	0.721158456843762,	0.746570917352214,	0.771513200653001,	0.795861280027415,	0.819487208160658,	0.842259728091535,	0.86404494358855,	0.884707046109044,	0.90410909475876,	0.922113844932701,	0.938584620593488,	0.95338622444111,	0.966385879557766,	0.977454195483903,	0.986466151106528,	0.993302086228303,	0.997848693245159,	1,	0.999658334608501,	0.996735262877173,	0.99115248885963,	0.982842709129457,	0.971750411490379,	0.957832609098869,	0.941059501341448,	0.921415053287494,	0.898897486126068,	0.873519671687591,	0.845309424942457,	0.814309689251451,	0.780578610108199,	0.744189494151644,	0.70523065132514,	0.663805119205529,	0.620030269706927,	0.57403729956559,	0.525970607219161,	0.475987059890558,	0.42425515585854,	0.370954088028045,	0.316272715989177,	0.260408454759576,	0.203566089327616,	0.145956524940774,	0.0877954838034294,	0.0293021594514236,	-0.0293021594514245,	-0.0877954838034304,	-0.145956524940774,	-0.203566089327617,	-0.260408454759578,	-0.316272715989178,	-0.370954088028045,	-0.424255155858542,	-0.47598705989056,	-0.525970607219161,	-0.57403729956559,	-0.620030269706927,	-0.663805119205531,	-0.705230651325142,	-0.744189494151644,	-0.780578610108199,	-0.814309689251451,	-0.845309424942458,	-0.873519671687591,	-0.898897486126069,	-0.921415053287494,	-0.941059501341448,	-0.957832609098869,	-0.971750411490379,	-0.982842709129457,	-0.99115248885963,	-0.996735262877174,	-0.999658334608503,	-1,	-0.997848693245161,	-0.993302086228305,	-0.986466151106529,	-0.977454195483903,	-0.966385879557766,	-0.95338622444111,	-0.938584620593488,	-0.922113844932701,	-0.904109094758762,	-0.884707046109044,	-0.86404494358855,	-0.842259728091535,	-0.81948720816066,	-0.795861280027415,	-0.771513200653001,	-0.746570917352214,	-0.721158456843762,	-0.695395375838056,	-0.669396274556466,	-0.643270373882431,	-0.617121156181683,	-0.591046069202778,	-0.565136291885586,	-0.539476560369099,	-0.514145052004557,	-0.489213324748381,	-0.464746308933637,	-0.440802348099881,	-0.417433285299343,	-0.394684591092216,	-0.372595529293766,	-0.351199356439436,	-0.330523550888371,	-0.310590067487949,	-0.291415613768501,	-0.273011943724686,	-0.255386165363778,	-0.23854105835747,	-0.222475398317956,	-0.207184284427028,	-0.192659467374028,	-0.178889674800692,	-0.165860931703912,	-0.153556873507335,	-0.141959049775888,	-0.131047216810128,	-0.120799617616832,	-0.111193248005461,	-0.10220410780464,	-0.0938074364264019,	-0.0859779322267875,	-0.0786899553180237,	-0.0719177136787086,	-0.0656354325833708,	-0.0598175075308785,	-0.0544386409921513,	-0.0494739634214455,	-0.0448991390823399,	-0.0406904573298344,	-0.0368249100642671,	-0.0332802561317895,	-0.0300350734907687,	-0.0270687999946623,	-0.0243617636606782,	-0.021895203300969,	-0.0196512803903429,	-0.0176130830326379,	-0.0157646228681272,	-0.0140908257377166,	-0.0125775168873355,	-0.0112114014588401,	-0.00998004097292293,	-0.00887182646586075,	-0.00787594889629661,	-0.00698236739139106,	-0.00618177585431327,	-0.0054655684077809,	-0.00482580410175279,	-0.00425517126789902,	-0.00374695185952026,	-0.00329498607349084,	-0.00289363751082587,	-0.00253775909482396,	-0.00222265993056489,	-0.00194407325694091,	-0.00169812561242549,	-0.00148130730844348,	-0.00129044427947807,	-0.00112267135687822,	-0.000975406993632486,	-0.000846329450052815,	-0.000733354435242717,	-0.000634614186279159,	-0.000548437956074162,	-0.000473333871754026,	-0.000407972117951058,	-0.000351169393493994,	-0.000301874585459822,	-0.000259155601265877,	-0.000222187297296324,	-0.00019024044133701,	-0.00016267164571019,	-0.000138914208336566,	-0.000118469799895754,	-0.000100900936705706,	-8.58241808037005E-05,	-7.29040109018349E-05,	-6.18473103330234E-05,	-5.23984207318817E-05,	-4.43347129494397E-05,	-3.74626295296604E-05,	-3.16141559347925E-05,	-2.66436805577782E-05,	-2.24252063715181E-05,	-1.88498798106082E-05,	-1.58238051400721E-05,	-1.32661151210241E-05,	-1.11072712225503E-05,	-9.28756894331343E-06,	-7.75582598955712E-06,	-6.46823310496068E-06,	-5.38734926114557E-06,	-4.48122469639552E-06,	-3.72263693672212E-06,	-3.08842645145369E-06,	-2.55891998970628E-06,	-2.11743091983458E-06,	-1.74982705723125E-06,	-1.44415752299972E-06,	-1.19033113364039E-06,	-9.79839686631661E-07,	-8.05520285308228E-07,	-6.61351545286602E-07,	-5.42279150237526E-07,	-4.44066783198031E-07,	-3.6316895672282E-07,	-2.96622706560844E-07,	-2.41955504449981E-07,	-1.97107090955884E-07,	-1.60363233604531E-07,	-1.30299683096659E-07,	-1.0573483503648E-07,	-8.56898099242044E-08,	-6.93548434071995E-08,	-5.60610349177208E-08,	-4.52566385254907E-08,	-3.648719751806E-08,	-2.93789260634371E-08,	-2.36248292540009E-08,	-1.89731286156694E-08,	-1.52176253379481E-08,	-1.21896894107027E-08,	-9.75161075580404E-09,	-7.79108938648158E-09,	-6.21667655570292E-09,	-4.95400859279792E-09,	-3.9427003971696E-09,	-3.13378698939575E-09,	-2.487619661706E-09,	-1.97213859772073E-09,	-1.56145675763611E-09,	-1.23470070578345E-09,	-9.75063199502915E-10,	-7.69030025356895E-10,	-6.05749987143255E-10,	-4.76522313622393E-10,	-3.74380227718723E-10,	-2.93753144112661E-10,	-2.30193058410446E-10,	-1.80153260112575E-10,	-1.40809629452325E-10,	-1.09916537547152E-10,	-8.56908214946309E-11,	-6.67185026251758E-11,	-5.18799004022028E-11,	-4.0289602727662E-11,	-3.12484160187545E-11,	-2.42049607052147E-11,	-1.87250208619481E-11,	-1.44671180885289E-11,	-1.11630740392114E-11,	-8.60256525847056E-12,	-6.62086816152344E-12,	-5.0891493474305E-12,	-3.90678372372004E-12,	-2.99528572995842E-12,	-2.29352184767451E-12,	-1.75393932279676E-12,	-1.33959966843097E-12,	-1.02184856148694E-12,	-7.78488234209659E-13,	-5.92346048567396E-13,	-4.50154973423158E-13,	-3.4167925216963E-13,	-2.590325345056E-13,	-1.96146863020266E-13,	-1.48359727908323E-13,	-1.12093394022884E-13,	-8.46062352644517E-14,	-6.38001801504864E-14,	-4.80718180700495E-14,	-3.61974290794863E-14,	-2.72443334808611E-14,	-2.05026322295396E-14,	-1.54327219975927E-14,	-1.16249961553804E-14,	-8.76894581913584E-15,	-6.62950156813477E-15,	-5.02894448766078E-15,	-3.83309484431811E-15,	-2.94078157886584E-15,	-2.2758244914328E-15,	-1.78093808437863E-15,	-1.41310295536566E-15,	-1.14005634730712E-15,	-9.37634955736819E-16,	-7.87765823023325E-16,	-6.76949361588171E-16,	-5.95115548904751E-16,	-5.34762690793086E-16,	-4.90309844892551E-16,	-4.57610572210875E-16,	-4.33588330582687E-16,	-4.15963457289721E-16,	-4.03049015830059E-16,	-3.93598347525974E-16,	-3.8669138982275E-16,	-3.81650019869991E-16,	-3.77975099364223E-16,	-3.75299721998947E-16,	-3.73354541190686E-16,	-3.71942091981665E-16,	-3.70917800063109E-16,	-3.7017595570036E-16,	-3.69639368748613E-16,	-3.69251749112448E-16,	-3.68972102290475E-16,	-3.687706127233E-16,	-3.68625624110026E-16,	-3.68521427405681E-16,	-3.68446642677248E-16,	-3.68393036997636E-16,	-3.68354662055488E-16,	-3.68327225865799E-16,	-3.68307635655955E-16,	-3.68293665743037E-16,	-3.68283716553145E-16,	-3.6827664000867E-16,	-3.68271613176926E-16,	-3.68268046965124E-16,	-3.68265520230226E-16,	-3.68263732293804E-16,	-3.68262468767182E-16,	-3.68261576989217E-16,	-3.68260948396746E-16,	-3.68260505888067E-16,	-3.68260194777565E-16,	-3.68259976329656E-16,	-3.68259823142789E-16,	-3.68259715858614E-16,	-3.68259640819036E-16,	-3.68259588400396E-16,	-3.6825955183066E-16,	-3.682595263507E-16,	-3.68259508620383E-16,	-3.68259496298553E-16,	-3.68259487746409E-16,	-3.68259481818306E-16,	-3.68259477714395E-16,	-3.68259474876992E-16,	-3.68259472917762E-16,	-3.68259471566649E-16,	-3.68259470636101E-16,	-3.6825946999603E-16,	-3.68259469556328E-16,	-3.6825946925466E-16,	-3.68259469047959E-16,	-3.68259468906509E-16,	-3.68259468809838E-16,	-3.68259468743854E-16,	-3.68259468698874E-16,	-3.68259468668252E-16,	-3.68259468647431E-16,	-3.68259468633291E-16,	-3.68259468623704E-16,	-3.68259468617208E-16,	-3.68259468612815E-16,	-3.68259468609847E-16,	-3.68259468607844E-16,	-3.68259468606495E-16,	-3.68259468605587E-16,	-3.68259468604977E-16,	-3.68259468604566E-16,	-3.68259468604293E-16,	-3.68259468604109E-16,	-3.68259468603986E-16,	-3.68259468603903E-16,	-3.68259468603848E-16,	-3.68259468603812E-16,	-3.68259468603788E-16,	-3.68259468603772E-16,	-3.6825946860376E-16,	-3.68259468603754E-16,	-3.68259468603749E-16,	-3.68259468603746E-16,	-3.68259468603744E-16,	-3.68259468603743E-16,	-3.68259468603741E-16,	-3.68259468603741E-16,	-3.68259468603741E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16};


    //double Force_2[1000]={5.91124150106578E-16,	1.44280906975299E-15,	2.66832859567292E-15,	4.42950339721609E-15,	6.95720362633263E-15,	1.05803799964225E-14,	1.57671156037512E-14,	2.31825978202313E-14,	3.3770873992334E-14,	4.88699934507724E-14,	7.03739513734053E-14,	1.00960151283069E-13,	1.44408462895194E-13,	2.06048155783132E-13,	2.93383115844471E-13,	4.1696528138808E-13,	5.91613197533099E-13,	8.38109757283117E-13,	1.18556436748207E-12,	1.67469512160385E-12,	2.3623831268982E-12,	3.32798350778865E-12,	4.68205881220292E-12,	6.57844824345608E-12,	9.23092424479162E-12,	1.29361488195462E-11,	1.81052692489066E-11,	2.53073454531499E-11,	3.53289584026793E-11,	4.92559172858191E-11,	6.85851056206094E-11,	9.53773748722783E-11,	1.32466264977056E-10,	1.837425473303E-10,	2.54541604740929E-10,	3.52170094068342E-10,	4.86620991535375E-10,	6.71543075023316E-10,	9.25553580381976E-10,	1.27401307151672E-09,	1.75142198796354E-09,	2.40465218114194E-09,	3.29729837483368E-09,	4.51552961684349E-09,	6.17594695818391E-09,	8.4361206859413E-09,	1.15087013945863E-08,	1.56802909472646E-08,	2.13366440799427E-08,	2.89962778943405E-08,	3.93552323474972E-08,	5.33465989328368E-08,	7.22195804818425E-08,	9.76443445992193E-08,	1.31850893054803E-07,	1.77812726938176E-07,	2.39489419846236E-07,	3.22146547194081E-07,	4.32777052596444E-07,	5.8065543922564E-07,	7.78065638466864E-07,	1.04125564868058E-06,	1.39168784391276E-06,	1.8576742344041E-06,	2.47651220218108E-06,	3.2972699692088E-06,	4.3844143562086E-06,	5.82252888568589E-06,	7.72244129906874E-06,	1.02291702972303E-05,	1.35322170761678E-05,	1.78788746892341E-05,	2.35914158141719E-05,	3.10892576672158E-05,	4.09175047811854E-05,	5.37836526527788E-05,	7.06047184741849E-05,	9.25676749861559E-05,	0.000121206831895919,	0.00015850277601547,	0.000207008695558934,	0.000270011436887736,	0.000351736548798603,	0.000457608953228707,	0.000594583856569693,	0.000771566223494502,	0.000999941748183324,	0.00129424798788206,	0.00167302142984865,	0.00215986506112053,	0.00278479188644022,	0.00358591326170448,	0.00461155744735637,	0.00592292412761312,	0.00759740562116338,	0.00973273613364896,	0.0124521678848854,	0.0159109187432343,	0.0203041918669714,	0.0258771358858224,	0.0329371968634842,	0.0418694136572463,	0.053155329898573,	0.0673963428955459,	0.0853424873296562,	0.107927865633713,	0.136314194420361,	0.171944245563024,	0.216607331252938,	0.272519425998269,	0.342421048480098,	0.429696658076092,	0.538520072939513,	0.674031309987081,	0.842551306613682,	1.05184223790488,	1.31142262446098,	1.63294817258245,	2.03067134405068,	2.52199506699422,	3.12813882942783,	3.87493870805877,	4.79380675085155,	5.92287963650692,	7.30839177226742,	9.00631407018703,	11.0843066812912,	13.6240421014166,	16.7239644421287,	20.5025614520315,	25.1022382632273,	30.6938960291939,	37.4823348389637,	45.712618784756,	55.6775620950508,	67.7265191145265,	82.2756879317698,	99.820167964287,	120.94804617037,	146.356825153512,	176.872549663045,	213.472036299025,	257.308665040779,	309.742250989405,	372.373580904929,	447.084272189429,	536.082692372735,	641.956765336372,	767.734586877227,	916.953877139732,	1093.74141125632,	1302.90369248379,	1550.03026437387,	1841.61120012818,	2185.17045818741,	2589.41695306468,	3064.41535903982,	3621.77884096221,	4274.88609019313,	5039.12523350939,	5932.16737712635,	6974.27274507297,	8188.63256876907,	9601.75008019802,	11243.8641514629,	13149.4193051887,	15357.5859890965,	17912.835159472,	20865.5713469575,	24272.8284782795,	28199.0327927517,	32716.837215613,	37908.0315238208,	43864.5325555901,	50689.4585639512,	58498.2915876252,	67420.1313998521,	77599.0441873972,	89195.508597499,	102387.961159588,	117374.442330825,	134374.34351973,	153630.25440062,	175409.908634207,	200008.224748259,	227749.437399662,	258989.312530015,	294117.438037155,	333559.579513128,	377780.088345791,	427284.347050108,	482621.234092925,	544385.587711997,	613220.645320271,	689820.432048408,	774932.068835006,	869357.967252603,	973957.87499095,	1089650.73264493,	1217416.29921591,	1358296.50058048,	1513396.45216296,	1683885.10422606,	1870995.45563104,	2076024.27968191,	2300331.30382927,	2545337.78364212,	2812524.41063813,	3103428.49337131,	3419640.3516903,	3762798.86537749,	4134586.1205343,	4536721.10016257,	4970952.36947233,	5439049.71158179,	5942794.67551554,	6483970.00579289,	7064347.93145672,	7685677.30213748,	8349669.56967545,	9057983.62591605,	9812209.52050938,	10613851.096827,	11464307.5993744,	12364854.3222239,	13316622.3848928,	14320577.7395962,	15377499.5317371,	16487957.9536683,	17652291.7499443,	18870585.5502503,	20142647.2236788,	21467985.464767,	22845787.8374042,	24274899.5170939,	25753802.9847852,	27280598.9362941,	28852988.6798907,	30468258.3006638,	32123264.8734931,	33814425.0066082,	35537705.9945551,	37288619.8527122,	39062220.4951317,	40853104.3032961,	42655414.315278,	44462848.2427536,	46268670.4973429,	48065728.3779264,	49846472.5370408,	51602981.807385,	53326992.4291355,	55009931.675495,	56642955.8280648,	58216992.4056869,	59722786.5008407,	61150951.0270571,	62492020.6297256,	63736508.9617484,	64874968.9754219,	65898055.8333695,	66796591.9950328,	67561633.9918474,	68184540.3644833,	68657040.2001004,	68971301.6771055,	69120000,	69096384.0881397,	68894341.3700702,	68508460.0299776,	67934088.0550281,	67167388.442215,	66205389.9409138,	65046032.7327209,	63688208.4832316,	62131794.2410338,	60377679.7070463,	58427787.4520226,	56285085.7210603,	53953593.5306786,	51438377.8357615,	48745542.6195937,	45882209.8394862,	42856492.2421427,	39677458.1459735,	36355088.3709883,	32900225.5796354,	29324516.3729422,	25640346.5644984,	21860770.1291718,	17999432.3929818,	14070488.0943248,	10088515.0039062,	6068423.84049295,	2025365.26128231,	-2025365.26128255,	-6068423.8404932,	-10088515.0039064,	-14070488.094325,	-17999432.3929821,	-21860770.1291721,	-25640346.5644986,	-29324516.3729424,	-32900225.5796356,	-36355088.3709885,	-39677458.1459737,	-42856492.2421429,	-45882209.8394863,	-48745542.6195938,	-51438377.8357617,	-53953593.5306788,	-56285085.7210604,	-58427787.4520227,	-60377679.7070464,	-62131794.241034,	-63688208.4832317,	-65046032.732721,	-66205389.9409139,	-67167388.442215,	-67934088.0550282,	-68508460.0299777,	-68894341.3700703,	-69096384.0881397,	-69120000,	-68971301.6771055,	-68657040.2001004,	-68184540.3644832,	-67561633.9918474,	-66796591.9950328,	-65898055.8333695,	-64874968.9754219,	-63736508.9617483,	-62492020.6297255,	-61150951.0270571,	-59722786.5008406,	-58216992.4056868,	-56642955.8280647,	-55009931.6754949,	-53326992.4291354,	-51602981.8073849,	-49846472.5370408,	-48065728.3779263,	-46268670.4973428,	-44462848.2427535,	-42655414.3152779,	-40853104.303296,	-39062220.4951316,	-37288619.8527121,	-35537705.994555,	-33814425.006608,	-32123264.873493,	-30468258.3006637,	-28852988.6798906,	-27280598.936294,	-25753802.9847851,	-24274899.5170938,	-22845787.8374042,	-21467985.4647669,	-20142647.2236788,	-18870585.5502502,	-17652291.7499443,	-16487957.9536682,	-15377499.5317371,	-14320577.7395961,	-13316622.3848927,	-12364854.3222238,	-11464307.5993743,	-10613851.0968269,	-9812209.52050933,	-9057983.625916,	-8349669.5696754,	-7685677.30213743,	-7064347.93145668,	-6483970.00579285,	-5942794.6755155,	-5439049.71158176,	-4970952.3694723,	-4536721.10016254,	-4134586.12053428,	-3762798.86537746,	-3419640.35169027,	-3103428.49337129,	-2812524.41063811,	-2545337.78364209,	-2300331.30382925,	-2076024.27968189,	-1870995.45563102,	-1683885.10422604,	-1513396.45216294,	-1358296.50058046,	-1217416.29921589,	-1089650.73264491,	-973957.874990935,	-869357.967252589,	-774932.068834993,	-689820.432048395,	-613220.645320258,	-544385.587711985,	-482621.234092913,	-427284.347050096,	-377780.088345779,	-333559.579513116,	-294117.438037143,	-258989.312530003,	-227749.43739965,	-200008.224748247,	-175409.908634195,	-153630.254400608,	-134374.343519719,	-117374.442330814,	-102387.961159576,	-89195.5085974878,	-77599.044187386,	-67420.1313998411,	-58498.2915876142,	-50689.4585639402,	-43864.5325555791,	-37908.0315238098,	-32716.837215602,	-28199.0327927408,	-24272.8284782686,	-20865.5713469466,	-17912.8351594611,	-15357.5859890856,	-13149.4193051778,	-11243.864151452,	-9601.75008018721,	-8188.63256875825,	-6974.27274506216,	-5932.16737711554,	-5039.1252334986,	-4274.88609018234,	-3621.77884095142,	-3064.41535902904,	-2589.4169530539,	-2185.17045817663,	-1841.6112001174,	-1550.0302643631,	-1302.90369247302,	-1093.74141124555,	-916.95387712896,	-767.734586866455,	-641.956765325601,	-536.082692361965,	-447.084272178659,	-372.373580894159,	-309.742250978635,	-257.308665030011,	-213.472036288257,	-176.872549652276,	-146.356825142744,	-120.948046159602,	-99.8201679535189,	-82.2756879210017,	-67.7265191037587,	-55.677562084283,	-45.7126187739882,	-37.482334828196,	-30.6938960184262,	-25.1022382524597,	-20.5025614412639,	-16.7239644313611,	-13.624042090649,	-11.0843066705236,	-9.00631405941953,	-7.30839176149992,	-5.92287962573943,	-4.79380674008408,	-3.87493869729131,	-3.12813881866037,	-2.52199505622676,	-2.03067133328322,	-1.63294816181499,	-1.31142261369353,	-1.05184222713742,	-0.842551295846233,	-0.674031299219633,	-0.538520062172066,	-0.429696647308645,	-0.342421037712652,	-0.272519415230823,	-0.216607320485494,	-0.171944234795579,	-0.136314183652917,	-0.107927854866269,	-0.0853424765622123,	-0.0673963321281021,	-0.0531553191311291,	-0.0418694028898026,	-0.0329371860960406,	-0.0258771251183789,	-0.0203041810995279,	-0.0159109079757908,	-0.012452157117442,	-0.0097327253662056,	-0.00759739485372003,	-0.0059229133601698,	-0.00461154667991306,	-0.00358590249426117,	-0.00278478111899691,	-0.00215985429367723,	-0.00167301066240536,	-0.00129423722043877,	-0.000999930980740036,	-0.000771555456051218,	-0.000594573089126411,	-0.000457598185785426,	-0.000351725781355323,	-0.000270000669444456,	-0.000206997928115654,	-0.000158492008572191,	-0.00012119606445264,	-9.25569075428769E-05,	-7.05939510309063E-05,	-5.37728852095004E-05,	-4.09067373379071E-05,	-3.10784902239375E-05,	-2.35806483708937E-05,	-0.000017868107245956,	-1.35214496328897E-05,	-1.02184028539523E-05,	-7.71167385579073E-06,	-5.8117614424079E-06,	-4.37364691293063E-06,	-3.28650252593082E-06,	-2.46574475890311E-06,	-1.84690679112614E-06,	-1.3809204006348E-06,	-1.03048820540263E-06,	-7.67298195188918E-07,	-5.69887995947696E-07,	-4.22009609318502E-07,	-3.1137910391614E-07,	-2.28721976568295E-07,	-1.67045283660235E-07,	-1.21083449776863E-07,	-8.68769013212792E-08,	-6.14521372039023E-08,	-4.25791556548968E-08,	-2.85877890695573E-08,	-1.82288346164007E-08,	-1.05692008020029E-08,	-4.91284766932494E-09,	-7.41258116646708E-10,	2.33132259199832E-09,	4.59149631975572E-09,	6.2519136610961E-09,	7.47014490310589E-09,	8.36279109679763E-09,	9.01602128997603E-09,	9.49343020642283E-09,	9.84188969755757E-09,	1.00959002029162E-08,	1.02808222864042E-08,	1.04152731838712E-08,	1.05129016731986E-08,	1.05837007306092E-08,	1.06349770129625E-08,	1.06720659030673E-08,	1.06988581723189E-08,	1.07181873606537E-08,	1.07321143195369E-08,	1.07421359324864E-08,	1.07493380086906E-08,	1.075450712912E-08,	1.07582123536948E-08,	1.07608648296961E-08,	1.07627612191273E-08,	1.07641152944318E-08,	1.07650808948126E-08,	1.07657685828179E-08,	1.07662577135721E-08,	1.07666051681823E-08,	1.0766851664742E-08,	1.07670263126582E-08,	1.07671498948237E-08,	1.07672372297838E-08,	1.07672988694766E-08,	1.07673423177883E-08,	1.07673729039882E-08,	1.07673944079461E-08,	1.07674095070656E-08,	1.07674200953417E-08,	1.07674275108239E-08,	1.07674326975595E-08,	1.07674363207359E-08,	1.07674388484361E-08,	1.07674406096109E-08,	1.07674418351305E-08,	1.07674426868154E-08,	1.07674432779395E-08,	1.07674436876901E-08,	1.07674439713523E-08,	1.07674441674737E-08,	1.07674443028957E-08,	1.07674443962846E-08,	1.07674444606041E-08,	1.07674445048458E-08,	1.07674445352379E-08,	1.07674445560893E-08,	1.07674445703765E-08,	1.07674445801535E-08,	1.07674445868354E-08,	1.07674445913962E-08,	1.07674445945052E-08,	1.07674445966218E-08,	1.0767444598061E-08,	1.07674445990382E-08,	1.0767444599701E-08,	1.07674446001499E-08,	1.07674446004536E-08,	1.07674446006587E-08,	1.07674446007972E-08,	1.07674446008904E-08,	1.07674446009532E-08,	1.07674446009954E-08,	1.07674446010237E-08,	1.07674446010427E-08,	1.07674446010554E-08,	1.07674446010639E-08,	1.07674446010696E-08,	1.07674446010734E-08,	1.07674446010759E-08,	1.07674446010776E-08,	1.07674446010787E-08,	1.07674446010794E-08,	1.07674446010799E-08,	1.07674446010802E-08,	1.07674446010804E-08,	1.07674446010806E-08,	1.07674446010807E-08,	1.07674446010807E-08,	1.07674446010808E-08,	1.07674446010808E-08,	1.07674446010808E-08,	1.07674446010808E-08,	1.07674446010808E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08,	1.07674446010809E-08};
/*
	double sum1=0,sum2=0;				 // comment this line for ricker pulse
	for(int i=0; i < timestep; i++)  // comment this line for ricker pulse
	{                  // comment this line for ricker pulse

		// first node
		double t=i*ri_dt;
		double alfa1 = ( PI * fc ) * ( PI * fc ) * ( t - zp / Vs - Ts) * ( t - zp / Vs - Ts);
		double alfa2 = ( PI * fc ) * ( PI * fc ) * ( t + zp / Vs - Ts) * ( t + zp / Vs - Ts);

		double uo1 = ( 2.0 * alfa1 - 1.0 ) * exp(-alfa1);
		double uo2 = ( 2.0 * alfa2 - 1.0 ) * exp(-alfa2);

		double force_11 = (uo1+uo2);

		sum1 = sum1 + force_11;
	}

	//double force_1 = sum1 * ri_dt * force_coef;
	 * */

	double force_1 = Force_1[timestep];
	double force_2 = -1*Force_1[timestep+5];

	double force_coefficient =    3538944;

	double force_1x = force_1 * force_coefficient;
	double force_1z = force_1 * force_coefficient/2;

	double force_2x = force_2 * force_coefficient;
	double force_2z = force_2 * force_coefficient/2;


	// second node
/*
	for(int i=0;i<timestep;i++)  // comment this line for ricker pulse
	{                  // comment this line for ricker pulse
		double t=(i+5)*ri_dt;
		double alfa1 = ( PI * fc ) * ( PI * fc ) * ( t - zp / Vs - Ts) * ( t - zp / Vs - Ts);
		double alfa2 = ( PI * fc ) * ( PI * fc ) * ( t + zp / Vs - Ts) * ( t + zp / Vs - Ts);

		double uo1 = ( 2.0 * alfa1 - 1.0 ) * exp(-alfa1);
		double uo2 = ( 2.0 * alfa2 - 1.0 ) * exp(-alfa2);

		double force_22 = -(uo1+uo2)*(1);
		sum2 = sum2 + force_22;
	}

	double force_2 = sum2*ri_dt* force_coef;
	*/
	/*
         	//first node
         	double mean=3;
         	double sigma=0.2;
         	double t=(timestep)*0.001;
         	double force_1=(1/sqrt(2*PI))/sigma*exp(-(t-mean)*(t-mean)/2/sigma/sigma);

            //second node
         	t=(timestep+64)*0.001;
         	double force_2=(1/sqrt(2*PI))/sigma*exp(-(t-mean)*(t-mean)/2/sigma/sigma);

	 */

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

	//          printf("Number of bottom nodes are : %i - %i",k1,k2);




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

