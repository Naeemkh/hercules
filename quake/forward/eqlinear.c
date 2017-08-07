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


    if ( z_m == depth) {
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

	  double thGDtable[11][3] = {{ 0.0001,	0.999997975588779, 0.24},
				                       { 0.0003,	0.999993221555393, 0.42},
				                       { 0.0010,	0.999974514770621, 0.80},
				                       { 0.0030,	0.999914671141506, 1.40},
				                       { 0.0100,	0.999679254704751, 2.80},
			                           { 0.0300,	0.998926834567793, 5.10},
				                       { 0.1000,	0.995977010223477, 9.80},
				                       { 0.3000,	0.986655614282645, 15.50},
				                       { 1.0000,	0.951609682956969, 21.00},
				                       { 3.0000,	0.85450513321844, 25.00},
				                       { 10.0000,	0.609690543239149, 28.00},
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

		    x_m = (myMesh->ticksize)*(double)myMesh->nodeTable[nindex].x;
		    y_m = (myMesh->ticksize)*(double)myMesh->nodeTable[nindex].y;
		    z_m = (myMesh->ticksize)*(double)myMesh->nodeTable[nindex].z;



// the following part of code is written to extract the strain.Uncomment it. Alos double check the isThisElementsAtTheBottom_eq function





			//For each element you need to define another strain tensor, to only keep the maximum value of strain.

            maxstrains->qp[i].xx = MAX(fabs(maxstrains->qp[i].xx),fabs(tstrains->qp[i].xx));
            maxstrains->qp[i].yy = MAX(fabs(maxstrains->qp[i].yy),fabs(tstrains->qp[i].yy));
            maxstrains->qp[i].zz = MAX(fabs(maxstrains->qp[i].zz),fabs(tstrains->qp[i].zz));
            maxstrains->qp[i].xy = MAX(fabs(maxstrains->qp[i].xy),fabs(tstrains->qp[i].xy));
            maxstrains->qp[i].yz = MAX(fabs(maxstrains->qp[i].yz),fabs(tstrains->qp[i].yz));
            maxstrains->qp[i].xz = MAX(fabs(maxstrains->qp[i].xz),fabs(tstrains->qp[i].xz));


//            printf("This is xy max strain: %.20f (node %i) \n", maxstrains->qp[i].xy, i);


            // the following part of code is written to extract the max strain.Uncomment it
//			if (  x_m == 4096 && y_m == 4096) {
//				printf("STST el_eindex = %i , node =%i , depth = %f , timestep = %i , maxxz = %.20f   \n",el_eindex,i,z_m,step,maxstrains->qp[i].xz);
//			 }

//			if (  x_m == 4096 && y_m == 4096) {
//				printf("STST el_eindex = %i,node =%i,depth = %f, timestep = %i ,xz = %.20f, maxxz = %.20f   \n",el_eindex,i,z_m,step,tstrains->qp[i].xz,maxstrains->qp[i].xz);
//			 }


			 //FILE *fp_st = hu_fopen( "element_strain.txt", "a" );
            // step % 10 == 0 &&

      //      if (x_m == 4096 && y_m == 4096 && ( z_m == 32 || z_m == 64 || z_m == 96 || z_m == 128 || z_m == 160 || z_m == 192 || z_m == 224 || z_m == 256 || z_m == 288 || z_m == 320 || z_m == 352 || z_m == 384 || z_m == 416 || z_m == 448 || z_m == 480 || z_m == 512)){
			    //printf("\n STST el_eindex = %i node =%i depth = %f  timestep = %i xx = %.20f  yy = %.20f  zz = %.20f  xy = %.20f  yz = %.20f  xz = %.20f \n",el_eindex,i,z_m,step,tstrains->qp[i].xx,tstrains->qp[i].yy,tstrains->qp[i].zz,tstrains->qp[i].xy,tstrains->qp[i].yz,tstrains->qp[i].xz);
	//		    printf("\n STST el_eindex = %i node = %i depth = %f  timestep = %i xx = %.10f  yy = %.10f  zz = %.10f  xy = %.10f  yz = %.10f  xz = %.10f \n",el_eindex,i,z_m,step,tstrains->qp[i].xx,tstrains->qp[i].yy,tstrains->qp[i].zz,tstrains->qp[i].xy,tstrains->qp[i].yz,tstrains->qp[i].xz);

      //      }


			 //hu_fclosep( &fp_st );

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
							   //double      *theQTABLE,
							   //int         QTable_Size,
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
		double         updated_Q;
        double         x_m,y_m,z_m;
        double         depth=512;

		/* Capture data from the element and mesh */
		eindex = myEqlinElementsMapping[el_eindex];

		elemp = &myMesh->elemTable[eindex];
		edata = (edata_t *)elemp->data;
		h     = edata->edgesize;
		ep    = &mySolver->eTable[eindex];

	    /* Capture the element's last node at the bottom */
	    elemp  = &myMesh->elemTable[eindex];
	    nindex = elemp->lnid[7];


	    x_m = (myMesh->ticksize)*(double)myMesh->nodeTable[nindex].x;
	    y_m = (myMesh->ticksize)*(double)myMesh->nodeTable[nindex].y;
		z_m = (myMesh->ticksize)*(double)myMesh->nodeTable[nindex].z;



		if (isThisElementsAtTheBottom_eq(myMesh, el_eindex, depth)==YES){
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
		// --------------------------------------
        // Here I need to add a function to pick the best strain combination
		double myteststrain = 0;


		// block of code for single strain
		// ---------------
		for (i=0; i < 8; i++){
		myteststrain = strain_mat[i][5] + myteststrain;
		}

		// ---------------


		// J2 strain (copied the formula from nonlinear.c)
		// ---------------
		//for (i=0; i < 8; i++){
		//myteststrain = ((strain_mat[i][1]*strain_mat[i][1]+strain_mat[i][2]*strain_mat[i][2]+strain_mat[i][3]*strain_mat[i][3])*0.5 +
		//		        (strain_mat[i][1]*strain_mat[i][2]+strain_mat[i][2]*strain_mat[i][3]+strain_mat[i][1]*strain_mat[i][3])) + myteststrain;
		//}
		// --------------------------------------

		myteststrain = 100* (myteststrain/8)*1; // I add a factor of 2 just to make the effective strain to be consistent with 1D. Later on we need to define a 3D curves.



//		printf("Here is the max strain of element %i : %.10f \n", el_eindex, myteststrain);

		 GD_t GD = search_GD_table(myteststrain);

		//printf("Results : strain = %f , G = %f, D = %f \n", myteststrain, GD.g, GD.d);


		 updated_mu = original_mu * GD.g;



			if (x_m == 4096 && y_m == 4096) {
				printf("STST el_eindex = %i , depth = %f , updated_mu_f = %f , myteststrain = %.20f   \n",el_eindex,z_m,GD.g,myteststrain);
			 }




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
         // Q set to be high value to run as without damping. In the real simulation I need to set the real values (like 100vs)
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




	    updated_Q = 1 / (2*total_requested_damping);


//	    /* regenerate the Q table */
//
//	    int i,k,j,l;
//	    double myQTABLE[26][6];
//
//	    int array_size=(QTable_Size*6)+1;
//
//	    // regenertate the Q_table
//	    k=0;
//	    j=0;
//		for(l = 1; l < array_size; l++)
//		{
//			i = l-1;
//			myQTABLE[k][j] = theQTABLE[i];
//	        if (l % 6 ==0){
//	        	k=k+1;
//	        }
//	        j=j+1;
//	        if (j==6){
//	        	j=0;
//	        }
//		}
//
//	    double index_Qs = Search_Quality_Table(new_Q,&(myQTABLE[0][0]), QTable_Size);

		update_Q_params(edata,updated_Q,0,theFreq);
	    control_correction_factor(edata,theFreq_Vel,theFreq);

	} /* for all nonlinear elements */

}


void    compute_addforce_bottom(int32_t timestep, mesh_t *myMesh, mysolver_t *mySolver)
{

    // dt = 0.01
    double Force_1[1000]={1.35111255205411E-24,	3.31034338608617E-24,	6.14774662326503E-24,	1.02516565323711E-23,	1.61797605160242E-23,	2.47319038914087E-23,	3.70537343261079E-23,	5.47840515413847E-23,	8.02640110544514E-23,	1.16833748165512E-22,	1.69252387779081E-22,	2.44291965431406E-22,	3.51576141017056E-22,	5.04763006932892E-22,	7.23210915764468E-22,	1.03432141811737E-21,	1.47683009736597E-21,	2.10542256860622E-21,	2.99720053284144E-21,	4.26072715409923E-21,	6.0486635768945E-21,	8.5753984748259E-21,	1.21416102767648E-20,	1.71684420215414E-20,	2.42449864966829E-20,	3.41941763880942E-20,	4.81640893060695E-20,	6.77542991500259E-20,	9.5190488838295E-20,	1.33565430987247E-19,	1.87171110598859E-19,	2.61955839032891E-19,	3.6615254337669E-19,	5.11141156651131E-19,	7.12630723826186E-19,	9.92277545799426E-19,	1.37989718196362E-18,	1.91648413371098E-18,	2.65832849646046E-18,	3.682620415017E-18,	5.09506962403629E-18,	7.04025043229777E-18,	9.7156277975747E-18,	1.33905483033423E-17,	1.8431918256102E-17,	2.53388759593255E-17,	3.47895442634112E-17,	4.77039857230719E-17,	6.53288590160379E-17,	8.93511006442267E-17,	1.22050373325902E-16,	1.66503219226438E-16,	2.26856077338103E-16,	3.08689890021523E-16,	4.19506351456677E-16,	5.69375484170172E-16,	7.71796875740479E-16,	1.04484348379901E-15,	1.41267861281199E-15,	1.90756501986615E-15,	2.5725221072992E-15,	3.46483537275146E-15,	4.66068501609413E-15,	6.26124209656813E-15,	8.40068634756923E-15,	1.12567366838138E-14,	1.50644625260262E-14,	2.01343727579731E-14,	2.68760740092946E-14,	3.582916960792E-14,	4.77035585984831E-14,	6.34319206789198E-14,	8.42379757928853E-14,	1.11725134551318E-13,	1.47991468436758E-13,	1.95778603548699E-13,	2.58664275034035E-13,	3.41310992698064E-13,	4.49786713951593E-13,	5.9197778909583E-13,	7.78119974738096E-13,	1.02148030201538E-12,	1.33923140895941E-12,	1.75357106332519E-12,	2.29315358820295E-12,	2.99491747048686E-12,	3.90641546424848E-12,	5.08878108795893E-12,	6.62049990205187E-12,	8.60219699899902E-12,	1.11627057797399E-11,	1.44667498290573E-11,	1.87246526024766E-11,	2.42045924457432E-11,	3.12480477592831E-11,	4.02892344681907E-11,	5.18795321427313E-11,	6.67181343657044E-11,	8.56904532351597E-11,	1.09916169287681E-10,	1.40809261192854E-10,	1.80152891853103E-10,	2.30192690150975E-10,	2.9375277585319E-10,	3.74379859459252E-10,	4.76521945362924E-10,	6.05749618883788E-10,	7.69029657097426E-10,	9.75062831243446E-10,	1.23470033752398E-09,	1.56145638937665E-09,	1.97213822946126E-09,	2.48761929344654E-09,	3.13378662113627E-09,	3.94270002891016E-09,	4.95400822453848E-09,	6.21667618744348E-09,	7.79108901822214E-09,	9.75161038754459E-09,	1.21896890424433E-08,	1.52176249696887E-08,	1.897312824741E-08,	2.36248288857415E-08,	2.93789256951777E-08,	3.64871971498006E-08,	4.52566381572313E-08,	5.60610345494614E-08,	6.93548430389403E-08,	8.5689809555945E-08,	1.05734834668221E-07,	1.302996827284E-07,	1.60363233236273E-07,	1.97107090587624E-07,	2.41955504081723E-07,	2.96622706192584E-07,	3.63168956354562E-07,	4.44066782829773E-07,	5.42279149869267E-07,	6.61351544918343E-07,	8.05520284939971E-07,	9.79839686263406E-07,	1.19033113327213E-06,	1.44415752263147E-06,	1.74982705686299E-06,	2.11743091946632E-06,	2.55891998933804E-06,	3.08842645108543E-06,	3.72263693635386E-06,	4.48122469602727E-06,	5.38734926077733E-06,	6.46823310459243E-06,	7.75582598918887E-06,	9.2875689429452E-06,	1.11072712221821E-05,	1.32661151206558E-05,	1.58238051397037E-05,	0.00001884987981024,	2.24252063711497E-05,	0.00002664368055741,	3.16141559344243E-05,	3.74626295292922E-05,	4.43347129490715E-05,	5.23984207315135E-05,	6.18473103326552E-05,	7.29040109014669E-05,	8.58241808033325E-05,	0.000100900936705338,	0.000118469799895386,	0.000138914208336198,	0.000162671645709822,	0.000190240441336642,	0.000222187297295956,	0.000259155601265509,	0.000301874585459455,	0.000351169393493627,	0.00040797211795069,	0.000473333871753659,	0.000548437956073796,	0.000634614186278792,	0.00073335443524235,	0.000846329450052448,	0.000975406993632119,	0.00112267135687785,	0.0012904442794777,	0.00148130730844311,	0.00169812561242513,	0.00194407325694054,	0.00222265993056453,	0.0025377590948236,	0.00289363751082551,	0.00329498607349048,	0.00374695185951989,	0.00425517126789866,	0.00482580410175243,	0.00546556840778054,	0.00618177585431291,	0.0069823673913907,	0.00787594889629625,	0.00887182646586039,	0.00998004097292257,	0.0112114014588398,	0.0125775168873351,	0.0140908257377163,	0.0157646228681269,	0.0176130830326376,	0.0196512803903425,	0.0218952033009686,	0.024361763660678,	0.027068799994662,	0.0300350734907684,	0.0332802561317892,	0.0368249100642666,	0.0406904573298339,	0.0448991390823394,	0.0494739634214452,	0.054438640992151,	0.0598175075308782,	0.0656354325833704,	0.0719177136787083,	0.0786899553180234,	0.0859779322267872,	0.0938074364264016,	0.10220410780464,	0.111193248005461,	0.120799617616832,	0.131047216810128,	0.141959049775888,	0.153556873507335,	0.165860931703912,	0.178889674800692,	0.192659467374028,	0.207184284427028,	0.222475398317956,	0.23854105835747,	0.255386165363778,	0.273011943724686,	0.291415613768501,	0.310590067487947,	0.330523550888371,	0.351199356439436,	0.372595529293766,	0.394684591092216,	0.417433285299343,	0.440802348099881,	0.464746308933637,	0.489213324748381,	0.514145052004557,	0.539476560369099,	0.565136291885586,	0.591046069202778,	0.617121156181683,	0.643270373882431,	0.669396274556464,	0.695395375838056,	0.721158456843762,	0.746570917352214,	0.771513200653001,	0.795861280027415,	0.819487208160658,	0.842259728091535,	0.86404494358855,	0.884707046109044,	0.90410909475876,	0.922113844932701,	0.938584620593488,	0.95338622444111,	0.966385879557766,	0.977454195483903,	0.986466151106528,	0.993302086228303,	0.997848693245159,	1,	0.999658334608501,	0.996735262877173,	0.99115248885963,	0.982842709129457,	0.971750411490379,	0.957832609098869,	0.941059501341448,	0.921415053287494,	0.898897486126068,	0.873519671687591,	0.845309424942457,	0.814309689251451,	0.780578610108199,	0.744189494151644,	0.70523065132514,	0.663805119205529,	0.620030269706927,	0.57403729956559,	0.525970607219161,	0.475987059890558,	0.42425515585854,	0.370954088028045,	0.316272715989177,	0.260408454759576,	0.203566089327616,	0.145956524940774,	0.0877954838034294,	0.0293021594514236,	-0.0293021594514245,	-0.0877954838034304,	-0.145956524940774,	-0.203566089327617,	-0.260408454759578,	-0.316272715989178,	-0.370954088028045,	-0.424255155858542,	-0.47598705989056,	-0.525970607219161,	-0.57403729956559,	-0.620030269706927,	-0.663805119205531,	-0.705230651325142,	-0.744189494151644,	-0.780578610108199,	-0.814309689251451,	-0.845309424942458,	-0.873519671687591,	-0.898897486126069,	-0.921415053287494,	-0.941059501341448,	-0.957832609098869,	-0.971750411490379,	-0.982842709129457,	-0.99115248885963,	-0.996735262877174,	-0.999658334608503,	-1,	-0.997848693245161,	-0.993302086228305,	-0.986466151106529,	-0.977454195483903,	-0.966385879557766,	-0.95338622444111,	-0.938584620593488,	-0.922113844932701,	-0.904109094758762,	-0.884707046109044,	-0.86404494358855,	-0.842259728091535,	-0.81948720816066,	-0.795861280027415,	-0.771513200653001,	-0.746570917352214,	-0.721158456843762,	-0.695395375838056,	-0.669396274556466,	-0.643270373882431,	-0.617121156181683,	-0.591046069202778,	-0.565136291885586,	-0.539476560369099,	-0.514145052004557,	-0.489213324748381,	-0.464746308933637,	-0.440802348099881,	-0.417433285299343,	-0.394684591092216,	-0.372595529293766,	-0.351199356439436,	-0.330523550888371,	-0.310590067487949,	-0.291415613768501,	-0.273011943724686,	-0.255386165363778,	-0.23854105835747,	-0.222475398317956,	-0.207184284427028,	-0.192659467374028,	-0.178889674800692,	-0.165860931703912,	-0.153556873507335,	-0.141959049775888,	-0.131047216810128,	-0.120799617616832,	-0.111193248005461,	-0.10220410780464,	-0.0938074364264019,	-0.0859779322267875,	-0.0786899553180237,	-0.0719177136787086,	-0.0656354325833708,	-0.0598175075308785,	-0.0544386409921513,	-0.0494739634214455,	-0.0448991390823399,	-0.0406904573298344,	-0.0368249100642671,	-0.0332802561317895,	-0.0300350734907687,	-0.0270687999946623,	-0.0243617636606782,	-0.021895203300969,	-0.0196512803903429,	-0.0176130830326379,	-0.0157646228681272,	-0.0140908257377166,	-0.0125775168873355,	-0.0112114014588401,	-0.00998004097292293,	-0.00887182646586075,	-0.00787594889629661,	-0.00698236739139106,	-0.00618177585431327,	-0.0054655684077809,	-0.00482580410175279,	-0.00425517126789902,	-0.00374695185952026,	-0.00329498607349084,	-0.00289363751082587,	-0.00253775909482396,	-0.00222265993056489,	-0.00194407325694091,	-0.00169812561242549,	-0.00148130730844348,	-0.00129044427947807,	-0.00112267135687822,	-0.000975406993632486,	-0.000846329450052815,	-0.000733354435242717,	-0.000634614186279159,	-0.000548437956074162,	-0.000473333871754026,	-0.000407972117951058,	-0.000351169393493994,	-0.000301874585459822,	-0.000259155601265877,	-0.000222187297296324,	-0.00019024044133701,	-0.00016267164571019,	-0.000138914208336566,	-0.000118469799895754,	-0.000100900936705706,	-8.58241808037005E-05,	-7.29040109018349E-05,	-6.18473103330234E-05,	-5.23984207318817E-05,	-4.43347129494397E-05,	-3.74626295296604E-05,	-3.16141559347925E-05,	-2.66436805577782E-05,	-2.24252063715181E-05,	-1.88498798106082E-05,	-1.58238051400721E-05,	-1.32661151210241E-05,	-1.11072712225503E-05,	-9.28756894331343E-06,	-7.75582598955712E-06,	-6.46823310496068E-06,	-5.38734926114557E-06,	-4.48122469639552E-06,	-3.72263693672212E-06,	-3.08842645145369E-06,	-2.55891998970628E-06,	-2.11743091983458E-06,	-1.74982705723125E-06,	-1.44415752299972E-06,	-1.19033113364039E-06,	-9.79839686631661E-07,	-8.05520285308228E-07,	-6.61351545286602E-07,	-5.42279150237526E-07,	-4.44066783198031E-07,	-3.6316895672282E-07,	-2.96622706560844E-07,	-2.41955504449981E-07,	-1.97107090955884E-07,	-1.60363233604531E-07,	-1.30299683096659E-07,	-1.0573483503648E-07,	-8.56898099242044E-08,	-6.93548434071995E-08,	-5.60610349177208E-08,	-4.52566385254907E-08,	-3.648719751806E-08,	-2.93789260634371E-08,	-2.36248292540009E-08,	-1.89731286156694E-08,	-1.52176253379481E-08,	-1.21896894107027E-08,	-9.75161075580404E-09,	-7.79108938648158E-09,	-6.21667655570292E-09,	-4.95400859279792E-09,	-3.9427003971696E-09,	-3.13378698939575E-09,	-2.487619661706E-09,	-1.97213859772073E-09,	-1.56145675763611E-09,	-1.23470070578345E-09,	-9.75063199502915E-10,	-7.69030025356895E-10,	-6.05749987143255E-10,	-4.76522313622393E-10,	-3.74380227718723E-10,	-2.93753144112661E-10,	-2.30193058410446E-10,	-1.80153260112575E-10,	-1.40809629452325E-10,	-1.09916537547152E-10,	-8.56908214946309E-11,	-6.67185026251758E-11,	-5.18799004022028E-11,	-4.0289602727662E-11,	-3.12484160187545E-11,	-2.42049607052147E-11,	-1.87250208619481E-11,	-1.44671180885289E-11,	-1.11630740392114E-11,	-8.60256525847056E-12,	-6.62086816152344E-12,	-5.0891493474305E-12,	-3.90678372372004E-12,	-2.99528572995842E-12,	-2.29352184767451E-12,	-1.75393932279676E-12,	-1.33959966843097E-12,	-1.02184856148694E-12,	-7.78488234209659E-13,	-5.92346048567396E-13,	-4.50154973423158E-13,	-3.4167925216963E-13,	-2.590325345056E-13,	-1.96146863020266E-13,	-1.48359727908323E-13,	-1.12093394022884E-13,	-8.46062352644517E-14,	-6.38001801504864E-14,	-4.80718180700495E-14,	-3.61974290794863E-14,	-2.72443334808611E-14,	-2.05026322295396E-14,	-1.54327219975927E-14,	-1.16249961553804E-14,	-8.76894581913584E-15,	-6.62950156813477E-15,	-5.02894448766078E-15,	-3.83309484431811E-15,	-2.94078157886584E-15,	-2.2758244914328E-15,	-1.78093808437863E-15,	-1.41310295536566E-15,	-1.14005634730712E-15,	-9.37634955736819E-16,	-7.87765823023325E-16,	-6.76949361588171E-16,	-5.95115548904751E-16,	-5.34762690793086E-16,	-4.90309844892551E-16,	-4.57610572210875E-16,	-4.33588330582687E-16,	-4.15963457289721E-16,	-4.03049015830059E-16,	-3.93598347525974E-16,	-3.8669138982275E-16,	-3.81650019869991E-16,	-3.77975099364223E-16,	-3.75299721998947E-16,	-3.73354541190686E-16,	-3.71942091981665E-16,	-3.70917800063109E-16,	-3.7017595570036E-16,	-3.69639368748613E-16,	-3.69251749112448E-16,	-3.68972102290475E-16,	-3.687706127233E-16,	-3.68625624110026E-16,	-3.68521427405681E-16,	-3.68446642677248E-16,	-3.68393036997636E-16,	-3.68354662055488E-16,	-3.68327225865799E-16,	-3.68307635655955E-16,	-3.68293665743037E-16,	-3.68283716553145E-16,	-3.6827664000867E-16,	-3.68271613176926E-16,	-3.68268046965124E-16,	-3.68265520230226E-16,	-3.68263732293804E-16,	-3.68262468767182E-16,	-3.68261576989217E-16,	-3.68260948396746E-16,	-3.68260505888067E-16,	-3.68260194777565E-16,	-3.68259976329656E-16,	-3.68259823142789E-16,	-3.68259715858614E-16,	-3.68259640819036E-16,	-3.68259588400396E-16,	-3.6825955183066E-16,	-3.682595263507E-16,	-3.68259508620383E-16,	-3.68259496298553E-16,	-3.68259487746409E-16,	-3.68259481818306E-16,	-3.68259477714395E-16,	-3.68259474876992E-16,	-3.68259472917762E-16,	-3.68259471566649E-16,	-3.68259470636101E-16,	-3.6825946999603E-16,	-3.68259469556328E-16,	-3.6825946925466E-16,	-3.68259469047959E-16,	-3.68259468906509E-16,	-3.68259468809838E-16,	-3.68259468743854E-16,	-3.68259468698874E-16,	-3.68259468668252E-16,	-3.68259468647431E-16,	-3.68259468633291E-16,	-3.68259468623704E-16,	-3.68259468617208E-16,	-3.68259468612815E-16,	-3.68259468609847E-16,	-3.68259468607844E-16,	-3.68259468606495E-16,	-3.68259468605587E-16,	-3.68259468604977E-16,	-3.68259468604566E-16,	-3.68259468604293E-16,	-3.68259468604109E-16,	-3.68259468603986E-16,	-3.68259468603903E-16,	-3.68259468603848E-16,	-3.68259468603812E-16,	-3.68259468603788E-16,	-3.68259468603772E-16,	-3.6825946860376E-16,	-3.68259468603754E-16,	-3.68259468603749E-16,	-3.68259468603746E-16,	-3.68259468603744E-16,	-3.68259468603743E-16,	-3.68259468603741E-16,	-3.68259468603741E-16,	-3.68259468603741E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16,	-3.6825946860374E-16};
    // dt = 0.001
    //double Force_1[8000] ={1.35097001840288E-25,	2.75317100130399E-25,	4.20852699153141E-25,	5.71903357559868E-25,	7.28676051684058E-25,	8.91385448444084E-25,	1.06025418818064E-24,	1.23551317778651E-24,	1.41740189449914E-24,	1.60616870073975E-24,	1.80207117039611E-24,	2.00537642696058E-24,	2.21636149394938E-24,	2.43531365804426E-24,	2.66253084541348E-24,	2.89832201168516E-24,	3.14300754606287E-24,	3.39691969009068E-24,	3.66040297159291E-24,	3.93381465433255E-24,	4.21752520395128E-24,	4.51191877077426E-24,	4.81739369008347E-24,	5.13436300048453E-24,	5.46325498101433E-24,	5.80451370765942E-24,	6.15859962997907E-24,	6.5259901685511E-24,	6.90718033398432E-24,	7.3026833682673E-24,	7.71303140925068E-24,	8.13877617908807E-24,	8.58048969749003E-24,	9.03876502067529E-24,	9.51421700693493E-24,	1.00074831097573E-23,	1.05192241994947E-23,	1.1050125414588E-23,	1.16008970433997E-23,	1.21722754377452E-23,	1.27650239592473E-23,	1.33799339596817E-23,	1.40178257965197E-23,	1.4679549884918E-23,	1.53659877874493E-23,	1.60780533429121E-23,	1.68166938356051E-23,	1.75828912065022E-23,	1.83776633078098E-23,	1.92020652024465E-23,	2.00571905100317E-23,	2.09441728010326E-23,	2.1864187040771E-23,	2.28184510850531E-23,	2.38082272292489E-23,	2.48348238127066E-23,	2.58995968804601E-23,	2.70039519042485E-23,	2.81493455649444E-23,	2.93372875985543E-23,	3.05693427080343E-23,	3.18471325432418E-23,	3.31723377514212E-23,	3.45467001007111E-23,	3.59720246792415E-23,	3.7450182172482E-23,	3.89831112215941E-23,	4.05728208656349E-23,	4.2221393070561E-23,	4.39309853480814E-23,	4.57038334675155E-23,	4.75422542639206E-23,	4.94486485458695E-23,	5.14255041063705E-23,	5.347539884055E-23,	5.56010039738374E-23,	5.78050874045255E-23,	6.00905171647109E-23,	6.24602650037603E-23,	6.4917410098589E-23,	6.74651428951888E-23,	7.01067690859936E-23,	7.28457137278321E-23,	7.56855255053773E-23,	7.86298811451768E-23,	8.16825899855184E-23,	8.48475987075712E-23,	8.81289962334252E-23,	9.15310187968526E-23,	9.50580551928064E-23,	9.87146522118878E-23,	1.0250552026622E-22,	1.06435539213398E-22,	1.10509764385392E-22,	1.14733432829563E-22,	1.19111969769125E-22,	1.23650995290717E-22,	1.28356331266958E-22,	1.33234008522141E-22,	1.38290274249521E-22,	1.43531599688916E-22,	1.48964688073647E-22,	1.54596482856172E-22,	1.60434176222048E-22,	1.66485217902237E-22,	1.72757324294063E-22,	1.79258487901516E-22,	1.85996987105951E-22,	1.92981396278611E-22,	2.00220596246783E-22,	2.07723785125818E-22,	2.15500489529649E-22,	2.23560576172872E-22,	2.31914263877916E-22,	2.40572136001262E-22,	2.49545153293186E-22,	2.58844667205941E-22,	2.68482433665862E-22,	2.78470627325347E-22,	2.88821856311256E-22,	2.99549177486801E-22,	3.10666112244601E-22,	3.22186662849159E-22,	3.34125329347653E-22,	3.46497127068572E-22,	3.59317604728372E-22,	3.72602863167051E-22,	3.86369574734196E-22,	4.00635003347844E-22,	4.1541702524921E-22,	4.30734150477136E-22,	4.46605545086929E-22,	4.63051054139065E-22,	4.80091225484128E-22,	4.97747334371226E-22,	5.16041408908048E-22,	5.34996256401681E-22,	5.54635490610295E-22,	5.74983559936797E-22,	5.96065776596649E-22,	6.17908346793064E-22,	6.40538401933993E-22,	6.63984030926397E-22,	6.88274313584558E-22,	7.13439355190379E-22,	7.39510322244916E-22,	7.66519479451701E-22,	7.94500227973794E-22,	8.23487145007867E-22,	8.5351602472015E-22,	8.84623920590492E-22,	9.1684918921241E-22,	9.50231535598564E-22,	9.84812060042781E-22,	1.02063330659143E-21,	1.05773931317877E-21,	1.09617566348265E-21,	1.13598954055892E-21,	1.17722978231481E-21,	1.21994693888346E-21,	1.26419333196413E-21,	1.31002311619443E-21,	1.35749234262342E-21,	1.40665902435656E-21,	1.45758320444586E-21,	1.51032702610123E-21,	1.56495480530117E-21,	1.62153310588403E-21,	1.68013081720337E-21,	1.74081923443391E-21,	1.80367214161752E-21,	1.86876589754144E-21,	1.9361795245442E-21,	2.00599480034786E-21,	2.07829635301828E-21,	2.15317175915874E-21,	2.23071164544562E-21,	2.31100979361845E-21,	2.39416324904039E-21,	2.48027243294904E-21,	2.5694412585215E-21,	2.66177725088156E-21,	2.75739167118135E-21,	2.85639964489403E-21,	2.95892029445854E-21,	3.06507687642232E-21,	3.17499692323255E-21,	3.28881238983148E-21,	3.40665980521656E-21,	3.5286804291315E-21,	3.65502041405957E-21,	3.78583097269654E-21,	3.9212685510861E-21,	4.06149500760679E-21,	4.20667779800588E-21,	4.35699016668168E-21,	4.51261134442275E-21,	4.67372675281919E-21,	4.84052821556819E-21,	5.01321417690354E-21,	5.19198992738612E-21,	5.37706783730029E-21,	5.56866759790919E-21,	5.76701647083012E-21,	5.97234954579985E-21,	6.1849100071085E-21,	6.40494940898987E-21,	6.63272796026529E-21,	6.86851481854806E-21,	7.1125883943255E-21,	7.36523666524585E-21,	7.62675750094823E-21,	7.89745899878474E-21,	8.17765983079517E-21,	8.46768960230682E-21,	8.76788922254357E-21,	9.07861128764159E-21,	9.40022047648131E-21,	9.73309395975918E-21,	1.00776218227363E-20,	1.0434207502115E-20,	1.080326823751E-20,	1.1185235537995E-20,	1.15805556642207E-20,	1.19896901266199E-20,	1.24131162002255E-20,	1.28513274566523E-20,	1.33048343138038E-20,	1.37741646038901E-20,	1.42598641603555E-20,	1.47624974243405E-20,	1.52826480713176E-20,	1.58209196585632E-20,	1.63779362941514E-20,	1.69543433281729E-20,	1.75508080669101E-20,	1.81680205107193E-20,	1.88066941163971E-20,	1.9467566584833E-20,	2.0151400674775E-20,	2.08589850435642E-20,	2.15911351157179E-20,	2.23486939802743E-20,	2.31325333178362E-20,	2.3943554358285E-20,	2.47826888701661E-20,	2.56509001827787E-20,	2.65491842420374E-20,	2.74785707012057E-20,	2.84401240476395E-20,	2.94349447667114E-20,	3.04641705441284E-20,	3.15289775078917E-20,	3.26305815111881E-20,	3.37702394575442E-20,	3.49492506696174E-20,	3.61689583030412E-20,	3.74307508067876E-20,	3.8736063431557E-20,	4.00863797877544E-20,	4.14832334546585E-20,	4.29282096424457E-20,	4.44229469087794E-20,	4.59691389317332E-20,	4.75685363408717E-20,	4.92229486083708E-20,	5.09342460021201E-20,	5.270436160281E-20,	5.45352933870736E-20,	5.64291063788155E-20,	5.838793487093E-20,	6.04139847196803E-20,	6.25095357140838E-20,	6.46769440227205E-20,	6.69186447204616E-20,	6.92371543976934E-20,	7.16350738546913E-20,	7.41150908838876E-20,	7.66799831428588E-20,	7.93326211209521E-20,	8.20759712025599E-20,	8.49130988301485E-20,	8.7847171770246E-20,	9.08814634856924E-20,	9.40193566175654E-20,	9.72643465802964E-20,	1.00620045273609E-19,	1.04090184915022E-19,	1.07678621996777E-19,	1.11389341371184E-19,	1.15226460468478E-19,	1.19194233651447E-19,	1.23297056711185E-19,	1.27539471508495E-19,	1.31926170765587E-19,	1.36462003012875E-19,	1.41151977695823E-19,	1.46001270446941E-19,	1.51015228528209E-19,	1.56199376449352E-19,	1.61559421767569E-19,	1.67101261074506E-19,	1.72830986176417E-19,	1.78754890473679E-19,	1.84879475545984E-19,	1.91211457949757E-19,	1.97757776234547E-19,	2.04525598185328E-19,	2.11522328297915E-19,	2.18755615494848E-19,	2.26233361089422E-19,	2.33963727005686E-19,	2.41955144262551E-19,	2.50216321730365E-19,	2.58756255168578E-19,	2.67584236553399E-19,	2.76709863704618E-19,	2.86143050221045E-19,	2.95894035734335E-19,	3.05973396491244E-19,	3.16392056274687E-19,	3.27161297674306E-19,	3.38292773717559E-19,	3.49798519872702E-19,	3.61690966435396E-19,	3.73982951311012E-19,	3.86687733205103E-19,	3.99819005234891E-19,	4.13390908975011E-19,	4.27418048951174E-19,	4.41915507595826E-19,	4.56898860680314E-19,	4.72384193238544E-19,	4.88388115997524E-19,	5.04927782330741E-19,	5.22020905750725E-19,	5.39685777957743E-19,	5.57941287462018E-19,	5.76806938797449E-19,	5.96302872345359E-19,	6.16449884787333E-19,	6.37269450206854E-19,	6.58783741859991E-19,	6.8101565463607E-19,	7.03988828229862E-19,	7.27727671047508E-19,	7.52257384869085E-19,	7.77603990291395E-19,	8.0379435297532E-19,	8.3085621072281E-19,	8.5881820140934E-19,	8.87709891798475E-19,	9.17561807265996E-19,	9.48405462461878E-19,	9.80273392939267E-19,	1.01319918778054E-18,	1.04721752325135E-18,	1.08236419751472E-18,	1.11867616643783E-18,	1.15619158052571E-18,	1.19494982301649E-18,	1.2349915491744E-18,	1.27635872681752E-18,	1.31909467811855E-18,	1.363244122718E-18,	1.40885322219036E-18,	1.45596962590514E-18,	1.50464251832583E-18,	1.55492266779127E-18,	1.60686247682509E-18,	1.66051603402051E-18,	1.71593916754901E-18,	1.77318950034293E-18,	1.83232650700369E-18,	1.89341157248865E-18,	1.9565080526315E-18,	2.02168133655246E-18,	2.08899891101664E-18,	2.15853042680022E-18,	2.23034776712636E-18,	2.30452511823423E-18,	2.38113904214689E-18,	2.46026855170519E-18,	2.54199518793758E-18,	2.62640309983707E-18,	2.7135791266194E-18,	2.80361288253828E-18,	2.89659684433603E-18,	2.99262644141028E-18,	3.0918001487798E-18,	3.19421958293503E-18,	3.29998960066153E-18,	3.40921840092706E-18,	3.52201762992595E-18,	3.63850248937697E-18,	3.75879184817412E-18,	3.88300835749238E-18,	4.01127856945382E-18,	4.14373305946253E-18,	4.28050655231999E-18,	4.42173805223601E-18,	4.56757097685362E-18,	4.71815329541011E-18,	4.87363767115974E-18,	5.03418160818766E-18,	5.19994760274835E-18,	5.37110329926585E-18,	5.54782165113719E-18,	5.73028108648459E-18,	5.91866567900653E-18,	6.11316532408191E-18,	6.3139759202865E-18,	6.52129955648542E-18,	6.73534470467017E-18,	6.95632641871406E-18,	7.18446653922475E-18,	7.41999390467799E-18,	7.66314456902238E-18,	7.91416202595011E-18,	8.17329744003498E-18,	8.44080988494456E-18,	8.71696658893962E-18,	9.00204318788041E-18,	9.29632398596556E-18,	9.60010222443646E-18,	9.91368035848649E-18,	1.02373703426218E-17,	1.05714939247275E-17,	1.09163829491009E-17,	1.12723796687201E-17,	1.16398370670266E-17,	1.20191191895055E-17,	1.24106014853581E-17,	1.28146711595682E-17,	1.32317275356747E-17,	1.36621824295699E-17,	1.41064605346528E-17,	1.45649998186792E-17,	1.50382519326548E-17,	1.55266826321345E-17,	1.60307722112954E-17,	1.65510159501664E-17,	1.70879245754068E-17,	1.76420247350364E-17,	1.82138594875346E-17,	1.88039888057354E-17,	1.94129900959588E-17,	2.00414587328326E-17,	2.06900086102705E-17,	2.13592727090874E-17,	2.20499036817454E-17,	2.27625744547405E-17,	2.34979788491519E-17,	2.42568322198947E-17,	2.50398721142289E-17,	2.58478589500966E-17,	2.66815767148743E-17,	2.75418336851448E-17,	2.84294631681111E-17,	2.93453242652922E-17,	3.02903026591597E-17,	3.12653114233934E-17,	3.22712918574531E-17,	3.33092143461848E-17,	3.43800792452004E-17,	3.54849177927902E-17,	3.66247930491514E-17,	3.78008008637366E-17,	3.90140708715513E-17,	4.0265767519252E-17,	4.15570911219218E-17,	4.28892789514261E-17,	4.42636063572761E-17,	4.56813879209553E-17,	4.71439786446923E-17,	4.86527751756899E-17,	5.02092170668511E-17,	5.18147880750726E-17,	5.34710174982059E-17,	5.51794815518196E-17,	5.69418047869276E-17,	5.87596615498817E-17,	6.06347774856636E-17,	6.25689310858426E-17,	6.45639552825058E-17,	6.66217390895037E-17,	6.87442292923909E-17,	7.09334321884842E-17,	7.31914153784984E-17,	7.55203096112646E-17,	7.79223106830762E-17,	8.03996813932537E-17,	8.2954753557565E-17,	8.55899300811839E-17,	8.83076870929171E-17,	9.11105761424825E-17,	9.4001226462667E-17,	9.69823472982508E-17,	1.00056730303633E-16,	1.03227252011151E-16,	1.06496876372147E-16,	1.09868657372886E-16,	1.13345741727488E-16,	1.16931371650117E-16,	1.20628887708706E-16,	1.24441731762586E-16,	1.28373449986434E-16,	1.32427695983045E-16,	1.36608233987482E-16,	1.4091894216525E-16,	1.45363816007196E-16,	1.49946971823933E-16,	1.54672650342657E-16,	1.59545220409285E-16,	1.64569182798984E-16,	1.69749174138165E-16,	1.75089970941186E-16,	1.8059649376504E-16,	1.86273811485414E-16,	1.92127145697625E-16,	1.98161875245996E-16,	2.04383540885365E-16,	2.10797850078515E-16,	2.17410681933421E-16,	2.24228092284313E-16,	2.31256318920677E-16,	2.38501786968427E-16,	2.45971114427597E-16,	2.53671117871025E-16,	2.61608818308634E-16,	2.69791447222033E-16,	2.78226452774306E-16,	2.86921506199971E-16,	2.9588450838027E-16,	3.05123596609043E-16,	3.14647151554634E-16,	3.2446380442341E-16,	3.34582444330605E-16,	3.45012225884423E-16,	3.55762576989423E-16,	3.66843206875444E-16,	3.7826411435846E-16,	3.90035596339945E-16,	4.02168256551527E-16,	4.14673014551867E-16,	4.27561114982927E-16,	4.40844137092957E-16,	4.54534004533768E-16,	4.68642995440033E-16,	4.83183752798604E-16,	4.98169295116041E-16,	5.13613027392764E-16,	5.29528752412506E-16,	5.45930682355945E-16,	5.62833450747666E-16,	5.80252124745854E-16,	5.98202217784371E-16,	6.16699702577132E-16,	6.35761024494994E-16,	6.55403115325624E-16,	6.75643407427115E-16,	6.96499848286429E-16,	7.17990915494015E-16,	7.40135632146311E-16,	7.62953582688109E-16,	7.86464929207146E-16,	8.10690428193571E-16,	8.35651447777335E-16,	8.61369985456862E-16,	8.87868686332773E-16,	9.15170861860777E-16,	9.43300509138238E-16,	9.7228233073935E-16,	1.00214175511422E-15,	1.0329049575676E-15,	1.06459888183347E-15,	1.09725126226207E-15,	1.13089064663635E-15,	1.16554641963564E-15,	1.20124882696418E-15,	1.23802900016338E-15,	1.27591898212656E-15,	1.31495175333581E-15,	1.35516125884102E-15,	1.39658243600171E-15,	1.43925124301276E-15,	1.48320468823583E-15,	1.52848086035881E-15,	1.57511895940626E-15,	1.62315932862432E-15,	1.67264348726442E-15,	1.72361416429054E-15,	1.77611533303562E-15,	1.83019224683339E-15,	1.8858914756524E-15,	1.94326094376013E-15,	2.00234996844543E-15,	2.06320929982858E-15,	2.12589116178887E-15,	2.19044929404057E-15,	2.25693899538882E-15,	2.32541716819796E-15,	2.39594236410569E-15,	2.46857483101714E-15,	2.54337656141426E-15,	2.62041134201641E-15,	2.69974480482942E-15,	2.78144447962104E-15,	2.86557984786209E-15,	2.95222239817323E-15,	3.04144568331879E-15,	3.13332537878986E-15,	3.22793934302029E-15,	3.32536767927999E-15,	3.4256927992916E-15,	3.52899948861756E-15,	3.63537497386571E-15,	3.7449089917634E-15,	3.85769386015065E-15,	3.97382455094502E-15,	4.09339876513163E-15,	4.21651700983369E-15,	4.34328267751993E-15,	4.47380212740725E-15,	4.60818476911811E-15,	4.74654314865398E-15,	4.88899303674773E-15,	5.0356535196595E-15,	5.18664709248233E-15,	5.34209975502564E-15,	5.50214111034624E-15,	5.66690446599878E-15,	5.83652693807904E-15,	6.01114955813578E-15,	6.1909173830285E-15,	6.37597960781095E-15,	6.56648968172182E-15,	6.76260542736672E-15,	6.96448916317744E-15,	7.17230782923677E-15,	7.38623311655971E-15,	7.60644159992404E-15,	7.83311487434591E-15,	8.06643969529833E-15,	8.30660812277337E-15,	8.55381766929128E-15,	8.80827145196255E-15,	9.07017834871162E-15,	9.33975315877412E-15,	9.61721676758191E-15,	9.90279631615382E-15,	1.01967253751125E-14,	1.04992441234515E-14,	1.08105995321794E-14,	1.11310455529713E-14,	1.14608433119625E-14,	1.18002613088197E-14,	1.21495756212332E-14,	1.25090701149723E-14,	1.28790366596538E-14,	1.32597753503755E-14,	1.36515947353704E-14,	1.40548120498428E-14,	1.44697534561498E-14,	1.48967542904981E-14,	1.53361593163274E-14,	1.57883229845609E-14,	1.62536097009018E-14,	1.67323941003658E-14,	1.72250613292391E-14,	1.77320073346609E-14,	1.82536391620302E-14,	1.87903752604456E-14,	1.93426457963902E-14,	1.99108929758794E-14,	2.04955713752952E-14,	2.10971482811374E-14,	2.1716104038926E-14,	2.23529324114963E-14,	2.30081409469359E-14,	2.36822513564152E-14,	2.43757999021748E-14,	2.50893377959348E-14,	2.58234316080015E-14,	2.65786636873534E-14,	2.73556325929928E-14,	2.81549535368612E-14,	2.89772588386214E-14,	2.98231983926171E-14,	3.06934401473295E-14,	3.15886705976595E-14,	3.25095952903697E-14,	3.34569393430322E-14,	3.44314479768345E-14,	3.54338870636071E-14,	3.64650436874432E-14,	3.7525726721292E-14,	3.86167674189172E-14,	3.97390200226209E-14,	4.08933623871439E-14,	4.20806966201645E-14,	4.33019497398278E-14,	4.45580743497484E-14,	4.58500493319424E-14,	4.7178880558153E-14,	4.85456016200496E-14,	4.99512745787892E-14,	5.13969907344434E-14,	5.28838714158055E-14,	5.4413068791108E-14,	5.59857667001895E-14,	5.76031815086692E-14,	5.9266562984697E-14,	6.09771951988635E-14,	6.27363974478693E-14,	6.4545525202566E-14,	6.64059710810005E-14,	6.83191658471058E-14,	7.02865794357005E-14,	7.23097220044761E-14,	7.43901450136659E-14,	7.65294423341095E-14,	7.87292513844429E-14,	8.09912542981638E-14,	8.33171791213394E-14,	8.57088010417448E-14,	8.81679436502376E-14,	9.06964802351979E-14,	9.32963351108791E-14,	9.59694849805411E-14,	9.87179603352554E-14,	1.01543846889295E-13,	1.04449287053048E-13,	1.07436481444407E-13,	1.10507690439632E-13,	1.13665235764674E-13,	1.16911502128016E-13,	1.20248938896067E-13,	1.23680061812213E-13,	1.2720745476063E-13,	1.30833771575992E-13,	1.34561737900255E-13,	1.38394153087706E-13,	1.42333892159507E-13,	1.46383907808988E-13,	1.50547232458986E-13,	1.54826980372538E-13,	1.59226349818292E-13,	1.63748625292022E-13,	1.6839717979566E-13,	1.73175477175307E-13,	1.78087074519722E-13,	1.83135624620799E-13,	1.88324878497624E-13,	1.93658687985694E-13,	1.99141008392962E-13,	2.04775901224376E-13,	2.10567536976658E-13,	2.16520198005076E-13,	2.22638281464033E-13,	2.28926302323327E-13,	2.3538889646198E-13,	2.42030823841601E-13,	2.48856971761262E-13,	2.55872358195947E-13,	2.63082135220668E-13,	2.70491592522395E-13,	2.78106161001999E-13,	2.85931416468469E-13,	2.93973083427704E-13,	3.02237038968256E-13,	3.10729316746443E-13,	3.19456111073308E-13,	3.28423781105986E-13,	3.37638855146063E-13,	3.47108035047616E-13,	3.56838200737648E-13,	3.66836414851739E-13,	3.77109927487764E-13,	3.87666181080628E-13,	3.98512815401017E-13,	4.09657672681263E-13,	4.21108802871458E-13,	4.32874469029072E-13,	4.44963152845375E-13,	4.57383560312056E-13,	4.70144627531514E-13,	4.83255526674385E-13,	4.96725672087933E-13,	5.10564726559059E-13,	5.24782607735725E-13,	5.39389494710727E-13,	5.54395834771816E-13,	5.69812350322268E-13,	5.85650045976117E-13,	6.01920215832341E-13,	6.18634450932426E-13,	6.358046469058E-13,	6.53443011807776E-13,	6.71562074154736E-13,	6.90174691161384E-13,	7.0929405718505E-13,	7.28933712382116E-13,	7.49107551581767E-13,	7.69829833382389E-13,	7.91115189476081E-13,	8.12978634206849E-13,	8.35435574368205E-13,	8.58501819246025E-13,	8.8219359091266E-13,	9.06527534778431E-13,	9.31520730406787E-13,	9.57190702599565E-13,	9.83555432758924E-13,	1.0106333705327E-12,	1.03844344575009E-12,	1.06700508065469E-12,	1.09633820244217E-12,	1.12646325610993E-12,	1.1574012176264E-12,	1.18917360742763E-12,	1.22180250424915E-12,	1.25531055930131E-12,	1.28972101079629E-12,	1.32505769883537E-12,	1.3613450806651E-12,	1.39860824631148E-12,	1.43687293460099E-12,	1.47616554957807E-12,	1.51651317732847E-12,	1.55794360321834E-12,	1.60048532955896E-12,	1.6441675937076E-12,	1.68902038661467E-12,	1.73507447182818E-12,	1.78236140496643E-12,	1.83091355367008E-12,	1.88076411804517E-12,	1.93194715160899E-12,	1.98449758275063E-12,	2.03845123671876E-12,	2.09384485814906E-12,	2.15071613414445E-12,	2.20910371792116E-12,	2.26904725303417E-12,	2.33058739819607E-12,	2.39376585270317E-12,	2.45862538248367E-12,	2.52520984678254E-12,	2.59356422549827E-12,	2.66373464718718E-12,	2.73576841775097E-12,	2.8097140498239E-12,	2.88562129287615E-12,	2.96354116405047E-12,	3.04352597974938E-12,	3.12562938799086E-12,	3.20990640155068E-12,	3.29641343191001E-12,	3.38520832402739E-12,	3.47635039195457E-12,	3.56990045531619E-12,	3.66592087667361E-12,	3.76447559979395E-12,	3.86563018884556E-12,	3.96945186854178E-12,	4.07600956525542E-12,	4.18537394912676E-12,	4.29761747718851E-12,	4.41281443753156E-12,	4.53104099453612E-12,	4.65237523519314E-12,	4.77689721654172E-12,	4.90468901424864E-12,	5.03583477235684E-12,	5.17042075423013E-12,	5.30853539472227E-12,	5.45026935359903E-12,	5.59571557024252E-12,	5.74496931966771E-12,	5.89812826988194E-12,	6.05529254061859E-12,	6.21656476347708E-12,	6.38205014350196E-12,	6.5518565222346E-12,	6.72609444227172E-12,	6.90487721336591E-12,	7.08832098010398E-12,	7.27654479119965E-12,	7.4696706704383E-12,	7.66782368931192E-12,	7.87113204138358E-12,	8.07972711842136E-12,	8.29374358834294E-12,	8.51331947501247E-12,	8.73859623993282E-12,	8.96971886587685E-12,	9.20683594250257E-12,	9.45009975399792E-12,	9.69966636880201E-12,	9.95569573145076E-12,	1.02183517565958E-11,	1.04878024252464E-11,	1.07642198832865E-11,	1.10477805423177E-11,	1.13386651828829E-11,	1.16370590601247E-11,	1.19431520119346E-11,	1.22571385696495E-11,	1.25792180713552E-11,	1.29095947778552E-11,	1.32484779913664E-11,	1.35960821770044E-11,	1.39526270871211E-11,	1.43183378885611E-11,	1.46934452929027E-11,	1.50781856897508E-11,	1.54728012831535E-11,	1.58775402312115E-11,	1.62926567889531E-11,	1.67184114545506E-11,	1.71550711189525E-11,	1.760290921901E-11,	1.80622058941761E-11,	1.85332481468596E-11,	1.90163300065165E-11,	1.95117526975621E-11,	2.00198248111924E-11,	2.05408624812021E-11,	2.10751895638888E-11,	2.16231378221381E-11,	2.21850471137808E-11,	2.27612655843214E-11,	2.33521498641351E-11,	2.39580652702333E-11,	2.45793860127028E-11,	2.52164954059206E-11,	2.58697860846551E-11,	2.65396602251607E-11,	2.72265297713789E-11,	2.79308166663619E-11,	2.86529530890327E-11,	2.93933816964045E-11,	3.01525558713789E-11,	3.09309399762504E-11,	3.1729009612042E-11,	3.25472518838043E-11,	3.33861656720106E-11,	3.42462619101834E-11,	3.51280638688919E-11,	3.60321074462629E-11,	3.69589414651486E-11,	3.79091279771014E-11,	3.88832425733046E-11,	3.98818747026168E-11,	4.09056279968844E-11,	4.1955120603686E-11,	4.30309855266732E-11,	4.41338709736739E-11,	4.52644407127339E-11,	4.6423374436268E-11,	4.76113681335039E-11,	4.88291344713991E-11,	5.00774031842199E-11,	5.13569214719724E-11,	5.26684544078813E-11,	5.40127853551159E-11,	5.53907163929661E-11,	5.68030687526772E-11,	5.8250683263156E-11,	5.97344208067638E-11,	6.12551627854193E-11,	6.28138115972364E-11,	6.44112911239286E-11,	6.60485472292159E-11,	6.77265482684741E-11,	6.9446285609874E-11,	7.12087741672604E-11,	7.30150529450279E-11,	7.48661855952556E-11,	7.67632609873673E-11,	7.87073937905915E-11,	8.06997250694989E-11,	8.27414228929029E-11,	8.4833682956413E-11,	8.69777292189384E-11,	8.91748145534452E-11,	9.14262214122753E-11,	9.3733262507344E-11,	9.6097281505539E-11,	9.85196537396484E-11,	1.01001786935156E-10,	1.03545121953244E-10,	1.06151133550359E-10,	1.08821331154692E-10,	1.1155725965994E-10,	1.14360500236729E-10,	1.17232671162062E-10,	1.20175428667204E-10,	1.23190467804373E-10,	1.26279523332669E-10,	1.29444370623634E-10,	1.32686826586873E-10,	1.36008750616163E-10,	1.394120455565E-10,	1.42898658692512E-10,	1.46470582758717E-10,	1.50129856972078E-10,	1.53878568087341E-10,	1.57718851475639E-10,	1.61652892226852E-10,	1.65682926276248E-10,	1.69811241555902E-10,	1.74040179171431E-10,	1.78372134604585E-10,	1.82809558942241E-10,	1.87354960132361E-10,	1.92010904267485E-10,	1.96780016896345E-10,	2.01664984364205E-10,	2.06668555182516E-10,	2.11793541428522E-10,	2.17042820175451E-10,	2.22419334953928E-10,	2.27926097245275E-10,	2.33566188007374E-10,	2.39342759233774E-10,	2.45259035546749E-10,	2.51318315825011E-10,	2.57523974866824E-10,	2.63879465089251E-10,	2.70388318264291E-10,	2.77054147292697E-10,	2.83880648016246E-10,	2.90871601069284E-10,	2.9803087377035E-10,	3.05362422054744E-10,	3.12870292448866E-10,	3.20558624087227E-10,	3.28431650773001E-10,	3.36493703083038E-10,	3.44749210518269E-10,	3.53202703700436E-10,	3.61858816616123E-10,	3.70722288909066E-10,	3.79797968221755E-10,	3.89090812587335E-10,	3.98605892872867E-10,	4.08348395275009E-10,	4.183236238692E-10,	4.2853700321346E-10,	4.38994081007931E-10,	4.49700530811319E-10,	4.60662154815411E-10,	4.71884886678854E-10,	4.83374794421447E-10,	4.95138083380157E-10,	5.07181099228163E-10,	5.19510331058205E-10,	5.32132414531574E-10,	5.45054135094079E-10,	5.5828243126038E-10,	5.71824397968078E-10,	5.85687290003007E-10,	5.99878525497167E-10,	6.14405689500807E-10,	6.29276537630156E-10,	6.44498999792365E-10,	6.60081183989225E-10,	6.76031380201269E-10,	6.92358064353915E-10,	7.0906990236729E-10,	7.26175754291474E-10,	7.43684678528878E-10,	7.61605936145532E-10,	7.79948995273101E-10,	7.98723535603454E-10,	8.17939452977677E-10,	8.37606864071432E-10,	8.57736111178627E-10,	8.78337767095371E-10,	8.99422640106256E-10,	9.21001779075024E-10,	9.43086478641725E-10,	9.65688284528525E-10,	9.88818998956343E-10,	1.01249068617455E-09,	1.03671567810601E-09,	1.06150658010978E-09,	1.08687627686384E-09,	1.11283793837025E-09,	1.13940502608516E-09,	1.1665912991763E-09,	1.19441082091029E-09,	1.22287796517262E-09,	1.2520074231227E-09,	1.28181420998691E-09,	1.31231367199223E-09,	1.34352149344333E-09,	1.37545370394598E-09,	1.40812668577969E-09,	1.44155718142253E-09,	1.47576230123116E-09,	1.51075953127917E-09,	1.5465667413569E-09,	1.58320219313582E-09,	1.6206845485009E-09,	1.65903287805417E-09,	1.69826666979291E-09,	1.73840583796591E-09,	1.77947073211132E-09,	1.82148214627968E-09,	1.86446132844576E-09,	1.90842999011295E-09,	1.953410316114E-09,	1.99942497461191E-09,	2.04649712730504E-09,	2.09465043984026E-09,	2.14390909243838E-09,	2.19429779073597E-09,	2.24584177684779E-09,	2.29856684065414E-09,	2.3524993313176E-09,	2.40766616903355E-09,	2.46409485701911E-09,	2.52181349374509E-09,	2.58085078541577E-09,	2.64123605870116E-09,	2.70299927372697E-09,	2.7661710373269E-09,	2.83078261656271E-09,	2.896865952517E-09,	2.96445367436418E-09,	3.0335791137249E-09,	3.10427631930946E-09,	3.17658007185582E-09,	3.25052589936785E-09,	3.32615009265973E-09,	3.40348972121221E-09,	3.48258264934693E-09,	3.5634675527249E-09,	3.64618393517522E-09,	3.73077214586059E-09,	3.81727339678602E-09,	3.90572978065721E-09,	3.99618428909547E-09,	4.08868083121599E-09,	4.18326425257624E-09,	4.27998035450186E-09,	4.37887591379701E-09,	4.47999870284669E-09,	4.58339751011839E-09,	4.68912216107081E-09,	4.79722353947722E-09,	4.90775360917161E-09,	5.02076543622538E-09,	5.13631321156302E-09,	5.2544522740249E-09,	5.37523913388579E-09,	5.49873149683767E-09,	5.62498828844566E-09,	5.75406967908602E-09,	5.88603710937529E-09,	6.02095331609991E-09,	6.15888235865569E-09,	6.29988964600683E-09,	6.44404196417416E-09,	6.59140750426272E-09,	6.74205589103863E-09,	6.89605821206571E-09,	7.05348704741232E-09,	7.21441649993912E-09,	7.37892222617855E-09,	7.54708146781725E-09,	7.71897308379259E-09,	7.89467758301482E-09,	8.07427715772654E-09,	8.2578557175113E-09,	8.44549892396358E-09,	8.63729422603229E-09,	8.8333308960504E-09,	9.03370006646353E-09,	9.2384947672703E-09,	9.44780996418785E-09,	9.66174259755578E-09,	9.88039162199232E-09,	1.01038580468166E-08,	1.03322449772513E-08,	1.05656576564197E-08,	1.08042035081522E-08,	1.10479921806174E-08,	1.12971355907924E-08,	1.15517479697876E-08,	1.18119459090427E-08,	1.20778484074094E-08,	1.23495769191363E-08,	1.26272554027736E-08,	1.29110103710144E-08,	1.32009709414892E-08,	1.34972688885307E-08,	1.38000386959281E-08,	1.41094176106875E-08,	1.44255456978173E-08,	1.47485658961576E-08,	1.50786240752718E-08,	1.54158690934205E-08,	1.57604528566371E-08,	1.61125303789248E-08,	1.64722598435957E-08,	1.68398026657726E-08,	1.72153235560746E-08,	1.75989905855077E-08,	1.79909752515823E-08,	1.83914525456796E-08,	1.88006010216903E-08,	1.92186028659469E-08,	1.9645643968474E-08,	2.00819139955807E-08,	2.05276064638177E-08,	2.09829188153253E-08,	2.14480524945956E-08,	2.19232130266761E-08,	2.24086100968384E-08,	2.29044576317402E-08,	2.34109738821055E-08,	2.39283815069517E-08,	2.445690765939E-08,	2.49967840740276E-08,	2.55482471560004E-08,	2.61115380716644E-08,	2.66869028409766E-08,	2.7274592431594E-08,	2.78748628547217E-08,	2.84879752627421E-08,	2.91141960486541E-08,	2.97537969473576E-08,	3.04070551388131E-08,	3.10742533531114E-08,	3.17556799774857E-08,	3.24516291653014E-08,	3.31624009470572E-08,	3.38883013434345E-08,	3.46296424804292E-08,	3.53867427066042E-08,	3.61599267124977E-08,	3.69495256522278E-08,	3.77558772673289E-08,	3.85793260128608E-08,	3.94202231858299E-08,	4.02789270559619E-08,	4.11558029988683E-08,	4.20512236316469E-08,	4.29655689509604E-08,	4.38992264736337E-08,	4.48525913798169E-08,	4.58260666587548E-08,	4.68200632572109E-08,	4.7835000230591E-08,	4.88713048968115E-08,	4.99294129929629E-08,	5.10097688348136E-08,	5.21128254792046E-08,	5.32390448893853E-08,	5.43888981033394E-08,	5.5562865405153E-08,	5.67614364994781E-08,	5.7985110689143E-08,	5.92343970559644E-08,	6.05098146448155E-08,	6.18118926510067E-08,	6.31411706110343E-08,	6.44981985967547E-08,	6.58835374130445E-08,	6.7297758799002E-08,	6.8741445632755E-08,	7.02151921399318E-08,	7.17196041058606E-08,	7.32552990915596E-08,	7.48229066535806E-08,	7.64230685677743E-08,	7.80564390570402E-08,	7.97236850231316E-08,	8.14254862825808E-08,	8.3162535806817E-08,	8.49355399665458E-08,	8.67452187804616E-08,	8.85923061683674E-08,	9.0477550208774E-08,	9.24017134010549E-08,	9.43655729322334E-08,	9.6369920948478E-08,	9.8415564831386E-08,	1.00503327479135E-07,	1.02634047592583E-07,	1.04808579966401E-07,	1.07027795785319E-07,	1.0929258292557E-07,	1.1160384626163E-07,	1.1396250797832E-07,	1.16369507888376E-07,	1.18825803755577E-07,	1.21332371623513E-07,	1.23890206150091E-07,	1.26500320947879E-07,	1.29163748930377E-07,	1.31881542664312E-07,	1.34654774728066E-07,	1.37484538076324E-07,	1.40371946411057E-07,	1.4331813455893E-07,	1.46324258855258E-07,	1.49391497534597E-07,	1.52521051128099E-07,	1.55714142867719E-07,	1.58972019097408E-07,	1.62295949691388E-07,	1.6568722847964E-07,	1.691471736807E-07,	1.72677128341912E-07,	1.76278460787231E-07,	1.7995256507271E-07,	1.83700861449803E-07,	1.87524796836592E-07,	1.91425845297083E-07,	1.95405508528693E-07,	1.99465316358056E-07,	2.03606827245298E-07,	2.07831628796892E-07,	2.12141338287258E-07,	2.16537603189224E-07,	2.21022101713509E-07,	2.25596543357361E-07,	2.30262669462499E-07,	2.35022253782508E-07,	2.39877103059842E-07,	2.44829057612577E-07,	2.49879991931079E-07,	2.55031815284746E-07,	2.6028647233897E-07,	2.65645943782502E-07,	2.71112246965359E-07,	2.76687436547469E-07,	2.82373605158196E-07,	2.88172884066936E-07,	2.94087443864952E-07,	3.00119495158622E-07,	3.06271289274278E-07,	3.12545118974838E-07,	3.18943319188378E-07,	3.25468267748873E-07,	3.32122386149267E-07,	3.38908140307081E-07,	3.45828041342742E-07,	3.52884646370851E-07,	3.60080559304568E-07,	3.67418431673334E-07,	3.74900963454139E-07,	3.82530903916528E-07,	3.9031105248158E-07,	3.98244259595066E-07,	4.06333427614999E-07,	4.14581511713815E-07,	4.22991520795395E-07,	4.31566518427163E-07,	4.40309623787498E-07,	4.49224012628684E-07,	4.58312918255644E-07,	4.67579632520698E-07,	4.77027506834596E-07,	4.8665995319406E-07,	4.96480445226112E-07,	5.06492519249413E-07,	5.16699775352903E-07,	5.27105878491986E-07,	5.37714559602533E-07,	5.48529616732973E-07,	5.59554916194746E-07,	5.70794393731403E-07,	5.8225205570662E-07,	5.93931980311429E-07,	6.05838318790948E-07,	6.17975296690896E-07,	6.30347215124213E-07,	6.42958452058063E-07,	6.55813463621542E-07,	6.689167854344E-07,	6.8227303395708E-07,	6.95886907862414E-07,	7.09763189429277E-07,	7.23906745958543E-07,	7.38322531211673E-07,	7.53015586872264E-07,	7.67991044030908E-07,	7.83254124693705E-07,	7.98810143314785E-07,	8.14664508353183E-07,	8.30822723854444E-07,	8.47290391057311E-07,	8.64073210025868E-07,	8.8117698130752E-07,	8.98607607617181E-07,	9.16371095548068E-07,	9.34473557309474E-07,	9.52921212491935E-07,	9.71720389860175E-07,	9.90877529174249E-07,	1.01039918303929E-06,	1.03029201878427E-06,	1.05056282037023E-06,	1.07121849032832E-06,	1.09226605172826E-06,	1.1137126501774E-06,	1.13556555585112E-06,	1.15783216555477E-06,	1.18052000481773E-06,	1.20363673002011E-06,	1.22719013055236E-06,	1.25118813100851E-06,	1.27563879341331E-06,	1.30055031948392E-06,	1.32593105292652E-06,	1.35178948176844E-06,	1.37813424072624E-06,	1.40497411361028E-06,	1.43231803576637E-06,	1.46017509655491E-06,	1.48855454186807E-06,	1.5174657766857E-06,	1.54691836767035E-06,	1.57692204580194E-06,	1.60748670905286E-06,	1.6386224251038E-06,	1.67033943410109E-06,	1.70264815145598E-06,	1.7355591706866E-06,	1.76908326630307E-06,	1.80323139673645E-06,	1.83801470731207E-06,	1.87344453326796E-06,	1.90953240281893E-06,	1.94629004026692E-06,	1.98372936915838E-06,	2.02186251548922E-06,	0.000002060701810958,	2.10025979626811E-06,	2.14054922447958E-06,	2.1815830644111E-06,	2.22337450409315E-06,	2.26593695427276E-06,	2.30928405197074E-06,	2.35342966409195E-06,	2.39838789108958E-06,	2.44417307068387E-06,	2.49079978163633E-06,	2.53828284757994E-06,	2.58663734090632E-06,	2.63587858671051E-06,	2.68602216679417E-06,	2.73708392372803E-06,	2.78907996497437E-06,	2.8420266670703E-06,	2.89594067987284E-06,	2.95083893086629E-06,	3.00673862953317E-06,	3.06365727178919E-06,	3.12161264448337E-06,	3.18062282996407E-06,	3.24070621071181E-06,	3.3018814740399E-06,	3.36416761686355E-06,	3.42758395053864E-06,	3.4921501057709E-06,	3.55788603759647E-06,	3.62481203043485E-06,	3.69294870321513E-06,	3.76231701457654E-06,	3.83293826814421E-06,	3.9048341178812E-06,	3.97802657351786E-06,	4.05253800605935E-06,	4.12839115337256E-06,	4.20560912585335E-06,	4.28421541217521E-06,	4.36423388512032E-06,	4.44568880749426E-06,	4.52860483812527E-06,	4.61300703794925E-06,	4.69892087618167E-06,	4.78637223657743E-06,	4.87538742377976E-06,	4.96599316975953E-06,	5.05821664034586E-06,	5.15208544184947E-06,	5.24762762777973E-06,	5.34487170565676E-06,	5.44384664391983E-06,	5.5445818789331E-06,	5.64710732209016E-06,	5.75145336701854E-06,	5.85765089688542E-06,	5.96573129180591E-06,	6.07572643635517E-06,	6.18766872718564E-06,	6.30159108075084E-06,	6.41752694113694E-06,	6.53551028800349E-06,	6.65557564463486E-06,	6.77775808610346E-06,	6.90209324754646E-06,	7.02861733255727E-06,	7.15736712169317E-06,	7.28837998110074E-06,	7.42169387126031E-06,	7.55734735585113E-06,	7.69537961073862E-06,	7.83583043308526E-06,	7.97874025058668E-06,	8.12415013083443E-06,	8.27210179080711E-06,	8.42263760649122E-06,	8.57580062263366E-06,	8.73163456262717E-06,	8.89018383853047E-06,	9.05149356122489E-06,	9.21560955070884E-06,	9.38257834653215E-06,	9.55244721837164E-06,	9.72526417675001E-06,	9.90107798389936E-06,	1.00799381647715E-05,	1.02618950181966E-05,	1.04469996281918E-05,	1.06353038754223E-05,	1.08268604488155E-05,	1.10217228573315E-05,	1.12199454418908E-05,	1.14215833874612E-05,	1.16266927353063E-05,	1.18353303953972E-05,	0.000012047554158989,	0.000012263422711365,	1.24829956447491E-05,	1.27063334713903E-05,	1.29334976368197E-05,	1.31645505332824E-05,	1.33995555133475E-05,	1.36385769036964E-05,	1.38816800190926E-05,	1.41289311765362E-05,	0.000014380397709602,	0.000014636147982967,	1.48962514071274E-05,	1.51607784533073E-05,	1.54298006685625E-05,	1.57033906910801E-05,	0.000015981622265678,	1.62645702595049E-05,	1.65523106779436E-05,	1.68449206807208E-05,	1.71424785982245E-05,	1.74450639480322E-05,	1.77527574516521E-05,	1.80656410514793E-05,	0.00001838379792797,	1.87073125170362E-05,	1.90362705276628E-05,	0.00001937075895975,	1.97108661221834E-05,	0.000020056681651135,	2.04082965285962E-05,	2.07658031011469E-05,	0.000021129295098963,	2.14988676550637E-05,	2.18746173248033E-05,	2.22566421056086E-05,	2.26450414569653E-05,	2.30399163206558E-05,	0.000023441369141252,	2.38495038868643E-05,	2.42644260701517E-05,	2.46862427695939E-05,	2.51150626510293E-05,	0.000025550995989462,	2.59941546911396E-05,	2.64446523159061E-05,	2.69026040998317E-05,	2.73681269781237E-05,	2.78413396083197E-05,	2.83223623937685E-05,	2.88113175073999E-05,	2.93083289157867E-05,	2.98135224035035E-05,	3.03270255977836E-05,	3.08489679934782E-05,	3.13794809783209E-05,	3.19186978585007E-05,	3.24667538845463E-05,	3.30237862775263E-05,	0.000033589934255567,	3.41653390606919E-05,	3.47501439859867E-05,	0.000035344494403092,	3.59485377900282E-05,	0.000036562423759356,	3.71863040866746E-05,	3.78203327394631E-05,	0.000038464665906267,	3.91194620262342E-05,	3.97848818190038E-05,	4.04610883149515E-05,	4.11482468857946E-05,	4.18465252755612E-05,	4.25560936319267E-05,	4.32771245379213E-05,	4.40097930440123E-05,	4.47542767005656E-05,	4.55107555906885E-05,	4.62794123634602E-05,	0.000047060432267551,	4.78540031852369E-05,	4.86603156668106E-05,	4.94795629653954E-05,	0.000050311941072165,	5.11576487519719E-05,	5.20168875793913E-05,	5.28898619751813E-05,	5.37767792431662E-05,	5.46778496075451E-05,	5.55932862506315E-05,	5.65233053510262E-05,	5.74681261222299E-05,	5.84279708516977E-05,	5.94030649403411E-05,	6.03936369424815E-05,	6.13999186062584E-05,	6.24221449144989E-05,	6.34605541260503E-05,	6.45153878175823E-05,	6.55868909258619E-05,	6.66753117905063E-05,	0.000067780902197217,	6.89039174215017E-05,	7.00446162728868E-05,	7.12032611396252E-05,	7.23801180339048E-05,	7.35754566375622E-05,	7.47895503483056E-05,	7.60226763264517E-05,	7.72751155421826E-05,	7.85471528233249E-05,	7.98390769036587E-05,	8.11511804717592E-05,	8.24837602203763E-05,	8.38371168963573E-05,	8.52115553511178E-05,	8.66073845916645E-05,	8.80249178321766E-05,	8.94644725461492E-05,	9.09263705191044E-05,	9.24109379018754E-05,	9.39185052644681E-05,	0.000095449407650505,	9.70039846322578E-05,	9.85825803662726E-05,	0.000100185543649593,	0.000101813227976586,	0.000103465991596381,	0.000105144197570911,	0.00010684821383359,	0.000108578413248599,	0.00011033517367081,	0.000112118878006341,	0.000113929914273759,	0.000115768675665916,	0.000117635560612447,	0.000119530972842921,	0.000121455321450654,	0.000123409020957186,	0.000125392491377441,	0.000127406158285551,	0.000129450452881381,	0.000131525812057728,	0.000133632678468227,	0.000135771500595953,	0.000137942732822727,	0.000140146835499143,	0.000142384275015304,	0.000144655523872286,	0.000146961060754332,	0.000149301370601782,	0.000151676944684739,	0.00015408828067749,	0.000156535882733672,	0.000159020261562201,	0.000161541934503966,	0.000164101425609292,	0.000166699265716178,	0.000169335992529318,	0.000172012150699915,	0.000174728291906283,	0.000177484974935249,	0.000180282765764368,	0.000183122237644945,	0.000186003971185874,	0.000188928554438307,	0.000191896582981147,	0.000194908660007384,	0.000197965396411268,	0.000201067410876335,	0.000204215329964285,	0.000207409788204722,	0.000210651428185765,	0.000213940900645523,	0.000217278864564463,	0.000220665987258649,	0.000224102944473888,	0.000227590420480761,	0.000231129108170567,	0.000234719709152177,	0.000238362933849799,	0.000242059501601674,	0.0002458101407597,	0.000249615588789983,	0.000253476592374347,	0.000257393907512772,	0.000261368299626809,	0.000265400543663936,	0.000269491424202893,	0.000273641735559988,	0.000277852281896384,	0.000282123877326369,	0.000286457346026622,	0.000290853522346477,	0.000295313250919192,	0.000299837386774229,	0.000304426795450556,	0.000309082353110968,	0.000313804946657438,	0.000318595473847517,	0.000323454843411756,	0.000328383975172202,	0.000333383800161928,	0.000338455260745641,	0.00034359931074135,	0.000348816915543109,	0.000354109052244851,	0.000359476709765292,	0.000364920888973944,	0.000370442602818222,	0.000376042876451654,	0.00038172274736321,	0.000387483265507743,	0.000393325493437558,	0.000399250506435114,	0.000405259392646856,	0.000411353253218199,	0.000417533202429655,	0.000423800367834114,	0.000430155890395297,	0.000436600924627362,	0.000443136638735693,	0.000449764214758869,	0.000456484848711812,	0.000463299750730132,	0.000470210145215666,	0.000477217270983228,	0.000484322381408552,	0.000491526744577469,	0.000498831643436293,	0.000506238375943436,	0.000513748255222257,	0.000521362609715149,	0.000529082783338867,	0.000536910135641108,	0.000544846041958349,	0.000552891893574934,	0.000561049097883444,	0.000569319078546322,	0.000577703275658788,	0.000586203145913026,	0.000594820162763667,	0.000603555816594561,	0.000612411614886845,	0.00062138908238832,	0.00063048976128413,	0.000639715211368763,	0.000649067010219361,	0.00065854675337036,	0.000668156054489456,	0.000677896545554908,	0.000687769877034172,	0.000697777718063884,	0.000707921756631189,	0.00071820369975642,	0.00072862527367714,	0.000739188224033531,	0.000749894316055169,	0.00076074533474915,	0.0007717430850896,	0.000782889392208564,	0.000794186101588269,	0.000805635079254782,	0.000817238211973058,	0.000828997407443371,	0.00084091459449916,	0.000852991723306263,	0.000865230765563561,	0.000877633714705032,	0.000890202586103218,	0.000902939417274099,	0.0009158462680834,	0.000928925220954307,	0.000942178381076615,	0.0009556078766173,	0.000969215858932527,	0.000983004502781084,	0.000996976006539259,	0.00101113259241715,	0.00102547650667643,	0.00104001001984953,	0.00105473542696027,	0.001069655047746,	0.00108477122688108,	0.00110008633420192,	0.00111560276493338,	0.00113132293991672,	0.00114724930583892,	0.00116338433546348,	0.00117973052786277,	0.00119629040865165,	0.00121306653022278,	0.00123006147198317,	0.0012472778405924,	0.0012647182702021,	0.00128238542269709,	0.00130028198793781,	0.00131841068400435,	0.00133677425744185,	0.00135537548350743,	0.00137421716641854,	0.00139330213960278,	0.00141263326594919,	0.00143221343806102,	0.00145204557850994,	0.0014721326400917,	0.00149247760608324,	0.00151308349050133,	0.00153395333836258,	0.00155509022594491,	0.00157649726105057,	0.00159817758327046,	0.00162013436425002,	0.00164237080795651,	0.00166489015094774,	0.00168769566264223,	0.00171079064559081,	0.00173417843574972,	0.00175786240275497,	0.00178184595019835,	0.00180613251590465,	0.00183072557221047,	0.00185562862624432,	0.00188084522020818,	0.00190637893166049,	0.00193223337380049,	0.00195841219575399,	0.0019849190828605,	0.00201175775696178,	0.00203893197669173,	0.00206644553776768,	0.00209430227328305,	0.00212250605400132,	0.00215106078865145,	0.00217997042422451,	0.0022092389462718,	0.0022388703792042,	0.00226886878659284,	0.00229923827147118,	0.00232998297663828,	0.00236110708496347,	0.00239261481969225,	0.00242451044475351,	0.00245679826506797,	0.00248948262685796,	0.0025225679179584,	0.00255605856812903,	0.00258995904936793,	0.00262427387622615,	0.00265900760612368,	0.00269416483966656,	0.00272975022096517,	0.00276576843795372,	0.00280222422271095,	0.00283912235178189,	0.00287646764650084,	0.00291426497331547,	0.00295251924411203,	0.00299123541654164,	0.00303041849434774,	0.0030700735276945,	0.00311020561349647,	0.00315081989574909,	0.00319192156586036,	0.0032335158629835,	0.00327560807435058,	0.00331820353560719,	0.00336130763114804,	0.00340492579445355,	0.00344906350842735,	0.00349372630573473,	0.00353891976914195,	0.00358464953185652,	0.00363092127786827,	0.00367774074229128,	0.00372511371170671,	0.00377304602450638,	0.00382154357123712,	0.00387061229494601,	0.00392025819152625,	0.00397048731006383,	0.0040213057531849,	0.00407271967740383,	0.00412473529347193,	0.00417735886672689,	0.00423059671744268,	0.00428445522118029,	0.00433894080913882,	0.00439405996850729,	0.00444981924281692,	0.0045062252322939,	0.00456328459421273,	0.0046210040432499,	0.0046793903518381,	0.00473845035052083,	0.00479819092830735,	0.00485861903302802,	0.00491974167168996,	0.00498156591083304,	0.00504409887688609,	0.00510734775652342,	0.00517131979702154,	0.00523602230661605,	0.00530146265485868,	0.00536764827297455,	0.00543458665421944,	0.00550228535423714,	0.00557075199141691,	0.0056399942472509,	0.00571001986669149,	0.00578083665850875,	0.00585245249564759,	0.00592487531558499,	0.00599811312068693,	0.00607217397856526,	0.00614706602243421,	0.00622279745146677,	0.00629937653115072,	0.0063768115936443,	0.00645511103813165,	0.00653428333117764,	0.00661433700708251,	0.0066952806682358,	0.00677712298546998,	0.00685987269841335,	0.0069435386158425,	0.00702812961603404,	0.00711365464711571,	0.00720012272741679,	0.00728754294581775,	0.00737592446209908,	0.00746527650728938,	0.00755560838401247,	0.00764692946683366,	0.00773924920260504,	0.0078325771108098,	0.00792692278390542,	0.00802229588766591,	0.00811870616152278,	0.00821616341890494,	0.00831467754757729,	0.00841425850997812,	0.00851491634355515,	0.00861666116110021,	0.00871950315108253,	0.00882345257798055,	0.00892851978261223,	0.00903471518246384,	0.00914204927201708,	0.00925053262307461,	0.00936017588508384,	0.00947098978545897,	0.00958298512990122,	0.00969617280271725,	0.00981056376713558,	0.00992616906562113,	0.0100429998201877,	0.0101610672327084,	0.0102803825852239,	0.0104009572402487,	0.0105228026410748,	0.0106459303120735,	0.0107703518589945,	0.0108960789692629,	0.0110231234122734,	0.0111514970396824,	0.0112812117856973,	0.0114122796673629,	0.0115447127848455,	0.0116785233217145,	0.0118137235452203,	0.0119503258065702,	0.0120883425412008,	0.0122277862690479,	0.0123686695948127,	0.0125110052082257,	0.0126548058843067,	0.0128000844836219,	0.0129468539525378,	0.0130951273234711,	0.0132449177151363,	0.0133962383327888,	0.0135491024684646,	0.0137035235012169,	0.0138595148973486,	0.0140170902106406,	0.0141762630825776,	0.0143370472425685,	0.0144994565081638,	0.0146635047852686,	0.0148292060683515,	0.0149965744406495,	0.0151656240743686,	0.0153363692308798,	0.0155088242609114,	0.0156830036047359,	0.0158589217923531,	0.0160365934436687,	0.0162160332686673,	0.0163972560675819,	0.0165802767310575,	0.0167651102403106,	0.0169517716672832,	0.0171402761747921,	0.0173306390166728,	0.0175228755379188,	0.0177170011748147,	0.0179130314550648,	0.0181109819979159,	0.0183108685142748,	0.0185127068068197,	0.018716512770107,	0.0189223023906712,	0.0191300917471199,	0.0193398970102223,	0.0195517344429922,	0.0197656204007649,	0.0199815713312674,	0.0201996037746832,	0.0204197343637106,	0.0206419798236142,	0.0208663569722706,	0.0210928827202068,	0.0213215740706333,	0.0215524481194687,	0.0217855220553594,	0.0220208131596909,	0.0222583388065933,	0.0224981164629394,	0.022740163688335,	0.0229844981351034,	0.0232311375482617,	0.0234800997654897,	0.0237314027170919,	0.0239850644259516,	0.0242411030074776,	0.0244995366695431,	0.0247603837124169,	0.0250236625286866,	0.0252893916031741,	0.0255575895128429,	0.0258282749266972,	0.0261014666056734,	0.0263771834025222,	0.0266554442616837,	0.0269362682191525,	0.0272196744023358,	0.0275056820299021,	0.0277943104116207,	0.0280855789481935,	0.0283795071310771,	0.0286761145422961,	0.0289754208542475,	0.0292774458294962,	0.02958220932056,	0.0298897312696873,	0.0302000317086232,	0.0305131307583681,	0.030829048628925,	0.0311478056190385,	0.0314694221159227,	0.0317939185949809,	0.0321213156195135,	0.0324516338404174,	0.0327848939958744,	0.0331211169110297,	0.0334603234976601,	0.0338025347538319,	0.0341477717635477,	0.0344960556963837,	0.0348474078071159,	0.0352018494353347,	0.0355594020050502,	0.0359200870242855,	0.0362839260846594,	0.0366509408609581,	0.0370211531106952,	0.0373945846736613,	0.0377712574714615,	0.0381511935070419,	0.0385344148642042,	0.0389209437071093,	0.0393108022797686,	0.039704012905524,	0.0401005979865153,	0.0405005800031369,	0.0409039815134814,	0.0413108251527711,	0.0417211336327784,	0.0421349297412323,	0.0425522363412142,	0.0429730763705397,	0.0433974728411292,	0.0438254488383651,	0.0442570275204362,	0.0446922321176701,	0.0451310859318521,	0.0455736123355309,	0.0460198347713124,	0.0464697767511392,	0.0469234618555577,	0.0473809137329714,	0.0478421560988812,	0.0483072127351124,	0.0487761074890277,	0.049248864272727,	0.0497255070622338,	0.0502060598966669,	0.0506905468773998,	0.0511789921672045,	0.0516714199893828,	0.0521678546268829,	0.0526683204214018,	0.0531728417724739,	0.053681443136545,	0.0541941490260325,	0.0547109840083707,	0.0552319727050416,	0.0557571397905917,	0.056286509991634,	0.0568201080858348,	0.0573579589008867,	0.057900087313466,	0.058446518248176,	0.0589972766764743,	0.0595523876155863,	0.0601118761274031,	0.0606757673173638,	0.0612440863333238,	0.0618168583644066,	0.0623941086398415,	0.0629758624277847,	0.0635621450341262,	0.0641529818012804,	0.0647483981069615,	0.0653484193629434,	0.0659530710138036,	0.0665623785356521,	0.0671763674348436,	0.067795063246675,	0.0684184915340661,	0.0690466778862253,	0.0696796479172984,	0.0703174272650024,	0.0709600415892424,	0.0716075165707135,	0.0722598779094854,	0.0729171513235714,	0.0735793625474817,	0.0742465373307599,	0.0749187014365029,	0.0755958806398653,	0.0762781007265472,	0.0769653874912655,	0.0776577667362086,	0.0783552642694752,	0.0790579059034961,	0.0797657174534402,	0.0804787247356031,	0.0811969535657798,	0.0819204297576206,	0.0826491791209708,	0.0833832274601926,	0.0841226005724723,	0.084867324246109,	0.0856174242587874,	0.0863729263758346,	0.0871338563484587,	0.0879002399119721,	0.0886721027839973,	0.0894494706626561,	0.0902323692247424,	0.0910208241238775,	0.0918148609886497,	0.0926145054207366,	0.09341978299301,	0.0942307192476257,	0.0950473396940946,	0.0958696698073389,	0.0966977350257299,	0.0975315607491106,	0.0983711723368005,	0.0992165951055839,	0.100067854327682,	0.100924975228708,	0.101787982985604,	0.102656902724567,	0.103531759518949,	0.104412578387147,	0.105299384290479,	0.106192202131033,	0.10709105674951,	0.107995972923046,	0.108906975363019,	0.109824088712833,	0.110747337545699,	0.111676746362387,	0.112612339588966,	0.113554141574531,	0.114502176588911,	0.115456468820358,	0.116417042373224,	0.117383921265622,	0.118357129427068,	0.11933669069611,	0.120322628817938,	0.12131496744198,	0.122313730119485,	0.123318940301084,	0.12433062133434,	0.125348796461281,	0.126373488815917,	0.127404721421743,	0.128442517189225,	0.129486898913274,	0.130537889270696,	0.131595510817641,	0.132659785987022,	0.13373073708593,	0.134808386293031,	0.135892755655943,	0.136983867088606,	0.138081742368634,	0.139186403134652,	0.140297870883618,	0.141416166968133,	0.142541312593736,	0.143673328816184,	0.144812236538718,	0.145958056509316,	0.147110809317934,	0.148270515393729,	0.149437195002269,	0.150610868242737,	0.151791555045111,	0.15297927516734,	0.154174048192498,	0.155375893525935,	0.156584830392408,	0.157800877833203,	0.159024054703241,	0.160254379668176,	0.16149187120148,	0.16273654758151,	0.163988426888572,	0.16524752700197,	0.166513865597038,	0.16778746014217,	0.16906832789583,	0.170356485903557,	0.171651950994957,	0.172954739780683,	0.174264868649406,	0.175582353764777,	0.176907211062372,	0.178239456246639,	0.17957910478782,	0.180926171918879,	0.182280672632405,	0.183642621677521,	0.185012033556771,	0.186388922523006,	0.187773302576258,	0.189165187460605,	0.190564590661031,	0.191971525400274,	0.193386004635669,	0.19480804105598,	0.196237647078229,	0.197674834844515,	0.199119616218823,	0.200572002783834,	0.202032005837722,	0.203499636390943,	0.204974905163029,	0.206457822579359,	0.207948398767942,	0.209446643556179,	0.210952566467633,	0.212466176718784,	0.213987483215785,	0.215516494551211,	0.217053219000805,	0.218597664520221,	0.220149838741757,	0.221709748971094,	0.223277402184028,	0.224852805023192,	0.226435963794789,	0.228026884465309,	0.229625572658253,	0.231232033650854,	0.232846272370789,	0.2344682933929,	0.236098100935908,	0.237735698859125,	0.239381090659172,	0.24103427946669,	0.242695268043055,	0.24436405877709,	0.246040653681783,	0.247725054391003,	0.249417262156215,	0.2511172778432,	0.252825101928775,	0.254540734497518,	0.256264175238493,	0.257995423441975,	0.259734477996188,	0.261481337384033,	0.26323599967983,	0.264998462546064,	0.266768723230125,	0.268546778561066,	0.27033262494636,	0.272126258368657,	0.273927674382559,	0.27573686811139,	0.277553834243978,	0.279378567031445,	0.281211060283995,	0.283051307367724,	0.284899301201425,	0.286755034253405,	0.288618498538317,	0.290489685613988,	0.292368586578271,	0.294255192065891,	0.296149492245315,	0.298051476815623,	0.299961135003395,	0.301878455559603,	0.303803426756522,	0.305736036384645,	0.307676271749617,	0.309624119669176,	0.31157956647011,	0.313542597985224,	0.315513199550327,	0.317491356001225,	0.319477051670731,	0.321470270385694,	0.323470995464038,	0.325479209711819,	0.327494895420294,	0.329518034363011,	0.331548607792916,	0.33358659643947,	0.335631980505791,	0.337684739665811,	0.33974485306145,	0.341812299299807,	0.343887056450378,	0.345969102042283,	0.34805841306152,	0.350154965948233,	0.35225873659401,	0.354369700339191,	0.356487831970202,	0.358613105716917,	0.360745495250027,	0.36288497367845,	0.365031513546747,	0.367185086832573,	0.369345664944146,	0.371513218717744,	0.373687718415219,	0.375869133721548,	0.378057433742397,	0.380252587001719,	0.382454561439375,	0.384663324408781,	0.386878842674584,	0.389101082410366,	0.391330009196369,	0.393565588017259,	0.39580778325991,	0.39805655871122,	0.400311877555958,	0.402573702374638,	0.404841995141425,	0.40711671722207,	0.40939782937188,	0.411685291733712,	0.41397906383601,	0.416279104590862,	0.4185853722921,	0.420897824613425,	0.423216418606571,	0.425541110699504,	0.427871856694647,	0.430208611767149,	0.432551330463183,	0.434899966698286,	0.437254473755727,	0.439614804284915,	0.441980910299847,	0.44435274317759,	0.446730253656797,	0.449113391836269,	0.451502107173549,	0.453896348483558,	0.456296063937267,	0.458701201060413,	0.461111706732251,	0.463527527184348,	0.46594860799942,	0.468374894110203,	0.470806329798376,	0.473242858693515,	0.4756844237721,	0.478130967356552,	0.480582431114329,	0.483038756057048,	0.485499882539666,	0.487965750259696,	0.490436298256471,	0.492911464910451,	0.495391187942579,	0.497875404413679,	0.5003640507239,	0.50285706261221,	0.505354375155934,	0.50785592277034,	0.510361639208272,	0.51287145755983,	0.515385310252102,	0.51790312904894,	0.520424845050786,	0.522950388694553,	0.525479689753545,	0.528012677337437,	0.530549279892298,	0.533089425200671,	0.535633040381699,	0.538180051891302,	0.540730385522413,	0.543283966405254,	0.545840719007678,	0.548400567135552,	0.550963433933196,	0.553529241883882,	0.556097912810377,	0.558669367875547,	0.561243527583008,	0.56382031177784,	0.566399639647349,	0.568981429721888,	0.571565599875729,	0.574152067327997,	0.576740748643652,	0.579331559734535,	0.581924415860465,	0.584519231630394,	0.587115921003621,	0.589714397291061,	0.59231457315657,	0.594916360618332,	0.597519671050304,	0.60012441518371,	0.602730503108608,	0.605337844275504,	0.60794634749703,	0.610555920949678,	0.6131664721756,	0.615777908084461,	0.61839013495535,	0.621003058438764,	0.623616583558632,	0.626230614714419,	0.628845055683278,	0.63145980962227,	0.634074779070636,	0.636689865952143,	0.639304971577479,	0.641919996646717,	0.644534841251838,	0.647149404879315,	0.649763586412763,	0.652377284135646,	0.654990395734045,	0.657602818299501,	0.660214448331901,	0.662825181742442,	0.665434913856656,	0.668043539417485,	0.670650952588437,	0.673257046956793,	0.67586171553688,	0.678464850773409,	0.681066344544875,	0.683666088167021,	0.686263972396363,	0.688859887433782,	0.691453722928181,	0.694045367980197,	0.696634711145987,	0.69922164044107,	0.701806043344242,	0.70438780680154,	0.706966817230286,	0.709542960523185,	0.71211612205249,	0.714686186674227,	0.717253038732493,	0.719816562063806,	0.72237664000153,	0.724933155380354,	0.727485990540844,	0.730035027334051,	0.732580147126183,	0.735121230803354,	0.737658158776372,	0.740190810985617,	0.742719066905962,	0.745242805551771,	0.74776190548195,	0.750276244805068,	0.752785701184542,	0.755290151843877,	0.757789473571978,	0.760283542728519,	0.762772235249377,	0.765255426652129,	0.767732992041609,	0.770204806115527,	0.772670743170159,	0.775130677106082,	0.777584481433988,	0.780032029280549,	0.782473193394345,	0.78490784615186,	0.787335859563527,	0.789757105279847,	0.792171454597556,	0.794578778465865,	0.796978947492749,	0.799371831951301,	0.801757301786149,	0.804135226619924,	0.806505475759797,	0.808867918204066,	0.811222422648807,	0.813568857494582,	0.815907090853207,	0.818236990554572,	0.820558424153526,	0.822871258936813,	0.825175361930071,	0.827470599904878,	0.82975683938587,	0.832033946657894,	0.834301787773237,	0.836560228558896,	0.838809134623909,	0.84104837136674,	0.843277803982719,	0.84549729747153,	0.847706716644759,	0.849905926133495,	0.852094790395976,	0.854273173725295,	0.856440940257157,	0.858597953977681,	0.86074407873126,	0.862879178228466,	0.86500311605401,	0.867115755674751,	0.869216960447743,	0.871306593628352,	0.873384518378399,	0.875450597774365,	0.877504694815637,	0.879546672432805,	0.881576393495997,	0.883593720823271,	0.885598517189042,	0.887590645332558,	0.889569967966422,	0.891536347785154,	0.893489647473797,	0.895429729716564,	0.897356457205532,	0.899269692649371,	0.901169298782117,	0.903055138371986,	0.904927074230224,	0.906784969219998,	0.908628686265328,	0.910458088360051,	0.912273038576826,	0.914073400076177,	0.915859036115567,	0.917629810058509,	0.919385585383714,	0.921126225694272,	0.922851594726863,	0.924561556361005,	0.926255974628331,	0.927934713721901,	0.929597638005536,	0.931244612023194,	0.932875500508365,	0.934490168393501,	0.93608848081947,	0.937670303145038,	0.939235500956379,	0.940783940076608,	0.942315486575341,	0.943830006778279,	0.945327367276816,	0.946807434937663,	0.948270076912507,	0.949715160647679,	0.951142553893849,	0.952552124715738,	0.95394374150185,	0.955317272974222,	0.95667258819819,	0.958009556592175,	0.95932804793748,	0.960627932388105,	0.961909080480575,	0.963171363143783,	0.964414651708844,	0.965638817918959,	0.966843733939293,	0.968029272366863,	0.969195306240431,	0.970341709050411,	0.971468354748781,	0.972575117758997,	0.973661872985924,	0.974728495825759,	0.975774862175968,	0.976800848445224,	0.97780633156334,	0.978791188991215,	0.979755298730772,	0.980698539334902,	0.981620789917399,	0.982521930162901,	0.983401840336825,	0.984260401295296,	0.985097494495073,	0.985913002003474,	0.986706806508282,	0.98747879132766,	0.988228840420045,	0.988956838394037,	0.989662670518282,	0.990346222731333,	0.991007381651514,	0.991646034586757,	0.992262069544434,	0.992855375241171,	0.993425841112648,	0.993973357323379,	0.994497814776481,	0.994999105123414,	0.995477120773714,	0.995931754904697,	0.996362901471143,	0.996770455214961,	0.997154311674832,	0.997514367195816,	0.997850518938954,	0.998162664890825,	0.998450703873085,	0.998714535551981,	0.99895406044783,	0.999169179944469,	0.999359796298678,	0.999525812649569,	0.999667133027942,	0.999783662365609,	0.999875306504681,	0.999941972206822,	0.999983567162467,	1,	0.999991180294894,	0.999957018578816,	0.999897426348684,	0.999812316075695,	0.999701601214296,	0.999565196211123,	0.999403016513891,	0.999214978580245,	0.999000999886554,	0.99876099893667,	0.998494895270633,	0.998202609473328,	0.997884063183094,	0.997539179100281,	0.997167880995755,	0.996770093719354,	0.996345743208287,	0.99589475649548,	0.995417061717869,	0.994912588124632,	0.994381266085369,	0.993823027098219,	0.993237803797929,	0.992625529963844,	0.99198614052786,	0.991319571582292,	0.990625760387701,	0.989904645380641,	0.989156166181349,	0.988380263601375,	0.987576879651132,	0.986745957547393,	0.985887441720714,	0.985001277822787,	0.984087412733729,	0.983145794569293,	0.982176372688015,	0.981179097698286,	0.980153921465349,	0.979100797118227,	0.978019679056571,	0.976910522957438,	0.975773285781986,	0.974607925782102,	0.973414402506939,	0.972192676809387,	0.970942710852455,	0.969664468115578,	0.968357913400842,	0.967023012839123,	0.965659733896151,	0.964268045378482,	0.962847917439393,	0.961399321584687,	0.959922230678414,	0.958416618948504,	0.956882461992318,	0.955319736782104,	0.953728421670368,	0.952108496395153,	0.950459942085232,	0.948782741265207,	0.947076877860516,	0.945342337202347,	0.943579106032466,	0.941787172507938,	0.939966526205768,	0.938117158127435,	0.93623906070334,	0.934332227797149,	0.932396654710047,	0.930432338184886,	0.928439276410242,	0.92641746902437,	0.924366917119058,	0.922287623243384,	0.920179591407369,	0.918042827085533,	0.915877337220341,	0.913683130225559,	0.911460215989495,	0.909208605878142,	0.906928312738217,	0.904619350900095,	0.902281736180636,	0.899915485885907,	0.897520618813802,	0.895097155256544,	0.892645117003098,	0.890164527341456,	0.88765541106083,	0.885117794453726,	0.882551705317915,	0.879957172958294,	0.87733422818863,	0.874682903333205,	0.872003232228343,	0.869295250223825,	0.866558994184195,	0.863794502489958,	0.861001815038657,	0.858180973245845,	0.855332020045938,	0.852454999892963,	0.849549958761178,	0.846616944145598,	0.843656005062389,	0.840667192049155,	0.837650557165112,	0.834606153991141,	0.831534037629734,	0.828434264704816,	0.825306893361454,	0.822151983265453,	0.818969595602834,	0.815759793079191,	0.81252263991894,	0.809258201864443,	0.805966546175021,	0.802647741625844,	0.79930185850671,	0.795928968620705,	0.792529145282737,	0.789102463317968,	0.785648999060112,	0.78216883034963,	0.778662036531791,	0.775128698454634,	0.771568898466791,	0.767982720415213,	0.764370249642759,	0.760731572985681,	0.757066778770981,	0.753375956813659,	0.749659198413833,	0.745916596353751,	0.742148244894674,	0.738354239773652,	0.734534678200176,	0.730689658852711,	0.726819281875116,	0.722923648872942,	0.719002862909616,	0.715057028502505,	0.711086251618864,	0.707090639671665,	0.703070301515315,	0.699025347441249,	0.69495588917341,	0.690862039863618,	0.686743914086811,	0.682601627836183,	0.678435298518198,	0.674245044947489,	0.670030987341645,	0.665793247315886,	0.66153194787761,	0.657247213420845,	0.652939169720569,	0.648607943926929,	0.644253664559339,	0.639876461500468,	0.635476465990119,	0.631053810618986,	0.626608629322309,	0.622141057373413,	0.617651231377135,	0.613139289263144,	0.608605370279145,	0.60404961498398,	0.599472165240614,	0.594873164209016,	0.590252756338931,	0.585611087362538,	0.580948304287013,	0.57626455538697,	0.571559990196811,	0.566834759502956,	0.562089015335978,	0.557322910962627,	0.552536600877755,	0.547730240796131,	0.542903987644161,	0.538057999551496,	0.533192435842548,	0.528307457027897,	0.523403224795603,	0.518479902002412,	0.51353765266487,	0.508576641950335,	0.503597036167888,	0.498599002759152,	0.49358271028901,	0.488548328436234,	0.483496027984008,	0.478425980810369,	0.473338359878543,	0.468233339227199,	0.463111093960601,	0.457971800238673,	0.452815635266978,	0.4476427772866,	0.442453405563936,	0.437247700380411,	0.432025843022094,	0.426788015769233,	0.421534401885701,	0.416265185608364,	0.410980552136356,	0.405680687620279,	0.400365779151314,	0.395036014750256,	0.389691583356468,	0.384332674816749,	0.378959479874134,	0.373572190156608,	0.368170998165745,	0.362756097265276,	0.357327681669576,	0.351885946432078,	0.346431087433618,	0.340963301370707,	0.335482785743727,	0.329989738845063,	0.324484359747162,	0.318966848290531,	0.313437405071654,	0.30789623143086,	0.302343529440115,	0.29677950189075,	0.291204352281134,	0.285618284804279,	0.280021504335385,	0.274414216419324,	0.268796627258074,	0.26316894369808,	0.257531373217575,	0.251884123913835,	0.246227404490381,	0.240561424244135,	0.234886393052513,	0.229202521360474,	0.223510020167523,	0.217809101014651,	0.212099975971246,	0.206382857621942,	0.200657959053429,	0.194925493841224,	0.189185676036388,	0.183438720152209,	0.177684841150843,	0.171924254429914,	0.166157175809076,	0.160383821516538,	0.154604408175555,	0.148819152790879,	0.143028272735185,	0.137231985735452,	0.131430509859328,	0.125624063501451,	0.119812865369748,	0.113997134471711,	0.108177090100632,	0.102352951821829,	0.0965249394588403,	0.0906932730795914,	0.0848581729825504,	0.0790198596828564,	0.0731785538984303,	0.067334476536068,	0.0614878486775168,	0.055638891565536,	0.0497878265899439,	0.0439348752736512,	0.0380802592586829,	0.0322242002921896,	0.0263669202124499,	0.0205086409348638,	0.0146495844379403,	0.00878997274927886,	0.00293002793154675,	-0.00293002793154699,	-0.0087899727492791,	-0.0146495844379405,	-0.020508640934864,	-0.0263669202124501,	-0.0322242002921899,	-0.0380802592586831,	-0.0439348752736515,	-0.0497878265899442,	-0.0556388915655363,	-0.061487848677517,	-0.0673344765360683,	-0.0731785538984306,	-0.0790198596828567,	-0.0848581729825507,	-0.0906932730795916,	-0.0965249394588406,	-0.10235295182183,	-0.108177090100632,	-0.113997134471711,	-0.119812865369749,	-0.125624063501451,	-0.131430509859328,	-0.137231985735453,	-0.143028272735185,	-0.148819152790879,	-0.154604408175555,	-0.160383821516538,	-0.166157175809076,	-0.171924254429914,	-0.177684841150843,	-0.183438720152209,	-0.189185676036388,	-0.194925493841224,	-0.200657959053429,	-0.206382857621942,	-0.212099975971246,	-0.217809101014652,	-0.223510020167523,	-0.229202521360475,	-0.234886393052513,	-0.240561424244136,	-0.246227404490382,	-0.251884123913835,	-0.257531373217575,	-0.26316894369808,	-0.268796627258074,	-0.274414216419325,	-0.280021504335385,	-0.285618284804279,	-0.291204352281134,	-0.29677950189075,	-0.302343529440115,	-0.307896231430861,	-0.313437405071654,	-0.318966848290531,	-0.324484359747163,	-0.329989738845063,	-0.335482785743727,	-0.340963301370707,	-0.346431087433618,	-0.351885946432078,	-0.357327681669576,	-0.362756097265276,	-0.368170998165745,	-0.373572190156608,	-0.378959479874134,	-0.384332674816749,	-0.389691583356468,	-0.395036014750256,	-0.400365779151314,	-0.405680687620279,	-0.410980552136356,	-0.416265185608364,	-0.421534401885701,	-0.426788015769233,	-0.432025843022094,	-0.437247700380411,	-0.442453405563936,	-0.447642777286599,	-0.452815635266978,	-0.457971800238673,	-0.463111093960601,	-0.468233339227199,	-0.473338359878543,	-0.478425980810369,	-0.483496027984008,	-0.488548328436234,	-0.49358271028901,	-0.498599002759152,	-0.503597036167888,	-0.508576641950335,	-0.51353765266487,	-0.518479902002412,	-0.523403224795603,	-0.528307457027897,	-0.533192435842548,	-0.538057999551496,	-0.542903987644161,	-0.547730240796131,	-0.552536600877754,	-0.557322910962627,	-0.562089015335978,	-0.566834759502956,	-0.571559990196811,	-0.57626455538697,	-0.580948304287012,	-0.585611087362538,	-0.59025275633893,	-0.594873164209016,	-0.599472165240613,	-0.604049614983979,	-0.608605370279145,	-0.613139289263144,	-0.617651231377135,	-0.622141057373413,	-0.626608629322308,	-0.631053810618985,	-0.635476465990118,	-0.639876461500468,	-0.644253664559338,	-0.648607943926929,	-0.652939169720569,	-0.657247213420845,	-0.66153194787761,	-0.665793247315885,	-0.670030987341645,	-0.674245044947488,	-0.678435298518198,	-0.682601627836183,	-0.686743914086811,	-0.690862039863617,	-0.69495588917341,	-0.699025347441248,	-0.703070301515315,	-0.707090639671665,	-0.711086251618863,	-0.715057028502505,	-0.719002862909616,	-0.722923648872942,	-0.726819281875115,	-0.730689658852711,	-0.734534678200176,	-0.738354239773651,	-0.742148244894673,	-0.74591659635375,	-0.749659198413833,	-0.753375956813658,	-0.75706677877098,	-0.76073157298568,	-0.764370249642759,	-0.767982720415212,	-0.771568898466791,	-0.775128698454633,	-0.778662036531791,	-0.782168830349629,	-0.785648999060112,	-0.789102463317967,	-0.792529145282736,	-0.795928968620704,	-0.79930185850671,	-0.802647741625843,	-0.80596654617502,	-0.809258201864442,	-0.81252263991894,	-0.815759793079191,	-0.818969595602833,	-0.822151983265453,	-0.825306893361453,	-0.828434264704816,	-0.831534037629734,	-0.834606153991141,	-0.837650557165111,	-0.840667192049155,	-0.843656005062389,	-0.846616944145598,	-0.849549958761177,	-0.852454999892962,	-0.855332020045938,	-0.858180973245844,	-0.861001815038656,	-0.863794502489957,	-0.866558994184194,	-0.869295250223824,	-0.872003232228343,	-0.874682903333205,	-0.877334228188629,	-0.879957172958293,	-0.882551705317915,	-0.885117794453725,	-0.887655411060829,	-0.890164527341455,	-0.892645117003097,	-0.895097155256544,	-0.897520618813801,	-0.899915485885907,	-0.902281736180636,	-0.904619350900095,	-0.906928312738217,	-0.909208605878141,	-0.911460215989494,	-0.913683130225559,	-0.915877337220341,	-0.918042827085532,	-0.920179591407369,	-0.922287623243383,	-0.924366917119057,	-0.926417469024369,	-0.928439276410241,	-0.930432338184885,	-0.932396654710046,	-0.934332227797149,	-0.936239060703339,	-0.938117158127435,	-0.939966526205767,	-0.941787172507938,	-0.943579106032465,	-0.945342337202347,	-0.947076877860515,	-0.948782741265206,	-0.950459942085231,	-0.952108496395152,	-0.953728421670367,	-0.955319736782104,	-0.956882461992317,	-0.958416618948503,	-0.959922230678413,	-0.961399321584687,	-0.962847917439393,	-0.964268045378481,	-0.96565973389615,	-0.967023012839122,	-0.968357913400841,	-0.969664468115577,	-0.970942710852454,	-0.972192676809386,	-0.973414402506938,	-0.974607925782101,	-0.975773285781985,	-0.976910522957437,	-0.97801967905657,	-0.979100797118226,	-0.980153921465348,	-0.981179097698285,	-0.982176372688014,	-0.983145794569292,	-0.984087412733728,	-0.985001277822786,	-0.985887441720713,	-0.986745957547392,	-0.987576879651131,	-0.988380263601374,	-0.989156166181349,	-0.98990464538064,	-0.9906257603877,	-0.991319571582291,	-0.991986140527859,	-0.992625529963844,	-0.993237803797928,	-0.993823027098219,	-0.994381266085368,	-0.994912588124632,	-0.995417061717869,	-0.995894756495479,	-0.996345743208286,	-0.996770093719353,	-0.997167880995754,	-0.99753917910028,	-0.997884063183093,	-0.998202609473327,	-0.998494895270632,	-0.998760998936669,	-0.999000999886553,	-0.999214978580244,	-0.999403016513891,	-0.999565196211122,	-0.999701601214295,	-0.999812316075695,	-0.999897426348684,	-0.999957018578815,	-0.999991180294893,	-0.999999999999999,	-0.999983567162466,	-0.999941972206821,	-0.99987530650468,	-0.999783662365608,	-0.999667133027941,	-0.999525812649568,	-0.999359796298677,	-0.999169179944468,	-0.998954060447829,	-0.998714535551981,	-0.998450703873084,	-0.998162664890824,	-0.997850518938954,	-0.997514367195816,	-0.997154311674831,	-0.99677045521496,	-0.996362901471142,	-0.995931754904696,	-0.995477120773713,	-0.994999105123413,	-0.99449781477648,	-0.993973357323378,	-0.993425841112647,	-0.99285537524117,	-0.992262069544433,	-0.991646034586756,	-0.991007381651513,	-0.990346222731333,	-0.989662670518281,	-0.988956838394037,	-0.988228840420044,	-0.987478791327659,	-0.986706806508281,	-0.985913002003473,	-0.985097494495073,	-0.984260401295295,	-0.983401840336824,	-0.9825219301629,	-0.981620789917398,	-0.980698539334901,	-0.979755298730771,	-0.978791188991214,	-0.977806331563339,	-0.976800848445223,	-0.975774862175968,	-0.974728495825758,	-0.973661872985923,	-0.972575117758996,	-0.97146835474878,	-0.97034170905041,	-0.96919530624043,	-0.968029272366861,	-0.966843733939292,	-0.965638817918958,	-0.964414651708843,	-0.963171363143782,	-0.961909080480574,	-0.960627932388104,	-0.959328047937479,	-0.958009556592174,	-0.956672588198189,	-0.955317272974221,	-0.953943741501849,	-0.952552124715737,	-0.951142553893848,	-0.949715160647678,	-0.948270076912506,	-0.946807434937662,	-0.945327367276815,	-0.943830006778278,	-0.94231548657534,	-0.940783940076607,	-0.939235500956378,	-0.937670303145037,	-0.936088480819469,	-0.9344901683935,	-0.932875500508364,	-0.931244612023193,	-0.929597638005535,	-0.9279347137219,	-0.92625597462833,	-0.924561556361004,	-0.922851594726862,	-0.921126225694271,	-0.919385585383713,	-0.917629810058508,	-0.915859036115566,	-0.914073400076177,	-0.912273038576825,	-0.91045808836005,	-0.908628686265327,	-0.906784969219997,	-0.904927074230223,	-0.903055138371985,	-0.901169298782116,	-0.89926969264937,	-0.897356457205531,	-0.895429729716563,	-0.893489647473796,	-0.891536347785153,	-0.889569967966421,	-0.887590645332557,	-0.885598517189041,	-0.88359372082327,	-0.881576393495996,	-0.879546672432804,	-0.877504694815636,	-0.875450597774364,	-0.873384518378398,	-0.871306593628351,	-0.869216960447742,	-0.86711575567475,	-0.865003116054009,	-0.862879178228465,	-0.860744078731259,	-0.85859795397768,	-0.856440940257156,	-0.854273173725294,	-0.852094790395975,	-0.849905926133494,	-0.847706716644758,	-0.845497297471529,	-0.843277803982718,	-0.841048371366739,	-0.838809134623908,	-0.836560228558895,	-0.834301787773236,	-0.832033946657893,	-0.829756839385869,	-0.827470599904877,	-0.82517536193007,	-0.822871258936812,	-0.820558424153525,	-0.818236990554572,	-0.815907090853206,	-0.813568857494581,	-0.811222422648806,	-0.808867918204065,	-0.806505475759796,	-0.804135226619923,	-0.801757301786148,	-0.7993718319513,	-0.796978947492748,	-0.794578778465864,	-0.792171454597555,	-0.789757105279846,	-0.787335859563526,	-0.784907846151859,	-0.782473193394344,	-0.780032029280548,	-0.777584481433987,	-0.775130677106081,	-0.772670743170158,	-0.770204806115526,	-0.767732992041607,	-0.765255426652128,	-0.762772235249376,	-0.760283542728518,	-0.757789473571977,	-0.755290151843876,	-0.752785701184541,	-0.750276244805067,	-0.747761905481948,	-0.74524280555177,	-0.742719066905961,	-0.740190810985616,	-0.737658158776371,	-0.735121230803352,	-0.732580147126182,	-0.730035027334049,	-0.727485990540843,	-0.724933155380353,	-0.722376640001529,	-0.719816562063805,	-0.717253038732492,	-0.714686186674226,	-0.712116122052489,	-0.709542960523184,	-0.706966817230285,	-0.704387806801539,	-0.701806043344241,	-0.699221640441069,	-0.696634711145985,	-0.694045367980196,	-0.69145372292818,	-0.688859887433781,	-0.686263972396361,	-0.68366608816702,	-0.681066344544874,	-0.678464850773408,	-0.675861715536879,	-0.673257046956792,	-0.670650952588436,	-0.668043539417484,	-0.665434913856655,	-0.662825181742441,	-0.660214448331899,	-0.6576028182995,	-0.654990395734044,	-0.652377284135644,	-0.649763586412762,	-0.647149404879314,	-0.644534841251836,	-0.641919996646716,	-0.639304971577478,	-0.636689865952142,	-0.634074779070635,	-0.631459809622268,	-0.628845055683277,	-0.626230614714417,	-0.62361658355863,	-0.621003058438763,	-0.618390134955349,	-0.615777908084459,	-0.613166472175599,	-0.610555920949677,	-0.607946347497028,	-0.605337844275503,	-0.602730503108607,	-0.600124415183709,	-0.597519671050303,	-0.594916360618331,	-0.592314573156568,	-0.58971439729106,	-0.58711592100362,	-0.584519231630393,	-0.581924415860464,	-0.579331559734534,	-0.576740748643651,	-0.574152067327996,	-0.571565599875728,	-0.568981429721887,	-0.566399639647348,	-0.563820311777838,	-0.561243527583006,	-0.558669367875545,	-0.556097912810376,	-0.553529241883881,	-0.550963433933195,	-0.54840056713555,	-0.545840719007677,	-0.543283966405253,	-0.540730385522411,	-0.538180051891301,	-0.535633040381698,	-0.53308942520067,	-0.530549279892297,	-0.528012677337436,	-0.525479689753544,	-0.522950388694552,	-0.520424845050785,	-0.517903129048938,	-0.515385310252101,	-0.512871457559829,	-0.510361639208271,	-0.507855922770339,	-0.505354375155933,	-0.502857062612209,	-0.500364050723899,	-0.497875404413678,	-0.495391187942578,	-0.49291146491045,	-0.49043629825647,	-0.487965750259695,	-0.485499882539665,	-0.483038756057047,	-0.480582431114328,	-0.478130967356551,	-0.475684423772099,	-0.473242858693514,	-0.470806329798375,	-0.468374894110202,	-0.465948607999418,	-0.463527527184347,	-0.46111170673225,	-0.458701201060412,	-0.456296063937266,	-0.453896348483557,	-0.451502107173548,	-0.449113391836268,	-0.446730253656796,	-0.444352743177588,	-0.441980910299846,	-0.439614804284913,	-0.437254473755725,	-0.434899966698285,	-0.432551330463182,	-0.430208611767148,	-0.427871856694646,	-0.425541110699503,	-0.42321641860657,	-0.420897824613423,	-0.418585372292099,	-0.416279104590861,	-0.413979063836009,	-0.411685291733711,	-0.409397829371878,	-0.407116717222069,	-0.404841995141424,	-0.402573702374637,	-0.400311877555957,	-0.398056558711219,	-0.395807783259908,	-0.393565588017258,	-0.391330009196368,	-0.389101082410365,	-0.386878842674583,	-0.38466332440878,	-0.382454561439374,	-0.380252587001718,	-0.378057433742396,	-0.375869133721547,	-0.373687718415218,	-0.371513218717743,	-0.369345664944145,	-0.367185086832571,	-0.365031513546746,	-0.362884973678449,	-0.360745495250026,	-0.358613105716915,	-0.356487831970201,	-0.35436970033919,	-0.352258736594009,	-0.350154965948232,	-0.348058413061519,	-0.345969102042282,	-0.343887056450377,	-0.341812299299806,	-0.339744853061449,	-0.33768473966581,	-0.33563198050579,	-0.333586596439469,	-0.331548607792915,	-0.32951803436301,	-0.327494895420292,	-0.325479209711818,	-0.323470995464037,	-0.321470270385693,	-0.319477051670729,	-0.317491356001223,	-0.315513199550326,	-0.313542597985223,	-0.311579566470108,	-0.309624119669175,	-0.307676271749616,	-0.305736036384644,	-0.30380342675652,	-0.301878455559602,	-0.299961135003394,	-0.298051476815622,	-0.296149492245314,	-0.29425519206589,	-0.29236858657827,	-0.290489685613987,	-0.288618498538316,	-0.286755034253404,	-0.284899301201423,	-0.283051307367723,	-0.281211060283994,	-0.279378567031443,	-0.277553834243977,	-0.275736868111388,	-0.273927674382557,	-0.272126258368656,	-0.270332624946358,	-0.268546778561065,	-0.266768723230123,	-0.264998462546062,	-0.263235999679829,	-0.261481337384032,	-0.259734477996187,	-0.257995423441974,	-0.256264175238492,	-0.254540734497517,	-0.252825101928774,	-0.251117277843198,	-0.249417262156214,	-0.247725054391002,	-0.246040653681782,	-0.244364058777089,	-0.242695268043054,	-0.241034279466689,	-0.239381090659171,	-0.237735698859124,	-0.236098100935906,	-0.234468293392899,	-0.232846272370788,	-0.231232033650853,	-0.229625572658252,	-0.228026884465308,	-0.226435963794788,	-0.224852805023191,	-0.223277402184027,	-0.221709748971093,	-0.220149838741756,	-0.21859766452022,	-0.217053219000804,	-0.21551649455121,	-0.213987483215784,	-0.212466176718782,	-0.210952566467632,	-0.209446643556178,	-0.207948398767941,	-0.206457822579358,	-0.204974905163028,	-0.203499636390942,	-0.202032005837721,	-0.200572002783833,	-0.199119616218822,	-0.197674834844514,	-0.196237647078228,	-0.194808041055979,	-0.193386004635668,	-0.191971525400273,	-0.19056459066103,	-0.189165187460604,	-0.187773302576257,	-0.186388922523005,	-0.18501203355677,	-0.18364262167752,	-0.182280672632404,	-0.180926171918877,	-0.179579104787819,	-0.178239456246638,	-0.176907211062371,	-0.175582353764776,	-0.174264868649405,	-0.172954739780682,	-0.171651950994956,	-0.170356485903556,	-0.169068327895829,	-0.167787460142169,	-0.166513865597037,	-0.165247527001969,	-0.163988426888571,	-0.162736547581509,	-0.161491871201479,	-0.160254379668175,	-0.15902405470324,	-0.157800877833202,	-0.156584830392407,	-0.155375893525934,	-0.154174048192497,	-0.152979275167339,	-0.15179155504511,	-0.150610868242736,	-0.149437195002268,	-0.148270515393728,	-0.147110809317933,	-0.145958056509315,	-0.144812236538717,	-0.143673328816183,	-0.142541312593735,	-0.141416166968132,	-0.140297870883617,	-0.139186403134651,	-0.138081742368633,	-0.136983867088605,	-0.135892755655942,	-0.13480838629303,	-0.133730737085929,	-0.132659785987021,	-0.13159551081764,	-0.130537889270695,	-0.129486898913273,	-0.128442517189225,	-0.127404721421742,	-0.126373488815916,	-0.12534879646128,	-0.124330621334339,	-0.123318940301083,	-0.122313730119484,	-0.121314967441979,	-0.120322628817937,	-0.119336690696109,	-0.118357129427067,	-0.117383921265621,	-0.116417042373223,	-0.115456468820357,	-0.11450217658891,	-0.11355414157453,	-0.112612339588965,	-0.111676746362386,	-0.110747337545698,	-0.109824088712832,	-0.108906975363018,	-0.107995972923046,	-0.107091056749509,	-0.106192202131032,	-0.105299384290478,	-0.104412578387146,	-0.103531759518948,	-0.102656902724566,	-0.101787982985603,	-0.100924975228707,	-0.100067854327681,	-0.0992165951055829,	-0.0983711723367996,	-0.0975315607491097,	-0.096697735025729,	-0.0958696698073379,	-0.0950473396940937,	-0.0942307192476248,	-0.0934197829930091,	-0.0926145054207357,	-0.0918148609886488,	-0.0910208241238766,	-0.0902323692247414,	-0.0894494706626552,	-0.0886721027839964,	-0.0879002399119712,	-0.0871338563484578,	-0.0863729263758336,	-0.0856174242587865,	-0.0848673242461081,	-0.0841226005724714,	-0.0833832274601917,	-0.0826491791209699,	-0.0819204297576197,	-0.0811969535657789,	-0.0804787247356022,	-0.0797657174534393,	-0.0790579059034952,	-0.0783552642694743,	-0.0776577667362077,	-0.0769653874912646,	-0.0762781007265464,	-0.0755958806398644,	-0.074918701436502,	-0.074246537330759,	-0.0735793625474809,	-0.0729171513235705,	-0.0722598779094845,	-0.0716075165707127,	-0.0709600415892416,	-0.0703174272650015,	-0.0696796479172975,	-0.0690466778862244,	-0.0684184915340653,	-0.0677950632466741,	-0.0671763674348427,	-0.0665623785356512,	-0.0659530710138028,	-0.0653484193629425,	-0.0647483981069607,	-0.0641529818012795,	-0.0635621450341253,	-0.0629758624277838,	-0.0623941086398406,	-0.0618168583644058,	-0.0612440863333229,	-0.060675767317363,	-0.0601118761274022,	-0.0595523876155855,	-0.0589972766764735,	-0.0584465182481752,	-0.0579000873134652,	-0.0573579589008858,	-0.0568201080858339,	-0.0562865099916332,	-0.0557571397905909,	-0.0552319727050408,	-0.0547109840083699,	-0.0541941490260317,	-0.0536814431365442,	-0.053172841772473,	-0.052668320421401,	-0.0521678546268821,	-0.051671419989382,	-0.0511789921672037,	-0.050690546877399,	-0.0502060598966661,	-0.049725507062233,	-0.0492488642727262,	-0.0487761074890269,	-0.0483072127351116,	-0.0478421560988804,	-0.0473809137329706,	-0.0469234618555569,	-0.0464697767511384,	-0.0460198347713116,	-0.0455736123355301,	-0.0451310859318513,	-0.0446922321176693,	-0.0442570275204354,	-0.0438254488383643,	-0.0433974728411284,	-0.0429730763705389,	-0.0425522363412133,	-0.0421349297412315,	-0.0417211336327776,	-0.0413108251527703,	-0.0409039815134806,	-0.0405005800031361,	-0.0401005979865145,	-0.0397040129055231,	-0.0393108022797678,	-0.0389209437071085,	-0.0385344148642034,	-0.0381511935070411,	-0.0377712574714607,	-0.0373945846736605,	-0.0370211531106944,	-0.0366509408609573,	-0.0362839260846586,	-0.0359200870242847,	-0.0355594020050494,	-0.0352018494353339,	-0.0348474078071151,	-0.0344960556963829,	-0.0341477717635468,	-0.0338025347538311,	-0.0334603234976593,	-0.0331211169110289,	-0.0327848939958736,	-0.0324516338404166,	-0.0321213156195127,	-0.0317939185949801,	-0.0314694221159219,	-0.0311478056190377,	-0.0308290486289242,	-0.0305131307583673,	-0.0302000317086224,	-0.0298897312696865,	-0.0295822093205592,	-0.0292774458294954,	-0.0289754208542467,	-0.0286761145422953,	-0.0283795071310763,	-0.0280855789481927,	-0.0277943104116199,	-0.0275056820299013,	-0.027219674402335,	-0.0269362682191517,	-0.0266554442616829,	-0.0263771834025215,	-0.0261014666056726,	-0.0258282749266964,	-0.0255575895128421,	-0.0252893916031733,	-0.0250236625286858,	-0.0247603837124161,	-0.0244995366695423,	-0.0242411030074768,	-0.0239850644259508,	-0.0237314027170911,	-0.0234800997654889,	-0.0232311375482609,	-0.0229844981351026,	-0.0227401636883342,	-0.0224981164629386,	-0.0222583388065926,	-0.0220208131596901,	-0.0217855220553586,	-0.0215524481194679,	-0.0213215740706325,	-0.0210928827202061,	-0.0208663569722698,	-0.0206419798236135,	-0.0204197343637098,	-0.0201996037746824,	-0.0199815713312666,	-0.0197656204007641,	-0.0195517344429915,	-0.0193398970102215,	-0.0191300917471191,	-0.0189223023906704,	-0.0187165127701062,	-0.0185127068068189,	-0.018310868514274,	-0.0181109819979152,	-0.017913031455064,	-0.0177170011748139,	-0.017522875537918,	-0.0173306390166721,	-0.0171402761747913,	-0.0169517716672824,	-0.0167651102403098,	-0.0165802767310568,	-0.0163972560675811,	-0.0162160332686665,	-0.0160365934436679,	-0.0158589217923523,	-0.0156830036047351,	-0.0155088242609106,	-0.0153363692308791,	-0.0151656240743678,	-0.0149965744406487,	-0.0148292060683507,	-0.0146635047852678,	-0.0144994565081631,	-0.0143370472425678,	-0.0141762630825768,	-0.0140170902106398,	-0.0138595148973478,	-0.0137035235012162,	-0.0135491024684638,	-0.013396238332788,	-0.0132449177151355,	-0.0130951273234703,	-0.012946853952537,	-0.0128000844836212,	-0.012654805884306,	-0.012511005208225,	-0.012368669594812,	-0.0122277862690471,	-0.0120883425412001,	-0.0119503258065695,	-0.0118137235452196,	-0.0116785233217138,	-0.0115447127848448,	-0.0114122796673621,	-0.0112812117856966,	-0.0111514970396817,	-0.0110231234122726,	-0.0108960789692621,	-0.0107703518589937,	-0.0106459303120727,	-0.010522802641074,	-0.0104009572402479,	-0.0102803825852231,	-0.0101610672327076,	-0.0100429998201869,	-0.00992616906562036,	-0.00981056376713482,	-0.00969617280271649,	-0.00958298512990046,	-0.0094709897854582,	-0.00936017588508308,	-0.00925053262307385,	-0.00914204927201632,	-0.00903471518246307,	-0.00892851978261147,	-0.00882345257797978,	-0.00871950315108177,	-0.00861666116109945,	-0.00851491634355438,	-0.00841425850997735,	-0.00831467754757652,	-0.00821616341890417,	-0.00811870616152201,	-0.00802229588766514,	-0.00792692278390465,	-0.00783257711080903,	-0.00773924920260427,	-0.00764692946683289,	-0.0075556083840117,	-0.00746527650728861,	-0.00737592446209831,	-0.00728754294581697,	-0.00720012272741602,	-0.00711365464711493,	-0.00702812961603326,	-0.00694353861584172,	-0.00685987269841258,	-0.0067771229854692,	-0.00669528066823502,	-0.00661433700708173,	-0.00653428333117686,	-0.00645511103813087,	-0.00637681159364352,	-0.00629937653114994,	-0.00622279745146599,	-0.00614706602243343,	-0.00607217397856448,	-0.00599811312068615,	-0.0059248753155842,	-0.00585245249564681,	-0.00578083665850797,	-0.00571001986669071,	-0.00563999424725011,	-0.00557075199141613,	-0.00550228535423636,	-0.00543458665421865,	-0.00536764827297377,	-0.00530146265485789,	-0.00523602230661526,	-0.00517131979702076,	-0.00510734775652264,	-0.0050440988768853,	-0.00498156591083225,	-0.00491974167168918,	-0.00485861903302724,	-0.00479819092830657,	-0.00473845035052005,	-0.00467939035183731,	-0.00462100404324911,	-0.00456328459421194,	-0.00450622523229312,	-0.00444981924281613,	-0.0043940599685065,	-0.00433894080913803,	-0.0042844552211795,	-0.00423059671744189,	-0.0041773588667261,	-0.00412473529347114,	-0.00407271967740303,	-0.00402130575318411,	-0.00397048731006304,	-0.00392025819152546,	-0.00387061229494522,	-0.00382154357123633,	-0.00377304602450558,	-0.00372511371170592,	-0.00367774074229049,	-0.00363092127786747,	-0.00358464953185573,	-0.00353891976914116,	-0.00349372630573393,	-0.00344906350842656,	-0.00340492579445276,	-0.00336130763114725,	-0.0033182035356064,	-0.00327560807434979,	-0.00323351586298271,	-0.00319192156585957,	-0.0031508198957483,	-0.00311020561349568,	-0.00307007352769371,	-0.00303041849434694,	-0.00299123541654085,	-0.00295251924411124,	-0.00291426497331468,	-0.00287646764650004,	-0.00283912235178109,	-0.00280222422271016,	-0.00276576843795293,	-0.00272975022096437,	-0.00269416483966576,	-0.00265900760612288,	-0.00262427387622535,	-0.00258995904936713,	-0.00255605856812824,	-0.0025225679179576,	-0.00248948262685716,	-0.00245679826506717,	-0.00242451044475271,	-0.00239261481969145,	-0.00236110708496267,	-0.00232998297663748,	-0.00229923827147038,	-0.00226886878659204,	-0.0022388703792034,	-0.00220923894627101,	-0.00217997042422371,	-0.00215106078865065,	-0.00212250605400052,	-0.00209430227328225,	-0.00206644553776688,	-0.00203893197669093,	-0.00201175775696098,	-0.0019849190828597,	-0.00195841219575319,	-0.00193223337379969,	-0.00190637893165969,	-0.00188084522020738,	-0.00185562862624352,	-0.00183072557220967,	-0.00180613251590385,	-0.00178184595019755,	-0.00175786240275417,	-0.00173417843574892,	-0.00171079064559001,	-0.00168769566264142,	-0.00166489015094694,	-0.00164237080795571,	-0.00162013436424922,	-0.00159817758326966,	-0.00157649726104977,	-0.00155509022594411,	-0.00153395333836178,	-0.00151308349050053,	-0.00149247760608244,	-0.0014721326400909,	-0.00145204557850914,	-0.00143221343806022,	-0.00141263326594838,	-0.00139330213960197,	-0.00137421716641774,	-0.00135537548350663,	-0.00133677425744105,	-0.00131841068400354,	-0.00130028198793701,	-0.00128238542269629,	-0.0012647182702013,	-0.0012472778405916,	-0.00123006147198237,	-0.00121306653022197,	-0.00119629040865085,	-0.00117973052786196,	-0.00116338433546268,	-0.00114724930583811,	-0.00113132293991592,	-0.00111560276493258,	-0.00110008633420112,	-0.00108477122688028,	-0.0010696550477452,	-0.00105473542695947,	-0.00104001001984872,	-0.00102547650667563,	-0.00101113259241635,	-0.000996976006538455,	-0.000983004502780279,	-0.000969215858931722,	-0.000955607876616496,	-0.000942178381075811,	-0.000928925220953503,	-0.000915846268082596,	-0.000902939417273294,	-0.000890202586102413,	-0.000877633714704227,	-0.000865230765562756,	-0.000852991723305458,	-0.000840914594498355,	-0.000828997407442566,	-0.000817238211972252,	-0.000805635079253977,	-0.000794186101587463,	-0.000782889392207758,	-0.000771743085088795,	-0.000760745334748344,	-0.000749894316054363,	-0.000739188224032726,	-0.000728625273676334,	-0.000718203699755615,	-0.000707921756630383,	-0.000697777718063078,	-0.000687769877033366,	-0.000677896545554102,	-0.00066815605448865,	-0.000658546753369554,	-0.000649067010218555,	-0.000639715211367957,	-0.000630489761283324,	-0.000621389082387513,	-0.000612411614886039,	-0.000603555816593755,	-0.000594820162762861,	-0.00058620314591222,	-0.000577703275657982,	-0.000569319078545516,	-0.000561049097882638,	-0.000552891893574128,	-0.000544846041957542,	-0.000536910135640302,	-0.00052908278333806,	-0.000521362609714342,	-0.00051374825522145,	-0.000506238375942629,	-0.000498831643435486,	-0.000491526744576662,	-0.000484322381407745,	-0.000477217270982421,	-0.00047021014521486,	-0.000463299750729325,	-0.000456484848711005,	-0.000449764214758062,	-0.000443136638734886,	-0.000436600924626555,	-0.00043015589039449,	-0.000423800367833307,	-0.000417533202428848,	-0.000411353253217392,	-0.000405259392646049,	-0.000399250506434306,	-0.000393325493436751,	-0.000387483265506936,	-0.000381722747362403,	-0.000376042876450847,	-0.000370442602817414,	-0.000364920888973137,	-0.000359476709764484,	-0.000354109052244043,	-0.000348816915542302,	-0.000343599310740542,	-0.000338455260744834,	-0.000333383800161121,	-0.000328383975171395,	-0.000323454843410949,	-0.000318595473846709,	-0.000313804946656631,	-0.00030908235311016,	-0.000304426795449748,	-0.000299837386773421,	-0.000295313250918384,	-0.000290853522345669,	-0.000286457346025814,	-0.000282123877325561,	-0.000277852281895576,	-0.00027364173555918,	-0.000269491424202085,	-0.000265400543663128,	-0.000261368299626001,	-0.000257393907511964,	-0.000253476592373538,	-0.000249615588789175,	-0.000245810140758892,	-0.000242059501600866,	-0.00023836293384899,	-0.000234719709151368,	-0.000231129108169759,	-0.000227590420479953,	-0.00022410294447308,	-0.000220665987257841,	-0.000217278864563655,	-0.000213940900644715,	-0.000210651428184956,	-0.000207409788203914,	-0.000204215329963476,	-0.000201067410875527,	-0.00019796539641046,	-0.000194908660006576,	-0.000191896582980339,	-0.000188928554437499,	-0.000186003971185066,	-0.000183122237644137,	-0.00018028276576356,	-0.00017748497493444,	-0.000174728291905474,	-0.000172012150699107,	-0.000169335992528509,	-0.000166699265715369,	-0.000164101425608484,	-0.000161541934503158,	-0.000159020261561392,	-0.000156535882732863,	-0.000154088280676681,	-0.00015167694468393,	-0.000149301370600973,	-0.000146961060753523,	-0.000144655523871477,	-0.000142384275014495,	-0.000140146835498334,	-0.000137942732821918,	-0.000135771500595144,	-0.000133632678467418,	-0.000131525812056919,	-0.000129450452880572,	-0.000127406158284742,	-0.000125392491376632,	-0.000123409020956377,	-0.000121455321449845,	-0.000119530972842112,	-0.000117635560611638,	-0.000115768675665106,	-0.00011392991427295,	-0.000112118878005532,	-0.000110335173670001,	-0.00010857841324779,	-0.000106848213832781,	-0.000105144197570102,	-0.000103465991595571,	-0.000101813227975777,	-0.000100185543648784,	-9.85825803654633E-05,	-9.70039846314486E-05,	-9.54494076496957E-05,	-9.39185052636588E-05,	-9.24109379010662E-05,	-9.09263705182951E-05,	-8.94644725453399E-05,	-8.80249178313673E-05,	-8.66073845908552E-05,	-8.52115553503085E-05,	-0.000083837116895548,	-0.000082483760219567,	-8.11511804709499E-05,	-7.98390769028494E-05,	-7.85471528225155E-05,	-7.72751155413732E-05,	-7.60226763256424E-05,	-7.47895503474962E-05,	-7.35754566367528E-05,	-7.23801180330954E-05,	-7.12032611388158E-05,	-7.00446162720775E-05,	-6.89039174206923E-05,	-6.77809021964075E-05,	-6.66753117896968E-05,	-6.55868909250525E-05,	-6.45153878167729E-05,	-6.34605541252409E-05,	-6.24221449136895E-05,	-0.000061399918605449,	-6.03936369416721E-05,	-5.94030649395317E-05,	-5.84279708508882E-05,	-5.74681261214204E-05,	-5.65233053502167E-05,	-0.000055593286249822,	-5.46778496067357E-05,	-5.37767792423567E-05,	-5.28898619743719E-05,	-5.20168875785819E-05,	-5.11576487511625E-05,	-5.03119410713555E-05,	-4.94795629645859E-05,	-4.86603156660011E-05,	-4.78540031844274E-05,	-4.70604322667416E-05,	-4.62794123626507E-05,	-0.000045510755589879,	-4.47542766997561E-05,	-4.40097930432028E-05,	-4.32771245371118E-05,	-4.25560936311172E-05,	-4.18465252747517E-05,	-4.11482468849851E-05,	-0.000040461088314142,	-3.97848818181943E-05,	-3.91194620254247E-05,	-3.84646659054575E-05,	-3.78203327386536E-05,	-3.71863040858651E-05,	-3.65624237585465E-05,	-3.59485377892187E-05,	-3.53444944022824E-05,	-3.47501439851772E-05,	-3.41653390598824E-05,	-3.35899342547574E-05,	-3.30237862767167E-05,	-3.24667538837367E-05,	-3.19186978576911E-05,	-3.13794809775114E-05,	-3.08489679926686E-05,	-0.000030327025596974,	-2.98135224026939E-05,	-2.93083289149771E-05,	-2.88113175065903E-05,	-0.000028322362392959,	-2.78413396075101E-05,	-2.73681269773141E-05,	-2.69026040990221E-05,	-2.64446523150965E-05,	-0.00002599415469033,	-2.55509959886524E-05,	-2.51150626502197E-05,	-2.46862427687843E-05,	-2.42644260693421E-05,	-2.38495038860547E-05,	-2.34413691404424E-05,	-2.30399163198462E-05,	-2.26450414561557E-05,	-0.000022256642104799,	-2.18746173239937E-05,	-2.14988676542541E-05,	-2.11292950981533E-05,	-2.07658031003373E-05,	-2.04082965277865E-05,	-2.00566816503254E-05,	-1.97108661213738E-05,	-1.93707589589403E-05,	-1.90362705268532E-05,	-1.87073125162265E-05,	-1.83837979271603E-05,	-1.80656410506696E-05,	-1.77527574508425E-05,	-1.74450639472226E-05,	-1.71424785974148E-05,	-1.68449206799111E-05,	-0.000016552310677134,	-1.62645702586952E-05,	-1.59816222648684E-05,	-1.57033906902704E-05,	-1.54298006677528E-05,	-1.51607784524977E-05,	-1.48962514063177E-05,	-1.46361479821573E-05,	-1.43803977087923E-05,	-1.41289311757265E-05,	-1.38816800182829E-05,	-1.36385769028867E-05,	-1.33995555125378E-05,	-1.31645505324727E-05,	-0.00001293349763601,	-1.27063334705806E-05,	-1.24829956439394E-05,	-1.22634227105553E-05,	-1.20475541581793E-05,	-1.18353303945875E-05,	-1.16266927344966E-05,	-1.14215833866515E-05,	-1.12199454410811E-05,	-1.10217228565218E-05,	-1.08268604480058E-05,	-1.06353038746126E-05,	-1.04469996273821E-05,	-1.02618950173869E-05,	-1.00799381639618E-05,	-9.90107798308964E-06,	-9.72526417594029E-06,	-9.55244721756192E-06,	-9.38257834572242E-06,	-9.21560954989912E-06,	-9.05149356041516E-06,	-8.89018383772075E-06,	-8.73163456181744E-06,	-8.57580062182394E-06,	-8.42263760568149E-06,	-8.27210178999738E-06,	-8.1241501300247E-06,	-7.97874024977695E-06,	-7.83583043227553E-06,	-7.69537960992889E-06,	-7.5573473550414E-06,	-7.42169387045058E-06,	-0.000007288379980291,	-7.15736712088343E-06,	-7.02861733174753E-06,	-6.90209324673672E-06,	-6.77775808529372E-06,	-6.65557564382512E-06,	-6.53551028719375E-06,	-6.41752694032719E-06,	-6.3015910799411E-06,	-6.1876687263759E-06,	-6.07572643554542E-06,	-5.96573129099617E-06,	-5.85765089607568E-06,	-5.7514533662088E-06,	-5.64710732128042E-06,	-5.54458187812335E-06,	-5.44384664311009E-06,	-5.34487170484702E-06,	-5.24762762696998E-06,	-5.15208544103972E-06,	-5.05821663953611E-06,	-4.96599316894978E-06,	-4.87538742297001E-06,	-4.78637223576768E-06,	-4.69892087537192E-06,	-4.61300703713949E-06,	-4.52860483731551E-06,	-4.44568880668451E-06,	-4.36423388431057E-06,	-4.28421541136546E-06,	-4.2056091250436E-06,	-4.12839115256281E-06,	-4.0525380052496E-06,	-3.97802657270811E-06,	-3.90483411707145E-06,	-3.83293826733445E-06,	-3.76231701376678E-06,	-3.69294870240537E-06,	-3.62481202962509E-06,	-3.55788603678671E-06,	-3.49215010496114E-06,	-3.42758394972888E-06,	-3.36416761605379E-06,	-3.30188147323014E-06,	-3.24070620990205E-06,	-3.18062282915431E-06,	-3.12161264367362E-06,	-3.06365727097943E-06,	-3.00673862872341E-06,	-2.95083893005653E-06,	-2.89594067906308E-06,	-2.84202666626054E-06,	-2.78907996416461E-06,	-2.73708392291827E-06,	-2.68602216598441E-06,	-2.63587858590075E-06,	-2.58663734009656E-06,	-2.53828284677017E-06,	-2.49079978082657E-06,	-2.44417306987411E-06,	-2.39838789027982E-06,	-2.35342966328219E-06,	-2.30928405116097E-06,	-0.000002265936953463,	-2.22337450328338E-06,	-2.18158306360134E-06,	-2.14054922366982E-06,	-2.10025979545835E-06,	-2.06070181014823E-06,	-2.02186251467945E-06,	-1.98372936834862E-06,	-1.94629003945715E-06,	-1.90953240200916E-06,	-1.8734445324582E-06,	-1.8380147065023E-06,	-1.80323139592668E-06,	-1.7690832654933E-06,	-1.73555916987683E-06,	-1.70264815064621E-06,	-1.67033943329132E-06,	-1.63862242429403E-06,	-1.60748670824309E-06,	-1.57692204499217E-06,	-1.54691836686058E-06,	-1.51746577587593E-06,	-1.48855454105829E-06,	-1.46017509574513E-06,	-1.4323180349566E-06,	-1.40497411280051E-06,	-1.37813423991647E-06,	-1.35178948095867E-06,	-1.32593105211675E-06,	-1.30055031867415E-06,	-1.27563879260354E-06,	-1.25118813019873E-06,	-1.22719012974259E-06,	-1.20363672921033E-06,	-1.18052000400796E-06,	-1.15783216474499E-06,	-1.13556555504135E-06,	-1.11371264936763E-06,	-1.09226605091848E-06,	-1.07121848951855E-06,	-1.05056281956045E-06,	-1.0302920179745E-06,	-1.01039918222952E-06,	-9.90877528364475E-07,	-9.717203890504E-07,	-9.5292121168216E-07,	-9.34473556499699E-07,	-9.16371094738292E-07,	-8.98607606807406E-07,	-8.81176980497744E-07,	-8.64073209216093E-07,	-8.47290390247536E-07,	-8.30822723044669E-07,	-8.14664507543407E-07,	-7.98810142505009E-07,	-7.83254123883929E-07,	-7.67991043221132E-07,	-7.53015586062488E-07,	-7.38322530401897E-07,	-7.23906745148766E-07,	-7.097631886195E-07,	-6.95886907052638E-07,	-6.82273033147304E-07,	-6.68916784624623E-07,	-6.55813462811765E-07,	-6.42958451248286E-07,	-6.30347214314436E-07,	-6.17975295881119E-07,	-6.0583831798117E-07,	-5.93931979501652E-07,	-5.82252054896842E-07,	-5.70794392921625E-07,	-5.59554915384968E-07,	-5.48529615923195E-07,	-5.37714558792756E-07,	-5.27105877682208E-07,	-5.16699774543125E-07,	-5.06492518439635E-07,	-4.96480444416334E-07,	-4.86659952384282E-07,	-4.77027506024817E-07,	-4.6757963171092E-07,	-4.58312917445865E-07,	-4.49224011818906E-07,	-4.4030962297772E-07,	-4.31566517617384E-07,	-4.22991519985616E-07,	-4.14581510904036E-07,	-4.0633342680522E-07,	-3.98244258785287E-07,	-3.90311051671801E-07,	-3.82530903106749E-07,	-3.7490096264436E-07,	-3.67418430863555E-07,	-3.60080558494789E-07,	-3.52884645561072E-07,	-3.45828040532963E-07,	-3.38908139497302E-07,	-3.32122385339488E-07,	-3.25468266939093E-07,	-3.18943318378598E-07,	-3.12545118165058E-07,	-3.06271288464499E-07,	-3.00119494348842E-07,	-2.94087443055173E-07,	-2.88172883257156E-07,	-2.82373604348416E-07,	-2.76687435737689E-07,	-2.71112246155579E-07,	-2.65645942972722E-07,	-2.60286471529191E-07,	-2.55031814474966E-07,	-2.49879991121299E-07,	-2.44829056802797E-07,	-2.39877102250063E-07,	-2.35022252972728E-07,	-2.30262668652719E-07,	-2.25596542547581E-07,	-2.21022100903729E-07,	-2.16537602379444E-07,	-2.12141337477478E-07,	-2.07831627987112E-07,	-2.03606826435518E-07,	-1.99465315548276E-07,	-1.95405507718913E-07,	-1.91425844487303E-07,	-1.87524796026812E-07,	-1.83700860640023E-07,	-1.7995256426293E-07,	-1.76278459977451E-07,	-1.72677127532132E-07,	-1.69147172870919E-07,	-1.65687227669859E-07,	-1.62295948881608E-07,	-1.58972018287627E-07,	-1.55714142057939E-07,	-1.52521050318319E-07,	-1.49391496724817E-07,	-1.46324258045477E-07,	-1.43318133749149E-07,	-1.40371945601276E-07,	-1.37484537266544E-07,	-1.34654773918286E-07,	-1.31881541854531E-07,	-1.29163748120596E-07,	-1.26500320138099E-07,	-1.2389020534031E-07,	-1.21332370813732E-07,	-1.18825802945796E-07,	-1.16369507078595E-07,	-1.13962507168539E-07,	-1.1160384545185E-07,	-1.09292582115789E-07,	-1.07027794975538E-07,	-1.04808579156621E-07,	-1.02634046782802E-07,	-1.00503326669354E-07,	-9.8415564021605E-08,	-9.6369920138697E-08,	-9.43655721224524E-08,	-9.24017125912739E-08,	-9.04775493989929E-08,	-8.85923053585864E-08,	-8.67452179706806E-08,	-8.49355391567648E-08,	-8.3162534997036E-08,	-8.14254854727997E-08,	-7.97236842133505E-08,	-7.80564382472591E-08,	-7.64230677579931E-08,	-7.48229058437995E-08,	-7.32552982817784E-08,	-7.17196032960795E-08,	-7.02151913301506E-08,	-6.87414448229738E-08,	-6.72977579892208E-08,	-6.58835366032632E-08,	-6.44981977869735E-08,	-6.3141169801253E-08,	-6.18118918412255E-08,	-6.05098138350342E-08,	-5.92343962461831E-08,	-5.79851098793617E-08,	-5.67614356896968E-08,	-5.55628645953717E-08,	-5.43888972935581E-08,	-5.3239044079604E-08,	-5.21128246694233E-08,	-5.10097680250322E-08,	-4.99294121831816E-08,	-4.88713040870302E-08,	-4.78349994208097E-08,	-4.68200624474296E-08,	-4.58260658489734E-08,	-4.48525905700355E-08,	-4.38992256638523E-08,	-4.2965568141179E-08,	-4.20512228218656E-08,	-4.11558021890869E-08,	-4.02789262461805E-08,	-3.94202223760485E-08,	-3.85793252030794E-08,	-3.77558764575475E-08,	-3.69495248424464E-08,	-3.61599259027163E-08,	-3.53867418968227E-08,	-3.46296416706478E-08,	-3.3888300533653E-08,	-3.31624001372757E-08,	-3.24516283555199E-08,	-3.17556791677043E-08,	-3.10742525433299E-08,	-3.04070543290316E-08,	-2.97537961375761E-08,	-2.91141952388726E-08,	-2.84879744529606E-08,	-2.78748620449403E-08,	-2.72745916218124E-08,	-2.66869020311951E-08,	-2.61115372618829E-08,	-2.55482463462188E-08,	-2.49967832642461E-08,	-2.44569068496085E-08,	-2.39283806971702E-08,	-2.3410973072324E-08,	-2.29044568219587E-08,	-2.24086092870569E-08,	-2.19232122168946E-08,	-2.14480516848141E-08,	-2.09829180055437E-08,	-2.05276056540362E-08,	-2.00819131857992E-08,	-1.96456431586925E-08,	-1.92186020561653E-08,	-1.88006002119088E-08,	-1.83914517358981E-08,	-1.79909744418007E-08,	-1.75989897757262E-08,	-1.72153227462931E-08,	-1.6839801855991E-08,	-1.64722590338141E-08,	-1.61125295691432E-08,	-1.57604520468555E-08,	-1.54158682836389E-08,	-1.50786232654902E-08,	-1.4748565086376E-08,	-1.44255448880357E-08,	-1.41094168009059E-08,	-1.38000378861465E-08,	-1.34972680787491E-08,	-1.32009701317076E-08,	-1.29110095612328E-08,	-1.2627254592992E-08,	-1.23495761093547E-08,	-1.20778475976278E-08,	-1.18119450992611E-08,	-1.1551747160006E-08,	-1.12971347810108E-08,	-1.10479913708358E-08,	-1.08042026983705E-08,	-1.0565656846638E-08,	-1.03322441674697E-08,	-1.0103857237035E-08,	-9.88039081221069E-09,	-9.66174178777415E-09,	-9.44780915440622E-09,	-9.23849395748867E-09,	-9.03369925668189E-09,	-8.83333008626876E-09,	-8.63729341625064E-09,	-8.44549811418193E-09,	-8.25785490772965E-09,	-8.07427634794488E-09,	-7.89467677323317E-09,	-7.71897227401094E-09,	-7.5470806580356E-09,	-7.37892141639689E-09,	-7.21441569015746E-09,	-7.05348623763066E-09,	-6.89605740228405E-09,	-6.74205508125697E-09,	-6.59140669448106E-09,	-6.44404115439249E-09,	-6.29988883622516E-09,	-6.15888154887402E-09,	-6.02095250631824E-09,	-5.88603629959362E-09,	-5.75406886930434E-09,	-5.62498747866398E-09,	-5.49873068705599E-09,	-5.37523832410411E-09,	-5.25445146424322E-09,	-5.13631240178134E-09,	-5.0207646264437E-09,	-4.90775279938992E-09,	-4.79722272969554E-09,	-4.68912135128913E-09,	-4.58339670033671E-09,	-4.479997893065E-09,	-4.37887510401532E-09,	-4.27997954472017E-09,	-4.18326344279455E-09,	-4.0886800214343E-09,	-3.99618347931378E-09,	-3.90572897087551E-09,	-3.81727258700433E-09,	-3.7307713360789E-09,	-3.64618312539352E-09,	-3.56346674294321E-09,	-3.48258183956524E-09,	-3.40348891143051E-09,	-3.32614928287804E-09,	-3.25052508958616E-09,	-3.17657926207412E-09,	-3.10427550952777E-09,	-3.03357830394321E-09,	-2.96445286458249E-09,	-2.8968651427353E-09,	-2.83078180678101E-09,	-2.76617022754521E-09,	-2.70299846394527E-09,	-2.64123524891946E-09,	-2.58084997563406E-09,	-2.52181268396339E-09,	-2.46409404723741E-09,	-2.40766535925185E-09,	-2.3524985215359E-09,	-2.29856603087244E-09,	-2.24584096706609E-09,	-2.19429698095427E-09,	-2.14390828265668E-09,	-2.09464963005855E-09,	-2.04649631752334E-09,	-1.99942416483021E-09,	-1.95340950633229E-09,	-1.90842918033125E-09,	-1.86446051866406E-09,	-1.82148133649798E-09,	-1.77946992232962E-09,	-1.73840502818421E-09,	-1.69826586001121E-09,	-1.65903206827247E-09,	-1.6206837387192E-09,	-1.58320138335412E-09,	-1.54656593157519E-09,	-1.51075872149746E-09,	-1.47576149144945E-09,	-1.44155637164082E-09,	-1.40812587599798E-09,	-1.37545289416427E-09,	-1.34352068366162E-09,	-1.31231286221052E-09,	-1.2818134002052E-09,	-1.25200661334099E-09,	-1.2228771553909E-09,	-1.19441001112858E-09,	-1.16659048939459E-09,	-1.13940421630345E-09,	-1.11283712858854E-09,	-1.08687546708213E-09,	-1.06150577032807E-09,	-1.03671486832429E-09,	-1.01248987639284E-09,	-9.8881818917463E-10,	-9.65687474746813E-10,	-9.43085668860012E-10,	-9.21000969293311E-10,	-8.99421830324543E-10,	-8.78336957313657E-10,	-8.57735301396913E-10,	-8.37606054289719E-10,	-8.17938643195963E-10,	-7.98722725821741E-10,	-7.79948185491387E-10,	-7.61605126363818E-10,	-7.43683868747164E-10,	-7.2617494450976E-10,	-7.09069092585576E-10,	-6.923572545722E-10,	-6.76030570419555E-10,	-6.6008037420751E-10,	-6.4449819001065E-10,	-6.29275727848441E-10,	-6.14404879719091E-10,	-5.99877715715452E-10,	-5.85686480221292E-10,	-5.71823588186363E-10,	-5.58281621478664E-10,	-5.45053325312364E-10,	-5.32131604749858E-10,	-5.19509521276489E-10,	-5.07180289446447E-10,	-4.95137273598441E-10,	-4.83373984639731E-10,	-4.71884076897138E-10,	-4.60661345033694E-10,	-4.49699721029602E-10,	-4.38993271226214E-10,	-4.28536193431743E-10,	-4.18322814087483E-10,	-4.08347585493291E-10,	-3.98605083091149E-10,	-3.89090002805617E-10,	-3.79797158440038E-10,	-3.70721479127349E-10,	-3.61858006834405E-10,	-3.53201893918718E-10,	-3.44748400736551E-10,	-3.3649289330132E-10,	-3.28430840991282E-10,	-3.20557814305509E-10,	-3.12869482667147E-10,	-3.05361612273025E-10,	-2.98030063988631E-10,	-2.90870791287565E-10,	-2.83879838234528E-10,	-2.77053337510978E-10,	-2.70387508482572E-10,	-2.63878655307532E-10,	-2.57523165085105E-10,	-2.51317506043292E-10,	-2.4525822576503E-10,	-2.39341949452055E-10,	-2.33565378225654E-10,	-2.27925287463555E-10,	-2.22418525172208E-10,	-2.17042010393732E-10,	-2.11792731646802E-10,	-2.06667745400796E-10,	-2.01664174582486E-10,	-1.96779207114625E-10,	-1.92010094485765E-10,	-1.87354150350642E-10,	-1.82808749160521E-10,	-1.78371324822865E-10,	-1.74039369389711E-10,	-1.69810431774182E-10,	-1.65682116494528E-10,	-1.61652082445132E-10,	-1.57718041693919E-10,	-1.53877758305621E-10,	-1.50129047190358E-10,	-1.46469772976997E-10,	-1.42897848910792E-10,	-1.3941123577478E-10,	-1.36007940834443E-10,	-1.32686016805153E-10,	-1.29443560841914E-10,	-1.26278713550949E-10,	-1.23189658022652E-10,	-1.20174618885483E-10,	-1.17231861380342E-10,	-1.14359690455009E-10,	-1.1155644987822E-10,	-1.08820521372971E-10,	-1.06150323768639E-10,	-1.03544312171523E-10,	-1.01000977153435E-10,	-9.85188439579279E-11,	-9.60964717238184E-11,	-9.37324527256234E-11,	-9.14254116305546E-11,	-8.91740047717246E-11,	-8.69769194372177E-11,	-8.48328731746923E-11,	-8.27406131111822E-11,	-8.06989152877781E-11,	-7.87065840088707E-11,	-7.67624512056465E-11,	-7.48653758135349E-11,	-7.30142431633071E-11,	-7.12079643855396E-11,	-6.94454758281532E-11,	-6.77257384867532E-11,	-6.6047737447495E-11,	-6.44104813422078E-11,	-6.28130018155155E-11,	-6.12543530036984E-11,	-5.97336110250429E-11,	-5.8249873481435E-11,	-5.68022589709562E-11,	-5.53899066112451E-11,	-5.4011975573395E-11,	-5.26676446261604E-11,	-5.13561116902514E-11,	-5.00765934024989E-11,	-4.88283246896781E-11,	-4.76105583517829E-11,	-4.6422564654547E-11,	-4.52636309310129E-11,	-4.41330611919529E-11,	-4.30301757449521E-11,	-4.1954310821965E-11,	-4.09048182151633E-11,	-3.98810649208957E-11,	-3.88824327915835E-11,	-3.79083181953803E-11,	-3.69581316834275E-11,	-3.60312976645417E-11,	-3.51272540871708E-11,	-3.42454521284623E-11,	-3.33853558902895E-11,	-3.25464421020832E-11,	-3.17281998303208E-11,	-3.09301301945293E-11,	-3.01517460896578E-11,	-2.93925719146833E-11,	-2.86521433073116E-11,	-2.79300068846407E-11,	-2.72257199896577E-11,	-2.65388504434395E-11,	-2.5868976302934E-11,	-2.52156856241994E-11,	-2.45785762309816E-11,	-2.39572554885121E-11,	-2.33513400824139E-11,	-2.27604558026002E-11,	-2.21842373320596E-11,	-2.16223280404169E-11,	-2.10743797821676E-11,	-2.05400526994808E-11,	-2.00190150294712E-11,	-1.95109429158409E-11,	-1.90155202247953E-11,	-1.85324383651384E-11,	-1.80613961124548E-11,	-1.76020994372887E-11,	-1.71542613372313E-11,	-1.67176016728293E-11,	-1.62918470072318E-11,	-1.58767304494902E-11,	-1.54719915014323E-11,	-1.50773759080295E-11,	-1.46926355111814E-11,	-1.43175281068398E-11,	-1.39518173053998E-11,	-1.35952723952831E-11,	-1.32476682096451E-11,	-1.29087849961338E-11,	-1.25784082896338E-11,	-1.22563287879281E-11,	-1.19423422302132E-11,	-1.16362492784034E-11,	-1.13378554011616E-11,	-1.10469707605964E-11,	-1.07634101015652E-11,	-1.04869926435251E-11,	-1.02175419748744E-11,	-9.95488594972942E-12,	-9.69885658708066E-12,	-9.44928997227657E-12,	-9.20602616078122E-12,	-8.9689090841555E-12,	-8.73778645821146E-12,	-8.51250969329111E-12,	-8.29293380662158E-12,	-8.0789173367E-12,	-7.87032225966221E-12,	-7.66701390759055E-12,	-7.46886088871692E-12,	-7.27573500947827E-12,	-7.08751119838261E-12,	-6.90406743164454E-12,	-6.72528466055034E-12,	-6.55104674051323E-12,	-6.38124036178059E-12,	-6.2157549817557E-12,	-6.0544827588972E-12,	-5.89731848816055E-12,	-5.74415953794632E-12,	-5.59490578852113E-12,	-5.44945957187764E-12,	-5.30772561300087E-12,	-5.16961097250873E-12,	-5.03502499063544E-12,	-4.90387923252724E-12,	-4.77608743482032E-12,	-4.65156545347174E-12,	-4.53023121281472E-12,	-4.41200465581016E-12,	-4.2968076954671E-12,	-4.18456416740535E-12,	-4.07519978353401E-12,	-3.96864208682037E-12,	-3.86482040712415E-12,	-3.76366581807254E-12,	-3.66511109495219E-12,	-3.56909067359477E-12,	-3.47554061023316E-12,	-3.38439854230598E-12,	-3.2956036501886E-12,	-3.20909661982927E-12,	-3.12481960626944E-12,	-3.04271619802796E-12,	-2.96273138232905E-12,	-2.88481151115473E-12,	-2.80890426810248E-12,	-2.73495863602955E-12,	-2.66292486546576E-12,	-2.59275444377685E-12,	-2.52440006506111E-12,	-2.45781560076225E-12,	-2.39295607098174E-12,	-2.32977761647464E-12,	-2.26823747131274E-12,	-2.20829393619973E-12,	-2.14990635242303E-12,	-2.09303507642763E-12,	-2.03764145499733E-12,	-1.98368780102921E-12,	-1.93113736988756E-12,	-1.87995433632374E-12,	-1.83010377194865E-12,	-1.781551623245E-12,	-1.73426469010675E-12,	-1.68821060489323E-12,	-1.64335781198617E-12,	-1.59967554783753E-12,	-1.5571338214969E-12,	-1.51570339560704E-12,	-1.47535576785664E-12,	-1.43606315287956E-12,	-1.39779846459005E-12,	-1.36053529894367E-12,	-1.32424791711393E-12,	-1.28891122907486E-12,	-1.25450077757987E-12,	-1.22099272252771E-12,	-1.18836382570619E-12,	-1.15659143590497E-12,	-1.12565347438849E-12,	-1.09552842072073E-12,	-1.06619529893326E-12,	-1.03763366402866E-12,	-1.00982358881127E-12,	-9.82745651037487E-13,	-9.56380920878127E-13,	-9.30710948685349E-13,	-9.05717753056993E-13,	-8.81383809191222E-13,	-8.57692037524586E-13,	-8.34625792646766E-13,	-8.12168852485411E-13,	-7.90305407754642E-13,	-7.69020051660949E-13,	-7.48297769860327E-13,	-7.28123930660676E-13,	-7.0848427546361E-13,	-6.89364909439944E-13,	-6.70752292433295E-13,	-6.52633230086335E-13,	-6.34994865184359E-13,	-6.17824669210985E-13,	-6.011104341109E-13,	-5.84840264254675E-13,	-5.69002568600826E-13,	-5.53586053050374E-13,	-5.38579712989286E-13,	-5.23972826014283E-13,	-5.09754944837617E-13,	-4.95915890366491E-13,	-4.82445744952942E-13,	-4.69334845810071E-13,	-4.56573778590613E-13,	-4.44153371123932E-13,	-4.32064687307629E-13,	-4.20299021150014E-13,	-4.0884789095982E-13,	-3.97703033679574E-13,	-3.86856399359184E-13,	-3.76300145766321E-13,	-3.66026633130295E-13,	-3.56028419016204E-13,	-3.46298253326172E-13,	-3.3682907342462E-13,	-3.27613999384542E-13,	-3.18646329351864E-13,	-3.09919535024999E-13,	-3.01427257246812E-13,	-2.93163301706259E-13,	-2.85121634747024E-13,	-2.77296379280555E-13,	-2.6968181080095E-13,	-2.62272353499223E-13,	-2.55062576474502E-13,	-2.48047190039817E-13,	-2.41221042120156E-13,	-2.34579114740535E-13,	-2.28116520601882E-13,	-2.21828499742588E-13,	-2.15710416283631E-13,	-2.09757755255213E-13,	-2.03966119502931E-13,	-1.98331226671517E-13,	-1.92848906264249E-13,	-1.87515096776178E-13,	-1.82325842899353E-13,	-1.77277292798276E-13,	-1.72365695453862E-13,	-1.67587398074214E-13,	-1.62938843570576E-13,	-1.58416568096846E-13,	-1.54017198651092E-13,	-1.4973745073754E-13,	-1.45574126087543E-13,	-1.41524110438061E-13,	-1.37584371366261E-13,	-1.33751956178809E-13,	-1.30023989854546E-13,	-1.26397673039184E-13,	-1.22870280090767E-13,	-1.19439157174621E-13,	-1.1610172040657E-13,	-1.12855454043228E-13,	-1.09697908718185E-13,	-1.06626699722961E-13,	-1.03639505331601E-13,	-1.00734065167849E-13,	-9.7908178613809E-14,	-9.51597032590948E-14,	-9.24865533894327E-14,	-8.98866985137515E-14,	-8.73581619287912E-14,	-8.48990193202983E-14,	-8.25073973998929E-14,	-8.01814725767172E-14,	-7.79194696629963E-14,	-7.57196606126628E-14,	-7.35803632922193E-14,	-7.14999402830294E-14,	-6.94767977142538E-14,	-6.75093841256591E-14,	-6.55961893595538E-14,	-6.37357434811192E-14,	-6.19266157264225E-14,	-6.01674134774168E-14,	-5.84567812632502E-14,	-5.67933997872224E-14,	-5.51759849787427E-14,	-5.36032870696611E-14,	-5.20740896943586E-14,	-5.05872090129965E-14,	-4.91414928573423E-14,	-4.77358198986027E-14,	-4.6369098836706E-14,	-4.50402676104954E-14,	-4.37482926283015E-14,	-4.24921680183808E-14,	-4.12709148987175E-14,	-4.00835806656969E-14,	-3.89292383011739E-14,	-3.78069856974702E-14,	-3.6715944999845E-14,	-3.56552619659962E-14,	-3.46241053421601E-14,	-3.36216662553875E-14,	-3.26471576215851E-14,	-3.16998135689227E-14,	-3.07788888762125E-14,	-2.98836584258824E-14,	-2.901341667117E-14,	-2.81674771171743E-14,	-2.73451718154141E-14,	-2.65458508715457E-14,	-2.57688819659064E-14,	-2.50136498865544E-14,	-2.42795560744877E-14,	-2.35660181807277E-14,	-2.28724696349681E-14,	-2.21983592254888E-14,	-2.15431506900492E-14,	-2.09063223174788E-14,	-2.02873665596902E-14,	-1.9685789653848E-14,	-1.91011112544322E-14,	-1.8532864074943E-14,	-1.79805935389984E-14,	-1.7443857440583E-14,	-1.69222256132137E-14,	-1.64152796077919E-14,	-1.59226123789186E-14,	-1.54438279794546E-14,	-1.49785412631137E-14,	-1.45263775948802E-14,	-1.40869725690508E-14,	-1.36599717347026E-14,	-1.32450303283955E-14,	-1.28418130139231E-14,	-1.24499936289282E-14,	-1.20692549382066E-14,	-1.1699288393525E-14,	-1.13397938997859E-14,	-1.09904795873724E-14,	-1.06510615905152E-14,	-1.0321263831524E-14,	-1.00008178107321E-14,	-9.68946240200426E-15,	-9.38694365366523E-15,	-9.09301459470653E-15,	-8.80743504613461E-15,	-8.52997143732681E-15,	-8.26039662726431E-15,	-7.99848973051524E-15,	-7.74403594784397E-15,	-7.49682640132605E-15,	-7.256657973851E-15,	-7.02333315289858E-15,	-6.79665987847671E-15,	-6.57645139511238E-15,	-6.36252610778944E-15,	-6.15470744173011E-15,	-5.95282370591939E-15,	-5.75670796027448E-15,	-5.56619788636361E-15,	-5.38113566158116E-15,	-5.20136783668843E-15,	-5.02674521663169E-15,	-4.85712274455143E-15,	-4.69235938889889E-15,	-4.53231803357829E-15,	-4.37686537103498E-15,	-4.22587179821214E-15,	-4.07921131530037E-15,	-3.93676142720662E-15,	-3.79840304767075E-15,	-3.66402040595989E-15,	-3.53350095607257E-15,	-3.40673528838633E-15,	-3.28361704368427E-15,	-3.16404282949766E-15,	-3.04791213870329E-15,	-2.93512727031604E-15,	-2.82559325241834E-15,	-2.71921776717019E-15,	-2.61591107784424E-15,	-2.51558595783262E-15,	-2.41815762157292E-15,	-2.32354365734249E-15,	-2.23166396187141E-15,	-2.14244067672586E-15,	-2.05579812641472E-15,	-1.97166275817367E-15,	-1.88996308338204E-15,	-1.81062962056904E-15,	-1.73359483996688E-15,	-1.65879310956976E-15,	-1.58616064265831E-15,	-1.51563544675058E-15,	-1.44715727394143E-15,	-1.38066757259319E-15,	-1.31610944034149E-15,	-1.25342757838119E-15,	-1.19256824699804E-15,	-1.13347922231274E-15,	-1.07610975420501E-15,	-1.020410525386E-15,	-9.66333611588232E-16,	-9.13832442843147E-16,	-8.62861765817032E-16,	-8.13377607176934E-16,	-7.65337237958871E-16,	-7.18699138911421E-16,	-6.73422966788439E-16,	-6.29469521565371E-16,	-5.8680071455432E-16,	-5.45379537393625E-16,	-5.05170031888413E-16,	-4.66137260679167E-16,	-4.28247278715987E-16,	-3.91467105516783E-16,	-3.55764698188241E-16,	-3.21108925188957E-16,	-2.87469540814671E-16,	-2.54817160386079E-16,	-2.23123236120202E-16,	-1.92360033666821E-16,	-1.62500609291954E-16,	-1.33518787690841E-16,	-1.0538914041338E-16,	-7.80869648853761E-17,	-5.1588264009464E-17,	-2.5869726329937E-17,	-9.08706746173374E-19,	2.33167922402527E-17,	4.68281387592897E-17,	6.9646089301088E-17,	9.17908059533839E-17,	1.1328187316097E-16,	1.34138314020284E-16,	1.54378606121775E-16,	1.74020696952405E-16,	1.93082018870268E-16,	2.11579503663029E-16,	2.29529596701546E-16,	2.46948270699735E-16,	2.63851039091455E-16,	2.80252969034894E-16,	2.96168694054637E-16,	3.11612426331361E-16,	3.26597968648797E-16,	3.41138726007369E-16,	3.55247716913633E-16,	3.68937584354444E-16,	3.82220606464475E-16,	3.95108706895534E-16,	4.07613464895875E-16,	4.19746125107457E-16,	4.31517607088942E-16,	4.42938514571958E-16,	4.5401914445798E-16,	4.6476949556298E-16,	4.75199277116798E-16,	4.85317917023994E-16,	4.95134569892769E-16,	5.04658124838361E-16,	5.13897213067133E-16,	5.22860215247433E-16,	5.31555268673098E-16,	5.39990274225371E-16,	5.4817290313877E-16,	5.56110603576379E-16,	5.63810607019807E-16,	5.71279934478977E-16,	5.78525402526728E-16,	5.85553629163091E-16,	5.92371039513983E-16,	5.98983871368889E-16,	6.0539818056204E-16,	6.11619846201409E-16,	6.1765457574978E-16,	6.23507909961991E-16,	6.29185227682365E-16,	6.34691750506218E-16,	6.4003254730924E-16,	6.45212538648421E-16,	6.5023650103812E-16,	6.55109071104748E-16,	6.59834749623471E-16,	6.64417905440209E-16,	6.68862779282155E-16,	6.73173487459923E-16,	6.7735402546436E-16,	6.81408271460972E-16,	6.85339989684819E-16,	6.89152833738699E-16,	6.92850349797289E-16,	6.96435979719918E-16,	6.99913064074519E-16,	7.03284845075259E-16,	7.06554469436255E-16,	7.09724991143773E-16,	7.12799374149155E-16,	7.15780494984739E-16,	7.18671145304923E-16,	7.21474034354489E-16,	7.24191791366222E-16,	7.26826967889841E-16,	7.29382040054152E-16,	7.3185941076433E-16,	7.34261411836141E-16,	7.36590306068908E-16,	7.38848289258922E-16,	7.41037492155015E-16,	7.43159982357902E-16,	7.452177661649E-16,	7.47212790361563E-16,	7.49146943961742E-16,	7.51022059897524E-16,	7.52839916660478E-16,	7.54602239895586E-16,	7.563107039492E-16,	7.57966933372333E-16,	7.59572504380555E-16,	7.61128946271716E-16,	7.62637742802714E-16,	7.64100333526451E-16,	7.6551811509013E-16,	7.6689244249598E-16,	7.68224630325484E-16,	7.69515953928154E-16,	7.70767650575855E-16,	7.7198092058367E-16,	7.73156928398255E-16,	7.74296803654616E-16,	7.75401642202206E-16,	7.76472507101221E-16,	7.77510429589953E-16,	7.78516410024013E-16,	7.79491418788246E-16,	7.80436397182114E-16,	7.81352258279295E-16,	7.82239887762261E-16,	7.83100144732532E-16,	7.8393386249731E-16,	7.84741849333177E-16,	7.85524889227512E-16,	7.86283742598254E-16,	7.87019146992666E-16,	7.87731817765661E-16,	7.88422448738319E-16,	7.89091712837136E-16,	7.89740262714574E-16,	7.90368731351448E-16,	7.90977732641671E-16,	7.91567861959872E-16,	7.9213969671237E-16,	7.92693796872E-16,	7.9323070549724E-16,	7.93750949236111E-16,	7.94255038815272E-16,	7.94743469514752E-16,	7.95216721628727E-16,	7.95675260912753E-16,	7.96119539017836E-16,	7.96549993911732E-16,	7.96967050287838E-16,	7.97371119962048E-16,	7.97762602257901E-16,	7.9814188438038E-16,	7.98509341778686E-16,	7.98865338498305E-16,	7.99210227522679E-16,	7.99544351104784E-16,	7.9986804108892E-16,	8.0018161922297E-16,	8.00485397461441E-16,	8.00779678259526E-16,	8.01064754858467E-16,	8.01340911562462E-16,	8.01608424007371E-16,	8.01867559421456E-16,	8.02118576878384E-16,	8.02361727542728E-16,	8.02597254908181E-16,	8.02825395028692E-16,	8.03046376742736E-16,	8.03260421890921E-16,	8.0346774552712E-16,	8.03668556123324E-16,	8.038630557684E-16,	8.04051440360922E-16,	8.04233899796269E-16,	8.0441061814814E-16,	8.04581773844658E-16,	8.04747539839218E-16,	8.04908083776246E-16,	8.05063568151996E-16,	8.05214150470552E-16,	8.0535998339517E-16,	8.05501214895086E-16,	8.05637988387944E-16,	8.05770442877952E-16,	8.05898713089914E-16,	8.06022929599232E-16,	8.06143218958029E-16,	8.0625970381748E-16,	8.06372503046479E-16,	8.06481731846745E-16,	8.06587501864471E-16,	8.06689921298626E-16,	8.06789095005996E-16,	8.0688512460307E-16,	8.06978108564868E-16,	8.07068142320787E-16,	8.07155318347569E-16,	8.07239726259469E-16,	8.07321452895701E-16,	8.07400582405259E-16,	8.07477196329172E-16,	8.0755137368028E-16,	8.07623191020606E-16,	8.0769272253639E-16,	8.07760040110854E-16,	8.07825213394775E-16,	8.07888309874917E-16,	8.07949394940402E-16,	8.08008531947063E-16,	8.08065782279857E-16,	8.08121205413386E-16,	8.08174858970581E-16,	8.08226798779615E-16,	8.0827707892908E-16,	8.08325751821501E-16,	8.08372868225216E-16,	8.08418477324688E-16,	8.08462626769288E-16,	8.08505362720588E-16,	8.08546729898232E-16,	8.08586771624389E-16,	8.0862552986688E-16,	8.08663045280968E-16,	8.08699357249891E-16,	8.08734503924155E-16,	8.08768522259625E-16,	8.08801448054467E-16,	8.08833315984944E-16,	8.0886415964014E-16,	8.08894011555607E-16,	8.08922903245997E-16,	8.08950865236683E-16,	8.08977927094431E-16,	8.09004117457115E-16,	8.09029464062537E-16,	8.09053993776358E-16,	8.09077732619176E-16,	8.0910070579277E-16,	8.09122937705546E-16,	8.09144451997199E-16,	8.09165271562619E-16,	8.09185418575061E-16,	8.09204914508608E-16,	8.09223780159944E-16,	8.09242035669448E-16,	8.09259700541655E-16,	8.09276793665075E-16,	8.09293333331408E-16,	8.09309337254167E-16,	8.09324822586725E-16,	8.0933980593981E-16,	8.09354303398454E-16,	8.09368330538431E-16,	8.09381902442171E-16,	8.09395033714201E-16,	8.09407738496095E-16,	8.0942003048097E-16,	8.09431922927533E-16,	8.09443428673688E-16,	8.09454560149732E-16,	8.09465329391131E-16,	8.09475748050915E-16,	8.09485827411672E-16,	8.09495578397185E-16,	8.09505011583701E-16,	8.09514137210853E-16,	8.09522965192237E-16,	8.09531505125676E-16,	8.09539766303143E-16,	8.095477577204E-16,	8.09555488086317E-16,	8.09562965831911E-16,	8.09570199119108E-16,	8.09577195849221E-16,	8.09583963671171E-16,	8.09590509989456E-16,	8.0959684197186E-16,	8.09602966556932E-16,	8.0960889046123E-16,	8.09614620186332E-16,	8.09620162025638E-16,	8.09625522070957E-16,	8.09630706218878E-16,	8.09635720176959E-16,	8.0964056946971E-16,	8.09645259444393E-16,	8.0964979527664E-16,	8.09654181975897E-16,	8.09658424390695E-16,	8.09662527213755E-16,	8.09666494986937E-16,	8.09670332106035E-16,	8.09674042825409E-16,	8.09677631262491E-16,	8.09681101402132E-16,	8.09684457100826E-16,	8.09687702090788E-16,	8.0969083998392E-16,	8.09693874275636E-16,	8.09696808348576E-16,	8.09699645476203E-16,	8.09702388826285E-16,	8.09705041464263E-16,	8.09707606356522E-16,	8.09710086373551E-16,	8.09712484293008E-16,	8.09714802802686E-16,	8.09717044503383E-16,	8.09719211911692E-16,	8.09721307462686E-16,	8.09723333512535E-16,	8.09725292341027E-16,	8.09727186154019E-16,	8.09729017085803E-16,	8.09730787201404E-16,	8.09732498498798E-16,	8.09734152911065E-16,	8.09735752308474E-16,	8.09737298500497E-16,	8.09738793237764E-16,	8.09740238213952E-16,	8.09741635067619E-16,	8.09742985383975E-16,	8.097442906966E-16,	8.09745552489103E-16,	8.09746772196737E-16,	8.09747951207949E-16,	8.09749090865895E-16,	8.09750192469899E-16,	8.09751257276862E-16,	8.0975228650264E-16,	8.09753281323359E-16,	8.09754242876705E-16,	8.09755172263164E-16,	8.09756070547224E-16,	8.09756938758536E-16,	8.09757777893048E-16,	8.09758588914089E-16,	8.09759372753426E-16,	8.09760130312291E-16,	8.09760862462363E-16,	8.09761570046732E-16,	8.09762253880822E-16,	8.0976291475329E-16,	8.09763553426896E-16,	8.09764170639339E-16,	8.09764767104078E-16,	8.09765343511112E-16,	8.09765900527748E-16,	8.09766438799335E-16,	8.09766958949982E-16,	8.09767461583246E-16,	8.09767947282802E-16,	8.09768416613093E-16,	8.0976887011995E-16,	8.09769308331206E-16,	8.0976973175728E-16,	8.09770140891742E-16,	8.09770536211868E-16,	8.09770918179169E-16,	8.09771287239904E-16,	8.09771643825584E-16,	8.09771988353447E-16,	8.0977232122693E-16,	8.09772642836119E-16,	8.09772953558184E-16,	8.09773253757804E-16,	8.09773543787576E-16,	8.09773823988408E-16,	8.09774094689906E-16,	8.09774356210742E-16,	8.09774608859013E-16,	8.09774852932588E-16,	8.09775088719446E-16,	8.09775316497998E-16,	8.097755365374E-16,	8.09775749097861E-16,	8.09775954430936E-16,	8.09776152779809E-16,	8.0977634437957E-16,	8.09776529457479E-16,	8.0977670823323E-16,	8.09776880919191E-16,	8.09777047720654E-16,	8.09777208836062E-16,	8.0977736445724E-16,	8.09777514769609E-16,	8.09777659952399E-16,	8.09777800178856E-16,	8.09777935616434E-16,	8.09778066426993E-16,	8.09778192766978E-16,	8.09778314787602E-16,	8.09778432635017E-16,	8.09778546450484E-16,	8.0977865637053E-16,	8.09778762527112E-16,	8.09778865047762E-16,	8.09778964055736E-16,	8.09779059670156E-16,	8.09779152006148E-16,	8.09779241174974E-16,	8.09779327284158E-16,	8.09779410437613E-16,	8.09779490735761E-16,	8.09779568275648E-16,	8.09779643151054E-16,	8.09779715452606E-16,	8.09779785267882E-16,	8.09779852681509E-16,	8.09779917775265E-16,	8.09779980628172E-16,	8.0978004131659E-16,	8.09780099914301E-16,	8.09780156492602E-16,	8.09780211120381E-16,	8.09780263864203E-16,	8.09780314788383E-16,	8.09780363955064E-16,	8.09780411424291E-16,	8.09780457254075E-16,	8.09780501500468E-16,	8.09780544217624E-16,	8.09780585457866E-16,	8.09780625271743E-16,	8.09780663708093E-16,	8.097807008141E-16,	8.09780736635347E-16,	8.09780771215871E-16,	8.09780804598217E-16,	8.09780836823486E-16,	8.09780867931382E-16,	8.09780897960262E-16,	8.09780926947179E-16,	8.09780954927927E-16,	8.09780981937084E-16,	8.09781008008051E-16,	8.09781033173093E-16,	8.09781057463376E-16,	8.09781080909005E-16,	8.0978110353906E-16,	8.0978112538163E-16,	8.09781146463847E-16,	8.09781166811916E-16,	8.0978118645115E-16,	8.09781205405998E-16,	8.09781223700072E-16,	8.09781241356181E-16,	8.09781258396353E-16,	8.09781274841862E-16,	8.09781290713256E-16,	8.09781306030381E-16,	8.09781320812403E-16,	8.09781335077832E-16,	8.09781348844544E-16,	8.09781362129802E-16,	8.0978137495028E-16,	8.09781387322077E-16,	8.09781399260744E-16,	8.09781410781295E-16,	8.09781421898229E-16,	8.0978143262555E-16,	8.09781442976779E-16,	8.09781452964973E-16,	8.0978146260274E-16,	8.09781471902253E-16,	8.09781480875271E-16,	8.09781489533143E-16,	8.09781497886831E-16,	8.09781505946917E-16,	8.09781513723621E-16,	8.0978152122681E-16,	8.0978152846601E-16,	8.0978153545042E-16,	8.09781542188919E-16,	8.09781548690082E-16,	8.09781554962189E-16,	8.0978156101323E-16,	8.09781566850924E-16,	8.09781572482718E-16,	8.09781577915807E-16,	8.09781583157132E-16,	8.09781588213398E-16,	8.09781593091075E-16,	8.09781597796411E-16,	8.09781602335437E-16,	8.09781606713974E-16,	8.09781610937642E-16,	8.09781615011867E-16,	8.09781618941886E-16,	8.09781622732754E-16,	8.09781626389351E-16,	8.09781629916388E-16,	8.0978163331841E-16,	8.09781636599808E-16,	8.09781639764817E-16,	8.09781642817525E-16,	8.09781645761881E-16,	8.09781648601693E-16,	8.09781651340638E-16,	8.09781653982264E-16,	8.09781656529996E-16,	8.09781658987141E-16,	8.09781661356889E-16,	8.09781663642319E-16,	8.09781665846402E-16,	8.09781667972007E-16,	8.09781670021902E-16,	8.09781671998758E-16,	8.09781673905152E-16,	8.09781675743573E-16,	8.09781677516421E-16,	8.09781679226013E-16,	8.09781680874585E-16,	8.09781682464295E-16,	8.09781683997224E-16,	8.09781685475382E-16,	8.09781686900706E-16,	8.09781688275069E-16,	8.09781689600274E-16,	8.09781690878064E-16,	8.09781692110119E-16,	8.09781693298061E-16,	8.09781694443455E-16,	8.0978169554781E-16,	8.09781696612583E-16,	8.09781697639179E-16,	8.09781698628956E-16,	8.0978169958322E-16,	8.09781700503234E-16,	8.09781701390216E-16,	8.09781702245342E-16,	8.09781703069743E-16,	8.09781703864516E-16,	8.09781704630713E-16,	8.09781705369353E-16,	8.09781706081419E-16,	8.09781706767857E-16,	8.09781707429581E-16,	8.09781708067473E-16,	8.09781708682383E-16,	8.09781709275131E-16,	8.09781709846509E-16,	8.09781710397281E-16,	8.09781710928182E-16,	8.09781711439923E-16,	8.09781711933189E-16,	8.09781712408642E-16,	8.09781712866917E-16,	8.0978171330863E-16,	8.09781713734375E-16,	8.09781714144723E-16,	8.09781714540226E-16,	8.09781714921416E-16,	8.09781715288807E-16,	8.09781715642893E-16,	8.09781715984152E-16,	8.09781716313043E-16,	8.09781716630013E-16,	8.09781716935488E-16,	8.09781717229881E-16,	8.09781717513592E-16,	8.09781717787004E-16,	8.09781718050487E-16,	8.09781718304399E-16,	8.09781718549085E-16,	8.09781718784876E-16,	8.09781719012093E-16,	8.09781719231045E-16,	8.0978171944203E-16,	8.09781719645335E-16,	8.09781719841238E-16,	8.09781720030005E-16,	8.09781720211893E-16,	8.09781720387152E-16,	8.09781720556021E-16,	8.09781720718731E-16,	8.09781720875503E-16,	8.09781721026554E-16,	8.09781721172089E-16,	8.09781721312309E-16,	8.09781721447407E-16,	8.09781721577566E-16,	8.09781721702967E-16,	8.09781721823781E-16,	8.09781721940175E-16,	8.0978172205231E-16,	8.09781722160339E-16,	8.09781722264412E-16,	8.09781722364672E-16,	8.09781722461258E-16,	8.09781722554303E-16,	8.09781722643935E-16,	8.0978172273028E-16,	8.09781722813456E-16,	8.09781722893579E-16,	8.09781722970759E-16,	8.09781723045104E-16,	8.09781723116718E-16,	8.09781723185699E-16,	8.09781723252143E-16,	8.09781723316143E-16,	8.09781723377788E-16,	8.09781723437165E-16,	8.09781723494354E-16,	8.09781723549438E-16,	8.09781723602492E-16,	8.0978172365359E-16,	8.09781723702805E-16,	8.09781723750204E-16,	8.09781723795854E-16,	8.0978172383982E-16,	8.09781723882162E-16,	8.0978172392294E-16,	8.09781723962212E-16,	8.09781724000033E-16,	8.09781724036455E-16,	8.09781724071531E-16,	8.09781724105309E-16,	8.09781724137836E-16,	8.0978172416916E-16,	8.09781724199323E-16,	8.0978172422837E-16,	8.0978172425634E-16,	8.09781724283273E-16,	8.09781724309207E-16,	8.09781724334179E-16,	8.09781724358225E-16,	8.09781724381378E-16,	8.09781724403672E-16,	8.09781724425137E-16,	8.09781724445806E-16,	8.09781724465705E-16,	8.09781724484865E-16,	8.09781724503312E-16,	8.09781724521074E-16,	8.09781724538173E-16,	8.09781724554637E-16,	8.09781724570487E-16,	8.09781724585747E-16,	8.09781724600438E-16,	8.09781724614581E-16,	8.09781724628197E-16,	8.09781724641305E-16,	8.09781724653924E-16,	8.09781724666072E-16,	8.09781724677766E-16,	8.09781724689023E-16,	8.0978172469986E-16,	8.09781724710292E-16,	8.09781724720334E-16,	8.0978172473E-16,	8.09781724739304E-16,	8.0978172474826E-16,	8.09781724756881E-16,	8.09781724765179E-16,	8.09781724773166E-16,	8.09781724780854E-16,	8.09781724788254E-16,	8.09781724795376E-16,	8.09781724802231E-16,	8.09781724808828E-16,	8.09781724815178E-16,	8.0978172482129E-16,	8.09781724827172E-16,	8.09781724832833E-16,	8.09781724838281E-16,	8.09781724843524E-16,	8.0978172484857E-16,	8.09781724853427E-16,	8.097817248581E-16,	8.09781724862598E-16,	8.09781724866926E-16,	8.09781724871091E-16,	8.09781724875099E-16,	8.09781724878956E-16,	8.09781724882667E-16,	8.09781724886239E-16,	8.09781724889675E-16,	8.09781724892982E-16,	8.09781724896164E-16,	8.09781724899226E-16,	8.09781724902173E-16,	8.09781724905008E-16,	8.09781724907735E-16,	8.0978172491036E-16,	8.09781724912885E-16,	8.09781724915314E-16,	8.09781724917652E-16,	8.09781724919901E-16,	8.09781724922065E-16,	8.09781724924147E-16,	8.0978172492615E-16,	8.09781724928077E-16,	8.09781724929931E-16,	8.09781724931715E-16,	8.09781724933431E-16,	8.09781724935081E-16,	8.09781724936669E-16,	8.09781724938197E-16,	8.09781724939667E-16,	8.09781724941081E-16,	8.0978172494244E-16,	8.09781724943749E-16,	8.09781724945007E-16,	8.09781724946218E-16,	8.09781724947382E-16,	8.09781724948502E-16,	8.0978172494958E-16,	8.09781724950616E-16,	8.09781724951612E-16,	8.09781724952571E-16,	8.09781724953493E-16,	8.09781724954381E-16,	8.09781724955234E-16,	8.09781724956054E-16,	8.09781724956843E-16,	8.09781724957602E-16,	8.09781724958332E-16,	8.09781724959034E-16,	8.0978172495971E-16,	8.09781724960359E-16,	8.09781724960984E-16,	8.09781724961584E-16,	8.09781724962162E-16,	8.09781724962717E-16,	8.09781724963252E-16,	8.09781724963765E-16,	8.09781724964259E-16,	8.09781724964734E-16,	8.09781724965191E-16,	8.0978172496563E-16,	8.09781724966053E-16,	8.09781724966459E-16,	8.0978172496685E-16,	8.09781724967225E-16,	8.09781724967586E-16,	8.09781724967933E-16,	8.09781724968267E-16,	8.09781724968588E-16,	8.09781724968897E-16,	8.09781724969194E-16,	8.09781724969479E-16,	8.09781724969754E-16,	8.09781724970017E-16,	8.09781724970271E-16,	8.09781724970515E-16,	8.09781724970749E-16,	8.09781724970975E-16,	8.09781724971191E-16,	8.097817249714E-16,	8.097817249716E-16,	8.09781724971792E-16,	8.09781724971977E-16,	8.09781724972155E-16,	8.09781724972326E-16,	8.09781724972491E-16,	8.09781724972648E-16,	8.09781724972801E-16,	8.09781724972947E-16,	8.09781724973087E-16,	8.09781724973222E-16,	8.09781724973351E-16,	8.09781724973476E-16,	8.09781724973596E-16,	8.09781724973711E-16,	8.09781724973822E-16,	8.09781724973928E-16,	8.09781724974031E-16,	8.09781724974129E-16,	8.09781724974223E-16,	8.09781724974314E-16,	8.09781724974402E-16,	8.09781724974485E-16,	8.09781724974566E-16,	8.09781724974643E-16,	8.09781724974718E-16,	8.0978172497479E-16,	8.09781724974858E-16,	8.09781724974925E-16,	8.09781724974988E-16,	8.09781724975049E-16,	8.09781724975108E-16,	8.09781724975164E-16,	8.09781724975218E-16,	8.0978172497527E-16,	8.0978172497532E-16,	8.09781724975368E-16,	8.09781724975414E-16,	8.09781724975459E-16,	8.09781724975501E-16,	8.09781724975542E-16,	8.09781724975582E-16,	8.0978172497562E-16,	8.09781724975656E-16,	8.09781724975691E-16,	8.09781724975724E-16,	8.09781724975757E-16,	8.09781724975788E-16,	8.09781724975817E-16,	8.09781724975846E-16,	8.09781724975873E-16,	8.097817249759E-16,	8.09781724975925E-16,	8.09781724975949E-16,	8.09781724975973E-16,	8.09781724975995E-16,	8.09781724976017E-16,	8.09781724976038E-16,	8.09781724976058E-16,	8.09781724976077E-16,	8.09781724976095E-16,	8.09781724976113E-16,	8.0978172497613E-16,	8.09781724976146E-16,	8.09781724976162E-16,	8.09781724976177E-16,	8.09781724976191E-16,	8.09781724976205E-16,	8.09781724976218E-16,	8.09781724976231E-16,	8.09781724976244E-16,	8.09781724976255E-16,	8.09781724976267E-16,	8.09781724976278E-16,	8.09781724976288E-16,	8.09781724976298E-16,	8.09781724976308E-16,	8.09781724976317E-16,	8.09781724976326E-16,	8.09781724976335E-16,	8.09781724976343E-16,	8.09781724976351E-16,	8.09781724976358E-16,	8.09781724976365E-16,	8.09781724976373E-16,	8.09781724976379E-16,	8.09781724976386E-16,	8.09781724976392E-16,	8.09781724976398E-16,	8.09781724976404E-16,	8.09781724976409E-16,	8.09781724976414E-16,	8.09781724976419E-16,	8.09781724976424E-16,	8.09781724976429E-16,	8.09781724976433E-16,	8.09781724976438E-16,	8.09781724976442E-16,	8.09781724976446E-16,	8.0978172497645E-16,	8.09781724976453E-16,	8.09781724976457E-16,	8.0978172497646E-16,	8.09781724976463E-16,	8.09781724976467E-16,	8.0978172497647E-16,	8.09781724976472E-16,	8.09781724976475E-16,	8.09781724976478E-16,	8.0978172497648E-16,	8.09781724976483E-16,	8.09781724976485E-16,	8.09781724976487E-16,	8.0978172497649E-16,	8.09781724976492E-16,	8.09781724976493E-16,	8.09781724976495E-16,	8.09781724976497E-16,	8.09781724976499E-16,	8.09781724976501E-16,	8.09781724976502E-16,	8.09781724976504E-16,	8.09781724976505E-16,	8.09781724976507E-16,	8.09781724976508E-16,	8.09781724976509E-16,	8.09781724976511E-16,	8.09781724976512E-16,	8.09781724976513E-16,	8.09781724976514E-16,	8.09781724976515E-16,	8.09781724976516E-16,	8.09781724976517E-16,	8.09781724976518E-16,	8.09781724976519E-16,	8.0978172497652E-16,	8.09781724976521E-16,	8.09781724976522E-16,	8.09781724976522E-16,	8.09781724976523E-16,	8.09781724976524E-16,	8.09781724976525E-16,	8.09781724976525E-16,	8.09781724976526E-16,	8.09781724976526E-16,	8.09781724976527E-16,	8.09781724976528E-16,	8.09781724976528E-16,	8.09781724976529E-16,	8.09781724976529E-16,	8.09781724976529E-16,	8.0978172497653E-16,	8.0978172497653E-16,	8.09781724976531E-16,	8.09781724976531E-16,	8.09781724976532E-16,	8.09781724976532E-16,	8.09781724976532E-16,	8.09781724976533E-16,	8.09781724976533E-16,	8.09781724976533E-16,	8.09781724976534E-16,	8.09781724976534E-16,	8.09781724976534E-16,	8.09781724976534E-16,	8.09781724976535E-16,	8.09781724976535E-16,	8.09781724976535E-16,	8.09781724976536E-16,	8.09781724976536E-16,	8.09781724976536E-16,	8.09781724976536E-16,	8.09781724976536E-16,	8.09781724976536E-16,	8.09781724976536E-16,	8.09781724976537E-16,	8.09781724976537E-16,	8.09781724976537E-16,	8.09781724976537E-16,	8.09781724976537E-16,	8.09781724976537E-16,	8.09781724976537E-16,	8.09781724976538E-16,	8.09781724976538E-16,	8.09781724976538E-16,	8.09781724976538E-16,	8.09781724976538E-16,	8.09781724976538E-16,	8.09781724976538E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16,	8.09781724976539E-16};
    
	double force_1 = Force_1[timestep];
	double force_2 = -1*Force_1[timestep+5];

	double force_coefficient =     3538944;

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

