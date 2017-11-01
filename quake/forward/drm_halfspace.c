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
#include <stdio.h>

#include "psolve.h"
#include "octor.h"
#include "util.h"
#include "stiffness.h"
#include "quake_util.h"
#include "cvm.h"
#include "drm_halfspace.h"
#include "topography.h"
#include "geometrics.h"

static planewavetype_t	thePlaneWaveType;
static int32_t	        theDRMBox_halfwidthElements = 0;
static int32_t	        theDRMBox_DepthElements = 0;
static double 	        thedrmbox_esize         = 0.0;

static double 	        theTs = 0.0;
static double 	        thefc = 0.0;
static double           theUo = 0.0;
static double 	        theplanewave_strike = 0.0;
static double 	        theXc  = 0.0;
static double 	        theYc  = 0.0;

static int32_t            *myDRMFace1ElementsMapping;
static int32_t            *myDRMFace2ElementsMapping;
static int32_t            *myDRMFace3ElementsMapping;
static int32_t            *myDRMFace4ElementsMapping;
static int32_t            *myDRMBottomElementsMapping;

static int32_t            *myDRMBorder1ElementsMapping;
static int32_t            *myDRMBorder2ElementsMapping;
static int32_t            *myDRMBorder3ElementsMapping;
static int32_t            *myDRMBorder4ElementsMapping;

static int32_t            *myDRMBorder5ElementsMapping;
static int32_t            *myDRMBorder6ElementsMapping;
static int32_t            *myDRMBorder7ElementsMapping;
static int32_t            *myDRMBorder8ElementsMapping;
static double              theDRMdepth;

static int32_t            myDRM_Face1Count  = 0;
static int32_t            myDRM_Face2Count  = 0;
static int32_t            myDRM_Face3Count  = 0;
static int32_t            myDRM_Face4Count  = 0;
static int32_t            myDRM_BottomCount = 0;
static int32_t            myDRM_Brd1 = 0;
static int32_t            myDRM_Brd2 = 0;
static int32_t            myDRM_Brd3 = 0;
static int32_t            myDRM_Brd4 = 0;
static int32_t            myDRM_Brd5 = 0;
static int32_t            myDRM_Brd6 = 0;
static int32_t            myDRM_Brd7 = 0;
static int32_t            myDRM_Brd8 = 0;


void drm_planewaves_init ( int32_t myID, const char *parametersin ) {

    int     int_message[3];
    double  double_message[7];

    /* Capturing data from file --- only done by PE0 */
    if (myID == 0) {
        if ( drm_planewaves_initparameters( parametersin ) != 0 ) {
            fprintf(stderr,"Thread %d: drm_planewaves_init: "
                    "incidentPlaneWaves_initparameters error\n",myID);
            MPI_Abort(MPI_COMM_WORLD, ERROR);
            exit(1);
        }

    }

    /* Broadcasting data */

    int_message   [0]    = (int)thePlaneWaveType;
    int_message   [1]    = theDRMBox_halfwidthElements;
    int_message   [2]    = theDRMBox_DepthElements;

    double_message[0] = theTs;
    double_message[1] = thefc;
    double_message[2] = theUo;
    double_message[3] = theplanewave_strike;
    double_message[4] = theXc;
    double_message[5] = theYc;
    double_message[6] = thedrmbox_esize;

    MPI_Bcast(double_message, 7, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(int_message,    3, MPI_INT,    0, comm_solver);

    thePlaneWaveType             = int_message[0];
    theDRMBox_halfwidthElements  = int_message[1];
    theDRMBox_DepthElements      = int_message[2];

    theTs               = double_message[0];
    thefc               = double_message[1];
    theUo               = double_message[2];
    theplanewave_strike = double_message[3];
    theXc               = double_message[4];
    theYc               = double_message[5];
    thedrmbox_esize     = double_message[6];

    return;

}



int32_t
drm_planewaves_initparameters ( const char *parametersin ) {
    FILE                *fp;

    double              drmbox_halfwidth_elements, drmbox_depth_elements, Ts, fc, Uo, planewave_strike, L_ew, L_ns, drmbox_esize;
    char                type_of_wave[64];

    planewavetype_t     planewave;


    /* Opens parametersin file */

   if ( ( fp = fopen(parametersin, "r" ) ) == NULL ) {
        fprintf( stderr,
                 "Error opening %s\n at drm_planewaves_initparameters",
                 parametersin );
        return -1;
    }


     /* Parses parametersin to capture drm_planewaves single-value parameters */
    if ( ( parsetext(fp, "type_of_wave",                     's', &type_of_wave                ) != 0) ||
         ( parsetext(fp, "DRMBox_Noelements_in_halfwidth",   'd', &drmbox_halfwidth_elements   ) != 0) ||
         ( parsetext(fp, "DRMBox_Noelements_in_depth",       'd', &drmbox_depth_elements       ) != 0) ||
         ( parsetext(fp, "DRMBox_element_size_m",            'd', &drmbox_esize                ) != 0) ||
         ( parsetext(fp, "Ts",                               'd', &Ts                          ) != 0) ||
         ( parsetext(fp, "region_length_east_m",             'd', &L_ew                        ) != 0) ||
         ( parsetext(fp, "region_length_north_m",            'd', &L_ns                        ) != 0) ||
         ( parsetext(fp, "fc",                               'd', &fc                          ) != 0) ||
         ( parsetext(fp, "Uo",                               'd', &Uo                          ) != 0) ||
         ( parsetext(fp, "planewave_strike",                 'd', &planewave_strike            ) != 0) )
    {
        fprintf( stderr,
                 "Error parsing planewaves parameters from %s\n",
                 parametersin );
        return -1;
    }

    if ( strcasecmp(type_of_wave, "SV") == 0 ) {
    	planewave = SV;
    } else if ( strcasecmp(type_of_wave, "P") == 0 ) {
    	planewave = P;
    } else if ( strcasecmp(type_of_wave, "SINEWAVE") == 0 ) {
    	planewave = SINEWAVE;
    }else {
        fprintf(stderr,
                "Illegal type_of_wave for incident plane wave analysis"
                "(SV, P, SINE): %s\n", type_of_wave);
        return -1;
    }

    /*  Initialize the static global variables */
	thePlaneWaveType                 = planewave;
	theDRMBox_halfwidthElements      = drmbox_halfwidth_elements;
	theDRMBox_DepthElements          = drmbox_depth_elements;
	theTs                            = Ts;
	thefc                            = fc;
    theUo                            = Uo;
	theplanewave_strike              = planewave_strike * PI / 180.00;
	theXc                            = L_ew / 2.0;
	theYc                            = L_ns / 2.0;
	thedrmbox_esize                  = drmbox_esize;

    fclose(fp);

    return 0;
}


void PlaneWaves_solver_init( int32_t myID, mesh_t *myMesh, mysolver_t *mySolver) {

	int32_t theFaceElem, theBaseElem;

	double DRM_D = theDRMBox_DepthElements * thedrmbox_esize;
	double DRM_B = theDRMBox_halfwidthElements * thedrmbox_esize;
	double thebase_zcoord = get_thebase_topo();

	theFaceElem = 2 * ( theDRMBox_halfwidthElements + 0 ) * theDRMBox_DepthElements;
	theBaseElem = 4 * theDRMBox_halfwidthElements * theDRMBox_halfwidthElements;
	theDRMdepth = DRM_D;

	/*  mapping of face1 elements */
	int32_t eindex;
	int32_t countf1 = 0, countf2 = 0, countf3 = 0, countf4 = 0, countbott=0;
	int32_t countb1 = 0, countb2 = 0, countb3 = 0, countb4 = 0;
	int32_t countb5 = 0, countb6 = 0, countb7 = 0, countb8 = 0;

	XMALLOC_VAR_N(myDRMFace1ElementsMapping, int32_t, theFaceElem);
	XMALLOC_VAR_N(myDRMFace2ElementsMapping, int32_t, theFaceElem);
	XMALLOC_VAR_N(myDRMFace3ElementsMapping, int32_t, theFaceElem);
	XMALLOC_VAR_N(myDRMFace4ElementsMapping, int32_t, theFaceElem);
	XMALLOC_VAR_N(myDRMBottomElementsMapping , int32_t, theBaseElem);

	/* border elements*/
	XMALLOC_VAR_N(myDRMBorder1ElementsMapping, int32_t, theDRMBox_DepthElements + 1);
	XMALLOC_VAR_N(myDRMBorder2ElementsMapping, int32_t, theDRMBox_DepthElements + 1);
	XMALLOC_VAR_N(myDRMBorder3ElementsMapping, int32_t, theDRMBox_DepthElements + 1);
	XMALLOC_VAR_N(myDRMBorder4ElementsMapping, int32_t, theDRMBox_DepthElements + 1);

	XMALLOC_VAR_N(myDRMBorder5ElementsMapping, int32_t, theFaceElem);
	XMALLOC_VAR_N(myDRMBorder6ElementsMapping, int32_t, theFaceElem);
	XMALLOC_VAR_N(myDRMBorder7ElementsMapping, int32_t, theFaceElem);
	XMALLOC_VAR_N(myDRMBorder8ElementsMapping, int32_t, theFaceElem);

	for (eindex = 0; eindex < myMesh->lenum; eindex++) {

		elem_t     *elemp;
		node_t     *node0dat;
		edata_t    *edata;
		double      xo, yo, zo;
		int32_t	    node0;

		elemp    = &myMesh->elemTable[eindex]; //Takes the information of the "eindex" element
		node0    = elemp->lnid[0];             //Takes the ID for the zero node in element eindex
		node0dat = &myMesh->nodeTable[node0];
		edata    = (edata_t *)elemp->data;

		/* get coordinates of element zero node */
		xo = (node0dat->x)*(myMesh->ticksize);
		yo = (node0dat->y)*(myMesh->ticksize);
		zo = (node0dat->z)*(myMesh->ticksize);


		if ( ( ( yo - theYc ) == DRM_B ) &&                            /*face 1*/
				( xo >= ( theXc - DRM_B ) ) &&
				( xo <  ( theXc + DRM_B ) ) &&
				( zo <  DRM_D + thebase_zcoord ) &&
				( zo >=  thebase_zcoord ) ) {

			myDRMFace1ElementsMapping[countf1] = eindex;
			countf1++;
		} else 	if ( ( ( theYc - ( yo + thedrmbox_esize ) ) == DRM_B ) &&        /* face 2*/
				( xo >= ( theXc - DRM_B ) ) &&
				( xo <  ( theXc + DRM_B ) ) &&
				( zo <  DRM_D + thebase_zcoord ) &&
				( zo >=  thebase_zcoord ) ) {

			myDRMFace2ElementsMapping[countf2] = eindex;
			countf2++;

		} else 	if ( ( ( theXc - ( xo + thedrmbox_esize ) ) == DRM_B ) &&      /*face 3*/
				( yo >= ( theYc - DRM_B ) ) &&
				( yo <  ( theYc + DRM_B ) ) &&
				( zo <  DRM_D + thebase_zcoord ) &&
				( zo >=  thebase_zcoord ) ) {

			myDRMFace3ElementsMapping[countf3] = eindex;
			countf3++;

		} else 	if ( ( ( xo - theXc ) == DRM_B ) &&             /* face 4*/
				( yo >= ( theYc - DRM_B ) ) &&
				( yo <  ( theYc + DRM_B ) ) &&
				( zo <  DRM_D + thebase_zcoord ) &&
				( zo >=  thebase_zcoord ) ) {

			myDRMFace4ElementsMapping[countf4] = eindex;
			countf4++;

		} else 	if ( ( yo >= ( theYc - DRM_B ) ) &&         /*bottom*/
				( yo <  ( theYc + DRM_B ) ) &&
				( xo >= ( theXc - DRM_B ) ) &&
				( xo <  ( theXc + DRM_B ) ) &&
				( zo ==  DRM_D + thebase_zcoord )          ) {

			myDRMBottomElementsMapping[countbott] = eindex;
			countbott++;

		} else 	if ( ( ( yo - theYc ) == DRM_B ) &&      /*border 1*/
				( ( xo + thedrmbox_esize ) == ( theXc - DRM_B ) ) &&
				( zo <=  DRM_D + thebase_zcoord ) &&
				( zo >=  thebase_zcoord ) ) {

			myDRMBorder1ElementsMapping[countb1] = eindex;
			countb1++;

		} else if ( ( ( yo - theYc ) == DRM_B  ) &&      /*border 2*/
				( xo  == ( theXc + DRM_B ) ) &&
				( zo <=  DRM_D + thebase_zcoord ) &&
				( zo >=  thebase_zcoord ) ) {

			myDRMBorder2ElementsMapping[countb2] = eindex;
			countb2++;

		} else if ( ( ( theYc - ( yo + thedrmbox_esize ) ) == DRM_B ) &&  /* border 3*/
				( ( xo + thedrmbox_esize)  == ( theXc - DRM_B ) ) &&
				( zo <=  DRM_D + thebase_zcoord ) &&
				( zo >=  thebase_zcoord ) ) {

			myDRMBorder3ElementsMapping[countb3] = eindex;
			countb3++;

		} else if ( ( ( theYc - ( yo + thedrmbox_esize ) ) == DRM_B ) &&  /* border 4*/
				(  xo  == ( theXc + DRM_B ) ) &&
				( zo <=  DRM_D + thebase_zcoord ) &&
				( zo >=  thebase_zcoord ) ) {

			myDRMBorder4ElementsMapping[countb4] = eindex;
			countb4++;

		} else 	if ( ( ( yo - theYc ) == DRM_B ) &&            /* border 5*/
				( xo >= ( theXc - DRM_B ) ) &&
				( xo <  ( theXc + DRM_B ) ) &&
				( zo ==  DRM_D + thebase_zcoord )          ) {

			myDRMBorder5ElementsMapping[countb5] = eindex;
			countb5++;

		} else if ( ( ( theYc - ( yo + thedrmbox_esize ) ) == DRM_B ) &&          /* border 6*/
				( xo >= ( theXc - DRM_B ) ) &&
				( xo <  ( theXc + DRM_B ) ) &&
				( zo ==  DRM_D + thebase_zcoord )          ) {

			myDRMBorder6ElementsMapping[countb6] = eindex;
			countb6++;

		} else if ( ( ( theXc - ( xo + thedrmbox_esize ) ) == DRM_B ) &&      /* border 7*/
				( yo >= ( theYc - DRM_B ) ) &&
				( yo <  ( theYc + DRM_B ) ) &&
				( zo ==  DRM_D + thebase_zcoord )          ) {

			myDRMBorder7ElementsMapping[countb7] = eindex;
			countb7++;

		} else if ( ( ( xo - theXc ) == DRM_B ) &&             /* border 8*/
				( yo >= ( theYc - DRM_B ) ) &&
				( yo <  ( theYc + DRM_B ) ) &&
				( zo ==  DRM_D + thebase_zcoord )          ) {

			myDRMBorder8ElementsMapping[countb8] = eindex;
			countb8++;
		}
	}


	myDRM_Face1Count  = countf1;
	myDRM_Face2Count  = countf2;
	myDRM_Face3Count  = countf3;
	myDRM_Face4Count  = countf4;
	myDRM_BottomCount = countbott;
	myDRM_Brd1        = countb1;
	myDRM_Brd2        = countb2;
	myDRM_Brd3        = countb3;
	myDRM_Brd4        = countb4;
	myDRM_Brd5        = countb5;
	myDRM_Brd6        = countb6;
	myDRM_Brd7        = countb7;
	myDRM_Brd8        = countb8;

//	fprintf(stdout,"myID = %d, myDRM_Face1Count= %d, myDRM_Face2Count= %d, myDRM_Face3Count= %d, myDRM_Face4Count= %d, myDRM_BottomCount=%d \n"
//			       "myDRM_Brd1=%d, myDRM_Brd2=%d, myDRM_Brd3=%d, myDRM_Brd4=%d, myDRM_Brd5=%d, myDRM_Brd6=%d, myDRM_Brd7=%d, myDRM_Brd8=%d \n\n",
//			       myID, countf1, countf2, countf3, countf4,countbott,countb1,countb2,countb3,countb4,countb5,countb6,countb7,countb8);

}


void compute_addforce_PlaneWaves ( mesh_t     *myMesh,
                                mysolver_t *mySolver,
                                double      theDeltaT,
                                int         step,
                                fmatrix_t (*theK1)[8], fmatrix_t (*theK2)[8])
{

    int32_t   eindex;
    int32_t   face_eindex;

    double thebase_zcoord = get_thebase_topo();

    int  f_nodes_face1[4] = { 0, 1, 4, 5 };
    int  e_nodes_face1[4] = { 2, 3, 6, 7 };

    int  f_nodes_face2[4] = { 2, 3, 6, 7 };
    int  e_nodes_face2[4] = { 0, 1, 4, 5 };

    int  f_nodes_face3[4] = { 1, 3, 5, 7 };
    int  e_nodes_face3[4] = { 0, 2, 4, 6 };

    int  f_nodes_face4[4] = { 0, 2, 4, 6 };
    int  e_nodes_face4[4] = { 1, 3, 5, 7 };

    int  f_nodes_bottom[4] = { 0, 1, 2, 3 };
    int  e_nodes_bottom[4] = { 4, 5, 6, 7 };

    int  f_nodes_border1[2] = { 1, 5 };
    int  e_nodes_border1[6] = { 0, 2, 3, 4, 6, 7 };

    int  f_nodes_border2[2] = { 0, 4 };
    int  e_nodes_border2[6] = { 1, 2, 3, 5, 6, 7 };

    int  f_nodes_border3[2] = { 3, 7 };
    int  e_nodes_border3[6] = { 0, 1, 2, 4, 5, 6 };

    int  f_nodes_border4[2] = { 2, 6 };
    int  e_nodes_border4[6] = { 0, 1, 3, 4, 5, 7 };

    int  f_nodes_border5[2] = { 0, 1 };
    int  e_nodes_border5[6] = { 2, 3, 4, 5, 6, 7 };

    int  f_nodes_border6[2] = { 2, 3 };
    int  e_nodes_border6[6] = { 0, 1, 4, 5, 6, 7 };

    int  f_nodes_border7[2] = { 1, 3 };
    int  e_nodes_border7[6] = { 0, 2, 4, 5, 6, 7 };

    int  f_nodes_border8[2] = { 0, 2 };
    int  e_nodes_border8[6] = { 1, 3, 4, 5, 6, 7 };

    double tt = theDeltaT * step;
    int  *f_nodes, *e_nodes;

    /* Loop over face1 elements */
    f_nodes = &f_nodes_face1[0];
    e_nodes = &e_nodes_face1[0];
    for ( face_eindex = 0; face_eindex < myDRM_Face1Count ; face_eindex++) {
    	eindex = myDRMFace1ElementsMapping[face_eindex];
    	DRM_ForcesinElement ( myMesh, mySolver, theK1, theK2, f_nodes, e_nodes, eindex, tt, 4, 4 );
    } /* all elements in face 1*/

    /* Loop over face2 elements */
    f_nodes = &f_nodes_face2[0];
    e_nodes = &e_nodes_face2[0];
    for ( face_eindex = 0; face_eindex < myDRM_Face2Count ; face_eindex++) {
    	eindex = myDRMFace2ElementsMapping[face_eindex];
    	DRM_ForcesinElement ( myMesh, mySolver, theK1, theK2, f_nodes, e_nodes, eindex, tt, 4, 4 );
    } /* all elements in face 2*/

    /* Loop over face3 elements */
    f_nodes = &f_nodes_face3[0];
    e_nodes = &e_nodes_face3[0];
    for ( face_eindex = 0; face_eindex < myDRM_Face3Count ; face_eindex++) {
    	eindex = myDRMFace3ElementsMapping[face_eindex];
    	DRM_ForcesinElement ( myMesh, mySolver, theK1, theK2, f_nodes, e_nodes, eindex, tt, 4, 4 );
    } /* all elements in face 3*/

    /* Loop over face4 elements */
    f_nodes = &f_nodes_face4[0];
    e_nodes = &e_nodes_face4[0];
    for ( face_eindex = 0; face_eindex < myDRM_Face4Count ; face_eindex++) {
    	eindex = myDRMFace4ElementsMapping[face_eindex];
    	DRM_ForcesinElement ( myMesh, mySolver, theK1, theK2, f_nodes, e_nodes, eindex, tt, 4, 4 );
    } /* all elements in face 4*/

    /* Loop over bottom elements */
    f_nodes = &f_nodes_bottom[0];
    e_nodes = &e_nodes_bottom[0];
    for ( face_eindex = 0; face_eindex < myDRM_BottomCount ; face_eindex++) {
    	eindex = myDRMBottomElementsMapping[face_eindex];
    	DRM_ForcesinElement ( myMesh, mySolver, theK1, theK2, f_nodes, e_nodes, eindex, tt, 4, 4 );
    } /* all elements in bottom*/

    /* Loop over border1 elements */
    f_nodes = &f_nodes_border1[0];
    e_nodes = &e_nodes_border1[0];
    for ( face_eindex = 0; face_eindex < myDRM_Brd1 ; face_eindex++) {
    	eindex = myDRMBorder1ElementsMapping[face_eindex];

    	/* check for bottom element */
    	elem_t        *elemp;
    	int32_t	      node0;
    	node_t        *node0dat;
    	double        zo;

    	elemp        = &myMesh->elemTable[eindex];
    	node0        = elemp->lnid[0];
    	node0dat     = &myMesh->nodeTable[node0];
    	zo           = (node0dat->z)*(myMesh->ticksize);

    	if ( zo != theDRMdepth + thebase_zcoord ) {
    		DRM_ForcesinElement ( myMesh, mySolver, theK1, theK2, f_nodes, e_nodes, eindex, tt, 6, 2 );
    	} else {
    	    int  f_corner[1] = { 1 };
    	    int  e_corner[7] = { 0, 2, 3, 4, 5, 6, 7 };
    	    int  *fcorner_nodes, *ecorner_nodes;
    	    fcorner_nodes = &f_corner[0];
    	    ecorner_nodes = &e_corner[0];
    	    DRM_ForcesinElement ( myMesh, mySolver, theK1, theK2, fcorner_nodes, ecorner_nodes, eindex, tt, 7, 1 );

    	}
    } /* all elements in border1*/

    /* Loop over border2 elements */
    f_nodes = &f_nodes_border2[0];
    e_nodes = &e_nodes_border2[0];
    for ( face_eindex = 0; face_eindex < myDRM_Brd2 ; face_eindex++) {
    	eindex = myDRMBorder2ElementsMapping[face_eindex];

    	/* check for bottom element */
    	elem_t        *elemp;
    	int32_t	      node0;
    	node_t        *node0dat;
    	double        zo;

    	elemp        = &myMesh->elemTable[eindex];
    	node0        = elemp->lnid[0];
    	node0dat     = &myMesh->nodeTable[node0];
    	zo           = (node0dat->z)*(myMesh->ticksize);

    	if ( zo != theDRMdepth + thebase_zcoord ) {
    		DRM_ForcesinElement ( myMesh, mySolver, theK1, theK2, f_nodes, e_nodes, eindex, tt, 6, 2 );
    	} else {
    	    int  f_corner[1] = { 0 };
    	    int  e_corner[7] = { 1, 2, 3, 4, 5, 6, 7 };
    	    int  *fcorner_nodes, *ecorner_nodes;
    	    fcorner_nodes = &f_corner[0];
    	    ecorner_nodes = &e_corner[0];
    	    DRM_ForcesinElement ( myMesh, mySolver, theK1, theK2, fcorner_nodes, ecorner_nodes, eindex, tt, 7, 1 );

    	}

    } /* all elements in border2*/

    /* Loop over border3 elements */
    f_nodes = &f_nodes_border3[0];
    e_nodes = &e_nodes_border3[0];
    for ( face_eindex = 0; face_eindex < myDRM_Brd3 ; face_eindex++) {
    	eindex = myDRMBorder3ElementsMapping[face_eindex];

    	/* check for bottom element */
    	elem_t        *elemp;
    	int32_t	      node0;
    	node_t        *node0dat;
    	double        zo;

    	elemp        = &myMesh->elemTable[eindex];
    	node0        = elemp->lnid[0];
    	node0dat     = &myMesh->nodeTable[node0];
    	zo           = (node0dat->z)*(myMesh->ticksize);

    	if ( zo != theDRMdepth + thebase_zcoord ) {
    		DRM_ForcesinElement ( myMesh, mySolver, theK1, theK2, f_nodes, e_nodes, eindex, tt, 6, 2 );
    	} else {
    	    int  f_corner[1] = { 3 };
    	    int  e_corner[7] = { 0, 1, 2, 4, 5, 6, 7 };
    	    int  *fcorner_nodes, *ecorner_nodes;
    	    fcorner_nodes = &f_corner[0];
    	    ecorner_nodes = &e_corner[0];
    	    DRM_ForcesinElement ( myMesh, mySolver, theK1, theK2, fcorner_nodes, ecorner_nodes, eindex, tt, 7, 1 );

    	}
    } /* all elements in border3*/

    /* Loop over border4 elements */
    f_nodes = &f_nodes_border4[0];
    e_nodes = &e_nodes_border4[0];
    for ( face_eindex = 0; face_eindex < myDRM_Brd4 ; face_eindex++) {
    	eindex = myDRMBorder4ElementsMapping[face_eindex];

    	/* check for bottom element */
    	elem_t        *elemp;
    	int32_t	      node0;
    	node_t        *node0dat;
    	double        zo;

    	elemp        = &myMesh->elemTable[eindex];
    	node0        = elemp->lnid[0];
    	node0dat     = &myMesh->nodeTable[node0];
    	zo           = (node0dat->z)*(myMesh->ticksize);

    	if ( zo != theDRMdepth + thebase_zcoord ) {
    		DRM_ForcesinElement ( myMesh, mySolver, theK1, theK2, f_nodes, e_nodes, eindex, tt, 6, 2 );
    	} else {
    	    int  f_corner[1] = { 3 };
    	    int  e_corner[7] = { 0, 1, 2, 4, 5, 6, 7 };
    	    int  *fcorner_nodes, *ecorner_nodes;
    	    fcorner_nodes = &f_corner[0];
    	    ecorner_nodes = &e_corner[0];
    	    DRM_ForcesinElement ( myMesh, mySolver, theK1, theK2, fcorner_nodes, ecorner_nodes, eindex, tt, 7, 1 );

    	}
    } /* all elements in border4*/

    /* Loop over border5 elements */
    f_nodes = &f_nodes_border5[0];
    e_nodes = &e_nodes_border5[0];
    for ( face_eindex = 0; face_eindex < myDRM_Brd5 ; face_eindex++) {
    	eindex = myDRMBorder5ElementsMapping[face_eindex];
    	DRM_ForcesinElement ( myMesh, mySolver, theK1, theK2, f_nodes, e_nodes, eindex, tt, 6, 2 );
    } /* all elements in border5*/

    /* Loop over border6 elements */
    f_nodes = &f_nodes_border6[0];
    e_nodes = &e_nodes_border6[0];
    for ( face_eindex = 0; face_eindex < myDRM_Brd6 ; face_eindex++) {
    	eindex = myDRMBorder6ElementsMapping[face_eindex];
    	DRM_ForcesinElement ( myMesh, mySolver, theK1, theK2, f_nodes, e_nodes, eindex, tt, 6, 2 );
    } /* all elements in border6*/

    /* Loop over border7 elements */
    f_nodes = &f_nodes_border7[0];
    e_nodes = &e_nodes_border7[0];
    for ( face_eindex = 0; face_eindex < myDRM_Brd7 ; face_eindex++) {
    	eindex = myDRMBorder7ElementsMapping[face_eindex];
    	DRM_ForcesinElement ( myMesh, mySolver, theK1, theK2, f_nodes, e_nodes, eindex, tt, 6, 2 );
    } /* all elements in border7*/

    /* Loop over border8 elements */
    f_nodes = &f_nodes_border8[0];
    e_nodes = &e_nodes_border8[0];
    for ( face_eindex = 0; face_eindex < myDRM_Brd8 ; face_eindex++) {
    	eindex = myDRMBorder8ElementsMapping[face_eindex];
    	DRM_ForcesinElement ( myMesh, mySolver, theK1, theK2, f_nodes, e_nodes, eindex, tt, 6, 2 );
    } /* all elements in border7*/


    return;
}


void DRM_ForcesinElement ( mesh_t     *myMesh,
		mysolver_t *mySolver,
		fmatrix_t (*theK1)[8], fmatrix_t (*theK2)[8],
		int *f_nodes, int *e_nodes, int32_t   eindex, double tt, int Nnodes_e, int Nnodes_f )
{

    int       i, j;
    int  CoordArr[8]      = { 0, 0, 0, 0, 1, 1, 1, 1 };

    double thebase_zcoord = get_thebase_topo();

	fvector_t localForce[8];

	elem_t        *elemp;
	edata_t       *edata;
	node_t        *node0dat;
	double        xo, yo, zo;
	int32_t	      node0;
	e_t*          ep;

	/* Capture the table of elements from the mesh and the size
	 * This is what gives me the connectivity to nodes */
	elemp        = &myMesh->elemTable[eindex];
	edata        = (edata_t *)elemp->data;
	node0        = elemp->lnid[0];
	node0dat     = &myMesh->nodeTable[node0];
	ep           = &mySolver->eTable[eindex];

	/* get coordinates of element zero node */
	xo = (node0dat->x)*(myMesh->ticksize);
	yo = (node0dat->y)*(myMesh->ticksize);
	zo = (node0dat->z)*(myMesh->ticksize) - thebase_zcoord;


	/* get material properties  */
	double  h, Vs;
	h    = (double)edata->edgesize;

	if ( thePlaneWaveType == SV  )
		Vs = edata->Vs;
	else
		Vs = edata->Vp;


	/* Force contribution from external nodes */
	/* -------------------------------
	 * Ku DONE IN THE CONVENTIONAL WAY
	 * ------------------------------- */
	memset( localForce, 0, 8 * sizeof(fvector_t) );

	fvector_t myDisp;
	/* forces over f nodes */
	for (i = 0; i < Nnodes_f; i++) {

		int  nodef = *(f_nodes + i);
		fvector_t* toForce = &localForce[ nodef ];

		/* incoming displacements over e nodes */
		for (j = 0; j < Nnodes_e; j++) {

			int  nodee = *(e_nodes + j);

			double z_ne = zo + h * CoordArr[ nodee ];   /* get zcoord */
			getRicker ( &myDisp, z_ne, tt, Vs ); /* get Displ */

			MultAddMatVec( &theK1[ nodef ][ nodee ], &myDisp, -ep->c1, toForce );
			MultAddMatVec( &theK2[ nodef ][ nodee ], &myDisp, -ep->c2, toForce );
		}
	}

	/* forces over e nodes */
	for (i = 0; i < Nnodes_e; i++) {

		int  nodee = *(e_nodes + i);
		fvector_t* toForce = &localForce[ nodee ];

		/* incoming displacements over f nodes */
		for (j = 0; j < Nnodes_f; j++) {

			int  nodef = *(f_nodes + j);

			double z_nf = zo + h * CoordArr[ nodef ];   /* get zcoord */
			getRicker ( &myDisp, z_nf, tt, Vs ); /* get Displ */

			MultAddMatVec( &theK1[ nodee ][ nodef ], &myDisp, ep->c1, toForce );
			MultAddMatVec( &theK2[ nodee ][ nodef ], &myDisp, ep->c2, toForce );
		}
	}
	/* end Ku */

	/* Loop over the 8 element nodes:
	 * Add the contribution calculated above to the node
	 * forces carried from the source and stiffness.
	 */
	for (i = 0; i < 8; i++) {

		int32_t    lnid;
		fvector_t *nodalForce;

		lnid = elemp->lnid[i];

		nodalForce = mySolver->force + lnid;

		nodalForce->f[0] += localForce[i].f[0] ;
		nodalForce->f[1] += localForce[i].f[1] ;
		nodalForce->f[2] += localForce[i].f[2] ;

	} /* element nodes */

}

void getRicker ( fvector_t *myDisp, double zp, double t, double Vs ) {

	double Rz = Ricker_displ ( zp, theTs, t, thefc, Vs  );


	if (thePlaneWaveType == P || thePlaneWaveType == SV ){
		if ( thePlaneWaveType == SV ) {
			myDisp->f[0] = Rz * theUo * cos (theplanewave_strike);
			myDisp->f[1] = Rz * theUo * sin (theplanewave_strike);
			myDisp->f[2] = 0.0;
		} else {
			myDisp->f[0] = 0.0;
			myDisp->f[1] = 0.0;
			myDisp->f[2] = Rz * theUo;
		}
	}else if (thePlaneWaveType == SINEWAVE){
        double f = 2;
		double amp = t*t*t*(0.000125)*sin(2*PI*t*f);
		myDisp->f[0] = 1 * theUo * cos (theplanewave_strike) * amp;
		myDisp->f[1] = 1 * theUo * sin (theplanewave_strike) * amp;
		myDisp->f[2] = 0.0;

	}
}

double Ricker_displ ( double zp, double Ts, double t, double fc, double Vs  ) {

	double alfa1 = ( PI * fc ) * ( PI * fc ) * ( t - zp / Vs - Ts) * ( t - zp / Vs - Ts);
	double alfa2 = ( PI * fc ) * ( PI * fc ) * ( t + zp / Vs - Ts) * ( t + zp / Vs - Ts);

	double uo1 = ( 2.0 * alfa1 - 1.0 ) * exp(-alfa1);
	double uo2 = ( 2.0 * alfa2 - 1.0 ) * exp(-alfa2);

	return (uo1+uo2);
}

