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

#ifndef Q_EQLINEAR_H
#define Q_EQLINEAR_H

/* -------------------------------------------------------------------------- */
/*                        Structures and definitions                          */
/* -------------------------------------------------------------------------- */

typedef struct eq_tensor_t {

    double xx;
    double yy;
    double zz;
    double xy;
    double yz;
    double xz;

} eq_tensor_t;

typedef struct eq_qpvectors_t {

    double qv[8];

} eq_qpvectors_t;

typedef struct eq_qptensors_t {

    eq_tensor_t qp[8];

} eq_qptensors_t;




typedef struct elconstants_t {

    double lambda;
    double mu;
    double original_mu;
    double Qs_value;
    double Qp_value;

//    double alpha;        /*  yield function constants in Drucker-Prager model*/
//    double beta;         /*  constant of the plastic potential flow law */
//    double gamma;        /*  constant for teh hardening function in Drucker-Prager model */
//
//    double c;            /* soil cohesion */
//    double phi;          /* angle of internal friction */
//    double dil_angle;    /* angle of dilatancy */
//
//    double k;
//    double h;             /*  variable used for the isotropic hardening function.
//                              vonMises H=0. Hardening is considered Kinematic in vonMises criterion
//                              In Drucker-Prager H = gamma(c + h*ep)  */
//                          /*  In MohrCoulomb    H =     2(c + h*ep)cos(phi)  */
//
//    double Sstrain0;      /* Defines the elastic range of the vonMises model. Sy=G*Sstrain0   */
//
//    double fs[8];         /* F(sigma) */
//    double dLambda[8];    /* yield control */
//    double strainrate;
//    double sensitivity;
//
//    double sigmaZ_st;    /* static vertical stress */
//
//    double maxFs;
//    double avgFs;

} elconstants_t;


typedef struct elsolver_t {

    elconstants_t    *constants;
    eq_qptensors_t   *stresses;
    eq_qptensors_t   *strains;
    eq_qptensors_t   *maxstrains;
    eq_qpvectors_t   *effectivestrain;
//    qptensors_t   *pstrains1;
//    qptensors_t   *pstrains2;
//    qptensors_t   *alphastress1;
//    qptensors_t   *alphastress2;
//    qpvectors_t   *ep1;         /* effective plastic strains */
//    qpvectors_t   *ep2;

} elsolver_t;


typedef struct GD_t {

    double g;
    double d;

} GD_t;

/* -------------------------------------------------------------------------- */
/*                                 Utilities                                  */
/* -------------------------------------------------------------------------- */
int isThisElementEqLinear(mesh_t *myMesh, int32_t eindex);


/* -------------------------------------------------------------------------- */
/*       Initialization of parameters, structures and memory allocations      */
/* -------------------------------------------------------------------------- */
void eqlinear_solver_init(int32_t myID, mesh_t *myMesh, double depth);
void eqlinear_stats(int32_t myID, int32_t theGroupSize);
void constract_GD_Table();
void compute_addforce_bottom(int32_t timestep, mesh_t *myMesh, mysolver_t *mySolver, double dt);
GD_t  search_GD_table(double strain);
void eqlinear_init ( int32_t     myID,
                      const char *parametersin,
                      double      theDeltaT,
                      double      theEndT );
/* -------------------------------------------------------------------------- */
/*                   Auxiliary tensor manipulation methods                    */
/* -------------------------------------------------------------------------- */

eq_tensor_t point_strain_eq    ( fvector_t *u, double lx, double ly, double lz,
                           double h);
eq_tensor_t init_tensor_eq     ( );
int get_displacements_eq ( mysolver_t *solver, elem_t *elemp, fvector_t *u );

/* -------------------------------------------------------------------------- */
/*                              Stability methods                             */
/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */
/*                   Nonlinear core computational methods                     */
/* -------------------------------------------------------------------------- */

void compute_eqlinear_state ( mesh_t     *myMesh,
                               mysolver_t *mySolver,
                               int32_t     theNumberOfStations,
                               int32_t     myNumberOfStations,
                               station_t  *myStations,
                               double      theDeltaT,
                               int         step,
							   double      theVSmaxeq,
							   double      theEQA,
							   double      theEQB,
							   double      theEQC,
							   double      theEQD,
							   double      theEQE,
							   double      theEQF,
							   double      theEQG,
							   double      theEQH);

void material_update_eq ( mesh_t     *myMesh,
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
							   double      theFreq,
							   double      theEQA,
							   double      theEQB,
							   double      theEQC,
							   double      theEQD,
							   double      theEQE,
							   double      theEQF,
							   double      theEQG,
							   double      theEQH,
							   double      theVSmaxeq
							   );
/* -------------------------------------------------------------------------- */
/*                        Nonlinear Output to Stations                        */
/* -------------------------------------------------------------------------- */
void eqlinear_stations_init ( mesh_t    *myMesh,
                               station_t *myStations,
                               int32_t    myNumberOfStations );

#endif /* Q_EQLINEAR_H */
