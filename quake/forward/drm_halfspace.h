/*
 * drm_halfspace.h
 *
 *  Created on: May 14, 2015
 *      Author: eafit
 */

#ifndef DRM_HALFSPACE_H_
#define DRM_HALFSPACE_H_


typedef enum {
  SV = 0,  P, SINEWAVE
} planewavetype_t;

void    drm_planewaves_init ( int32_t myID, const char *parametersin );
int32_t drm_planewaves_initparameters ( const char *parametersin );
void    PlaneWaves_solver_init( int32_t myID, mesh_t *myMesh, mysolver_t *mySolver);
void    compute_addforce_PlaneWaves ( mesh_t     *myMesh,
                                mysolver_t *mySolver,
                                double      theDeltaT,
                                int         step,
                                fmatrix_t (*theK1)[8], fmatrix_t (*theK2)[8]);

void DRM_ForcesinElement ( mesh_t     *myMesh,
		                   mysolver_t *mySolver,
		                   fmatrix_t (*theK1)[8], fmatrix_t (*theK2)[8],
		                   int *f_nodes, int *e_nodes, int32_t   eindex, double tt, int Nnodes_e, int Nnodes_f );

void   getRicker    ( fvector_t *myDisp, double zp, double t, double Vs );
double Ricker_displ ( double zp, double Ts, double t, double fc, double Vs  );


#endif /* DRM_HALFSPACE_H_ */
