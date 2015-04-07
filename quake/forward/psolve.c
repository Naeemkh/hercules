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

/*
 * Generate an unstructured mesh and solve the linear system thereof
 * derived.
 */
#define _GNU_SOURCE
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <float.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include <stdarg.h>

#include "util.h"
#include "commutil.h"
#include "timers.h"
#include "etree.h"
#include "octor.h"
#include "psolve.h"
#include "cvm.h"
#include "nrutila.h"
#include "quakesource.h"
#include "output.h"
#include "nonlinear.h"
#include "io_planes.h"
#include "io_checkpoint.h"
#include "stiffness.h"
#include "damping.h"
#include "quake_util.h"
#include "buildings.h"
#include "drm.h"
#include "meshformatlab.h"


/* ONLY GLOBAL VARIABLES ALLOWED OUTSIDE OF PARAM. and GLOBAL. IN ALL OF PSOLVE!! */
MPI_Comm comm_solver;
MPI_Comm comm_IO;

#ifndef PROCPERNODE
#define PROCPERNODE     4
#endif

#define PI		3.14159265358979323846264338327

#define GOAHEAD_MSG     100
#define MESH_MSG	101
#define STAT_MSG	102
#define OUT4D_MSG       103
#define DN_MASS_MSG     104
#define AN_MASS_MSG     105
#define DN_FORCE_MSG    106
#define AN_FORCE_MSG    107
#define DN_DISP_MSG     108
#define AN_DISP_MSG     109
#define CVMRECORD_MSG   110

#define BATCH		(1 << 20)
#define LINESIZE	512
#define FILEBUFSIZE     (1 << 25)


#define CONTRIBUTION    901  /**< Harboring processors to owner processors */
#define SHARING		902  /**< Owner processors to harboring processors */

#define DISTRIBUTION    903  /**< Dangling nodes to anchored nodes */
#define ASSIGNMENT      904  /**< Anchored nodes to dangling nodes */


/*---------------Initialization and cleanup routines----------------------*/
static void    read_parameters(int argc, char **argv);
static int32_t parse_parameters(const char *numericalin);
static void    local_finalize(void);

/*---------------- Mesh generation data structures -----------------------*/
#ifdef USECVMDB

#ifndef CVMBUFSIZE
#define CVMBUFSIZE      100
#endif

static void     replicateDB(const char *dbname);
static void     open_cvmdb(void);

#else

static int32_t zsearch(void *base, int32_t count, int32_t recordsize,
		       const point_t *searchpt);
static cvmrecord_t *sliceCVM(const char *cvm_flatfile);

#endif
/** cvmrecord_t: cvm record.  Onlye used if USECVMDB not defined **/
typedef struct cvmrecord_t {
    char key[12];
    float Vp, Vs, density;
} cvmrecord_t;

/**
 * mrecord_t: Complete mesh database record
 *
 */
typedef struct mrecord_t {
    etree_addr_t addr;
    mdata_t mdata;
} mrecord_t;

/* Mesh generation related routines */
static int32_t toexpand(octant_t *leaf, double ticksize, const void *data);
static void    setrec(octant_t *leaf, double ticksize, void *data);
static void    mesh_generate(void);
static int32_t bulkload(etree_t *mep, mrecord_t *partTable, int32_t count);
static void    mesh_output(void);
static void    mesh_correct_properties( etree_t* cvm );

static void    solver_init(void);
static void    solver_printstat(mysolver_t* solver);
static void    solver_delete(void);
static void    solver_run(void);
       void    solver_output_seq(void);
static int     solver_print_schedules(mysolver_t* solver);

static schedule_t* schedule_new(void);
static void    schedule_build(mesh_t *mesh, schedule_t *dnsched,
			      schedule_t *ansched);
static void    schedule_allocMPIctl(schedule_t *sched);
static void    schedule_allocmapping(schedule_t *sched);
static void    schedule_delete(schedule_t *sched);
static void    schedule_prepare(schedule_t *sched, int32_t c_outsize,
				int32_t c_insize, int32_t s_outsize,
				int32_t s_insize);
static void    schedule_senddata(schedule_t *sched, void *valuetable,
				 int32_t itemsperentry, int32_t direction,
				 int32_t msgtag);
static int     schedule_print( schedule_t *sched, char type, FILE* out );
static int     schedule_print_detail(schedule_t* sched, char type, FILE* out);
static int     schedule_print_messenger_list( schedule_t* sched,
					      messenger_t* msg, int count,
					      char type, char cs, FILE* out );

static messenger_t *messenger_new(int32_t procid);
static void    messenger_delete(messenger_t *messenger);
static void    messenger_set(messenger_t *messenger, int32_t outsize,
			     int32_t insize);
static int32_t messenger_countnodes(messenger_t *first);

static void    compute_K(void);
static void constract_Quality_Factor_Table(void);

static void constract_1DVelocity_Model_Table(void);


#ifdef BOUNDARY
static char    compute_setflag(tick_t ldb[3], tick_t ruf[3],
			       tick_t nearendp[3], tick_t farendp[3]);
static void    compute_setboundary(float size, float Vp, float Vs,
				   float rho, int flag, double dashpot[8][3]);
#endif /* BOUNDARY */

static void    compute_setab(double freq, double *aBasePtr, double *bBasePtr);
static void    compute_addforce_s(int32_t timestep);
static void    compute_adjust(void *valuetable, int32_t itemsperentry,
			      int32_t how);

static int     interpolate_station_displacements(int32_t step);


/* ---------- Static global variables ------------------------------------ */

/* These are all of the input parameters - add new ones here */
static struct Param_t {
    char  FourDOutFile[256]; 
    FILE*  FourDOutFp;
    FILE*  theMonitorFileFp;
    char*  theMonitorFileName;
    char  parameters_input_file[256];
    char  cvmdb_input_file[256];
    char  mesh_etree_output_file[256];
    char  planes_input_file[256];
    double  theVsCut;
    double  theFactor;
    double  theFreq;
    double  theFreq_Vel;
    double  theDeltaT;
    double  theDeltaTSquared;
    double  theEndT;
    double  theStartT;
    double  theDomainAzimuth;
    int    monitor_stats_rate;
    double  theSofteningFactor;
    int     theStepMeshingFactor;
    int32_t  theTotalSteps;
    int32_t  theRate;
    damping_type_t  theTypeOfDamping;
    double	theThresholdDamping;
    double	theThresholdVpVs;
    int	   theDampingStatisticsFlag;
    int	   theSchedulePrintErrorCheckFlag;
    int	   theSchedulePrintToStdout;
    int	   theSchedulePrintToFile;
    char*  theSchedulePrintFilename;
    char*	theScheduleStatFilename;
    char*       theMeshStatFilename;
    noyesflag_t  printStationVelocities;
    noyesflag_t  printK;
    noyesflag_t  printStationAccelerations;
    noyesflag_t  includeBuildings;
    noyesflag_t  includeNonlinearAnalysis;
    noyesflag_t  useInfQk;
    int  theTimingBarriersFlag;
    stiffness_type_t   theStiffness;
    int      theStationsPrintRate;
    double*  theStationX;
    double*  theStationY;
    double*  theStationZ;
    int32_t  theNumberOfStations;
    int32_t  myNumberOfStations;
    int      IO_pool_pe_count;
    int32_t  thePlanePrintRate;
    int      theNumberOfPlanes;
    char     theStationsDirOut[256];
    station_t*  myStations;
    int  theCheckPointingRate;
    int    theUseCheckPoint;
    char   theCheckPointingDirOut[256];
    noyesflag_t  storeMeshCoordinatesForMatlab;
    double  the4DOutSize;
    int    theMeshOutFlag;
    char  theCVMFlatFile[128];
    output_parameters_t  theOutputParameters;
    double  theRegionLong;
    double  theRegionLat; 
    double  theRegionDepth;
    double  region_depth_deep_m;
    double  theSurfaceCornersLong[4];
    double  theSurfaceCornersLat[4];
    double  theDomainX;
    double  theDomainY;
    double  theDomainZ;
    noyesflag_t  drmImplement;
    drm_part_t   theDrmPart;

} Param = {
    .FourDOutFp = NULL,
    .theMonitorFileFp = NULL,
    .theMonitorFileName = NULL,
    .theFreq_Vel = 0,
    .monitor_stats_rate = 50,
    .theSchedulePrintErrorCheckFlag = 0,
    .theSchedulePrintToStdout = 0,
    .theSchedulePrintToFile = 0,
    .theSchedulePrintFilename = "schedule_info.txt",
    .theScheduleStatFilename = NULL,
    .theMeshStatFilename = NULL,
    .theSofteningFactor = 0,
    .theStepMeshingFactor = 0,
    .myNumberOfStations = 0,
    .theUseCheckPoint =0,
    .theTimingBarriersFlag = 0,
    .the4DOutSize   = 0,
    .theMeshOutFlag = DO_OUTPUT
};

/* These are all of the remaining global variables - this list should not grow */
static struct Global_t {
    int32_t  myID;
    int32_t  theGroupSize;
    octree_t*  myOctree;
    mesh_t*  myMesh;
    int64_t  theETotal;
    int64_t  theNTotal;
    mysolver_t*  mySolver;
    fvector_t*  myVelocityTable;
    fmatrix_t  theK1[8][8];
    fmatrix_t  theK2[8][8];
    fmatrix_t  theK3[8][8];
    double  theQTABLE[26][6];
    double  the1DVModel[1551][6];
    double  theABase;
    double  theBBase;
    double  theCriticalT;
    double  fastestTimeSteps;
    double  slowestTimeSteps;
    numerics_info_t  theNumericsInformation;
    mpi_info_t  theMPIInformation;
    int32_t  theNodesLoaded;
    int32_t*  theNodesLoadedList;
    vector3D_t*  myForces;
    FILE*  fpsource;
    etree_t*  theCVMEp;
    int32_t  theCVMQueryStage;
    double  theXForMeshOrigin;
    double  theYForMeshOrigin;
    double  theZForMeshOrigin;
    cvmrecord_t*  theCVMRecord;
    int  theCVMRecordSize;
    int  theCVMRecordCount;

} Global = {
    .myID = -1,
    .theGroupSize = -1,
    .myOctree = NULL,
    .myMesh = NULL,
    .theETotal = 0,
    .theNTotal = 0,
    .mySolver = NULL,
    .myVelocityTable = NULL,
    .theCriticalT = 0,
    .fastestTimeSteps = 10000000,
    .slowestTimeSteps = 0,
    .theNodesLoaded = -1,
    .theNodesLoadedList = NULL,
    .myForces = NULL,
    .fpsource = NULL,
    .theCVMRecordSize = sizeof(cvmrecord_t)
};

/* ------------------------------End of declarations------------------------------ */

static inline int
monitor_print( const char* format, ... )
{
    int ret = 0;

    if (format != NULL) {
	va_list args;

	va_start( args, format );

	if (Param.theMonitorFileFp == NULL) {
	    ret = vfprintf( stderr, format, args );
	} else {
	    ret = hu_print_tee_va( Param.theMonitorFileFp, format, args );
	}

	va_end( args );
    }

    return ret;
}



/*-----------Parameter input routines---------------------------------*/

static void read_parameters( int argc, char** argv ){

#define LOCAL_INIT_DOUBLE_MESSAGE_LENGTH 18  /* Must adjust this if adding double params */
#define LOCAL_INIT_INT_MESSAGE_LENGTH 20     /* Must adjust this if adding int params */

    double  double_message[LOCAL_INIT_DOUBLE_MESSAGE_LENGTH];
    int     int_message[LOCAL_INIT_INT_MESSAGE_LENGTH];

    strcpy(Param.parameters_input_file, argv[1]);

    /* PE 0 reads all params from disk */
    if (Global.myID == 0) {
        if (parse_parameters(Param.parameters_input_file) != 0) {
            fprintf(stderr, "Thread 0: Problem reading parameters!\n");
            MPI_Abort(MPI_COMM_WORLD, ERROR);
            exit(1);
        }
    }

    /*Broadcast all double params*/
    double_message[0]  = Param.theVsCut;
    double_message[1]  = Param.theFactor;
    double_message[2]  = Param.theFreq;
    double_message[3]  = Param.theDeltaT;
    double_message[4]  = Param.theDeltaTSquared;
    double_message[5]  = Param.theEndT;
    double_message[6]  = Param.theStartT;
    double_message[7]  = Param.theDomainX;
    double_message[8]  = Param.theDomainY;
    double_message[9]  = Param.theDomainZ;
    double_message[10] = Param.theDomainAzimuth;
    double_message[11] = Param.theThresholdDamping;
    double_message[12] = Param.theThresholdVpVs;
    double_message[13] = Param.theSofteningFactor;
    double_message[14] = Param.theFreq_Vel;
    double_message[15] = Param.theRegionLat;
    double_message[16] = Param.theRegionLong;
    double_message[17] = Param.theRegionDepth;


    MPI_Bcast(double_message, LOCAL_INIT_DOUBLE_MESSAGE_LENGTH, MPI_DOUBLE, 0, comm_solver);

    Param.theVsCut            = double_message[0];
    Param.theFactor           = double_message[1];
    Param.theFreq             = double_message[2];
    Param.theDeltaT           = double_message[3];
    Param.theDeltaTSquared    = double_message[4];
    Param.theEndT             = double_message[5];
    Param.theStartT           = double_message[6];
    Param.theDomainX          = double_message[7];
    Param.theDomainY          = double_message[8];
    Param.theDomainZ          = double_message[9];
    Param.theDomainAzimuth	= double_message[10];
    Param.theThresholdDamping = double_message[11];
    Param.theThresholdVpVs    = double_message[12];
    Param.theSofteningFactor  = double_message[13];
    Param.theFreq_Vel		= double_message[14];
    Param.theRegionLat		= double_message[15];
    Param.theRegionLong		= double_message[16];
    Param.theRegionDepth    = double_message[17];

    /*Broadcast all integer params*/
    int_message[0]  = Param.theTotalSteps;
    int_message[1]  = Param.theRate;
    int_message[2]  = Param.theNumberOfPlanes;
    int_message[3]  = Param.theNumberOfStations;
    int_message[4]  = (int)Param.theTypeOfDamping;
    int_message[5]  = Param.theDampingStatisticsFlag;
    int_message[6]  = Param.theMeshOutFlag;
    int_message[7]  = Param.theCheckPointingRate;
    int_message[8]  = Param.theUseCheckPoint;
    int_message[9]  = (int)Param.includeNonlinearAnalysis;
    int_message[10] = (int)Param.theStiffness;
    int_message[11] = (int)Param.printK;
    int_message[12] = (int)Param.printStationVelocities;
    int_message[13] = (int)Param.printStationAccelerations;
    int_message[14] = Param.theTimingBarriersFlag;
    int_message[15] = (int)Param.includeBuildings;
    int_message[16] = (int)Param.storeMeshCoordinatesForMatlab;
    int_message[17] = (int)Param.drmImplement;
    int_message[18] = (int)Param.useInfQk;
    int_message[19] = Param.theStepMeshingFactor;


    MPI_Bcast(int_message, LOCAL_INIT_INT_MESSAGE_LENGTH, MPI_INT, 0, comm_solver);

    Param.theTotalSteps            = int_message[0];
    Param.theRate                  = int_message[1];
    Param.theNumberOfPlanes              = int_message[2];
    Param.theNumberOfStations            = int_message[3];
    Param.theTypeOfDamping         = int_message[4];
    Param.theDampingStatisticsFlag = int_message[5];
    Param.theMeshOutFlag                 = int_message[6];
    Param.theCheckPointingRate           = int_message[7];
    Param.theUseCheckPoint               = int_message[8];
    Param.includeNonlinearAnalysis       = int_message[9];
    Param.theStiffness                   = int_message[10];
    Param.printK                         = int_message[11];
    Param.printStationVelocities         = int_message[12];
    Param.printStationAccelerations      = int_message[13];
    Param.theTimingBarriersFlag          = int_message[14];
    Param.includeBuildings               = int_message[15];
    Param.storeMeshCoordinatesForMatlab  = int_message[16];
    Param.drmImplement                   = int_message[17];
    Param.useInfQk                       = int_message[18];
    Param.theStepMeshingFactor           = int_message[19];

    /*Broadcast all string params*/
    MPI_Bcast (Param.parameters_input_file,  256, MPI_CHAR, 0, comm_solver);
    MPI_Bcast (Param.theCheckPointingDirOut, 256, MPI_CHAR, 0, comm_solver);
    MPI_Bcast (Param.FourDOutFile,           256, MPI_CHAR, 0, comm_solver);
    MPI_Bcast (Param.cvmdb_input_file,       256, MPI_CHAR, 0, comm_solver);
    MPI_Bcast (Param.mesh_etree_output_file, 256, MPI_CHAR, 0, comm_solver);
    MPI_Bcast (Param.planes_input_file,      256, MPI_CHAR, 0, comm_solver);

    return;
}



static void
local_finalize()
{

    /* Free memory associated with the octree and mesh */
    octor_deletetree(Global.myOctree);
    octor_deletemesh(Global.myMesh);

    /* Free the memory used by the source */
    free(Global.myForces);

    if (Global.myVelocityTable != NULL)
	free(Global.myVelocityTable);

    /* Free memory associated with the solver */
    solver_delete();

    return;
}


/**
 * Parse a text file and return the value of a match string.
 *
 * \return 0 if OK, -1 on error.
 */
int
parsetext (FILE* fp, const char* querystring, const char type, void* result)
{
    const static char delimiters[] = " =\n\t";

    int32_t res = 0, found = 0;

    /* Start from the beginning */
    rewind(fp);


    /* Look for the string until found */
    while (!found) {
	char line[LINESIZE];
	char *name, *value;

	/* Read in one line */
	if (fgets(line, LINESIZE, fp) == NULL)
	    break;

	name = strtok(line, delimiters);
	if ((name != NULL) && (strcmp(name, querystring) == 0)) {
	    found = 1;
	    value = strtok(NULL, delimiters);

	    switch (type) {
	    case 'i':
		res = sscanf(value, "%d", (int *)result);
		break;
	    case 'f':
		res = sscanf(value, "%f", (float *)result);
		break;
	    case 'd':
		res = sscanf(value, "%lf", (double *)result);
		break;
	    case 's':
		res = 1;
		strcpy((char *)result, value);
		break;
	    case 'u':
		res = sscanf(value, "%u", (uint32_t *)result);
		break;
	    default:
		fprintf(stderr, "parsetext: unknown type %c\n", type);
		return -1;
	    }
	}

    }

    return (res == 1) ? 0 : -1;
}




/**
 * This is like \c parsetext with the following differences:
 * - works only for strings;
 * - avoids writting past the end of the string;
 * - return convention is different, it distinguishes between "key not found"
 *   and other type of errors.
 *
 * \return
 *	1 if the key name is found and the value is stored in \c value;
 *	0 if the key name was not found; -1 on error.
 */
static int
read_config_string (FILE* fp, const char* key, char* value_ptr, size_t size)
{
    static const char delimiters[] = " =\n\t";

    int  ret;
    char line[LINESIZE];
    char state[LINESIZE];
    char *name, *value, *state_ptr;


    HU_ASSERT_PTR_ALWAYS( fp );
    HU_ASSERT_PTR_ALWAYS( value_ptr );
    HU_ASSERT_PTR_ALWAYS( key );

    rewind (fp);
    ret   = 0;
    *value_ptr = '\0';

    while (0 == ret && !ferror (fp)) {

	if (fgets (line, LINESIZE, fp) == NULL) {
	    if (!feof (fp)) {
		ret = -1;	/* input error */
	    }
	    break;
	}

	state_ptr = state;
	name      = strtok_r (line, delimiters, &state_ptr);

	if ((name != NULL) && (strcmp (name, key) == 0)) {
	    size_t value_len;

	    value = strtok_r (NULL, delimiters, &state_ptr);

	    if (NULL != value) {
		value_len = strlen (value);

		if (value_len >= size) {
		    ret = -2;	/* return buffer is too short */
		} else {
		    strncpy (value_ptr, value, size);
		    ret = 1;
		}
	    }

	    break;
	}
    }

    return ret;
}


/**
 * Open material database and initialize various static global variables.
 *
 * \return 0 if OK, -1 on error.
 */
static int32_t parse_parameters( const char* numericalin )
{
    FILE* fp;

    /*
     * This used to be a seperate file, now aliased to numericalin.
     * This fakes a seperate file as per old code
     */
    const char* physicsin = numericalin;

    int32_t   samples, rate;
    int       number_output_planes, number_output_stations,
              damping_statistics, use_checkpoint, checkpointing_rate,
              step_meshing;

    double    freq, vscut,
              region_origin_latitude_deg, region_origin_longitude_deg,
              region_azimuth_leftface_deg,
              region_depth_shallow_m, region_length_east_m,
              region_length_north_m, region_depth_deep_m,
              startT, endT, deltaT, softening_factor,
              threshold_damping, threshold_VpVs, freq_vel;
    char      type_of_damping[64],
	      	  checkpoint_path[256],
              include_buildings[64],
              include_nonlinear_analysis[64],
              stiffness_calculation_method[64],
              print_matrix_k[64],
              print_station_velocities[64],
              print_station_accelerations[64],
	      	  mesh_coordinates_for_matlab[64],
    		  implement_drm[64],
    		  use_infinite_qk[64];

    damping_type_t   typeOfDamping     = -1;
    stiffness_type_t stiffness_method  = -1;
    noyesflag_t      have_buildings    = -1;
    noyesflag_t      includeNonlinear  = -1;
    noyesflag_t      printMatrixK      = -1;
    noyesflag_t      printStationVels  = -1;
    noyesflag_t      printStationAccs  = -1;
    noyesflag_t      useInfQk          = -1;

    noyesflag_t      meshCoordinatesForMatlab  = -1;
    noyesflag_t      implementdrm  = -1;


    /* Obtain the specification of the simulation */
    if ((fp = fopen(physicsin, "r")) == NULL)
    {
	fprintf(stderr, "Error opening %s\n", physicsin);
	return -1;
    }

    /* Julio, I know this violates the 80 chars printing limit, but we never
     * print this code and is a heck of a lot easier to read it like this
     * --Ricardo
     */
    if ((parsetext(fp, "region_origin_latitude_deg",  'd', &region_origin_latitude_deg  ) != 0) ||
        (parsetext(fp, "region_origin_longitude_deg", 'd', &region_origin_longitude_deg ) != 0) ||
        (parsetext(fp, "region_depth_shallow_m",      'd', &region_depth_shallow_m      ) != 0) ||
        (parsetext(fp, "region_length_east_m",        'd', &region_length_east_m        ) != 0) ||
        (parsetext(fp, "region_length_north_m",       'd', &region_length_north_m       ) != 0) ||
        (parsetext(fp, "region_depth_deep_m",         'd', &region_depth_deep_m         ) != 0) ||
        (parsetext(fp, "region_azimuth_leftface_deg", 'd', &region_azimuth_leftface_deg ) != 0) ||
        (parsetext(fp, "type_of_damping",             's', &type_of_damping             ) != 0) )
    {
        fprintf(stderr, "Error reading region origin from %s\n", physicsin);
        return -1;
    }

    if ( strcasecmp(type_of_damping, "rayleigh") == 0) {
    	typeOfDamping = RAYLEIGH;
    } else if (strcasecmp(type_of_damping, "mass") == 0) {
    	typeOfDamping = MASS;
    } else if (strcasecmp(type_of_damping, "none") == 0) {
    	typeOfDamping = NONE;
    } else if (strcasecmp(type_of_damping, "bkt") == 0) {
    	typeOfDamping = BKT;
    } else {
    	solver_abort( __FUNCTION_NAME, NULL,
    			"Unknown damping type: %s\n",
    			type_of_damping );
    }


    fclose(fp); /* physics.in */

    if ((fp = fopen(numericalin, "r")) == NULL) {
	fprintf(stderr, "Error opening %s\n", numericalin);
	return -1;
    }

    size_t monitor_name_len = 0;
    hu_config_get_string_def( fp, "monitor_file", &Param.theMonitorFileName,
			      &monitor_name_len, "monitor.txt" );

    /* open the monitor file of the simulation in pe 0 */
    Param.theMonitorFileFp = fopen( Param.theMonitorFileName, "w" );
    if (Param.theMonitorFileFp == NULL) {
	fprintf( stderr,"\n Err opening the monitor file" );
    } else {
	setlinebuf ( Param.theMonitorFileFp );
    }

    xfree_char( &Param.theMonitorFileName );

    /* numerical.in parse texts */
    if ((parsetext(fp, "simulation_wave_max_freq_hz",    'd', &freq                        ) != 0) ||
        (parsetext(fp, "simulation_node_per_wavelength", 'i', &samples                     ) != 0) ||
        (parsetext(fp, "simulation_shear_velocity_min",  'd', &vscut                       ) != 0) ||
        (parsetext(fp, "simulation_start_time_sec",      'd', &startT                      ) != 0) ||
        (parsetext(fp, "simulation_end_time_sec",        'd', &endT                        ) != 0) ||
        (parsetext(fp, "simulation_delta_time_sec",      'd', &deltaT                      ) != 0) ||
        (parsetext(fp, "softening_factor",               'd', &softening_factor            ) != 0) ||
        (parsetext(fp, "use_progressive_meshing",        'i', &step_meshing                ) != 0) ||
        (parsetext(fp, "simulation_output_rate",         'i', &rate                        ) != 0) ||
        (parsetext(fp, "number_output_planes",           'i', &number_output_planes        ) != 0) ||
        (parsetext(fp, "number_output_stations",         'i', &number_output_stations      ) != 0) ||
        (parsetext(fp, "the_threshold_damping",          'd', &threshold_damping           ) != 0) ||
        (parsetext(fp, "the_threshold_Vp_over_Vs",       'd', &threshold_VpVs              ) != 0) ||
        (parsetext(fp, "do_damping_statistics",          'i', &damping_statistics          ) != 0) ||
        (parsetext(fp, "use_checkpoint",                 'i', &use_checkpoint              ) != 0) ||
        (parsetext(fp, "checkpointing_rate",             'i', &checkpointing_rate          ) != 0) ||
        (parsetext(fp, "checkpoint_path",                's', &checkpoint_path             ) != 0) ||
        (parsetext(fp, "4D_output_file",                 's', &Param.FourDOutFile          ) != 0) ||
        (parsetext(fp, "cvmdb_input_file",               's', &Param.cvmdb_input_file      ) != 0) ||
        (parsetext(fp, "mesh_etree_output_file",         's', &Param.mesh_etree_output_file) != 0) ||
        (parsetext(fp, "planes_input_file",              's', &Param.planes_input_file     ) != 0) ||
        (parsetext(fp, "include_nonlinear_analysis",     's', &include_nonlinear_analysis  ) != 0) ||
        (parsetext(fp, "stiffness_calculation_method",   's', &stiffness_calculation_method) != 0) ||
        (parsetext(fp, "print_matrix_k",                 's', &print_matrix_k              ) != 0) ||
        (parsetext(fp, "print_station_velocities",       's', &print_station_velocities    ) != 0) ||
        (parsetext(fp, "print_station_accelerations",    's', &print_station_accelerations ) != 0) ||
        (parsetext(fp, "include_buildings",              's', &include_buildings           ) != 0) ||
        (parsetext(fp, "mesh_coordinates_for_matlab",    's', &mesh_coordinates_for_matlab ) != 0) ||
        (parsetext(fp, "implement_drm",    				 's', &implement_drm               ) != 0) ||
        (parsetext(fp, "simulation_velocity_profile_freq_hz",'d', &freq_vel                ) != 0) ||
        (parsetext(fp, "use_infinite_qk",                's', &use_infinite_qk             ) != 0) )
    {
        fprintf( stderr, "Error parsing simulation parameters from %s\n",
                numericalin );
        return -1;
    }

    hu_config_get_int_opt(fp, "output_mesh", &Param.theMeshOutFlag );
    hu_config_get_int_opt(fp, "enable_timing_barriers",&Param.theTimingBarriersFlag);
    hu_config_get_int_opt(fp, "forces_buffer_size", &theForcesBufferSize );
    hu_config_get_int_opt(fp, "schedule_print_file", &Param.theSchedulePrintToFile );

    hu_config_get_int_opt(fp, "schedule_print_error_check",
			  &Param.theSchedulePrintErrorCheckFlag);
    hu_config_get_int_opt(fp, "schedule_print_stdout",
			  &Param.theSchedulePrintToStdout );

    size_t schedule_stat_len = 0;
    hu_config_get_string_def( fp, "stat_schedule_filename",
			      &Param.theScheduleStatFilename,
			      &schedule_stat_len, "stat-sched.txt" );
    size_t mesh_stat_len = 0;
    hu_config_get_string_def( fp, "stat_mesh_filename", &Param.theMeshStatFilename,
			      &mesh_stat_len, "stat-mesh.txt" );

    fclose( fp );

    /* sanity check */
    if (freq <= 0) {
        fprintf(stderr, "Illegal frequency value %f\n", freq);
        return -1;
    }

    if (freq_vel < 0 || freq_vel > freq) {
        fprintf(stderr, "Illegal frequency value, velocity profile frequency can not be smaller than zero or bigger than max freq %f\n", freq_vel);
        return -1;
    }

    if (samples <= 0) {
        fprintf(stderr, "Illegal samples value %d\n", samples);
        return -1;
    }

    if (vscut <= 0) {
        fprintf(stderr, "Illegal vscut value %f\n", vscut);
        return -1;
    }

    if ((startT < 0) || (endT < 0) || (startT > endT)) {
        fprintf(stderr, "Illegal startT %f or endT %f\n", startT, endT);
        return -1;
    }

    if (deltaT <= 0) {
        fprintf(stderr, "Illegal deltaT %f\n", deltaT);
        return -1;
    }

    if ( (softening_factor <= 1) && (softening_factor != 0) ) {
        fprintf(stderr, "Illegal softening factor %f\n", softening_factor);
        return -1;
    }

    if (step_meshing < 0) {
        fprintf(stderr, "Illegal progressive meshing factor %d\n", step_meshing);
        return -1;
    }

    if (rate <= 0) {
        fprintf(stderr, "Illegal output rate %d\n", rate);
        return -1;
    }

    if (number_output_planes < 0) {
        fprintf(stderr, "Illegal number of output planes %d\n",
                number_output_planes);
        return -1;
    }

    if (number_output_stations < 0) {
        fprintf(stderr, "Illegal number of output stations %d\n",
                number_output_planes);
        return -1;
    }

    if (threshold_damping < 0) {
        fprintf(stderr, "Illegal threshold damping %f\n",
                threshold_damping);
        return -1;
    }

    if (threshold_VpVs < 0) {
        fprintf(stderr, "Illegal threshold Vp over Vs %f\n",
                threshold_VpVs);
        return -1;
    }

    if ( (damping_statistics < 0) || (damping_statistics > 1) ) {
        fprintf(stderr, "Illegal do damping statistics flag %d\n",
                damping_statistics);
        return -1;
    }

    if ( (use_checkpoint < 0) || (use_checkpoint > 1) ) {
        fprintf(stderr, "Illegal use checkpoint flag %d\n",
                use_checkpoint);
        return -1;
    }

    if ( checkpointing_rate < 0 ) {
        fprintf(stderr, "Illegal checkpointing rate %d\n",
                use_checkpoint);
        return -1;
    }

    if ( strcasecmp(include_nonlinear_analysis, "yes") == 0 ) {
        includeNonlinear = YES;
    } else if ( strcasecmp(include_nonlinear_analysis, "no") == 0 ) {
        includeNonlinear = NO;
    } else {
        solver_abort( __FUNCTION_NAME, NULL,
        	"Unknown response for including"
                "nonlinear analysis (yes or no): %s\n",
                include_nonlinear_analysis );
    }

    if ( strcasecmp(stiffness_calculation_method, "effective") == 0 ) {
        stiffness_method = EFFECTIVE;
    } else if ( strcasecmp(stiffness_calculation_method, "conventional") == 0 ) {
        stiffness_method = CONVENTIONAL;
    } else {
        solver_abort( __FUNCTION_NAME, NULL, "Unknown response for stiffness"
                "calculation method (effective or conventional): %s\n",
                stiffness_calculation_method );
    }

    if ( strcasecmp(print_matrix_k, "yes") == 0 ) {
        printMatrixK = YES;
    } else if ( strcasecmp(print_matrix_k, "no") == 0 ) {
        printMatrixK = NO;
    } else {
        solver_abort( __FUNCTION_NAME, NULL,
                "Unknown response for printing K matrix (yes or no): %s\n",
                print_matrix_k );
    }

    if ( strcasecmp(print_station_velocities, "yes") == 0 ) {
        printStationVels = YES;
    } else if ( strcasecmp(print_station_velocities, "no") == 0 ) {
        printStationVels = NO;
    } else {
        solver_abort( __FUNCTION_NAME, NULL,
                "Unknown response for printing"
                "station velocities (yes or no): %s\n",
                print_station_velocities );
    }

    if ( strcasecmp(print_station_accelerations, "yes") == 0 ) {
        printStationAccs = YES;
    } else if ( strcasecmp(print_station_accelerations, "no") == 0 ) {
        printStationAccs = NO;
    } else {
        solver_abort( __FUNCTION_NAME, NULL,
                "Unknown response for printing"
                "station accelerations (yes or no): %s\n",
                print_station_accelerations );
    }

    if ( strcasecmp(mesh_coordinates_for_matlab, "yes") == 0 ) {
    	meshCoordinatesForMatlab = YES;
    } else if ( strcasecmp(mesh_coordinates_for_matlab, "no") == 0 ) {
    	meshCoordinatesForMatlab = NO;
    } else {
    	solver_abort( __FUNCTION_NAME, NULL,
    			"Unknown response for mesh coordinates"
    			"for matlab (yes or no): %s\n",
    			mesh_coordinates_for_matlab );
    }

    if ( strcasecmp(include_buildings, "yes") == 0 ) {
        have_buildings = YES;
    } else if ( strcasecmp(include_buildings, "no") == 0 ) {
        have_buildings = NO;
    } else {
        solver_abort( __FUNCTION_NAME, NULL,
                "Unknown response for including buildings (yes or no): %s\n",
                include_buildings );
    }

    if ( strcasecmp(implement_drm, "yes") == 0 ) {
        implementdrm = YES;
    } else if ( strcasecmp(implement_drm, "no") == 0 ) {
        implementdrm = NO;
    } else {
        solver_abort( __FUNCTION_NAME, NULL,
                "Unknown response for impelement_drm (yes or no): %s\n",
                implement_drm );
    }

    if ( strcasecmp(use_infinite_qk, "yes") == 0 ) {
        useInfQk = YES;
    } else if ( strcasecmp(use_infinite_qk, "no") == 0 ) {
        useInfQk = NO;
    } else {
        solver_abort( __FUNCTION_NAME, NULL,
            "Unknown response using infinite Qk (yes or no): %s\n",
                use_infinite_qk);
    }

    /* Init the static global variables */

    Param.theRegionLat      = region_origin_latitude_deg;
    Param.theRegionLong     = region_origin_longitude_deg ;
    Param.theRegionDepth    = region_depth_shallow_m ;

    Param.theVsCut	      = vscut;
    Param.theFactor	      = freq * samples;
    Param.theFreq         = freq;
    Param.theFreq_Vel	  = freq_vel;
    Param.theDeltaT	      = deltaT;
    Param.theDeltaTSquared  = deltaT * deltaT;
    Param.theStartT	      = startT;
    Param.theEndT           = endT;
    Param.theTotalSteps     = (int)(((endT - startT) / deltaT));

    Param.theDomainX	      = region_length_north_m;
    Param.theDomainY	      = region_length_east_m;
    Param.region_depth_deep_m = region_depth_deep_m;
    Param.theDomainZ	      = region_depth_deep_m - region_depth_shallow_m;
    Param.theDomainAzimuth  = region_azimuth_leftface_deg;
    Param.theTypeOfDamping  = typeOfDamping;
    Param.useInfQk          = useInfQk;

    Param.theRate           = rate;

    Param.theNumberOfPlanes	      = number_output_planes;
    Param.theNumberOfStations	      = number_output_stations;

    Param.theSofteningFactor        = softening_factor;
    Param.theStepMeshingFactor     = step_meshing;
    Param.theThresholdDamping	      = threshold_damping;
    Param.theThresholdVpVs	      = threshold_VpVs;
    Param.theDampingStatisticsFlag  = damping_statistics;

    Param.theCheckPointingRate      = checkpointing_rate;
    Param.theUseCheckPoint	      = use_checkpoint;

    Param.includeNonlinearAnalysis  = includeNonlinear;
    Param.theStiffness              = stiffness_method;

    Param.printK                    = printMatrixK;
    Param.printStationVelocities    = printStationVels;
    Param.printStationAccelerations = printStationAccs;

    Param.includeBuildings          = have_buildings;

    Param.storeMeshCoordinatesForMatlab  = meshCoordinatesForMatlab;

    Param.drmImplement              = implementdrm;

    strcpy( Param.theCheckPointingDirOut, checkpoint_path );

    monitor_print("\n\n---------------- Some Input Data ----------------\n\n");
    monitor_print("Vs cut:                             %f\n", Param.theVsCut);
    monitor_print("Softening factor:                   %f\n", Param.theSofteningFactor);
    monitor_print("Number of stations:                 %d\n", Param.theNumberOfStations);
    monitor_print("Number of planes:                   %d\n", Param.theNumberOfPlanes);
    monitor_print("Stiffness calculation method:       %s\n", stiffness_calculation_method);
    monitor_print("Include buildings:                  %s\n", include_buildings);
    monitor_print("Include nonlinear analysis:         %s\n", include_nonlinear_analysis);
    monitor_print("Printing velocities on stations:    %s\n", print_station_velocities);
    monitor_print("Printing accelerations on stations: %s\n", print_station_accelerations);
    monitor_print("Mesh Coordinates For Matlab:        %s\n", mesh_coordinates_for_matlab);
    monitor_print("cvmdb_input_file:                   %s\n", Param.cvmdb_input_file);
    monitor_print("Implement drm:      	               %s\n", implement_drm);
    monitor_print("\n-------------------------------------------------\n\n");

    fflush(Param.theMonitorFileFp);

    return 0;
}



/*-----------Mesh generation related routines------------------------------*/

static void  open_cvmdb(void){

#ifdef USECVMDB

    MPI_Barrier(comm_solver);
    replicateDB(Param.cvmdb_input_file);

    MPI_Barrier(comm_solver);
    Global.theCVMEp = etree_open(Param.cvmdb_input_file, O_RDONLY, CVMBUFSIZE, 0, 0 );

    if (Global.theCVMEp == NULL) {
	fprintf( stderr, "Thread %d: open_cvmdb: error opening CVM etree %s\n",
		 Global.myID, Param.cvmdb_input_file );
	MPI_Abort(MPI_COMM_WORLD, ERROR );
	exit( 1 );
    }

    dbctl_t  *myctl;
    /* Obtain the material database application control/meta data */
    if ((myctl = cvm_getdbctl(Global.theCVMEp)) == NULL) {
    	fprintf(stderr, "Error reading CVM etree control data\n");
    	MPI_Abort(MPI_COMM_WORLD, ERROR );
    	exit( 1 );
    }

    /* Check the ranges of the mesh and the scope of the CVM etree */
    if ((Param.theRegionLat < myctl->region_origin_latitude_deg) ||
        (Param.theRegionLong < myctl->region_origin_longitude_deg) ||
        (Param.theRegionDepth < myctl->region_depth_shallow_m) ||
        (Param.region_depth_deep_m > myctl->region_depth_deep_m) ||
        (Param.theRegionLat + Param.theDomainX / DIST1LAT
         > myctl->region_origin_latitude_deg
         + myctl->region_length_north_m / DIST1LAT) ||
        (Param.theRegionLong + Param.theDomainY / DIST1LON
         > myctl->region_origin_longitude_deg +
         myctl->region_length_east_m / DIST1LON)) {
        fprintf(stderr, "Mesh area out of the CVM etree\n");
        MPI_Abort(MPI_COMM_WORLD, ERROR );
    	exit( 1 );
    }

    /* Compute the coordinates of the origin of the mesh coordinate
       system in the CVM etree domain coordinate system */
    Global.theXForMeshOrigin = (Param.theRegionLat
				- myctl->region_origin_latitude_deg) * DIST1LAT;
    Global.theYForMeshOrigin = (Param.theRegionLong
				- myctl->region_origin_longitude_deg) * DIST1LON;
    Global.theZForMeshOrigin = Param.theRegionDepth - myctl->region_depth_shallow_m;

    /* Free memory used by myctl */
    cvm_freedbctl(myctl);

    double  double_message_extra[3];

    double_message_extra[0] = Global.theXForMeshOrigin;
    double_message_extra[1] = Global.theYForMeshOrigin;
    double_message_extra[2] = Global.theZForMeshOrigin;

    MPI_Bcast(double_message_extra, 3, MPI_DOUBLE, 0, comm_solver);

    Global.theXForMeshOrigin = double_message_extra[0];
    Global.theYForMeshOrigin = double_message_extra[1];
    Global.theZForMeshOrigin = double_message_extra[2];

#else
    strcpy(Param.theCVMFlatFile, cvmdb_input_file);
#endif

}



#ifdef USECVMDB

/**
 * replicateDB: Copy the material database to local disks.
 *
 */
static void
replicateDB(const char *dbname)
{
    char* destdir;

#ifndef SCEC
    char* srcpath;
    MPI_Comm replica_comm;
#endif /* SCEC */


    /* Change the working directory to $LOCAL */
#ifndef CVM_DESTDIR
    char  curdir[256];
    destdir = getenv( "CVM_DESTDIR" );
    if (destdir == NULL) { /* then use current directory */
	destdir = getcwd( curdir, 256 );
    }
#else
    destdir = CVM_DESTDIR;
#endif

    /* Clean the disks:
     * NOTE: Guys! cleanup the disk in your job script before launching
     * psolve, using rm -rf.
     * E.g., on Lemieux add the following line to your PBS job script,
     * before launching psolve.
     *   prun -m cyclic -n ${RMS_NODES} -N ${RMS_NODES} rm -rf <dirname>/
     * where dirname is the directory you want to wipe out.
     * This will take care of this issue.
     *
     * On BigBen it is even easier, it can be done from the front end
     * since everything is a shared parallel file system.
     *
     * Unfortunately the 'system' function is not supported on all platforms,
     * e.g., Cray's XT3 catamount platform (BigBen) does not have it.
     */
#ifndef SCEC
    if (chdir(destdir) != 0) {
	fprintf(stderr, "Thread %d: replicateDB: cannot chdir to %s\n",
		Global.myID, destdir);
	MPI_Abort(MPI_COMM_WORLD, ERROR);
	exit(1);
    }

    MPI_Barrier(comm_solver);

    /* Replicate the material database among the processors */
    if (Global.myID % PROCPERNODE != 0) {
	MPI_Comm_split(comm_solver, MPI_UNDEFINED, Global.myID, &replica_comm);

    } else {
	int replica_id;
	off_t filesize, remains, batchsize;
	void *filebuf;
	int src_fd = -1, dest_fd;

	MPI_Comm_split(comm_solver, 0, Global.myID, &replica_comm);
	MPI_Comm_rank(replica_comm, &replica_id);

	if (replica_id == 0) {
	    struct stat statbuf;

#ifndef CVM_SRCPATH
	    srcpath = getenv("CVM_SRCPATH");
#else
	    srcpath = CVM_SRCPATH;
#endif

	    if (stat(srcpath, &statbuf) != 0) {
		fprintf(stderr,
			"Thread 0: replicateDB: Cannot get stat of %s\n",
			srcpath);
		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	    }

	    filesize = statbuf.st_size;
	    src_fd = open(srcpath, O_RDONLY);
	    if (src_fd == -1) {
		fprintf(stderr,
			"Thread 0: replicateDB: Cannot open cvm source db\n");
		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	    }
	}


	MPI_Bcast(&filesize, sizeof(off_t), MPI_CHAR, 0, replica_comm);
	theDBSize = filesize;

	if ((filebuf = malloc(FILEBUFSIZE)) == NULL) {
	    fprintf(stderr, "Thread %d: replicateDB: ", Global.myID);
	    fprintf(stderr, "run out of memory while ");
	    fprintf(stderr, "preparing to receive material database\n");
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	/* Everyone opens a replicate db */
	dest_fd = open(dbname, O_CREAT|O_TRUNC|O_WRONLY, S_IRUSR|S_IWUSR);
	if (dest_fd == -1) {
	    fprintf(stderr, "Thread %d: replicateDB: ", Global.myID);
	    fprintf(stderr, "cannot create replica database\n");
	    perror("open");
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	remains = filesize;
	while (remains > 0) {
	    batchsize = (remains > FILEBUFSIZE) ? FILEBUFSIZE : remains;

	    if (replica_id == 0) {
		if (read(src_fd, filebuf, batchsize) !=	 batchsize) {
		    fprintf(stderr, "Thread 0: replicateDB: ");
		    fprintf(stderr, "Cannot read database\n");
		    perror("read");
		    MPI_Abort(MPI_COMM_WORLD, ERROR);
		    exit(1);
		}
	    }

	    MPI_Bcast(filebuf, batchsize, MPI_CHAR, 0, replica_comm);

	    if (write(dest_fd, filebuf, batchsize) != batchsize) {
		fprintf(stderr, "Thread %d: replicateDB: ", Global.myID);
		fprintf(stderr, "Cannot write replica database\n");
		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	    }

	    remains -= batchsize;
	}

	free(filebuf);

	if (close(dest_fd) != 0) {
	    fprintf(stderr, "Thread %d: replicateDB: ", Global.myID);
	    fprintf(stderr, "cannot close replica database\n");
	    perror("close");
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	if (replica_id == 0) {
	    close(src_fd);
	}

	MPI_Comm_free(&replica_comm);

    } /* processors participating in the replication */

#endif /* SCEC */

    return ;
}



/**
 * Assign values (material properties) to a leaf octant specified by
 * octleaf.  In order to refine the mesh, select the minimum Vs of 27
 * sample points: 8 near the corners and 19 midpoints.
 */
static void
setrec( octant_t* leaf, double ticksize, void* data )
{
	double x_m, y_m, z_m;	/* x:south-north, y:east-west, z:depth */
	tick_t halfticks;
	cvmpayload_t g_props;	/* cvm record with ground properties */
	cvmpayload_t g_props_min;	/* cvm record with the min Vs found */

	int i_x, i_y, i_z, n_points = 3;
	double points[3];

	int res = 0;
	edata_t* edata = (edata_t*)data;

	points[0] = 0.01;
	points[1] = 1;
	points[2] = 1.99;

	halfticks = (tick_t)1 << (PIXELLEVEL - leaf->level - 1);
	edata->edgesize = ticksize * halfticks * 2;

	/* Check for buildings and proceed according to the buildings setrec */
	if ( Param.includeBuildings == YES ) {
		if ( bldgs_setrec( leaf, ticksize, edata, Global.theCVMEp,Global.theXForMeshOrigin,Global.theYForMeshOrigin,Global.theZForMeshOrigin ) ) {
			return;
		}
	}

	g_props_min.Vs  = FLT_MAX;
	g_props_min.Vp  = NAN;
	g_props_min.rho = NAN;

	for ( i_x = 0; i_x < n_points; i_x++ ) {

		x_m = (Global.theXForMeshOrigin
				+ (leaf->lx + points[i_x] * halfticks) * ticksize);

		for ( i_y = 0; i_y < n_points; i_y++ ) {

			y_m  = Global.theYForMeshOrigin
					+ (leaf->ly + points[i_y] * halfticks) * ticksize;

			for ( i_z = 0; i_z < n_points; i_z++) {

				z_m = Global.theZForMeshOrigin
						+ (leaf->lz +  points[i_z] * halfticks) * ticksize;

				/* Shift the domain if buildings are considered */
				if ( Param.includeBuildings == YES ) {
					z_m -= get_surface_shift();
				}

				res = cvm_query( Global.theCVMEp, y_m, x_m, z_m, &g_props );

				if (res != 0) {
					continue;
				}

				if ( g_props.Vs < g_props_min.Vs ) {
					/* assign minimum value of vs to produce elements
					 * that are small enough to rightly represent the model */
					g_props_min = g_props;
				}

				if (g_props.Vs <= Param.theVsCut) {
					/* stop early if needed, completely break out of all
					 * the loops, the label is just outside the loop */
					goto outer_loop_label;
				}
			}
		}
	}
	outer_loop_label: /* in order to completely break out from the inner loop */

	edata->Vp  = g_props_min.Vp;
	edata->Vs  = g_props_min.Vs;
	edata->rho = g_props_min.rho;

	if (res != 0 && g_props_min.Vs == DBL_MAX) {
		/* all the queries failed, then center out of bound point. Set Vs
		 * to force split */
		edata->Vs = Param.theFactor * edata->edgesize / 2;
	} else if (edata->Vs <= Param.theVsCut) {	/* adjust Vs and Vp */
		double VpVsRatio = edata->Vp / edata->Vs;

		edata->Vs = Param.theVsCut;
		edata->Vp = Param.theVsCut * VpVsRatio;
	}

	return;
}

#else /* USECVMDB */

static int32_t
zsearch(void *base, int32_t count, int32_t recordsize, const point_t *searchpt)
{
    int32_t start, end, offset, found;

    start = 0;
    end = count - 1;
    offset = (start + end ) / 2;

    found = 0;
    do {
	if (end < start) {
	    /* the two pointer crossed each other */
	    offset = end;
	    found = 1;
	} else {
	    const void *pivot = (char *)base + offset * recordsize;

	    switch (octor_zcompare(searchpt, pivot)) {
	    case (0): /* equal */
		found = 1;
		break;
	    case (1): /* searchpoint larger than the pivot */
		start = offset + 1;
		offset = (start + end) / 2;
		break;
	    case (-1): /* searchpoint smaller than the pivot */
		end = offset - 1;
		offset = (start + end) / 2;
		break;
	    }
	}
    } while (!found);

    return offset;
}


static cvmrecord_t *sliceCVM(const char *cvm_flatfile)
{
    cvmrecord_t *cvmrecord;
    int32_t bufferedbytes, bytecount, recordcount;
    if (Global.myID == Global.theGroupSize - 1) {
	/* the last processor reads data and
	   distribute to other processors*/

	struct timeval starttime, endtime;
	float iotime = 0, memmovetime = 0;
	MPI_Request *isendreqs;
	MPI_Status *isendstats;
	FILE *fp;
	int fd, procid;
	struct stat statbuf;
	void *maxbuf;
	const point_t *intervaltable;
	off_t bytesent;
	int32_t offset;
	const int maxcount =  (1 << 29) / sizeof(cvmrecord_t);
	const int maxbufsize = maxcount * sizeof(cvmrecord_t);

	fp = fopen(cvm_flatfile, "r");
	if (fp == NULL) {
	    fprintf(stderr, "Thread %d: Cannot open flat CVM file\n", Global.myID);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	fd = fileno(fp);
	if (fstat(fd, &statbuf) != 0) {
	    fprintf(stderr, "Thread %d: Cannot get the status of CVM file\n",
		    Global.myID);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	intervaltable = octor_getintervaltable(Global.myOctree);

	/*
	for (procid = 0; procid <= Global.myID; procid++) {
	    fprintf(stderr, "interval[%d] = {%d, %d, %d}\n", procid,
		    intervaltable[procid].x << 1, intervaltable[procid].y << 1,
		    intervaltable[procid].z << 1);
	}
	*/

	bytesent = 0;
	maxbuf = malloc(maxbufsize) ;
	if (maxbuf == NULL) {
	    fprintf(stderr, "Thread %d: Cannot allocate send buffer\n", Global.myID);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	isendreqs = (MPI_Request *)malloc(sizeof(MPI_Request) * Global.theGroupSize);
	isendstats = (MPI_Status *)malloc(sizeof(MPI_Status) * Global.theGroupSize);
	if ((isendreqs == NULL) || (isendstats == NULL)) {
	    fprintf(stderr, "Thread %d: Cannot allocate isend controls\n",
		    Global.myID);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	/* Try to read max number of CVM records as allowed */
	//	gettimeofday(&starttime, NULL);
	recordcount = fread(maxbuf, sizeof(cvmrecord_t),
			    maxbufsize / sizeof(cvmrecord_t), fp);
	//	gettimeofday(&endtime, NULL);

	iotime += (endtime.tv_sec - starttime.tv_sec) * 1000.0
	    + (endtime.tv_usec - starttime.tv_usec) / 1000.0;

	if (recordcount != maxbufsize / sizeof(cvmrecord_t)) {
	    fprintf(stderr, "Thread %d: Cannot read-init buffer\n", Global.myID);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	/* start with proc 0 */
	procid = 0;

	while (procid < Global.myID) { /* repeatedly fill the buffer */
	    point_t searchpoint, *point;
	    int newreads;
	    int isendcount = 0;

	    /* we have recordcount to work with */
	    cvmrecord = (cvmrecord_t *)maxbuf;

	    while (procid < Global.myID) { /* repeatedly send out data */

		searchpoint.x = intervaltable[procid + 1].x << 1;
		searchpoint.y = intervaltable[procid + 1].y << 1;
		searchpoint.z = intervaltable[procid + 1].z << 1;

		offset = zsearch(cvmrecord, recordcount, Global.theCVMRecordSize,
				 &searchpoint);

		point = (point_t *)(cvmrecord + offset);

		if ((point->x != searchpoint.x) ||
		    (point->y != searchpoint.y) ||
		    (point->z != searchpoint.z)) {
		    break;
		} else {
		    bytecount = offset * sizeof(cvmrecord_t);
		    MPI_Isend(cvmrecord, bytecount, MPI_CHAR, procid,
			      CVMRECORD_MSG, comm_solver,
			      &isendreqs[isendcount]);
		    isendcount++;

		    /*
		      fprintf(stderr,
		      "Procid = %d offset = %qd bytecount = %d\n",
		      procid, (int64_t)bytesent, bytecount);
		    */

		    bytesent += bytecount;

		    /* prepare for the next processor */
		    recordcount -= offset;
		    cvmrecord = (cvmrecord_t *)point;
		    procid++;
		}
	    }

	    /* Wait till the data in the buffer has been sent */
	    MPI_Waitall(isendcount, isendreqs, isendstats);

	    /* Move residual data to the beginning of the buffer
	       and try to fill the newly free space */
	    bufferedbytes = sizeof(cvmrecord_t) * recordcount;

	    // gettimeofday(&starttime, NULL);
	    memmove(maxbuf, cvmrecord, bufferedbytes);
	    // gettimeofday(&endtime, NULL);
	    memmovetime += (endtime.tv_sec - starttime.tv_sec) * 1000.0
		+ (endtime.tv_usec - starttime.tv_usec) / 1000.0;

	    // gettimeofday(&starttime, NULL);
	    newreads = fread((char *)maxbuf + bufferedbytes,
			     sizeof(cvmrecord_t), maxcount - recordcount, fp);
	    // gettimeofday(&endtime, NULL);
	    iotime += (endtime.tv_sec - starttime.tv_sec) * 1000.0
		+ (endtime.tv_usec - starttime.tv_usec) / 1000.0;

	    recordcount += newreads;

	    if (newreads == 0)
		break;
	}

	free(maxbuf);
	free(isendreqs);
	free(isendstats);

	/* I am supposed to accomodate the remaining octants */
	bytecount = statbuf.st_size - bytesent;

	cvmrecord = (cvmrecord_t *)malloc(bytecount);
	if (cvmrecord == NULL) {
	    fprintf(stderr, "Thread %d: out of memory\n", Global.myID);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	/* fseek exiting the for loop has file cursor propertly */
	if (fseeko(fp, bytesent, SEEK_SET) != 0) {
	    fprintf(stderr, "Thread %d: fseeko failed\n", Global.myID);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	//	gettimeofday(&starttime, NULL);
	if (fread(cvmrecord, 1, bytecount, fp) != (size_t)bytecount) {
	    fprintf(stderr, "Thread %d: fail to read the last chunk\n",
		    Global.myID);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}
	//	gettimeofday(&endtime, NULL);
	iotime += (endtime.tv_sec - starttime.tv_sec) * 1000.0
	    + (endtime.tv_usec - starttime.tv_usec) / 1000.0;

	/*
	fprintf(stderr, "Procid = %d offset = %qd bytecount = %d\n",
		Global.myID, (int64_t)bytesent, bytecount);
	*/

	fclose(fp);

	fprintf(stdout, "Read %s (%.2fMB) in %.2f seconds (%.2fMB/sec)\n",
		cvm_flatfile, (float)statbuf.st_size / (1 << 20),
		iotime / 1000,
		(float)statbuf.st_size / (1 << 20) / (iotime / 1000));

	fprintf(stdout, "Memmove takes %.2f seconds\n",
		(float)memmovetime / 1000);

    } else {
	/* wait for my turn till PE(n - 1) tells me to go ahead */

	MPI_Status status;

	MPI_Probe(Global.theGroupSize - 1, CVMRECORD_MSG, comm_solver, &status);
	MPI_Get_count(&status, MPI_CHAR, &bytecount);

	cvmrecord = (cvmrecord_t *)malloc(bytecount);
	if (cvmrecord == NULL) {
	    fprintf(stderr, "Thread %d: out of memory\n", Global.myID);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	MPI_Recv(cvmrecord, bytecount, MPI_CHAR, Global.theGroupSize - 1,
		 CVMRECORD_MSG, comm_solver,	 &status);

    }

    /* Every processor should set these parameters correctly */
    Global.theCVMRecordCount = bytecount / sizeof(cvmrecord_t);
    if (Global.theCVMRecordCount * sizeof(cvmrecord_t) != (size_t)bytecount) {
	fprintf(stderr, "Thread %d: received corrupted CVM data\n",
		Global.myID);
	MPI_Abort(MPI_COMM_WORLD, ERROR);
	exit(1);
    }

    return cvmrecord;
}


static cvmrecord_t *sliceCVM_old(const char *cvm_flatfile)
{
    cvmrecord_t *cvmrecord;
    int32_t bufferedbytes, bytecount, recordcount;

    if (Global.myID == Global.theGroupSize - 1) {
	/* the last processor reads data and
	   distribute to other processors*/

	FILE *fp;
	int fd, procid;
	struct stat statbuf;
	void *maxbuf;
	const point_t *intervaltable;
	off_t bytesent;
	int32_t offset;
	const int maxcount =  (1 << 29) / sizeof(cvmrecord_t);
	const int maxbufsize = maxcount * sizeof(cvmrecord_t);

	fp = fopen(cvm_flatfile, "r");
	if (fp == NULL) {
	    fprintf(stderr, "Thread %d: Cannot open flat CVM file\n", Global.myID);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	fd = fileno(fp);
	if (fstat(fd, &statbuf) != 0) {
	    fprintf(stderr, "Thread %d: Cannot get the status of CVM file\n",
		    Global.myID);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	intervaltable = octor_getintervaltable(Global.myOctree);
	/*
	for (procid = 0; procid <= Global.myID; procid++) {
	    fprintf(stderr, "interval[%d] = {%d, %d, %d}\n", procid,
		    intervaltable[procid].x << 1, intervaltable[procid].y << 1,
		    intervaltable[procid].z << 1);
	}
	*/

	bytesent = 0;
	maxbuf = malloc(maxbufsize) ;
	if (maxbuf == NULL) {
	    fprintf(stderr, "Thread %d: Cannot allocate send buffer\n", Global.myID);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	/* Try to read max number of CVM records as allowed */
	recordcount = fread(maxbuf, sizeof(cvmrecord_t),
			    maxbufsize / sizeof(cvmrecord_t), fp);

	if (recordcount != maxbufsize / sizeof(cvmrecord_t)) {
	    fprintf(stderr, "Thread %d: Cannot read-init buffer\n", Global.myID);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	/* start with proc 0 */
	procid = 0;

	while (procid < Global.myID) { /* repeatedly fill the buffer */
	    point_t searchpoint, *point;
	    int newreads;

	    /* we have recordcount to work with */
	    cvmrecord = (cvmrecord_t *)maxbuf;

	    while (procid < Global.myID) { /* repeatedly send out data */
		searchpoint.x = intervaltable[procid + 1].x << 1;
		searchpoint.y = intervaltable[procid + 1].y << 1;
		searchpoint.z = intervaltable[procid + 1].z << 1;

		offset = zsearch(cvmrecord, recordcount, Global.theCVMRecordSize,
				 &searchpoint);

		point = (point_t *)(cvmrecord + offset);

		if ((point->x != searchpoint.x) ||
		    (point->y != searchpoint.y) ||
		    (point->z != searchpoint.z)) {
		    break;
		} else {
		    bytecount = offset * sizeof(cvmrecord_t);
		    MPI_Send(cvmrecord, bytecount, MPI_CHAR, procid,
			     CVMRECORD_MSG, comm_solver);
		    /*
		    fprintf(stderr,
			    "Procid = %d offset = %qd bytecount = %d\n",
			    procid, (int64_t)bytesent, bytecount);
		    */

		    bytesent += bytecount;

		    /* prepare for the next processor */
		    recordcount -= offset;
		    cvmrecord = (cvmrecord_t *)point;
		    procid++;
		}
	    }

	    /* Move residual data to the beginning of the buffer
	       and try to fill the newly free space */
	    bufferedbytes = sizeof(cvmrecord_t) * recordcount;
	    memmove(maxbuf, cvmrecord, bufferedbytes);
	    newreads = fread((char *)maxbuf + bufferedbytes,
			     sizeof(cvmrecord_t), maxcount - recordcount, fp);
	    recordcount += newreads;

	    if (newreads == 0)
		break;
	}

	free(maxbuf);

	/* I am supposed to accomodate the remaining octants */
	bytecount = statbuf.st_size - bytesent;

	cvmrecord = (cvmrecord_t *)malloc(bytecount);
	if (cvmrecord == NULL) {
	    fprintf(stderr, "Thread %d: out of memory for %d bytes\n",
		    Global.myID, bytecount);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	/* fseek exiting the for loop has file cursor propertly */
	if (fseeko(fp, bytesent, SEEK_SET) != 0) {
	    fprintf(stderr, "Thread %d: fseeko failed\n", Global.myID);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	if (fread(cvmrecord, 1, bytecount, fp) != (size_t)bytecount) {
	    fprintf(stderr, "Thread %d: fail to read the last chunk\n",
		    Global.myID);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	/*
	  fprintf(stderr, "Procid = %d offset = %qd bytecount = %d\n",
	  Global.myID, (int64_t)bytesent, bytecount);
	*/

	fclose(fp);

    } else {
	/* wait for my turn till PE(n - 1) tells me to go ahead */

	MPI_Status status;

	MPI_Probe(Global.theGroupSize - 1, CVMRECORD_MSG, comm_solver, &status);
	MPI_Get_count(&status, MPI_CHAR, &bytecount);

	cvmrecord = (cvmrecord_t *)malloc(bytecount);
	if (cvmrecord == NULL) {
	    fprintf(stderr, "Thread %d: out of memory\n", Global.myID);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	MPI_Recv(cvmrecord, bytecount, MPI_CHAR, Global.theGroupSize - 1,
		 CVMRECORD_MSG, comm_solver,	 &status);

    }

    /* Every processor should set these parameters correctly */
    Global.theCVMRecordCount = bytecount / sizeof(cvmrecord_t);
    if (Global.theCVMRecordCount * sizeof(cvmrecord_t) != (size_t)bytecount) {
	fprintf(stderr, "Thread %d: received corrupted CVM data\n",
		Global.myID);
	MPI_Abort(MPI_COMM_WORLD, ERROR);
	exit(1);
    }

    return cvmrecord;
}



/**
 * setrec: Search the CVM record array to obtain the material property of
 *	   a leaf octant.
 *
 */
void setrec(octant_t *leaf, double ticksize, void *data)
{
    cvmrecord_t *agghit;
    edata_t *edata;
    etree_tick_t x, y, z;
    etree_tick_t halfticks;
    point_t searchpoint;

    edata = (edata_t *)data;

    halfticks = (tick_t)1 << (PIXELLEVEL - leaf->level - 1);

    edata->edgesize = ticksize * halfticks * 2;

    searchpoint.x = x = leaf->lx + halfticks;
    searchpoint.y = y = leaf->ly + halfticks;
    searchpoint.z = z = leaf->lz + halfticks;

    if ((x * ticksize >= Param.theDomainX) ||
    		(y * ticksize >= Param.theDomainY) ||
			(z * ticksize >= Param.theDomainZ)) {
    	/* Center point out the bound. Set Vs to force split */
    	edata->Vs = Param.theFactor * edata->edgesize / 2;
    } else {
    	int offset;

    	/* map the coordinate from the octor address space to the
	   etree address space */
    	searchpoint.x = x << 1;
    	searchpoint.y = y << 1;
    	searchpoint.z = z << 1;

    	/* Inbound */
    	offset = zsearch(Global.theCVMRecord, Global.theCVMRecordCount, Global.theCVMRecordSize,
    			&searchpoint);
    	if (offset < 0) {
    		fprintf(stderr, "setrec: fatal error\n");
    		MPI_Abort(MPI_COMM_WORLD, ERROR);
    		exit(1);
    	}

    	agghit = Global.theCVMRecord + offset;
    	edata->Vs = agghit->Vs;
    	edata->Vp = agghit->Vp;
    	edata->rho = agghit->density;

    	/* Adjust the Vs */
    	edata->Vs = (edata->Vs < Param.theVsCut) ? Param.theVsCut : edata->Vs;
    }

    return;
}
#endif	/* USECVMDB */


/**
 * mesh_generate: Generate and partition an unstructured octree mesh.
 *
 */
static void
mesh_generate()
{

    int mstep, step = 1;
    double originalFactor = Param.theFactor;
    double ppwl = Param.theFactor / Param.theFreq;
    double prevtref = 0, prevtbal = 0, prevtpar = 0;
    int64_t tote, mine, maxe;

    if (Global.myID == 0) {
        fprintf(stdout, "Meshing: ");
        if (Param.theStepMeshingFactor == 0) {
            fprintf(stdout, "Conventional\n\n");
        } else {
            fprintf(stdout, "Progressive\n\n");
        }
        fprintf(stdout, "Stage %14s Min %7s Max %5s Total    Time(s)","","","");
        if (Param.theStepMeshingFactor == 0) {
            fprintf(stdout, "\n\n");
        } else {
            fprintf(stdout, "   Step  f(Hz)\n\n");
        }
    }

    /*----  Generate and partition an unstructured octree mesh ----*/
    MPI_Barrier(comm_solver);
    Timer_Start("Octor Newtree");
    if (Global.myID == 0) {
        fprintf(stdout, "New tree %41s","");
    }
    Global.myOctree = octor_newtree( Param.theDomainX, Param.theDomainY, Param.theDomainZ,
            sizeof(edata_t), Global.myID, Global.theGroupSize,
            comm_solver, get_surface_shift());

    /* NOTE:
     * If you want to see the carving process, replace by:
     *     Global.myOctree = octor_newtree(
     *             Param.theDomainX, Param.theDomainY, Param.theDomainZ,
     *             sizeof(edata_t), Global.myID, Global.theGroupSize,
     *             comm_solver, 0);
     */

    if (Global.myOctree == NULL) {
        fprintf(stderr, "Thread %d: mesh_generate: fail to create octree\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }
    MPI_Barrier(comm_solver);
    Timer_Stop("Octor Newtree");
    if (Global.myID == 0) {
        fprintf(stdout, "%9.2f\n\n", Timer_Value("Octor Newtree", 0) );
    }

    /* Essential for DRM implementation */
    if (Param.drmImplement == YES) {
        drm_fix_coordinates(Global.myOctree->ticksize);
    }

#ifdef USECVMDB
    Global.theCVMQueryStage = 0; /* Query CVM database to refine the mesh */
#else
     /* Use flat data record file and distibute the data in memories */
    if (Global.myID == 0) {
	fprintf(stdout, "slicing CVMDB ...");
    }
    Timer_Start("Slice CVM");
    Global.theCVMRecord = sliceCVM(Param.theCVMFlatFile);
    MPI_Barrier(comm_solver);
    Timer_Stop("Slice CVM");
    if (Global.theCVMRecord == NULL) {
	fprintf(stderr, "Thread %d: Error obtaining the CVM records from %s\n",
		Global.myID, Param.theCVMFlatFile);
	MPI_Abort(MPI_COMM_WORLD, ERROR);
	exit(1);
    };
    if (Global.myID == 0) {
	fprintf(stdout, "done : %9.2f seconds\n", Timer_Value("Slice CVM", 0));
    }
#endif

    for ( mstep = Param.theStepMeshingFactor; mstep >= 0; mstep-- ) {

        double myFactor = (double)(1 << mstep); // 2^mstep
        Param.theFactor = originalFactor / myFactor;

        /* Refinement */
        Timer_Start("Octor Refinetree");
        if (Global.myID == 0) {
            fprintf(stdout, "Refining     ");
            fflush(stdout);
        }
        if (octor_refinetree(Global.myOctree, toexpand, setrec) != 0) {
            fprintf(stderr, "Thread %d: mesh_generate: fail to refine octree\n",Global.myID);
            MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
        }
        MPI_Barrier(comm_solver);
        tote = octor_getleavescount(Global.myOctree, GLOBAL);
        mine = octor_getminleavescount(Global.myOctree, GLOBAL);
        maxe = octor_getmaxleavescount(Global.myOctree, GLOBAL);
        if (Global.myID == 0) {
            fprintf(stdout, "%11"INT64_FMT" %11"INT64_FMT" %11"INT64_FMT, mine, maxe, tote);
            fflush(stdout);
        }
        Timer_Stop("Octor Refinetree");
        if (Global.myID == 0) {
            fprintf(stdout, "%11.2f", Timer_Value("Octor Refinetree", 0) - prevtref);
            if (Param.theStepMeshingFactor == 0 ) {
                fprintf(stdout, "\n");
            } else {
                fprintf(stdout, "   %4d %6.2f\n", step, Param.theFactor/ppwl);
            }
            prevtref = Timer_Value("Octor Refinetree", 0);
            fflush(stdout);
        }

        /* Balancing */
        Timer_Start("Octor Balancetree");
        if (Global.myID == 0) {
            fprintf(stdout, "Balancing    ");
            fflush(stdout);
        }
        if (octor_balancetree(Global.myOctree, setrec, Param.theStepMeshingFactor) != 0) {
            fprintf(stderr, "Thread %d: mesh_generate: fail to balance octree\n",Global.myID);
            MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
        }
        MPI_Barrier(comm_solver);
        tote = octor_getleavescount(Global.myOctree, GLOBAL);
        mine = octor_getminleavescount(Global.myOctree, GLOBAL);
        maxe = octor_getmaxleavescount(Global.myOctree, GLOBAL);
        if (Global.myID == 0) {
            fprintf(stdout, "%11"INT64_FMT" %11"INT64_FMT" %11"INT64_FMT, mine, maxe, tote);
            fflush(stdout);
        }
        Timer_Stop("Octor Balancetree");
        if (Global.myID == 0) {
            fprintf(stdout, "%11.2f\n", Timer_Value("Octor Balancetree", 0) - prevtbal);
            prevtbal = Timer_Value("Octor Balancetree", 0);
            fflush(stdout);
        }

        /* Partitioning */
        Timer_Start("Octor Partitiontree");
        if (Global.myID == 0) {
            fprintf(stdout, "Partitioning ");
            fflush(stdout);
        }
        if (octor_partitiontree(Global.myOctree, bldgs_nodesearch) != 0) {
            fprintf(stderr, "Thread %d: mesh_generate: fail to balance load\n",Global.myID);
            MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
        }
        MPI_Barrier(comm_solver);
        tote = octor_getleavescount(Global.myOctree, GLOBAL);
        mine = octor_getminleavescount(Global.myOctree, GLOBAL);
        maxe = octor_getmaxleavescount(Global.myOctree, GLOBAL);
        if (Global.myID == 0) {
            fprintf(stdout, "%11"INT64_FMT" %11"INT64_FMT" %11"INT64_FMT, mine, maxe, tote);
            fflush(stdout);
        }
        Timer_Stop("Octor Partitiontree");
        if (Global.myID == 0) {
            fprintf(stdout, "%11.2f\n\n", Timer_Value("Octor Partitiontree", 0) - prevtpar);
            prevtpar = Timer_Value("Octor Partitiontree", 0);
            fflush(stdout);
        }

        step++;
        fflush(stdout);
        MPI_Barrier(comm_solver);
    }

    /* Buildings Carving */
    if ( Param.includeBuildings == YES ) {

        Timer_Start("Carve Buildings");
        if (Global.myID == 0) {
            fprintf(stdout, "Carving buildings");
            fflush(stdout);
        }

        /* NOTE: If you want to see the carving process, comment next line */
        octor_carvebuildings(Global.myOctree, 1, bldgs_nodesearch);
        MPI_Barrier(comm_solver);
        Timer_Stop("Carve Buildings");
        if (Global.myID == 0) {
            fprintf(stdout, "%9.2f\n", Timer_Value("Carve Buildings", 0) );
            fflush(stdout);
        }

        Timer_Start("Octor Partitiontree");
        if (Global.myID == 0) {
            fprintf(stdout, "Repartitioning");
            fflush(stdout);
        }
        if (octor_partitiontree(Global.myOctree, bldgs_nodesearch) != 0) {
            fprintf(stderr, "Thread %d: mesh_generate: fail to balance load\n",
                    Global.myID);
            MPI_Abort(MPI_COMM_WORLD, ERROR);
            exit(1);
        }
        MPI_Barrier(comm_solver);
        Timer_Stop("Octor Partitiontree");
        if (Global.myID == 0) {
            fprintf(stdout, "%9.2f\n", Timer_Value("Octor Partitiontree", 0));
            fflush(stdout);
        }
    }

    if ( Global.myID == 0 && Param.theStepMeshingFactor !=0 ) {
        fprintf(stdout, "Total refine    %33s %9.2f\n", "", Timer_Value("Octor Refinetree", 0));
        fprintf(stdout, "Total balance   %33s %9.2f\n", "", Timer_Value("Octor Balancetree", 0));
        fprintf(stdout, "Total partition %33s %9.2f\n\n", "", Timer_Value("Octor Partitiontree", 0));
        fflush(stdout);
    }

    Timer_Start("Octor Extractmesh");
    if (Global.myID == 0) {
        fprintf(stdout, "Extracting the mesh %30s","");
        fflush(stdout);
    }
    Global.myMesh = octor_extractmesh(Global.myOctree, bldgs_nodesearch);
    if (Global.myMesh == NULL) {
        fprintf(stderr, "Thread %d: mesh_generate: fail to extract mesh\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }
    MPI_Barrier(comm_solver);
    Timer_Stop("Octor Extractmesh");
    if (Global.myID == 0) {
        fprintf(stdout, "%9.2f\n", Timer_Value("Octor Partitiontree", 0));
    }

    Timer_Start( "Mesh correct properties" );
    /* Re-populates the mesh with actual values from the CVM-etree */
    if (Global.myID == 0) {
        fprintf(stdout,"Correcting mesh properties %23s","");
        fflush(stdout);
    }

    mesh_correct_properties( Global.theCVMEp );

    MPI_Barrier( comm_solver );
    Timer_Stop( "Mesh correct properties" );
    if (Global.myID == 0) {
        fprintf(stdout, "%9.2f\n\n",Timer_Value( "Mesh correct properties", 0 ) );
        fflush(stdout);
    }

#ifdef USECVMDB
    /* Close the material database */
    etree_close(Global.theCVMEp);
#else
    free(Global.theCVMRecord);
#endif /* USECVMDB */
}


/**
 * toexpand: Instruct the Octor library whether a leaf octant needs to
 *	     be expanded or not. Return 1 if true, 0 otherwise.
 *
 */
static int32_t
toexpand(octant_t *leaf, double ticksize, const void *data) {

	if ( data == NULL ) {
		return 1;
	}

	int      res;
	edata_t *edata = (edata_t *)data;

	if ( Param.includeBuildings == YES ) {
		res = bldgs_toexpand( leaf, ticksize, edata, Param.theFactor );
		if ( res != -1 ) {
			return res;
		}
	}

	if ( Param.drmImplement == YES) {
		//if( Param.drmImplement == YES && Param.theDrmPart != PART1 ) {
		res = drm_toexpand( leaf, ticksize, edata );
		if ( res != -1 ) {
			return res;
		}
	}

	return vsrule( edata, Param.theFactor );
}

/**
 * bulkload: Append the data to the end of the mesh database. Return 0 if OK,
 *	     -1 on error.
 *
 */
static int32_t
bulkload(etree_t *mep, mrecord_t *partTable, int32_t count)
{
    int index;

    for (index = 0; index < count; index++) {
	void *payload = &partTable[index].mdata;

	if (etree_append(mep, partTable[index].addr, payload) != 0) {
	    /* Append error */
	    return -1;
	}
    }

    return 0;
}


/** Enumeration of the counters used in the mesh statistics */
enum mesh_stat_count_t {
    ELEMENT_COUNT, NODE_COUNT, DANGLING_COUNT, HARBORED_COUNT, MESH_COUNT_LENGTH
};



static void
mesh_print_stat_imp( int32_t* st, int group_size, FILE* out )
{
    int pid;
    global_id_t total_count[MESH_COUNT_LENGTH] = { 0, 0, 0, 0 };

    fputs( "\n"
	   "# ------------------------------------------------------------\n"
	   "# Mesh statistics:\n"
	   "# ------------------------------------------------------------\n"
	   "# Rank    Elements       Nodes     D-nodes     H-nodes\n", out );

    for (pid = 0; pid < group_size; pid++) {
	fprintf( out, "%06d %11d %11d %11d %11d\n", pid, st[ELEMENT_COUNT],
		 st[NODE_COUNT], st[DANGLING_COUNT], st[HARBORED_COUNT] );

	/* add to total count */
	total_count[ELEMENT_COUNT]  += st[ELEMENT_COUNT];
	total_count[NODE_COUNT]	    += st[NODE_COUNT];
	total_count[DANGLING_COUNT] += st[DANGLING_COUNT];
	total_count[HARBORED_COUNT] += st[HARBORED_COUNT];

	st += MESH_COUNT_LENGTH;	/* move to next row */
    }

    fputs( "\n\n"
	    "# ------------------------------------------------------------\n"
	   "# Total\n"
	   "# ------------------------------------------------------------\n",
	   out );

    fprintf( out, "       %11"INT64_FMT" %11"INT64_FMT" %11"INT64_FMT
	     " %11"INT64_FMT"\n\n",
	     total_count[ELEMENT_COUNT], total_count[NODE_COUNT],
	     total_count[DANGLING_COUNT], total_count[HARBORED_COUNT] );

    /* copy totals to static globals */
    /* TODO this should be computed through different means */
    Global.theETotal  = total_count[ELEMENT_COUNT];
    Global.theNTotal  = total_count[NODE_COUNT];


    /* output aggregate information to the monitor file / stdout */
    monitor_print(
		   "Total elements:                      %11"INT64_FMT"\n"
		   "Total nodes:                         %11"INT64_FMT"\n"
		   "Total dangling nodes:                %11"INT64_FMT"\n\n",
		   total_count[ELEMENT_COUNT], total_count[NODE_COUNT],
		   total_count[DANGLING_COUNT] );
}


static int
mesh_collect_print_stats( local_id_t mesh_stat[MESH_COUNT_LENGTH], int my_id,
			  int group_size, const char* fname )
{
    local_id_t* st = NULL;

    if (0 == my_id) { /* only the root node allocates memory */
	XMALLOC_VAR_N( st, local_id_t, (group_size * MESH_COUNT_LENGTH) );
    }

    MPI_Gather( mesh_stat, MESH_COUNT_LENGTH, MPI_INT,
		st,        MESH_COUNT_LENGTH, MPI_INT, 0, comm_solver );

    if (0 == my_id) { /* the root node prints the stats */
        const size_t bufsize = 1048576;  // 1MB
	FILE* out = hu_fopen( Param.theMeshStatFilename, "w" );

	setvbuf( out, NULL, _IOFBF, bufsize );
	mesh_print_stat_imp( st, group_size, out );
	hu_fclosep( &out );
	xfree_int32_t( &st );
    }

    return 0;
}


static void
mesh_printstat_imp( const mesh_t* mesh, int my_id, int group_size,
		    const char* fname )
{
    local_id_t mesh_stat[MESH_COUNT_LENGTH];

    mesh_stat[ ELEMENT_COUNT  ] = mesh->lenum;
    mesh_stat[ NODE_COUNT     ] = mesh->lnnum;
    mesh_stat[ DANGLING_COUNT ] = mesh->ldnnum;
    mesh_stat[ HARBORED_COUNT ] = mesh->nharbored;

    mesh_collect_print_stats( mesh_stat, my_id, group_size, fname );
}


/**
 * Gather and print mesh statistics to a file with the given name.
 *
 * \param fname Name of the file where the statistics should be stored.
 */
static void
mesh_print_stat( const octree_t* oct, const mesh_t* mesh, int my_id,
		 int group_size, const char* fname )
{
    /* collective function calls */
    int32_t gmin = octor_getminleaflevel( oct, GLOBAL );
    int32_t gmax = octor_getmaxleaflevel( oct, GLOBAL );

    mesh_printstat_imp( mesh, my_id, group_size, fname );

    if (Global.myID == 0) {
	monitor_print( "Maximum leaf level = %d\n", gmax );
	monitor_print( "Minimum leaf level = %d\n", gmin );
    }
}


/**
 * Join elements and nodes, and send to Thread 0 for output.
 */
static void
mesh_output()
{
    int32_t eindex;
    int32_t remains, batch, batchlimit, idx;
    mrecord_t *partTable;

    Timer_Start("Mesh Out");

    batchlimit = BATCH;

    /* Allocate a fixed size buffer space to store the join results */
    partTable = (mrecord_t *)calloc(batchlimit, sizeof(mrecord_t));
    if (partTable == NULL) {
	fprintf(stderr,	 "Thread %d: mesh_output: out of memory\n", Global.myID);
	MPI_Abort(MPI_COMM_WORLD, ERROR);
	exit(1);
    }

    if (Global.myID == 0) {
	etree_t *mep;
	int32_t procid;

	printf("mesh_output ... ");

	mep = etree_open(Param.mesh_etree_output_file, O_CREAT|O_RDWR|O_TRUNC, 0,
			 sizeof(mdata_t),3);
	if (mep == NULL) {
	    fprintf(stderr, "Thread 0: mesh_output: ");
	    fprintf(stderr, "cannot create mesh etree\n");
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	/* Begin an appending operation */
	if (etree_beginappend(mep, 1) != 0) {
	    fprintf(stderr, "Thread 0: mesh_output: \n");
	    fprintf(stderr, "cannot begin an append operation\n");
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	eindex = 0;
	while (eindex < Global.myMesh->lenum) {
	    remains = Global.myMesh->lenum - eindex;
	    batch = (remains < batchlimit) ? remains : batchlimit;

	    for (idx = 0; idx < batch; idx++) {
		mrecord_t *mrecord;
		int32_t whichnode;
		int32_t localnid0;

		mrecord = &partTable[idx];

		/* Fill the address field */
		localnid0 = Global.myMesh->elemTable[eindex].lnid[0];

		mrecord->addr.x = Global.myMesh->nodeTable[localnid0].x;
		mrecord->addr.y = Global.myMesh->nodeTable[localnid0].y;
		mrecord->addr.z = Global.myMesh->nodeTable[localnid0].z;
		mrecord->addr.level = Global.myMesh->elemTable[eindex].level;
		mrecord->addr.type = ETREE_LEAF;

		/* Find the global node ids for the vertices */
		for (whichnode = 0; whichnode < 8; whichnode++) {
		    int32_t localnid;
		    int64_t globalnid;

		    localnid = Global.myMesh->elemTable[eindex].lnid[whichnode];
		    globalnid = Global.myMesh->nodeTable[localnid].gnid;

		    mrecord->mdata.nid[whichnode] = globalnid;
		}

		/* data points to mdata_t type */
		memcpy(&mrecord->mdata.edgesize,
		       Global.myMesh->elemTable[eindex].data,
		       sizeof(edata_t));

		eindex++;
	    } /* for a batch */

	    if (bulkload(mep, partTable, batch) != 0) {
		fprintf(stderr, "Thread 0: mesh_output: ");
		fprintf(stderr, "error bulk-loading data\n");
		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	    }
	} /* for all the elements Thread 0 has */

	/* Receive data from other processors */
	for (procid = 1; procid < Global.theGroupSize; procid++) {
	    MPI_Status status;
	    int32_t rcvbytecount;

	    /* Signal the next processor to go ahead */
	    MPI_Send(NULL, 0, MPI_CHAR, procid, GOAHEAD_MSG, comm_solver);

	    while (1) {
		MPI_Probe(procid, MESH_MSG, comm_solver, &status);
		MPI_Get_count(&status, MPI_CHAR, &rcvbytecount);

		batch = rcvbytecount / sizeof(mrecord_t);

		MPI_Recv(partTable, rcvbytecount, MPI_CHAR, procid,
			 MESH_MSG, comm_solver, &status);

		if (batch == 0) {
		    /* Done */
		    break;
		}

		if (bulkload(mep, partTable, batch) != 0) {
		    fprintf(stderr, "Thread 0: mesh_output: ");
		    fprintf(stderr, "cannot bulkloading data from ");
		    fprintf(stderr, "Thread %d\n", procid);
		    MPI_Abort(MPI_COMM_WORLD, ERROR);
		    exit(1);
		}
	    } /* while there is more data to be received from procid */
	} /* for all the processors */

	/* End the appending operation */
	etree_endappend(mep);

	/* Close the mep to ensure the data is on disk */
	if (etree_close(mep) != 0) {
	    fprintf(stderr, "Thread 0: mesh_output ");
	    fprintf(stderr, "error closing the etree database\n");
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

    } else {
	/* Processors other than 0 needs to send data to 0 */
	int32_t sndbytecount;
	MPI_Status status;

	/* Wait for me turn */
	MPI_Recv(NULL, 0, MPI_CHAR, 0, GOAHEAD_MSG, comm_solver, &status);

	eindex = 0;
	while (eindex < Global.myMesh->lenum) {
	    remains = Global.myMesh->lenum - eindex;
	    batch = (remains < batchlimit) ? remains : batchlimit;

	    for (idx = 0; idx < batch; idx++) {
		mrecord_t *mrecord;
		int32_t whichnode;
		int32_t localnid0;

		mrecord = &partTable[idx];

		/* Fill the address field */
		localnid0 = Global.myMesh->elemTable[eindex].lnid[0];

		mrecord->addr.x = Global.myMesh->nodeTable[localnid0].x;
		mrecord->addr.y = Global.myMesh->nodeTable[localnid0].y;
		mrecord->addr.z = Global.myMesh->nodeTable[localnid0].z;
		mrecord->addr.level = Global.myMesh->elemTable[eindex].level;
		mrecord->addr.type = ETREE_LEAF;

		/* Find the global node ids for the vertices */
		for (whichnode = 0; whichnode < 8; whichnode++) {
		    int32_t localnid;
		    int64_t globalnid;

		    localnid = Global.myMesh->elemTable[eindex].lnid[whichnode];
		    globalnid = Global.myMesh->nodeTable[localnid].gnid;

		    mrecord->mdata.nid[whichnode] = globalnid;
		}

		memcpy(&mrecord->mdata.edgesize,
		       Global.myMesh->elemTable[eindex].data,
		       sizeof(edata_t));

		eindex++;
	    } /* for a batch */


	    /* Send data to proc 0 */
	    sndbytecount = batch * sizeof(mrecord_t);
	    MPI_Send(partTable, sndbytecount, MPI_CHAR, 0, MESH_MSG,
		     comm_solver);
	} /* While there is data left to be sent */

	/* Send an empty message to indicate the end of my transfer */
	MPI_Send(NULL, 0, MPI_CHAR, 0, MESH_MSG, comm_solver);
    }

    /* Free the memory for the partial join results */
    free(partTable);

    Timer_Stop("Mesh Out");

    if (Global.myID == 0) {
	printf("done : %9.2f seconds\n", Timer_Value("Mesh Out", 0) );
    }

    return;
}


/*-----------Computation routines ------------------------------*/

/*
   Macros to facitilate computation

   INTEGRAL_1: Integral_Delfixl_Delfjxl()
   INTEGRAL_2: Integral_Delfixl_Delfjxm()
 */

#define INTEGRAL_1(xki, xkj, xli, xlj, xmi, xmj) \
(4.5 * xki * xkj * (1 + xli * xlj / 3) * (1 + xmi * xmj / 3) / 8)

#define INTEGRAL_2(xki, xlj, xmi, xmj) \
(4.5 * xki * xlj * (1 + xmi * xmj / 3) / 8)

#define DS_TOTAL_INTERVALS	    40
#define DS_TOTAL_PARAMETERS	     6

/**
 * Compute histograms for xi, zeta and associated values to understand
 * what is happenning with those parameters involved the damping and
 * delta T.
 */
static void
damping_statistics (
	double min_xi,
	double max_xi,
	double min_zeta,
	double max_zeta,
	double min_VsVp,
	double max_VsVp,
	double min_VpVsZ,
	double max_VpVsZ,
	double min_Vs,
	double max_Vs
	)
{
    static const int totalintervals  = DS_TOTAL_INTERVALS;
    static const int totalparameters = DS_TOTAL_PARAMETERS;

    int interval, parameter, row, col, matrixelements;

    double  min_VpVs, max_VpVs;

    double  themins[DS_TOTAL_PARAMETERS],
	    themaxs[DS_TOTAL_PARAMETERS],
	    spacing[DS_TOTAL_PARAMETERS];

    int32_t counters[DS_TOTAL_PARAMETERS][DS_TOTAL_INTERVALS],
	    global_counter[DS_TOTAL_PARAMETERS][DS_TOTAL_INTERVALS],
	    global_total[DS_TOTAL_PARAMETERS],
	    eindex;

    /* Initializing clue values and variables */

    min_VpVs = 1 / max_VsVp;
    max_VpVs = 1 / min_VsVp;

    themins[0] = min_zeta;
    themins[1] = min_xi;
    themins[2] = min_VsVp;
    themins[3] = min_VpVs;
    themins[4] = min_VpVsZ;
    themins[5] = min_Vs;

    themaxs[0] = max_zeta;
    themaxs[1] = max_xi;
    themaxs[2] = max_VsVp;
    themaxs[3] = max_VpVs;
    themaxs[4] = max_VpVsZ;
    themaxs[5] = max_Vs;

    for ( row = 0; row < totalparameters; row++ )
    {
	for ( col = 0; col < totalintervals; col++ )
	{
	    counters[row][col] = 0;
	}
	global_total[row] = 0;
    }

    for ( row = 0; row < totalparameters; row++ )
    {
	spacing[row] = ( themaxs[row] - themins[row] ) / totalintervals;
    }

    /* loop over the elements */
    for ( eindex = 0; eindex < Global.myMesh->lenum; eindex++)
    {
	/* loop variables */
	elem_t	*elemp;
	edata_t *edata;
	double	 a, b,
		 omega,
		 elementvalues[6];

	/* capturing the elements */
	elemp = &Global.myMesh->elemTable[eindex];
	edata = (edata_t *)elemp->data;

	/* the parameteres */
	elementvalues[0] = 10 / edata->Vs;
	  /* (edata->Vs < 1500) ? 25 / edata->Vs : 5 / edata->Vs; */
	  /* zeta */
	omega = 3.46410161514 / ( edata->edgesize / edata->Vp );     /* freq in rad */
	a = elementvalues[0] * Global.theABase;			    /* a     */
	b = elementvalues[0] * Global.theBBase;			    /* b     */
	elementvalues[1] = ( a / (2 * omega)) + ( b * omega / 2 );  /* xi    */
	elementvalues[2] = ( edata->Vs / edata->Vp);		    /* Vs/Vp */
	elementvalues[3] = edata->Vp / edata->Vs;		    /* Vp/Vs */
	elementvalues[4] = elementvalues[0] * ( edata->Vp / edata->Vs );
	/* Vp/Vs*zeta  */
	elementvalues[5] = edata->Vs;				    /* Vs    */

	/* loop over the parameters */
	for ( parameter = 0; parameter < totalparameters; parameter++ )
	{
	    /* loop over each interval */
	    for ( interval = 0; interval < totalintervals; interval++)
	    {
		/* loop variables */
		double liminf, limsup;

		/* histogram limits */
		liminf = themins[parameter] + (interval * spacing[parameter]);
		limsup = liminf + spacing[parameter];

		/* for the last interval adjust to the max value */
		if ( interval == totalintervals-1 )
		{
		    limsup = themaxs[parameter];
		}

		/* counting elements within the interval */
		if ( ( elementvalues[parameter] >  liminf ) &&
		     ( elementvalues[parameter] <= limsup ) )
		{
		    counters[parameter][interval]++;
		}
	    } /* ends loop on intervals */
	} /* ends loop on parameters */
    } /*ends loop on elements */

    /* add all counting results from each processor */
    matrixelements = totalparameters * totalintervals;
    MPI_Reduce (&counters[0][0], &global_counter[0][0], matrixelements,
		MPI_INT, MPI_SUM, 0, comm_solver);
    MPI_Bcast (&global_counter[0][0], matrixelements, MPI_INT,0,comm_solver);

    /* sums the total of elements for each histogram */
    if (Global.myID == 0)
    {
	for ( parameter = 0; parameter < totalparameters; parameter++)
	{
	    global_counter[parameter][0]++;
	    for ( interval = 0; interval < totalintervals; interval++)
	    {
		global_total[parameter] = global_total[parameter]
		    + global_counter[parameter][interval];
	    }
	}
    }

    /* MPI Barrier */
    MPI_Barrier( comm_solver );

    /* prints to the terminal the histograms */
    if (Global.myID == 0)
    {
	/* header to identify each column */
	printf("\n\n\tThe histograms of the following parameters: \n\n");
	printf("\t 1. Zeta\n");
	printf("\t 2. Xi\n");
	printf("\t 3. Vs/Vp\n");
	printf("\t 4. Vp/Vs\n");
	printf("\t 5. Vp/Vs*zeta\n");
	printf("\t 6. Vs\n\n");
	printf("\tAre given in the following table\n");
	printf("\t(each column is one of the parameters)\n\n\t");

	/* printing the histograms */
	for ( interval = 0; interval < totalintervals; interval++)
	{
	    for ( parameter = 0; parameter < totalparameters; parameter++)
	    {
		printf("%12d",global_counter[parameter][interval]);
	    }
	    printf("\n\t");
	}

	/* prints the total of elements for each column */
	printf("\n\tTotals:\n\n\t");
	for ( parameter = 0; parameter < totalparameters; parameter++)
	{
	    printf("%12d",global_total[parameter]);
	}

	/* prints the interval witdth */
	printf("\n\n\tAnd the intervals width is:");
	for ( parameter = 0; parameter < totalparameters; parameter++)
	{
	    printf("\n\t %2d. %.6f ",parameter+1,spacing[parameter]);
	}
	printf ("\n\n");
	fflush (stdout);
    }

    return;
} /* end damping_statistics */


/**
 * Determine the limit values associated with the damping and delta_t problem.
 */
static void solver_set_critical_T()
{
    int32_t eindex;			/* element index	 */

    double  min_h_over_Vp = 1e32;	/* the min h/Vp group	 */
    double  min_h_over_Vp_global;
    int32_t min_h_over_Vp_elem_index = -1;

    double  min_dt_factor_X = 1e32;	/* the min delta_t group */
    double  min_dt_factor_Z = 1e32,
	    min_dt_factor_X_global,
	    min_dt_factor_Z_global;
    int32_t min_dt_factor_X_elem_index = -1,
	    min_dt_factor_Z_elem_index = -1;

    double  min_zeta = 1e32;		/* the zeta group	 */
    double  max_zeta = 0,
	    min_zeta_global,
	    max_zeta_global;
    int32_t min_zeta_elem_index = -1,
	    max_zeta_elem_index = -1;

    double  min_xi = 1e32;		/* the xi group		 */
    double  max_xi = 0,
	    min_xi_global,
	    max_xi_global;
    int32_t min_xi_elem_index = -1,
	    max_xi_elem_index = -1;

    double  min_VsVp = 1e32;		/* the Vs/Vp group	 */
    double  min_VsVp_global,
	    max_VsVp = 0,
	    max_VsVp_global;
    int32_t min_VsVp_elem_index = -1,
	    max_VsVp_elem_index = -1;

    double  min_VpVsZ = 1e32;		/* the Vp/Vs group	 */
    double  min_VpVsZ_global,
	    max_VpVsZ = 0,
	    max_VpVsZ_global;
    int32_t min_VpVsZ_elem_index = -1,
	    max_VpVsZ_elem_index = -1;

    double  min_Vs = 1e32;		/* the Vs group		 */
    double  min_Vs_global,
	    max_Vs = 0,
	    max_Vs_global;
    int32_t min_Vs_elem_index = -1,
	    max_Vs_elem_index = -1;

    /* Find the minima and maxima for all needed coefficients */
    /* Start loop over the mesh elements */
    for (eindex = 0; eindex < Global.myMesh->lenum; eindex++)
    {
	/* Loop local variables */

	elem_t	*elemp;	      /* pointer to the mesh database		      */
	edata_t *edata;	      /* pointer to the element data		      */

	double	 ratio;	      /* the h/Vp ratio				      */
	double	 zeta;	      /* the time domain zeta-damping		      */
	double	 xi;	      /* the freq domain xi-damping		      */
	double	 omega;	      /* the element associated freq from w=3.46*Vp/h */
	double	 a, b;	      /* the same constants we use for C = aM + bK    */
	double	 dt_factor_X; /* the factor of 0.577(1-xi)*h/Vp		      */
	double	 dt_factor_Z; /* the factor of 0.577(1-zeta)*h/Vp	      */
	double	 VsVp;	      /* the quotient Vs/Vp			      */
	double	 VpVsZ;	      /* the result of Vp / Vs * zeta		      */
	double	 Vs;	      /* the Vs					      */

	/* Captures the element */

	elemp = &Global.myMesh->elemTable[eindex];
	edata = (edata_t *)elemp->data;

	/* Calculate the clue quantities */

	ratio	    = edata->edgesize / edata->Vp;

	/* Old formula for damping */
	/* zeta	       = (edata->Vs < 1500) ? 25 / edata->Vs : 5 / edata->Vs; */
	/* New formula acording to Graves */

	zeta	    = 10 / edata->Vs;

	omega	    = 3.46410161514 / ratio;
	a	    = zeta * Global.theABase;
	b	    = zeta * Global.theBBase;
	xi	    = ( a / (2 * omega)) + ( b * omega / 2 );
	dt_factor_X = 0.57735026919 * ( 1 - xi ) * ratio;
	dt_factor_Z = 0.57735026919 * ( 1 - zeta ) * ratio;
	VsVp	    = edata->Vs / edata->Vp;
	VpVsZ	    = zeta * ( edata->Vp / edata->Vs );
	Vs	    = edata->Vs;

	/* Updating for extreme values */

	/* ratio */
	if ( ratio < min_h_over_Vp )
	{
	    min_h_over_Vp = ratio;
	    min_h_over_Vp_elem_index = eindex;
	}

	/* dt_factors */
	if ( dt_factor_X < min_dt_factor_X )
	{
	    min_dt_factor_X = dt_factor_X;
	    min_dt_factor_X_elem_index = eindex;
	}
	if ( dt_factor_Z < min_dt_factor_Z )
	{
	    min_dt_factor_Z = dt_factor_Z;
	    min_dt_factor_Z_elem_index = eindex;
	}

	/* min_zeta and max_zeta */
	if ( zeta < min_zeta )
	{
	    min_zeta = zeta;
	    min_zeta_elem_index = eindex;
	}
	if ( zeta > max_zeta )
	{
	    max_zeta = zeta;
	    max_zeta_elem_index = eindex;
	}

	/* min_xi and max_xi */
	if ( xi < min_xi )
	{
	    min_xi = xi;
	    min_xi_elem_index = eindex;
	}
	if ( xi > max_xi )
	{
	    max_xi = xi;
	    max_xi_elem_index = eindex;
	}

	/* min Vs/Vp */
	if ( VsVp < min_VsVp )
	{
	    min_VsVp = VsVp;
	    min_VsVp_elem_index = eindex;
	}
	if ( VsVp > max_VsVp )
	{
	    max_VsVp = VsVp;
	    max_VsVp_elem_index = eindex;
	}

	/* min and max VpVsZ */
	if ( VpVsZ < min_VpVsZ )
	{
	    min_VpVsZ = VpVsZ;
	    min_VpVsZ_elem_index = eindex;
	}
	if ( VpVsZ > max_VpVsZ )
	{
	    max_VpVsZ = VpVsZ;
	    max_VpVsZ_elem_index = eindex;
	}

	/* min Vs */
	if ( Vs < min_Vs )
	{
	    min_Vs = Vs;
	    min_Vs_elem_index = eindex;
	}
	if ( Vs > max_Vs )
	{
	    max_Vs = Vs;
	    max_Vs_elem_index = eindex;
	}

    } /* End of the loop over the mesh elements */

    /* Reducing to global values */
    MPI_Reduce(&min_h_over_Vp,	 &min_h_over_Vp_global,	  1, MPI_DOUBLE, MPI_MIN, 0, comm_solver);
    if ( Param.theDampingStatisticsFlag == 1 )
    {
	MPI_Reduce(&min_dt_factor_X, &min_dt_factor_X_global, 1, MPI_DOUBLE, MPI_MIN, 0, comm_solver);
	MPI_Reduce(&min_dt_factor_Z, &min_dt_factor_Z_global, 1, MPI_DOUBLE, MPI_MIN, 0, comm_solver);
	MPI_Reduce(&min_zeta,	     &min_zeta_global,	      1, MPI_DOUBLE, MPI_MIN, 0, comm_solver);
	MPI_Reduce(&max_zeta,	     &max_zeta_global,	      1, MPI_DOUBLE, MPI_MAX, 0, comm_solver);
	MPI_Reduce(&min_xi,	     &min_xi_global,	      1, MPI_DOUBLE, MPI_MIN, 0, comm_solver);
	MPI_Reduce(&max_xi,	     &max_xi_global,	      1, MPI_DOUBLE, MPI_MAX, 0, comm_solver);
	MPI_Reduce(&min_VsVp,	     &min_VsVp_global,	      1, MPI_DOUBLE, MPI_MIN, 0, comm_solver);
	MPI_Reduce(&max_VsVp,	     &max_VsVp_global,	      1, MPI_DOUBLE, MPI_MAX, 0, comm_solver);
	MPI_Reduce(&min_VpVsZ,	     &min_VpVsZ_global,	      1, MPI_DOUBLE, MPI_MIN, 0, comm_solver);
	MPI_Reduce(&max_VpVsZ,	     &max_VpVsZ_global,	      1, MPI_DOUBLE, MPI_MAX, 0, comm_solver);
	MPI_Reduce(&min_Vs,	     &min_Vs_global,	      1, MPI_DOUBLE, MPI_MIN, 0, comm_solver);
	MPI_Reduce(&max_Vs,	     &max_Vs_global,	      1, MPI_DOUBLE, MPI_MAX, 0, comm_solver);
    }

    /* Inform everyone about the global values */
    MPI_Bcast(&min_h_over_Vp_global,   1, MPI_DOUBLE, 0, comm_solver);
    if ( Param.theDampingStatisticsFlag == 1 )
    {
	MPI_Bcast(&min_dt_factor_X_global, 1, MPI_DOUBLE, 0, comm_solver);
	MPI_Bcast(&min_dt_factor_Z_global, 1, MPI_DOUBLE, 0, comm_solver);
	MPI_Bcast(&min_zeta_global,	   1, MPI_DOUBLE, 0, comm_solver);
	MPI_Bcast(&max_zeta_global,	   1, MPI_DOUBLE, 0, comm_solver);
	MPI_Bcast(&min_xi_global,	   1, MPI_DOUBLE, 0, comm_solver);
	MPI_Bcast(&max_xi_global,	   1, MPI_DOUBLE, 0, comm_solver);
	MPI_Bcast(&min_VsVp_global,	   1, MPI_DOUBLE, 0, comm_solver);
	MPI_Bcast(&max_VsVp_global,	   1, MPI_DOUBLE, 0, comm_solver);
	MPI_Bcast(&min_VpVsZ_global,	   1, MPI_DOUBLE, 0, comm_solver);
	MPI_Bcast(&max_VpVsZ_global,	   1, MPI_DOUBLE, 0, comm_solver);
	MPI_Bcast(&min_Vs_global,	   1, MPI_DOUBLE, 0, comm_solver);
	MPI_Bcast(&max_Vs_global,	   1, MPI_DOUBLE, 0, comm_solver);
    }

    /* go for damping statistics */
    if ( Param.theDampingStatisticsFlag == 1 )
    {
	damping_statistics(min_xi_global,   max_xi_global,   min_zeta_global,  max_zeta_global,
			   min_VsVp_global, max_VsVp_global, min_VpVsZ_global, max_VpVsZ_global,
			   min_Vs_global,   max_Vs_global);
    }

    /* Static global variable for the critical delta t */
    Global.theCriticalT = min_h_over_Vp_global;

    /* Printing of information */
    MPI_Barrier( comm_solver );
    if (Global.myID == 0)
    {
	if ( Param.theDampingStatisticsFlag == 1 )
	{
	    printf("\n\n Critical delta t related information: \n\n");
	    printf("\t 1. The minimum h/Vp	   = %.6f \n", min_h_over_Vp_global);
	    printf("\t 2. The minimum dt X	   = %.6f \n", min_dt_factor_X_global);
	    printf("\t 3. The minimum dt Z	   = %.6f \n", min_dt_factor_Z_global);
	    printf("\t 4. The minimum zeta	   = %.6f \n", min_zeta_global);
	    printf("\t 5. The maximum zeta	   = %.6f \n", max_zeta_global);
	    printf("\t 6. The minimum xi	   = %.6f \n", min_xi_global);
	    printf("\t 7. The maximum xi	   = %.6f \n", max_xi_global);
	    printf("\t 8. The minimum Vs/Vp	   = %.6f \n", min_VsVp_global);
	    printf("\t 9. The maximum Vs/Vp	   = %.6f \n", max_VsVp_global);
	    printf("\t10. The minimum (Vp/Vs)*zeta = %.6f \n", min_VpVsZ_global);
	    printf("\t11. The maximum (Vp/Vs)*zeta = %.6f \n", max_VpVsZ_global);
	    printf("\t12. The minimum Vs	   = %.6f \n", min_Vs_global);
	    printf("\t13. The maximum Vs	   = %.6f \n", max_Vs_global);
	}
	else
	{
	    printf("\n\n Critical delta t related information: \n\n");
	    printf("\t The minimum h/Vp = %.6f \n\n", min_h_over_Vp_global);
	}
    }

#ifdef AUTO_DELTA_T
    /* Set the critical delta T */
    Param.theDeltaT	     = Global.theCriticalT;
    Param.theDeltaTSquared = Param.theDeltaT * Param.theDeltaT;

    /* Set the total steps */
    Param.theTotalSteps    = (int)(((Param.theEndT - Param.theStartT) / Param.theDeltaT));
#endif /* AUTO_DELTA_T */

    /* Printing location and element properties of the maximum values */
    if ( Param.theDampingStatisticsFlag == 1 )
    {
	/* Local variables */

	double	local_extremes[13],
		global_extremes[13];
	int32_t element_indices[13];
	int32_t extreme_index;

	local_extremes[0]  = min_h_over_Vp;
	local_extremes[1]  = min_dt_factor_X;
	local_extremes[2]  = min_dt_factor_Z;
	local_extremes[3]  = min_zeta;
	local_extremes[4]  = max_zeta;
	local_extremes[5]  = min_xi;
	local_extremes[6]  = max_xi;
	local_extremes[7]  = min_VsVp;
	local_extremes[8]  = max_VsVp;
	local_extremes[9]  = min_VpVsZ;
	local_extremes[10] = max_VpVsZ;
	local_extremes[11] = min_Vs;
	local_extremes[12] = max_Vs;

	global_extremes[0]  = min_h_over_Vp_global;
	global_extremes[1]  = min_dt_factor_X_global;
	global_extremes[2]  = min_dt_factor_Z_global;
	global_extremes[3]  = min_zeta_global;
	global_extremes[4]  = max_zeta_global;
	global_extremes[5]  = min_xi_global;
	global_extremes[6]  = max_xi_global;
	global_extremes[7]  = min_VsVp_global;
	global_extremes[8]  = max_VsVp_global;
	global_extremes[9]  = min_VpVsZ_global;
	global_extremes[10] = max_VpVsZ_global;
	global_extremes[11] = min_Vs_global;
	global_extremes[12] = max_Vs_global;

	element_indices[0]  = min_h_over_Vp_elem_index;
	element_indices[1]  = min_dt_factor_X_elem_index;
	element_indices[2]  = min_dt_factor_Z_elem_index;
	element_indices[3]  = min_zeta_elem_index;
	element_indices[4]  = max_zeta_elem_index;
	element_indices[5]  = min_xi_elem_index;
	element_indices[6]  = max_xi_elem_index;
	element_indices[7]  = min_VsVp_elem_index;
	element_indices[8]  = max_VsVp_elem_index;
	element_indices[9]  = min_VpVsZ_elem_index;
	element_indices[10] = max_VpVsZ_elem_index;
	element_indices[11] = min_Vs_elem_index;
	element_indices[12] = max_Vs_elem_index;

	/* Printing section title */
	MPI_Barrier( comm_solver );
	if (Global.myID == 0)
	{
	    printf("\n\t Their corresponding element properties and coordinates are: \n\n");
	}

	/* Loop over the six extreme values */
	MPI_Barrier( comm_solver );
	for ( extreme_index = 0; extreme_index < 13; extreme_index++ )
	{
	    MPI_Barrier( comm_solver );
	    if ( local_extremes[extreme_index] == global_extremes[extreme_index] )
	    {
		tick_t	 ldb[3];
		elem_t	*elemp;
		edata_t *edata;
		int lnid0 = Global.myMesh->elemTable[element_indices[extreme_index]].lnid[0];

		ldb[0] = Global.myMesh->nodeTable[lnid0].x;
		ldb[1] = Global.myMesh->nodeTable[lnid0].y;
		ldb[2] = Global.myMesh->nodeTable[lnid0].z;

		elemp  = &Global.myMesh->elemTable[element_indices[extreme_index]];
		edata  = (edata_t *)elemp->data;

		printf("\t For extreme value No. %d:", extreme_index + 1);
		printf("\n\t\t h = %.6f, Vp = %.6f Vs = %.6f rho = %.6f",
		       edata->edgesize, edata->Vp , edata->Vs, edata->rho);
		printf("\n\t\t by thread %d, element_coord = (%.6f, %.6f, %.6f)\n\n",
		       Global.myID, ldb[0] * Global.myMesh->ticksize, ldb[1] * Global.myMesh->ticksize,
		       ldb[2] * Global.myMesh->ticksize);
	    }
	    MPI_Barrier( comm_solver );
	} /* End of loop over the extreme values */

	if (Global.myID == 0) {
	    fflush (stdout);
	}
    } /* end if damping statistics */

    return;
} /* End solver_set_critical_T */


/**
 * Iterates through all processor to obtain the minimum * edgesize of the
 * mesh.
 */
static void get_minimum_edge()
{
    int32_t eindex;
    double min_h = 1e32, min_h_global;
    int32_t min_h_elem_index;

    /* Find the minimal h/Vp in the domain */
    for (eindex = 0; eindex < Global.myMesh->lenum; eindex++) {
        elem_t *elemp;	 /* pointer to the mesh database */
        edata_t *edata;
        double h;

        elemp = &Global.myMesh->elemTable[eindex];
        edata = (edata_t *)elemp->data;

        /* Update the min_h value. (h = edgesize)  */
        h = edata->edgesize;
        if (h < min_h) {
            min_h = h;
            min_h_elem_index = eindex;
        }
    }

    MPI_Reduce(&min_h, &min_h_global, 1, MPI_DOUBLE,
            MPI_MIN, 0, comm_solver);
    /* Inform everyone of this value */
    MPI_Bcast(&min_h_global, 1, MPI_DOUBLE, 0, comm_solver);

    if (Global.myID == 0) {
        printf("\nThe minimal h	 = %.6f\n\n\n", min_h_global);
    }

    Global.theNumericsInformation.minimumh=min_h_global;

    return;
}


/**
 * Print the stiffness matrix K1, K2 and K3 to the given output stream.
 */
static void
print_K_stdoutput()
{
    int i, j, iloc, jloc;

    fprintf(stdout, "\n\nStiffness Matrix K1 \n\n");
    for ( i = 0 ; i < 8 ; i++){
        for ( iloc = 0; iloc < 3; iloc++){
            for ( j = 0; j < 8; j++) {
                for (jloc = 0; jloc < 3; jloc++) {
                    fprintf(stdout, "%10.2e", Global.theK1[i][j].f[iloc][jloc]);
                }
            }
            fprintf(stdout, "\n");
        }
    }

    fprintf(stdout, "\n\nStiffness Matrix K2 \n\n");
    for ( i = 0 ; i < 8 ; i++){
        for ( iloc = 0; iloc < 3; iloc++){
            for ( j = 0; j < 8; j++) {
                for (jloc = 0; jloc < 3; jloc++) {
                    fprintf(stdout, "%10.2e", Global.theK2[i][j].f[iloc][jloc]);
                }
            }
            fprintf(stdout, "\n");
        }
    }

    fprintf(stdout, "\n\nStiffness Matrix K3 \n\n");
    for ( i = 0 ; i < 8 ; i++){
        for ( iloc = 0; iloc < 3; iloc++){
            for ( j = 0; j < 8; j++) {
                for (jloc = 0; jloc < 3; jloc++) {
                    fprintf(stdout, "%10.2e", Global.theK3[i][j].f[iloc][jloc]);
                }
            }
            fprintf(stdout, "\n");
        }
    }

    fprintf(stdout, "\n\n");
}



/**
 * mu_and_lambda: Calculates mu and lambda according to the element values
 *                of Vp, Vs, and Rho and verifies/applies some rules defined
 *                by Jacobo, Leonardo and Ricardo.  It was originally within
 *                solver_init but was moved out because it needs to be used
 *                in other places as well (nonlinear)
 */
void mu_and_lambda(double *theMu, double *theLambda,
                   edata_t *edata, int32_t eindex)
{

    double mu, lambda;

    mu = edata->rho * edata->Vs * edata->Vs;

    if ( edata->Vp > (edata->Vs * Param.theThresholdVpVs) ) {
        lambda = edata->rho * edata->Vs * edata->Vs * Param.theThresholdVpVs
               * Param.theThresholdVpVs - 2 * mu;
    } else {
        lambda = edata->rho * edata->Vp * edata->Vp - 2 * mu;
    }

    /* Adjust Vs, Vp to fix Poisson ratio problem, formula provided by Jacobo */
    if ( lambda < 0 ) {
        if ( edata->Vs < 500 )
            edata->Vp = 2.45 * edata->Vs;
        else if ( edata->Vs < 1200 )
            edata->Vp = 2 * edata->Vs;
        else
            edata->Vp = 1.87 * edata->Vs;

        lambda = edata->rho * edata->Vp * edata->Vp;
    }

    if ( lambda < 0) {
        fprintf(stderr, "\nThread %d: %d element produces negative lambda = %.6f; Vp = %f; Vs = %f; Rho = %f",
                Global.myID, eindex, lambda, edata->Vp, edata->Vs, edata->rho);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
    }

    /* assign results to returning values */
    *theMu = mu;
    *theLambda = lambda;
}



/**
 * Init matrices and constants, build comm schedule, allocate/init space
 * for the solver.
 */
static void solver_init()
{
    /* local variables */
    int32_t eindex;
    int32_t c_outsize, c_insize, s_outsize, s_insize;

    /* compute the damping parameters a/zeta and b/zeta */
    compute_setab(Param.theFreq, &Global.theABase, &Global.theBBase);

    /* find out the critical delta T of the current simulation */
    /* and goes for the damping statistics if falg is == 1     */
    MPI_Barrier( comm_solver );
    solver_set_critical_T();

    /* find the minimum edge size */

    get_minimum_edge();

    /* Init stiffness matrices and other constants */
    compute_K();

    /* For debugging */
    if ( ( Param.printK == YES ) && ( Global.myID == 0 ) ) {
        print_K_stdoutput();
    }

    compute_setab(Param.theFreq, &Global.theABase, &Global.theBBase);

    /* allocation of memory */
    Global.mySolver = (mysolver_t *)malloc(sizeof(mysolver_t));
    if (Global.mySolver == NULL) {
        fprintf(stderr, "Thread %d: solver_init: out of memory\n", Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    /* Allocate memory */
    Global.mySolver->eTable = (e_t *)calloc(Global.myMesh->lenum, sizeof(e_t));
    Global.mySolver->nTable = (n_t *)calloc(Global.myMesh->nharbored, sizeof(n_t));
    Global.mySolver->tm1    = (fvector_t *)calloc(Global.myMesh->nharbored, sizeof(fvector_t));
    Global.mySolver->tm2    = (fvector_t *)calloc(Global.myMesh->nharbored, sizeof(fvector_t));
    Global.mySolver->force  = (fvector_t *)calloc(Global.myMesh->nharbored, sizeof(fvector_t));
    Global.mySolver->conv_shear_1 = (fvector_t *)calloc(8 * Global.myMesh->lenum, sizeof(fvector_t));
    Global.mySolver->conv_shear_2 = (fvector_t *)calloc(8 * Global.myMesh->lenum, sizeof(fvector_t));
    Global.mySolver->conv_kappa_1 = (fvector_t *)calloc(8 * Global.myMesh->lenum, sizeof(fvector_t));
    Global.mySolver->conv_kappa_2 = (fvector_t *)calloc(8 * Global.myMesh->lenum, sizeof(fvector_t));

    Global.mySolver->dn_sched = schedule_new();
    Global.mySolver->an_sched = schedule_new();

    if ( (Global.mySolver->eTable == NULL) ||
         (Global.mySolver->nTable == NULL) ||
         (Global.mySolver->tm1    == NULL) ||
         (Global.mySolver->tm2    == NULL) ||
         (Global.mySolver->force  == NULL) ||
         (Global.mySolver->conv_shear_1 == NULL) ||
         (Global.mySolver->conv_shear_2 == NULL) ||
         (Global.mySolver->conv_kappa_1 == NULL) ||
         (Global.mySolver->conv_kappa_2 == NULL) ) {

        fprintf(stderr, "Thread %d: solver_init: out of memory\n", Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    if ( Param.printStationAccelerations == YES ) {

        Global.mySolver->tm3 = (fvector_t *)calloc(Global.myMesh->nharbored, sizeof(fvector_t));

        if ( Global.mySolver->tm3 == NULL ) {

            fprintf(stderr, "Thread %d: solver_init: out of memory for accs\n", Global.myID);
            MPI_Abort(MPI_COMM_WORLD, ERROR);
            exit(1);
        }
    }

    /* For each element:
     * Initialize the data structures. tm1, tm2 and force have
     * already been initialized to 0 by calloc(). */
    for (eindex = 0; eindex < Global.myMesh->lenum; eindex++)
    {
        elem_t  *elemp; /* pointer to the mesh database */
        edata_t *edata;
        e_t     *ep;    /* pointer to the element constant table */
        double   mass, M, mu, lambda;
        int j;

#ifdef BOUNDARY
        tick_t  edgeticks;
        tick_t  ldb[3], ruf[3];
        int32_t lnid0;
        char    flag;
        double  dashpot[8][3];
#endif

        double zeta, a, b;

        /* Note the difference between the two tables */
        elemp = &Global.myMesh->elemTable[eindex];
        edata = (edata_t *)elemp->data;
        ep    = &Global.mySolver->eTable[eindex];

        /* Calculate the Lame constants */
        mu_and_lambda(&mu, &lambda, edata, eindex);

        /* coefficients for term (deltaT_squared * Ke * Ut) */
        ep->c1 = Param.theDeltaTSquared * edata->edgesize * mu / 9;
        ep->c2 = Param.theDeltaTSquared * edata->edgesize * lambda / 9;

        /* coefficients for term (b * deltaT * Ke_off * (Ut-1 - Ut)) */
        /* Anelastic attenuation (material damping) */

        /* Old formula for damping */
        /* zeta = (edata->Vs < 1500) ? 25 / edata->Vs : 5 / edata->Vs; */

        /* New formula for damping according to Graves */
        	zeta = 10 / edata->Vs;

        if ( zeta > Param.theThresholdDamping ) {
        	zeta = Param.theThresholdDamping;
        }

        /* the a,b coefficients */
        a = zeta * Global.theABase;
        b = zeta * Global.theBBase;

        /* coefficients for term (b * deltaT * Ke_off * (Ut-1 - Ut)) */
        ep->c3 = b * Param.theDeltaT * edata->edgesize * mu / 9;
        ep->c4 = b * Param.theDeltaT * edata->edgesize * lambda / 9;

#ifdef BOUNDARY

        /* Set the flag for the element */
        lnid0 = elemp->lnid[0];

        ldb[0] = Global.myMesh->nodeTable[lnid0].x;
        ldb[1] = Global.myMesh->nodeTable[lnid0].y;
        ldb[2] = Global.myMesh->nodeTable[lnid0].z;

        edgeticks = (tick_t)1 << (PIXELLEVEL - elemp->level);
        ruf[0] = ldb[0] + edgeticks;
        ruf[1] = ldb[1] + edgeticks;
        ruf[2] = ldb[2] + edgeticks;

        flag = compute_setflag(ldb, ruf, Global.myOctree->nearendp,
                Global.myOctree->farendp);
        if (flag != 13) {
            compute_setboundary(edata->edgesize, edata->Vp, edata->Vs,
                                edata->rho, flag, dashpot);
        }
#endif /* BOUNDARY */

        /* Assign the element mass to its vertices */
        /* mass is the total mass of the element   */
        /* and M is the mass assigned to each node */
        mass = edata->rho * edata->edgesize * edata->edgesize *edata->edgesize;
        M    = mass / 8;

        /* For each node */
        for (j = 0; j < 8; j++)
        {
            int32_t lnid;
            int axis;
            n_t *np;

            lnid = elemp->lnid[j];
            np   = &Global.mySolver->nTable[lnid];

            np->mass_simple += M;

            /* loop for each axis */
            for (axis = 0; axis < 3; axis++ )
            {
                np->mass_minusaM[axis]	-= (Param.theDeltaT * a * M);
                np->mass2_minusaM[axis] -= (Param.theDeltaT * a * M);

#ifdef BOUNDARY
                if (flag != 13)
                {
                    /* boundary impact */
                    np->mass_minusaM[axis]  -= (Param.theDeltaT * dashpot[j][axis]);
                    np->mass2_minusaM[axis] -= (Param.theDeltaT * dashpot[j][axis]);
                }
#endif /* BOUNDARY */

                np->mass_minusaM[axis]	+= M;
                np->mass2_minusaM[axis] += (M * 2);

            } /* end loop for each axis */

        } /* end loop for each node */

    } /* eindex for elements */

    /* Build the communication schedules */
    schedule_build(Global.myMesh, Global.mySolver->dn_sched, Global.mySolver->an_sched);

#ifdef DEBUG
    /* For debug purpose, add gnid into the data field. */
    c_outsize = sizeof(n_t) + sizeof(int64_t);
    c_insize  = sizeof(n_t) + sizeof(int64_t);
    s_outsize = sizeof(n_t) + sizeof(int64_t);
    s_insize  = sizeof(n_t) + sizeof(int64_t);
#else
    c_outsize = sizeof(n_t);
    c_insize  = sizeof(n_t);
    s_outsize = sizeof(n_t);
    s_insize  = sizeof(n_t);
#endif /* DEBUG */

    schedule_prepare(Global.mySolver->dn_sched, c_outsize, c_insize,
		     s_outsize, s_insize);

    schedule_prepare(Global.mySolver->an_sched, c_outsize, c_insize,
		     s_outsize, s_insize);

    /* Send mass information of dangling nodes to their owners */
    schedule_senddata(Global.mySolver->dn_sched, Global.mySolver->nTable,
	      sizeof(n_t) / sizeof(solver_float), CONTRIBUTION, DN_MASS_MSG);

    /* Distribute the mass from dangling nodes to anchored nodes. (local) */
    compute_adjust(Global.mySolver->nTable, sizeof(n_t) / sizeof(solver_float),
		   DISTRIBUTION);

    /* Send mass information of anchored nodes to their owner processors*/
    schedule_senddata(Global.mySolver->an_sched, Global.mySolver->nTable,
	      sizeof(n_t) / sizeof(solver_float), CONTRIBUTION, AN_MASS_MSG);

    return;
}


/**
 * Print the communication schedule of each processor to the given output
 * stream.
 */
static void
solver_printstat_to_stream( mysolver_t* solver, FILE* out )
{
    /**
     * Enumeration for the indices in the counts array
     * the first 4 values are for the dangling nodes:
     * the last 4 values are for the anchored nodes.
     */
    enum solver_msg_idx_t {
	I_DN_C_P,	/**< 0: Dangling contribution processor count */
	I_DN_C_N,	/**< 1: Dangling contribution node count */
	I_DN_S_P,	/**< 2: Dangling shared processor count */
	I_DN_S_N,	/**< 3: Dangling shared node count */
	I_AN_C_P,	/**< 4: Shared contribution processor count */
	I_AN_C_N,	/**< 5: Shared contribution node count */
	I_AN_S_P,	/**< 6: Shared shared processor count */
	I_AN_S_N,	/**< 7: Shared shared node count */
	I_SCHED_LAST	/**< 8: Number of counters */
    };

    int32_t* recv_counts = NULL;
    int32_t  send_counts[I_SCHED_LAST];

    /* number of processors this PE communicates with */
    send_counts[I_DN_C_P] = solver->dn_sched->c_count;
    send_counts[I_DN_S_P] = solver->dn_sched->s_count;
    send_counts[I_AN_C_P] = solver->an_sched->c_count;
    send_counts[I_AN_S_P] = solver->an_sched->s_count;

    /* number of nodes this PE communicates */
    send_counts[I_DN_C_N] = messenger_countnodes( solver->dn_sched->first_c );
    send_counts[I_DN_S_N] = messenger_countnodes( solver->dn_sched->first_s );
    send_counts[I_AN_C_N] = messenger_countnodes( solver->an_sched->first_c );
    send_counts[I_AN_S_N] = messenger_countnodes( solver->an_sched->first_s );


    if (Global.myID == 0) {
	/* allocate a buffer in PE 0 to receive stats from all other PEs */
	XMALLOC_VAR_N( recv_counts, int32_t, I_SCHED_LAST * Global.theGroupSize );
    }

    MPI_Gather( send_counts, I_SCHED_LAST, MPI_INT,
		recv_counts, I_SCHED_LAST, MPI_INT, 0, comm_solver );

    if (Global.myID == 0) {
	int      pe_id;
	int32_t* pe_counts = recv_counts;

	fprintf( out, "# Solver communication schedule summary:\n"
		 "#   PE dc_p     dc_n ds_p     ds_n ac_p     ac_n as_p"
		 "     as_n total_ncnt\n" );

	for (pe_id = 0; pe_id < Global.theGroupSize; pe_id++) {
	    /* print a row for a PE, each row has I_SCHED_LAST (8) items
	     * besides the pe_id and the total node sum at the end */
	    int pe_comm_count = pe_counts[I_DN_C_N] +
				pe_counts[I_DN_S_N] +
				pe_counts[I_AN_C_N] +
				pe_counts[I_AN_S_N];

	    fprintf( out, "%6d %4d %8d %4d %8d %4d %8d %4d %8d %10d\n",
		     pe_id,
		     pe_counts[I_DN_C_P],
		     pe_counts[I_DN_C_N],
		     pe_counts[I_DN_S_P],
		     pe_counts[I_DN_S_N],
		     pe_counts[I_AN_C_P],
		     pe_counts[I_AN_C_N],
		     pe_counts[I_AN_S_P],
		     pe_counts[I_AN_S_N],
		     pe_comm_count );

	    pe_counts += I_SCHED_LAST; /* advance a row of size I_SCHED_LAST */
	}

	fprintf( out, "\n\n" );
	fflush( out );

	free( recv_counts );
    }

    return;
}

/**
 * Print the communication schedule of each processor to the given output
 * stream.
 */
static void
solver_printstat( mysolver_t* solver )
{
    FILE* stat_out = NULL;

    if (Global.myID == 0) {
	stat_out = hu_fopen( Param.theScheduleStatFilename, "w" );
    }

    solver_printstat_to_stream( solver, stat_out );

    if (Global.myID == 0) {
	hu_fclose( stat_out );
	xfree_char( & Param.theScheduleStatFilename );
    }
}


/**
 * solver_delete: Release all the memory associate with Global.mySolver.
 *
 */
static void solver_delete()
{
    if (Global.mySolver == NULL) {
	return;
    }

    free(Global.mySolver->eTable);
    free(Global.mySolver->nTable);

    free(Global.mySolver->tm1);
    free(Global.mySolver->tm2);
    free(Global.mySolver->force);

    free(Global.mySolver->conv_shear_1);
    free(Global.mySolver->conv_shear_2);
    free(Global.mySolver->conv_kappa_1);
    free(Global.mySolver->conv_kappa_2);

    schedule_delete(Global.mySolver->dn_sched);
    schedule_delete(Global.mySolver->an_sched);

    free(Global.mySolver);
}

static int
read_myForces( int32_t timestep )
{
    off_t   whereToRead;
    size_t  to_read, read_count;

    whereToRead = ((off_t)sizeof(int32_t))
		+ Global.theNodesLoaded * sizeof(int32_t)
		+ Global.theNodesLoaded * timestep * sizeof(double) * 3;

    hu_fseeko( Global.fpsource, whereToRead, SEEK_SET );

    to_read    = Global.theNodesLoaded * 3;
    read_count = hu_fread( Global.myForces, sizeof(double), to_read, Global.fpsource );

    return 0;	/* if we got here everything went OK */
}


/**
 * check the max and min value of displacement.
 */
static void
solver_debug_overflow( mysolver_t* solver, mesh_t* mesh, int step )
{
    int nindex;
    double max_disp, min_disp, global_max_disp, global_min_disp;

    max_disp = DBL_MIN;
    min_disp = DBL_MAX;

    /* find the min and max X displacement components */
    for (nindex = 0; nindex < mesh->nharbored; nindex++) {
	fvector_t* tm2Disp;
	n_t* np;

	np      = &solver->nTable[nindex];
	tm2Disp = solver->tm2 + nindex;

	max_disp = (max_disp > tm2Disp->f[0]) ? max_disp : tm2Disp->f[0];
	min_disp = (min_disp < tm2Disp->f[0]) ? min_disp : tm2Disp->f[0];
    }

    /* get global min and max values */
    MPI_Reduce( &min_disp, &global_min_disp, 1, MPI_DOUBLE,
		MPI_MIN, 0, comm_solver );

    MPI_Reduce( &max_disp, &global_max_disp, 1, MPI_DOUBLE,
		MPI_MAX, 0, comm_solver );

    if (Global.myID == 0) {
	printf("Timestep %d: max_dx = %.6f min_dx = %.6f\n",
	       step, global_max_disp, global_min_disp);
    }
}


int
darray_has_nan_nd( const double* v, int dim, const size_t len[] )
{
    HU_ASSERT_PTR_ALWAYS( v );
    HU_ASSERT_PTR_ALWAYS( len );
    HU_ASSERT_ALWAYS( dim > 0 );

    size_t* idx = XMALLOC_N( size_t, dim );
    int     ret = hu_darray_has_nan_nd( v, dim, len, idx );

    if (ret != 0) {
	int i;

	fputs( "WARNING!: Found NAN value at index", stderr );

	for (i = 0; i < dim; i++) {
	    fprintf( stderr, " %zu", idx[i] );
	}

	fputc( '\n', stderr );
    }

    free( idx );
    idx = NULL;

    return ret;
}



/**
 * Check a fvector_t array for NAN values.
 */
static int
fvector_array_has_nan( const fvector_t* f, size_t len, const char* varname )
{
    size_t lengths[2];
    const double* v = (double*)f;

    lengths[0] = len;
    lengths[1] = 3;

    if (darray_has_nan_nd( v, 2, lengths ) != 0) {
        if (NULL == varname) {
	    varname = "(unknown)";
	}
	fputs( "fvector_t varname=", stderr );
	fputs( varname, stderr );
	fputc( '\n', stderr );

	return -1;
    }

    return 0;
}

/**
 * Debug function to check the solver structures for NAN.
 * The checked structures are the displacement arrays (tm1 & tm2) and
 * the forces array.
 */
void
solver_check_nan( mysolver_t* solver, int node_count, int time_step )
{
    HU_ASSERT_PTR_ALWAYS( solver );

    int ret1 = fvector_array_has_nan( solver->tm1, node_count, "tm1" );
    int ret2 = fvector_array_has_nan( solver->tm2, node_count, "tm2" );
    int ret3 = fvector_array_has_nan( solver->force, node_count, "force" );

    if ((ret1 | ret2 | ret3) != 0) {
	hu_solver_abort( __FUNCTION_NAME, NULL,
			 "Found a NAN value at timestep %d", time_step );
    }
}


static void solver_run_init_comm( mysolver_t* solver )
{
    /* The int64_t (global node id) is for debug purpose */
    static const int debug_size = (DO_DEBUG) ? sizeof(int64_t) : 0;

    /* properly setup the communication schedule */
    int c_outsize = sizeof(fvector_t) + debug_size;     /* force */
    int c_insize  = sizeof(fvector_t) + debug_size;     /* displacement */
    int s_outsize = sizeof(fvector_t) + debug_size;     /* displacement */
    int s_insize  = sizeof(fvector_t) + debug_size;     /* force */

    schedule_prepare(solver->dn_sched,c_outsize,c_insize,s_outsize,s_insize);
    schedule_prepare(solver->an_sched,c_outsize,c_insize,s_outsize,s_insize);

    solver_print_schedules( solver );
}


/**
 * \note Globals used:
 * - Global.myID (read)
 * - Param.theDeltaT (read)
 * - startingStep (read)
 */
/* show a progress bar to make the process less anxious! */
static void solver_update_status( int step, const int start_step ){

    static double lastCheckedTime = 0;
    double interval = 0;
    double CurrTime;

    if (Global.myID == 0) {

	CurrTime = Timer_Value( "Total Wall Clock", 0 );

	if (lastCheckedTime==0) {
	    lastCheckedTime = CurrTime;
	    return;
	}

	monitor_print( "*" );

        if (step % Param.monitor_stats_rate == 0) {
	    interval = CurrTime - lastCheckedTime;
	    if (interval > Global.slowestTimeSteps) Global.slowestTimeSteps = interval;
	    if (interval < Global.fastestTimeSteps) Global.fastestTimeSteps = interval;
            monitor_print( "     Sim=% 12.6f     ETA=% 6.1f    WC=% 6.1f\n",
                           step * Param.theDeltaT,
			   ((Param.theTotalSteps - step) / Param.monitor_stats_rate) * interval,
                           CurrTime);
	    lastCheckedTime = CurrTime;
        }

    }
}


static void solver_write_checkpoint( int step, int start_step ){

    if ((Param.theCheckPointingRate != 0) && (step != start_step) &&
        ((step % Param.theCheckPointingRate) == 0)) {

        checkpoint_write( step, Global.myID, Global.myMesh, Param.theCheckPointingDirOut, Global.theGroupSize,
			     Global.mySolver, comm_solver );
    }

}


/**
 * \note Globals used:
 * - Param.theRate (read)
 */
static void solver_output_wavefield( int step )
{
    if (DO_OUTPUT && (step % Param.theRate == 0)) {
        /* output the current timestep */
        do_solver_output();
    }
}


/**
 * \note Globals used:
 * - thePlanePrintRate (read)
 */
static void solver_output_planes( mysolver_t* solver, int my_id, int step )
{
    if (Param.theNumberOfPlanes != 0) {
        if (step % Param.thePlanePrintRate == 0) {
            Timer_Start( "Print Planes" );
            planes_print( my_id, Param.IO_pool_pe_count, Param.theNumberOfPlanes, solver );
            Timer_Stop( "Print Planes" );
        }
    }
}


static void solver_output_stations( int step )
{
    if (Param.theNumberOfStations !=0) {
        if (step % Param.theStationsPrintRate == 0) {
            Timer_Start( "Print Stations" );
            interpolate_station_displacements( step );
            Timer_Stop( "Print Stations" );
        }
    }
}


/**
 * Calculate the nonlinear entities necessary for the next step computation
 * of force correction.
 *
 * \note Globals used:
 * - Param.theNumberOfStations
 * - Param.myNumberOfStations
 * - Param.myStations
 * - Param.theDeltaT
 */
static void solver_nonlinear_state( mysolver_t *solver,
                                    mesh_t     *mesh,
                                    fmatrix_t   k1[8][8],
                                    fmatrix_t   k2[8][8],
                                    int step )
{
    if ( Param.includeNonlinearAnalysis == YES ) {
        Timer_Start( "Compute Non-linear Entities" );
        compute_nonlinear_state ( mesh, solver, Param.theNumberOfStations,
                                  Param.myNumberOfStations, Param.myStations, Param.theDeltaT );
        if ( get_geostatic_total_time() > 0 ) {
            compute_bottom_reactions( mesh, solver, k1, k2, step, Param.theDeltaT );
        }
        Timer_Stop( "Compute Non-linear Entities" );
        if (Param.theNumberOfStations != 0) {
            Timer_Start( "Print Stations" );
            print_nonlinear_stations( mesh, solver, Param.myStations,
                                      Param.myNumberOfStations, Param.theDeltaT,
                                      step, Param.theStationsPrintRate);
            Timer_Stop( "Print Stations" );
        }
    }
}


static void solver_read_source_forces( int step )
{
    Timer_Start( "Read My Forces" );
    if (Global.theNodesLoaded > 0) {
        read_myForces( step );
    }
    Timer_Stop( "Read My Forces" );
}

/**
 * TODO: This method uses global variable Param.theDeltaT
 */
static void solver_load_fixedbase_displacements( mysolver_t* solver, int step )
{
    Timer_Start( "Load Fixedbase Disps" );
    if ( get_fixedbase_flag() == YES ) {
        bldgs_load_fixedbase_disps ( solver, Param.theDeltaT, step);
    }
    Timer_Stop( "Load Fixedbase Disps" );
}


/** Compute the force due to the earthquake source. */
static void solver_compute_force_source( int step )
{
    Timer_Start( "Compute addforces s" );
    compute_addforce_s( step );
    Timer_Stop( "Compute addforces s" );
}


/** Compute the force due to element stiffness matrices. */
static void
solver_compute_force_stiffness( mysolver_t *solver,
                                mesh_t     *mesh,
                                fmatrix_t   k1[8][8],
                                fmatrix_t   k2[8][8] )
{
	Timer_Start( "Compute addforces e" );
	if(Param.theTypeOfDamping != BKT)
	{
		if (Param.theStiffness == EFFECTIVE) {
			compute_addforce_effective( mesh, solver );
		}
		else if (Param.theStiffness == CONVENTIONAL) {
			compute_addforce_conventional( mesh, solver, k1, k2 );
		}
	}
	Timer_Stop( "Compute addforces e" );
}


/** Compute contribution of damping to the force vector */
static void
solver_compute_force_damping( mysolver_t *solver,
                              mesh_t     *mesh,
                              fmatrix_t   k1[8][8],
                              fmatrix_t   k2[8][8] )
{
	Timer_Start( "Damping addforce" );

	if(Param.theTypeOfDamping == RAYLEIGH  || Param.theTypeOfDamping == MASS)
	{
		damping_addforce(Global.myMesh, Global.mySolver, Global.theK1, Global.theK2);
	}
	else if(Param.theTypeOfDamping == BKT)
	{
		calc_conv(Global.myMesh, Global.mySolver, Param.theFreq, Param.theDeltaT, Param.theDeltaTSquared);
		//addforce_conv(myMesh, mySolver, theFreq, theDeltaT, theDeltaTSquared);
		constant_Q_addforce(Global.myMesh, Global.mySolver, Param.theFreq, Param.theDeltaT, Param.theDeltaTSquared);
	}
	else
	{}

	Timer_Stop( "Damping addforce" );
}

/**
 * Compute the nonlinear contribution to the force.
 * \param deltaT2 Delta t^2 (i.e., squared).
 */
static void
solver_compute_force_nonlinear( mysolver_t *solver,
                                mesh_t     *mesh,
                                double      deltaT2 )
{
    if ( Param.includeNonlinearAnalysis == YES ) {
        Timer_Start( "Compute addforces Non-linear" );
        compute_addforce_nl( mesh, solver, deltaT2 );
        Timer_Stop( "Compute addforces Non-linear" );
    }
}

static void
solver_compute_force_gravity( mysolver_t *solver, mesh_t *mesh, int step )
{
    if ( Param.includeNonlinearAnalysis == YES ) {
        Timer_Start( "Compute addforces gravity" );
        if ( get_geostatic_total_time() > 0 ) {
            compute_addforce_gravity( mesh, solver, step, Param.theDeltaT );
        }
        Timer_Stop( "Compute addforces gravity" );
    }
}

/** Send the forces on dangling nodes to their owner processors */
static void solver_send_force_dangling( mysolver_t* solver )
{
//    HU_COND_GLOBAL_BARRIER( Param.theTimingBarriersFlag );

    Timer_Start( "1st schedule send data (contribution)" );
    schedule_senddata( solver->dn_sched, solver->force,
                       sizeof(fvector_t) / sizeof(solver_float), CONTRIBUTION,
                       DN_FORCE_MSG);
    Timer_Stop( "1st schedule send data (contribution)" );
}


static void solver_adjust_forces(  mysolver_t* solver )
{
    Timer_Start( "1st compute adjust (distribution)" );
    /* Distribute the forces to LOCAL anchored nodes */
    compute_adjust( solver->force, sizeof(fvector_t) / sizeof(solver_float),
                    DISTRIBUTION);
    Timer_Stop( "1st compute adjust (distribution)" );
}


static void solver_send_force_anchored( mysolver_t* solver )
{
//    HU_COND_GLOBAL_BARRIER( Param.theTimingBarriersFlag );

    Timer_Start( "2nd schedule send data (contribution)" );
    /* Send the forces on anchored nodes to their owner processors */
    schedule_senddata( solver->an_sched, solver->force,
                       sizeof(fvector_t) / sizeof(solver_float), CONTRIBUTION,
                       AN_FORCE_MSG );
    Timer_Stop( "2nd schedule send data (contribution)" );
}


/** Compute new displacements of my harbored nodes */
static void
solver_compute_displacement( mysolver_t* solver, mesh_t* mesh )
{
    lnid_t nindex;

    Timer_Start( "Compute new displacement" );
    for (nindex = 0; nindex < mesh->nharbored; nindex++) {

        const n_t*       np         = &solver->nTable[nindex];
        fvector_t        nodalForce = solver->force[nindex];
        const fvector_t* tm1Disp    = solver->tm1 + nindex;
        fvector_t*       tm2Disp    = solver->tm2 + nindex;

        /* total nodal forces */
        nodalForce.f[0] += np->mass2_minusaM[0] * tm1Disp->f[0]
                         - np->mass_minusaM[0]  * tm2Disp->f[0];
        nodalForce.f[1] += np->mass2_minusaM[1] * tm1Disp->f[1]
                         - np->mass_minusaM[1]  * tm2Disp->f[1];
        nodalForce.f[2] += np->mass2_minusaM[2] * tm1Disp->f[2]
                         - np->mass_minusaM[2]  * tm2Disp->f[2];

        /* Save tm3 for accelerations */
        if ( Param.printStationAccelerations == YES ) {

            fvector_t* tm3Disp = solver->tm3 + nindex;

            tm3Disp->f[0] = tm2Disp->f[0];
            tm3Disp->f[1] = tm2Disp->f[1];
            tm3Disp->f[2] = tm2Disp->f[2];
        }

        /* overwrite tm2 */
        tm2Disp->f[0] = nodalForce.f[0] / np->mass_simple;
        tm2Disp->f[1] = nodalForce.f[1] / np->mass_simple;
        tm2Disp->f[2] = nodalForce.f[2] / np->mass_simple;

    } /* for (nindex ...): all my harbored nodes */

    /* zero out the force vector for all nodes */
    memset( solver->force, 0, sizeof(fvector_t) * mesh->nharbored );

    Timer_Stop( "Compute new displacement" );
}

static void
solver_geostatic_fix(int step)
{
    if ( Param.includeNonlinearAnalysis == YES ) {
        Timer_Start( "Compute addforces gravity" );
        if ( get_geostatic_total_time() > 0 ) {
            geostatic_displacements_fix( Global.myMesh, Global.mySolver, Param.theDomainZ,
                                         Param.theDeltaT, step );
        }
        Timer_Stop( "Compute addforces gravity" );
    }
}

/** Share the displacement of anchored nodes with other processors. */
static void solver_send_displacement_anchored( mysolver_t* solver )
{
//    HU_COND_GLOBAL_BARRIER( Param.theTimingBarriersFlag );

    Timer_Start("3rd schedule send data (sharing)");
    schedule_senddata(Global.mySolver->an_sched, Global.mySolver->tm2,
                      sizeof(fvector_t) / sizeof(solver_float), SHARING,
                      AN_DISP_MSG);
    Timer_Stop("3rd schedule send data (sharing)");

}


/** Adjust the displacements of my LOCAL dangling nodes. */
static void solver_adjust_displacement(  mysolver_t* solver )
{
    Timer_Start( "2nd compute adjust (assignment)" );
    compute_adjust( Global.mySolver->tm2, sizeof(fvector_t) / sizeof(solver_float),
                    ASSIGNMENT );
    Timer_Stop( "2nd compute adjust (assignment)" );
}


/** Share the displacement of dangling nodes with other processors. */
static void solver_send_displacement_dangling( mysolver_t* solver )
{
//    HU_COND_GLOBAL_BARRIER( Param.theTimingBarriersFlag );

    Timer_Start( "4th schadule send data (sharing)" );
    schedule_senddata( Global.mySolver->dn_sched, Global.mySolver->tm2,
                       sizeof(fvector_t) / sizeof(solver_float), SHARING,
                       DN_DISP_MSG );
    Timer_Stop( "4th schadule send data (sharing)" );
}


/**
 * Hook for instrumentation or other functionality in the solver_run loop.
 */
static void
solver_loop_hook_bottom( mysolver_t* solver, mesh_t* mesh, int step )
{
    if (0) {
        solver_debug_overflow( solver, mesh, step );
    }
}


static void solver_output_wavefield_close( void )
{
    if (DO_OUTPUT && (Param.FourDOutFp != NULL)) {
        fclose( Param.FourDOutFp );     /* close the output file */
    }
}


static void solver_run_collect_timers( void )
{
    /* Get min and max of individual cumulative times.
     * this is NOT best and worst single step times
     */
    if( Timer_Exists("Print Planes") ) {
        Timer_Reduce("Print Planes",   MAX | MIN, comm_solver);
    }

    if( Timer_Exists("Print Stations") ) {
        Timer_Reduce("Print Stations", MAX | MIN, comm_solver);
    }

    if ( Timer_Exists("Compute Non-linear Entities") ) {
        Timer_Reduce("Compute Non-linear Entities", MAX | MIN | AVERAGE, comm_solver);
    }

    if ( Timer_Exists("Compute addforces Non-linear") ) {
        Timer_Reduce("Compute addforces Non-linear", MAX | MIN | AVERAGE, comm_solver);
        Timer_Reduce("Compute addforces gravity",    MAX | MIN | AVERAGE, comm_solver);
    }

    if ( Timer_Exists("Solver drm output") ) {
        Timer_Reduce("Solver drm output", MAX | MIN | AVERAGE, comm_solver);
    }

    if ( Timer_Exists("Solver drm read displacements") ) {
        Timer_Reduce("Solver drm read displacements", MAX | MIN | AVERAGE, comm_solver);
    }

    if ( Timer_Exists("Solver drm force compute") ) {
        Timer_Reduce("Solver drm force compute", MAX | MIN | AVERAGE, comm_solver);
    }

    Timer_Reduce("Read My Forces",                        MAX | MIN | AVERAGE, comm_solver);
    Timer_Reduce("Compute addforces s",                   MAX | MIN | AVERAGE, comm_solver);
    Timer_Reduce("Compute addforces e",                   MAX | MIN | AVERAGE, comm_solver);
    Timer_Reduce("Damping addforce",                      MAX | MIN | AVERAGE, comm_solver);
    Timer_Reduce("1st schedule send data (contribution)", MAX | MIN | AVERAGE, comm_solver);
    Timer_Reduce("1st compute adjust (distribution)",     MAX | MIN | AVERAGE, comm_solver);
    Timer_Reduce("2nd schedule send data (contribution)", MAX | MIN | AVERAGE, comm_solver);
    Timer_Reduce("Compute new displacement",              MAX | MIN | AVERAGE, comm_solver);
    Timer_Reduce("3rd schedule send data (sharing)",      MAX | MIN | AVERAGE, comm_solver);
    Timer_Reduce("2nd compute adjust (assignment)",       MAX | MIN | AVERAGE, comm_solver);
    Timer_Reduce("4th schadule send data (sharing)",      MAX | MIN | AVERAGE, comm_solver);

    Timer_Reduce("Solver I/O",      MAX | MIN | AVERAGE, comm_solver);
    Timer_Reduce("Compute Physics", MAX | MIN | AVERAGE, comm_solver);
    Timer_Reduce("Communication",   MAX | MIN | AVERAGE, comm_solver);
}


/**
 * March forward in time and output the result whenever necessary.
 */
static void solver_run()
{
    int32_t step, startingStep;

    solver_run_init_comm( Global.mySolver );

    /* sets new starting step if loading checkpoint */
    if (Param.theUseCheckPoint == 1) {
        startingStep = checkpoint_read(Global.myID, Global.myMesh, Param.theCheckPointingDirOut,
				       Global.theGroupSize, Global.mySolver,comm_solver);
    } else {
        startingStep = 0;
    }

    if (Global.myID == 0) {
        /* print header for monitor file */
        monitor_print( "solver_run() start\nStarting time step = %d\n\n",
                       startingStep );
        monitor_print( "Sim = Simulation time (s), Sol = Solver time (s), WC = Wall Clock Time (s)\n");
    }

    MPI_Barrier( comm_solver );

    /* march forward in time */
    for (step = startingStep; step < Param.theTotalSteps; step++) {

        fvector_t* tmpvector;

        /* prepare for a new iteration
         * swap displacement vectors for t(n) and t(n-1) */
        tmpvector     = Global.mySolver->tm2;
        Global.mySolver->tm2 = Global.mySolver->tm1;
        Global.mySolver->tm1 = tmpvector;

        Timer_Start( "Solver I/O" );
        solver_write_checkpoint( step, startingStep );
        solver_update_status( step, startingStep );
        solver_output_wavefield( step );
        solver_output_planes( Global.mySolver, Global.myID, step );
        solver_output_stations( step );
        solver_output_drm_nodes( Global.mySolver, step, Param.theTotalSteps );
        solver_read_source_forces( step );
        solver_read_drm_displacements( step , Param.theDeltaT ,Param.theTotalSteps );
        Timer_Stop( "Solver I/O" );

        Timer_Start( "Compute Physics" );
        solver_nonlinear_state( Global.mySolver, Global.myMesh, Global.theK1, Global.theK2, step );
        solver_compute_force_source( step );
        solver_compute_effective_drm_force( Global.mySolver, Global.myMesh,Global.theK1, Global.theK2, step, Param.theDeltaT );
        solver_compute_force_stiffness( Global.mySolver, Global.myMesh, Global.theK1, Global.theK2 );
        solver_compute_force_damping( Global.mySolver, Global.myMesh, Global.theK1, Global.theK2 );
        solver_compute_force_gravity( Global.mySolver, Global.myMesh, step );
        solver_compute_force_nonlinear( Global.mySolver, Global.myMesh, Param.theDeltaTSquared );
        Timer_Stop( "Compute Physics" );

        Timer_Start( "Communication" );
        HU_COND_GLOBAL_BARRIER( Param.theTimingBarriersFlag );
        solver_send_force_dangling( Global.mySolver );
        solver_adjust_forces( Global.mySolver );
        HU_COND_GLOBAL_BARRIER( Param.theTimingBarriersFlag );
        solver_send_force_anchored( Global.mySolver );
        Timer_Stop( "Communication" );

        Timer_Start( "Compute Physics" );
        solver_compute_displacement( Global.mySolver, Global.myMesh );
        solver_geostatic_fix( step );
        solver_load_fixedbase_displacements( Global.mySolver, step );
        Timer_Stop( "Compute Physics" );

        Timer_Start( "Communication" );
        HU_COND_GLOBAL_BARRIER( Param.theTimingBarriersFlag );
        solver_send_displacement_anchored( Global.mySolver );
        solver_adjust_displacement( Global.mySolver );
        HU_COND_GLOBAL_BARRIER( Param.theTimingBarriersFlag );
        solver_send_displacement_dangling( Global.mySolver );
        Timer_Stop( "Communication" );

        solver_loop_hook_bottom( Global.mySolver, Global.myMesh, step );
    } /* for (step = ....): all steps */

    solver_drm_close();
    solver_output_wavefield_close();
    solver_run_collect_timers();
}


/**
 * Output the velocity of the mesh nodes for the current timestep. Send
 * the data to Thread 0 who is responsible for dumping the data to disk.
 *
 * \note This is only useful for very small distributed setups where:
 * - There is no parallel file system.
 * - The distributed file system does not have POSIX semantics.
 * - The distributed file system performs poorly with just a few clients
 *   and the effective overall performance is better when only one client
 *   (i.e., PE 0) is writing to the FS.
 */
void
solver_output_seq()
{
    int32_t nindex;
    int32_t batchlimit, idx;

#ifdef DEBUG
    int64_t gnid_prev, gnid_current;
    int32_t first_counted;
#endif /* DEBUG */

    batchlimit = BATCH * 10;

    /* Allocate a fixed size buffer space if not initiazlied */
    if (Global.myVelocityTable == NULL) {
	Global.myVelocityTable = (fvector_t *)calloc(batchlimit, sizeof(fvector_t));
	if (Global.myVelocityTable == NULL) {
	    fprintf(stderr,  "Thread %d: solver_output_seq: out of memory\n",
		    Global.myID);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}
    }

    if (Global.myID == 0) {
	int32_t procid;
#ifdef DEBUG
	first_counted = 0;
#endif

	if (Param.FourDOutFp == NULL) {
	    out_hdr_t out_hdr;

	    /* First output, create the output file */
	    Param.FourDOutFp = fopen(Param.FourDOutFile, "w+");
	    if (Param.FourDOutFp == NULL) {
		fprintf(stderr, "Thread 0: solver_output_seq: ");
		fprintf(stderr, "cannot create %s\n", Param.FourDOutFile);
		perror("fopen");
		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	    }

	    /* Write the header that contains the metadata */
	    out_hdr.domain_x = Param.theDomainX;
	    out_hdr.domain_y = Param.theDomainY;
	    out_hdr.domain_z = Param.theDomainZ;
	    out_hdr.total_nodes = Global.theNTotal;
	    out_hdr.total_elements = Global.theETotal;
	    out_hdr.mesh_ticksize = Global.myMesh->ticksize;
	    out_hdr.output_steps = (Param.theTotalSteps - 1) / Param.theRate + 1;

	    if (fwrite(&out_hdr, sizeof(out_hdr_t), 1, Param.FourDOutFp) != 1){
		fprintf(stderr, "Thread 0: solver_output_seq: ");
		fprintf(stderr, "fail to write 4D-out header info\n");
		perror("fwrite");
		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	    }
	}

	/* Output my owned nodes' velocities */
	nindex = 0;
	while (nindex < Global.myMesh->nharbored) {
	    fvector_t vel;

	    if (Global.myMesh->nodeTable[nindex].ismine) {
		vel.f[0] =
		    (Global.mySolver->tm1[nindex].f[0]
		     - Global.mySolver->tm2[nindex].f[0])  / Param.theDeltaT;
		vel.f[1] =
		    (Global.mySolver->tm1[nindex].f[1]
		     - Global.mySolver->tm2[nindex].f[1])  / Param.theDeltaT;
		vel.f[2] =
		    (Global.mySolver->tm1[nindex].f[2]
		     - Global.mySolver->tm2[nindex].f[2])  / Param.theDeltaT;


		if (fwrite(&vel, sizeof(fvector_t), 1, Param.FourDOutFp) != 1) {
		    fprintf(stderr, "Thread 0: solver_output_seq: error\n");
		    MPI_Abort(MPI_COMM_WORLD, ERROR);
		    exit(1);
		}

		Param.the4DOutSize += sizeof(fvector_t);

#ifdef DEBUG
		gnid_current = Global.myMesh->nodeTable[nindex].gnid;

		if (first_counted) {
		    if (gnid_prev != (gnid_current - 1)) {
			fprintf( stderr, "PE 0: uncontinuous gnid\n"
				 "   gnid_prev = %" INT64_FMT
				 ", gnid_current = %" INT64_FMT "\n",
				 gnid_prev, gnid_current);
		    }
		} else {
		    first_counted = 1;
		}

		gnid_prev = gnid_current;

		if ((vel.f[0] != 0) ||
		    (vel.f[1] != 0) ||
		    (vel.f[2] != 0)) {
		    /*
		    fprintf(stderr, "Thread 0: Node %ld	 non-zero values\n",
			    gnid_current);
		    */
		}
#endif /* DEBUG */

	    }

	    nindex++;
	}

	/* Receive data from other processors */
	for (procid = 1; procid < Global.theGroupSize; procid++) {
	    MPI_Status status;
	    int32_t rcvbytecount;

	    /* Signal the next processor to go ahead */
	    MPI_Send(NULL, 0, MPI_CHAR, procid, GOAHEAD_MSG, comm_solver);

	    while (1) {
		MPI_Probe(procid, OUT4D_MSG, comm_solver, &status);
		MPI_Get_count(&status, MPI_CHAR, &rcvbytecount);

		/* Receive the data even if rcvbytecount == 0. Otherwise
		   the 0-byte message would get stuck in the message queue */
		MPI_Recv(Global.myVelocityTable, rcvbytecount, MPI_CHAR, procid,
			 OUT4D_MSG, comm_solver, &status);

		if (rcvbytecount == 0) {
		    /* Done */
		    break;
		}

		if (fwrite(Global.myVelocityTable, rcvbytecount, 1, Param.FourDOutFp) != 1) {
		    fprintf(stderr, "Thread 0: solver_output_seq: error\n");
		    MPI_Abort(MPI_COMM_WORLD, ERROR);
		    exit(1);
		}

		Param.the4DOutSize += rcvbytecount;

	    } /* while there is more data to be received from procid */
	} /* for all the processors */

    } else {
	/* Processors other than 0 needs to send data to 0 */
	int32_t sndbytecount;
	MPI_Status status;

	/* Wait for me turn */
	MPI_Recv(NULL, 0, MPI_CHAR, 0, GOAHEAD_MSG, comm_solver, &status);

#ifdef DEBUG
	first_counted = 0;
#endif


	nindex = 0;
	while (nindex < Global.myMesh->nharbored) {
	    fvector_t *velp;

	    idx = 0;
	    while ((idx < batchlimit) &&
		   (nindex < Global.myMesh->nharbored)) {

		if (Global.myMesh->nodeTable[nindex].ismine) {

		    velp = &Global.myVelocityTable[idx];

		    velp->f[0] =
			(Global.mySolver->tm1[nindex].f[0]
			 - Global.mySolver->tm2[nindex].f[0])	/ Param.theDeltaT;
		    velp->f[1] =
			(Global.mySolver->tm1[nindex].f[1]
			 - Global.mySolver->tm2[nindex].f[1])	/ Param.theDeltaT;
		    velp->f[2] =
			(Global.mySolver->tm1[nindex].f[2]
			 - Global.mySolver->tm2[nindex].f[2])	/ Param.theDeltaT;


		    idx++;

#ifdef DEBUG
		    gnid_current = Global.myMesh->nodeTable[nindex].gnid;

		    if (first_counted) {
			if (gnid_prev != (gnid_current - 1)) {
			    fprintf( stderr, "PE %d uncontinuous gnid\n"
				     "	gnid_prev = %" INT64_FMT
				     ", gnid_current = %" INT64_FMT "\n",
				     Global.myID, gnid_prev, gnid_current );
			}
		    } else {
			first_counted = 1;
		    }

		    gnid_prev = gnid_current;

		    /* debug */
		    /*
		    if ((velp->f[0] != 0) ||
			(velp->f[1] != 0) ||
			(velp->f[2] != 0)) {
			fprintf(stderr,
				"Thread %d: there are non-zero values\n",
				Global.myID);
				}
		    */
#endif /* DEBUG */

		}

		nindex++;
	    }

	    /* Send data to proc 0 */

	    if (idx > 0) {
		/* I have some real data to send */
		sndbytecount = idx * sizeof(fvector_t);
		MPI_Send(Global.myVelocityTable, sndbytecount, MPI_CHAR, 0, OUT4D_MSG,
			 comm_solver);
	    }
	} /* While there is data left to be sent */

	/* Send an empty message to indicate the end of my transfer */
	MPI_Send(NULL, 0, MPI_CHAR, 0, OUT4D_MSG, comm_solver);
    }

    return;
}


/**
 * Allocate and initialize a scheduler.
 */
static schedule_t*
schedule_new()
{
    schedule_t *sched;

    sched = (schedule_t *)malloc(sizeof(schedule_t));
    if (sched == NULL)
	return NULL;

    sched->c_count = 0;
    sched->first_c = NULL;
    sched->messenger_c = (messenger_t **)
	calloc(Global.theGroupSize, sizeof(messenger_t *));
    if (sched->messenger_c == NULL)
	return NULL;

    sched->s_count = 0;
    sched->first_s = NULL;
    sched->messenger_s = (messenger_t **)
	calloc(Global.theGroupSize, sizeof(messenger_t *));
    if (sched->messenger_s == NULL)
	return NULL;

    return sched;
}



/**
 * Allocate the mapping table for each of my messenger.
 */
static void
schedule_allocmapping( schedule_t *sched )
{
    messenger_t *messenger;
    int32_t nodecount;

    messenger = sched->first_c;
    while (messenger != NULL) {
	nodecount = messenger->nodecount;

	messenger->mapping =
	    (int32_t *)calloc(nodecount, sizeof(int32_t));

	if (messenger->mapping == NULL) {
	    fprintf(stderr, "Thread %d: schedule_allocamapping: ", Global.myID);
	    fprintf(stderr, " out of memory\n");
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	messenger = messenger->next;
    }

    messenger = sched->first_s;
    while (messenger != NULL) {
	nodecount = messenger->nodecount;

	messenger->mapping =
	    (int32_t *)calloc(nodecount, sizeof(int32_t));

	if (messenger->mapping == NULL) {
	    fprintf(stderr, "Thread %d: schedule_allocamapping: ", Global.myID);
	    fprintf(stderr, " out of memory\n");
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}

	messenger = messenger->next;
    }

    return;
}


/**
 * Allocate MPI controls for non-blocing receives.
 */
static void
schedule_allocMPIctl( schedule_t* sched )
{
    if (sched->c_count != 0) {
	sched->irecvreqs_c =
	    (MPI_Request *)malloc(sizeof(MPI_Request) * sched->c_count);
	sched->irecvstats_c =
	    (MPI_Status *)malloc(sizeof(MPI_Status) * sched->c_count);

	if ((sched->irecvreqs_c == NULL) ||
	    (sched->irecvstats_c == NULL)) {
	    fprintf(stderr, "Thread %d: schedule_allocMPIctl: ", Global.myID);
	    fprintf(stderr, "out of memory\n");
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}
    } else {
	sched->irecvreqs_c = NULL;
	sched->irecvstats_c = NULL;
    }

    if (sched->s_count != 0) {
	sched->irecvreqs_s =
	    (MPI_Request *)malloc(sizeof(MPI_Request) * sched->s_count);
	sched->irecvstats_s =
	    (MPI_Status *)malloc(sizeof(MPI_Status) * sched->s_count);

	if ((sched->irecvreqs_s == NULL) ||
	    (sched->irecvstats_s == NULL)) {
	    fprintf(stderr, "Thread %d: schedule_allocMPIctl: ", Global.myID);
	    fprintf(stderr, "out of memory\n");
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}
    } else {
	sched->irecvreqs_s = NULL;
	sched->irecvstats_s = NULL;
    }

    return;
}


/**
 * Build a communication schedule using local information.
 */
static void
schedule_build( mesh_t* mesh, schedule_t* dnsched, schedule_t* ansched )
{
    int32_t nindex;
    node_t* nodep;
    messenger_t* messenger;

    for (nindex = 0; nindex < mesh->nharbored; nindex++) {
	nodep = &mesh->nodeTable[nindex];
	int32_t owner, sharer;

	if (!nodep->ismine) {
	    /* I do not own this node. Add its owner processor to my c-list */

	    owner = nodep->proc.ownerid;

	    if (nodep->isanchored) {
		messenger = ansched->messenger_c[owner];
	    } else {
		messenger = dnsched->messenger_c[owner];
	    }

	    if (messenger == NULL) {
		messenger = messenger_new(owner);
		if (messenger == NULL) {
		    fprintf(stderr, "Thread %d: schedule_build: ", Global.myID);
		    fprintf(stderr, "out of memory.\n");
		    MPI_Abort(MPI_COMM_WORLD, ERROR);
		    exit(1);
		}

		if (nodep->isanchored) {
		    ansched->c_count++;
		    ansched->messenger_c[owner] = messenger;
		    messenger->next = ansched->first_c;
		    ansched->first_c = messenger;
		} else {
		    dnsched->c_count++;
		    dnsched->messenger_c[owner] = messenger;
		    messenger->next = dnsched->first_c;
		    dnsched->first_c = messenger;
		}
	    }

	    /* Update the number of nodecount for the messenger */
	    messenger->nodecount++;

	} else {
	    /* I own this node. Add any sharing processor to my s-list */

	    int32link_t *int32link;

	    int32link = nodep->proc.share;
	    while (int32link != NULL) {
		sharer = int32link->id;

		if (nodep->isanchored) {
		    messenger = ansched->messenger_s[sharer];
		} else {
		    messenger = dnsched->messenger_s[sharer];
		}

		if (messenger == NULL) {
		    messenger = messenger_new(sharer);
		    if (messenger == NULL) {
			fprintf(stderr, "Thread %d: schedule_build: ", Global.myID);
			fprintf(stderr, "out of memory.\n");
			MPI_Abort(MPI_COMM_WORLD, ERROR);
			exit(1);
		    }

		    if (nodep->isanchored) {
			ansched->s_count++;
			ansched->messenger_s[sharer] = messenger;
			messenger->next = ansched->first_s;
			ansched->first_s = messenger;
		    } else {
			dnsched->s_count++;
			dnsched->messenger_s[sharer] = messenger;
			messenger->next = dnsched->first_s;
			dnsched->first_s = messenger;
		    }
		}

		/* Update the nodecount */
		messenger->nodecount++;

		/* Move to the next sharing processor */
		int32link = int32link->next;
	    }
	}
    }

    /* Allocate MPI controls */
    schedule_allocMPIctl(ansched);
    schedule_allocMPIctl(dnsched);

    /* Allocate localnode table for each of the messegners I have */
    schedule_allocmapping(ansched);
    schedule_allocmapping(dnsched);

    /* Go through the nodes again and fill out the mapping table */
    for (nindex = 0; nindex < mesh->nharbored; nindex++) {
	nodep = &mesh->nodeTable[nindex];
	int32_t owner, sharer;

	if (!nodep->ismine) {
	    /* I do not own this node. Add its owner processor to my c-list */
	    owner = nodep->proc.ownerid;

	    if (nodep->isanchored) {
		messenger = ansched->messenger_c[owner];
	    } else {
		messenger = dnsched->messenger_c[owner];
	    }

	    if (messenger == NULL) {
		fprintf(stderr, "Thread %d: schedule_build: ", Global.myID);
		fprintf(stderr, "encounter NULL messenger.\n");
		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	    }

	    /* Fill in the mapping */
	    messenger->mapping[messenger->nidx] = nindex;
	    messenger->nidx++;

	} else {
	    /* I own this node. Add any sharing processor to my s-list */

	    int32link_t *int32link;

	    int32link = nodep->proc.share;
	    while (int32link != NULL) {
		sharer = int32link->id;

		if (nodep->isanchored) {
		    messenger = ansched->messenger_s[sharer];
		} else {
		    messenger = dnsched->messenger_s[sharer];
		}

		if (messenger == NULL) {
		    fprintf(stderr, "Thread %d: schedule_build: ", Global.myID);
		    fprintf(stderr, "encounter NULL messenger.\n");
		    MPI_Abort(MPI_COMM_WORLD, ERROR);
		    exit(1);
		}

		messenger->mapping[messenger->nidx] = nindex;
		messenger->nidx++;

		/* Move to the next sharing processor */
		int32link = int32link->next;
	    }
	}
    }

    return ;
}



/**
 * Release memory used by a scheduler.
 */
static void
schedule_delete( schedule_t* sched )
{
    messenger_t *current, *next;

    /* Release messengers overseeing my contributions */
    current = sched->first_c;
    while (current != NULL) {
	next = current->next;

	messenger_delete(current);
	current = next;
    }

    /* Release messengers overseeing shareing with others */
    current = sched->first_s;
    while (current != NULL) {
	next = current->next;

	messenger_delete(current);
	current = next;
    }

    if (sched->irecvreqs_c != NULL)
	free(sched->irecvreqs_c);
    if (sched->irecvstats_c != NULL)
	free(sched->irecvstats_c);
    free(sched->messenger_c);

    if (sched->irecvreqs_s != NULL)
	free(sched->irecvreqs_s);
    if (sched->irecvstats_s != NULL)
	free(sched->irecvstats_s);
    free(sched->messenger_s);

    free(sched);

    return;
}



/**
 * Allocate the memory for data exchange.
 */
static void
schedule_prepare( schedule_t* sched, int32_t c_outsize, int32_t c_insize,
		  int32_t s_outsize, int32_t s_insize )
{
    messenger_t* messenger;

    messenger = sched->first_c;
    while (messenger != NULL) {
	messenger_set(messenger, c_outsize, c_insize);
	messenger = messenger->next;
    }

    messenger = sched->first_s;
    while (messenger != NULL) {
	messenger_set(messenger, s_outsize, s_insize);
	messenger = messenger->next;
    }

    return;

}



/**
 * Assemble the proper information for the group of messengers managed by
 * a scheduler and send the data.
 *
 * \param direction: CONTRIBUTION or SHARING.
 */
static void
schedule_senddata(schedule_t *sched, void *valuetable, int32_t itemsperentry,
		  int32_t direction, int32_t msgtag)
{
    messenger_t *send_messenger, *recv_messenger;
    int32_t irecvcount, irecvnum, bytesize;
    MPI_Request *irecvreqs;
    MPI_Status *irecvstats;

#ifdef DEBUG
    int64_t *gnidp;
#endif /* DEBUG */

    if (direction == CONTRIBUTION) {
	send_messenger = sched->first_c;
	recv_messenger = sched->first_s;
	irecvcount = sched->s_count;
	irecvreqs = sched->irecvreqs_s;
	irecvstats = sched->irecvstats_s;
    } else {
	send_messenger = sched->first_s;
	recv_messenger = sched->first_c;
	irecvcount = sched->c_count;
	irecvreqs = sched->irecvreqs_c;
	irecvstats = sched->irecvstats_c;
    }

    /* Post receives */
    irecvnum = 0;
    while (recv_messenger != NULL) {
	bytesize = recv_messenger->nodecount * recv_messenger->insize;
	MPI_Irecv(recv_messenger->indata, bytesize, MPI_CHAR,
		  recv_messenger->procid, msgtag, comm_solver,
		  &irecvreqs[irecvnum]);

	irecvnum++;
	recv_messenger = recv_messenger->next;
    }

    /* Asssemble outgoing messages */
    while (send_messenger != NULL) {
	int32_t lnid, idx, entry;
	solver_float *dvalue;
	solver_float *out;

	for (idx = 0; idx < send_messenger->nidx; idx++) {

	    lnid = send_messenger->mapping[idx];

	    out = (solver_float *)((char *)send_messenger->outdata +
			     send_messenger->outsize * idx);

	    dvalue = (solver_float *)valuetable + itemsperentry * lnid;

	    for (entry = 0; entry < itemsperentry; entry++)
		*(out + entry) = *(dvalue + entry);

#ifdef DEBUG
	    /* For debug, carry the global node id */
	    gnidp = (int64_t *)
		((char *)out + itemsperentry * sizeof(solver_float));
	    *gnidp = Global.myMesh->nodeTable[lnid].gnid;
#endif /* DEBUG */
	}

	send_messenger = send_messenger->next;
    }

    /* Revisit messengers */
    if (direction == CONTRIBUTION) {
	send_messenger = sched->first_c;
	recv_messenger = sched->first_s;
    } else {
	send_messenger = sched->first_s;
	recv_messenger = sched->first_c;
    }

    /* Send the data */
    while (send_messenger != NULL) {
	bytesize = send_messenger->nodecount * send_messenger->outsize;
	MPI_Send(send_messenger->outdata, bytesize, MPI_CHAR,
		 send_messenger->procid, msgtag, comm_solver);
	send_messenger = send_messenger->next;
    }

    /* Wait till I receive all the data I want */
    if (irecvcount != 0) {
	MPI_Waitall(irecvcount, irecvreqs, irecvstats);
    }

    while (recv_messenger != NULL) {
	int32_t lnid, idx, entry;
	solver_float *dvalue;
	solver_float *in;

	for (idx = 0; idx < recv_messenger->nidx; idx++) {

	    lnid = recv_messenger->mapping[idx];

	    in = (solver_float *)((char *)recv_messenger->indata +
			    recv_messenger->insize * idx);

	    dvalue = (solver_float *)valuetable + itemsperentry * lnid;

	    for (entry = 0; entry < itemsperentry; entry++) {
		if (direction == CONTRIBUTION) {
		    *(dvalue + entry) += *(in + entry);
		} else {
		    /* SHARING, overwrite my local value */
		    *(dvalue + entry) = *(in + entry);
		}
	    }

#ifdef DEBUG
	    /* For debug, check the global node id */
	    gnidp = (int64_t *)
		((char *)in + itemsperentry * sizeof(solver_float));

	    if (*gnidp != Global.myMesh->nodeTable[lnid].gnid) {
		fprintf(stderr, "Thread %d: solver_init: gnids do not match\n",
			Global.myID);
		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	    }
#endif /* DEBUG */
	} /* for all the incoming data */

	recv_messenger = recv_messenger->next;
    }


    /* MPI_Barrier(comm_solver);  */

    return;
}


/**
 * Print messenger entry data to the given output stream.
 * The output line contains the following information:
 * "pe_id type c/s rpe nodecount outsize insize"
 *
 * \note Printing is not synchronized accross PEs, so lines can be out of
 * order / interleaved
 */
static int
messenger_print( messenger_t* msg, char type, char cs, FILE* out )
{
    int ret;

    assert( NULL != msg );

    ret = fprintf( out, "msg_info: %d %c %c %d %d %d %d %d\n", Global.myID, type, cs,
		   msg->procid, msg->nodecount, msg->outsize, msg->insize,
		   msg->nodecount * msg->outsize );

    if (ret > 0 || !Param.theSchedulePrintErrorCheckFlag) {
	ret = 0;
    } else {  /* ret <= 0 && Param.theSchedulePrintErrorCheckFlag */
	perror( "Warning! fprintf(...) failed" );
	fprintf(stderr, "PE id = %d, ret = %d, out = 0x%p\n", Global.myID, ret, out);
	fflush( stderr );
    }

    return ret;
}


/**
 * Print the messenger list for a scheduler
 *
 * \note Printing is not synchronized accross PEs, so lines can be out of
 * order / interleaved
 */
static int
schedule_print_messenger_list( schedule_t* sched, messenger_t* msg, int count,
			       char type, char cs, FILE* out )
{
    messenger_t* m_p;
    int my_count;
    int ret = 0;

    assert( NULL != sched );

    m_p = msg;


    for (my_count = 0; m_p != NULL && ret >= 0; my_count++, m_p = m_p->next) {
	ret += messenger_print( msg, type, cs, out );
    }


    if (ret < 0) {
	return -1;
    }

    if (my_count != count) {
	fprintf( stderr, "Warning! schedule list count does not match: "
		 "%u %c %c %d %d\n", Global.myID, type, cs, count, my_count );
    }

    return ret;
}


/**
 * Print scheduler communication to a given output stream.
 *
 * Print a line per entry in the schedule_t structure, i.e., an line
 * per messenger_t entry in the list.
 *
 * Each line has the following format:
 * "pe_id type c/s rpe nodecount outsize insize"
 */
static int
schedule_print_detail( schedule_t* sched, char type, FILE* out )
{
    int ret = 0;

    ret += schedule_print_messenger_list( sched, sched->first_c,
					  sched->c_count, type, 'c', out );
    ret += schedule_print_messenger_list( sched, sched->first_s,
					  sched->s_count, type, 's', out );

    return ret;
}


/**
 * Print scheduler communication to a given output stream.
 *
 * Print a line per entry in the schedule_t structure, i.e., an line
 * per messenger_t entry in the list.
 *
 * Each line has the following format:
 * "pe_id type s_count c_count"
 *
 * \note Printing is not synchronized accross PEs, so lines can be out of
 * order / interleaved
 */
static int
schedule_print( schedule_t* sched, char type, FILE* out )
{
    int ret;

    assert( NULL != sched );
    assert( NULL != out );

    ret = fprintf( out, "sch_info: %u %c %d %d\n", Global.myID, type, sched->c_count,
		   sched->s_count );


    if (ret > 0 || !Param.theSchedulePrintErrorCheckFlag) {
	ret = 0;
    }

    return ret;
}


/**
 * \note Printing is not synchronized accross PEs, so lines can be out of
 * order / interleaved
 */
static int
solver_print_schedules_imp( mysolver_t* solver, FILE* out )
{
    int ret = 0;

    assert( solver != NULL );
    assert( out != NULL );

    MPI_Barrier(comm_solver);

    /* print the high-level schedule per PE */
    if (Global.myID == 0) { /* print some header information */
	fputs( "# ----------------------------------------------------------\n"
	       "# Content: Solver schedule information\n"
	       "# pe_id d/a s_count c_count\n", out );
    }
    fflush( out );
    fdatasync( fileno( out ) );

    MPI_Barrier( comm_solver );

    /* this is not synchronized, so it might be out of order */
    ret += schedule_print( solver->an_sched, 'a', out );
    ret += schedule_print( solver->dn_sched, 'd', out );

    fflush( out );
    fdatasync( fileno( out ) );
    MPI_Barrier( comm_solver );

    /* print the schedule detail */
    if (Global.myID == 0) { /* print some header information */
	fputs( "\n\n"
	       "# ----------------------------------------------------------\n"
	       "# Content: Solver schedule detail\n"
	       "# pe_id d/a c/s rpe nodecount outsize insize msgsize\n", out );
	fflush( out );
	fdatasync( fileno( out ) );
    }

    MPI_Barrier( comm_solver );

    /* this is not synchronized, so it might be out of order */
    ret += schedule_print_detail( solver->an_sched, 'a', out );
    ret += schedule_print_detail( solver->dn_sched, 'd', out );

    fflush( out );
    fdatasync( fileno( out ) );
    MPI_Barrier( comm_solver );

    if (Global.myID == 0) {
	fputs( "# ----------------------------------------------------------\n"
	       "\n", out );
	fflush( out );
	fdatasync( fileno( out ) );
    }

    MPI_Barrier( comm_solver );

    return ret;
}


/**
 * Wrapper to print the solver schedules to stdout and a file with
 * the given name.
 *
 * \note Printing is not synchronized accross PEs, so lines can be out of
 * order / interleaved
 */
static int
solver_print_schedules( mysolver_t* solver )
{
    FILE* out;
    int ret = 0;

    if (Param.theSchedulePrintToStdout) {
	/* print schedules to standard output */
	ret += solver_print_schedules_imp( solver, stdout );

	if (ret < 0) {
	    fprintf( stderr, "Warning! printing schedules to standard output "
		     "failed for PE %d\n", Global.myID );
	}
    }

    if (Param.theSchedulePrintToFile && (Param.theSchedulePrintFilename != NULL)) {
	/* print schedules to the given file */
	out = fopen( Param.theSchedulePrintFilename, "a" );

	if (NULL == out) { /* this is not fatal */
	    fprintf( stderr, "Warning!, PE# %d failed to open output file for "
		     "printing the communication schedule\n", Global.myID );
	    return -1;
	}

	ret = solver_print_schedules_imp( solver, out );

	if (ret < 0) {
	    fprintf( stderr, "Warning! PE %d could not print schedules "
		     "to file\n", Global.myID );
	}

	fclose( out );
    }

    return ret;
}

/**
 * messenger_new: Allocate and initialize a messenger.
 *
 */
static messenger_t *messenger_new(int32_t procid)
{
    messenger_t *messenger;

    messenger = (messenger_t *)calloc(1, sizeof(messenger_t));

    if (messenger == NULL)
	return NULL;

    messenger->procid = procid;

    return messenger;
}



/**
 * messenger_delete: Release memory used by a messenger.
 *
 */
static void messenger_delete(messenger_t *messenger)
{
    if (messenger == NULL)
	return;

    if (messenger->outdata != NULL)
	free(messenger->outdata);

    if (messenger->indata != NULL)
	free(messenger->indata);

    free(messenger->mapping);

    free(messenger);

    return;
}



/**
 * messenger_set: Free any data memory used by the messenger for
 *                previous communication. And allocate new memory
 *                for the new round of communication.
 *
 */
static void
messenger_set(messenger_t *messenger, int32_t outsize, int32_t insize)
{
    if (messenger->outdata != NULL) {
	free(messenger->outdata);
	messenger->outdata = NULL;
    }

    if (messenger->indata != NULL) {
	free(messenger->indata);
	messenger->indata = NULL;
    }

    messenger->outsize = outsize;
    messenger->insize = insize;

    if (outsize != 0) {
	messenger->outdata = calloc(messenger->nodecount, outsize);
	if (messenger->outdata == NULL) {
	    fprintf(stderr, "Thread %d: messenger_set: out of memory\n",
		    Global.myID);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}
    }

    if (insize != 0) {
	messenger->indata = calloc(messenger->nodecount, insize);
	if (messenger->indata == NULL) {
	    fprintf(stderr, "Thread %d: messenger_set: out of memory\n",
		    Global.myID);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}
    }

    return;
}



/**
 * messenger_countnodes: Count the total number of nodes (that need to
 *                       be communicated) on the messenger link list.
 *
 */
static int32_t
messenger_countnodes(messenger_t *first)
{
    messenger_t *messenger;
    int32_t total;

    total = 0;
    messenger = first;
    while (messenger != NULL) {
	total += messenger->nodecount;
	messenger = messenger->next;
    }

    return total;
}


/**
 * Compute K1, K2 and K3 for a  linear  elastic,  homogeneous
 * and isotropic cube; and their off-diagonal matrices.
 *
 *            K1 = Integral(   Del(Phi(j))  * T(Del(Phi(i)))     )
 *            K2 = Integral(   Del(Phi(i))  * T(Del(Phi(j)))     )
 *            K3 = Integral( T(Del(Phi(i))) *   Del(Phi(j))  * I )
 *
 *            where:
 *
 *            T(W) = transpose of W.
 *
 * \note K offdiagonal has been commented out because is no longer needed
 * after the changes made on the solution algorith due to the new
 * form to handle the damping involved terms.
 */
static void compute_K()
{
    int i, j; /* indices of the block matrices, i for rows, j for columns */
    int k, l; /* indices of 3 x 3 matrices, k for rows, l for columns     */

    double x[3][8]={ {-1,  1, -1,  1, -1,  1, -1, 1} ,
		     {-1, -1,  1,  1, -1, -1,  1, 1} ,
		     {-1, -1, -1, -1,  1,  1,  1, 1} };

    /* compute K3 first */
    memset(Global.theK3, 0, 8 * 8 * sizeof(fmatrix_t));
    for (i = 0 ;  i < 8; i++) {
	for (j = 0; j <8; j++){
	    /* for  each 3 x 3 matrix representing (node i, node j),
               each of these matrices is diagonal, off-diagonal elements
               have been set to 0 in memset() */

	    fmatrix_t *matPtr;
	    matPtr = &Global.theK3[i][j];

	    for (k = 0; k < 3; k++) {
		/* set the diagonal values for the 3 x 3 matrix.
		   Note the rotation of the index */
		double I1, I2, I3;

		I1 = INTEGRAL_1
		    (x[(k + 0) % 3][i], x[(k + 0) % 3][j],
		     x[(k + 1) % 3][i], x[(k + 1) % 3][j],
		     x[(k + 2) % 3][i], x[(k + 2) % 3][j]);

		I2 = INTEGRAL_1
		    (x[(k + 1) % 3][i], x[(k + 1) % 3][j],
		     x[(k + 2) % 3][i], x[(k + 2) % 3][j],
		     x[(k + 0) % 3][i], x[(k + 0) % 3][j]);


		I3 = INTEGRAL_1
		    (x[(k + 2) % 3][i], x[(k + 2) % 3][j],
		     x[(k + 0) % 3][i], x[(k + 0) % 3][j],
		     x[(k + 1) % 3][i], x[(k + 1) % 3][j]);

		matPtr->f[k][k] = I1 + I2 + I3;
	    } /* k*/
	} /* j*/
    } /* i */

    /* compute K1 and K2. They are not diagonal either, but they are
       indeed symmetric globablly */
    for (i = 0; i < 8; i++) {
	for (j = 0; j < 8; j++) {
	    /* for  each 3 x 3 matrix representing (node i, node j)*/

	    fmatrix_t *matPtr0, *matPtr1;
	    matPtr0 = &Global.theK1[i][j];
	    matPtr1 = &Global.theK2[i][j];

	    for (k = 0; k <  3; k++)
		for (l = 0; l < 3; l++) {
		    if (k == l) {
			/* Initialize the diagnoal elements */
			matPtr0->f[k][k] =
			    INTEGRAL_1
			    (x[(k + 0) % 3][i], x[(k + 0) % 3][j],
			     x[(k + 1) % 3][i], x[(k + 1) % 3][j],
			     x[(k + 2) % 3][i], x[(k + 2) % 3][j]);


			matPtr1->f[k][k] =
			    INTEGRAL_1
			    (x[(k + 0) % 3][j], x[(k + 0) % 3][i],
			     x[(k + 1) % 3][j], x[(k + 1) % 3][i],
			     x[(k + 2) % 3][j], x[(k + 2) % 3][i]);

		    } else {
			/* Initialize off-diagonal elements */

			matPtr0->f[k][l] =
			    INTEGRAL_2
			    (x[k][j], x[l][i],
			     x[3 - (k + l)][j], x[3 - (k + l)][i]);

			matPtr1->f[k][l] =
			    INTEGRAL_2
			    (x[k][i], x[l][j],
			     x[3 - (k + l)][i], x[3 - (k + l)][j]);

		    } /* else */
		} /* for l */
	} /* for j */
    } /* i */

    /* Old code for damping                                    */
    /* Create the off-diagonal matrices                        */
    /* memcpy(theK1_offdiag, Global.theK1, sizeof(double) * 24 * 24); */
    /* memcpy(theK2_offdiag, Global.theK2, sizeof(double) * 24 * 24); */
    /* memcpy(theK3_offdiag, Global.theK3, sizeof(double) * 24 * 24); */
    /* for (i = 0; i <8; i++) {                                */
    /*     theK1_offdiag[i][i].f[0][0] = 0;                    */
    /*     theK1_offdiag[i][i].f[1][1] = 0;                    */
    /*     theK1_offdiag[i][i].f[2][2] = 0;                    */
    /*     theK2_offdiag[i][i].f[0][0] = 0;                    */
    /*     theK2_offdiag[i][i].f[1][1] = 0;                    */
    /*     theK2_offdiag[i][i].f[2][2] = 0;                    */
    /*     theK3_offdiag[i][i].f[0][0] = 0;                    */
    /*     theK3_offdiag[i][i].f[1][1] = 0;                    */
    /*     theK3_offdiag[i][i].f[2][2] = 0;                    */
    /* }                                                       */

    /* New code:
     * First option to solve double K1 K3 computation in the
     * time steps is to merge both of them in K1 only
     */
    for ( i = 0; i < 8; i++ )
    {
	for ( j = 0; j < 8; j++ )
	{
	    for ( k = 0; k < 3; k++ )
	    {
		for ( l = 0; l < 3; l++ )
		{
		    Global.theK1[i][j].f[k][l] += Global.theK3[i][j].f[k][l];
		}
	    }
	}
    }

    return;
}

static void constract_Quality_Factor_Table()
{
int i,j;
double local_QTABLE[26][6] = {{ 5.,	0.211111102, 0.236842104, 0.032142857, 0.271428571,	0.14},
		{6.25,	0.188888889,	0.184210526,	0.039893617,	0.336879433,	0.10152},
		{8.33,	0.157777778,	0.139473684,	0.045,	0.38,	0.07},
		{10., 0.137777765, 0.12105263, 0.032942899, 0.27818448, 0.0683},
		{15., 0.097777765,	0.08105263,	0.032942899,	0.27818448,	0.045},
		{20., 0.078139527, 0.060526314,	0.031409788, 0.277574872, 0.034225},
		{25., 0.064285708, 0.049999999,	0.031578947, 0.285714286, 0.0266},
		{30.,	0.053658537,	0.044736842,	0.026640676,	0.24691358,	0.023085},
		{35.,	0.046341463,	0.038157895,	0.02709848,	0.251156642,	0.019669},
		{40.,	0.040487805,	0.034210526,	0.025949367,	0.240506329,	0.01738},
		{45.,	0.036585366,	0.028947368,	0.031393568,	0.290964778,	0.014366},
		{50.,	0.032926829,	0.026315789,	0.032488114,	0.30110935,	0.01262},
		{60.,     0.0279,    0.0223,    0.0275,    0.2545,    0.0114},
		{70.,   0.024,			0.019,			0.032488114,    0.30110935, 0.0083},
		{80.,  0.0207,    0.0174,    0.0251,    0.2326,    0.0088},
		{90.,    0.0187,    0.0154,    0.0244,    0.2256,    0.0079},
		{100.,	0.017,	0.014,	0.028021016,	0.288966725,	0.006281},
		{120.,     0.0142,    0.0115,    0.0280,   0.2700,    0.0052},
		{150.,  0.0114,    0.0094,    0.0240,    0.2316,    0.0047},
		{200.,	0.0085,	0.00705,	0.022603978,	0.226039783,	0.0035392},
		{250., 0.0069,    0.0055,    0.0269,    0.2596,    0.0027},
		{300.,	0.0057,	0.0047,	0.027072758,	0.279187817,	0.0021276},
		{350,  0.0048,    0.0040,    0.0242,    0.2339,    0.0020},
		{400.,	0.0043,	0.0036,	0.021425572,	0.214255718,	0.0017935},
		{450., 0.0039,    0.0030,   0.0280,    0.2710,    0.0015},
		{500.,	0.0035,	0.00285,	0.023408925,	0.241404535,	0.001367}
};

for(i = 0; i < 18; i++)
{
	for(j = 0; j < 6; j++)
	{
		Global.theQTABLE[i][j] = local_QTABLE[i][j];
//		printf("%f ",theQTABLE[i][j]);
	}
//	printf("\n");
}
return;
}


int VModel_query(double depth, cvmpayload_t * g_props)

{
	int i;
	double local_1DVelocityModel[1551][6] = {
			{	0	,	1.7	,	0.5	,	2	,	50	,	100	}	,
			{	0.02	,	1.98	,	0.78	,	2.08	,	78	,	156	}	,
			{	0.04	,	2.26	,	1.06	,	2.16	,	106	,	212	}	,
			{	0.06	,	2.425	,	1.2125	,	2.20625	,	121.25	,	242.5	}	,
			{	0.08	,	2.475	,	1.2375	,	2.21875	,	123.75	,	247.5	}	,
			{	0.1	,	2.525	,	1.2625	,	2.23125	,	126.25	,	252.5	}	,
			{	0.12	,	2.575	,	1.2875	,	2.24375	,	128.75	,	257.5	}	,
			{	0.14	,	2.625	,	1.3125	,	2.25625	,	131.25	,	262.5	}	,
			{	0.16	,	2.675	,	1.3375	,	2.26875	,	133.75	,	267.5	}	,
			{	0.18	,	2.725	,	1.3625	,	2.28125	,	136.25	,	272.5	}	,
			{	0.2	,	2.775	,	1.3875	,	2.29375	,	138.75	,	277.5	}	,
			{	0.22	,	2.81364	,	1.40909	,	2.30455	,	140.90909	,	281.81818	}	,
			{	0.24	,	2.84091	,	1.42727	,	2.31364	,	142.72727	,	285.45455	}	,
			{	0.26	,	2.86818	,	1.44545	,	2.32273	,	144.54545	,	289.09091	}	,
			{	0.28	,	2.89545	,	1.46364	,	2.33182	,	146.36364	,	292.72727	}	,
			{	0.3	,	2.92273	,	1.48182	,	2.34091	,	148.18182	,	296.36364	}	,
			{	0.32	,	2.95	,	1.5	,	2.35	,	150	,	300	}	,
			{	0.34	,	2.97727	,	1.51818	,	2.35909	,	151.81818	,	303.63636	}	,
			{	0.36	,	3.00455	,	1.53636	,	2.36818	,	153.63636	,	307.27273	}	,
			{	0.38	,	3.03182	,	1.55455	,	2.37727	,	155.45455	,	310.90909	}	,
			{	0.4	,	3.05909	,	1.57273	,	2.38636	,	157.27273	,	314.54545	}	,
			{	0.42	,	3.08636	,	1.59091	,	2.39545	,	159.09091	,	318.18182	}	,
			{	0.44	,	3.11364	,	1.60909	,	2.40227	,	160.90909	,	321.81818	}	,
			{	0.46	,	3.14091	,	1.62727	,	2.40682	,	162.72727	,	325.45455	}	,
			{	0.48	,	3.16818	,	1.64545	,	2.41136	,	164.54545	,	329.09091	}	,
			{	0.5	,	3.19545	,	1.66364	,	2.41591	,	166.36364	,	332.72727	}	,
			{	0.52	,	3.22273	,	1.68182	,	2.42045	,	168.18182	,	336.36364	}	,
			{	0.54	,	3.25	,	1.7	,	2.425	,	170	,	340	}	,
			{	0.56	,	3.27727	,	1.71818	,	2.42955	,	171.81818	,	343.63636	}	,
			{	0.58	,	3.30455	,	1.73636	,	2.43409	,	173.63636	,	347.27273	}	,
			{	0.6	,	3.33182	,	1.75455	,	2.43864	,	175.45455	,	350.90909	}	,
			{	0.62	,	3.35909	,	1.77273	,	2.44318	,	177.27273	,	354.54545	}	,
			{	0.64	,	3.38636	,	1.79091	,	2.44773	,	179.09091	,	358.18182	}	,
			{	0.66	,	3.41111	,	1.81111	,	2.45185	,	181.11111	,	362.22222	}	,
			{	0.68	,	3.43333	,	1.83333	,	2.45556	,	183.33333	,	366.66667	}	,
			{	0.7	,	3.45556	,	1.85556	,	2.45926	,	185.55556	,	371.11111	}	,
			{	0.72	,	3.47778	,	1.87778	,	2.46296	,	187.77778	,	375.55556	}	,
			{	0.74	,	3.5	,	1.9	,	2.46667	,	190	,	380	}	,
			{	0.76	,	3.52222	,	1.92222	,	2.47037	,	192.22222	,	384.44444	}	,
			{	0.78	,	3.54444	,	1.94444	,	2.47407	,	194.44444	,	388.88889	}	,
			{	0.8	,	3.56667	,	1.96667	,	2.47778	,	196.66667	,	393.33333	}	,
			{	0.82	,	3.58889	,	1.98889	,	2.48148	,	198.88889	,	397.77778	}	,
			{	0.84	,	3.61111	,	2.01111	,	2.48519	,	201.11111	,	402.22222	}	,
			{	0.86	,	3.63333	,	2.03333	,	2.48889	,	203.33333	,	406.66667	}	,
			{	0.88	,	3.65556	,	2.05556	,	2.49259	,	205.55556	,	411.11111	}	,
			{	0.9	,	3.67778	,	2.07778	,	2.4963	,	207.77778	,	415.55556	}	,
			{	0.92	,	3.7	,	2.1	,	2.5	,	210	,	420	}	,
			{	0.94	,	3.71197	,	2.10513	,	2.50171	,	210.51282	,	421.02564	}	,
			{	0.96	,	3.72393	,	2.11026	,	2.50342	,	211.02564	,	422.05128	}	,
			{	0.98	,	3.7359	,	2.11538	,	2.50513	,	211.53846	,	423.07692	}	,
			{	1	,	3.74786	,	2.12051	,	2.50684	,	212.05128	,	424.10256	}	,
			{	1.02	,	3.75983	,	2.12564	,	2.50855	,	212.5641	,	425.12821	}	,
			{	1.04	,	3.77179	,	2.13077	,	2.51026	,	213.07692	,	426.15385	}	,
			{	1.06	,	3.78376	,	2.1359	,	2.51197	,	213.58974	,	427.17949	}	,
			{	1.08	,	3.79573	,	2.14103	,	2.51368	,	214.10256	,	428.20513	}	,
			{	1.1	,	3.80769	,	2.14615	,	2.51538	,	214.61538	,	429.23077	}	,
			{	1.12	,	3.81966	,	2.15128	,	2.51709	,	215.12821	,	430.25641	}	,
			{	1.14	,	3.83162	,	2.15641	,	2.5188	,	215.64103	,	431.28205	}	,
			{	1.16	,	3.84359	,	2.16154	,	2.52051	,	216.15385	,	432.30769	}	,
			{	1.18	,	3.85556	,	2.16667	,	2.52222	,	216.66667	,	433.33333	}	,
			{	1.2	,	3.86752	,	2.17179	,	2.52393	,	217.17949	,	434.35897	}	,
			{	1.22	,	3.87949	,	2.17692	,	2.52564	,	217.69231	,	435.38462	}	,
			{	1.24	,	3.89145	,	2.18205	,	2.52735	,	218.20513	,	436.41026	}	,
			{	1.26	,	3.90342	,	2.18718	,	2.52906	,	218.71795	,	437.4359	}	,
			{	1.28	,	3.91538	,	2.19231	,	2.53077	,	219.23077	,	438.46154	}	,
			{	1.3	,	3.92735	,	2.19744	,	2.53248	,	219.74359	,	439.48718	}	,
			{	1.32	,	3.93932	,	2.20256	,	2.53419	,	220.25641	,	440.51282	}	,
			{	1.34	,	3.95128	,	2.20769	,	2.5359	,	220.76923	,	441.53846	}	,
			{	1.36	,	3.96325	,	2.21282	,	2.53761	,	221.28205	,	442.5641	}	,
			{	1.38	,	3.97521	,	2.21795	,	2.53932	,	221.79487	,	443.58974	}	,
			{	1.4	,	3.98718	,	2.22308	,	2.54103	,	222.30769	,	444.61538	}	,
			{	1.42	,	3.99915	,	2.22821	,	2.54274	,	222.82051	,	445.64103	}	,
			{	1.44	,	4.01111	,	2.23333	,	2.54444	,	223.33333	,	446.66667	}	,
			{	1.46	,	4.02308	,	2.23846	,	2.54615	,	223.84615	,	447.69231	}	,
			{	1.48	,	4.03504	,	2.24359	,	2.54786	,	224.35897	,	448.71795	}	,
			{	1.5	,	4.04701	,	2.24872	,	2.54957	,	224.87179	,	449.74359	}	,
			{	1.52	,	4.05897	,	2.25385	,	2.55128	,	225.38462	,	450.76923	}	,
			{	1.54	,	4.07094	,	2.25897	,	2.55299	,	225.89744	,	451.79487	}	,
			{	1.56	,	4.08291	,	2.2641	,	2.5547	,	226.41026	,	452.82051	}	,
			{	1.58	,	4.09487	,	2.26923	,	2.55641	,	226.92308	,	453.84615	}	,
			{	1.6	,	4.10684	,	2.27436	,	2.55812	,	227.4359	,	454.87179	}	,
			{	1.62	,	4.1188	,	2.27949	,	2.55983	,	227.94872	,	455.89744	}	,
			{	1.64	,	4.13077	,	2.28462	,	2.56154	,	228.46154	,	456.92308	}	,
			{	1.66	,	4.14274	,	2.28974	,	2.56325	,	228.97436	,	457.94872	}	,
			{	1.68	,	4.1547	,	2.29487	,	2.56496	,	229.48718	,	458.97436	}	,
			{	1.7	,	4.16667	,	2.3	,	2.56667	,	230	,	460	}	,
			{	1.72	,	4.17863	,	2.30513	,	2.56838	,	230.51282	,	461.02564	}	,
			{	1.74	,	4.1906	,	2.31026	,	2.57009	,	231.02564	,	462.05128	}	,
			{	1.76	,	4.20256	,	2.31538	,	2.57179	,	231.53846	,	463.07692	}	,
			{	1.78	,	4.21453	,	2.32051	,	2.5735	,	232.05128	,	464.10256	}	,
			{	1.8	,	4.2265	,	2.32564	,	2.57521	,	232.5641	,	465.12821	}	,
			{	1.82	,	4.23846	,	2.33077	,	2.57692	,	233.07692	,	466.15385	}	,
			{	1.84	,	4.25043	,	2.3359	,	2.57863	,	233.58974	,	467.17949	}	,
			{	1.86	,	4.26239	,	2.34103	,	2.58034	,	234.10256	,	468.20513	}	,
			{	1.88	,	4.27436	,	2.34615	,	2.58205	,	234.61538	,	469.23077	}	,
			{	1.9	,	4.28632	,	2.35128	,	2.58376	,	235.12821	,	470.25641	}	,
			{	1.92	,	4.29829	,	2.35641	,	2.58547	,	235.64103	,	471.28205	}	,
			{	1.94	,	4.31026	,	2.36154	,	2.58718	,	236.15385	,	472.30769	}	,
			{	1.96	,	4.32222	,	2.36667	,	2.58889	,	236.66667	,	473.33333	}	,
			{	1.98	,	4.33419	,	2.37179	,	2.5906	,	237.17949	,	474.35897	}	,
			{	2	,	4.34615	,	2.37692	,	2.59231	,	237.69231	,	475.38462	}	,
			{	2.02	,	4.35812	,	2.38205	,	2.59402	,	238.20513	,	476.41026	}	,
			{	2.04	,	4.37009	,	2.38718	,	2.59573	,	238.71795	,	477.4359	}	,
			{	2.06	,	4.38205	,	2.39231	,	2.59744	,	239.23077	,	478.46154	}	,
			{	2.08	,	4.39402	,	2.39744	,	2.59915	,	239.74359	,	479.48718	}	,
			{	2.1	,	4.40347	,	2.40198	,	2.6005	,	240.19802	,	480.39604	}	,
			{	2.12	,	4.4104	,	2.40594	,	2.60149	,	240.59406	,	481.18812	}	,
			{	2.14	,	4.41733	,	2.4099	,	2.60248	,	240.9901	,	481.9802	}	,
			{	2.16	,	4.42426	,	2.41386	,	2.60347	,	241.38614	,	482.77228	}	,
			{	2.18	,	4.43119	,	2.41782	,	2.60446	,	241.78218	,	483.56436	}	,
			{	2.2	,	4.43812	,	2.42178	,	2.60545	,	242.17822	,	484.35644	}	,
			{	2.22	,	4.44505	,	2.42574	,	2.60644	,	242.57426	,	485.14851	}	,
			{	2.24	,	4.45198	,	2.4297	,	2.60743	,	242.9703	,	485.94059	}	,
			{	2.26	,	4.45891	,	2.43366	,	2.60842	,	243.36634	,	486.73267	}	,
			{	2.28	,	4.46584	,	2.43762	,	2.60941	,	243.76238	,	487.52475	}	,
			{	2.3	,	4.47277	,	2.44158	,	2.6104	,	244.15842	,	488.31683	}	,
			{	2.32	,	4.4797	,	2.44554	,	2.61139	,	244.55446	,	489.10891	}	,
			{	2.34	,	4.48663	,	2.4495	,	2.61238	,	244.9505	,	489.90099	}	,
			{	2.36	,	4.49356	,	2.45347	,	2.61337	,	245.34653	,	490.69307	}	,
			{	2.38	,	4.5005	,	2.45743	,	2.61436	,	245.74257	,	491.48515	}	,
			{	2.4	,	4.50743	,	2.46139	,	2.61535	,	246.13861	,	492.27723	}	,
			{	2.42	,	4.51436	,	2.46535	,	2.61634	,	246.53465	,	493.06931	}	,
			{	2.44	,	4.52129	,	2.46931	,	2.61733	,	246.93069	,	493.86139	}	,
			{	2.46	,	4.52822	,	2.47327	,	2.61832	,	247.32673	,	494.65347	}	,
			{	2.48	,	4.53515	,	2.47723	,	2.61931	,	247.72277	,	495.44554	}	,
			{	2.5	,	4.54208	,	2.48119	,	2.6203	,	248.11881	,	496.23762	}	,
			{	2.52	,	4.54901	,	2.48515	,	2.62129	,	248.51485	,	497.0297	}	,
			{	2.54	,	4.55594	,	2.48911	,	2.62228	,	248.91089	,	497.82178	}	,
			{	2.56	,	4.56287	,	2.49307	,	2.62327	,	249.30693	,	498.61386	}	,
			{	2.58	,	4.5698	,	2.49703	,	2.62426	,	249.70297	,	499.40594	}	,
			{	2.6	,	4.57673	,	2.50099	,	2.62525	,	250.09901	,	500.19802	}	,
			{	2.62	,	4.58366	,	2.50495	,	2.62624	,	250.49505	,	500.9901	}	,
			{	2.64	,	4.59059	,	2.50891	,	2.62723	,	250.89109	,	501.78218	}	,
			{	2.66	,	4.59752	,	2.51287	,	2.62822	,	251.28713	,	502.57426	}	,
			{	2.68	,	4.60446	,	2.51683	,	2.62921	,	251.68317	,	503.36634	}	,
			{	2.7	,	4.61139	,	2.52079	,	2.6302	,	252.07921	,	504.15842	}	,
			{	2.72	,	4.61832	,	2.52475	,	2.63119	,	252.47525	,	504.9505	}	,
			{	2.74	,	4.62525	,	2.52871	,	2.63218	,	252.87129	,	505.74257	}	,
			{	2.76	,	4.63218	,	2.53267	,	2.63317	,	253.26733	,	506.53465	}	,
			{	2.78	,	4.63911	,	2.53663	,	2.63416	,	253.66337	,	507.32673	}	,
			{	2.8	,	4.64604	,	2.54059	,	2.63515	,	254.05941	,	508.11881	}	,
			{	2.82	,	4.65297	,	2.54455	,	2.63614	,	254.45545	,	508.91089	}	,
			{	2.84	,	4.6599	,	2.54851	,	2.63713	,	254.85149	,	509.70297	}	,
			{	2.86	,	4.66683	,	2.55248	,	2.63812	,	255.24752	,	510.49505	}	,
			{	2.88	,	4.67376	,	2.55644	,	2.63911	,	255.64356	,	511.28713	}	,
			{	2.9	,	4.68069	,	2.5604	,	2.6401	,	256.0396	,	512.07921	}	,
			{	2.92	,	4.68762	,	2.56436	,	2.64109	,	256.43564	,	512.87129	}	,
			{	2.94	,	4.69455	,	2.56832	,	2.64208	,	256.83168	,	513.66337	}	,
			{	2.96	,	4.70149	,	2.57228	,	2.64307	,	257.22772	,	514.45545	}	,
			{	2.98	,	4.70842	,	2.57624	,	2.64406	,	257.62376	,	515.24752	}	,
			{	3	,	4.71535	,	2.5802	,	2.64505	,	258.0198	,	516.0396	}	,
			{	3.02	,	4.72228	,	2.58416	,	2.64604	,	258.41584	,	516.83168	}	,
			{	3.04	,	4.72921	,	2.58812	,	2.64703	,	258.81188	,	517.62376	}	,
			{	3.06	,	4.73614	,	2.59208	,	2.64802	,	259.20792	,	518.41584	}	,
			{	3.08	,	4.74307	,	2.59604	,	2.64901	,	259.60396	,	519.20792	}	,
			{	3.1	,	4.75	,	2.6	,	2.65	,	260	,	520	}	,
			{	3.12	,	4.75693	,	2.60396	,	2.65099	,	260.39604	,	520.79208	}	,
			{	3.14	,	4.76386	,	2.60792	,	2.65198	,	260.79208	,	521.58416	}	,
			{	3.16	,	4.77079	,	2.61188	,	2.65297	,	261.18812	,	522.37624	}	,
			{	3.18	,	4.77772	,	2.61584	,	2.65396	,	261.58416	,	523.16832	}	,
			{	3.2	,	4.78465	,	2.6198	,	2.65495	,	261.9802	,	523.9604	}	,
			{	3.22	,	4.79158	,	2.62376	,	2.65594	,	262.37624	,	524.75248	}	,
			{	3.24	,	4.79851	,	2.62772	,	2.65693	,	262.77228	,	525.54455	}	,
			{	3.26	,	4.80545	,	2.63168	,	2.65792	,	263.16832	,	526.33663	}	,
			{	3.28	,	4.81238	,	2.63564	,	2.65891	,	263.56436	,	527.12871	}	,
			{	3.3	,	4.81931	,	2.6396	,	2.6599	,	263.9604	,	527.92079	}	,
			{	3.32	,	4.82624	,	2.64356	,	2.66089	,	264.35644	,	528.71287	}	,
			{	3.34	,	4.83317	,	2.64752	,	2.66188	,	264.75248	,	529.50495	}	,
			{	3.36	,	4.8401	,	2.65149	,	2.66287	,	265.14851	,	530.29703	}	,
			{	3.38	,	4.84703	,	2.65545	,	2.66386	,	265.54455	,	531.08911	}	,
			{	3.4	,	4.85396	,	2.65941	,	2.66485	,	265.94059	,	531.88119	}	,
			{	3.42	,	4.86089	,	2.66337	,	2.66584	,	266.33663	,	532.67327	}	,
			{	3.44	,	4.86782	,	2.66733	,	2.66683	,	266.73267	,	533.46535	}	,
			{	3.46	,	4.87475	,	2.67129	,	2.66782	,	267.12871	,	534.25743	}	,
			{	3.48	,	4.88168	,	2.67525	,	2.66881	,	267.52475	,	535.0495	}	,
			{	3.5	,	4.88861	,	2.67921	,	2.6698	,	267.92079	,	535.84158	}	,
			{	3.52	,	4.89554	,	2.68317	,	2.67079	,	268.31683	,	536.63366	}	,
			{	3.54	,	4.90248	,	2.68713	,	2.67178	,	268.71287	,	537.42574	}	,
			{	3.56	,	4.90941	,	2.69109	,	2.67277	,	269.10891	,	538.21782	}	,
			{	3.58	,	4.91634	,	2.69505	,	2.67376	,	269.50495	,	539.0099	}	,
			{	3.6	,	4.92327	,	2.69901	,	2.67475	,	269.90099	,	539.80198	}	,
			{	3.62	,	4.9302	,	2.70297	,	2.67574	,	270.29703	,	540.59406	}	,
			{	3.64	,	4.93713	,	2.70693	,	2.67673	,	270.69307	,	541.38614	}	,
			{	3.66	,	4.94406	,	2.71089	,	2.67772	,	271.08911	,	542.17822	}	,
			{	3.68	,	4.95099	,	2.71485	,	2.67871	,	271.48515	,	542.9703	}	,
			{	3.7	,	4.95792	,	2.71881	,	2.6797	,	271.88119	,	543.76238	}	,
			{	3.72	,	4.96485	,	2.72277	,	2.68069	,	272.27723	,	544.55446	}	,
			{	3.74	,	4.97178	,	2.72673	,	2.68168	,	272.67327	,	545.34653	}	,
			{	3.76	,	4.97871	,	2.73069	,	2.68267	,	273.06931	,	546.13861	}	,
			{	3.78	,	4.98564	,	2.73465	,	2.68366	,	273.46535	,	546.93069	}	,
			{	3.8	,	4.99257	,	2.73861	,	2.68465	,	273.86139	,	547.72277	}	,
			{	3.82	,	4.9995	,	2.74257	,	2.68564	,	274.25743	,	548.51485	}	,
			{	3.84	,	5.00644	,	2.74653	,	2.68663	,	274.65347	,	549.30693	}	,
			{	3.86	,	5.01337	,	2.7505	,	2.68762	,	275.0495	,	550.09901	}	,
			{	3.88	,	5.0203	,	2.75446	,	2.68861	,	275.44554	,	550.89109	}	,
			{	3.9	,	5.02723	,	2.75842	,	2.6896	,	275.84158	,	551.68317	}	,
			{	3.92	,	5.03416	,	2.76238	,	2.69059	,	276.23762	,	552.47525	}	,
			{	3.94	,	5.04109	,	2.76634	,	2.69158	,	276.63366	,	553.26733	}	,
			{	3.96	,	5.04802	,	2.7703	,	2.69257	,	277.0297	,	554.05941	}	,
			{	3.98	,	5.05495	,	2.77426	,	2.69356	,	277.42574	,	554.85149	}	,
			{	4	,	5.06188	,	2.77822	,	2.69455	,	277.82178	,	555.64356	}	,
			{	4.02	,	5.06881	,	2.78218	,	2.69554	,	278.21782	,	556.43564	}	,
			{	4.04	,	5.07574	,	2.78614	,	2.69653	,	278.61386	,	557.22772	}	,
			{	4.06	,	5.08267	,	2.7901	,	2.69752	,	279.0099	,	558.0198	}	,
			{	4.08	,	5.0896	,	2.79406	,	2.69851	,	279.40594	,	558.81188	}	,
			{	4.1	,	5.09653	,	2.79802	,	2.6995	,	279.80198	,	559.60396	}	,
			{	4.12	,	5.10329	,	2.8023	,	2.70033	,	280.23026	,	560.46053	}	,
			{	4.14	,	5.10987	,	2.80691	,	2.70099	,	280.69079	,	561.38158	}	,
			{	4.16	,	5.11645	,	2.81151	,	2.70164	,	281.15132	,	562.30263	}	,
			{	4.18	,	5.12303	,	2.81612	,	2.7023	,	281.61184	,	563.22368	}	,
			{	4.2	,	5.12961	,	2.82072	,	2.70296	,	282.07237	,	564.14474	}	,
			{	4.22	,	5.13618	,	2.82533	,	2.70362	,	282.53289	,	565.06579	}	,
			{	4.24	,	5.14276	,	2.82993	,	2.70428	,	282.99342	,	565.98684	}	,
			{	4.26	,	5.14934	,	2.83454	,	2.70493	,	283.45395	,	566.90789	}	,
			{	4.28	,	5.15592	,	2.83914	,	2.70559	,	283.91447	,	567.82895	}	,
			{	4.3	,	5.1625	,	2.84375	,	2.70625	,	284.375	,	568.75	}	,
			{	4.32	,	5.16908	,	2.84836	,	2.70691	,	284.83553	,	569.67105	}	,
			{	4.34	,	5.17566	,	2.85296	,	2.70757	,	285.29605	,	570.59211	}	,
			{	4.36	,	5.18224	,	2.85757	,	2.70822	,	285.75658	,	571.51316	}	,
			{	4.38	,	5.18882	,	2.86217	,	2.70888	,	286.21711	,	572.43421	}	,
			{	4.4	,	5.19539	,	2.86678	,	2.70954	,	286.67763	,	573.35526	}	,
			{	4.42	,	5.20197	,	2.87138	,	2.7102	,	287.13816	,	574.27632	}	,
			{	4.44	,	5.20855	,	2.87599	,	2.71086	,	287.59868	,	575.19737	}	,
			{	4.46	,	5.21513	,	2.88059	,	2.71151	,	288.05921	,	576.11842	}	,
			{	4.48	,	5.22171	,	2.8852	,	2.71217	,	288.51974	,	577.03947	}	,
			{	4.5	,	5.22829	,	2.8898	,	2.71283	,	288.98026	,	577.96053	}	,
			{	4.52	,	5.23487	,	2.89441	,	2.71349	,	289.44079	,	578.88158	}	,
			{	4.54	,	5.24145	,	2.89901	,	2.71414	,	289.90132	,	579.80263	}	,
			{	4.56	,	5.24803	,	2.90362	,	2.7148	,	290.36184	,	580.72368	}	,
			{	4.58	,	5.25461	,	2.90822	,	2.71546	,	290.82237	,	581.64474	}	,
			{	4.6	,	5.26118	,	2.91283	,	2.71612	,	291.28289	,	582.56579	}	,
			{	4.62	,	5.26776	,	2.91743	,	2.71678	,	291.74342	,	583.48684	}	,
			{	4.64	,	5.27434	,	2.92204	,	2.71743	,	292.20395	,	584.40789	}	,
			{	4.66	,	5.28092	,	2.92664	,	2.71809	,	292.66447	,	585.32895	}	,
			{	4.68	,	5.2875	,	2.93125	,	2.71875	,	293.125	,	586.25	}	,
			{	4.7	,	5.29408	,	2.93586	,	2.71941	,	293.58553	,	587.17105	}	,
			{	4.72	,	5.30066	,	2.94046	,	2.72007	,	294.04605	,	588.09211	}	,
			{	4.74	,	5.30724	,	2.94507	,	2.72072	,	294.50658	,	589.01316	}	,
			{	4.76	,	5.31382	,	2.94967	,	2.72138	,	294.96711	,	589.93421	}	,
			{	4.78	,	5.32039	,	2.95428	,	2.72204	,	295.42763	,	590.85526	}	,
			{	4.8	,	5.32697	,	2.95888	,	2.7227	,	295.88816	,	591.77632	}	,
			{	4.82	,	5.33355	,	2.96349	,	2.72336	,	296.34868	,	592.69737	}	,
			{	4.84	,	5.34013	,	2.96809	,	2.72401	,	296.80921	,	593.61842	}	,
			{	4.86	,	5.34671	,	2.9727	,	2.72467	,	297.26974	,	594.53947	}	,
			{	4.88	,	5.35329	,	2.9773	,	2.72533	,	297.73026	,	595.46053	}	,
			{	4.9	,	5.35987	,	2.98191	,	2.72599	,	298.19079	,	596.38158	}	,
			{	4.92	,	5.36645	,	2.98651	,	2.72664	,	298.65132	,	597.30263	}	,
			{	4.94	,	5.37303	,	2.99112	,	2.7273	,	299.11184	,	598.22368	}	,
			{	4.96	,	5.37961	,	2.99572	,	2.72796	,	299.57237	,	599.14474	}	,
			{	4.98	,	5.38618	,	3.00033	,	2.72862	,	300.03289	,	600.06579	}	,
			{	5	,	5.39276	,	3.00493	,	2.72928	,	300.49342	,	600.98684	}	,
			{	5.02	,	5.39934	,	3.00954	,	2.72993	,	300.95395	,	601.90789	}	,
			{	5.04	,	5.40592	,	3.01414	,	2.73059	,	301.41447	,	602.82895	}	,
			{	5.06	,	5.4125	,	3.01875	,	2.73125	,	301.875	,	603.75	}	,
			{	5.08	,	5.41908	,	3.02336	,	2.73191	,	302.33553	,	604.67105	}	,
			{	5.1	,	5.42566	,	3.02796	,	2.73257	,	302.79605	,	605.59211	}	,
			{	5.12	,	5.43224	,	3.03257	,	2.73322	,	303.25658	,	606.51316	}	,
			{	5.14	,	5.43882	,	3.03717	,	2.73388	,	303.71711	,	607.43421	}	,
			{	5.16	,	5.44539	,	3.04178	,	2.73454	,	304.17763	,	608.35526	}	,
			{	5.18	,	5.45197	,	3.04638	,	2.7352	,	304.63816	,	609.27632	}	,
			{	5.2	,	5.45855	,	3.05099	,	2.73586	,	305.09868	,	610.19737	}	,
			{	5.22	,	5.46513	,	3.05559	,	2.73651	,	305.55921	,	611.11842	}	,
			{	5.24	,	5.47171	,	3.0602	,	2.73717	,	306.01974	,	612.03947	}	,
			{	5.26	,	5.47829	,	3.0648	,	2.73783	,	306.48026	,	612.96053	}	,
			{	5.28	,	5.48487	,	3.06941	,	2.73849	,	306.94079	,	613.88158	}	,
			{	5.3	,	5.49145	,	3.07401	,	2.73914	,	307.40132	,	614.80263	}	,
			{	5.32	,	5.49803	,	3.07862	,	2.7398	,	307.86184	,	615.72368	}	,
			{	5.34	,	5.50461	,	3.08322	,	2.74046	,	308.32237	,	616.64474	}	,
			{	5.36	,	5.51118	,	3.08783	,	2.74112	,	308.78289	,	617.56579	}	,
			{	5.38	,	5.51776	,	3.09243	,	2.74178	,	309.24342	,	618.48684	}	,
			{	5.4	,	5.52434	,	3.09704	,	2.74243	,	309.70395	,	619.40789	}	,
			{	5.42	,	5.53092	,	3.10164	,	2.74309	,	310.16447	,	620.32895	}	,
			{	5.44	,	5.5375	,	3.10625	,	2.74375	,	310.625	,	621.25	}	,
			{	5.46	,	5.54408	,	3.11086	,	2.74441	,	311.08553	,	622.17105	}	,
			{	5.48	,	5.55066	,	3.11546	,	2.74507	,	311.54605	,	623.09211	}	,
			{	5.5	,	5.55724	,	3.12007	,	2.74572	,	312.00658	,	624.01316	}	,
			{	5.52	,	5.56382	,	3.12467	,	2.74638	,	312.46711	,	624.93421	}	,
			{	5.54	,	5.57039	,	3.12928	,	2.74704	,	312.92763	,	625.85526	}	,
			{	5.56	,	5.57697	,	3.13388	,	2.7477	,	313.38816	,	626.77632	}	,
			{	5.58	,	5.58355	,	3.13849	,	2.74836	,	313.84868	,	627.69737	}	,
			{	5.6	,	5.59013	,	3.14309	,	2.74901	,	314.30921	,	628.61842	}	,
			{	5.62	,	5.59671	,	3.1477	,	2.74967	,	314.76974	,	629.53947	}	,
			{	5.64	,	5.60182	,	3.15149	,	2.75025	,	315.14901	,	630.29801	}	,
			{	5.66	,	5.60546	,	3.15447	,	2.75075	,	315.44702	,	630.89404	}	,
			{	5.68	,	5.60911	,	3.15745	,	2.75124	,	315.74503	,	631.49007	}	,
			{	5.7	,	5.61275	,	3.16043	,	2.75174	,	316.04305	,	632.08609	}	,
			{	5.72	,	5.61639	,	3.16341	,	2.75224	,	316.34106	,	632.68212	}	,
			{	5.74	,	5.62003	,	3.16639	,	2.75273	,	316.63907	,	633.27815	}	,
			{	5.76	,	5.62368	,	3.16937	,	2.75323	,	316.93709	,	633.87417	}	,
			{	5.78	,	5.62732	,	3.17235	,	2.75373	,	317.2351	,	634.4702	}	,
			{	5.8	,	5.63096	,	3.17533	,	2.75422	,	317.53311	,	635.06623	}	,
			{	5.82	,	5.6346	,	3.17831	,	2.75472	,	317.83113	,	635.66225	}	,
			{	5.84	,	5.63825	,	3.18129	,	2.75522	,	318.12914	,	636.25828	}	,
			{	5.86	,	5.64189	,	3.18427	,	2.75571	,	318.42715	,	636.8543	}	,
			{	5.88	,	5.64553	,	3.18725	,	2.75621	,	318.72517	,	637.45033	}	,
			{	5.9	,	5.64917	,	3.19023	,	2.75671	,	319.02318	,	638.04636	}	,
			{	5.92	,	5.65281	,	3.19321	,	2.7572	,	319.32119	,	638.64238	}	,
			{	5.94	,	5.65646	,	3.19619	,	2.7577	,	319.61921	,	639.23841	}	,
			{	5.96	,	5.6601	,	3.19917	,	2.7582	,	319.91722	,	639.83444	}	,
			{	5.98	,	5.66374	,	3.20215	,	2.75869	,	320.21523	,	640.43046	}	,
			{	6	,	5.66738	,	3.20513	,	2.75919	,	320.51325	,	641.02649	}	,
			{	6.02	,	5.67103	,	3.20811	,	2.75969	,	320.81126	,	641.62252	}	,
			{	6.04	,	5.67467	,	3.21109	,	2.76018	,	321.10927	,	642.21854	}	,
			{	6.06	,	5.67831	,	3.21407	,	2.76068	,	321.40728	,	642.81457	}	,
			{	6.08	,	5.68195	,	3.21705	,	2.76118	,	321.7053	,	643.4106	}	,
			{	6.1	,	5.6856	,	3.22003	,	2.76167	,	322.00331	,	644.00662	}	,
			{	6.12	,	5.68924	,	3.22301	,	2.76217	,	322.30132	,	644.60265	}	,
			{	6.14	,	5.69288	,	3.22599	,	2.76267	,	322.59934	,	645.19868	}	,
			{	6.16	,	5.69652	,	3.22897	,	2.76316	,	322.89735	,	645.7947	}	,
			{	6.18	,	5.70017	,	3.23195	,	2.76366	,	323.19536	,	646.39073	}	,
			{	6.2	,	5.70381	,	3.23493	,	2.76416	,	323.49338	,	646.98675	}	,
			{	6.22	,	5.70745	,	3.23791	,	2.76465	,	323.79139	,	647.58278	}	,
			{	6.24	,	5.71109	,	3.24089	,	2.76515	,	324.0894	,	648.17881	}	,
			{	6.26	,	5.71474	,	3.24387	,	2.76565	,	324.38742	,	648.77483	}	,
			{	6.28	,	5.71838	,	3.24685	,	2.76614	,	324.68543	,	649.37086	}	,
			{	6.3	,	5.72202	,	3.24983	,	2.76664	,	324.98344	,	649.96689	}	,
			{	6.32	,	5.72566	,	3.25281	,	2.76714	,	325.28146	,	650.56291	}	,
			{	6.34	,	5.7293	,	3.25579	,	2.76763	,	325.57947	,	651.15894	}	,
			{	6.36	,	5.73295	,	3.25877	,	2.76813	,	325.87748	,	651.75497	}	,
			{	6.38	,	5.73659	,	3.26175	,	2.76863	,	326.1755	,	652.35099	}	,
			{	6.4	,	5.74023	,	3.26474	,	2.76912	,	326.47351	,	652.94702	}	,
			{	6.42	,	5.74387	,	3.26772	,	2.76962	,	326.77152	,	653.54305	}	,
			{	6.44	,	5.74752	,	3.2707	,	2.77012	,	327.06954	,	654.13907	}	,
			{	6.46	,	5.75116	,	3.27368	,	2.77061	,	327.36755	,	654.7351	}	,
			{	6.48	,	5.7548	,	3.27666	,	2.77111	,	327.66556	,	655.33113	}	,
			{	6.5	,	5.75844	,	3.27964	,	2.77161	,	327.96358	,	655.92715	}	,
			{	6.52	,	5.76209	,	3.28262	,	2.7721	,	328.26159	,	656.52318	}	,
			{	6.54	,	5.76573	,	3.2856	,	2.7726	,	328.5596	,	657.11921	}	,
			{	6.56	,	5.76937	,	3.28858	,	2.7731	,	328.85762	,	657.71523	}	,
			{	6.58	,	5.77301	,	3.29156	,	2.77359	,	329.15563	,	658.31126	}	,
			{	6.6	,	5.77666	,	3.29454	,	2.77409	,	329.45364	,	658.90728	}	,
			{	6.62	,	5.7803	,	3.29752	,	2.77459	,	329.75166	,	659.50331	}	,
			{	6.64	,	5.78394	,	3.3005	,	2.77508	,	330.04967	,	660.09934	}	,
			{	6.66	,	5.78758	,	3.30348	,	2.77558	,	330.34768	,	660.69536	}	,
			{	6.68	,	5.79123	,	3.30646	,	2.77608	,	330.6457	,	661.29139	}	,
			{	6.7	,	5.79487	,	3.30944	,	2.77657	,	330.94371	,	661.88742	}	,
			{	6.72	,	5.79851	,	3.31242	,	2.77707	,	331.24172	,	662.48344	}	,
			{	6.74	,	5.80215	,	3.3154	,	2.77757	,	331.53974	,	663.07947	}	,
			{	6.76	,	5.80579	,	3.31838	,	2.77806	,	331.83775	,	663.6755	}	,
			{	6.78	,	5.80944	,	3.32136	,	2.77856	,	332.13576	,	664.27152	}	,
			{	6.8	,	5.81308	,	3.32434	,	2.77906	,	332.43377	,	664.86755	}	,
			{	6.82	,	5.81672	,	3.32732	,	2.77955	,	332.73179	,	665.46358	}	,
			{	6.84	,	5.82036	,	3.3303	,	2.78005	,	333.0298	,	666.0596	}	,
			{	6.86	,	5.82401	,	3.33328	,	2.78055	,	333.32781	,	666.65563	}	,
			{	6.88	,	5.82765	,	3.33626	,	2.78104	,	333.62583	,	667.25166	}	,
			{	6.9	,	5.83129	,	3.33924	,	2.78154	,	333.92384	,	667.84768	}	,
			{	6.92	,	5.83493	,	3.34222	,	2.78204	,	334.22185	,	668.44371	}	,
			{	6.94	,	5.83858	,	3.3452	,	2.78253	,	334.51987	,	669.03974	}	,
			{	6.96	,	5.84222	,	3.34818	,	2.78303	,	334.81788	,	669.63576	}	,
			{	6.98	,	5.84586	,	3.35116	,	2.78353	,	335.11589	,	670.23179	}	,
			{	7	,	5.8495	,	3.35414	,	2.78402	,	335.41391	,	670.82781	}	,
			{	7.02	,	5.85315	,	3.35712	,	2.78452	,	335.71192	,	671.42384	}	,
			{	7.04	,	5.85679	,	3.3601	,	2.78502	,	336.00993	,	672.01987	}	,
			{	7.06	,	5.86043	,	3.36308	,	2.78551	,	336.30795	,	672.61589	}	,
			{	7.08	,	5.86407	,	3.36606	,	2.78601	,	336.60596	,	673.21192	}	,
			{	7.1	,	5.86772	,	3.36904	,	2.78651	,	336.90397	,	673.80795	}	,
			{	7.12	,	5.87136	,	3.37202	,	2.787	,	337.20199	,	674.40397	}	,
			{	7.14	,	5.875	,	3.375	,	2.7875	,	337.5	,	675	}	,
			{	7.16	,	5.87864	,	3.37798	,	2.788	,	337.79801	,	675.59603	}	,
			{	7.18	,	5.88228	,	3.38096	,	2.78849	,	338.09603	,	676.19205	}	,
			{	7.2	,	5.88593	,	3.38394	,	2.78899	,	338.39404	,	676.78808	}	,
			{	7.22	,	5.88957	,	3.38692	,	2.78949	,	338.69205	,	677.38411	}	,
			{	7.24	,	5.89321	,	3.3899	,	2.78998	,	338.99007	,	677.98013	}	,
			{	7.26	,	5.89685	,	3.39288	,	2.79048	,	339.28808	,	678.57616	}	,
			{	7.28	,	5.9005	,	3.39586	,	2.79098	,	339.58609	,	679.17219	}	,
			{	7.3	,	5.90414	,	3.39884	,	2.79147	,	339.88411	,	679.76821	}	,
			{	7.32	,	5.90778	,	3.40182	,	2.79197	,	340.18212	,	680.36424	}	,
			{	7.34	,	5.91142	,	3.4048	,	2.79247	,	340.48013	,	680.96026	}	,
			{	7.36	,	5.91507	,	3.40778	,	2.79296	,	340.77815	,	681.55629	}	,
			{	7.38	,	5.91871	,	3.41076	,	2.79346	,	341.07616	,	682.15232	}	,
			{	7.4	,	5.92235	,	3.41374	,	2.79396	,	341.37417	,	682.74834	}	,
			{	7.42	,	5.92599	,	3.41672	,	2.79445	,	341.67219	,	683.34437	}	,
			{	7.44	,	5.92964	,	3.4197	,	2.79495	,	341.9702	,	683.9404	}	,
			{	7.46	,	5.93328	,	3.42268	,	2.79545	,	342.26821	,	684.53642	}	,
			{	7.48	,	5.93692	,	3.42566	,	2.79594	,	342.56623	,	685.13245	}	,
			{	7.5	,	5.94056	,	3.42864	,	2.79644	,	342.86424	,	685.72848	}	,
			{	7.52	,	5.94421	,	3.43162	,	2.79694	,	343.16225	,	686.3245	}	,
			{	7.54	,	5.94785	,	3.4346	,	2.79743	,	343.46026	,	686.92053	}	,
			{	7.56	,	5.95149	,	3.43758	,	2.79793	,	343.75828	,	687.51656	}	,
			{	7.58	,	5.95513	,	3.44056	,	2.79843	,	344.05629	,	688.11258	}	,
			{	7.6	,	5.95877	,	3.44354	,	2.79892	,	344.3543	,	688.70861	}	,
			{	7.62	,	5.96242	,	3.44652	,	2.79942	,	344.65232	,	689.30464	}	,
			{	7.64	,	5.96606	,	3.4495	,	2.79992	,	344.95033	,	689.90066	}	,
			{	7.66	,	5.9697	,	3.45248	,	2.80041	,	345.24834	,	690.49669	}	,
			{	7.68	,	5.97334	,	3.45546	,	2.80091	,	345.54636	,	691.09272	}	,
			{	7.7	,	5.97699	,	3.45844	,	2.80141	,	345.84437	,	691.68874	}	,
			{	7.72	,	5.98063	,	3.46142	,	2.8019	,	346.14238	,	692.28477	}	,
			{	7.74	,	5.98427	,	3.4644	,	2.8024	,	346.4404	,	692.88079	}	,
			{	7.76	,	5.98791	,	3.46738	,	2.8029	,	346.73841	,	693.47682	}	,
			{	7.78	,	5.99156	,	3.47036	,	2.80339	,	347.03642	,	694.07285	}	,
			{	7.8	,	5.9952	,	3.47334	,	2.80389	,	347.33444	,	694.66887	}	,
			{	7.82	,	5.99884	,	3.47632	,	2.80439	,	347.63245	,	695.2649	}	,
			{	7.84	,	6.00248	,	3.4793	,	2.80488	,	347.93046	,	695.86093	}	,
			{	7.86	,	6.00613	,	3.48228	,	2.80538	,	348.22848	,	696.45695	}	,
			{	7.88	,	6.00977	,	3.48526	,	2.80588	,	348.52649	,	697.05298	}	,
			{	7.9	,	6.01341	,	3.48825	,	2.80637	,	348.8245	,	697.64901	}	,
			{	7.92	,	6.01705	,	3.49123	,	2.80687	,	349.12252	,	698.24503	}	,
			{	7.94	,	6.0207	,	3.49421	,	2.80737	,	349.42053	,	698.84106	}	,
			{	7.96	,	6.02434	,	3.49719	,	2.80786	,	349.71854	,	699.43709	}	,
			{	7.98	,	6.02798	,	3.50017	,	2.80836	,	350.01656	,	700.03311	}	,
			{	8	,	6.03162	,	3.50315	,	2.80886	,	350.31457	,	700.62914	}	,
			{	8.02	,	6.03526	,	3.50613	,	2.80935	,	350.61258	,	701.22517	}	,
			{	8.04	,	6.03891	,	3.50911	,	2.80985	,	350.9106	,	701.82119	}	,
			{	8.06	,	6.04255	,	3.51209	,	2.81035	,	351.20861	,	702.41722	}	,
			{	8.08	,	6.04619	,	3.51507	,	2.81084	,	351.50662	,	703.01325	}	,
			{	8.1	,	6.04983	,	3.51805	,	2.81134	,	351.80464	,	703.60927	}	,
			{	8.12	,	6.05348	,	3.52103	,	2.81184	,	352.10265	,	704.2053	}	,
			{	8.14	,	6.05712	,	3.52401	,	2.81233	,	352.40066	,	704.80132	}	,
			{	8.16	,	6.06076	,	3.52699	,	2.81283	,	352.69868	,	705.39735	}	,
			{	8.18	,	6.0644	,	3.52997	,	2.81333	,	352.99669	,	705.99338	}	,
			{	8.2	,	6.06805	,	3.53295	,	2.81382	,	353.2947	,	706.5894	}	,
			{	8.22	,	6.07169	,	3.53593	,	2.81432	,	353.59272	,	707.18543	}	,
			{	8.24	,	6.07533	,	3.53891	,	2.81482	,	353.89073	,	707.78146	}	,
			{	8.26	,	6.07897	,	3.54189	,	2.81531	,	354.18874	,	708.37748	}	,
			{	8.28	,	6.08262	,	3.54487	,	2.81581	,	354.48675	,	708.97351	}	,
			{	8.3	,	6.08626	,	3.54785	,	2.81631	,	354.78477	,	709.56954	}	,
			{	8.32	,	6.0899	,	3.55083	,	2.8168	,	355.08278	,	710.16556	}	,
			{	8.34	,	6.09354	,	3.55381	,	2.8173	,	355.38079	,	710.76159	}	,
			{	8.36	,	6.09719	,	3.55679	,	2.8178	,	355.67881	,	711.35762	}	,
			{	8.38	,	6.10083	,	3.55977	,	2.81829	,	355.97682	,	711.95364	}	,
			{	8.4	,	6.10447	,	3.56275	,	2.81879	,	356.27483	,	712.54967	}	,
			{	8.42	,	6.10811	,	3.56573	,	2.81929	,	356.57285	,	713.1457	}	,
			{	8.44	,	6.11175	,	3.56871	,	2.81978	,	356.87086	,	713.74172	}	,
			{	8.46	,	6.1154	,	3.57169	,	2.82028	,	357.16887	,	714.33775	}	,
			{	8.48	,	6.11904	,	3.57467	,	2.82078	,	357.46689	,	714.93377	}	,
			{	8.5	,	6.12268	,	3.57765	,	2.82127	,	357.7649	,	715.5298	}	,
			{	8.52	,	6.12632	,	3.58063	,	2.82177	,	358.06291	,	716.12583	}	,
			{	8.54	,	6.12997	,	3.58361	,	2.82227	,	358.36093	,	716.72185	}	,
			{	8.56	,	6.13361	,	3.58659	,	2.82276	,	358.65894	,	717.31788	}	,
			{	8.58	,	6.13725	,	3.58957	,	2.82326	,	358.95695	,	717.91391	}	,
			{	8.6	,	6.14089	,	3.59255	,	2.82376	,	359.25497	,	718.50993	}	,
			{	8.62	,	6.14454	,	3.59553	,	2.82425	,	359.55298	,	719.10596	}	,
			{	8.64	,	6.14818	,	3.59851	,	2.82475	,	359.85099	,	719.70199	}	,
			{	8.66	,	6.15034	,	3.6001	,	2.82505	,	360.00996	,	720.01992	}	,
			{	8.68	,	6.15102	,	3.6003	,	2.82515	,	360.02988	,	720.05976	}	,
			{	8.7	,	6.15169	,	3.6005	,	2.82525	,	360.0498	,	720.0996	}	,
			{	8.72	,	6.15237	,	3.6007	,	2.82535	,	360.06972	,	720.13944	}	,
			{	8.74	,	6.15305	,	3.6009	,	2.82545	,	360.08964	,	720.17928	}	,
			{	8.76	,	6.15373	,	3.6011	,	2.82555	,	360.10956	,	720.21912	}	,
			{	8.78	,	6.1544	,	3.60129	,	2.82565	,	360.12948	,	720.25896	}	,
			{	8.8	,	6.15508	,	3.60149	,	2.82575	,	360.1494	,	720.2988	}	,
			{	8.82	,	6.15576	,	3.60169	,	2.82585	,	360.16932	,	720.33865	}	,
			{	8.84	,	6.15643	,	3.60189	,	2.82595	,	360.18924	,	720.37849	}	,
			{	8.86	,	6.15711	,	3.60209	,	2.82605	,	360.20916	,	720.41833	}	,
			{	8.88	,	6.15779	,	3.60229	,	2.82615	,	360.22908	,	720.45817	}	,
			{	8.9	,	6.15847	,	3.60249	,	2.82625	,	360.249	,	720.49801	}	,
			{	8.92	,	6.15914	,	3.60269	,	2.82634	,	360.26892	,	720.53785	}	,
			{	8.94	,	6.15982	,	3.60289	,	2.82644	,	360.28884	,	720.57769	}	,
			{	8.96	,	6.1605	,	3.60309	,	2.82654	,	360.30876	,	720.61753	}	,
			{	8.98	,	6.16118	,	3.60329	,	2.82664	,	360.32869	,	720.65737	}	,
			{	9	,	6.16185	,	3.60349	,	2.82674	,	360.34861	,	720.69721	}	,
			{	9.02	,	6.16253	,	3.60369	,	2.82684	,	360.36853	,	720.73705	}	,
			{	9.04	,	6.16321	,	3.60388	,	2.82694	,	360.38845	,	720.77689	}	,
			{	9.06	,	6.16388	,	3.60408	,	2.82704	,	360.40837	,	720.81673	}	,
			{	9.08	,	6.16456	,	3.60428	,	2.82714	,	360.42829	,	720.85657	}	,
			{	9.1	,	6.16524	,	3.60448	,	2.82724	,	360.44821	,	720.89641	}	,
			{	9.12	,	6.16592	,	3.60468	,	2.82734	,	360.46813	,	720.93625	}	,
			{	9.14	,	6.16659	,	3.60488	,	2.82744	,	360.48805	,	720.9761	}	,
			{	9.16	,	6.16727	,	3.60508	,	2.82754	,	360.50797	,	721.01594	}	,
			{	9.18	,	6.16795	,	3.60528	,	2.82764	,	360.52789	,	721.05578	}	,
			{	9.2	,	6.16863	,	3.60548	,	2.82774	,	360.54781	,	721.09562	}	,
			{	9.22	,	6.1693	,	3.60568	,	2.82784	,	360.56773	,	721.13546	}	,
			{	9.24	,	6.16998	,	3.60588	,	2.82794	,	360.58765	,	721.1753	}	,
			{	9.26	,	6.17066	,	3.60608	,	2.82804	,	360.60757	,	721.21514	}	,
			{	9.28	,	6.17133	,	3.60627	,	2.82814	,	360.62749	,	721.25498	}	,
			{	9.3	,	6.17201	,	3.60647	,	2.82824	,	360.64741	,	721.29482	}	,
			{	9.32	,	6.17269	,	3.60667	,	2.82834	,	360.66733	,	721.33466	}	,
			{	9.34	,	6.17337	,	3.60687	,	2.82844	,	360.68725	,	721.3745	}	,
			{	9.36	,	6.17404	,	3.60707	,	2.82854	,	360.70717	,	721.41434	}	,
			{	9.38	,	6.17472	,	3.60727	,	2.82864	,	360.72709	,	721.45418	}	,
			{	9.4	,	6.1754	,	3.60747	,	2.82874	,	360.74701	,	721.49402	}	,
			{	9.42	,	6.17608	,	3.60767	,	2.82883	,	360.76693	,	721.53386	}	,
			{	9.44	,	6.17675	,	3.60787	,	2.82893	,	360.78685	,	721.57371	}	,
			{	9.46	,	6.17743	,	3.60807	,	2.82903	,	360.80677	,	721.61355	}	,
			{	9.48	,	6.17811	,	3.60827	,	2.82913	,	360.82669	,	721.65339	}	,
			{	9.5	,	6.17878	,	3.60847	,	2.82923	,	360.84661	,	721.69323	}	,
			{	9.52	,	6.17946	,	3.60867	,	2.82933	,	360.86653	,	721.73307	}	,
			{	9.54	,	6.18014	,	3.60886	,	2.82943	,	360.88645	,	721.77291	}	,
			{	9.56	,	6.18082	,	3.60906	,	2.82953	,	360.90637	,	721.81275	}	,
			{	9.58	,	6.18149	,	3.60926	,	2.82963	,	360.92629	,	721.85259	}	,
			{	9.6	,	6.18217	,	3.60946	,	2.82973	,	360.94622	,	721.89243	}	,
			{	9.62	,	6.18285	,	3.60966	,	2.82983	,	360.96614	,	721.93227	}	,
			{	9.64	,	6.18353	,	3.60986	,	2.82993	,	360.98606	,	721.97211	}	,
			{	9.66	,	6.1842	,	3.61006	,	2.83003	,	361.00598	,	722.01195	}	,
			{	9.68	,	6.18488	,	3.61026	,	2.83013	,	361.0259	,	722.05179	}	,
			{	9.7	,	6.18556	,	3.61046	,	2.83023	,	361.04582	,	722.09163	}	,
			{	9.72	,	6.18624	,	3.61066	,	2.83033	,	361.06574	,	722.13147	}	,
			{	9.74	,	6.18691	,	3.61086	,	2.83043	,	361.08566	,	722.17131	}	,
			{	9.76	,	6.18759	,	3.61106	,	2.83053	,	361.10558	,	722.21116	}	,
			{	9.78	,	6.18827	,	3.61125	,	2.83063	,	361.1255	,	722.251	}	,
			{	9.8	,	6.18894	,	3.61145	,	2.83073	,	361.14542	,	722.29084	}	,
			{	9.82	,	6.18962	,	3.61165	,	2.83083	,	361.16534	,	722.33068	}	,
			{	9.84	,	6.1903	,	3.61185	,	2.83093	,	361.18526	,	722.37052	}	,
			{	9.86	,	6.19098	,	3.61205	,	2.83103	,	361.20518	,	722.41036	}	,
			{	9.88	,	6.19165	,	3.61225	,	2.83113	,	361.2251	,	722.4502	}	,
			{	9.9	,	6.19233	,	3.61245	,	2.83123	,	361.24502	,	722.49004	}	,
			{	9.92	,	6.19301	,	3.61265	,	2.83132	,	361.26494	,	722.52988	}	,
			{	9.94	,	6.19369	,	3.61285	,	2.83142	,	361.28486	,	722.56972	}	,
			{	9.96	,	6.19436	,	3.61305	,	2.83152	,	361.30478	,	722.60956	}	,
			{	9.98	,	6.19504	,	3.61325	,	2.83162	,	361.3247	,	722.6494	}	,
			{	10	,	6.19572	,	3.61345	,	2.83172	,	361.34462	,	722.68924	}	,
			{	10.02	,	6.19639	,	3.61365	,	2.83182	,	361.36454	,	722.72908	}	,
			{	10.04	,	6.19707	,	3.61384	,	2.83192	,	361.38446	,	722.76892	}	,
			{	10.06	,	6.19775	,	3.61404	,	2.83202	,	361.40438	,	722.80876	}	,
			{	10.08	,	6.19843	,	3.61424	,	2.83212	,	361.4243	,	722.84861	}	,
			{	10.1	,	6.1991	,	3.61444	,	2.83222	,	361.44422	,	722.88845	}	,
			{	10.12	,	6.19978	,	3.61464	,	2.83232	,	361.46414	,	722.92829	}	,
			{	10.14	,	6.20046	,	3.61484	,	2.83242	,	361.48406	,	722.96813	}	,
			{	10.16	,	6.20114	,	3.61504	,	2.83252	,	361.50398	,	723.00797	}	,
			{	10.18	,	6.20181	,	3.61524	,	2.83262	,	361.5239	,	723.04781	}	,
			{	10.2	,	6.20249	,	3.61544	,	2.83272	,	361.54382	,	723.08765	}	,
			{	10.22	,	6.20317	,	3.61564	,	2.83282	,	361.56375	,	723.12749	}	,
			{	10.24	,	6.20384	,	3.61584	,	2.83292	,	361.58367	,	723.16733	}	,
			{	10.26	,	6.20452	,	3.61604	,	2.83302	,	361.60359	,	723.20717	}	,
			{	10.28	,	6.2052	,	3.61624	,	2.83312	,	361.62351	,	723.24701	}	,
			{	10.3	,	6.20588	,	3.61643	,	2.83322	,	361.64343	,	723.28685	}	,
			{	10.32	,	6.20655	,	3.61663	,	2.83332	,	361.66335	,	723.32669	}	,
			{	10.34	,	6.20723	,	3.61683	,	2.83342	,	361.68327	,	723.36653	}	,
			{	10.36	,	6.20791	,	3.61703	,	2.83352	,	361.70319	,	723.40637	}	,
			{	10.38	,	6.20859	,	3.61723	,	2.83362	,	361.72311	,	723.44622	}	,
			{	10.4	,	6.20926	,	3.61743	,	2.83372	,	361.74303	,	723.48606	}	,
			{	10.42	,	6.20994	,	3.61763	,	2.83381	,	361.76295	,	723.5259	}	,
			{	10.44	,	6.21062	,	3.61783	,	2.83391	,	361.78287	,	723.56574	}	,
			{	10.46	,	6.21129	,	3.61803	,	2.83401	,	361.80279	,	723.60558	}	,
			{	10.48	,	6.21197	,	3.61823	,	2.83411	,	361.82271	,	723.64542	}	,
			{	10.5	,	6.21265	,	3.61843	,	2.83421	,	361.84263	,	723.68526	}	,
			{	10.52	,	6.21333	,	3.61863	,	2.83431	,	361.86255	,	723.7251	}	,
			{	10.54	,	6.214	,	3.61882	,	2.83441	,	361.88247	,	723.76494	}	,
			{	10.56	,	6.21468	,	3.61902	,	2.83451	,	361.90239	,	723.80478	}	,
			{	10.58	,	6.21536	,	3.61922	,	2.83461	,	361.92231	,	723.84462	}	,
			{	10.6	,	6.21604	,	3.61942	,	2.83471	,	361.94223	,	723.88446	}	,
			{	10.62	,	6.21671	,	3.61962	,	2.83481	,	361.96215	,	723.9243	}	,
			{	10.64	,	6.21739	,	3.61982	,	2.83491	,	361.98207	,	723.96414	}	,
			{	10.66	,	6.21807	,	3.62002	,	2.83501	,	362.00199	,	724.00398	}	,
			{	10.68	,	6.21875	,	3.62022	,	2.83511	,	362.02191	,	724.04382	}	,
			{	10.7	,	6.21942	,	3.62042	,	2.83521	,	362.04183	,	724.08367	}	,
			{	10.72	,	6.2201	,	3.62062	,	2.83531	,	362.06175	,	724.12351	}	,
			{	10.74	,	6.22078	,	3.62082	,	2.83541	,	362.08167	,	724.16335	}	,
			{	10.76	,	6.22145	,	3.62102	,	2.83551	,	362.10159	,	724.20319	}	,
			{	10.78	,	6.22213	,	3.62122	,	2.83561	,	362.12151	,	724.24303	}	,
			{	10.8	,	6.22281	,	3.62141	,	2.83571	,	362.14143	,	724.28287	}	,
			{	10.82	,	6.22349	,	3.62161	,	2.83581	,	362.16135	,	724.32271	}	,
			{	10.84	,	6.22416	,	3.62181	,	2.83591	,	362.18127	,	724.36255	}	,
			{	10.86	,	6.22484	,	3.62201	,	2.83601	,	362.2012	,	724.40239	}	,
			{	10.88	,	6.22552	,	3.62221	,	2.83611	,	362.22112	,	724.44223	}	,
			{	10.9	,	6.2262	,	3.62241	,	2.83621	,	362.24104	,	724.48207	}	,
			{	10.92	,	6.22687	,	3.62261	,	2.8363	,	362.26096	,	724.52191	}	,
			{	10.94	,	6.22755	,	3.62281	,	2.8364	,	362.28088	,	724.56175	}	,
			{	10.96	,	6.22823	,	3.62301	,	2.8365	,	362.3008	,	724.60159	}	,
			{	10.98	,	6.2289	,	3.62321	,	2.8366	,	362.32072	,	724.64143	}	,
			{	11	,	6.22958	,	3.62341	,	2.8367	,	362.34064	,	724.68127	}	,
			{	11.02	,	6.23026	,	3.62361	,	2.8368	,	362.36056	,	724.72112	}	,
			{	11.04	,	6.23094	,	3.6238	,	2.8369	,	362.38048	,	724.76096	}	,
			{	11.06	,	6.23161	,	3.624	,	2.837	,	362.4004	,	724.8008	}	,
			{	11.08	,	6.23229	,	3.6242	,	2.8371	,	362.42032	,	724.84064	}	,
			{	11.1	,	6.23297	,	3.6244	,	2.8372	,	362.44024	,	724.88048	}	,
			{	11.12	,	6.23365	,	3.6246	,	2.8373	,	362.46016	,	724.92032	}	,
			{	11.14	,	6.23432	,	3.6248	,	2.8374	,	362.48008	,	724.96016	}	,
			{	11.16	,	6.235	,	3.625	,	2.8375	,	362.5	,	725	}	,
			{	11.18	,	6.23568	,	3.6252	,	2.8376	,	362.51992	,	725.03984	}	,
			{	11.2	,	6.23635	,	3.6254	,	2.8377	,	362.53984	,	725.07968	}	,
			{	11.22	,	6.23703	,	3.6256	,	2.8378	,	362.55976	,	725.11952	}	,
			{	11.24	,	6.23771	,	3.6258	,	2.8379	,	362.57968	,	725.15936	}	,
			{	11.26	,	6.23839	,	3.626	,	2.838	,	362.5996	,	725.1992	}	,
			{	11.28	,	6.23906	,	3.6262	,	2.8381	,	362.61952	,	725.23904	}	,
			{	11.3	,	6.23974	,	3.62639	,	2.8382	,	362.63944	,	725.27888	}	,
			{	11.32	,	6.24042	,	3.62659	,	2.8383	,	362.65936	,	725.31873	}	,
			{	11.34	,	6.2411	,	3.62679	,	2.8384	,	362.67928	,	725.35857	}	,
			{	11.36	,	6.24177	,	3.62699	,	2.8385	,	362.6992	,	725.39841	}	,
			{	11.38	,	6.24245	,	3.62719	,	2.8386	,	362.71912	,	725.43825	}	,
			{	11.4	,	6.24313	,	3.62739	,	2.8387	,	362.73904	,	725.47809	}	,
			{	11.42	,	6.2438	,	3.62759	,	2.83879	,	362.75896	,	725.51793	}	,
			{	11.44	,	6.24448	,	3.62779	,	2.83889	,	362.77888	,	725.55777	}	,
			{	11.46	,	6.24516	,	3.62799	,	2.83899	,	362.7988	,	725.59761	}	,
			{	11.48	,	6.24584	,	3.62819	,	2.83909	,	362.81873	,	725.63745	}	,
			{	11.5	,	6.24651	,	3.62839	,	2.83919	,	362.83865	,	725.67729	}	,
			{	11.52	,	6.24719	,	3.62859	,	2.83929	,	362.85857	,	725.71713	}	,
			{	11.54	,	6.24787	,	3.62878	,	2.83939	,	362.87849	,	725.75697	}	,
			{	11.56	,	6.24855	,	3.62898	,	2.83949	,	362.89841	,	725.79681	}	,
			{	11.58	,	6.24922	,	3.62918	,	2.83959	,	362.91833	,	725.83665	}	,
			{	11.6	,	6.2499	,	3.62938	,	2.83969	,	362.93825	,	725.87649	}	,
			{	11.62	,	6.25058	,	3.62958	,	2.83979	,	362.95817	,	725.91633	}	,
			{	11.64	,	6.25125	,	3.62978	,	2.83989	,	362.97809	,	725.95618	}	,
			{	11.66	,	6.25193	,	3.62998	,	2.83999	,	362.99801	,	725.99602	}	,
			{	11.68	,	6.25261	,	3.63018	,	2.84009	,	363.01793	,	726.03586	}	,
			{	11.7	,	6.25329	,	3.63038	,	2.84019	,	363.03785	,	726.0757	}	,
			{	11.72	,	6.25396	,	3.63058	,	2.84029	,	363.05777	,	726.11554	}	,
			{	11.74	,	6.25464	,	3.63078	,	2.84039	,	363.07769	,	726.15538	}	,
			{	11.76	,	6.25532	,	3.63098	,	2.84049	,	363.09761	,	726.19522	}	,
			{	11.78	,	6.256	,	3.63118	,	2.84059	,	363.11753	,	726.23506	}	,
			{	11.8	,	6.25667	,	3.63137	,	2.84069	,	363.13745	,	726.2749	}	,
			{	11.82	,	6.25735	,	3.63157	,	2.84079	,	363.15737	,	726.31474	}	,
			{	11.84	,	6.25803	,	3.63177	,	2.84089	,	363.17729	,	726.35458	}	,
			{	11.86	,	6.25871	,	3.63197	,	2.84099	,	363.19721	,	726.39442	}	,
			{	11.88	,	6.25938	,	3.63217	,	2.84109	,	363.21713	,	726.43426	}	,
			{	11.9	,	6.26006	,	3.63237	,	2.84119	,	363.23705	,	726.4741	}	,
			{	11.92	,	6.26074	,	3.63257	,	2.84128	,	363.25697	,	726.51394	}	,
			{	11.94	,	6.26141	,	3.63277	,	2.84138	,	363.27689	,	726.55378	}	,
			{	11.96	,	6.26209	,	3.63297	,	2.84148	,	363.29681	,	726.59363	}	,
			{	11.98	,	6.26277	,	3.63317	,	2.84158	,	363.31673	,	726.63347	}	,
			{	12	,	6.26345	,	3.63337	,	2.84168	,	363.33665	,	726.67331	}	,
			{	12.02	,	6.26412	,	3.63357	,	2.84178	,	363.35657	,	726.71315	}	,
			{	12.04	,	6.2648	,	3.63376	,	2.84188	,	363.37649	,	726.75299	}	,
			{	12.06	,	6.26548	,	3.63396	,	2.84198	,	363.39641	,	726.79283	}	,
			{	12.08	,	6.26616	,	3.63416	,	2.84208	,	363.41633	,	726.83267	}	,
			{	12.1	,	6.26683	,	3.63436	,	2.84218	,	363.43625	,	726.87251	}	,
			{	12.12	,	6.26751	,	3.63456	,	2.84228	,	363.45618	,	726.91235	}	,
			{	12.14	,	6.26819	,	3.63476	,	2.84238	,	363.4761	,	726.95219	}	,
			{	12.16	,	6.26886	,	3.63496	,	2.84248	,	363.49602	,	726.99203	}	,
			{	12.18	,	6.26954	,	3.63516	,	2.84258	,	363.51594	,	727.03187	}	,
			{	12.2	,	6.27022	,	3.63536	,	2.84268	,	363.53586	,	727.07171	}	,
			{	12.22	,	6.2709	,	3.63556	,	2.84278	,	363.55578	,	727.11155	}	,
			{	12.24	,	6.27157	,	3.63576	,	2.84288	,	363.5757	,	727.15139	}	,
			{	12.26	,	6.27225	,	3.63596	,	2.84298	,	363.59562	,	727.19124	}	,
			{	12.28	,	6.27293	,	3.63616	,	2.84308	,	363.61554	,	727.23108	}	,
			{	12.3	,	6.27361	,	3.63635	,	2.84318	,	363.63546	,	727.27092	}	,
			{	12.32	,	6.27428	,	3.63655	,	2.84328	,	363.65538	,	727.31076	}	,
			{	12.34	,	6.27496	,	3.63675	,	2.84338	,	363.6753	,	727.3506	}	,
			{	12.36	,	6.27564	,	3.63695	,	2.84348	,	363.69522	,	727.39044	}	,
			{	12.38	,	6.27631	,	3.63715	,	2.84358	,	363.71514	,	727.43028	}	,
			{	12.4	,	6.27699	,	3.63735	,	2.84368	,	363.73506	,	727.47012	}	,
			{	12.42	,	6.27767	,	3.63755	,	2.84377	,	363.75498	,	727.50996	}	,
			{	12.44	,	6.27835	,	3.63775	,	2.84387	,	363.7749	,	727.5498	}	,
			{	12.46	,	6.27902	,	3.63795	,	2.84397	,	363.79482	,	727.58964	}	,
			{	12.48	,	6.2797	,	3.63815	,	2.84407	,	363.81474	,	727.62948	}	,
			{	12.5	,	6.28038	,	3.63835	,	2.84417	,	363.83466	,	727.66932	}	,
			{	12.52	,	6.28106	,	3.63855	,	2.84427	,	363.85458	,	727.70916	}	,
			{	12.54	,	6.28173	,	3.63875	,	2.84437	,	363.8745	,	727.749	}	,
			{	12.56	,	6.28241	,	3.63894	,	2.84447	,	363.89442	,	727.78884	}	,
			{	12.58	,	6.28309	,	3.63914	,	2.84457	,	363.91434	,	727.82869	}	,
			{	12.6	,	6.28376	,	3.63934	,	2.84467	,	363.93426	,	727.86853	}	,
			{	12.62	,	6.28444	,	3.63954	,	2.84477	,	363.95418	,	727.90837	}	,
			{	12.64	,	6.28512	,	3.63974	,	2.84487	,	363.9741	,	727.94821	}	,
			{	12.66	,	6.2858	,	3.63994	,	2.84497	,	363.99402	,	727.98805	}	,
			{	12.68	,	6.28647	,	3.64014	,	2.84507	,	364.01394	,	728.02789	}	,
			{	12.7	,	6.28715	,	3.64034	,	2.84517	,	364.03386	,	728.06773	}	,
			{	12.72	,	6.28783	,	3.64054	,	2.84527	,	364.05378	,	728.10757	}	,
			{	12.74	,	6.28851	,	3.64074	,	2.84537	,	364.07371	,	728.14741	}	,
			{	12.76	,	6.28918	,	3.64094	,	2.84547	,	364.09363	,	728.18725	}	,
			{	12.78	,	6.28986	,	3.64114	,	2.84557	,	364.11355	,	728.22709	}	,
			{	12.8	,	6.29054	,	3.64133	,	2.84567	,	364.13347	,	728.26693	}	,
			{	12.82	,	6.29122	,	3.64153	,	2.84577	,	364.15339	,	728.30677	}	,
			{	12.84	,	6.29189	,	3.64173	,	2.84587	,	364.17331	,	728.34661	}	,
			{	12.86	,	6.29257	,	3.64193	,	2.84597	,	364.19323	,	728.38645	}	,
			{	12.88	,	6.29325	,	3.64213	,	2.84607	,	364.21315	,	728.42629	}	,
			{	12.9	,	6.29392	,	3.64233	,	2.84617	,	364.23307	,	728.46614	}	,
			{	12.92	,	6.2946	,	3.64253	,	2.84626	,	364.25299	,	728.50598	}	,
			{	12.94	,	6.29528	,	3.64273	,	2.84636	,	364.27291	,	728.54582	}	,
			{	12.96	,	6.29596	,	3.64293	,	2.84646	,	364.29283	,	728.58566	}	,
			{	12.98	,	6.29663	,	3.64313	,	2.84656	,	364.31275	,	728.6255	}	,
			{	13	,	6.29731	,	3.64333	,	2.84666	,	364.33267	,	728.66534	}	,
			{	13.02	,	6.29799	,	3.64353	,	2.84676	,	364.35259	,	728.70518	}	,
			{	13.04	,	6.29867	,	3.64373	,	2.84686	,	364.37251	,	728.74502	}	,
			{	13.06	,	6.29934	,	3.64392	,	2.84696	,	364.39243	,	728.78486	}	,
			{	13.08	,	6.30002	,	3.64412	,	2.84706	,	364.41235	,	728.8247	}	,
			{	13.1	,	6.3007	,	3.64432	,	2.84716	,	364.43227	,	728.86454	}	,
			{	13.12	,	6.30137	,	3.64452	,	2.84726	,	364.45219	,	728.90438	}	,
			{	13.14	,	6.30205	,	3.64472	,	2.84736	,	364.47211	,	728.94422	}	,
			{	13.16	,	6.30273	,	3.64492	,	2.84746	,	364.49203	,	728.98406	}	,
			{	13.18	,	6.30341	,	3.64512	,	2.84756	,	364.51195	,	729.0239	}	,
			{	13.2	,	6.30408	,	3.64532	,	2.84766	,	364.53187	,	729.06375	}	,
			{	13.22	,	6.30476	,	3.64552	,	2.84776	,	364.55179	,	729.10359	}	,
			{	13.24	,	6.30544	,	3.64572	,	2.84786	,	364.57171	,	729.14343	}	,
			{	13.26	,	6.30612	,	3.64592	,	2.84796	,	364.59163	,	729.18327	}	,
			{	13.28	,	6.30679	,	3.64612	,	2.84806	,	364.61155	,	729.22311	}	,
			{	13.3	,	6.30747	,	3.64631	,	2.84816	,	364.63147	,	729.26295	}	,
			{	13.32	,	6.30815	,	3.64651	,	2.84826	,	364.65139	,	729.30279	}	,
			{	13.34	,	6.30882	,	3.64671	,	2.84836	,	364.67131	,	729.34263	}	,
			{	13.36	,	6.3095	,	3.64691	,	2.84846	,	364.69124	,	729.38247	}	,
			{	13.38	,	6.31018	,	3.64711	,	2.84856	,	364.71116	,	729.42231	}	,
			{	13.4	,	6.31086	,	3.64731	,	2.84866	,	364.73108	,	729.46215	}	,
			{	13.42	,	6.31153	,	3.64751	,	2.84875	,	364.751	,	729.50199	}	,
			{	13.44	,	6.31221	,	3.64771	,	2.84885	,	364.77092	,	729.54183	}	,
			{	13.46	,	6.31289	,	3.64791	,	2.84895	,	364.79084	,	729.58167	}	,
			{	13.48	,	6.31357	,	3.64811	,	2.84905	,	364.81076	,	729.62151	}	,
			{	13.5	,	6.31424	,	3.64831	,	2.84915	,	364.83068	,	729.66135	}	,
			{	13.52	,	6.31492	,	3.64851	,	2.84925	,	364.8506	,	729.7012	}	,
			{	13.54	,	6.3156	,	3.64871	,	2.84935	,	364.87052	,	729.74104	}	,
			{	13.56	,	6.31627	,	3.6489	,	2.84945	,	364.89044	,	729.78088	}	,
			{	13.58	,	6.31695	,	3.6491	,	2.84955	,	364.91036	,	729.82072	}	,
			{	13.6	,	6.31763	,	3.6493	,	2.84965	,	364.93028	,	729.86056	}	,
			{	13.62	,	6.31831	,	3.6495	,	2.84975	,	364.9502	,	729.9004	}	,
			{	13.64	,	6.31898	,	3.6497	,	2.84985	,	364.97012	,	729.94024	}	,
			{	13.66	,	6.31966	,	3.6499	,	2.84995	,	364.99004	,	729.98008	}	,
			{	13.68	,	6.32046	,	3.6501	,	2.8501	,	365.00996	,	730.01992	}	,
			{	13.7	,	6.32137	,	3.6503	,	2.8503	,	365.02988	,	730.05976	}	,
			{	13.72	,	6.32229	,	3.6505	,	2.8505	,	365.0498	,	730.0996	}	,
			{	13.74	,	6.32321	,	3.6507	,	2.8507	,	365.06972	,	730.13944	}	,
			{	13.76	,	6.32412	,	3.6509	,	2.8509	,	365.08964	,	730.17928	}	,
			{	13.78	,	6.32504	,	3.6511	,	2.8511	,	365.10956	,	730.21912	}	,
			{	13.8	,	6.32596	,	3.65129	,	2.85129	,	365.12948	,	730.25896	}	,
			{	13.82	,	6.32687	,	3.65149	,	2.85149	,	365.1494	,	730.2988	}	,
			{	13.84	,	6.32779	,	3.65169	,	2.85169	,	365.16932	,	730.33865	}	,
			{	13.86	,	6.32871	,	3.65189	,	2.85189	,	365.18924	,	730.37849	}	,
			{	13.88	,	6.32962	,	3.65209	,	2.85209	,	365.20916	,	730.41833	}	,
			{	13.9	,	6.33054	,	3.65229	,	2.85229	,	365.22908	,	730.45817	}	,
			{	13.92	,	6.33145	,	3.65249	,	2.85249	,	365.249	,	730.49801	}	,
			{	13.94	,	6.33237	,	3.65269	,	2.85269	,	365.26892	,	730.53785	}	,
			{	13.96	,	6.33329	,	3.65289	,	2.85289	,	365.28884	,	730.57769	}	,
			{	13.98	,	6.3342	,	3.65309	,	2.85309	,	365.30876	,	730.61753	}	,
			{	14	,	6.33512	,	3.65329	,	2.85329	,	365.32869	,	730.65737	}	,
			{	14.02	,	6.33604	,	3.65349	,	2.85349	,	365.34861	,	730.69721	}	,
			{	14.04	,	6.33695	,	3.65369	,	2.85369	,	365.36853	,	730.73705	}	,
			{	14.06	,	6.33787	,	3.65388	,	2.85388	,	365.38845	,	730.77689	}	,
			{	14.08	,	6.33878	,	3.65408	,	2.85408	,	365.40837	,	730.81673	}	,
			{	14.1	,	6.3397	,	3.65428	,	2.85428	,	365.42829	,	730.85657	}	,
			{	14.12	,	6.34062	,	3.65448	,	2.85448	,	365.44821	,	730.89641	}	,
			{	14.14	,	6.34153	,	3.65468	,	2.85468	,	365.46813	,	730.93625	}	,
			{	14.16	,	6.34245	,	3.65488	,	2.85488	,	365.48805	,	730.9761	}	,
			{	14.18	,	6.34337	,	3.65508	,	2.85508	,	365.50797	,	731.01594	}	,
			{	14.2	,	6.34428	,	3.65528	,	2.85528	,	365.52789	,	731.05578	}	,
			{	14.22	,	6.3452	,	3.65548	,	2.85548	,	365.54781	,	731.09562	}	,
			{	14.24	,	6.34612	,	3.65568	,	2.85568	,	365.56773	,	731.13546	}	,
			{	14.26	,	6.34703	,	3.65588	,	2.85588	,	365.58765	,	731.1753	}	,
			{	14.28	,	6.34795	,	3.65608	,	2.85608	,	365.60757	,	731.21514	}	,
			{	14.3	,	6.34886	,	3.65627	,	2.85627	,	365.62749	,	731.25498	}	,
			{	14.32	,	6.34978	,	3.65647	,	2.85647	,	365.64741	,	731.29482	}	,
			{	14.34	,	6.3507	,	3.65667	,	2.85667	,	365.66733	,	731.33466	}	,
			{	14.36	,	6.35161	,	3.65687	,	2.85687	,	365.68725	,	731.3745	}	,
			{	14.38	,	6.35253	,	3.65707	,	2.85707	,	365.70717	,	731.41434	}	,
			{	14.4	,	6.35345	,	3.65727	,	2.85727	,	365.72709	,	731.45418	}	,
			{	14.42	,	6.35436	,	3.65747	,	2.85747	,	365.74701	,	731.49402	}	,
			{	14.44	,	6.35528	,	3.65767	,	2.85767	,	365.76693	,	731.53386	}	,
			{	14.46	,	6.3562	,	3.65787	,	2.85787	,	365.78685	,	731.57371	}	,
			{	14.48	,	6.35711	,	3.65807	,	2.85807	,	365.80677	,	731.61355	}	,
			{	14.5	,	6.35803	,	3.65827	,	2.85827	,	365.82669	,	731.65339	}	,
			{	14.52	,	6.35894	,	3.65847	,	2.85847	,	365.84661	,	731.69323	}	,
			{	14.54	,	6.35986	,	3.65867	,	2.85867	,	365.86653	,	731.73307	}	,
			{	14.56	,	6.36078	,	3.65886	,	2.85886	,	365.88645	,	731.77291	}	,
			{	14.58	,	6.36169	,	3.65906	,	2.85906	,	365.90637	,	731.81275	}	,
			{	14.6	,	6.36261	,	3.65926	,	2.85926	,	365.92629	,	731.85259	}	,
			{	14.62	,	6.36353	,	3.65946	,	2.85946	,	365.94622	,	731.89243	}	,
			{	14.64	,	6.36444	,	3.65966	,	2.85966	,	365.96614	,	731.93227	}	,
			{	14.66	,	6.36536	,	3.65986	,	2.85986	,	365.98606	,	731.97211	}	,
			{	14.68	,	6.36627	,	3.66006	,	2.86006	,	366.00598	,	732.01195	}	,
			{	14.7	,	6.36719	,	3.66026	,	2.86026	,	366.0259	,	732.05179	}	,
			{	14.72	,	6.36811	,	3.66046	,	2.86046	,	366.04582	,	732.09163	}	,
			{	14.74	,	6.36902	,	3.66066	,	2.86066	,	366.06574	,	732.13147	}	,
			{	14.76	,	6.36994	,	3.66086	,	2.86086	,	366.08566	,	732.17131	}	,
			{	14.78	,	6.37086	,	3.66106	,	2.86106	,	366.10558	,	732.21116	}	,
			{	14.8	,	6.37177	,	3.66125	,	2.86125	,	366.1255	,	732.251	}	,
			{	14.82	,	6.37269	,	3.66145	,	2.86145	,	366.14542	,	732.29084	}	,
			{	14.84	,	6.37361	,	3.66165	,	2.86165	,	366.16534	,	732.33068	}	,
			{	14.86	,	6.37452	,	3.66185	,	2.86185	,	366.18526	,	732.37052	}	,
			{	14.88	,	6.37544	,	3.66205	,	2.86205	,	366.20518	,	732.41036	}	,
			{	14.9	,	6.37635	,	3.66225	,	2.86225	,	366.2251	,	732.4502	}	,
			{	14.92	,	6.37727	,	3.66245	,	2.86245	,	366.24502	,	732.49004	}	,
			{	14.94	,	6.37819	,	3.66265	,	2.86265	,	366.26494	,	732.52988	}	,
			{	14.96	,	6.3791	,	3.66285	,	2.86285	,	366.28486	,	732.56972	}	,
			{	14.98	,	6.38002	,	3.66305	,	2.86305	,	366.30478	,	732.60956	}	,
			{	15	,	6.38094	,	3.66325	,	2.86325	,	366.3247	,	732.6494	}	,
			{	15.02	,	6.38185	,	3.66345	,	2.86345	,	366.34462	,	732.68924	}	,
			{	15.04	,	6.38277	,	3.66365	,	2.86365	,	366.36454	,	732.72908	}	,
			{	15.06	,	6.38369	,	3.66384	,	2.86384	,	366.38446	,	732.76892	}	,
			{	15.08	,	6.3846	,	3.66404	,	2.86404	,	366.40438	,	732.80876	}	,
			{	15.1	,	6.38552	,	3.66424	,	2.86424	,	366.4243	,	732.84861	}	,
			{	15.12	,	6.38643	,	3.66444	,	2.86444	,	366.44422	,	732.88845	}	,
			{	15.14	,	6.38735	,	3.66464	,	2.86464	,	366.46414	,	732.92829	}	,
			{	15.16	,	6.38827	,	3.66484	,	2.86484	,	366.48406	,	732.96813	}	,
			{	15.18	,	6.38918	,	3.66504	,	2.86504	,	366.50398	,	733.00797	}	,
			{	15.2	,	6.3901	,	3.66524	,	2.86524	,	366.5239	,	733.04781	}	,
			{	15.22	,	6.39102	,	3.66544	,	2.86544	,	366.54382	,	733.08765	}	,
			{	15.24	,	6.39193	,	3.66564	,	2.86564	,	366.56375	,	733.12749	}	,
			{	15.26	,	6.39285	,	3.66584	,	2.86584	,	366.58367	,	733.16733	}	,
			{	15.28	,	6.39376	,	3.66604	,	2.86604	,	366.60359	,	733.20717	}	,
			{	15.3	,	6.39468	,	3.66624	,	2.86624	,	366.62351	,	733.24701	}	,
			{	15.32	,	6.3956	,	3.66643	,	2.86643	,	366.64343	,	733.28685	}	,
			{	15.34	,	6.39651	,	3.66663	,	2.86663	,	366.66335	,	733.32669	}	,
			{	15.36	,	6.39743	,	3.66683	,	2.86683	,	366.68327	,	733.36653	}	,
			{	15.38	,	6.39835	,	3.66703	,	2.86703	,	366.70319	,	733.40637	}	,
			{	15.4	,	6.39926	,	3.66723	,	2.86723	,	366.72311	,	733.44622	}	,
			{	15.42	,	6.40018	,	3.66743	,	2.86743	,	366.74303	,	733.48606	}	,
			{	15.44	,	6.4011	,	3.66763	,	2.86763	,	366.76295	,	733.5259	}	,
			{	15.46	,	6.40201	,	3.66783	,	2.86783	,	366.78287	,	733.56574	}	,
			{	15.48	,	6.40293	,	3.66803	,	2.86803	,	366.80279	,	733.60558	}	,
			{	15.5	,	6.40384	,	3.66823	,	2.86823	,	366.82271	,	733.64542	}	,
			{	15.52	,	6.40476	,	3.66843	,	2.86843	,	366.84263	,	733.68526	}	,
			{	15.54	,	6.40568	,	3.66863	,	2.86863	,	366.86255	,	733.7251	}	,
			{	15.56	,	6.40659	,	3.66882	,	2.86882	,	366.88247	,	733.76494	}	,
			{	15.58	,	6.40751	,	3.66902	,	2.86902	,	366.90239	,	733.80478	}	,
			{	15.6	,	6.40843	,	3.66922	,	2.86922	,	366.92231	,	733.84462	}	,
			{	15.62	,	6.40934	,	3.66942	,	2.86942	,	366.94223	,	733.88446	}	,
			{	15.64	,	6.41026	,	3.66962	,	2.86962	,	366.96215	,	733.9243	}	,
			{	15.66	,	6.41118	,	3.66982	,	2.86982	,	366.98207	,	733.96414	}	,
			{	15.68	,	6.41209	,	3.67002	,	2.87002	,	367.00199	,	734.00398	}	,
			{	15.7	,	6.41301	,	3.67022	,	2.87022	,	367.02191	,	734.04382	}	,
			{	15.72	,	6.41392	,	3.67042	,	2.87042	,	367.04183	,	734.08367	}	,
			{	15.74	,	6.41484	,	3.67062	,	2.87062	,	367.06175	,	734.12351	}	,
			{	15.76	,	6.41576	,	3.67082	,	2.87082	,	367.08167	,	734.16335	}	,
			{	15.78	,	6.41667	,	3.67102	,	2.87102	,	367.10159	,	734.20319	}	,
			{	15.8	,	6.41759	,	3.67122	,	2.87122	,	367.12151	,	734.24303	}	,
			{	15.82	,	6.41851	,	3.67141	,	2.87141	,	367.14143	,	734.28287	}	,
			{	15.84	,	6.41942	,	3.67161	,	2.87161	,	367.16135	,	734.32271	}	,
			{	15.86	,	6.42034	,	3.67181	,	2.87181	,	367.18127	,	734.36255	}	,
			{	15.88	,	6.42125	,	3.67201	,	2.87201	,	367.2012	,	734.40239	}	,
			{	15.9	,	6.42217	,	3.67221	,	2.87221	,	367.22112	,	734.44223	}	,
			{	15.92	,	6.42309	,	3.67241	,	2.87241	,	367.24104	,	734.48207	}	,
			{	15.94	,	6.424	,	3.67261	,	2.87261	,	367.26096	,	734.52191	}	,
			{	15.96	,	6.42492	,	3.67281	,	2.87281	,	367.28088	,	734.56175	}	,
			{	15.98	,	6.42584	,	3.67301	,	2.87301	,	367.3008	,	734.60159	}	,
			{	16	,	6.42675	,	3.67321	,	2.87321	,	367.32072	,	734.64143	}	,
			{	16.02	,	6.42767	,	3.67341	,	2.87341	,	367.34064	,	734.68127	}	,
			{	16.04	,	6.42859	,	3.67361	,	2.87361	,	367.36056	,	734.72112	}	,
			{	16.06	,	6.4295	,	3.6738	,	2.8738	,	367.38048	,	734.76096	}	,
			{	16.08	,	6.43042	,	3.674	,	2.874	,	367.4004	,	734.8008	}	,
			{	16.1	,	6.43133	,	3.6742	,	2.8742	,	367.42032	,	734.84064	}	,
			{	16.12	,	6.43225	,	3.6744	,	2.8744	,	367.44024	,	734.88048	}	,
			{	16.14	,	6.43317	,	3.6746	,	2.8746	,	367.46016	,	734.92032	}	,
			{	16.16	,	6.43408	,	3.6748	,	2.8748	,	367.48008	,	734.96016	}	,
			{	16.18	,	6.435	,	3.675	,	2.875	,	367.5	,	735	}	,
			{	16.2	,	6.43592	,	3.6752	,	2.8752	,	367.51992	,	735.03984	}	,
			{	16.22	,	6.43683	,	3.6754	,	2.8754	,	367.53984	,	735.07968	}	,
			{	16.24	,	6.43775	,	3.6756	,	2.8756	,	367.55976	,	735.11952	}	,
			{	16.26	,	6.43867	,	3.6758	,	2.8758	,	367.57968	,	735.15936	}	,
			{	16.28	,	6.43958	,	3.676	,	2.876	,	367.5996	,	735.1992	}	,
			{	16.3	,	6.4405	,	3.6762	,	2.8762	,	367.61952	,	735.23904	}	,
			{	16.32	,	6.44141	,	3.67639	,	2.87639	,	367.63944	,	735.27888	}	,
			{	16.34	,	6.44233	,	3.67659	,	2.87659	,	367.65936	,	735.31873	}	,
			{	16.36	,	6.44325	,	3.67679	,	2.87679	,	367.67928	,	735.35857	}	,
			{	16.38	,	6.44416	,	3.67699	,	2.87699	,	367.6992	,	735.39841	}	,
			{	16.4	,	6.44508	,	3.67719	,	2.87719	,	367.71912	,	735.43825	}	,
			{	16.42	,	6.446	,	3.67739	,	2.87739	,	367.73904	,	735.47809	}	,
			{	16.44	,	6.44691	,	3.67759	,	2.87759	,	367.75896	,	735.51793	}	,
			{	16.46	,	6.44783	,	3.67779	,	2.87779	,	367.77888	,	735.55777	}	,
			{	16.48	,	6.44875	,	3.67799	,	2.87799	,	367.7988	,	735.59761	}	,
			{	16.5	,	6.44966	,	3.67819	,	2.87819	,	367.81873	,	735.63745	}	,
			{	16.52	,	6.45058	,	3.67839	,	2.87839	,	367.83865	,	735.67729	}	,
			{	16.54	,	6.45149	,	3.67859	,	2.87859	,	367.85857	,	735.71713	}	,
			{	16.56	,	6.45241	,	3.67878	,	2.87878	,	367.87849	,	735.75697	}	,
			{	16.58	,	6.45333	,	3.67898	,	2.87898	,	367.89841	,	735.79681	}	,
			{	16.6	,	6.45424	,	3.67918	,	2.87918	,	367.91833	,	735.83665	}	,
			{	16.62	,	6.45516	,	3.67938	,	2.87938	,	367.93825	,	735.87649	}	,
			{	16.64	,	6.45608	,	3.67958	,	2.87958	,	367.95817	,	735.91633	}	,
			{	16.66	,	6.45699	,	3.67978	,	2.87978	,	367.97809	,	735.95618	}	,
			{	16.68	,	6.45791	,	3.67998	,	2.87998	,	367.99801	,	735.99602	}	,
			{	16.7	,	6.45882	,	3.68018	,	2.88018	,	368.01793	,	736.03586	}	,
			{	16.72	,	6.45974	,	3.68038	,	2.88038	,	368.03785	,	736.0757	}	,
			{	16.74	,	6.46066	,	3.68058	,	2.88058	,	368.05777	,	736.11554	}	,
			{	16.76	,	6.46157	,	3.68078	,	2.88078	,	368.07769	,	736.15538	}	,
			{	16.78	,	6.46249	,	3.68098	,	2.88098	,	368.09761	,	736.19522	}	,
			{	16.8	,	6.46341	,	3.68118	,	2.88118	,	368.11753	,	736.23506	}	,
			{	16.82	,	6.46432	,	3.68137	,	2.88137	,	368.13745	,	736.2749	}	,
			{	16.84	,	6.46524	,	3.68157	,	2.88157	,	368.15737	,	736.31474	}	,
			{	16.86	,	6.46616	,	3.68177	,	2.88177	,	368.17729	,	736.35458	}	,
			{	16.88	,	6.46707	,	3.68197	,	2.88197	,	368.19721	,	736.39442	}	,
			{	16.9	,	6.46799	,	3.68217	,	2.88217	,	368.21713	,	736.43426	}	,
			{	16.92	,	6.4689	,	3.68237	,	2.88237	,	368.23705	,	736.4741	}	,
			{	16.94	,	6.46982	,	3.68257	,	2.88257	,	368.25697	,	736.51394	}	,
			{	16.96	,	6.47074	,	3.68277	,	2.88277	,	368.27689	,	736.55378	}	,
			{	16.98	,	6.47165	,	3.68297	,	2.88297	,	368.29681	,	736.59363	}	,
			{	17	,	6.47257	,	3.68317	,	2.88317	,	368.31673	,	736.63347	}	,
			{	17.02	,	6.47349	,	3.68337	,	2.88337	,	368.33665	,	736.67331	}	,
			{	17.04	,	6.4744	,	3.68357	,	2.88357	,	368.35657	,	736.71315	}	,
			{	17.06	,	6.47532	,	3.68376	,	2.88376	,	368.37649	,	736.75299	}	,
			{	17.08	,	6.47624	,	3.68396	,	2.88396	,	368.39641	,	736.79283	}	,
			{	17.1	,	6.47715	,	3.68416	,	2.88416	,	368.41633	,	736.83267	}	,
			{	17.12	,	6.47807	,	3.68436	,	2.88436	,	368.43625	,	736.87251	}	,
			{	17.14	,	6.47898	,	3.68456	,	2.88456	,	368.45618	,	736.91235	}	,
			{	17.16	,	6.4799	,	3.68476	,	2.88476	,	368.4761	,	736.95219	}	,
			{	17.18	,	6.48082	,	3.68496	,	2.88496	,	368.49602	,	736.99203	}	,
			{	17.2	,	6.48173	,	3.68516	,	2.88516	,	368.51594	,	737.03187	}	,
			{	17.22	,	6.48265	,	3.68536	,	2.88536	,	368.53586	,	737.07171	}	,
			{	17.24	,	6.48357	,	3.68556	,	2.88556	,	368.55578	,	737.11155	}	,
			{	17.26	,	6.48448	,	3.68576	,	2.88576	,	368.5757	,	737.15139	}	,
			{	17.28	,	6.4854	,	3.68596	,	2.88596	,	368.59562	,	737.19124	}	,
			{	17.3	,	6.48631	,	3.68616	,	2.88616	,	368.61554	,	737.23108	}	,
			{	17.32	,	6.48723	,	3.68635	,	2.88635	,	368.63546	,	737.27092	}	,
			{	17.34	,	6.48815	,	3.68655	,	2.88655	,	368.65538	,	737.31076	}	,
			{	17.36	,	6.48906	,	3.68675	,	2.88675	,	368.6753	,	737.3506	}	,
			{	17.38	,	6.48998	,	3.68695	,	2.88695	,	368.69522	,	737.39044	}	,
			{	17.4	,	6.4909	,	3.68715	,	2.88715	,	368.71514	,	737.43028	}	,
			{	17.42	,	6.49181	,	3.68735	,	2.88735	,	368.73506	,	737.47012	}	,
			{	17.44	,	6.49273	,	3.68755	,	2.88755	,	368.75498	,	737.50996	}	,
			{	17.46	,	6.49365	,	3.68775	,	2.88775	,	368.7749	,	737.5498	}	,
			{	17.48	,	6.49456	,	3.68795	,	2.88795	,	368.79482	,	737.58964	}	,
			{	17.5	,	6.49548	,	3.68815	,	2.88815	,	368.81474	,	737.62948	}	,
			{	17.52	,	6.49639	,	3.68835	,	2.88835	,	368.83466	,	737.66932	}	,
			{	17.54	,	6.49731	,	3.68855	,	2.88855	,	368.85458	,	737.70916	}	,
			{	17.56	,	6.49823	,	3.68875	,	2.88875	,	368.8745	,	737.749	}	,
			{	17.58	,	6.49914	,	3.68894	,	2.88894	,	368.89442	,	737.78884	}	,
			{	17.6	,	6.50006	,	3.68914	,	2.88914	,	368.91434	,	737.82869	}	,
			{	17.62	,	6.50098	,	3.68934	,	2.88934	,	368.93426	,	737.86853	}	,
			{	17.64	,	6.50189	,	3.68954	,	2.88954	,	368.95418	,	737.90837	}	,
			{	17.66	,	6.50281	,	3.68974	,	2.88974	,	368.9741	,	737.94821	}	,
			{	17.68	,	6.50373	,	3.68994	,	2.88994	,	368.99402	,	737.98805	}	,
			{	17.7	,	6.50464	,	3.69014	,	2.89014	,	369.01394	,	738.02789	}	,
			{	17.72	,	6.50556	,	3.69034	,	2.89034	,	369.03386	,	738.06773	}	,
			{	17.74	,	6.50647	,	3.69054	,	2.89054	,	369.05378	,	738.10757	}	,
			{	17.76	,	6.50739	,	3.69074	,	2.89074	,	369.07371	,	738.14741	}	,
			{	17.78	,	6.50831	,	3.69094	,	2.89094	,	369.09363	,	738.18725	}	,
			{	17.8	,	6.50922	,	3.69114	,	2.89114	,	369.11355	,	738.22709	}	,
			{	17.82	,	6.51014	,	3.69133	,	2.89133	,	369.13347	,	738.26693	}	,
			{	17.84	,	6.51106	,	3.69153	,	2.89153	,	369.15339	,	738.30677	}	,
			{	17.86	,	6.51197	,	3.69173	,	2.89173	,	369.17331	,	738.34661	}	,
			{	17.88	,	6.51289	,	3.69193	,	2.89193	,	369.19323	,	738.38645	}	,
			{	17.9	,	6.5138	,	3.69213	,	2.89213	,	369.21315	,	738.42629	}	,
			{	17.92	,	6.51472	,	3.69233	,	2.89233	,	369.23307	,	738.46614	}	,
			{	17.94	,	6.51564	,	3.69253	,	2.89253	,	369.25299	,	738.50598	}	,
			{	17.96	,	6.51655	,	3.69273	,	2.89273	,	369.27291	,	738.54582	}	,
			{	17.98	,	6.51747	,	3.69293	,	2.89293	,	369.29283	,	738.58566	}	,
			{	18	,	6.51839	,	3.69313	,	2.89313	,	369.31275	,	738.6255	}	,
			{	18.02	,	6.5193	,	3.69333	,	2.89333	,	369.33267	,	738.66534	}	,
			{	18.04	,	6.52022	,	3.69353	,	2.89353	,	369.35259	,	738.70518	}	,
			{	18.06	,	6.52114	,	3.69373	,	2.89373	,	369.37251	,	738.74502	}	,
			{	18.08	,	6.52205	,	3.69392	,	2.89392	,	369.39243	,	738.78486	}	,
			{	18.1	,	6.52297	,	3.69412	,	2.89412	,	369.41235	,	738.8247	}	,
			{	18.12	,	6.52388	,	3.69432	,	2.89432	,	369.43227	,	738.86454	}	,
			{	18.14	,	6.5248	,	3.69452	,	2.89452	,	369.45219	,	738.90438	}	,
			{	18.16	,	6.52572	,	3.69472	,	2.89472	,	369.47211	,	738.94422	}	,
			{	18.18	,	6.52663	,	3.69492	,	2.89492	,	369.49203	,	738.98406	}	,
			{	18.2	,	6.52755	,	3.69512	,	2.89512	,	369.51195	,	739.0239	}	,
			{	18.22	,	6.52847	,	3.69532	,	2.89532	,	369.53187	,	739.06375	}	,
			{	18.24	,	6.52938	,	3.69552	,	2.89552	,	369.55179	,	739.10359	}	,
			{	18.26	,	6.5303	,	3.69572	,	2.89572	,	369.57171	,	739.14343	}	,
			{	18.28	,	6.53122	,	3.69592	,	2.89592	,	369.59163	,	739.18327	}	,
			{	18.3	,	6.53213	,	3.69612	,	2.89612	,	369.61155	,	739.22311	}	,
			{	18.32	,	6.53305	,	3.69631	,	2.89631	,	369.63147	,	739.26295	}	,
			{	18.34	,	6.53396	,	3.69651	,	2.89651	,	369.65139	,	739.30279	}	,
			{	18.36	,	6.53488	,	3.69671	,	2.89671	,	369.67131	,	739.34263	}	,
			{	18.38	,	6.5358	,	3.69691	,	2.89691	,	369.69124	,	739.38247	}	,
			{	18.4	,	6.53671	,	3.69711	,	2.89711	,	369.71116	,	739.42231	}	,
			{	18.42	,	6.53763	,	3.69731	,	2.89731	,	369.73108	,	739.46215	}	,
			{	18.44	,	6.53855	,	3.69751	,	2.89751	,	369.751	,	739.50199	}	,
			{	18.46	,	6.53946	,	3.69771	,	2.89771	,	369.77092	,	739.54183	}	,
			{	18.48	,	6.54038	,	3.69791	,	2.89791	,	369.79084	,	739.58167	}	,
			{	18.5	,	6.54129	,	3.69811	,	2.89811	,	369.81076	,	739.62151	}	,
			{	18.52	,	6.54221	,	3.69831	,	2.89831	,	369.83068	,	739.66135	}	,
			{	18.54	,	6.54313	,	3.69851	,	2.89851	,	369.8506	,	739.7012	}	,
			{	18.56	,	6.54404	,	3.69871	,	2.89871	,	369.87052	,	739.74104	}	,
			{	18.58	,	6.54496	,	3.6989	,	2.8989	,	369.89044	,	739.78088	}	,
			{	18.6	,	6.54588	,	3.6991	,	2.8991	,	369.91036	,	739.82072	}	,
			{	18.62	,	6.54679	,	3.6993	,	2.8993	,	369.93028	,	739.86056	}	,
			{	18.64	,	6.54771	,	3.6995	,	2.8995	,	369.9502	,	739.9004	}	,
			{	18.66	,	6.54863	,	3.6997	,	2.8997	,	369.97012	,	739.94024	}	,
			{	18.68	,	6.54954	,	3.6999	,	2.8999	,	369.99004	,	739.98008	}	,
			{	18.7	,	6.55034	,	3.70013	,	2.90007	,	370.0135	,	740.02699	}	,
			{	18.72	,	6.55101	,	3.7004	,	2.9002	,	370.04049	,	740.08097	}	,
			{	18.74	,	6.55169	,	3.70067	,	2.90034	,	370.06748	,	740.13495	}	,
			{	18.76	,	6.55236	,	3.70094	,	2.90047	,	370.09447	,	740.18893	}	,
			{	18.78	,	6.55304	,	3.70121	,	2.90061	,	370.12146	,	740.24291	}	,
			{	18.8	,	6.55371	,	3.70148	,	2.90074	,	370.14845	,	740.2969	}	,
			{	18.82	,	6.55439	,	3.70175	,	2.90088	,	370.17544	,	740.35088	}	,
			{	18.84	,	6.55506	,	3.70202	,	2.90101	,	370.20243	,	740.40486	}	,
			{	18.86	,	6.55574	,	3.70229	,	2.90115	,	370.22942	,	740.45884	}	,
			{	18.88	,	6.55641	,	3.70256	,	2.90128	,	370.25641	,	740.51282	}	,
			{	18.9	,	6.55709	,	3.70283	,	2.90142	,	370.2834	,	740.5668	}	,
			{	18.92	,	6.55776	,	3.7031	,	2.90155	,	370.31039	,	740.62078	}	,
			{	18.94	,	6.55843	,	3.70337	,	2.90169	,	370.33738	,	740.67476	}	,
			{	18.96	,	6.55911	,	3.70364	,	2.90182	,	370.36437	,	740.72874	}	,
			{	18.98	,	6.55978	,	3.70391	,	2.90196	,	370.39136	,	740.78273	}	,
			{	19	,	6.56046	,	3.70418	,	2.90209	,	370.41835	,	740.83671	}	,
			{	19.02	,	6.56113	,	3.70445	,	2.90223	,	370.44534	,	740.89069	}	,
			{	19.04	,	6.56181	,	3.70472	,	2.90236	,	370.47233	,	740.94467	}	,
			{	19.06	,	6.56248	,	3.70499	,	2.9025	,	370.49933	,	740.99865	}	,
			{	19.08	,	6.56316	,	3.70526	,	2.90263	,	370.52632	,	741.05263	}	,
			{	19.1	,	6.56383	,	3.70553	,	2.90277	,	370.55331	,	741.10661	}	,
			{	19.12	,	6.56451	,	3.7058	,	2.9029	,	370.5803	,	741.16059	}	,
			{	19.14	,	6.56518	,	3.70607	,	2.90304	,	370.60729	,	741.21457	}	,
			{	19.16	,	6.56586	,	3.70634	,	2.90317	,	370.63428	,	741.26856	}	,
			{	19.18	,	6.56653	,	3.70661	,	2.90331	,	370.66127	,	741.32254	}	,
			{	19.2	,	6.56721	,	3.70688	,	2.90344	,	370.68826	,	741.37652	}	,
			{	19.22	,	6.56788	,	3.70715	,	2.90358	,	370.71525	,	741.4305	}	,
			{	19.24	,	6.56856	,	3.70742	,	2.90371	,	370.74224	,	741.48448	}	,
			{	19.26	,	6.56923	,	3.70769	,	2.90385	,	370.76923	,	741.53846	}	,
			{	19.28	,	6.56991	,	3.70796	,	2.90398	,	370.79622	,	741.59244	}	,
			{	19.3	,	6.57058	,	3.70823	,	2.90412	,	370.82321	,	741.64642	}	,
			{	19.32	,	6.57126	,	3.7085	,	2.90425	,	370.8502	,	741.7004	}	,
			{	19.34	,	6.57193	,	3.70877	,	2.90439	,	370.87719	,	741.75439	}	,
			{	19.36	,	6.5726	,	3.70904	,	2.90452	,	370.90418	,	741.80837	}	,
			{	19.38	,	6.57328	,	3.70931	,	2.90466	,	370.93117	,	741.86235	}	,
			{	19.4	,	6.57395	,	3.70958	,	2.90479	,	370.95816	,	741.91633	}	,
			{	19.42	,	6.57463	,	3.70985	,	2.90493	,	370.98516	,	741.97031	}	,
			{	19.44	,	6.5753	,	3.71012	,	2.90506	,	371.01215	,	742.02429	}	,
			{	19.46	,	6.57598	,	3.71039	,	2.9052	,	371.03914	,	742.07827	}	,
			{	19.48	,	6.57665	,	3.71066	,	2.90533	,	371.06613	,	742.13225	}	,
			{	19.5	,	6.57733	,	3.71093	,	2.90547	,	371.09312	,	742.18623	}	,
			{	19.52	,	6.578	,	3.7112	,	2.9056	,	371.12011	,	742.24022	}	,
			{	19.54	,	6.57868	,	3.71147	,	2.90574	,	371.1471	,	742.2942	}	,
			{	19.56	,	6.57935	,	3.71174	,	2.90587	,	371.17409	,	742.34818	}	,
			{	19.58	,	6.58003	,	3.71201	,	2.90601	,	371.20108	,	742.40216	}	,
			{	19.6	,	6.5807	,	3.71228	,	2.90614	,	371.22807	,	742.45614	}	,
			{	19.62	,	6.58138	,	3.71255	,	2.90628	,	371.25506	,	742.51012	}	,
			{	19.64	,	6.58205	,	3.71282	,	2.90641	,	371.28205	,	742.5641	}	,
			{	19.66	,	6.58273	,	3.71309	,	2.90655	,	371.30904	,	742.61808	}	,
			{	19.68	,	6.5834	,	3.71336	,	2.90668	,	371.33603	,	742.67206	}	,
			{	19.7	,	6.58408	,	3.71363	,	2.90682	,	371.36302	,	742.72605	}	,
			{	19.72	,	6.58475	,	3.7139	,	2.90695	,	371.39001	,	742.78003	}	,
			{	19.74	,	6.58543	,	3.71417	,	2.90709	,	371.417	,	742.83401	}	,
			{	19.76	,	6.5861	,	3.71444	,	2.90722	,	371.44399	,	742.88799	}	,
			{	19.78	,	6.58677	,	3.71471	,	2.90735	,	371.47099	,	742.94197	}	,
			{	19.8	,	6.58745	,	3.71498	,	2.90749	,	371.49798	,	742.99595	}	,
			{	19.82	,	6.58812	,	3.71525	,	2.90762	,	371.52497	,	743.04993	}	,
			{	19.84	,	6.5888	,	3.71552	,	2.90776	,	371.55196	,	743.10391	}	,
			{	19.86	,	6.58947	,	3.71579	,	2.90789	,	371.57895	,	743.15789	}	,
			{	19.88	,	6.59015	,	3.71606	,	2.90803	,	371.60594	,	743.21188	}	,
			{	19.9	,	6.59082	,	3.71633	,	2.90816	,	371.63293	,	743.26586	}	,
			{	19.92	,	6.5915	,	3.7166	,	2.9083	,	371.65992	,	743.31984	}	,
			{	19.94	,	6.59217	,	3.71687	,	2.90843	,	371.68691	,	743.37382	}	,
			{	19.96	,	6.59285	,	3.71714	,	2.90857	,	371.7139	,	743.4278	}	,
			{	19.98	,	6.59352	,	3.71741	,	2.9087	,	371.74089	,	743.48178	}	,
			{	20	,	6.5942	,	3.71768	,	2.90884	,	371.76788	,	743.53576	}	,
			{	20.02	,	6.59487	,	3.71795	,	2.90897	,	371.79487	,	743.58974	}	,
			{	20.04	,	6.59555	,	3.71822	,	2.90911	,	371.82186	,	743.64372	}	,
			{	20.06	,	6.59622	,	3.71849	,	2.90924	,	371.84885	,	743.69771	}	,
			{	20.08	,	6.5969	,	3.71876	,	2.90938	,	371.87584	,	743.75169	}	,
			{	20.1	,	6.59757	,	3.71903	,	2.90951	,	371.90283	,	743.80567	}	,
			{	20.12	,	6.59825	,	3.7193	,	2.90965	,	371.92982	,	743.85965	}	,
			{	20.14	,	6.59892	,	3.71957	,	2.90978	,	371.95682	,	743.91363	}	,
			{	20.16	,	6.5996	,	3.71984	,	2.90992	,	371.98381	,	743.96761	}	,
			{	20.18	,	6.60027	,	3.72011	,	2.91005	,	372.0108	,	744.02159	}	,
			{	20.2	,	6.60094	,	3.72038	,	2.91019	,	372.03779	,	744.07557	}	,
			{	20.22	,	6.60162	,	3.72065	,	2.91032	,	372.06478	,	744.12955	}	,
			{	20.24	,	6.60229	,	3.72092	,	2.91046	,	372.09177	,	744.18354	}	,
			{	20.26	,	6.60297	,	3.72119	,	2.91059	,	372.11876	,	744.23752	}	,
			{	20.28	,	6.60364	,	3.72146	,	2.91073	,	372.14575	,	744.2915	}	,
			{	20.3	,	6.60432	,	3.72173	,	2.91086	,	372.17274	,	744.34548	}	,
			{	20.32	,	6.60499	,	3.722	,	2.911	,	372.19973	,	744.39946	}	,
			{	20.34	,	6.60567	,	3.72227	,	2.91113	,	372.22672	,	744.45344	}	,
			{	20.36	,	6.60634	,	3.72254	,	2.91127	,	372.25371	,	744.50742	}	,
			{	20.38	,	6.60702	,	3.72281	,	2.9114	,	372.2807	,	744.5614	}	,
			{	20.4	,	6.60769	,	3.72308	,	2.91154	,	372.30769	,	744.61538	}	,
			{	20.42	,	6.60837	,	3.72335	,	2.91167	,	372.33468	,	744.66937	}	,
			{	20.44	,	6.60904	,	3.72362	,	2.91181	,	372.36167	,	744.72335	}	,
			{	20.46	,	6.60972	,	3.72389	,	2.91194	,	372.38866	,	744.77733	}	,
			{	20.48	,	6.61039	,	3.72416	,	2.91208	,	372.41565	,	744.83131	}	,
			{	20.5	,	6.61107	,	3.72443	,	2.91221	,	372.44265	,	744.88529	}	,
			{	20.52	,	6.61174	,	3.7247	,	2.91235	,	372.46964	,	744.93927	}	,
			{	20.54	,	6.61242	,	3.72497	,	2.91248	,	372.49663	,	744.99325	}	,
			{	20.56	,	6.61309	,	3.72524	,	2.91262	,	372.52362	,	745.04723	}	,
			{	20.58	,	6.61377	,	3.72551	,	2.91275	,	372.55061	,	745.10121	}	,
			{	20.6	,	6.61444	,	3.72578	,	2.91289	,	372.5776	,	745.1552	}	,
			{	20.62	,	6.61511	,	3.72605	,	2.91302	,	372.60459	,	745.20918	}	,
			{	20.64	,	6.61579	,	3.72632	,	2.91316	,	372.63158	,	745.26316	}	,
			{	20.66	,	6.61646	,	3.72659	,	2.91329	,	372.65857	,	745.31714	}	,
			{	20.68	,	6.61714	,	3.72686	,	2.91343	,	372.68556	,	745.37112	}	,
			{	20.7	,	6.61781	,	3.72713	,	2.91356	,	372.71255	,	745.4251	}	,
			{	20.72	,	6.61849	,	3.7274	,	2.9137	,	372.73954	,	745.47908	}	,
			{	20.74	,	6.61916	,	3.72767	,	2.91383	,	372.76653	,	745.53306	}	,
			{	20.76	,	6.61984	,	3.72794	,	2.91397	,	372.79352	,	745.58704	}	,
			{	20.78	,	6.62051	,	3.72821	,	2.9141	,	372.82051	,	745.64103	}	,
			{	20.8	,	6.62119	,	3.72848	,	2.91424	,	372.8475	,	745.69501	}	,
			{	20.82	,	6.62186	,	3.72874	,	2.91437	,	372.87449	,	745.74899	}	,
			{	20.84	,	6.62254	,	3.72901	,	2.91451	,	372.90148	,	745.80297	}	,
			{	20.86	,	6.62321	,	3.72928	,	2.91464	,	372.92848	,	745.85695	}	,
			{	20.88	,	6.62389	,	3.72955	,	2.91478	,	372.95547	,	745.91093	}	,
			{	20.9	,	6.62456	,	3.72982	,	2.91491	,	372.98246	,	745.96491	}	,
			{	20.92	,	6.62524	,	3.73009	,	2.91505	,	373.00945	,	746.01889	}	,
			{	20.94	,	6.62591	,	3.73036	,	2.91518	,	373.03644	,	746.07287	}	,
			{	20.96	,	6.62659	,	3.73063	,	2.91532	,	373.06343	,	746.12686	}	,
			{	20.98	,	6.62726	,	3.7309	,	2.91545	,	373.09042	,	746.18084	}	,
			{	21	,	6.62794	,	3.73117	,	2.91559	,	373.11741	,	746.23482	}	,
			{	21.02	,	6.62861	,	3.73144	,	2.91572	,	373.1444	,	746.2888	}	,
			{	21.04	,	6.62928	,	3.73171	,	2.91586	,	373.17139	,	746.34278	}	,
			{	21.06	,	6.62996	,	3.73198	,	2.91599	,	373.19838	,	746.39676	}	,
			{	21.08	,	6.63063	,	3.73225	,	2.91613	,	373.22537	,	746.45074	}	,
			{	21.1	,	6.63131	,	3.73252	,	2.91626	,	373.25236	,	746.50472	}	,
			{	21.12	,	6.63198	,	3.73279	,	2.9164	,	373.27935	,	746.5587	}	,
			{	21.14	,	6.63266	,	3.73306	,	2.91653	,	373.30634	,	746.61269	}	,
			{	21.16	,	6.63333	,	3.73333	,	2.91667	,	373.33333	,	746.66667	}	,
			{	21.18	,	6.63401	,	3.7336	,	2.9168	,	373.36032	,	746.72065	}	,
			{	21.2	,	6.63468	,	3.73387	,	2.91694	,	373.38731	,	746.77463	}	,
			{	21.22	,	6.63536	,	3.73414	,	2.91707	,	373.4143	,	746.82861	}	,
			{	21.24	,	6.63603	,	3.73441	,	2.91721	,	373.4413	,	746.88259	}	,
			{	21.26	,	6.63671	,	3.73468	,	2.91734	,	373.46829	,	746.93657	}	,
			{	21.28	,	6.63738	,	3.73495	,	2.91748	,	373.49528	,	746.99055	}	,
			{	21.3	,	6.63806	,	3.73522	,	2.91761	,	373.52227	,	747.04453	}	,
			{	21.32	,	6.63873	,	3.73549	,	2.91775	,	373.54926	,	747.09852	}	,
			{	21.34	,	6.63941	,	3.73576	,	2.91788	,	373.57625	,	747.1525	}	,
			{	21.36	,	6.64008	,	3.73603	,	2.91802	,	373.60324	,	747.20648	}	,
			{	21.38	,	6.64076	,	3.7363	,	2.91815	,	373.63023	,	747.26046	}	,
			{	21.4	,	6.64143	,	3.73657	,	2.91829	,	373.65722	,	747.31444	}	,
			{	21.42	,	6.64211	,	3.73684	,	2.91842	,	373.68421	,	747.36842	}	,
			{	21.44	,	6.64278	,	3.73711	,	2.91856	,	373.7112	,	747.4224	}	,
			{	21.46	,	6.64345	,	3.73738	,	2.91869	,	373.73819	,	747.47638	}	,
			{	21.48	,	6.64413	,	3.73765	,	2.91883	,	373.76518	,	747.53036	}	,
			{	21.5	,	6.6448	,	3.73792	,	2.91896	,	373.79217	,	747.58435	}	,
			{	21.52	,	6.64548	,	3.73819	,	2.9191	,	373.81916	,	747.63833	}	,
			{	21.54	,	6.64615	,	3.73846	,	2.91923	,	373.84615	,	747.69231	}	,
			{	21.56	,	6.64683	,	3.73873	,	2.91937	,	373.87314	,	747.74629	}	,
			{	21.58	,	6.6475	,	3.739	,	2.9195	,	373.90013	,	747.80027	}	,
			{	21.6	,	6.64818	,	3.73927	,	2.91964	,	373.92713	,	747.85425	}	,
			{	21.62	,	6.64885	,	3.73954	,	2.91977	,	373.95412	,	747.90823	}	,
			{	21.64	,	6.64953	,	3.73981	,	2.91991	,	373.98111	,	747.96221	}	,
			{	21.66	,	6.6502	,	3.74008	,	2.92004	,	374.0081	,	748.01619	}	,
			{	21.68	,	6.65088	,	3.74035	,	2.92018	,	374.03509	,	748.07018	}	,
			{	21.7	,	6.65155	,	3.74062	,	2.92031	,	374.06208	,	748.12416	}	,
			{	21.72	,	6.65223	,	3.74089	,	2.92045	,	374.08907	,	748.17814	}	,
			{	21.74	,	6.6529	,	3.74116	,	2.92058	,	374.11606	,	748.23212	}	,
			{	21.76	,	6.65358	,	3.74143	,	2.92072	,	374.14305	,	748.2861	}	,
			{	21.78	,	6.65425	,	3.7417	,	2.92085	,	374.17004	,	748.34008	}	,
			{	21.8	,	6.65493	,	3.74197	,	2.92099	,	374.19703	,	748.39406	}	,
			{	21.82	,	6.6556	,	3.74224	,	2.92112	,	374.22402	,	748.44804	}	,
			{	21.84	,	6.65628	,	3.74251	,	2.92126	,	374.25101	,	748.50202	}	,
			{	21.86	,	6.65695	,	3.74278	,	2.92139	,	374.278	,	748.55601	}	,
			{	21.88	,	6.65762	,	3.74305	,	2.92152	,	374.30499	,	748.60999	}	,
			{	21.9	,	6.6583	,	3.74332	,	2.92166	,	374.33198	,	748.66397	}	,
			{	21.92	,	6.65897	,	3.74359	,	2.92179	,	374.35897	,	748.71795	}	,
			{	21.94	,	6.65965	,	3.74386	,	2.92193	,	374.38596	,	748.77193	}	,
			{	21.96	,	6.66032	,	3.74413	,	2.92206	,	374.41296	,	748.82591	}	,
			{	21.98	,	6.661	,	3.7444	,	2.9222	,	374.43995	,	748.87989	}	,
			{	22	,	6.66167	,	3.74467	,	2.92233	,	374.46694	,	748.93387	}	,
			{	22.02	,	6.66235	,	3.74494	,	2.92247	,	374.49393	,	748.98785	}	,
			{	22.04	,	6.66302	,	3.74521	,	2.9226	,	374.52092	,	749.04184	}	,
			{	22.06	,	6.6637	,	3.74548	,	2.92274	,	374.54791	,	749.09582	}	,
			{	22.08	,	6.66437	,	3.74575	,	2.92287	,	374.5749	,	749.1498	}	,
			{	22.1	,	6.66505	,	3.74602	,	2.92301	,	374.60189	,	749.20378	}	,
			{	22.12	,	6.66572	,	3.74629	,	2.92314	,	374.62888	,	749.25776	}	,
			{	22.14	,	6.6664	,	3.74656	,	2.92328	,	374.65587	,	749.31174	}	,
			{	22.16	,	6.66707	,	3.74683	,	2.92341	,	374.68286	,	749.36572	}	,
			{	22.18	,	6.66775	,	3.7471	,	2.92355	,	374.70985	,	749.4197	}	,
			{	22.2	,	6.66842	,	3.74737	,	2.92368	,	374.73684	,	749.47368	}	,
			{	22.22	,	6.6691	,	3.74764	,	2.92382	,	374.76383	,	749.52767	}	,
			{	22.24	,	6.66977	,	3.74791	,	2.92395	,	374.79082	,	749.58165	}	,
			{	22.26	,	6.67045	,	3.74818	,	2.92409	,	374.81781	,	749.63563	}	,
			{	22.28	,	6.67112	,	3.74845	,	2.92422	,	374.8448	,	749.68961	}	,
			{	22.3	,	6.67179	,	3.74872	,	2.92436	,	374.87179	,	749.74359	}	,
			{	22.32	,	6.67247	,	3.74899	,	2.92449	,	374.89879	,	749.79757	}	,
			{	22.34	,	6.67314	,	3.74926	,	2.92463	,	374.92578	,	749.85155	}	,
			{	22.36	,	6.67382	,	3.74953	,	2.92476	,	374.95277	,	749.90553	}	,
			{	22.38	,	6.67449	,	3.7498	,	2.9249	,	374.97976	,	749.95951	}	,
			{	22.4	,	6.67517	,	3.75007	,	2.92503	,	375.00675	,	750.0135	}	,
			{	22.42	,	6.67584	,	3.75034	,	2.92517	,	375.03374	,	750.06748	}	,
			{	22.44	,	6.67652	,	3.75061	,	2.9253	,	375.06073	,	750.12146	}	,
			{	22.46	,	6.67719	,	3.75088	,	2.92544	,	375.08772	,	750.17544	}	,
			{	22.48	,	6.67787	,	3.75115	,	2.92557	,	375.11471	,	750.22942	}	,
			{	22.5	,	6.67854	,	3.75142	,	2.92571	,	375.1417	,	750.2834	}	,
			{	22.52	,	6.67922	,	3.75169	,	2.92584	,	375.16869	,	750.33738	}	,
			{	22.54	,	6.67989	,	3.75196	,	2.92598	,	375.19568	,	750.39136	}	,
			{	22.56	,	6.68057	,	3.75223	,	2.92611	,	375.22267	,	750.44534	}	,
			{	22.58	,	6.68124	,	3.7525	,	2.92625	,	375.24966	,	750.49933	}	,
			{	22.6	,	6.68192	,	3.75277	,	2.92638	,	375.27665	,	750.55331	}	,
			{	22.62	,	6.68259	,	3.75304	,	2.92652	,	375.30364	,	750.60729	}	,
			{	22.64	,	6.68327	,	3.75331	,	2.92665	,	375.33063	,	750.66127	}	,
			{	22.66	,	6.68394	,	3.75358	,	2.92679	,	375.35762	,	750.71525	}	,
			{	22.68	,	6.68462	,	3.75385	,	2.92692	,	375.38462	,	750.76923	}	,
			{	22.7	,	6.68529	,	3.75412	,	2.92706	,	375.41161	,	750.82321	}	,
			{	22.72	,	6.68596	,	3.75439	,	2.92719	,	375.4386	,	750.87719	}	,
			{	22.74	,	6.68664	,	3.75466	,	2.92733	,	375.46559	,	750.93117	}	,
			{	22.76	,	6.68731	,	3.75493	,	2.92746	,	375.49258	,	750.98516	}	,
			{	22.78	,	6.68799	,	3.7552	,	2.9276	,	375.51957	,	751.03914	}	,
			{	22.8	,	6.68866	,	3.75547	,	2.92773	,	375.54656	,	751.09312	}	,
			{	22.82	,	6.68934	,	3.75574	,	2.92787	,	375.57355	,	751.1471	}	,
			{	22.84	,	6.69001	,	3.75601	,	2.928	,	375.60054	,	751.20108	}	,
			{	22.86	,	6.69069	,	3.75628	,	2.92814	,	375.62753	,	751.25506	}	,
			{	22.88	,	6.69136	,	3.75655	,	2.92827	,	375.65452	,	751.30904	}	,
			{	22.9	,	6.69204	,	3.75682	,	2.92841	,	375.68151	,	751.36302	}	,
			{	22.92	,	6.69271	,	3.75709	,	2.92854	,	375.7085	,	751.417	}	,
			{	22.94	,	6.69339	,	3.75735	,	2.92868	,	375.73549	,	751.47099	}	,
			{	22.96	,	6.69406	,	3.75762	,	2.92881	,	375.76248	,	751.52497	}	,
			{	22.98	,	6.69474	,	3.75789	,	2.92895	,	375.78947	,	751.57895	}	,
			{	23	,	6.69541	,	3.75816	,	2.92908	,	375.81646	,	751.63293	}	,
			{	23.02	,	6.69609	,	3.75843	,	2.92922	,	375.84345	,	751.68691	}	,
			{	23.04	,	6.69676	,	3.7587	,	2.92935	,	375.87045	,	751.74089	}	,
			{	23.06	,	6.69744	,	3.75897	,	2.92949	,	375.89744	,	751.79487	}	,
			{	23.08	,	6.69811	,	3.75924	,	2.92962	,	375.92443	,	751.84885	}	,
			{	23.1	,	6.69879	,	3.75951	,	2.92976	,	375.95142	,	751.90283	}	,
			{	23.12	,	6.69946	,	3.75978	,	2.92989	,	375.97841	,	751.95682	}	,
			{	23.14	,	6.70013	,	3.76005	,	2.93003	,	376.0054	,	752.0108	}	,
			{	23.16	,	6.70081	,	3.76032	,	2.93016	,	376.03239	,	752.06478	}	,
			{	23.18	,	6.70148	,	3.76059	,	2.9303	,	376.05938	,	752.11876	}	,
			{	23.2	,	6.70216	,	3.76086	,	2.93043	,	376.08637	,	752.17274	}	,
			{	23.22	,	6.70283	,	3.76113	,	2.93057	,	376.11336	,	752.22672	}	,
			{	23.24	,	6.70351	,	3.7614	,	2.9307	,	376.14035	,	752.2807	}	,
			{	23.26	,	6.70418	,	3.76167	,	2.93084	,	376.16734	,	752.33468	}	,
			{	23.28	,	6.70486	,	3.76194	,	2.93097	,	376.19433	,	752.38866	}	,
			{	23.3	,	6.70553	,	3.76221	,	2.93111	,	376.22132	,	752.44265	}	,
			{	23.32	,	6.70621	,	3.76248	,	2.93124	,	376.24831	,	752.49663	}	,
			{	23.34	,	6.70688	,	3.76275	,	2.93138	,	376.2753	,	752.55061	}	,
			{	23.36	,	6.70756	,	3.76302	,	2.93151	,	376.30229	,	752.60459	}	,
			{	23.38	,	6.70823	,	3.76329	,	2.93165	,	376.32928	,	752.65857	}	,
			{	23.4	,	6.70891	,	3.76356	,	2.93178	,	376.35628	,	752.71255	}	,
			{	23.42	,	6.70958	,	3.76383	,	2.93192	,	376.38327	,	752.76653	}	,
			{	23.44	,	6.71026	,	3.7641	,	2.93205	,	376.41026	,	752.82051	}	,
			{	23.46	,	6.71093	,	3.76437	,	2.93219	,	376.43725	,	752.87449	}	,
			{	23.48	,	6.71161	,	3.76464	,	2.93232	,	376.46424	,	752.92848	}	,
			{	23.5	,	6.71228	,	3.76491	,	2.93246	,	376.49123	,	752.98246	}	,
			{	23.52	,	6.71296	,	3.76518	,	2.93259	,	376.51822	,	753.03644	}	,
			{	23.54	,	6.71363	,	3.76545	,	2.93273	,	376.54521	,	753.09042	}	,
			{	23.56	,	6.7143	,	3.76572	,	2.93286	,	376.5722	,	753.1444	}	,
			{	23.58	,	6.71498	,	3.76599	,	2.933	,	376.59919	,	753.19838	}	,
			{	23.6	,	6.71565	,	3.76626	,	2.93313	,	376.62618	,	753.25236	}	,
			{	23.62	,	6.71633	,	3.76653	,	2.93327	,	376.65317	,	753.30634	}	,
			{	23.64	,	6.717	,	3.7668	,	2.9334	,	376.68016	,	753.36032	}	,
			{	23.66	,	6.71768	,	3.76707	,	2.93354	,	376.70715	,	753.4143	}	,
			{	23.68	,	6.71835	,	3.76734	,	2.93367	,	376.73414	,	753.46829	}	,
			{	23.7	,	6.71903	,	3.76761	,	2.93381	,	376.76113	,	753.52227	}	,
			{	23.72	,	6.7197	,	3.76788	,	2.93394	,	376.78812	,	753.57625	}	,
			{	23.74	,	6.72038	,	3.76815	,	2.93408	,	376.81511	,	753.63023	}	,
			{	23.76	,	6.72105	,	3.76842	,	2.93421	,	376.84211	,	753.68421	}	,
			{	23.78	,	6.72173	,	3.76869	,	2.93435	,	376.8691	,	753.73819	}	,
			{	23.8	,	6.7224	,	3.76896	,	2.93448	,	376.89609	,	753.79217	}	,
			{	23.82	,	6.72308	,	3.76923	,	2.93462	,	376.92308	,	753.84615	}	,
			{	23.84	,	6.72375	,	3.7695	,	2.93475	,	376.95007	,	753.90013	}	,
			{	23.86	,	6.72443	,	3.76977	,	2.93489	,	376.97706	,	753.95412	}	,
			{	23.88	,	6.7251	,	3.77004	,	2.93502	,	377.00405	,	754.0081	}	,
			{	23.9	,	6.72578	,	3.77031	,	2.93516	,	377.03104	,	754.06208	}	,
			{	23.92	,	6.72645	,	3.77058	,	2.93529	,	377.05803	,	754.11606	}	,
			{	23.94	,	6.72713	,	3.77085	,	2.93543	,	377.08502	,	754.17004	}	,
			{	23.96	,	6.7278	,	3.77112	,	2.93556	,	377.11201	,	754.22402	}	,
			{	23.98	,	6.72848	,	3.77139	,	2.9357	,	377.139	,	754.278	}	,
			{	24	,	6.72915	,	3.77166	,	2.93583	,	377.16599	,	754.33198	}	,
			{	24.02	,	6.72982	,	3.77193	,	2.93596	,	377.19298	,	754.38596	}	,
			{	24.04	,	6.7305	,	3.7722	,	2.9361	,	377.21997	,	754.43995	}	,
			{	24.06	,	6.73117	,	3.77247	,	2.93623	,	377.24696	,	754.49393	}	,
			{	24.08	,	6.73185	,	3.77274	,	2.93637	,	377.27395	,	754.54791	}	,
			{	24.1	,	6.73252	,	3.77301	,	2.9365	,	377.30094	,	754.60189	}	,
			{	24.12	,	6.7332	,	3.77328	,	2.93664	,	377.32794	,	754.65587	}	,
			{	24.14	,	6.73387	,	3.77355	,	2.93677	,	377.35493	,	754.70985	}	,
			{	24.16	,	6.73455	,	3.77382	,	2.93691	,	377.38192	,	754.76383	}	,
			{	24.18	,	6.73522	,	3.77409	,	2.93704	,	377.40891	,	754.81781	}	,
			{	24.2	,	6.7359	,	3.77436	,	2.93718	,	377.4359	,	754.87179	}	,
			{	24.22	,	6.73657	,	3.77463	,	2.93731	,	377.46289	,	754.92578	}	,
			{	24.24	,	6.73725	,	3.7749	,	2.93745	,	377.48988	,	754.97976	}	,
			{	24.26	,	6.73792	,	3.77517	,	2.93758	,	377.51687	,	755.03374	}	,
			{	24.28	,	6.7386	,	3.77544	,	2.93772	,	377.54386	,	755.08772	}	,
			{	24.3	,	6.73927	,	3.77571	,	2.93785	,	377.57085	,	755.1417	}	,
			{	24.32	,	6.73995	,	3.77598	,	2.93799	,	377.59784	,	755.19568	}	,
			{	24.34	,	6.74062	,	3.77625	,	2.93812	,	377.62483	,	755.24966	}	,
			{	24.36	,	6.7413	,	3.77652	,	2.93826	,	377.65182	,	755.30364	}	,
			{	24.38	,	6.74197	,	3.77679	,	2.93839	,	377.67881	,	755.35762	}	,
			{	24.4	,	6.74265	,	3.77706	,	2.93853	,	377.7058	,	755.41161	}	,
			{	24.42	,	6.74332	,	3.77733	,	2.93866	,	377.73279	,	755.46559	}	,
			{	24.44	,	6.74399	,	3.7776	,	2.9388	,	377.75978	,	755.51957	}	,
			{	24.46	,	6.74467	,	3.77787	,	2.93893	,	377.78677	,	755.57355	}	,
			{	24.48	,	6.74534	,	3.77814	,	2.93907	,	377.81377	,	755.62753	}	,
			{	24.5	,	6.74602	,	3.77841	,	2.9392	,	377.84076	,	755.68151	}	,
			{	24.52	,	6.74669	,	3.77868	,	2.93934	,	377.86775	,	755.73549	}	,
			{	24.54	,	6.74737	,	3.77895	,	2.93947	,	377.89474	,	755.78947	}	,
			{	24.56	,	6.74804	,	3.77922	,	2.93961	,	377.92173	,	755.84345	}	,
			{	24.58	,	6.74872	,	3.77949	,	2.93974	,	377.94872	,	755.89744	}	,
			{	24.6	,	6.74939	,	3.77976	,	2.93988	,	377.97571	,	755.95142	}	,
			{	24.62	,	6.75007	,	3.78003	,	2.94001	,	378.0027	,	756.0054	}	,
			{	24.64	,	6.75074	,	3.7803	,	2.94015	,	378.02969	,	756.05938	}	,
			{	24.66	,	6.75142	,	3.78057	,	2.94028	,	378.05668	,	756.11336	}	,
			{	24.68	,	6.75209	,	3.78084	,	2.94042	,	378.08367	,	756.16734	}	,
			{	24.7	,	6.75277	,	3.78111	,	2.94055	,	378.11066	,	756.22132	}	,
			{	24.72	,	6.75344	,	3.78138	,	2.94069	,	378.13765	,	756.2753	}	,
			{	24.74	,	6.75412	,	3.78165	,	2.94082	,	378.16464	,	756.32928	}	,
			{	24.76	,	6.75479	,	3.78192	,	2.94096	,	378.19163	,	756.38327	}	,
			{	24.78	,	6.75547	,	3.78219	,	2.94109	,	378.21862	,	756.43725	}	,
			{	24.8	,	6.75614	,	3.78246	,	2.94123	,	378.24561	,	756.49123	}	,
			{	24.82	,	6.75682	,	3.78273	,	2.94136	,	378.2726	,	756.54521	}	,
			{	24.84	,	6.75749	,	3.783	,	2.9415	,	378.2996	,	756.59919	}	,
			{	24.86	,	6.75816	,	3.78327	,	2.94163	,	378.32659	,	756.65317	}	,
			{	24.88	,	6.75884	,	3.78354	,	2.94177	,	378.35358	,	756.70715	}	,
			{	24.9	,	6.75951	,	3.78381	,	2.9419	,	378.38057	,	756.76113	}	,
			{	24.92	,	6.76019	,	3.78408	,	2.94204	,	378.40756	,	756.81511	}	,
			{	24.94	,	6.76086	,	3.78435	,	2.94217	,	378.43455	,	756.8691	}	,
			{	24.96	,	6.76154	,	3.78462	,	2.94231	,	378.46154	,	756.92308	}	,
			{	24.98	,	6.76221	,	3.78489	,	2.94244	,	378.48853	,	756.97706	}	,
			{	25	,	6.76289	,	3.78516	,	2.94258	,	378.51552	,	757.03104	}	,
			{	25.02	,	6.76356	,	3.78543	,	2.94271	,	378.54251	,	757.08502	}	,
			{	25.04	,	6.76424	,	3.7857	,	2.94285	,	378.5695	,	757.139	}	,
			{	25.06	,	6.76491	,	3.78596	,	2.94298	,	378.59649	,	757.19298	}	,
			{	25.08	,	6.76559	,	3.78623	,	2.94312	,	378.62348	,	757.24696	}	,
			{	25.1	,	6.76626	,	3.7865	,	2.94325	,	378.65047	,	757.30094	}	,
			{	25.12	,	6.76694	,	3.78677	,	2.94339	,	378.67746	,	757.35493	}	,
			{	25.14	,	6.76761	,	3.78704	,	2.94352	,	378.70445	,	757.40891	}	,
			{	25.16	,	6.76829	,	3.78731	,	2.94366	,	378.73144	,	757.46289	}	,
			{	25.18	,	6.76896	,	3.78758	,	2.94379	,	378.75843	,	757.51687	}	,
			{	25.2	,	6.76964	,	3.78785	,	2.94393	,	378.78543	,	757.57085	}	,
			{	25.22	,	6.77031	,	3.78812	,	2.94406	,	378.81242	,	757.62483	}	,
			{	25.24	,	6.77099	,	3.78839	,	2.9442	,	378.83941	,	757.67881	}	,
			{	25.26	,	6.77166	,	3.78866	,	2.94433	,	378.8664	,	757.73279	}	,
			{	25.28	,	6.77233	,	3.78893	,	2.94447	,	378.89339	,	757.78677	}	,
			{	25.3	,	6.77301	,	3.7892	,	2.9446	,	378.92038	,	757.84076	}	,
			{	25.32	,	6.77368	,	3.78947	,	2.94474	,	378.94737	,	757.89474	}	,
			{	25.34	,	6.77436	,	3.78974	,	2.94487	,	378.97436	,	757.94872	}	,
			{	25.36	,	6.77503	,	3.79001	,	2.94501	,	379.00135	,	758.0027	}	,
			{	25.38	,	6.77571	,	3.79028	,	2.94514	,	379.02834	,	758.05668	}	,
			{	25.4	,	6.77638	,	3.79055	,	2.94528	,	379.05533	,	758.11066	}	,
			{	25.42	,	6.77706	,	3.79082	,	2.94541	,	379.08232	,	758.16464	}	,
			{	25.44	,	6.77773	,	3.79109	,	2.94555	,	379.10931	,	758.21862	}	,
			{	25.46	,	6.77841	,	3.79136	,	2.94568	,	379.1363	,	758.2726	}	,
			{	25.48	,	6.77908	,	3.79163	,	2.94582	,	379.16329	,	758.32659	}	,
			{	25.5	,	6.77976	,	3.7919	,	2.94595	,	379.19028	,	758.38057	}	,
			{	25.52	,	6.78043	,	3.79217	,	2.94609	,	379.21727	,	758.43455	}	,
			{	25.54	,	6.78111	,	3.79244	,	2.94622	,	379.24426	,	758.48853	}	,
			{	25.56	,	6.78178	,	3.79271	,	2.94636	,	379.27126	,	758.54251	}	,
			{	25.58	,	6.78246	,	3.79298	,	2.94649	,	379.29825	,	758.59649	}	,
			{	25.6	,	6.78313	,	3.79325	,	2.94663	,	379.32524	,	758.65047	}	,
			{	25.62	,	6.78381	,	3.79352	,	2.94676	,	379.35223	,	758.70445	}	,
			{	25.64	,	6.78448	,	3.79379	,	2.9469	,	379.37922	,	758.75843	}	,
			{	25.66	,	6.78516	,	3.79406	,	2.94703	,	379.40621	,	758.81242	}	,
			{	25.68	,	6.78583	,	3.79433	,	2.94717	,	379.4332	,	758.8664	}	,
			{	25.7	,	6.7865	,	3.7946	,	2.9473	,	379.46019	,	758.92038	}	,
			{	25.72	,	6.78718	,	3.79487	,	2.94744	,	379.48718	,	758.97436	}	,
			{	25.74	,	6.78785	,	3.79514	,	2.94757	,	379.51417	,	759.02834	}	,
			{	25.76	,	6.78853	,	3.79541	,	2.94771	,	379.54116	,	759.08232	}	,
			{	25.78	,	6.7892	,	3.79568	,	2.94784	,	379.56815	,	759.1363	}	,
			{	25.8	,	6.78988	,	3.79595	,	2.94798	,	379.59514	,	759.19028	}	,
			{	25.82	,	6.79055	,	3.79622	,	2.94811	,	379.62213	,	759.24426	}	,
			{	25.84	,	6.79123	,	3.79649	,	2.94825	,	379.64912	,	759.29825	}	,
			{	25.86	,	6.7919	,	3.79676	,	2.94838	,	379.67611	,	759.35223	}	,
			{	25.88	,	6.79258	,	3.79703	,	2.94852	,	379.7031	,	759.40621	}	,
			{	25.9	,	6.79325	,	3.7973	,	2.94865	,	379.73009	,	759.46019	}	,
			{	25.92	,	6.79393	,	3.79757	,	2.94879	,	379.75709	,	759.51417	}	,
			{	25.94	,	6.7946	,	3.79784	,	2.94892	,	379.78408	,	759.56815	}	,
			{	25.96	,	6.79528	,	3.79811	,	2.94906	,	379.81107	,	759.62213	}	,
			{	25.98	,	6.79595	,	3.79838	,	2.94919	,	379.83806	,	759.67611	}	,
			{	26	,	6.79663	,	3.79865	,	2.94933	,	379.86505	,	759.73009	}	,
			{	26.02	,	6.7973	,	3.79892	,	2.94946	,	379.89204	,	759.78408	}	,
			{	26.04	,	6.79798	,	3.79919	,	2.9496	,	379.91903	,	759.83806	}	,
			{	26.06	,	6.79865	,	3.79946	,	2.94973	,	379.94602	,	759.89204	}	,
			{	26.08	,	6.79933	,	3.79973	,	2.94987	,	379.97301	,	759.94602	}	,
			{	26.1	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.12	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.14	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.16	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.18	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.2	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.22	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.24	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.26	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.28	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.3	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.32	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.34	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.36	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.38	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.4	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.42	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.44	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.46	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.48	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.5	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.52	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.54	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.56	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.58	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.6	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.62	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.64	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.66	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.68	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.7	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.72	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.74	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.76	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.78	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.8	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.82	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.84	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.86	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.88	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.9	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.92	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.94	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.96	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	26.98	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.02	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.04	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.06	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.08	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.1	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.12	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.14	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.16	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.18	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.2	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.22	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.24	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.26	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.28	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.3	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.32	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.34	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.36	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.38	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.4	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.42	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.44	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.46	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.48	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.5	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.52	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.54	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.56	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.58	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.6	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.62	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.64	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.66	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.68	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.7	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.72	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.74	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.76	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.78	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.8	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.82	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.84	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.86	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.88	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.9	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.92	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.94	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.96	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	27.98	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.02	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.04	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.06	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.08	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.1	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.12	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.14	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.16	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.18	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.2	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.22	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.24	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.26	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.28	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.3	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.32	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.34	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.36	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.38	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.4	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.42	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.44	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.46	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.48	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.5	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.52	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.54	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.56	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.58	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.6	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.62	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.64	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.66	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.68	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.7	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.72	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.74	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.76	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.78	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.8	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.82	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.84	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.86	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.88	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.9	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.92	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.94	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.96	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	28.98	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.02	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.04	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.06	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.08	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.1	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.12	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.14	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.16	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.18	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.2	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.22	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.24	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.26	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.28	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.3	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.32	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.34	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.36	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.38	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.4	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.42	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.44	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.46	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.48	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.5	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.52	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.54	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.56	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.58	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.6	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.62	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.64	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.66	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.68	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.7	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.72	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.74	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.76	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.78	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.8	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.82	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.84	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.86	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.88	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.9	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.92	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.94	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.96	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	29.98	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.02	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.04	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.06	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.08	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.1	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.12	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.14	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.16	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.18	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.2	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.22	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.24	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.26	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.28	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.3	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.32	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.34	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.36	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.38	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.4	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.42	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.44	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.46	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.48	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.5	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.52	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.54	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.56	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.58	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.6	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.62	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.64	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.66	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.68	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.7	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.72	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.74	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.76	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.78	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.8	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.82	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.84	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.86	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.88	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.9	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.92	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.94	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.96	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	30.98	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
			{	31	,	6.8	,	3.8	,	2.95	,	380	,	760	}	,
	};

	if (depth < 0){
		return 1;
	}else if(depth >= local_1DVelocityModel[1550][1]){

		g_props->Vp  = local_1DVelocityModel[1550][2];
		g_props->Vs  = local_1DVelocityModel[1550][3];
		g_props->rho = local_1DVelocityModel[1550][4];

		return 0;
	}else{
		for(i = 0; i < 1550; i++)
		{

			if (depth == local_1DVelocityModel[i][1]){
				g_props->Vp  = local_1DVelocityModel[i][2];
				g_props->Vs  = local_1DVelocityModel[i][3];
				g_props->rho = local_1DVelocityModel[i][4];
				return 0;
			}else if (depth > local_1DVelocityModel[i][1] && depth < local_1DVelocityModel[i+1][1]){

				g_props->Vp  = local_1DVelocityModel[i+1][2]-
						((local_1DVelocityModel[i+1][1]-depth)*(local_1DVelocityModel[i+1][2]-local_1DVelocityModel[i][2]))/
						( local_1DVelocityModel[i+1][1]- local_1DVelocityModel[i][1]);

				g_props->Vs  = local_1DVelocityModel[i+1][3]-
						((local_1DVelocityModel[i+1][1]-depth)*(local_1DVelocityModel[i+1][3]-local_1DVelocityModel[i][3]))/
						( local_1DVelocityModel[i+1][1]- local_1DVelocityModel[i][1]);

				g_props->rho = local_1DVelocityModel[i+1][4]-
						((local_1DVelocityModel[i+1][1]-depth)*(local_1DVelocityModel[i+1][4]-local_1DVelocityModel[i][4]))/
						( local_1DVelocityModel[i+1][1]- local_1DVelocityModel[i][1]);
				return 0;
			}else{

				return 1;

			}
			return 1;
		}
	}
}


/**
 * compute_setflag:
 *
 * - results from the discussion with Leo
 * - set flags as if in the full space.
 * - the main() routine will set the flag properly in case half-space
 *   is desired.
 *
 */
#ifdef BOUNDARY
static char
compute_setflag(tick_t ldb[3], tick_t ruf[3], tick_t p1[3], tick_t p2[3])
{
    char flag;

    flag = 13; /* indicate internal element */

    if (ldb[0] == p1[0])
	flag = 12;
    if (ldb[1] == p1[1])
	flag = 10;
    if (ldb[2] == p1[2])
	flag = 4;

    if (ruf[0] == p2[0])
	flag = 14;
    if (ruf[1] == p2[1])
	flag = 16;
    if (ruf[2] == p2[2])
	flag = 22;


    if(ldb[0] == p1[0] && ldb[1] == p1[1])
	flag = 9;

    if(ruf[0] == p2[0] && ldb[1] == p1[1])
	flag = 11;

    if(ldb[0] == p1[0] && ruf[1] == p2[1])
	flag = 15;

    if(ruf[0] == p2[0] &&   ruf[1] == p2[1])
	flag = 17;


    if (ldb[0] == p1[0] && ldb[2] == p1[2])
	flag = 3;

    if (ruf[0] == p2[0] && ldb[2] == p1[2])
	flag = 5;

    if (ldb[0] == p1[0] && ruf[2] == p2[2])
	flag = 21;

    if (ruf[0] == p2[0] && ruf[2] == p2[2])
	flag = 23;


    if (ldb[1] == p1[1] && ldb[2] == p1[2])
	flag = 1;

    if (ruf[1] == p2[1] && ldb[2] == p1[2])
	flag = 7;

    if (ldb[1] == p1[1] && ruf[2] == p2[2])
	flag = 19;

    if (ruf[1] == p2[1] && ruf[2] == p2[2])
	flag = 25;

    if (ldb[0] == p1[0] && ldb[1] == p1[1] && ldb[2] == p1[2])
	flag = 0;

    if (ruf[0] == p2[0] && (ldb[1] == p1[1]) && ldb[2] == p1[2])
	flag = 2;

    if (ldb[0] == p1[0] && ruf[1] == p2[1] && ldb[2] == p1[2])
	flag = 6;

    if (ruf[0] == p2[0] && ruf[1] == p2[1] && ldb[2] == p1[2])
	flag = 8;

    if (ldb[0] == p1[0] && ldb[1] == p1[1] && ruf[2] == p2[2])
	flag = 18;

    if (ruf[0] == p2[0] && ldb[1] == p1[1] && ruf[2] == p2[2])
	flag = 20;

    if (ldb[0] == p1[0] && ruf[1] == p2[1] && ruf[2] == p2[2])
	flag = 24;

    if (ruf[0] == p2[0] && ruf[1] == p2[1] && ruf[2] == p2[2])
	flag = 26;

    return flag;
}



static const int theIDBoundaryMatrix[27][8] = {
    { 7, 6, 5, 4, 3, 2, 1, 0},
    { 6, 6, 4, 4, 2, 2, 0, 0},
    { 6, 7, 4, 5, 2, 3, 0, 1},
    { 5, 4, 5, 4, 1, 0, 1, 0},
    { 4, 4, 4, 4, 0, 0, 0, 0},
    { 4, 5, 4, 5, 0, 1, 0, 1},
    { 5, 4, 7, 6, 1, 0, 3, 2},
    { 4, 4, 6, 6, 0, 0, 2, 2},
    { 4, 5, 6, 7, 0, 1, 2, 3},
    { 3, 2, 1, 0, 3, 2, 1, 0},
    { 2, 2, 0, 0, 2, 2, 0, 0},
    { 2, 3, 0, 1, 2, 3, 0, 1},
    { 1, 0, 1, 0, 1, 0, 1, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0}, /* 13: internal elements */
    { 0, 1, 0, 1, 0, 1, 0, 1},
    { 1, 0, 3, 2, 1, 0, 3, 2},
    { 0, 0, 2, 2, 0, 0, 2, 2},
    { 0, 1, 2, 3, 0, 1, 2, 3},
    { 3, 2, 1, 0, 7, 6, 5, 4},
    { 2, 2, 0, 0, 6, 6, 4, 4},
    { 2, 3, 0, 1, 6, 7, 4, 5},
    { 1, 0, 1, 0, 5, 4, 5, 4},
    { 0, 0, 0, 0, 4, 4, 4, 4},
    { 0, 1, 0, 1, 4, 5, 4, 5},
    { 1, 0, 3, 2, 5, 4, 7, 6},
    { 0, 0, 2, 2, 4, 4, 6, 6},
    { 0, 1, 2, 3, 4, 5, 6, 7},
};





static void
compute_setboundary(float size, float Vp, float Vs, float rho, int flag,
		    double dashpot[8][3])
{
    int whichNode;
    double scale;

    /* init the damping vector to all zeroes */
    memset(dashpot, 0, sizeof(double) * 8 * 3);

#ifdef HALFSPACE
    flag = (flag < 9) ? flag + 9 : flag;
#endif /* HALFSPACE */

    scale = rho * (size / 2) * (size / 2);

    for (whichNode = 0 ; whichNode < 8; whichNode++) {
	int bitmark, component;

	bitmark = theIDBoundaryMatrix[flag][whichNode];

	switch (bitmark) {
	case 0:
	    break;
	case 7:
	    /* Three contributing faces */
	    dashpot[whichNode][0] = dashpot[whichNode][1]
		= dashpot[whichNode][2] = (Vp + 2 * Vs) * scale;
	    break;
	case 3:
	case 5:
	case 6:
	    /* Two contributing faces */
	    for (component = 0; component < 3; component++)
		dashpot[whichNode][component] =
		    (Vs + ((bitmark & (1<< component)) ? Vp : Vs)) * scale;
	    break;
	case 1:
	case 2:
	case 4:
	    /* One contributing face */
	    for (component = 0; component < 3; component++)
		dashpot[whichNode][component] =
		    ((bitmark & (1<<component)) ? Vp : Vs) * scale;
	    break;
	default:
	    fprintf(stderr, "SetBoundary: Unknown bitmark. Panic!\n");
	    exit(1);
	}
    }

    return;
}
#endif /* BOUNDARY */



/**
 * compute_setab: the base a and b values will be scaled by zeta
 *                specific to each element.
 */
static void compute_setab(double freq, double *aBasePtr, double *bBasePtr)
{
    /* old version which caused overflow because of the aproximation in
     * the derivative */

    double w1, w2, lw1, lw2, sw1, sw2, cw1, cw2;
    double numer, denom;

    if (Param.theTypeOfDamping == RAYLEIGH)
    {
	/* the factors 0.2 and 1 were calibrated heuristically by LEO */
	w1 = 2 * PI * freq *.2;
	w2 = 2 * PI * freq * 1;

	/* logs */
	lw1 = log(w1);
	lw2 = log(w2);

	/* squares */
	sw1 = w1 * w1;
	sw2 = w2 * w2;

	/* cubes */
	cw1 = w1 * w1 * w1;
	cw2 = w2 * w2 * w2;

	/* numerator */
	numer = w1 * w2 *
	    ( -2 * sw1 * lw2 + 2 * sw1 * lw1 - 2 * w1 * w2 * lw2
	      + 2 * w1 * w2 * lw1 + 3 * sw2 - 3 * sw1
		  - 2 * sw2 * lw2 + 2 * sw2 * lw1);

	/* denominator */
	denom = (cw1 - cw2 + 3 * sw2 * w1 - 3 * sw1 * w2);

	/* the a over zeta target is... */
	*aBasePtr = numer / denom;

	/* new numerator */
	numer = 3 * (2 * w1 * w2 * lw2 - 2 * w1 * w2 * lw1 + sw1 - sw2);

	/* the b over zeta target is... */
	*bBasePtr = numer / denom;

    }
    else if ( Param.theTypeOfDamping == MASS )
    {
	w1 = 2 * PI * freq * .1;  /* these .1 and 8 heuristics */
	w2 = 2 * PI * freq * 8;

	numer = 2 * w2 * w1 * log(w2 / w1);
	denom = w2 - w1;

	*aBasePtr = 1.3*numer / denom;  /* this 1.3 comes out from heuristics */
	*bBasePtr = 0;
    }
    else if ( Param.theTypeOfDamping == NONE || Param.theTypeOfDamping == BKT )
    {
	*aBasePtr = 0;
	*bBasePtr = 0;
    }

    return;
}



int
is_nodeloaded( int32_t iNode, char* onoff )
{
    /* \todo use a general bitmap data structure for this */
    /* \todo move the routine declaration to the source generation file */

    int32_t whichByte, whichBit;

    char mask,test;

    whichByte = iNode/8;
    whichBit = 7 - iNode % 8;

    mask = ( char )pow(2,whichBit);

    test = onoff[whichByte] & mask;

    return (test == mask);  /* 1 if equal, 0 otherwise */
}


/**
 * Add the force due to earthquake source.
 *
 * Globals accessed:
 * - Global.mySolver->force (W)
 * - Param.theDeltaTSquared (R)
 * - Global.myForces (R)
 *
 * Iterates over the nodes that are loaded by the source and set the
 * respective forces for those nodes.
 */
static void
compute_addforce_s( int32_t timestep )
{
    int i;	/* index for loaded nodes (from the source) */

    for (i = 0; i <  Global.theNodesLoaded; i++) {
	int lnid = Global.theNodesLoadedList[i];	/* local node id */

	/* node's force vector */
	fvector_t* nodalForce =	Global.mySolver->force + lnid;

	/* vector-scalar multiply */
	nodalForce->f[0] = ( Global.myForces [ i ].x [0] ) * Param.theDeltaTSquared;
	nodalForce->f[1] = ( Global.myForces [ i ].x [1] ) * Param.theDeltaTSquared;
	nodalForce->f[2] = ( Global.myForces [ i ].x [2] ) * Param.theDeltaTSquared;
    }
}

/**
 * compute_adjust: Either distribute the values from LOCAL dangling nodes
 *                 to LOCAL anchored nodes, or assign values from LOCAL
 *                 anchored nodes to LOCAL dangling nodes.
 *
 */
static void
compute_adjust(void *valuetable, int32_t itemsperentry, int32_t how)
{
    solver_float *vtable = (solver_float *)valuetable;
    int32_t dnindex;

    if (how == DISTRIBUTION) {
	for (dnindex = 0; dnindex < Global.myMesh->ldnnum; dnindex++) {
	    dnode_t *dnode;
	    solver_float *myvalue, *parentvalue;
	    solver_float darray[7]; /* A hack to avoid memory allocation */
	    int32link_t *int32link;
	    int32_t idx, parentlnid;
#ifdef DEBUG
	    int32_t deps = 0;
#endif /* DEBUG */

	    dnode = &Global.myMesh->dnodeTable[dnindex];
	    myvalue = vtable + dnode->ldnid * itemsperentry;

	    for (idx = 0; idx < itemsperentry; idx++) {
		darray[idx] = (*(myvalue + idx)) / dnode->deps;
	    }

	    /* Distribute my darray value to my anchors */
	    int32link = dnode->lanid;
	    while (int32link != NULL) {

#ifdef DEBUG
		deps++;
#endif

		parentlnid = int32link->id;
		parentvalue = vtable + parentlnid * itemsperentry;

		for (idx = 0; idx < itemsperentry; idx++) {
		    /* Accumulation the distributed values */
		    *(parentvalue + idx) += darray[idx];
		}

		int32link = int32link->next;
	    }

#ifdef DEBUG
	    if (deps != (int)dnode->deps) {
		fprintf(stderr, "Thread %d: compute_adjust distri: ", Global.myID);
		fprintf(stderr, "deps don't match\n");
		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	    }
#endif /* DEBUG */
	} /* for all my LOCAL dangling nodes */

    } else {
	/* Assign the value of the anchored parents to the dangling nodes*/

	for (dnindex = 0; dnindex < Global.myMesh->ldnnum; dnindex++) {
	    dnode_t *dnode;
	    solver_float *myvalue, *parentvalue;
	    int32link_t *int32link;
	    int32_t idx, parentlnid;
#ifdef DEBUG
	    int32_t deps = 0;
#endif /* DEBUG */

	    dnode = &Global.myMesh->dnodeTable[dnindex];
	    myvalue = vtable + dnode->ldnid * itemsperentry;

	    /* Zero out the residual values the dangling node might
	       still hold */
	    memset(myvalue, 0, sizeof(solver_float) * itemsperentry);

	    /* Assign prorated anchored values to a dangling node */
	    int32link = dnode->lanid;
	    while (int32link != NULL) {

#ifdef DEBUG
		deps++;
#endif

		parentlnid = int32link->id;
		parentvalue = vtable + parentlnid * itemsperentry;

		for (idx = 0; idx < itemsperentry; idx++) {
		    *(myvalue + idx) += (*(parentvalue + idx) / dnode->deps);
		}

		int32link = int32link->next;
	    }

#ifdef DEBUG
	    if (deps != (int)dnode->deps) {
		fprintf(stderr, "Thread %d: compute_adjust assign: ", Global.myID);
		fprintf(stderr, "deps don't match\n");
		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	    }
#endif /* DEBUG */

	} /* for all my LOCAL dangling nodes */
    }

    return;
}

static void print_timing_stat()
{

    double TotalMeshingTime;
    TotalMeshingTime = Timer_Value("Octor Newtree"          , 0)
                     + Timer_Value("Octor Refinetree"       , 0)
                     + Timer_Value("Octor Balancetree"      , 0)
                     + Timer_Value("Octor Partitiontree"    , 0)
                     + Timer_Value("Octor Extractmesh"      , 0)
                     + Timer_Value("Mesh correct properties", 0)
                     + Timer_Value("Mesh Stats Print"       , 0);

    if ( Param.includeBuildings == YES ) {
        TotalMeshingTime += Timer_Value("Carve Buildings", 0);
    }

    printf("\n\n__________________________Raw Timers__________________________\n\n");
    Timer_PrintAll(comm_solver);
    printf("\n\n\n\n\n");


    printf("\n\n__________________________Timer Statistics__________________________\n\n");


    printf("\n_____________Summary_____________\n");
    printf("Max Frequency             : %.2f\n", Param.theFreq);
    printf("Vs                        : %.2f\n", Param.theVsCut);
    printf("Total elements            : %" INT64_FMT "\n", Global.theETotal);
    printf("Elements/PE               : %" INT64_FMT "\n", Global.theETotal/Global.theGroupSize);
    printf("Simulation duration       : %.2f seconds\n", Param.theEndT - Param.theStartT);
    printf("Total steps               : %d\n", Param.theTotalSteps);
    printf("DeltaT used               : %.6f seconds\n", Param.theDeltaT);
    printf("Critical deltaT           : %.6f seconds\n", Global.theCriticalT);
    printf("\n");
    printf("Total Wall Clock          : %.2f seconds\n", Timer_Value("Total Wall Clock",0));
    printf("Time/step                 : %.6f seconds\n", Timer_Value("Solver",0)/Param.theTotalSteps);
    printf("Time/step/(elem/PE)       : %.6f millisec\n",Timer_Value("Solver",0) * 1000.0 / Param.theTotalSteps /
	   (Global.theETotal * 1.0 / Global.theGroupSize));
    printf("Simulation Rate Variation : %.3f (Average)   %.3f (Min)   %.3f (Max)  (sec/%d timesteps)\n",
	   (Timer_Value("Solver",0)/Param.theTotalSteps)*Param.monitor_stats_rate,
	   Global.fastestTimeSteps, Global.slowestTimeSteps, Param.monitor_stats_rate);
    printf("\n");


    printf("\n____________Breakdown____________\n");
    printf("TOTAL MESHING                       : %.2f seconds\n", TotalMeshingTime);
    printf("    Octor Newtree                   : %.2f seconds\n", Timer_Value("Octor Newtree",0) );
    printf("    Octor Refinetree                : %.2f seconds\n", Timer_Value("Octor Refinetree",0));
    printf("    Octor Balancetree               : %.2f seconds\n", Timer_Value("Octor Balancetree",0));
    if ( Timer_Exists("Carve Buildings") )
        printf("    Octor Carve Buildings           : %.2f seconds\n", Timer_Value("Carve Buildings",0));
    printf("    Octor Partitiontree             : %.2f seconds\n", Timer_Value("Octor Partitiontree",0));
    printf("    Octor Extractmesh               : %.2f seconds\n", Timer_Value("Octor Extractmesh",0));
    printf("    Mesh correct properties         : %.2f seconds\n", Timer_Value("Mesh correct properties",0));
    printf("    Mesh Stats Print                : %.2f seconds\n", Timer_Value("Mesh Stats Print",0));
    printf("\n");

    if(Param.drmImplement == YES) {

    	printf("DRM INIT PARAMETERS                 : %.2f (Max) %.2f (Min) seconds\n",
    			Timer_Value("Init Drm Parameters",MAX), Timer_Value("Init Drm Parameters",MIN));
    	printf("\n");

    	printf("DRM INITIALIZATION                  : %.2f (Max) %.2f (Min) seconds\n",
    			Timer_Value("Drm Init",MAX), Timer_Value("Drm Init",MIN));

    	if(Param.theDrmPart == PART2) {

    		printf("    Find Drm File To Readjust       : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
    				Timer_Value("Find Drm File To Readjust",AVERAGE),
    				Timer_Value("Find Drm File To Readjust",MAX),
    				Timer_Value("Find Drm File To Readjust",MIN));

    		printf("    Fill Drm Struct                 : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
    				Timer_Value("Fill Drm Struct",AVERAGE),
    				Timer_Value("Fill Drm Struct",MAX),
    				Timer_Value("Fill Drm Struct",MIN));

    		printf("    Comm of Drm Coordinates         : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
    				Timer_Value("Comm of Drm Coordinates",AVERAGE),
    				Timer_Value("Comm of Drm Coordinates",MAX),
    				Timer_Value("Comm of Drm Coordinates",MIN));

    		printf("    Find Which Drm Files To Print   : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
    				Timer_Value("Find Which Drm Files To Print",AVERAGE),
    				Timer_Value("Find Which Drm Files To Print",MAX),
    				Timer_Value("Find Which Drm Files To Print",MIN));

    		printf("    Read And Rearrange Drm Files    : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
    				Timer_Value("Read And Rearrange Drm Files",AVERAGE),
    				Timer_Value("Read And Rearrange Drm Files",MAX),
    				Timer_Value("Read And Rearrange Drm Files",MIN));

    		printf("    Find Which Drm Files To Read    : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
    				Timer_Value("Find Which Drm Files To Read",AVERAGE),
    				Timer_Value("Find Which Drm Files To Read",MAX),
    				Timer_Value("Find Which Drm Files To Read",MIN));

    		printf("    Locate where I am in file       : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
    				Timer_Value("Locate where I am in file",AVERAGE),
    				Timer_Value("Locate where I am in file",MAX),
    				Timer_Value("Locate where I am in file",MIN));

    	}
    	printf("\n");
    }

    printf("SOURCE INITIALIZATION               : %.2f (Max) %.2f (Min) seconds\n",
	   Timer_Value("Source Init",MAX), Timer_Value("Source Init",MIN));
    printf("\n");
    printf("TOTAL SOLVER                        : %.2f seconds\n", Timer_Value("Solver",0));
    printf("    Read My Forces                  : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
	   Timer_Value("Read My Forces",AVERAGE),
	   Timer_Value("Read My Forces",MAX),
	   Timer_Value("Read My Forces",MIN));
    printf("    Compute addforces s             : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
	   Timer_Value("Compute addforces s",AVERAGE),
	   Timer_Value("Compute addforces s",MAX),
	   Timer_Value("Compute addforces s",MIN));
    printf("    Compute addforces e             : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
	   Timer_Value("Compute addforces e",AVERAGE),
	   Timer_Value("Compute addforces e",MAX),
	   Timer_Value("Compute addforces e",MIN));
    printf("    Compute Damping addforce        : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
	   Timer_Value("Damping addforce",AVERAGE),
	   Timer_Value("Damping addforce",MAX),
	   Timer_Value("Damping addforce",MIN));
    if ( Timer_Exists("Compute Non-linear Entities") )
        printf("    Compute Non-linear Entities     : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
                Timer_Value("Compute Non-linear Entities",AVERAGE),
                Timer_Value("Compute Non-linear Entities",MAX),
                Timer_Value("Compute Non-linear Entities",MIN));
    if ( Timer_Exists("Compute addforces Non-linear") ) {
        printf("    Compute addforces Non-linear    : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
                Timer_Value("Compute addforces Non-linear",AVERAGE),
                Timer_Value("Compute addforces Non-linear",MAX),
                Timer_Value("Compute addforces Non-linear",MIN));
        printf("    Compute addforces gravity       : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
                Timer_Value("Compute addforces gravity",AVERAGE),
                Timer_Value("Compute addforces gravity",MAX),
                Timer_Value("Compute addforces gravity",MIN));
    }

    if ( Timer_Exists("Solver drm force compute") ) {
    	printf("    Solver drm force compute        : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
    			Timer_Value("Solver drm force compute",AVERAGE),
    			Timer_Value("Solver drm force compute",MAX),
    			Timer_Value("Solver drm force compute",MIN));
    }
    printf("    1st schedule send data          : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
	   Timer_Value("1st schedule send data (contribution)",AVERAGE),
	   Timer_Value("1st schedule send data (contribution)",MAX),
	   Timer_Value("1st schedule send data (contribution)",MIN));
    printf("    1st compute adjust              : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
	   Timer_Value("1st compute adjust (distribution)",AVERAGE),
	   Timer_Value("1st compute adjust (distribution)",MAX),
	   Timer_Value("1st compute adjust (distribution)",MIN));
    printf("    2nd schedule send data          : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
	   Timer_Value("2nd schedule send data (contribution)",AVERAGE),
	   Timer_Value("2nd schedule send data (contribution)",MAX),
	   Timer_Value("2nd schedule send data (contribution)",MIN));
    printf("    Compute new displacement        : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
	   Timer_Value("Compute new displacement",AVERAGE),
	   Timer_Value("Compute new displacement",MAX),
	   Timer_Value("Compute new displacement",MIN));
    printf("    3rd schedule send data          : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
	   Timer_Value("3rd schedule send data (sharing)",AVERAGE),
	   Timer_Value("3rd schedule send data (sharing)",MAX),
	   Timer_Value("3rd schedule send data (sharing)",MIN));
    printf("    2nd compute adjust              : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
	   Timer_Value("2nd compute adjust (assignment)",AVERAGE),
	   Timer_Value("2nd compute adjust (assignment)",MAX),
	   Timer_Value("2nd compute adjust (assignment)",MIN));
    printf("    4th schadule send data          : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
	   Timer_Value("4th schadule send data (sharing)",AVERAGE),
	   Timer_Value("4th schadule send data (sharing)",MAX),
	   Timer_Value("4th schadule send data (sharing)",MIN));
    printf("    IO\n");

    if ( Timer_Exists("Solver drm output") ) {
    	printf("        Drm Output                  : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
    			Timer_Value("Solver drm output",AVERAGE),
    			Timer_Value("Solver drm output",MAX),
    			Timer_Value("Solver drm output",MIN));
    }
    if ( Timer_Exists("Solver drm read displacements") ) {
    	printf("        Solver drm read disp        : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
    			Timer_Value("Solver drm read displacements",AVERAGE),
    			Timer_Value("Solver drm read displacements",MAX),
    			Timer_Value("Solver drm read displacements",MIN));
    }
    printf("        Solver Stats Print          : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
	   Timer_Value("Solver Stats Print",AVERAGE),
	   Timer_Value("Solver Stats Print",MAX),
	   Timer_Value("Solver Stats Print",MIN));
    if( Timer_Exists("Print Planes") )
	printf("        Planes                      : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
	       Timer_Value("Print Planes",AVERAGE),
	       Timer_Value("Print Planes",MAX),
	       Timer_Value("Print Planes",MIN));
    if( Timer_Exists("Print Stations") )
	printf("        Stations                    : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
	       Timer_Value("Print Stations",AVERAGE),
	       Timer_Value("Print Stations",MAX),
	       Timer_Value("Print Stations",MIN));
    if( Timer_Exists("Checkpoint") )
	printf("        Checkpoint                  : %.2f\n",
	       Timer_Value("Checkpoint",0));
    printf("\n");
    printf("TOTAL WALL CLOCK                    : %.2f seconds\n", Timer_Value("Total Wall Clock",0));
    printf("\n");
    printf("\n____________Analysis_____________\n");
    printf("Solver I/O                 : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
	   Timer_Value("Solver I/O",      AVERAGE),Timer_Value("Solver I/O",      MAX), Timer_Value("Solver I/O",      MIN));
    printf("Solver Compute             : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
	   Timer_Value("Compute Physics", AVERAGE),Timer_Value("Compute Physics", MAX), Timer_Value("Compute Physics", MIN));
    printf("Solver Communicate         : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
	   Timer_Value("Communication",   AVERAGE),Timer_Value("Communication",   MAX), Timer_Value("Communication",   MIN));
    printf("Compute/Communicate Ratio  : %.2f \n",
	   Timer_Value("Compute Physics", AVERAGE) / Timer_Value("Communication",   AVERAGE) );
    printf("\n\n\n\n");

    fflush (stdout);

    return;
}




/**
 * Prepare data to compute the myForce vector.  It calls compute_force.
 */
static void
source_init( const char* physicsin )
{
	if(( Param.drmImplement == NO )||( Param.drmImplement == YES && Param.theDrmPart == PART1 )){

		double globalDelayT = 0;
		double surfaceShift = 0;

		/* Load to Global.theMPIInformation */
		Global.theMPIInformation.myid      = Global.myID;
		Global.theMPIInformation.groupsize = Global.theGroupSize;

		/* Load to theNumericsInforamation */
		Global.theNumericsInformation.numberoftimesteps = Param.theTotalSteps;
		Global.theNumericsInformation.deltat	     = Param.theDeltaT;
		Global.theNumericsInformation.validfrequency    = Param.theFreq;

		Global.theNumericsInformation.xlength = Param.theDomainX;
		Global.theNumericsInformation.ylength = Param.theDomainY;
		Global.theNumericsInformation.zlength = Param.theDomainZ;

		if ( Param.includeNonlinearAnalysis == YES ) {
			globalDelayT = get_geostatic_total_time();
		}

		if ( Param.includeBuildings == YES ) {
			surfaceShift = get_surface_shift();
		}

		/* it will create the files to be read each time step to
       load (force) the mesh */
		if ( compute_print_source( physicsin, Global.myOctree, Global.myMesh,
				Global.theNumericsInformation, Global.theMPIInformation,
				globalDelayT, surfaceShift ) )
		{
			fprintf(stdout,"Err:cannot create source forces");
			MPI_Abort(MPI_COMM_WORLD, ERROR);
			exit(1);
		}

		Global.theNodesLoaded = source_get_local_loaded_nodes_count();

		if (Global.theNodesLoaded != 0) {
			size_t ret;
			int32_t node_count;

			Global.fpsource = source_open_forces_file( "r" );

			hu_fread( &node_count, sizeof(int32_t), 1, Global.fpsource );

			/* \todo add assertion for node_count == Global.theNodesLoaded */

			Global.theNodesLoadedList = malloc( sizeof(int32_t) * Global.theNodesLoaded );
			Global.myForces	   = calloc( Global.theNodesLoaded, sizeof(vector3D_t) );

			if (Global.myForces == NULL || Global.theNodesLoadedList == NULL) {
				solver_abort( "source_init", "memory allocation failed",
						"Cannot allocate memory for Global.myForces or "
						"loaded nodes list arrays\n" );
			}

			ret = hu_fread( Global.theNodesLoadedList, sizeof(int32_t), Global.theNodesLoaded,
					Global.fpsource );

			if (ret != Global.theNodesLoaded) {
				solver_abort( "source_init(", "fread failed",
						"Could not read nodal force file");
			}
		}
	}
}



/**
 * Search a point in the domain of the local mesh.
 *
 *   input: coordinates
 *  output: 0 fail 1 success
 */
int32_t search_point( vector3D_t point, octant_t **octant )
{
    tick_t  xTick, yTick, zTick;

    xTick = point.x[0] / Global.myMesh->ticksize;
    yTick = point.x[1] / Global.myMesh->ticksize;
    zTick = point.x[2] / Global.myMesh->ticksize;

    *octant = octor_searchoctant( Global.myOctree, xTick, yTick, zTick,
            PIXELLEVEL, AGGREGATE_SEARCH );

    if ( (*octant == NULL) || ((*octant)->where == REMOTE) ) {
        return 0;
    }

    return 1;
}

/**
 * \param octant where the point is located.
 * \param pointcoords coordinates of the point.
 * \param localcoords[out] the displacment.
 */
extern int
compute_csi_eta_dzeta( octant_t* octant, vector3D_t pointcoords,
		       vector3D_t* localcoords, int32_t* localNodeID )
{
    tick_t  edgeticks;
    int32_t eindex;
    double  center_x, center_y, center_z;

    /* various convienient variables */
    double xGlobal = pointcoords.x[0];
    double yGlobal = pointcoords.x[1];
    double zGlobal = pointcoords.x[2];
    double h;

    edgeticks = (tick_t)1 << (PIXELLEVEL - octant->level);
    h =  Global.myMesh->ticksize * edgeticks;

    /* Calculate the center coordinate of the element */
    center_x = Global.myMesh->ticksize * (octant->lx + edgeticks / 2);
    center_y = Global.myMesh->ticksize * (octant->ly + edgeticks / 2);
    center_z = Global.myMesh->ticksize * (octant->lz + edgeticks / 2);


    /* Go through my local elements to find which one matches the
     * containing octant. I should have a better solution than this.
     */
    for (eindex = 0; eindex < Global.myMesh->lenum; eindex++) {
        int32_t lnid0 = Global.myMesh->elemTable[eindex].lnid[0];

        if ((Global.myMesh->nodeTable[lnid0].x == octant->lx) &&
                (Global.myMesh->nodeTable[lnid0].y == octant->ly) &&
                (Global.myMesh->nodeTable[lnid0].z == octant->lz)) {

            /* Sanity check */
            if (Global.myMesh->elemTable[eindex].level != octant->level) {
                fprintf(stderr, "Thread %d: source_init: internal error\n",
                        Global.myID);
                MPI_Abort(MPI_COMM_WORLD, ERROR);
                exit(1);
            }

            /* Fill in the local node ids of the containing element */
            memcpy( localNodeID, Global.myMesh->elemTable[eindex].lnid,
                    sizeof(int32_t) * 8 );

            break;
        }
    }  /* for all the local elements */


    if (eindex == Global.myMesh->lenum) {
        fprintf(stderr, "Thread %d: source_init: ", Global.myID);
        fprintf(stderr, "No element matches the containing octant.\n");
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    /* Derive the local coordinate of the source inside the element */
    localcoords->x[0] =  2*(xGlobal- center_x)/h;
    localcoords->x[1] =  2*(yGlobal- center_y)/h;
    localcoords->x[2] =  2*(zGlobal- center_z)/h;

    return 1;
}


/**
 * Read stations info.  This is called by PE 0.
 */
static void
read_stations_info( const char* numericalin )
{
    static const char* fname = __FUNCTION_NAME;

    int    iStation, iCorner;
    double lon, lat, depth, *auxiliar;
    FILE*  fp;

    vector3D_t coords;

    /* obtain the stations specifications */
    if ( (fp = fopen ( numericalin, "r")) == NULL ) {
	solver_abort (fname, numericalin,
		      "Error opening numerical.in configuration file");
    }

    auxiliar = (double *)malloc(sizeof(double)*8);

    if ( parsedarray( fp, "domain_surface_corners", 8, auxiliar ) != 0) {
	solver_abort( fname, NULL,
		      "Error parsing domain_surface_corners field from %s\n",
		      numericalin);
    }

    for ( iCorner = 0; iCorner < 4; iCorner++){
	Param.theSurfaceCornersLong[ iCorner ] = auxiliar [ iCorner * 2 ];
	Param.theSurfaceCornersLat [ iCorner ] = auxiliar [ iCorner * 2 +1 ];
    }
    free(auxiliar);


    if (parsetext( fp, "output_stations_print_rate", 'i',
		   &Param.theStationsPrintRate ) != 0) {
	solver_abort( fname, NULL,
		      "Error parsing output_planes_print_rate field from %s\n",
		      numericalin );
    }

    auxiliar    = (double*)malloc( sizeof(double) * Param.theNumberOfStations * 3 );
    Param.theStationX = (double*)malloc( sizeof(double) * Param.theNumberOfStations );
    Param.theStationY = (double*)malloc( sizeof(double) * Param.theNumberOfStations );
    Param.theStationZ = (double*)malloc( sizeof(double) * Param.theNumberOfStations );

    if (Param.theStationX == NULL || Param.theStationY == NULL || Param.theStationZ == NULL) {
    	fprintf( stdout,
    			"Err alloc theStations arrays in output_stations_init" );
    	fflush( stdout );
    	MPI_Abort(MPI_COMM_WORLD, ERROR );
    	exit( 1 );
    }

    if (parsedarray( fp, "output_stations", Param.theNumberOfStations * 3,
    		auxiliar ) != 0) {
    	solver_abort (fname, NULL,
    			"Err parsing output_stations from %s\n", numericalin);
    }

    for (iStation = 0; iStation < Param.theNumberOfStations; iStation++) {
    	lat    = auxiliar [ iStation * 3 ];
    	lon    = auxiliar [ iStation * 3 +1 ];
    	depth  = auxiliar [ iStation * 3 +2 ];
    	coords = compute_domain_coords_linearinterp(lon,lat,
    			Param.theSurfaceCornersLong,
				Param.theSurfaceCornersLat,
				Param.theDomainY,Param.theDomainX);
    	Param.theStationX [ iStation ] = coords.x[0];
    	Param.theStationY [ iStation ] = coords.x[1];
    	Param.theStationZ [ iStation ] = depth;

    	if ( Param.includeBuildings == YES ) {
    		Param.theStationZ [ iStation ] += get_surface_shift();
    	}
    }

    free( auxiliar );

    if ( parsetext(fp, "output_stations_directory",'s',Param.theStationsDirOut)!= 0)
	solver_abort (fname, NULL, "Error parsing fields from %s\n",
		      numericalin);

    return;
}


/**
 * Broadcast info about the stations.
 */
void
broadcast_stations_info()
{
    /*initialize the local structures */
    if ( Global.myID != 0 ) {
	Param.theStationX = (double*)malloc( sizeof(double) * Param.theNumberOfStations);
	Param.theStationY = (double*)malloc( sizeof(double) * Param.theNumberOfStations);
	Param.theStationZ = (double*)malloc( sizeof(double) * Param.theNumberOfStations);

	if (Param.theStationX == NULL ||  Param.theStationY == NULL || Param.theStationZ==NULL) {
	    solver_abort( "broadcast_stations_info", NULL,
			  "Error: Unable to create stations arrays" );
	}
    }

    MPI_Barrier( comm_solver );

    MPI_Bcast(Param.theStationX, Param.theNumberOfStations, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(Param.theStationY, Param.theNumberOfStations, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(Param.theStationZ, Param.theNumberOfStations, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(&Param.theStationsPrintRate, 1, MPI_INT, 0, comm_solver);

    broadcast_char_array( Param.theStationsDirOut, sizeof(Param.theStationsDirOut), 0,
			  comm_solver );

    return;
}



/**
 * Prepare all the info every station needs once it is located in a processor
 */
void setup_stations_data()
{
    static char stationFile[256];

    int32_t    iStation, iLCStation = 0; /* LC local count */
    vector3D_t stationCoords;
    octant_t*  octant;

    /* look for the stations in the domain each processor has */
    for (iStation = 0; iStation < Param.theNumberOfStations; iStation++) {
        stationCoords.x[0] = Param.theStationX[iStation];
        stationCoords.x[1] = Param.theStationY[iStation];
        stationCoords.x[2] = Param.theStationZ[iStation];

        if (search_point( stationCoords, &octant ) == 1) {
            Param.myNumberOfStations++;
        }
    }

    /* allocate memory if necessary and generate the list of stations per
     * processor */
    if (Param.myNumberOfStations != 0) {
	monitor_print( "PE=%d local_stations=%d total_station_count=%d\n",
		       Global.myID, Param.myNumberOfStations, Param.theNumberOfStations );

	XMALLOC_VAR_N( Param.myStations, station_t, Param.myNumberOfStations );
    }

    for (iStation = 0; iStation < Param.theNumberOfStations; iStation++) {
	stationCoords.x[0] = Param.theStationX[iStation];
	stationCoords.x[1] = Param.theStationY[iStation];
	stationCoords.x[2] = Param.theStationZ[iStation];

	if (search_point( stationCoords, &octant ) == 1) {
	    Param.myStations[iLCStation].id = iStation;
	    Param.myStations[iLCStation].coords=stationCoords;
	    sprintf(stationFile, "%s/station.%d",Param.theStationsDirOut,iStation);
	    Param.myStations[iLCStation].fpoutputfile = hu_fopen( stationFile,"w" );
	    compute_csi_eta_dzeta( octant, Param.myStations[iLCStation].coords,
				   &(Param.myStations[iLCStation].localcoords),
				   Param.myStations[iLCStation].nodestointerpolate);

	    /*
	     * This section now enclosed in a DEBUG def because its information
	     * is only useful for debugging purposes and otherwise it just
	     * introduce noise to the post-processing.
	     */
#ifdef DEBUG
	    fprintf( Param.myStations[iLCStation].fpoutputfile,
		     "# Node identifiers:\n" );

	    int iNode;

	    for (iNode = 0; iNode < 8; iNode++) {
		fprintf( Param.myStations[iLCStation].fpoutputfile,
			 "#  %13" INT64_FMT "\n",
			 Global.myMesh->nodeTable[Param.myStations[iLCStation].nodestointerpolate[iNode]].gnid );
	    }
#endif /* DEBUG */

	    /*
	     * This is a one line heading for the station files. Spacing is
	     * such that it aligns with data.
	     *
	     * The lines after X and Y represent the actual orientation of the
	     * axes if one looks at the domain's surface on a screen. Zeta's dot
	     * means the Z axis goes in(+) and out(-) of the screen.
	     */

	    fputs( "#  Time(s)         X|(m)         Y-(m)         Z.(m)",
	            Param.myStations[iLCStation].fpoutputfile );

	    if ( ( Param.printStationVelocities    == YES ) ||
	         ( Param.printStationAccelerations == YES ) ) {
	        fputs( "       X|(m/s)       Y-(m/s)       Z.(m/s)",
	                Param.myStations[iLCStation].fpoutputfile );
	    }

            if ( Param.printStationAccelerations == YES ) {
                fputs( "      X|(m/s2)      Y-(m/s2)      Z.(m/s2)",
                        Param.myStations[iLCStation].fpoutputfile );
            }

            /*
	     * Additional headings for nonlinear data.
	     */
	    if ( Param.includeNonlinearAnalysis == YES ) {
	        fputs(  "    Epsilon_XX      Sigma_XX"
	                "    Epsilon_YY      Sigma_YY"
                        "    Epsilon_ZZ      Sigma_ZZ"
                        "    Epsilon_KK      Sigma_KK"
	                "    Epsilon_XY      Sigma_XY"
	                "    Epsilon_YZ      Sigma_YZ"
	                "    Epsilon_XZ      Sigma_XZ"
	                "        lambda            Fs"
	                "             k",
	                Param.myStations[iLCStation].fpoutputfile );
	    }

	    iLCStation += 1;
	}
    }

    free( Param.theStationX );
    free( Param.theStationY );
    free( Param.theStationZ );
}


/**
 * Interpolate the displacements for the stations.
 */
static int
interpolate_station_displacements( int32_t step )
{
    int iPhi;

    /* Auxiliary array to handle shape functions in a loop */
    double  xi[3][8]={ {-1,  1, -1,  1, -1,  1, -1, 1} ,
                       {-1, -1,  1,  1, -1, -1,  1, 1} ,
                       {-1, -1, -1, -1,  1,  1,  1, 1} };

    double     phi[8];
    double     dis_x, dis_y, dis_z;
    double     vel_x, vel_y, vel_z;
    double     acc_x, acc_y, acc_z;
    int32_t    iStation,nodesToInterpolate[8];;
    vector3D_t localCoords; /* convenient renaming */

    for (iStation=0;iStation<Param.myNumberOfStations; iStation++) {

        localCoords = Param.myStations[iStation].localcoords;

        for (iPhi=0; iPhi<8; iPhi++) {
            nodesToInterpolate[iPhi]
                    = Param.myStations[iStation].nodestointerpolate[iPhi];
        }

        /* Compute interpolation function (phi) for each node and
         * load the displacements
         */
        dis_x = 0;
        dis_y = 0;
        dis_z = 0;

        for (iPhi = 0; iPhi < 8; iPhi++) {
            phi[ iPhi ] = ( 1 + xi[0][iPhi]*localCoords.x[0] )
		                * ( 1 + xi[1][iPhi]*localCoords.x[1] )
		                * ( 1 + xi[2][iPhi]*localCoords.x[2] ) / 8;

            dis_x += phi[iPhi] * Global.mySolver->tm1[ nodesToInterpolate[iPhi] ].f[0];
            dis_y += phi[iPhi] * Global.mySolver->tm1[ nodesToInterpolate[iPhi] ].f[1];
            dis_z += phi[iPhi] * Global.mySolver->tm1[ nodesToInterpolate[iPhi] ].f[2];
        }

        double time = Param.theDeltaT * step;

        /*
         * Please DO NOT CHANGE the format for printing the displacements.
         * It has to be *this* one because it goes in hand with the printing
         * format for the nonlinear information.
         */
        fprintf( Param.myStations[iStation].fpoutputfile,
                 "\n%10.6f % 8e % 8e % 8e",
                 time, dis_x, dis_y, dis_z );

        /*
         * Addition for printing velocities on the fly
         */

        if ( ( Param.printStationVelocities    == YES ) ||
             ( Param.printStationAccelerations == YES ) ) {

            for (iPhi = 0; iPhi < 8; iPhi++) {

                phi[ iPhi ] = ( 1 + xi[0][iPhi]*localCoords.x[0] )
                                    * ( 1 + xi[1][iPhi]*localCoords.x[1] )
                                    * ( 1 + xi[2][iPhi]*localCoords.x[2] ) / 8;

                dis_x -= phi[iPhi] * Global.mySolver->tm2[ nodesToInterpolate[iPhi] ].f[0];
                dis_y -= phi[iPhi] * Global.mySolver->tm2[ nodesToInterpolate[iPhi] ].f[1];
                dis_z -= phi[iPhi] * Global.mySolver->tm2[ nodesToInterpolate[iPhi] ].f[2];
            }

            vel_x = dis_x / Param.theDeltaT;
            vel_y = dis_y / Param.theDeltaT;
            vel_z = dis_z / Param.theDeltaT;

            fprintf( Param.myStations[iStation].fpoutputfile,
                     " % 8e % 8e % 8e", vel_x, vel_y, vel_z );
        }

        /*
         * Addition for printing accelerations on the fly
         */

        if ( Param.printStationAccelerations == YES ) {

            for (iPhi = 0; iPhi < 8; iPhi++) {

                phi[ iPhi ] = ( 1 + xi[0][iPhi]*localCoords.x[0] )
                                            * ( 1 + xi[1][iPhi]*localCoords.x[1] )
                                            * ( 1 + xi[2][iPhi]*localCoords.x[2] ) / 8;

                dis_x -= phi[iPhi] * Global.mySolver->tm2[ nodesToInterpolate[iPhi] ].f[0];
                dis_y -= phi[iPhi] * Global.mySolver->tm2[ nodesToInterpolate[iPhi] ].f[1];
                dis_z -= phi[iPhi] * Global.mySolver->tm2[ nodesToInterpolate[iPhi] ].f[2];

                dis_x += phi[iPhi] * Global.mySolver->tm3[ nodesToInterpolate[iPhi] ].f[0];
                dis_y += phi[iPhi] * Global.mySolver->tm3[ nodesToInterpolate[iPhi] ].f[1];
                dis_z += phi[iPhi] * Global.mySolver->tm3[ nodesToInterpolate[iPhi] ].f[2];
            }

            acc_x = dis_x / Param.theDeltaTSquared;
            acc_y = dis_y / Param.theDeltaTSquared;
            acc_z = dis_z / Param.theDeltaTSquared;

            fprintf( Param.myStations[iStation].fpoutputfile,
                     " % 8e % 8e % 8e", acc_x, acc_y, acc_z );
        }

	/* TODO: Have this 10 as a parameter with a default value */
        if ( (step % (Param.theStationsPrintRate*10)) == 0 ) {
            fflush(Param.myStations[iStation].fpoutputfile);
        }
    }

    return 1;
}


/**
 * Init stations info and data structures
 */
void output_stations_init( const char* numericalin )
{
    if (Global.myID == 0) {
	read_stations_info( numericalin );
    }

    broadcast_stations_info();
    setup_stations_data();

    MPI_Barrier( comm_solver );

    return;
}



/**
 * \note This function should only be called by PE with rank 0.
 */
static int
load_output_parameters (const char* numericalin, output_parameters_t* params)
{
    FILE* fp;
    int   ret, value;
    char  filename[LINESIZE];

    assert (NULL != numericalin);
    assert (NULL != params);
    assert (Global.myID == 0);

     /* Read output parameters from numerical.in */
    fp = fopen (numericalin, "r");

    if (NULL == fp) {
	solver_abort ("load_output_parameters", "fopen", "numerical.in=\"%s\"",
		      numericalin);
    }

    params->do_output		= 0;
    params->parallel_output     = 0;
    params->output_displacement = 0;
    params->output_velocity     = 0;
    params->output_debug	= 0;

    params->displacement_filename = NULL;
    params->velocity_filename     = NULL;


    /* read parameters from configuration file */

    value = 0;
    ret = parsetext (fp, "output_parallel", 'i', &value);

    if (0 == ret && 0 != value) {
	value = 0;
	ret = parsetext (fp, "output_displacement", 'i', &value);

	if (0 == ret && 0 != value) { /* output_displacement = 1 in config */
	    ret = read_config_string (fp, "output_displacement_file",
				      filename, LINESIZE);

	    if (1 == ret && filename[0] != '\0') {
		params->displacement_filename = strdup (filename);
		params->output_displacement = 1;
	    } else {
		solver_abort ("load_output_parameters", NULL,
			      "Output displacement file name not specified in "
			      "numerical.in=\"%s\"",
			      numericalin);
	    }
	}

	value = 0;
	ret   = parsetext (fp, "output_velocity", 'i', &value);

	if (0 == ret && 0 != value) { /* output_displacement = 1 in config */
	    ret = read_config_string (fp, "output_velocity_file",
				      filename, LINESIZE);

	    if (1 == ret && filename[0] != '\0') {
		params->velocity_filename = strdup (filename);
		params->output_velocity = 1;
	    } else {
		solver_abort ("load_output_parameters", NULL,
			      "Output velocity file name not specified in "
			      "numerical.in=\"%s\"",
			      numericalin);
	    }
	}

	params->stats_filename = "output-stats.txt"; /* default value */

	ret = read_config_string (fp, "output_stats_file", filename, LINESIZE);

	if (1 == ret && filename[0] != '\0') {
	    params->stats_filename = strdup (filename);
	}

	params->debug_filename = "output-debug.txt"; /* default value */
	ret = read_config_string (fp, "output_debug_file", filename, LINESIZE);

	if (1 == ret && filename[0] != '\0') {
	    params->debug_filename = strdup (filename);
	}


	if (params->output_velocity || params->output_displacement) {
	    params->do_output = 1;
	}

	params->parallel_output = 1;

	ret = parsetext (fp, "output_debug", 'i', &value);
	if (0 == ret && 0 != value) { /* output_debug = 1 in config */
	    params->output_debug = 1;
	}
    }

    ret = 0;

    fclose (fp);

    return ret;
}



/**
 * Initialize output structures, including the opening of 4D output files.
 *
 * \param numericsin Name of the file with the solver and output parameters
 *	i.e., "numerical.in".
 * \param[out] params output parameters, including filenames.
 *
 * \pre The following global variables should be initialized:
 *	- Global.myID.
 *	- Global.theGroupSize.
 *
 * \post On a successful return, the output argument \c params will be
 *	properly initialized.  If the routine fails, the state of the
 *	\c param struct is undefined.
 *
 * \return 0 on success, -1 on error.
 */
static int
output_init_parameters (const char* numericalin, output_parameters_t* params)
{
    /* jc: this ugly #define here is because the catamount compiler does not
     * like static const  */
#define VALUES_COUNT    4
    int ret;
    int32_t values[VALUES_COUNT];
    off_t output_steps;

    /* sanity cleanup */
    memset (params, 0, sizeof (output_parameters_t));

    params->do_output		  = 0;
    params->displacement_filename = NULL;
    params->velocity_filename     = NULL;
    params->stats_filename	  = NULL;

    if (Global.myID == 0) {
	ret = load_output_parameters (numericalin, params);

	if (0 != ret) {
	    solver_abort ("output_init_parameters", NULL, NULL);
	    return -1;
	}
    }

    /* parameters that can be initialized from global variables */
    params->pe_id	   = Global.myID;
    params->pe_count       = Global.theGroupSize;
    params->total_nodes    = Global.theNTotal;
    params->total_elements = Global.theETotal;
    params->output_rate	   = Param.theRate;
    params->domain_x	   = Param.theDomainX;
    params->domain_y	   = Param.theDomainY;
    params->domain_z	   = Param.theDomainZ;
    params->mesh	   = Global.myMesh;
    params->solver	   = (solver_t*)Global.mySolver;
    params->delta_t	   = Param.theDeltaT;
    params->total_time_steps = Param.theTotalSteps;

    output_steps	   = (Param.theTotalSteps - 1) / Param.theRate + 1;


    values[0] = params->parallel_output;
    values[1] = params->output_displacement;
    values[2] = params->output_velocity;
    values[3] = params->output_debug;


    MPI_Bcast (values, VALUES_COUNT, MPI_INT, 0, comm_solver);
    MPI_Bcast (&params->total_nodes, 1, MPI_INT64, 0, comm_solver);

    params->parallel_output     = values[0];
    params->output_displacement = values[1];
    params->output_velocity	= values[2];
    params->output_debug	= values[3];
    params->output_size		= (output_steps * params->total_nodes
				   * sizeof(fvector_t) + sizeof(out_hdr_t));
    Global.theNTotal = params->total_nodes;


    if (params->parallel_output) {
	if (params->output_displacement) {
	    broadcast_string(&params->displacement_filename, 0,comm_solver);
	}

	if (params->output_velocity) {
	    broadcast_string( &params->velocity_filename, 0, comm_solver );
	}

	if (params->output_debug) {
	    broadcast_string( &params->debug_filename, 0, comm_solver );
	}

	/* set the expected file size */
	Param.the4DOutSize = params->output_size;
    }

    return 0;
#undef VALUES_COUNT
}


/**
 * Intialize parallel output structures according to config file
 * parameters.
 *
 * \note params parameter could be a local var instead, it doesn't need
 * to be global, since this is not really used after the initialization
 * except for the stats file name.
 *
 * \return 0 on success, -1 on error (probably aborts instead).
 */
static int
output_init( const char* numericalin, output_parameters_t* params )
{
    int ret = -1;

    assert (NULL != numericalin);
    assert (NULL != params);

    ret = output_init_parameters( numericalin, params );

    if (ret != 0) {
	return -1;
    }

    /* initialize structures, open files, etc */
    ret = output_init_state( params );

    /* these aren't needed anymore after this point */
    xfree_char( &params->displacement_filename );
    xfree_char( &params->velocity_filename );

    return ret;
}


static int
output_get_stats( void )
{
    output_stats_t disp_stats, vel_stats;
    int    ret      = 0;
    double avg_tput = 0;


    if (Param.theOutputParameters.parallel_output) {
	ret = output_collect_io_stats( Param.theOutputParameters.stats_filename,
				       &disp_stats, &vel_stats, Timer_Value("Total Wall Clock",0) );

	/* magic trick so print_timing_stat() prints something sensible
	 * for the 4D output time in the parallel case.
	 *
	 * if both displacement and velocity where written out, prefer
	 * the displacement stat.
	 */
	if (Param.theOutputParameters.output_velocity) {
	    avg_tput = vel_stats.tput_avg;
	}

	if (Param.theOutputParameters.output_displacement) {
	    avg_tput = disp_stats.tput_avg;
	}
    }

    return ret;
}



/**
 * Adjust the values of the properties (Vp, Vs, Rho) for the mesh elements.
 * Initial implementation (by Leo) querying the middle point.1D.
 * Modified implementation (Ricardo) querying the 27 points and averaging.
 *
 * \param cvm the material model database (CVM etree).
 */
static void
mesh_correct_properties( etree_t* cvm )
{
    elem_t*  elemp;
    edata_t* edata;
    int32_t  eindex;
    double   east_m, north_m, depth_m, VpVsRatio, RhoVpRatio;
    int	     res, iNorth, iEast, iDepth, numPoints = 3;
    double   vs, vp, rho;
    double   points[3];
    int32_t  lnid0;

    // INTRODUCE BKT MODEL

    double Qs, Qp, Qk, L, vs_vp_Ratio, vksquared, w;
    int index_Qs, index_Qk;
    int QTable_Size = (int)(sizeof(Global.theQTABLE)/( 6 * sizeof(double)));
    int VModel_Size = (int)(sizeof(Global.the1DVModel)/( 6 * sizeof(double)));

    points[0] = 0.005;
    points[1] = 0.5;
    points[2] = 0.995;

//    if (Global.myID == 0) {
//        fprintf( stdout,"mesh_correct_properties  ... " );
//        fflush( stdout );
//    }

    /* iterate over mesh elements */
    for (eindex = 0; eindex < Global.myMesh->lenum; eindex++) {

        elemp = &Global.myMesh->elemTable[eindex];
        edata = (edata_t*)elemp->data;
        lnid0 = Global.myMesh->elemTable[eindex].lnid[0];

        if ( Param.includeBuildings == YES ) {
            if( bldgs_correctproperties( Global.myMesh, edata, lnid0) ) {
                continue;
            }
        }

        vp  = 0;
        vs  = 0;
        rho = 0;

        for (iNorth = 0; iNorth < numPoints; iNorth++) {

        	north_m = (Global.myMesh->ticksize) * (double)Global.myMesh->nodeTable[lnid0].x
        			+ edata->edgesize * points[iNorth] + Global.theXForMeshOrigin ;

        	for (iEast = 0; iEast < numPoints; iEast++) {

        		east_m = ( (Global.myMesh->ticksize)
        				* (double)Global.myMesh->nodeTable[lnid0].y
        				+ edata->edgesize * points[iEast] + Global.theYForMeshOrigin  );

        		for (iDepth = 0; iDepth < numPoints; iDepth++) {
        			cvmpayload_t g_props; /* ground properties record */

        			depth_m = ( (Global.myMesh->ticksize)
        					* (double)Global.myMesh->nodeTable[lnid0].z
        					+ edata->edgesize * points[iDepth] + Global.theZForMeshOrigin );

        			/* NOTE: If you want to see the carving process,
        			 *       activate this and comment the query below */
        			if ( Param.includeBuildings == YES ) {
        				//                        if ( depth_m < get_surface_shift() ) {
        				//                            g_props.Vp  = NAN;
        				//                            g_props.Vs  = NAN;
        				//                            g_props.rho = NAN;
        				//                        } else {
        				depth_m -= get_surface_shift();
        				//                            res = cvm_query( Global.theCVMEp, east_m, north_m,
        				//                                             depth_m, &g_props );
        				//                        }
        				//
        			}

//        			res = cvm_query( Global.theCVMEp, east_m, north_m, depth_m, &g_props );
        			res = VModel_query(depth_m, &g_props);

        			if (res != 0) {
        				fprintf(stderr, "Cannot find the query point\n");
        				exit(1);
        			}

        			vp  += g_props.Vp;
        			vs  += g_props.Vs;
        			rho += g_props.rho;
                }
            }
        }

        edata->Vp  =  vp / 27;
        edata->Vs  =  vs / 27;
        edata->rho = rho / 27;

        /* Auxiliary ratios for adjustments */
        VpVsRatio  = edata->Vp  / edata->Vs;
        RhoVpRatio = edata->rho / edata->Vp;

        /* Adjust material properties according to the element size and
         * softening factor.
         *
         * A factor of 1 means perfect compliance between the mesh and the
         * elements' material properties resulting in strong changes to the
         * results. A factor of 4 tends to double the simulation delta_t
         * without affecting too much the results. Testing is needed but for
         * now I recommend factors > 4.
         */
        if ( Param.theSofteningFactor > 0 ) {

            double idealVs, factoredVs;

            idealVs    = edata->edgesize * Param.theFactor;
            factoredVs = idealVs * Param.theSofteningFactor;

            if ( edata->Vs > factoredVs ) {
                edata->Vs  = factoredVs;
                edata->Vp  = factoredVs * VpVsRatio;
                edata->rho = edata->Vp  * RhoVpRatio;
            }
        }

        /* Readjust Vs, Vp and Density according to VsCut */
        if ( edata->Vs < Param.theVsCut ) {
            edata->Vs  = Param.theVsCut;
            edata->Vp  = Param.theVsCut  * VpVsRatio;
            /* edata->rho = edata->Vp * RhoVpRatio; */ /* Discuss with Jacobo */
        }


        // IMPLEMENT BKT MODEL

        /* CALCULATE QUALITY FACTOR VALUES AND READ CORRESPONDING VISCOELASTICITY COEFFICIENTS FROM THE TABLE */

        	/* L IS THE COEFFICIENT DEFINED BY "SHEARER-2009" TO RELATE QK, QS AND QP */

        if(Param.theTypeOfDamping == BKT)
        {

            vksquared = edata->Vp * edata->Vp - 4. / 3. * edata->Vs * edata->Vs;
        	vs_vp_Ratio = edata->Vs / edata->Vp;
        	vs = edata->Vs * 0.001;
        	L = 4. / 3. * vs_vp_Ratio * vs_vp_Ratio;

          	//Qs = 0.02 * edata->Vs;

        	// Ricardo's Formula based on Brocher's paper (2008) on the subject. In the paper Qp = 2*Qs is given.
        	//TODO : Make sure Qp Qs relation is correct...

        	Qs = 10.5 + vs * (-16. + vs * (153. + vs * (-103. + vs * (34.7 + vs * (-5.29 + vs * 0.31)))));
        	Qp = 2. * Qs;

        	if (Param.useInfQk == YES) {
        	    Qk = 1000;
        	} else {
                Qk = (1. - L) / (1. / Qp - L / Qs);
        	}

        	index_Qs = Search_Quality_Table(Qs, &(Global.theQTABLE[0][0]), QTable_Size);

//        	printf("Quality Factor Table\n Qs : %lf \n Vs : %lf\n",Qs,edata->Vs);

        	if(index_Qs == -2 || index_Qs >= QTable_Size)
        	{
        		fprintf(stderr,"Problem with the Quality Factor Table\n Qs : %lf \n Vs : %lf\n",Qs,edata->Vs);
        		exit(1);
        	}
        	else if(index_Qs == -1)
        	{
        		edata->a0_shear = 0;
        		edata->a1_shear = 0;
        		edata->g0_shear = 0;
        		edata->g1_shear = 0;
        		edata->b_shear  = 0;
        	}
        	else
        	{
        		edata->a0_shear = Global.theQTABLE[index_Qs][1];
        		edata->a1_shear = Global.theQTABLE[index_Qs][2];
        		edata->g0_shear = Global.theQTABLE[index_Qs][3];
        		edata->g1_shear = Global.theQTABLE[index_Qs][4];
        		edata->b_shear  = Global.theQTABLE[index_Qs][5];
        	}

        	index_Qk = Search_Quality_Table(Qk, &(Global.theQTABLE[0][0]), QTable_Size);

//        	printf("Quality Factor Table\n Qs : %lf \n Vs : %lf\n",Qs,edata->Vs);

        	if(index_Qk == -2 || index_Qk >= QTable_Size)
        	{
        		fprintf(stderr,"Problem with the Quality Factor Table\n Qk : %lf \n Vs : %lf\n",Qk,edata->Vs);
        		exit(1);
        	}
        	else if(index_Qk == -1)
        	{
        		edata->a0_kappa = 0;
        		edata->a1_kappa = 0;
        		edata->g0_kappa = 0;
        		edata->g1_kappa = 0;
        		edata->b_kappa  = 0;
        	}
        	else
        	{
        		edata->a0_kappa = Global.theQTABLE[index_Qk][1];
        		edata->a1_kappa = Global.theQTABLE[index_Qk][2];
        		edata->g0_kappa = Global.theQTABLE[index_Qk][3];
        		edata->g1_kappa = Global.theQTABLE[index_Qk][4];
        		edata->b_kappa  = Global.theQTABLE[index_Qk][5];
        	}

        	if(Param.theFreq_Vel != 0.)
        	{
        		w = Param.theFreq_Vel / Param.theFreq;

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
    }
}


/*** Program's standard entry point. ***/
int main( int argc, char** argv )
{

#ifdef DEBUG
    int32_t flag;
    MPI_Status status;
#endif /* DEBUG */

    /* MPI initialization */
    MPI_Init(&argc, &argv);
    Timer_Start("Total Wall Clock");
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm_rank(MPI_COMM_WORLD, &Global.myID);
    MPI_Comm_size(MPI_COMM_WORLD, &Global.theGroupSize);

    /* Make sure using correct input arguments */
    if (argc != 2) {
        if (Global.myID == 0) {
            fputs ( "Usage: psolve <parameters.in>\n", stderr);
        }
        MPI_Finalize();
        exit(1);
    }

    /*Read in and verify IO Pool Configuration */
    char * IO_PES_ENV_VAR;
    IO_PES_ENV_VAR = getenv("IO_PES");
    if (IO_PES_ENV_VAR==NULL)
        Param.IO_pool_pe_count = 0;
    else
        Param.IO_pool_pe_count = atoi(IO_PES_ENV_VAR);
    if (Param.IO_pool_pe_count >= Global.theGroupSize)
    {
        Param.IO_pool_pe_count = 0;
        if (Global.myID==0)
            printf("Warning: IO_PES too large.  Set to 0.\n");
    }

    /* Split off PEs into IO Pool */
    MPI_Comm_dup(MPI_COMM_WORLD, &comm_IO);
    int in_io_pool, color;
    Global.theGroupSize -= Param.IO_pool_pe_count;
    if (Global.myID >= Global.theGroupSize)
        in_io_pool = 1;
    else
        in_io_pool=0;
    if (in_io_pool)
        color = MPI_UNDEFINED;
    else
        color = 0;
    MPI_Comm_split(MPI_COMM_WORLD, color, Global.myID, &comm_solver);
    if (in_io_pool) {
        planes_IO_PES_main(Global.myID);
        goto IO_PES_REJOIN;
    }

    if (Global.myID == 0) {
        printf( "Starting psolve $Revision: 1.166 $ on %d PEs (%d IO Pool PEs).\n\n",
                Global.theGroupSize, Param.IO_pool_pe_count );
        fflush( stdout );
    }

    /* Read input parameters from file */
    read_parameters(argc, argv);

    /* Create and open database */
    open_cvmdb();

    /* Initialize nonlinear parameters */
    if ( Param.includeNonlinearAnalysis == YES ) {
        nonlinear_init(Global.myID, Param.parameters_input_file, Param.theDeltaT, Param.theEndT);
    }

    if ( Param.includeBuildings == YES ){
        bldgs_init( Global.myID, Param.parameters_input_file );
    }

    if ( Param.drmImplement == YES ){
    	Timer_Start("Init Drm Parameters");
    	Param.theDrmPart = drm_init(Global.myID, Param.parameters_input_file , Param.includeBuildings);
    	Timer_Stop("Init Drm Parameters");
    	Timer_Reduce("Init Drm Parameters", MAX | MIN | AVERAGE , comm_solver);
    }

    // INTRODUCE BKT MODEL
    /* Init Quality Factor Table */
    constract_Quality_Factor_Table();

    /* Generate, partition and output unstructured octree mesh */
    mesh_generate();

    if ( Param.includeBuildings == YES ){
        if ( get_fixedbase_flag() == YES ) {
            bldgs_fixedbase_init( Global.myMesh, Param.theEndT-Param.theStartT );
        }
        bldgs_finalize();
    }

    if ( Param.drmImplement == YES ) {
    	Timer_Start("Drm Init");
    	if ( Param.theDrmPart == PART0 ) {
    		find_drm_nodes(Global.myMesh, Global.myID, Param.parameters_input_file,
    				Global.myOctree->ticksize, Global.theGroupSize);
    	}
    	if (Param.theDrmPart == PART1) {
    		setup_drm_data(Global.myMesh, Global.myID, Global.theGroupSize);
    	}
    	if (Param.theDrmPart == PART2) {
    		proc_drm_elems(Global.myMesh, Global.myID, Global.theGroupSize, Param.theTotalSteps);
    	}
    	drm_stats(Global.myID, Global.theGroupSize, Global.theXForMeshOrigin,
    			 Global.theYForMeshOrigin, Global.theZForMeshOrigin);
    	Timer_Stop("Drm Init");
    	Timer_Reduce("Drm Init", MAX | MIN | AVERAGE , comm_solver);
    }

    if (Param.theMeshOutFlag && DO_OUTPUT) {
        mesh_output();
    }

    if ( Param.storeMeshCoordinatesForMatlab == YES ) {
        saveMeshCoordinatesForMatlab( Global.myMesh, Global.myID, Param.parameters_input_file,
				      Global.myOctree->ticksize,Param.theTypeOfDamping,Global.theXForMeshOrigin,
				      Global.theYForMeshOrigin,Global.theZForMeshOrigin, Param.includeBuildings);
    }

    Timer_Start("Mesh Stats Print");
    mesh_print_stat(Global.myOctree, Global.myMesh, Global.myID, Global.theGroupSize,
		    Param.theMeshStatFilename);
    Timer_Stop("Mesh Stats Print");
    Timer_Reduce("Mesh Stats Print", MAX | MIN, comm_solver);

    /* Initialize the output planes */
    if ( Param.theNumberOfPlanes != 0 ) {
        planes_setup(Global.myID, &Param.thePlanePrintRate, Param.IO_pool_pe_count,
		     Param.theNumberOfPlanes, Param.parameters_input_file, get_surface_shift(),
		     Param.theSurfaceCornersLong, Param.theSurfaceCornersLat,
		     Param.theDomainX, Param.theDomainY, Param.theDomainZ,
		     Param.planes_input_file);
    }

    if ( Param.theNumberOfStations !=0 ){
        output_stations_init(Param.parameters_input_file);
    }

    /* Initialize the solver, source and output structures */
    solver_init();
    Timer_Start("Solver Stats Print");
    solver_printstat( Global.mySolver );
    Timer_Stop("Solver Stats Print");
    Timer_Reduce("Solver Stats Print", MAX | MIN, comm_solver);

    /* Initialize nonlinear solver analysis structures */
    if ( Param.includeNonlinearAnalysis == YES ) {
        nonlinear_solver_init(Global.myID, Global.myMesh, Param.theDomainZ);
        if ( Param.theNumberOfStations !=0 ){
            nonlinear_stations_init(Global.myMesh, Param.myStations, Param.myNumberOfStations);
        }
        nonlinear_stats(Global.myID, Global.theGroupSize);
    }

    Timer_Start("Source Init");
    source_init(Param.parameters_input_file);
    Timer_Stop("Source Init");
    Timer_Reduce("Source Init", MAX | MIN, comm_solver);

    /* Mapping element indices for stiffness
     * This is for compatibility with nonlinear
     * \TODO a more clever way should be possible
     */
    stiffness_init(Global.myID, Global.myMesh);

    /* this is a little too late to check for output parameters,
     * but let's do this in the mean time
     */
    output_init (Param.parameters_input_file, &Param.theOutputParameters);

    /* Run the solver and output the results */
    MPI_Barrier(comm_solver);
    Timer_Start("Solver");
    solver_run();
    Timer_Stop("Solver");
    MPI_Barrier(comm_solver);

    if ( Param.includeNonlinearAnalysis == YES ) {
        nonlinear_yield_stats( Global.myMesh, Global.myID, Param.theTotalSteps,
			       Global.theGroupSize );
    }

    output_fini();

#ifdef DEBUG
    /* Does the OS page out my resident set ? */
    if ((Global.myID % PROCPERNODE) == 0) {
        /* system("ps xl"); */
    }

    /* Are there pending messages that I haven't processed */
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm_solver, &flag, &status);
    if (flag != 0) {
        fprintf(stderr, "Thread %d: MPI error: unreceived incoming message\n",
                Global.myID);
    }
#endif /* DEBUG */

    local_finalize();

    output_get_stats();

    MPI_Barrier(comm_solver);
    Timer_Stop("Total Wall Clock");

    /* Print out the timing stat */
    if (Global.myID == 0) {
        print_timing_stat();
    }

    /* TODO: Think of a better place for this */
    if ( Param.includeNonlinearAnalysis == YES ) {
        if ( get_geostatic_total_time() > 0 ) {
            check_balance(Global.myID);
        }
    }

    /* Send a message to IO pool PEs to close output files and goto here */
    if (Param.theNumberOfPlanes != 0) {
        planes_close(Global.myID, Param.IO_pool_pe_count, Param.theNumberOfPlanes);
    }
    IO_PES_REJOIN:

    MPI_Finalize();

    return 0;
}
