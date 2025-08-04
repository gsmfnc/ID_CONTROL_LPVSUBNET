/*  File    : gyro_simulation.cpp
 *  Abstract:
 *
 *      mex file S-function for simulation and visualisation of the Gyroscope setup.
 *
 *      Implements Nonlinear state-space equations for combinations of locked frames of the form:
 *      Ex' = A,  where :
 *          E is a Matrix containing nonlinear elements, with size between 0x0 and 4x4 depending on locking;
 *          A is a Vector containing nonlinear elements, with size between 0x1 and 4x1 depending on locking;
 *      The function solves x' through LU factorization with partial pivoting.
 *
 *      Implements VRML angle calculations through quaternions.
 *
 *      Inputs:     <name>      <unit>  <size>
 *                  pos_disk    [rad]   [1x1]
 *                  pos_bg      [rad]   [1x1]
 *                  pos_rg      [rad]   [1x1]
 *                  pos_sg      [rad]   [1x1]
 *
 *      Outputs:    <name>      <unit>  <size>
 *                  States      rad/(s) [8x1]   where, States[1:4] are positions, States[5:8] are velocities
 *                  VRML        -       [16x1]  where, VRMLq1 = [1:4,1]; VRMLq2 = [5:8,1]; etc.
 *
 *      Ver: 1.0    15-05-2015
 *      Ver: 1.3    03-12-2019      Added inertias as parameter
 *      Ver. 1.4    05-12-2019      Added motor constants as parameter
 */

#define S_FUNCTION_NAME gyro_simulation
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include "math.h"
#include "Eigen/Dense"

using namespace Eigen;

// Initial condition vector
#define X0_IDX 0
#define X0_PARAM(S) ssGetSFcnParam(S,X0_IDX)

// Viscous friction vector
#define f_IDX 1
#define f_PARAM(S) ssGetSFcnParam(S,f_IDX)

// disk lock
#define w1l_IDX 2
#define w1l_PARAM(S) ssGetSFcnParam(S,w1l_IDX)

// blue gimbal lock
#define q2l_IDX 3
#define q2l_PARAM(S) ssGetSFcnParam(S,q2l_IDX)

// red gimbal lock
#define q3l_IDX 4
#define q3l_PARAM(S) ssGetSFcnParam(S,q3l_IDX)

// silver gimbal lock
#define q4l_IDX 5
#define q4l_PARAM(S) ssGetSFcnParam(S,q4l_IDX)

// Inertia
#define inertia_IDX 6
#define inertia_PARAM(S) ssGetSFcnParam(S,inertia_IDX)

// Motor constants
#define Km_IDX 7
#define Km_PARAM(S) ssGetSFcnParam(S,Km_IDX)

#define NPARAMS 8       // Define the number of parameters

#define NSTATES 8       // Define the number of state variables
#define NINPUTS 4       // Define the number of input arguments
#define NOUTPUTS 2      // Define the number of output arguments
#define NINERTIAS 9     // Define the number of inertia values

#define C_PI 3.14159265358979323846 /* pi */

#define OK_EMPTY_DOUBLE_PARAM(pVal) (mxIsNumeric(pVal) && !mxIsLogical(pVal) &&\
!mxIsSparse(pVal) && !mxIsComplex(pVal) && mxIsDouble(pVal))

/* Function Declarations */
void VRML2Quaternion(double axis[4], double angle, double quaternion[4]);
void multiplicateQuaternions(double q1[4], double q2[4], double result[4]);
void Quaternion2VRML(double quaternion[4], double vrml[4]);
void gyro_eq(int lockCombination, const double parameters[13], double x[8], double u[4], double Edss[16], double A[4], double dxdt[8]);


/*====================*
 * S-function methods *
 *====================*/

#define MDL_CHECK_PARAMETERS
#if defined(MDL_CHECK_PARAMETERS) && defined(MATLAB_MEX_FILE)
/* Function: mdlCheckParameters =============================================
 * Abstract:
 *    Validate our parameters to verify they are okay.
 */
static void mdlCheckParameters(SimStruct *S)
{
    /* Check 1st parameter: X0 */
    {
        if ( ((mxGetM(X0_PARAM(S)) != 0) &&
              (mxGetM(X0_PARAM(S)) != NSTATES)) || !OK_EMPTY_DOUBLE_PARAM(X0_PARAM(S)) ) {
            ssSetErrorStatus(S,"1st parameter to S-function "
                             "\"X0-Matrix\" is not dimensioned "
                             "correctly");
            return;
        }
    }
    
    /* Check 2nd parameter: f (friction) */
    {
        if ( ((mxGetM(f_PARAM(S)) != 0) &&
              (mxGetM(f_PARAM(S)) != NINPUTS)) || !OK_EMPTY_DOUBLE_PARAM(f_PARAM(S)) ) {
            ssSetErrorStatus(S,"2nd parameter to S-function "
                             "\"f-Matrix\" is not dimensioned "
                             "correctly");
            return;
        }
    }

    /* Check 7th parameter: I (inertia) */
    {
        if ( ((mxGetM(inertia_PARAM(S)) != 0) &&
              (mxGetM(inertia_PARAM(S)) != NINERTIAS)) || !OK_EMPTY_DOUBLE_PARAM(inertia_PARAM(S)) ) {
            ssSetErrorStatus(S,"7th parameter to S-function "
                             "\"inertia-matrix\" is not dimensioned "
                             "correctly");
            return;
        }
    }

    /* Check 8th parameter: Km (motor constant) */
    {
        if ( ((mxGetM(Km_PARAM(S)) != 0) &&
              (mxGetM(Km_PARAM(S)) != NINPUTS)) || !OK_EMPTY_DOUBLE_PARAM(Km_PARAM(S)) ) {
            ssSetErrorStatus(S,"8th parameter to S-function "
                             "\"Km-matrix\" is not dimensioned "
                             "correctly");
            return;
        }
    }
}
#endif /* MDL_CHECK_PARAMETERS */



/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, NPARAMS);  /* Number of expected parameters */
#if defined(MATLAB_MEX_FILE)
    if (ssGetNumSFcnParams(S) == ssGetSFcnParamsCount(S)) {
        mdlCheckParameters(S);
        if (ssGetErrorStatus(S) != NULL) {
            return;
        }
    } else {
        return; /* Parameter mismatch will be reported by Simulink */
    }
#endif
    
    {
        int iParam = 0;
        int nParam = ssGetNumSFcnParams(S);
        
        for ( iParam = 0; iParam < nParam; iParam++ )
        {
            ssSetSFcnParamTunable( S, iParam, SS_PRM_SIM_ONLY_TUNABLE );
        }
    }
    
    ssSetNumContStates(S, (int_T)NSTATES);
    ssSetNumDiscStates(S, 0);
    
    /* Input-port configuration */
    if (!ssSetNumInputPorts(S, (int_T)NINPUTS)) return;     // # of inputs
    for (int i = 0; i<(int_T)NINPUTS; i++) {
        ssSetInputPortWidth(S, i, 1);
        ssSetInputPortDirectFeedThrough(S, i, 1);
    }
    
    /* Output-port configuration */
    if (!ssSetNumOutputPorts(S, (int_T)NOUTPUTS)) return;   // # of outputs
    ssSetOutputPortWidth(S, 0, 8);
    ssSetOutputPortWidth(S, 1, 16);
    
    ssSetNumSampleTimes(S, 1);      // Number of sample times
    
    /*
     * Set size of the work vectors.
     */
    ssSetNumRWork(S, 0);            // Number of real work vector elements
    ssSetNumIWork(S, 0);            // Number of integer work vector elements
    ssSetNumPWork(S, 0);            // Number of pointer work vector elements
    ssSetNumModes(S, 0);            // Number of mode work vector elements
    ssSetNumNonsampledZCs(S, 0);    // Number of nonsampled zero crossings
    
    /* Specify the sim state compliance to be same as a built-in block */
    /* see sfun_simstate.c for example of other possible settings */
    ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);
    
    /*
     * All options have the form SS_OPTION_<name> and are documented in
     * matlabroot/simulink/include/simstruc.h. The options should be
     * bitwise or'd together as in
     *   ssSetOptions(S, (SS_OPTION_name1 | SS_OPTION_name2))
     */
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}



/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    S-function is comprised of only continuous sample time elements
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S);
    
}



#define MDL_INITIALIZE_CONDITIONS
/* Function: mdlInitializeConditions ========================================
 * Abstract:
 *    If the initial condition parameter (X0) is not an empty matrix,
 *    then use it to set up the initial conditions, otherwise,
 *    set the initial conditions to all 0.0
 */
static void mdlInitializeConditions(SimStruct *S)
{
    real_T *x0 = ssGetContStates(S);
    int_T  i, nStates;
    
    nStates = ssGetNumContStates(S);
    if (mxGetM(X0_PARAM(S)) != 0) {
        const real_T *pr = mxGetPr(X0_PARAM(S));
        
        for (i = 0; i < nStates; i++) {
            *x0++ = *pr++;
        }
    } else {
        for (i = 0; i < nStates; i++) {
            *x0++ = 0.0;
        }
    }
}



/* Function: mdlOutputs =======================================================
 * Abstract:
 *      y = Cx + Du
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
    real_T            *y       = ssGetOutputPortRealSignal(S,0);
    real_T            *y2      = ssGetOutputPortRealSignal(S,1);
    real_T            *x       = ssGetContStates(S);
    
    UNUSED_ARG(tid); /* not used in single tasking mode */

    /* Assign output 1 - states */
    y[0] = x[0];
    y[1] = x[1];
    y[2] = x[2];
    y[3] = x[3];
    y[4] = x[4];
    y[5] = x[5];
    y[6] = x[6];
    y[7] = x[7];
    
    /* Gyroscope VRML quaternion angle transformations ====================
     */
    
    double q1 = (-1) * x[0];
    double q2 = (-1) * x[1];
    double q3 = (-1) * x[2];
    double q4 = x[3];
    
    /* Default rotation axis q1, q2 */
    double q10 = C_PI/2;
    double q20 = C_PI/2;
    
    /* Declare rotation axis for VRML2Quaternions() */
    double rotationq1[3]    = {1,0,0};
    double rotationq10[3]   = {0,1,0};
    double rotationq2[3]    = {0,0,1};
    double rotationq20[3]   = {1,0,0};
    double rotationq3[3]    = {1,0,0};
    double rotationq4[3]    = {0,1,0};
    
    /* Declare outputs for VRML2Quaternion(), multiplicateQuaternions() */
    double Q1[4]    = {0,0,0,0};
    double Q2[4]    = {0,0,0,0};
    double Q3[4]    = {0,0,0,0};
    double Q4[4]    = {0,0,0,0};
    double Q34[4]   = {0,0,0,0};
    double Q20[4]   = {0,0,0,0};
    double Q234[4]  = {0,0,0,0};
    double Q2034[4] = {0,0,0,0};
    double Q10[4]   = {0,0,0,0};
    double Q1234[4] = {0,0,0,0};
    double Q10234[4]= {0,0,0,0};
    
    /* Declare outputs for Quaternion2VRML() */
    double VRMLq1[4] = {0,0,0,0};
    double VRMLq2[4] = {0,0,0,0};
    double VRMLq3[4] = {0,0,0,0};
    double VRMLq4[4] = {0,0,0,0};
    
    /* Calculate rotation of q4 */
    VRML2Quaternion(rotationq4, q4, Q4);
    Quaternion2VRML(Q4, VRMLq4);
    
    /* Calculate rotation of q3 */
    VRML2Quaternion(rotationq3, q3, Q3);
    multiplicateQuaternions(Q4, Q3, Q34);
    Quaternion2VRML(Q34, VRMLq3);
    
    /* Calculate rotation of q2 */
    VRML2Quaternion(rotationq2, q2, Q2);
    VRML2Quaternion(rotationq20, q20, Q20);
    multiplicateQuaternions(Q34, Q2, Q234);
    multiplicateQuaternions(Q234, Q20, Q2034);
    Quaternion2VRML(Q2034, VRMLq2);
    
    /* Calculate rotation of q1 */
    VRML2Quaternion(rotationq1, q1, Q1);
    VRML2Quaternion(rotationq10, q10, Q10);
    multiplicateQuaternions(Q234, Q1, Q1234);
    multiplicateQuaternions(Q1234, Q10, Q10234);
    Quaternion2VRML(Q10234, VRMLq1);
    
    /* End VRML quaternion angle transformations ==========================
     */
    
    /* Declare 2nd outputs */
    // q1
    y2[0] = VRMLq1[0];
    y2[1] = VRMLq1[1];
    y2[2] = VRMLq1[2];
    y2[3] = VRMLq1[3];
    //q2
    y2[4] = VRMLq2[0];
    y2[5] = VRMLq2[1];
    y2[6] = VRMLq2[2];
    y2[7] = VRMLq2[3];
    //q3
    y2[8] = VRMLq3[0];
    y2[9] = VRMLq3[1];
    y2[10] = VRMLq3[2];
    y2[11] = VRMLq3[3];
    //q4
    y2[12] = VRMLq4[0];
    y2[13] = VRMLq4[1];
    y2[14] = VRMLq4[2];
    y2[15] = VRMLq4[3];
    
}


#define MDL_DERIVATIVES
/* Function: mdlDerivatives =================================================
 * Abstract:
 *      xdot = Ax + Bu
 */
static void mdlDerivatives(SimStruct *S)
{
    real_T            *dx      = ssGetdX(S);
    real_T            *x       = ssGetContStates(S);
    InputRealPtrsType u0       = ssGetInputPortRealSignalPtrs(S,0);
    InputRealPtrsType u1       = ssGetInputPortRealSignalPtrs(S,1);
    InputRealPtrsType u2       = ssGetInputPortRealSignalPtrs(S,2);
    InputRealPtrsType u3       = ssGetInputPortRealSignalPtrs(S,3);
    real_T            *w1lock  = (real_T*)mxGetData(w1l_PARAM(S));
    real_T            *q2lock  = (real_T*)mxGetData(q2l_PARAM(S));
    real_T            *q3lock  = (real_T*)mxGetData(q3l_PARAM(S));
    real_T            *q4lock  = (real_T*)mxGetData(q4l_PARAM(S));
    real_T            *f       = (real_T*)mxGetData(f_PARAM(S));
    real_T            *Inertia = (real_T*)mxGetData(inertia_PARAM(S));
    real_T            *Km      = (real_T*)mxGetData(Km_PARAM(S));
    double            A[4]     = {0,0,0,0};
    double            Edss[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double            u[4]     = {*u0[0],*u1[0],*u2[0],*u3[0]};
    
    /* Define gyroscope parameters */
    // par = {iB,iC,iD,jB,jC,jD,kA,kB,kC,fv1,fv2,fv3,fv4,Km1,Km2,Km3,Km4} (initialize everything as 0)
    double par[17] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    par[0]  = Inertia[0];   // iB
    par[1]  = Inertia[1];   // iC
    par[2]  = Inertia[2];   // iD
    par[3]  = Inertia[3];   // jB
    par[4]  = Inertia[4];   // jC
    par[5]  = Inertia[5];   // jD
    par[6]  = Inertia[6];   // kA
    par[7]  = Inertia[7];   // kB
    par[8]  = Inertia[8];   // kC
    par[9]  = f[0];         // Disk viscous friction
    par[10] = f[1];         // BG viscous friction
    par[11] = f[2];         // RG viscous friction
    par[12] = f[3];         // SG viscous friction
    par[13] = Km[0];        // Motor constant 1
    par[14] = Km[1];        // Motor constant 2
    par[15] = Km[2];        // Motor constant 3
    par[16] = Km[3];        // Motor constant 4
    
    // Bitwise shift the numbers for bitwise addition to create the lookupCombination number
    int disk = (int)*w1lock << 3;
    int bg = (int)*q2lock << 2;
    int rg = (int)*q3lock << 1;
    int sg = (int)*q4lock;
    /* Create lookupCombination for gyroscope lock combinations based on the following lookupTable:
     
     Each body has the states 0 (free) or 1 (locked), 2^4 combinations
     
     #    | 0     1     2     3     4     5     6     7     8     9     10    11    12    13    14    15
     -----|---------------------------------------------------------------------------------------------
     disk | 0     0     0     0     0     0     0     0     1     1     1     1     1     1     1     1
     bg   | 0     0     0     0     1     1     1     1     0     0     0     0     1     1     1     1
     rg   | 0     0     1     1     0     0     1     1     0     0     1     1     0     0     1     1
     sg   | 0     1     0     1     0     1     0     1     0     1     0     1     0     1     0     1
     
     */
    // Lock combination is bitwise "or" addition of the lock states
    int lockCombination = disk | bg | rg | sg;
    
    /* Computate gyro dynamics for locking configuration */
    gyro_eq(lockCombination, par, x, u, Edss, A, dx);
    
}
/* Function Definitions */
void VRML2Quaternion(double axis[4], double angle, double quaternion[4]) {
    /*
     * Arguments    : double axis[3]
     *                double angle
     *                double quaternion[4]
     * Return Type  : void
     *                double quaternion[4]
     */
    
    double s;
    s = sin(angle/2);
    for (int i = 0; i<3; i++) {
        quaternion[i] = axis[i]*s;
    }
    quaternion[3] = cos(angle/2);
}

void multiplicateQuaternions(double q1[4], double q2[4], double result[4]) {
    /*
     * Arguments    : double q1[4]
     *                double q2[4]
     *                double result[4]
     * Return Type  : void
     *                double result[4]
     */
    
    result[0] = q1[3]*q2[0] + q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1];
    result[1] = q1[3]*q2[1] + q1[1]*q2[3] + q1[2]*q2[0] - q1[0]*q2[2];
    result[2] = q1[3]*q2[2] + q1[2]*q2[3] + q1[0]*q2[1] - q1[1]*q2[0];
    result[3] = q1[3]*q2[3] - q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2];
}

void Quaternion2VRML(double quaternion[4], double vrml[4]) {
    /*
     * Arguments    : double quaternion[4]
     *                double vrml[4]
     * Return Type  : void
     *                double vrml[4]
     */
    
    double s;
    s = sqrt(quaternion[0]*quaternion[0]+quaternion[1]*quaternion[1]+quaternion[2]*quaternion[2]);
    if (s<0.01) {
        vrml[0] = 1;
        vrml[1] = 0;
        vrml[2] = 0;
        vrml[3] = 0;
    } else {
        for (int i = 0; i<3; i++) {
            vrml[i] = quaternion[i]/s;
        }
        vrml[3] = acos(quaternion[3])*2;
    }
}

void gyro_eq(int lockCombination, const double parameters[17], double x[8], double u[4], double Edss[16], double A[4], double dxdt[8]) {
    /*
     * Arguments    : int LockCombination
     *                const double parameters[13]
     *                double x[8]
     *                double u[4]
     *                double Edss[16]
     *                double A[4]
     *                double dxdt[8]
     * Return Type  : void
     *                double dxdt[8]
     */

    // Multiply inputs (A) with motor constants (Nm/A) to obtain motor torques (Nm);
    u[0] *= parameters[13];
    u[1] *= parameters[14];
    u[2] *= parameters[15];
    u[3] *= parameters[16];

    switch (lockCombination) {
        case 0:
        {
            // q1=free; q2=free; q3=free; q4=free;
            
            /* Gyroscope -E matrix */
            Edss[0] = (-1)*(-parameters[5]);
            Edss[1] = 0.0;
            Edss[2] = (-1)*(-parameters[5] * cos(x[1]));
            Edss[3] = (-1)*(-parameters[5] * cos(x[2]) * sin(x[1]));
            
            Edss[4] = 0.0;
            Edss[5] = (-1)*(-parameters[1] + -parameters[2]);
            Edss[6] = 0.0;
            Edss[7] = (-1)*((parameters[1] + parameters[2]) * sin(x[2]));
            
            Edss[8] = (-1)*(-parameters[5] * cos(x[1]));
            Edss[9] = 0.0;
            Edss[10] = (-1)*(((-parameters[3] + -parameters[4]) + -parameters[5]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * (sin(x[1]) * sin(x[1])));
            Edss[11] = (-1)*((((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(x[1]) * cos(x[2]) * sin(x[1]));
            
            Edss[12] = (-1)*(-parameters[5] * cos(x[2]) * sin(x[1]));
            Edss[13] = (-1)*((parameters[1] + parameters[2]) * sin(x[2]));
            Edss[14] = (-1)*((((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(x[1]) * cos(x[2]) * sin(x[1]));
            Edss[15] = (-1)*(((((-parameters[2] + -parameters[6]) + -parameters[7]) + -parameters[8]) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * (cos(x[2]) * cos(x[2])) * (sin(x[1]) * sin(x[1]))) + (((-parameters[0] + -parameters[1]) + parameters[7]) + parameters[8]) * (sin(x[2]) * sin(x[2])));
            
            /* Gyroscope A matrix */
            A[0] = ((u[0] + -x[4] * parameters[9]) + -parameters[5] * x[5] * x[7] * cos(x[1]) * cos(x[2])) + parameters[5] * x[6] * sin(x[1]) * (x[5] + x[7] * sin(x[2]));
            A[1] = ((u[1] + -x[5] * parameters[10]) + x[7] * (parameters[5] * x[4] * cos(x[1]) + x[6] * ((parameters[1] + parameters[2]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * cos(2.0 * x[1]))) * cos(x[2])) + -(parameters[5] * x[4] * x[6] + -(((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(x[1]) * (x[6] * x[6] + -(x[7] * x[7]) * (cos(x[2]) * cos(x[2])))) * sin(x[1]);
            A[2] = (((u[2] + -x[6] * parameters[11]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * x[5] * x[6] * sin(2.0 * x[1])) + parameters[5] * x[4] * sin(x[1]) * (x[5] + -x[7] * sin(x[2]))) + 0.5 * x[7] * cos(x[2]) * ((-2.0 * (parameters[1] + parameters[2]) * x[5] + ((((((2.0 * parameters[0] + 2.0 * parameters[1]) + parameters[2]) + -parameters[4]) + -parameters[5]) + -2.0 * parameters[7]) + -parameters[8]) * x[7] * sin(x[2])) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(2.0 * x[1]) * (2.0 * x[5] + -x[7] * sin(x[2])));
            A[3] = (((u[3] + -x[7] * parameters[12]) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * x[5] * x[7] * (cos(x[2]) * cos(x[2])) * sin(2.0 * x[1])) + x[6] * (parameters[5] * x[4] + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * x[6] * cos(x[1])) * sin(x[1]) * sin(x[2])) + cos(x[2]) * ((((parameters[1] + parameters [2]) * x[5] * x[6] + -parameters[5] * x[4] * x[5] * cos(x[1])) + ((((((-2.0 * parameters[0] + -2.0 * parameters[1]) + -parameters[2]) + parameters[4]) + parameters[5]) + 2.0 * parameters[7]) + parameters[8]) * x[6] * x[7] * sin(x[2])) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * x[6] * cos(2.0 * x[1]) * (x[5] + x[7] * sin(x[2])));
            
            
            /* Solve system of linear equations Ex = A */
            Matrix4d E0breve;
            Vector4d A0breve;
            Vector4d l0Sol;
            E0breve <<  Edss[0],Edss[1],Edss[2],Edss[3],
            Edss[4],Edss[5],Edss[6],Edss[7],
            Edss[8],Edss[9],Edss[10],Edss[11],
            Edss[12],Edss[13],Edss[14],Edss[15];
            A0breve <<  A[0],A[1],A[2],A[3];
            l0Sol = E0breve.partialPivLu().solve(A0breve);
            
            /* Declare Outputs */
            dxdt[0] = x[4];                 //vel1
            dxdt[1] = x[5];                 //vel2
            dxdt[2] = x[6];                 //vel3
            dxdt[3] = x[7];                 //vel4
            dxdt[4] = *(l0Sol.data());      //acc1
            dxdt[5] = *(l0Sol.data()+1);    //acc2
            dxdt[6] = *(l0Sol.data()+2);    //acc3
            dxdt[7] = *(l0Sol.data()+3);    //acc4
            break;
        }
        case 1:
        {
            // q1=free; q2=free; q3=free; q4=locked;
            
            /* Gyroscope -E matrix */
            Edss[0] = (-1)*(-parameters[5]);
            Edss[1] = 0.0;
            Edss[2] = (-1)*(-parameters[5] * cos(x[1]));
            Edss[3] = 0;
            
            Edss[4] = 0.0;
            Edss[5] = (-1)*(-parameters[1] + -parameters[2]);
            Edss[6] = 0.0;
            Edss[7] = 0;
            
            Edss[8] = (-1)*(-parameters[5] * cos(x[1]));
            Edss[9] = 0.0;
            Edss[10] = (-1)*(((-parameters[3] + -parameters[4]) + -parameters[5]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * (sin(x[1]) * sin(x[1])));
            Edss[11] = 0;
            
            Edss[12] = 0;
            Edss[13] = 0;
            Edss[14] = 0;
            Edss[15] = 0;
            
            /* Gyroscope A matrix */
            A[0] = ((u[0] + -x[4] * parameters[9]) + -parameters[5] * x[5] * x[7] * cos(x[1]) * cos(x[2])) + parameters[5] * x[6] * sin(x[1]) * (x[5] + x[7] * sin(x[2]));
            A[1] = ((u[1] + -x[5] * parameters[10]) + x[7] * (parameters[5] * x[4] * cos(x[1]) + x[6] * ((parameters[1] + parameters[2]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * cos(2.0 * x[1]))) * cos(x[2])) + -(parameters[5] * x[4] * x[6] + -(((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(x[1]) * (x[6] * x[6] + -(x[7] * x[7]) * (cos(x[2]) * cos(x[2])))) * sin(x[1]);
            A[2] = (((u[2] + -x[6] * parameters[11]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * x[5] * x[6] * sin(2.0 * x[1])) + parameters[5] * x[4] * sin(x[1]) * (x[5] + -x[7] * sin(x[2]))) + 0.5 * x[7] * cos(x[2]) * ((-2.0 * (parameters[1] + parameters[2]) * x[5] + ((((((2.0 * parameters[0] + 2.0 * parameters[1]) + parameters[2]) + -parameters[4]) + -parameters[5]) + -2.0 * parameters[7]) + -parameters[8]) * x[7] * sin(x[2])) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(2.0 * x[1]) * (2.0 * x[5] + -x[7] * sin(x[2])));
            A[3] = 0;
            
            /* Solve system of linear equations Ex = A */
            Matrix3d E1breve;
            Vector3d A1breve;
            Vector3d l1Sol;
            E1breve <<  Edss[0],Edss[1],Edss[2],
            Edss[4],Edss[5],Edss[6],
            Edss[8],Edss[9],Edss[10];
            A1breve <<  A[0],A[1],A[2];
            l1Sol = E1breve.partialPivLu().solve(A1breve);
            
            /* Declare Outputs */
            dxdt[0] = x[4];                 //vel1
            dxdt[1] = x[5];                 //vel2
            dxdt[2] = x[6];                 //vel3
            dxdt[3] = x[7];                 //vel4
            dxdt[4] = *(l1Sol.data());      //acc1
            dxdt[5] = *(l1Sol.data()+1);    //acc2
            dxdt[6] = *(l1Sol.data()+2);    //acc3
            dxdt[7] = 0;                    //acc4
            break;
        }
        case 2:
        {
            // q1=free; q2=free; q3=locked; q4=free;
            
            /* Gyroscope -E matrix */
            Edss[0] = (-1)*(-parameters[5]);
            Edss[1] = 0.0;
            Edss[2] = 0;
            Edss[3] = (-1)*(-parameters[5] * cos(x[2]) * sin(x[1]));
            
            Edss[4] = 0.0;
            Edss[5] = (-1)*(-parameters[1] + -parameters[2]);
            Edss[6] = 0.0;
            Edss[7] = (-1)*((parameters[1] + parameters[2]) * sin(x[2]));
            
            Edss[8] = 0;
            Edss[9] = 0.0;
            Edss[10] = 0;
            Edss[11] = 0;
            
            Edss[12] = (-1)*(-parameters[5] * cos(x[2]) * sin(x[1]));
            Edss[13] = (-1)*((parameters[1] + parameters[2]) * sin(x[2]));
            Edss[14] = 0;
            Edss[15] = (-1)*(((((-parameters[2] + -parameters[6]) + -parameters[7]) + -parameters[8]) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * (cos(x[2]) * cos(x[2])) * (sin(x[1]) * sin(x[1]))) + (((-parameters[0] + -parameters[1]) + parameters[7]) + parameters[8]) * (sin(x[2]) * sin(x[2])));
            
            /* Gyroscope A matrix */
            A[0] = ((u[0] + -x[4] * parameters[9]) + -parameters[5] * x[5] * x[7] * cos(x[1]) * cos(x[2])) + parameters[5] * x[6] * sin(x[1]) * (x[5] + x[7] * sin(x[2]));
            A[1] = ((u[1] + -x[5] * parameters[10]) + x[7] * (parameters[5] * x[4] * cos(x[1]) + x[6] * ((parameters[1] + parameters[2]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * cos(2.0 * x[1]))) * cos(x[2])) + -(parameters[5] * x[4] * x[6] + -(((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(x[1]) * (x[6] * x[6] + -(x[7] * x[7]) * (cos(x[2]) * cos(x[2])))) * sin(x[1]);
            A[2] = 0;
            A[3] = (((u[3] + -x[7] * parameters[12]) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * x[5] * x[7] * (cos(x[2]) * cos(x[2])) * sin(2.0 * x[1])) + x[6] * (parameters[5] * x[4] + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * x[6] * cos(x[1])) * sin(x[1]) * sin(x[2])) + cos(x[2]) * ((((parameters[1] + parameters [2]) * x[5] * x[6] + -parameters[5] * x[4] * x[5] * cos(x[1])) + ((((((-2.0 * parameters[0] + -2.0 * parameters[1]) + -parameters[2]) + parameters[4]) + parameters[5]) + 2.0 * parameters[7]) + parameters[8]) * x[6] * x[7] * sin(x[2])) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * x[6] * cos(2.0 * x[1]) * (x[5] + x[7] * sin(x[2])));
            
            
            /* Solve system of linear equations Ex = A */
            Matrix3d E2breve;
            Vector3d A2breve;
            Vector3d l2Sol;
            E2breve <<  Edss[0],Edss[1],Edss[3],
                        Edss[4],Edss[5],Edss[7],
                        Edss[12],Edss[13],Edss[15];
            A2breve <<  A[0],A[1],A[3];
            l2Sol = E2breve.partialPivLu().solve(A2breve);
            
            /* Declare Outputs */
            dxdt[0] = x[4];                 //vel1
            dxdt[1] = x[5];                 //vel2
            dxdt[2] = x[6];                 //vel3
            dxdt[3] = x[7];                 //vel4
            dxdt[4] = *(l2Sol.data());      //acc1
            dxdt[5] = *(l2Sol.data()+1);    //acc2
            dxdt[6] = 0;                    //acc3
            dxdt[7] = *(l2Sol.data()+2);    //acc4
            break;
        }
        case 3:
        {
            // q1=free; q2=free; q3=locked; q4=locked;
            
            /* Gyroscope -E matrix */
            Edss[0] = (-1)*(-parameters[5]);
            Edss[1] = 0.0;
            Edss[2] = 0;
            Edss[3] = 0;
            
            Edss[4] = 0.0;
            Edss[5] = (-1)*(-parameters[1] + -parameters[2]);
            Edss[6] = 0.0;
            Edss[7] = 0;
            
            Edss[8] = 0;
            Edss[9] = 0.0;
            Edss[10] = 0;
            Edss[11] = 0;
            
            Edss[12] = 0;
            Edss[13] = 0;
            Edss[14] = 0;
            Edss[15] = 0;
            
            /* Gyroscope A matrix */
            A[0] = ((u[0] + -x[4] * parameters[9]) + -parameters[5] * x[5] * x[7] * cos(x[1]) * cos(x[2])) + parameters[5] * x[6] * sin(x[1]) * (x[5] + x[7] * sin(x[2]));
            A[1] = ((u[1] + -x[5] * parameters[10]) + x[7] * (parameters[5] * x[4] * cos(x[1]) + x[6] * ((parameters[1] + parameters[2]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * cos(2.0 * x[1]))) * cos(x[2])) + -(parameters[5] * x[4] * x[6] + -(((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(x[1]) * (x[6] * x[6] + -(x[7] * x[7]) * (cos(x[2]) * cos(x[2])))) * sin(x[1]);
            A[2] = 0;
            A[3] = 0;
            
            
            /* Solve system of linear equations Ex = A */
            Matrix2d E3breve;
            Vector2d A3breve;
            Vector2d l3Sol;
            E3breve <<  Edss[0],Edss[1],
            Edss[4],Edss[5];
            A3breve <<  A[0],A[1];
            l3Sol = E3breve.partialPivLu().solve(A3breve);
            
            /* Declare Outputs */
            dxdt[0] = x[4];                 //vel1
            dxdt[1] = x[5];                 //vel2
            dxdt[2] = x[6];                 //vel3
            dxdt[3] = x[7];                 //vel4
            dxdt[4] = *(l3Sol.data());      //acc1
            dxdt[5] = *(l3Sol.data()+1);    //acc2
            dxdt[6] = 0;                    //acc3
            dxdt[7] = 0;                    //acc4
            break;
        }
        case 4:
        {
            // q1=free; q2=locked; q3=free; q4=free;
            
            /* Gyroscope -E matrix */
            Edss[0] = (-1)*(-parameters[5]);
            Edss[1] = 0.0;
            Edss[2] = (-1)*(-parameters[5] * cos(x[1]));
            Edss[3] = (-1)*(-parameters[5] * cos(x[2]) * sin(x[1]));
            
            Edss[4] = 0.0;
            Edss[5] = 0;
            Edss[6] = 0.0;
            Edss[7] = 0;
            
            Edss[8] = (-1)*(-parameters[5] * cos(x[1]));
            Edss[9] = 0.0;
            Edss[10] = (-1)*(((-parameters[3] + -parameters[4]) + -parameters[5]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * (sin(x[1]) * sin(x[1])));
            Edss[11] = (-1)*((((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(x[1]) * cos(x[2]) * sin(x[1]));
            
            Edss[12] = (-1)*(-parameters[5] * cos(x[2]) * sin(x[1]));
            Edss[13] = 0;
            Edss[14] = (-1)*((((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(x[1]) * cos(x[2]) * sin(x[1]));
            Edss[15] = (-1)*(((((-parameters[2] + -parameters[6]) + -parameters[7]) + -parameters[8]) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * (cos(x[2]) * cos(x[2])) * (sin(x[1]) * sin(x[1]))) + (((-parameters[0] + -parameters[1]) + parameters[7]) + parameters[8]) * (sin(x[2]) * sin(x[2])));
            
            /* Gyroscope A matrix */
            A[0] = ((u[0] + -x[4] * parameters[9]) + -parameters[5] * x[5] * x[7] * cos(x[1]) * cos(x[2])) + parameters[5] * x[6] * sin(x[1]) * (x[5] + x[7] * sin(x[2]));
            A[1] = 0;
            A[2] = (((u[2] + -x[6] * parameters[11]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * x[5] * x[6] * sin(2.0 * x[1])) + parameters[5] * x[4] * sin(x[1]) * (x[5] + -x[7] * sin(x[2]))) + 0.5 * x[7] * cos(x[2]) * ((-2.0 * (parameters[1] + parameters[2]) * x[5] + ((((((2.0 * parameters[0] + 2.0 * parameters[1]) + parameters[2]) + -parameters[4]) + -parameters[5]) + -2.0 * parameters[7]) + -parameters[8]) * x[7] * sin(x[2])) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(2.0 * x[1]) * (2.0 * x[5] + -x[7] * sin(x[2])));
            A[3] = (((u[3] + -x[7] * parameters[12]) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * x[5] * x[7] * (cos(x[2]) * cos(x[2])) * sin(2.0 * x[1])) + x[6] * (parameters[5] * x[4] + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * x[6] * cos(x[1])) * sin(x[1]) * sin(x[2])) + cos(x[2]) * ((((parameters[1] + parameters [2]) * x[5] * x[6] + -parameters[5] * x[4] * x[5] * cos(x[1])) + ((((((-2.0 * parameters[0] + -2.0 * parameters[1]) + -parameters[2]) + parameters[4]) + parameters[5]) + 2.0 * parameters[7]) + parameters[8]) * x[6] * x[7] * sin(x[2])) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * x[6] * cos(2.0 * x[1]) * (x[5] + x[7] * sin(x[2])));
            
            
            /* Solve system of linear equations Ex = A */
            Matrix3d E4breve;
            Vector3d A4breve;
            Vector3d l4Sol;
            E4breve <<  Edss[0],Edss[2],Edss[3],
            Edss[8],Edss[10],Edss[11],
            Edss[12],Edss[14],Edss[15];
            A4breve <<  A[0],A[2],A[3];
            l4Sol = E4breve.partialPivLu().solve(A4breve);
            
            /* Declare Outputs */
            dxdt[0] = x[4];                 //vel1
            dxdt[1] = x[5];                 //vel2
            dxdt[2] = x[6];                 //vel3
            dxdt[3] = x[7];                 //vel4
            dxdt[4] = *(l4Sol.data());      //acc1
            dxdt[5] = 0;                    //acc2
            dxdt[6] = *(l4Sol.data()+1);    //acc3
            dxdt[7] = *(l4Sol.data()+2);    //acc4
            break;
        }
        case 5:
        {
            // q1=free; q2=locked; q3=free; q4=locked;
            
            /* Gyroscope -E matrix */
            Edss[0] = (-1)*(-parameters[5]);
            Edss[1] = 0.0;
            Edss[2] = (-1)*(-parameters[5] * cos(x[1]));
            Edss[3] = 0;
            
            Edss[4] = 0.0;
            Edss[5] = 0;
            Edss[6] = 0.0;
            Edss[7] = 0;
            
            Edss[8] = (-1)*(-parameters[5] * cos(x[1]));
            Edss[9] = 0.0;
            Edss[10] = (-1)*(((-parameters[3] + -parameters[4]) + -parameters[5]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * (sin(x[1]) * sin(x[1])));
            Edss[11] = 0;
            
            Edss[12] = 0;
            Edss[13] = 0;
            Edss[14] = 0;
            Edss[15] = 0;
            
            /* Gyroscope A matrix */
            A[0] = ((u[0] + -x[4] * parameters[9]) + -parameters[5] * x[5] * x[7] * cos(x[1]) * cos(x[2])) + parameters[5] * x[6] * sin(x[1]) * (x[5] + x[7] * sin(x[2]));
            A[1] = 0;
            A[2] = (((u[2] + -x[6] * parameters[11]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * x[5] * x[6] * sin(2.0 * x[1])) + parameters[5] * x[4] * sin(x[1]) * (x[5] + -x[7] * sin(x[2]))) + 0.5 * x[7] * cos(x[2]) * ((-2.0 * (parameters[1] + parameters[2]) * x[5] + ((((((2.0 * parameters[0] + 2.0 * parameters[1]) + parameters[2]) + -parameters[4]) + -parameters[5]) + -2.0 * parameters[7]) + -parameters[8]) * x[7] * sin(x[2])) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(2.0 * x[1]) * (2.0 * x[5] + -x[7] * sin(x[2])));
            A[3] = 0;
            
            
            /* Solve system of linear equations Ex = A */
            Matrix2d E5breve;
            Vector2d A5breve;
            Vector2d l5Sol;
            E5breve <<  Edss[0],Edss[2],
            Edss[8],Edss[10];
            A5breve <<  A[0],A[2];
            l5Sol = E5breve.partialPivLu().solve(A5breve);
            
            /* Declare Outputs */
            dxdt[0] = x[4];                 //vel1
            dxdt[1] = x[5];                 //vel2
            dxdt[2] = x[6];                 //vel3
            dxdt[3] = x[7];                 //vel4
            dxdt[4] = *(l5Sol.data());      //acc1
            dxdt[5] = 0;                    //acc2
            dxdt[6] = *(l5Sol.data()+1);    //acc3
            dxdt[7] = 0;                    //acc4
            break;
        }
        case 6:
        {
            // q1=free; q2=locked; q3=locked; q4=free;
            
            /* Gyroscope -E matrix */
            Edss[0] = (-1)*(-parameters[5]);
            Edss[1] = 0.0;
            Edss[2] = 0;
            Edss[3] = (-1)*(-parameters[5] * cos(x[2]) * sin(x[1]));
            
            Edss[4] = 0.0;
            Edss[5] = 0;
            Edss[6] = 0.0;
            Edss[7] = 0;
            
            Edss[8] = 0;
            Edss[9] = 0.0;
            Edss[10] = 0;
            Edss[11] = 0;
            
            Edss[12] = (-1)*(-parameters[5] * cos(x[2]) * sin(x[1]));
            Edss[13] = 0;
            Edss[14] = 0;
            Edss[15] = (-1)*(((((-parameters[2] + -parameters[6]) + -parameters[7]) + -parameters[8]) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * (cos(x[2]) * cos(x[2])) * (sin(x[1]) * sin(x[1]))) + (((-parameters[0] + -parameters[1]) + parameters[7]) + parameters[8]) * (sin(x[2]) * sin(x[2])));
            
            /* Gyroscope A matrix */
            A[0] = ((u[0] + -x[4] * parameters[9]) + -parameters[5] * x[5] * x[7] * cos(x[1]) * cos(x[2])) + parameters[5] * x[6] * sin(x[1]) * (x[5] + x[7] * sin(x[2]));
            A[1] = 0;
            A[2] = 0;
            A[3] = (((u[3] + -x[7] * parameters[12]) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * x[5] * x[7] * (cos(x[2]) * cos(x[2])) * sin(2.0 * x[1])) + x[6] * (parameters[5] * x[4] + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * x[6] * cos(x[1])) * sin(x[1]) * sin(x[2])) + cos(x[2]) * ((((parameters[1] + parameters [2]) * x[5] * x[6] + -parameters[5] * x[4] * x[5] * cos(x[1])) + ((((((-2.0 * parameters[0] + -2.0 * parameters[1]) + -parameters[2]) + parameters[4]) + parameters[5]) + 2.0 * parameters[7]) + parameters[8]) * x[6] * x[7] * sin(x[2])) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * x[6] * cos(2.0 * x[1]) * (x[5] + x[7] * sin(x[2])));
            
            
            /* Solve system of linear equations Ex = A */
            Matrix2d E6breve;
            Vector2d A6breve;
            Vector2d l6Sol;
            E6breve <<  Edss[0],Edss[3],
            Edss[12],Edss[15];
            A6breve <<  A[0],A[3];
            l6Sol = E6breve.partialPivLu().solve(A6breve);
            
            /* Declare Outputs */
            dxdt[0] = x[4];                 //vel1
            dxdt[1] = x[5];                 //vel2
            dxdt[2] = x[6];                 //vel3
            dxdt[3] = x[7];                 //vel4
            dxdt[4] = *(l6Sol.data());      //acc1
            dxdt[5] = 0;                    //acc2
            dxdt[6] = 0;                    //acc3
            dxdt[7] = *(l6Sol.data()+1);    //acc4
            break;
        }
        case 7:
        {
            // q1=free; q2=locked; q3=locked; q4=locked;
            
            /* Gyroscope -E matrix */
            Edss[0] = (-1)*(-parameters[5]);
            Edss[1] = 0.0;
            Edss[2] = 0;
            Edss[3] = 0;
            
            Edss[4] = 0.0;
            Edss[5] = 0;
            Edss[6] = 0.0;
            Edss[7] = 0;
            
            Edss[8] = 0;
            Edss[9] = 0.0;
            Edss[10] = 0;
            Edss[11] = 0;
            
            Edss[12] = 0;
            Edss[13] = 0;
            Edss[14] = 0;
            Edss[15] = 0;
            
            /* Gyroscope A matrix */
            A[0] = ((u[0] + -x[4] * parameters[9]) + -parameters[5] * x[5] * x[7] * cos(x[1]) * cos(x[2])) + parameters[5] * x[6] * sin(x[1]) * (x[5] + x[7] * sin(x[2]));
            A[1] = 0;
            A[2] = 0;
            A[3] = 0;
            
            
            /* Solve system of linear equations Ex = A */
            double l7Sol;
            l7Sol = (1/Edss[0])*A[0];
            
            /* Declare Outputs */
            dxdt[0] = x[4];                 //vel1
            dxdt[1] = x[5];                 //vel2
            dxdt[2] = x[6];                 //vel3
            dxdt[3] = x[7];                 //vel4
            dxdt[4] = l7Sol;                //acc1
            dxdt[5] = 0;                    //acc2
            dxdt[6] = 0;                    //acc3
            dxdt[7] = 0;                    //acc4
            break;
        }
        case 8:
        {
            // q1=locked; q2=free; q3=free; q4=free;
            
            /* Gyroscope -E matrix */
            Edss[0] = 0;
            Edss[1] = 0.0;
            Edss[2] = 0;
            Edss[3] = 0;
            
            Edss[4] = 0.0;
            Edss[5] = (-1)*(-parameters[1] + -parameters[2]);
            Edss[6] = 0.0;
            Edss[7] = (-1)*((parameters[1] + parameters[2]) * sin(x[2]));
            
            Edss[8] = 0;
            Edss[9] = 0.0;
            Edss[10] = (-1)*(((-parameters[3] + -parameters[4]) + -parameters[5]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * (sin(x[1]) * sin(x[1])));
            Edss[11] = (-1)*((((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(x[1]) * cos(x[2]) * sin(x[1]));
            
            Edss[12] = 0;
            Edss[13] = (-1)*((parameters[1] + parameters[2]) * sin(x[2]));
            Edss[14] = (-1)*((((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(x[1]) * cos(x[2]) * sin(x[1]));
            Edss[15] = (-1)*(((((-parameters[2] + -parameters[6]) + -parameters[7]) + -parameters[8]) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * (cos(x[2]) * cos(x[2])) * (sin(x[1]) * sin(x[1]))) + (((-parameters[0] + -parameters[1]) + parameters[7]) + parameters[8]) * (sin(x[2]) * sin(x[2])));
            
            /* Gyroscope A matrix */
            A[0] = 0;
            A[1] = ((u[1] + -x[5] * parameters[10]) + x[7] * (parameters[5] * x[4] * cos(x[1]) + x[6] * ((parameters[1] + parameters[2]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * cos(2.0 * x[1]))) * cos(x[2])) + -(parameters[5] * x[4] * x[6] + -(((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(x[1]) * (x[6] * x[6] + -(x[7] * x[7]) * (cos(x[2]) * cos(x[2])))) * sin(x[1]);
            A[2] = (((u[2] + -x[6] * parameters[11]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * x[5] * x[6] * sin(2.0 * x[1])) + parameters[5] * x[4] * sin(x[1]) * (x[5] + -x[7] * sin(x[2]))) + 0.5 * x[7] * cos(x[2]) * ((-2.0 * (parameters[1] + parameters[2]) * x[5] + ((((((2.0 * parameters[0] + 2.0 * parameters[1]) + parameters[2]) + -parameters[4]) + -parameters[5]) + -2.0 * parameters[7]) + -parameters[8]) * x[7] * sin(x[2])) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(2.0 * x[1]) * (2.0 * x[5] + -x[7] * sin(x[2])));
            A[3] = (((u[3] + -x[7] * parameters[12]) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * x[5] * x[7] * (cos(x[2]) * cos(x[2])) * sin(2.0 * x[1])) + x[6] * (parameters[5] * x[4] + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * x[6] * cos(x[1])) * sin(x[1]) * sin(x[2])) + cos(x[2]) * ((((parameters[1] + parameters [2]) * x[5] * x[6] + -parameters[5] * x[4] * x[5] * cos(x[1])) + ((((((-2.0 * parameters[0] + -2.0 * parameters[1]) + -parameters[2]) + parameters[4]) + parameters[5]) + 2.0 * parameters[7]) + parameters[8]) * x[6] * x[7] * sin(x[2])) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * x[6] * cos(2.0 * x[1]) * (x[5] + x[7] * sin(x[2])));
            
            
            /* Solve system of linear equations Ex = A */
            Matrix3d E8breve;
            Vector3d A8breve;
            Vector3d l8Sol;
            E8breve <<  Edss[5],Edss[6],Edss[7],
            Edss[9],Edss[10],Edss[11],
            Edss[13],Edss[14],Edss[15];
            A8breve <<  A[1],A[2],A[3];
            l8Sol = E8breve.partialPivLu().solve(A8breve);
            
            /* Declare Outputs */
            dxdt[0] = x[4];                 //vel1
            dxdt[1] = x[5];                 //vel2
            dxdt[2] = x[6];                 //vel3
            dxdt[3] = x[7];                 //vel4
            dxdt[4] = 0;                    //acc1
            dxdt[5] = *(l8Sol.data());      //acc2
            dxdt[6] = *(l8Sol.data()+1);    //acc3
            dxdt[7] = *(l8Sol.data()+2);    //acc4
            break;
        }
        case 9:
        {
            // q1=locked; q2=free; q3=free; q4=locked;
            
            /* Gyroscope -E matrix */
            Edss[0] = 0;
            Edss[1] = 0.0;
            Edss[2] = 0;
            Edss[3] = 0;
            
            Edss[4] = 0.0;
            Edss[5] = (-1)*(-parameters[1] + -parameters[2]);
            Edss[6] = 0.0;
            Edss[7] = 0;
            
            Edss[8] = 0;
            Edss[9] = 0.0;
            Edss[10] = (-1)*(((-parameters[3] + -parameters[4]) + -parameters[5]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * (sin(x[1]) * sin(x[1])));
            Edss[11] = 0;
            
            Edss[12] = 0;
            Edss[13] = 0;
            Edss[14] = 0;
            Edss[15] = 0;
            
            /* Gyroscope A matrix */
            A[0] = 0;
            A[1] = ((u[1] + -x[5] * parameters[10]) + x[7] * (parameters[5] * x[4] * cos(x[1]) + x[6] * ((parameters[1] + parameters[2]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * cos(2.0 * x[1]))) * cos(x[2])) + -(parameters[5] * x[4] * x[6] + -(((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(x[1]) * (x[6] * x[6] + -(x[7] * x[7]) * (cos(x[2]) * cos(x[2])))) * sin(x[1]);
            A[2] = (((u[2] + -x[6] * parameters[11]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * x[5] * x[6] * sin(2.0 * x[1])) + parameters[5] * x[4] * sin(x[1]) * (x[5] + -x[7] * sin(x[2]))) + 0.5 * x[7] * cos(x[2]) * ((-2.0 * (parameters[1] + parameters[2]) * x[5] + ((((((2.0 * parameters[0] + 2.0 * parameters[1]) + parameters[2]) + -parameters[4]) + -parameters[5]) + -2.0 * parameters[7]) + -parameters[8]) * x[7] * sin(x[2])) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(2.0 * x[1]) * (2.0 * x[5] + -x[7] * sin(x[2])));
            A[3] = 0;
            
            
            /* Solve system of linear equations Ex = A */
            Matrix2d E9breve;
            Vector2d A9breve;
            Vector2d l9Sol;
            E9breve <<  Edss[5],Edss[6],
            Edss[9],Edss[10];
            A9breve <<  A[1],A[2];
            l9Sol = E9breve.partialPivLu().solve(A9breve);
            
            /* Declare Outputs */
            dxdt[0] = x[4];                 //vel1
            dxdt[1] = x[5];                 //vel2
            dxdt[2] = x[6];                 //vel3
            dxdt[3] = x[7];                 //vel4
            dxdt[4] = 0;                    //acc1
            dxdt[5] = *(l9Sol.data());      //acc2
            dxdt[6] = *(l9Sol.data()+1);    //acc3
            dxdt[7] = 0;                    //acc4
            break;
        }
        case 10:
        {
            // q1=locked; q2=free; q3=locked; q4=free;
            
            /* Gyroscope -E matrix */
            Edss[0] = 0;
            Edss[1] = 0.0;
            Edss[2] = 0;
            Edss[3] = 0;
            
            Edss[4] = 0.0;
            Edss[5] = (-1)*(-parameters[1] + -parameters[2]);
            Edss[6] = 0.0;
            Edss[7] = (-1)*((parameters[1] + parameters[2]) * sin(x[2]));
            
            Edss[8] = 0;
            Edss[9] = 0.0;
            Edss[10] = 0;
            Edss[11] = 0;
            
            Edss[12] = 0;
            Edss[13] = (-1)*((parameters[1] + parameters[2]) * sin(x[2]));
            Edss[14] = 0;
            Edss[15] = (-1)*(((((-parameters[2] + -parameters[6]) + -parameters[7]) + -parameters[8]) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * (cos(x[2]) * cos(x[2])) * (sin(x[1]) * sin(x[1]))) + (((-parameters[0] + -parameters[1]) + parameters[7]) + parameters[8]) * (sin(x[2]) * sin(x[2])));
            
            /* Gyroscope A matrix */
            A[0] = 0;
            A[1] = ((u[1] + -x[5] * parameters[10]) + x[7] * (parameters[5] * x[4] * cos(x[1]) + x[6] * ((parameters[1] + parameters[2]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * cos(2.0 * x[1]))) * cos(x[2])) + -(parameters[5] * x[4] * x[6] + -(((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(x[1]) * (x[6] * x[6] + -(x[7] * x[7]) * (cos(x[2]) * cos(x[2])))) * sin(x[1]);
            A[2] = 0;
            A[3] = (((u[3] + -x[7] * parameters[12]) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * x[5] * x[7] * (cos(x[2]) * cos(x[2])) * sin(2.0 * x[1])) + x[6] * (parameters[5] * x[4] + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * x[6] * cos(x[1])) * sin(x[1]) * sin(x[2])) + cos(x[2]) * ((((parameters[1] + parameters [2]) * x[5] * x[6] + -parameters[5] * x[4] * x[5] * cos(x[1])) + ((((((-2.0 * parameters[0] + -2.0 * parameters[1]) + -parameters[2]) + parameters[4]) + parameters[5]) + 2.0 * parameters[7]) + parameters[8]) * x[6] * x[7] * sin(x[2])) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * x[6] * cos(2.0 * x[1]) * (x[5] + x[7] * sin(x[2])));
            
            
            /* Solve system of linear equations Ex = A */
            Matrix2d E10breve;
            Vector2d A10breve;
            Vector2d l10Sol;
            E10breve << Edss[5],Edss[7],
            Edss[13],Edss[15];
            A10breve << A[1],A[3];
            l10Sol = E10breve.partialPivLu().solve(A10breve);
            
            /* Declare Outputs */
            dxdt[0] = x[4];                 //vel1
            dxdt[1] = x[5];                 //vel2
            dxdt[2] = x[6];                 //vel3
            dxdt[3] = x[7];                 //vel4
            dxdt[4] = 0;                    //acc1
            dxdt[5] = *(l10Sol.data());   //acc2
            dxdt[6] = 0;                    //acc3
            dxdt[7] = *(l10Sol.data()+1);   //acc4
            break;
        }
        case 11:
        {
            // q1=locked; q2=free; q3=locked; q4=locked;
            
            /* Gyroscope -E matrix */
            Edss[0] = 0;
            Edss[1] = 0.0;
            Edss[2] = 0;
            Edss[3] = 0;
            
            Edss[4] = 0.0;
            Edss[5] = (-1)*(-parameters[1] + -parameters[2]);
            Edss[6] = 0.0;
            Edss[7] = 0;
            
            Edss[8] = 0;
            Edss[9] = 0.0;
            Edss[10] = 0;
            Edss[11] = 0;
            
            Edss[12] = 0;
            Edss[13] = 0;
            Edss[14] = 0;
            Edss[15] = 0;
            
            /* Gyroscope A matrix */
            A[0] = 0;
            A[1] = ((u[1] + -x[5] * parameters[10]) + x[7] * (parameters[5] * x[4] * cos(x[1]) + x[6] * ((parameters[1] + parameters[2]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * cos(2.0 * x[1]))) * cos(x[2])) + -(parameters[5] * x[4] * x[6] + -(((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(x[1]) * (x[6] * x[6] + -(x[7] * x[7]) * (cos(x[2]) * cos(x[2])))) * sin(x[1]);
            A[2] = 0;
            A[3] = 0;
            
            
            /* Solve system of linear equations Ex = A */
            double l11Sol;
            l11Sol = (1/Edss[5])*A[1];
            
            /* Declare Outputs */
            dxdt[0] = x[4];                 //vel1
            dxdt[1] = x[5];                 //vel2
            dxdt[2] = x[6];                 //vel3
            dxdt[3] = x[7];                 //vel4
            dxdt[4] = 0;                    //acc1
            dxdt[5] = l11Sol;               //acc2
            dxdt[6] = 0;                    //acc3
            dxdt[7] = 0;                    //acc4
            break;
        }
        case 12:
        {
            // q1=locked; q2=locked; q3=free; q4=free;
            
            /* Gyroscope -E matrix */
            Edss[0] = 0;
            Edss[1] = 0.0;
            Edss[2] = 0;
            Edss[3] = 0;
            
            Edss[4] = 0.0;
            Edss[5] = 0;
            Edss[6] = 0.0;
            Edss[7] = 0;
            
            Edss[8] = 0;
            Edss[9] = 0.0;
            Edss[10] = (-1)*(((-parameters[3] + -parameters[4]) + -parameters[5]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * (sin(x[1]) * sin(x[1])));
            Edss[11] = (-1)*((((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(x[1]) * cos(x[2]) * sin(x[1]));
            
            Edss[12] = 0;
            Edss[13] = 0;
            Edss[14] = (-1)*((((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(x[1]) * cos(x[2]) * sin(x[1]));
            Edss[15] = (-1)*(((((-parameters[2] + -parameters[6]) + -parameters[7]) + -parameters[8]) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * (cos(x[2]) * cos(x[2])) * (sin(x[1]) * sin(x[1]))) + (((-parameters[0] + -parameters[1]) + parameters[7]) + parameters[8]) * (sin(x[2]) * sin(x[2])));
            
            /* Gyroscope A matrix */
            A[0] = 0;
            A[1] = 0;
            A[2] = (((u[2] + -x[6] * parameters[11]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * x[5] * x[6] * sin(2.0 * x[1])) + parameters[5] * x[4] * sin(x[1]) * (x[5] + -x[7] * sin(x[2]))) + 0.5 * x[7] * cos(x[2]) * ((-2.0 * (parameters[1] + parameters[2]) * x[5] + ((((((2.0 * parameters[0] + 2.0 * parameters[1]) + parameters[2]) + -parameters[4]) + -parameters[5]) + -2.0 * parameters[7]) + -parameters[8]) * x[7] * sin(x[2])) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(2.0 * x[1]) * (2.0 * x[5] + -x[7] * sin(x[2])));
            A[3] = (((u[3] + -x[7] * parameters[12]) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * x[5] * x[7] * (cos(x[2]) * cos(x[2])) * sin(2.0 * x[1])) + x[6] * (parameters[5] * x[4] + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * x[6] * cos(x[1])) * sin(x[1]) * sin(x[2])) + cos(x[2]) * ((((parameters[1] + parameters [2]) * x[5] * x[6] + -parameters[5] * x[4] * x[5] * cos(x[1])) + ((((((-2.0 * parameters[0] + -2.0 * parameters[1]) + -parameters[2]) + parameters[4]) + parameters[5]) + 2.0 * parameters[7]) + parameters[8]) * x[6] * x[7] * sin(x[2])) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * x[6] * cos(2.0 * x[1]) * (x[5] + x[7] * sin(x[2])));
            
            
            /* Solve system of linear equations Ex = A */
            Matrix2d E12breve;
            Vector2d A12breve;
            Vector2d l12Sol;
            E12breve << Edss[10],Edss[11],
            Edss[14],Edss[15];
            A12breve << A[2],A[3];
            l12Sol = E12breve.partialPivLu().solve(A12breve);
            
            /* Declare Outputs */
            dxdt[0] = x[4];                 //vel1
            dxdt[1] = x[5];                 //vel2
            dxdt[2] = x[6];                 //vel3
            dxdt[3] = x[7];                 //vel4
            dxdt[4] = 0;                    //acc1
            dxdt[5] = 0;                    //acc2
            dxdt[6] = *(l12Sol.data());     //acc3
            dxdt[7] = *(l12Sol.data()+1);   //acc4
            break;
        }
        case 13:
        {
            // q1=locked; q2=locked; q3=free; q4=locked;
            
            /* Gyroscope -E matrix */
            Edss[0] = 0;
            Edss[1] = 0.0;
            Edss[2] = 0;
            Edss[3] = 0;
            
            Edss[4] = 0.0;
            Edss[5] = 0;
            Edss[6] = 0.0;
            Edss[7] = 0;
            
            Edss[8] = 0;
            Edss[9] = 0.0;
            Edss[10] = (-1)*(((-parameters[3] + -parameters[4]) + -parameters[5]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * (sin(x[1]) * sin(x[1])));
            Edss[11] = 0;
            
            Edss[12] = 0;
            Edss[13] = 0;
            Edss[14] = 0;
            Edss[15] = 0;
            
            /* Gyroscope A matrix */
            A[0] = 0;
            A[1] = 0;
            A[2] = (((u[2] + -x[6] * parameters[11]) + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * x[5] * x[6] * sin(2.0 * x[1])) + parameters[5] * x[4] * sin(x[1]) * (x[5] + -x[7] * sin(x[2]))) + 0.5 * x[7] * cos(x[2]) * ((-2.0 * (parameters[1] + parameters[2]) * x[5] + ((((((2.0 * parameters[0] + 2.0 * parameters[1]) + parameters[2]) + -parameters[4]) + -parameters[5]) + -2.0 * parameters[7]) + -parameters[8]) * x[7] * sin(x[2])) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * cos(2.0 * x[1]) * (2.0 * x[5] + -x[7] * sin(x[2])));
            A[3] = 0;
            
            
            /* Solve system of linear equations Ex = A */
            double l13Sol;
            l13Sol = (1/Edss[10])*A[2];
            
            /* Declare Outputs */
            dxdt[0] = x[4];                 //vel1
            dxdt[1] = x[5];                 //vel2
            dxdt[2] = x[6];                 //vel3
            dxdt[3] = x[7];                 //vel4
            dxdt[4] = 0;                    //acc1
            dxdt[5] = 0;                    //acc2
            dxdt[6] = l13Sol;               //acc3
            dxdt[7] = 0;                    //acc4
            break;
        }
        case 14:
        {
            // q1=locked; q2=locked; q3=locked; q4=free;
            
            /* Gyroscope -E matrix */
            Edss[0] = 0;
            Edss[1] = 0.0;
            Edss[2] = 0;
            Edss[3] = 0;
            
            Edss[4] = 0.0;
            Edss[5] = 0;
            Edss[6] = 0.0;
            Edss[7] = 0;
            
            Edss[8] = 0;
            Edss[9] = 0.0;
            Edss[10] = 0;
            Edss[11] = 0;
            
            Edss[12] = 0;
            Edss[13] = 0;
            Edss[14] = 0;
            Edss[15] = (-1)*(((((-parameters[2] + -parameters[6]) + -parameters[7]) + -parameters[8]) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * (cos(x[2]) * cos(x[2])) * (sin(x[1]) * sin(x[1]))) + (((-parameters[0] + -parameters[1]) + parameters[7]) + parameters[8]) * (sin(x[2]) * sin(x[2])));
            
            /* Gyroscope A matrix */
            A[0] = 0;
            A[1] = 0;
            A[2] = 0;
            A[3] = (((u[3] + -x[7] * parameters[12]) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * x[5] * x[7] * (cos(x[2]) * cos(x[2])) * sin(2.0 * x[1])) + x[6] * (parameters[5] * x[4] + (((-parameters[2] + parameters[4]) + parameters[5]) + -parameters[8]) * x[6] * cos(x[1])) * sin(x[1]) * sin(x[2])) + cos(x[2]) * ((((parameters[1] + parameters [2]) * x[5] * x[6] + -parameters[5] * x[4] * x[5] * cos(x[1])) + ((((((-2.0 * parameters[0] + -2.0 * parameters[1]) + -parameters[2]) + parameters[4]) + parameters[5]) + 2.0 * parameters[7]) + parameters[8]) * x[6] * x[7] * sin(x[2])) + (((parameters[2] + -parameters[4]) + -parameters[5]) + parameters[8]) * x[6] * cos(2.0 * x[1]) * (x[5] + x[7] * sin(x[2])));
            
            
            /* Solve system of linear equations Ex = A */
            double l14Sol;
            l14Sol = (1/Edss[15])*A[3];
            
            /* Declare Outputs */
            dxdt[0] = x[4];                 //vel1
            dxdt[1] = x[5];                 //vel2
            dxdt[2] = x[6];                 //vel3
            dxdt[3] = x[7];                 //vel4
            dxdt[4] = 0;                    //acc1
            dxdt[5] = 0;                    //acc2
            dxdt[6] = 0;                    //acc3
            dxdt[7] = l14Sol;               //acc4
            break;
        }
        case 15:
        {
            // q1=locked; q2=locked; q3=locked; q4=locked;
            /* Declare Outputs */
            dxdt[0] = x[4];                 //vel1
            dxdt[1] = x[5];                 //vel2
            dxdt[2] = x[6];                 //vel3
            dxdt[3] = x[7];                 //vel4
            dxdt[4] = 0;                    //acc1
            dxdt[5] = 0;                    //acc2
            dxdt[6] = 0;                    //acc3
            dxdt[7] = 0;                    //acc4
            break;
        }
        default:
        {
            // default case ??
            /* Declare Outputs */
            dxdt[0] = x[4];                 //vel1
            dxdt[1] = x[5];                 //vel2
            dxdt[2] = x[6];                 //vel3
            dxdt[3] = x[7];                 //vel4
            dxdt[4] = 0;                    //acc1
            dxdt[5] = 0;                    //acc2
            dxdt[6] = 0;                    //acc3
            dxdt[7] = 0;                    //acc4
            break;
        }
    }
}

/* Function: mdlTerminate =====================================================
 * Abstract:
 *    No termination needed, but we are required to have this routine.
 */
static void mdlTerminate(SimStruct *S)
{
    UNUSED_ARG(S); /* unused input argument */
}

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
