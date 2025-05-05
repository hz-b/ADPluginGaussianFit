/**
 * NDPluginGaussianFit.h
 * 
 * Header file for EPICS AreaDetector 2D Gaussian Fit Plugin
 * Main author: Luca Porzio <luca.porzio@helmholtz-berlin.de>
 * 
 * Created on: 16/08/2023
*/

#ifndef NDPluginGaussianFit_H
#define NDPluginGaussianFit_H

/* Standard includes */
#include <vector>

/* AreaDetector includes */
#include "NDPluginDriver.h"


/* Ceres include */
#include "ceres/ceres.h"

typedef struct NDGaussianFit{

    std::vector<double> x_data;
    std::vector<double> y_data;
    double a;
    double x_0;
    double y_0;
    double s_x;
    double s_y;
    int max_num_iterations;
    int num_threads;
    int termination_type;
    double total_time_in_seconds;

} NDGaussianFit_t;

//version numbers
#define GAUSSIAN_FIT_VERSION 2
#define GAUSSIAN_FIT_REVISION 1
#define GAUSSIAN_FIT_MODIFICATION 0

// Initial guesses
#define NDPluginGaussianFitInitAmplitudeString  "INITIAL_AMPLITUDE"     // (asynFloat64)
#define NDPluginGaussianFitInitXString          "INITIAL_X"             // (asynFloat64)
#define NDPluginGaussianFitInitYString          "INITIAL_Y"             // (asynFloat64)
#define NDPluginGaussianFitInitSxString         "INITIAL_SX"            // (asynFloat64)
#define NDPluginGaussianFitInitSyString         "INITIAL_SY"            // (asynFloat64)

// Solver configuration
#define NDPluginGaussianFitOptMaxNumIterationsString    "MAX_ITERATIONS"        // (asynInt32)
#define NDPluginGaussianFitOptNumThreadsString          "CERES_THREADS"           // (asynInt32)

// Fitted Parameters
#define NDPluginGaussianFitAmplitudeString      "FITTED_AMPLITUDE"      // (asynFloat64)
#define NDPluginGaussianFitXString              "FITTED_X"              // (asynFloat64)
#define NDPluginGaussianFitYString              "FITTED_Y"              // (asynFloat64)
#define NDPluginGaussianFitSxString             "FITTED_SX"             // (asynFloat64)
#define NDPluginGaussianFitSyString             "FITTED_SY"             // (asynFloat64)

// Summary
#define NDPluginGaussianFitTerminationString    "TERMINATION"       // (asynInt32)
#define NDPluginGaussianFitTotalTimeString      "TOTAL_TIME"        // (asynFloat64)


class NDPLUGIN_API NDPluginGaussianFit : public NDPluginDriver
{
private:
    /* data */
public:
    NDPluginGaussianFit(const char *portName, int queueSize, int blockingCallbacks,
                        const char *NDArrayPort, int NDArrayAddr, int maxBuffers,
                        size_t maxMemory, int priority, int stackSize);

    /* Override the virtual method in the base class */
    void processCallbacks(NDArray *pArray);

    template <typename epicsType> asynStatus doComputeParamsT(NDArray *pArray, NDGaussianFit_t * pFittedParams);
    asynStatus doComputeParams(NDArray *pArray, NDGaussianFit_t * pFittedParams);


protected:
    // Initial guess for Amplitude
    int NDPluginGaussianFitInitAmplitude;
    
    // Initial guess for Centroid X-axis
    int NDPluginGaussianFitInitX;
    
    // Initial guess for Centrois Y-axis
    int NDPluginGaussianFitInitY;
    
    // Initial guess for SigmaX
    int NDPluginGaussianFitInitSx;
    
    // Initial guess for SigmaY
    int NDPluginGaussianFitInitSy;

    // Ceres Solver Max Number of Iterations 
    int NDPluginGaussianFitOptMaxNumIterations;

    // Ceres Solver Number of Threads
    int NDPluginGaussianFitOptNumThreads;

    // Fitted value for Amplitude
    int NDPluginGaussianFitAmplitude;

    // Fitted value for Centroid X-axis
    int NDPluginGaussianFitX;

    // Fitted value for Centroid Y-axis
    int NDPluginGaussianFitY;

    // Fitted value for SigmaX
    int NDPluginGaussianFitSx;

    // Fitted value for SigmaY
    int NDPluginGaussianFitSy;

    // Ceres Solver Termination Type
    int NDPluginGaussianFitTermination;

    // Total time spent in Ceres Solver
    int NDPluginGaussianFitTotalTime;
};



#endif