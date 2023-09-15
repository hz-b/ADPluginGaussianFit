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

/* AreaDetector includes */
#include "NDPluginDriver.h"


/* Ceres include */
#include "ceres/ceres.h"

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

protected:
    // Initial guess for Amplitude
    int NDPluginGaussianFitInitAmplitude;
    
    // Initial guess for X
    int NDPluginGaussianFitInitX;
    
    // Initial guess for Y
    int NDPluginGaussianFitInitY;
    
    // Initial guess for SigmaX
    int NDPluginGaussianFitInitSx;
    
    // Initial guess for SigmaY
    int NDPluginGaussianFitInitSy;
};



#endif