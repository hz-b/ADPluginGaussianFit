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


class NDPluginGaussianFit : public NDPluginDriver
{
private:
    /* data */
public:
    NDPluginGaussianFit(const char *portName, int queueSize, int blockingCallbacks,
                        const char *NDArrayPort, int NDArrayAddr, int maxBuffers,
                        size_t maxMemory, int priority, int stackSize);
    ~NDPluginGaussianFit();
};

void processCallbacks(NDArray *pArray);



#endif