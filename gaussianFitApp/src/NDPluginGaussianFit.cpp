/**
 * NDPluginGaussianFit.cpp
 *
 * EPICS AreaDetector 2D Gaussian Fit Plugin
 * Main author: Luca Porzio <luca.porzio@helmholtz-berlin.de>
 *
 * Created on: 16/08/2023
 */

/* Standard Includes */
#include <iostream>
#include <vector>

/* EPICS Includes */
#include <epicsMutex.h>
#include <iocsh.h>

/* User Includes */
#include <NDPluginGaussianFit.h>

using namespace std;

static const char *driverName = "NDPluginGaussianFit";

/* Define the residual class for fitting a 2D Gaussian */
static struct GaussianResidual
{
    GaussianResidual(double x, double y, double z)
        : x_(x), y_(y), z_(z) {}

    template <typename T>
    bool operator()(const T *const params, T *residual) const
    {
        T x_diff = T(x_) - params[0];
        T y_diff = T(y_) - params[1];
        T predicted = params[2] * exp(-x_diff * x_diff / (params[3] * params[3] * 2.0) -
                                      y_diff * y_diff / (params[4] * params[4] * 2.0)) +
                      params[5];
        residual[0] = T(z_) - predicted;
        return true;
    }

private:
    double x_;
    double y_;
    double z_;
};

NDPluginGaussianFit::NDPluginGaussianFit(const char *portName, int queueSize, int blockingCallbacks,
                                         const char *NDArrayPort, int NDArrayAddr, int maxBuffers,
                                         size_t maxMemory, int priority, int stackSize)
    : NDPluginDriver(portName, queueSize, blockingCallbacks,
                     NDArrayPort, NDArrayAddr, 1, maxBuffers, maxMemory,
                     asynInt32ArrayMask | asynFloat64ArrayMask | asynGenericPointerMask,
                     asynInt32ArrayMask | asynFloat64ArrayMask | asynGenericPointerMask,
                     ASYN_MULTIDEVICE, 1, priority, stackSize, 1)
{

    char versionString[25];

    setStringParam(NDPluginDriverPluginType, driverName);
    epicsSnprintf(versionString, sizeof(versionString), "%d.%d.%d", GAUSSIAN_FIT_VERSION, GAUSSIAN_FIT_REVISION, GAUSSIAN_FIT_MODIFICATION);
    setStringParam(NDDriverVersion, versionString);
    connectToArrayPort();
}