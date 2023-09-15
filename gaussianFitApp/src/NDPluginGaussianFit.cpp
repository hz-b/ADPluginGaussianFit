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
#include <epicsExport.h>
#include <epicsString.h>

/* User Includes */
#include "NDPluginGaussianFit.h"

using namespace std;

static const char *driverName = "NDPluginGaussianFit";

/* Define the residual struct for fitting a 2D Gaussian */
struct GaussianResidual
{
    GaussianResidual(double x, double y, double z)
        : x_(x), y_(y), z_(z) {}

    /*  Params:
        a = amplitude
        x_0 = x position of center
        s_x = std deviation on x
        y_0 = y position of center
        s_y = std deviation on y
    */
    template <typename T>
    bool operator()(const T *const a, const T *const x_0, const T *const s_x, const T *const y_0, const T *const s_y, T *residual) const
    {
        residual[0] = z_ - a[0] * exp(-((x_ - x_0[0]) * (x_ - x_0[0]) / (2.0 * s_x[0] * s_x[0]) + (y_ - y_0[0]) * (y_ - y_0[0]) / (2.0 * s_y[0] * s_y[0])));
        return true;
    }

private:
    double x_;
    double y_;
    double z_;
};

void NDPluginGaussianFit::processCallbacks(NDArray *pArray)
{
    static const char *functionName = "processCallbacks";

    epicsFloat64 *pData = (epicsFloat64 *)pArray->pData;

    // Call base class and get information about the frame
    NDPluginDriver::beginProcessCallbacks(pArray);

    // Retrieve image size
    size_t sizeX = 0, sizeY = 0;
    if (pArray->ndims > 0)
        sizeX = pArray->dims[0].size;
    if (pArray->ndims == 1)
        sizeY = 1;
    if (pArray->ndims > 1)
        sizeY = pArray->dims[1].size;

    // Unlock before computing the Residual
    this->unlock();

    // Create dimension arrays
    vector<double> x_data, y_data;
    for (size_t i = 0; i < sizeY; i++)
    {
        for (size_t j = 0; j < sizeX; j++)
        {
            x_data.push_back((double)j);
            y_data.push_back((double)i);
        }
    }

    // Initial guesses for parameters
    double a = 200.0;
    double x_0 = 200.0;
    double y_0 = 200.0;
    double s_x = 1.0;
    double s_y = 1.0;

    // Add Residual to ceres problem
    ceres::Problem problem;

    for (size_t i = 0; i < x_data.size(); i++)
    {
        problem.AddResidualBlock(
            new ceres::AutoDiffCostFunction<GaussianResidual, 1, 1, 1, 1, 1, 1>(
                new GaussianResidual(x_data[i], y_data[i], (double)pData[i])),
            nullptr,
            &a, &x_0, &s_x, &y_0, &s_y);
    }

    // Configure solver options
    ceres::Solver::Options opts;
    opts.max_num_iterations = 100;
    opts.num_threads = 8;
    opts.linear_solver_type = ceres::DENSE_QR;
    opts.use_nonmonotonic_steps = false;
    opts.minimizer_progress_to_stdout = false;

    ceres::Solver::Summary summary;
    ceres::Solve(opts, &problem, &summary);

    // Return the lock after computation
    this->lock();

    asynPrint(this->pasynUserSelf, ASYN_TRACE_ERROR, "%s::%s Fitted Params:\nAmplitude: %f\nX0: %f\n Y0: %f\nSx: %f\nSy: %f\n", driverName, functionName,
              a, x_0, y_0, s_x, s_y);

    callParamCallbacks();
}

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

    createParam(NDPluginGaussianFitInitAmplitudeString, asynParamFloat64, &NDPluginGaussianFitInitAmplitude);
    createParam(NDPluginGaussianFitInitXString, asynParamFloat64, &NDPluginGaussianFitInitX);
    createParam(NDPluginGaussianFitInitYString, asynParamFloat64, &NDPluginGaussianFitInitY);
    createParam(NDPluginGaussianFitInitSxString, asynParamFloat64, &NDPluginGaussianFitInitSx);
    createParam(NDPluginGaussianFitInitSyString, asynParamFloat64, &NDPluginGaussianFitInitSy);

    connectToArrayPort();
}

/* External function that is called in the IOC shell to create the plugin object */
extern "C" int NDGaussianFitConfigure(const char *portName, int queueSize, int blockingCallbacks, const char *NDArrayPort,
                                      int NDArrayAddr, int maxBuffers, size_t maxMemory, int priority, int stackSize)
{

    // calls the plugin constructor
    NDPluginGaussianFit *pPlugin = new NDPluginGaussianFit(portName, queueSize, blockingCallbacks, NDArrayPort, NDArrayAddr,
                                                           maxBuffers, maxMemory, priority, stackSize);

    // starts the plugin
    return pPlugin->start();
}

/* IOC shell argument initialization here */
static const iocshArg initArg0 = {"portName", iocshArgString};
static const iocshArg initArg1 = {"frame queue size", iocshArgInt};
static const iocshArg initArg2 = {"blocking callbacks", iocshArgInt};
static const iocshArg initArg3 = {"NDArrayPort", iocshArgString};
static const iocshArg initArg4 = {"NDArrayAddr", iocshArgInt};
static const iocshArg initArg5 = {"maxBuffers", iocshArgInt};
static const iocshArg initArg6 = {"maxMemory", iocshArgInt};
static const iocshArg initArg7 = {"priority", iocshArgInt};
static const iocshArg initArg8 = {"stackSize", iocshArgInt};
static const iocshArg *const initArgs[] = {&initArg0,
                                           &initArg1,
                                           &initArg2,
                                           &initArg3,
                                           &initArg4,
                                           &initArg5,
                                           &initArg6,
                                           &initArg7,
                                           &initArg8};

/* defines the configuration function for Initializing the plugin */
static const iocshFuncDef initFuncDef = {"NDGaussianFitConfigure", 9, initArgs};

/* Init call function for the IOC shell */
static void initCallFunc(const iocshArgBuf *args)
{
    NDGaussianFitConfigure(args[0].sval, args[1].ival, args[2].ival,
                           args[3].sval, args[4].ival, args[5].ival,
                           args[6].ival, args[7].ival, args[8].ival);
}

/* Registration of NDPluginCV for PV autosaving and into the IOC shell command set*/
extern "C" void NDGaussianFitRegister(void)
{
    iocshRegister(&initFuncDef, initCallFunc);
}

/* Extern C function called from IOC startup script */
extern "C"
{
    epicsExportRegistrar(NDGaussianFitRegister);
}