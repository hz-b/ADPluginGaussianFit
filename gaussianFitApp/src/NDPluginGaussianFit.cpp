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

template <typename epicsType> asynStatus NDPluginGaussianFit::doComputeParamsT(NDArray *pArray, NDGaussianFit_t * pFittedParams)
{
    epicsType *pData = (epicsType *)pArray->pData;

    // Add Residual to ceres problem
    ceres::Problem problem;

    for (size_t i = 0; i < pFittedParams->x_data.size(); i++)
    {
        problem.AddResidualBlock(
            new ceres::AutoDiffCostFunction<GaussianResidual, 1, 1, 1, 1, 1, 1>(
                new GaussianResidual(pFittedParams->x_data[i], pFittedParams->y_data[i], (double)pData[i])),
            nullptr,
            &(pFittedParams->a), &(pFittedParams->x_0), &(pFittedParams->s_x), &(pFittedParams->y_0), &(pFittedParams->s_y));
    }

    // Configure solver options
    ceres::Solver::Options opts;
    opts.linear_solver_type = ceres::DENSE_QR;
    opts.use_nonmonotonic_steps = false;
    opts.minimizer_progress_to_stdout = false;

    opts.max_num_iterations = pFittedParams->max_num_iterations;
    opts.num_threads = pFittedParams->num_threads;

    this->unlock();

    // Ceres Solver execution
    ceres::Solver::Summary summary;
    ceres::Solve(opts, &problem, &summary);

    this->lock();

    pFittedParams->termination_type = summary.termination_type;
    pFittedParams->total_time_in_seconds = summary.total_time_in_seconds;

    return(asynSuccess);

}

asynStatus NDPluginGaussianFit::doComputeParams(NDArray *pArray, NDGaussianFit_t * pFittedParams)
{
    asynStatus status;

    switch(pArray->dataType) {
        case NDInt8:
            status = doComputeParamsT<epicsInt8>(pArray, pFittedParams);
            break;
        case NDUInt8:
            status = doComputeParamsT<epicsUInt8>(pArray, pFittedParams);
            break;
        case NDInt16:
            status = doComputeParamsT<epicsInt16>(pArray, pFittedParams);
            break;
        case NDUInt16:
            status = doComputeParamsT<epicsUInt16>(pArray, pFittedParams);
            break;
        case NDInt32:
            status = doComputeParamsT<epicsInt32>(pArray, pFittedParams);
            break;
        case NDUInt32:
            status = doComputeParamsT<epicsUInt32>(pArray, pFittedParams);
            break;
        case NDInt64:
            status = doComputeParamsT<epicsInt64>(pArray, pFittedParams);
            break;
        case NDUInt64:
            status = doComputeParamsT<epicsUInt64>(pArray, pFittedParams);
            break;
        case NDFloat32:
            status = doComputeParamsT<epicsFloat32>(pArray, pFittedParams);
            break;
        case NDFloat64:
            status = doComputeParamsT<epicsFloat64>(pArray, pFittedParams);
            break;
        default:
            status = asynError;
        break;
    }
    return(status);
}

void NDPluginGaussianFit::processCallbacks(NDArray *pArray)
{

    // Instantiate new struct
    NDGaussianFit_t fittedParams;

    // Create pointer to the struct
    NDGaussianFit_t * pFittedParams = &fittedParams;

    static const char *functionName = "processCallbacks";

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

    // Create linspace arrays for axis x,y
    for (size_t i = 0; i < sizeY; i++)
    {
        for (size_t j = 0; j < sizeX; j++)
        {
            pFittedParams->x_data.push_back((double)j);
            pFittedParams->y_data.push_back((double)i);
        }
    }

    epicsInt32 overwriteInitials = 0;
    getIntegerParam(NDPluginGaussianFitOverwriteInit, &overwriteInitials);

    const bool useInit = justEnabled || overwriteInitials == 0;

    if(useInit)
    {
        if(justEnabled) asynPrint(this->pasynUserSelf, ASYN_TRACE_FLOW, "%s::%s Just Enabled\n", driverName, functionName);
        // Initial guesses for parameters
        getDoubleParam(NDPluginGaussianFitInitAmplitude, &pFittedParams->a);
        getDoubleParam(NDPluginGaussianFitInitX, &pFittedParams->x_0);
        getDoubleParam(NDPluginGaussianFitInitY, &pFittedParams->y_0);
        getDoubleParam(NDPluginGaussianFitInitSx, &pFittedParams->s_x);
        getDoubleParam(NDPluginGaussianFitInitSy, &pFittedParams->s_y);
        asynPrint(this->pasynUserSelf, ASYN_TRACE_FLOW, "%s::%s Overwrite Off\n", driverName, functionName);
    }
    else 
    {
        getDoubleParam(NDPluginGaussianFitAmplitude, &pFittedParams->a);
        getDoubleParam(NDPluginGaussianFitX, &pFittedParams->x_0);
        getDoubleParam(NDPluginGaussianFitY, &pFittedParams->y_0);
        getDoubleParam(NDPluginGaussianFitSx, &pFittedParams->s_x);
        getDoubleParam(NDPluginGaussianFitSy, &pFittedParams->s_y);
        asynPrint(this->pasynUserSelf, ASYN_TRACE_FLOW, "%s::%s Overwrite On\n", driverName, functionName);
    }

    asynPrint(this->pasynUserSelf, ASYN_TRACE_FLOW, "%s::%s Initial Params:\nAmplitude: %f\nX0: %f\n Y0: %f\nSx: %f\nSy: %f\n", driverName, functionName,
        pFittedParams->a, pFittedParams->x_0, pFittedParams->y_0, pFittedParams->s_x, pFittedParams->s_y);

    // Reset after first use
    justEnabled = false;

    // options for ceres solver
    getIntegerParam(NDPluginGaussianFitOptMaxNumIterations, &pFittedParams->max_num_iterations);
    getIntegerParam(NDPluginGaussianFitOptNumThreads, &pFittedParams->num_threads);

    doComputeParams(pArray, pFittedParams);

    // Writing fitted params to PVs
    setDoubleParam(NDPluginGaussianFitAmplitude, pFittedParams->a);
    setDoubleParam(NDPluginGaussianFitX, pFittedParams->x_0);
    setDoubleParam(NDPluginGaussianFitY, pFittedParams->y_0);
    setDoubleParam(NDPluginGaussianFitSx, pFittedParams->s_x);
    setDoubleParam(NDPluginGaussianFitSy, pFittedParams->s_y);

    // Writing summary results to PVs
    setIntegerParam(NDPluginGaussianFitTermination, pFittedParams->termination_type);
    setDoubleParam(NDPluginGaussianFitTotalTime, pFittedParams->total_time_in_seconds);

    asynPrint(this->pasynUserSelf, ASYN_TRACE_FLOW, "%s::%s Fitted Params:\nAmplitude: %f\nX0: %f\n Y0: %f\nSx: %f\nSy: %f\n", driverName, functionName,
        pFittedParams->a, pFittedParams->x_0, pFittedParams->y_0, pFittedParams->s_x, pFittedParams->s_y);

    callParamCallbacks();

    NDPluginDriver::endProcessCallbacks(pArray, true, true);
}

asynStatus NDPluginGaussianFit::writeInt32(asynUser* pasynUser, epicsInt32 value)
{
    int function = pasynUser->reason;

    if (function == NDPluginDriverEnableCallbacks) {
        int oldVal = 0;
        getIntegerParam(NDPluginDriverEnableCallbacks, &oldVal);

        // Let base class do the actual enable/disable bookkeeping
        asynStatus status = NDPluginDriver::writeInt32(pasynUser, value);

        if (oldVal == 0 && value == 1) {
            // plugin has JUST been enabled
            justEnabled = true;
        }
        return status;
    }

    return NDPluginDriver::writeInt32(pasynUser, value);
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

    justEnabled = false;

    char versionString[25];

    setStringParam(NDPluginDriverPluginType, driverName);
    epicsSnprintf(versionString, sizeof(versionString), "%d.%d.%d", GAUSSIAN_FIT_VERSION, GAUSSIAN_FIT_REVISION, GAUSSIAN_FIT_MODIFICATION);
    setStringParam(NDDriverVersion, versionString);

    createParam(NDPluginGaussianFitInitAmplitudeString, asynParamFloat64, &NDPluginGaussianFitInitAmplitude);
    createParam(NDPluginGaussianFitInitXString, asynParamFloat64, &NDPluginGaussianFitInitX);
    createParam(NDPluginGaussianFitInitYString, asynParamFloat64, &NDPluginGaussianFitInitY);
    createParam(NDPluginGaussianFitInitSxString, asynParamFloat64, &NDPluginGaussianFitInitSx);
    createParam(NDPluginGaussianFitInitSyString, asynParamFloat64, &NDPluginGaussianFitInitSy);

    createParam(NDPluginGaussianFitOptMaxNumIterationsString, asynParamInt32, &NDPluginGaussianFitOptMaxNumIterations);
    createParam(NDPluginGaussianFitOptNumThreadsString, asynParamInt32, &NDPluginGaussianFitOptNumThreads);

    createParam(NDPluginGaussianFitAmplitudeString, asynParamFloat64, &NDPluginGaussianFitAmplitude);
    createParam(NDPluginGaussianFitXString, asynParamFloat64, &NDPluginGaussianFitX);
    createParam(NDPluginGaussianFitYString, asynParamFloat64, &NDPluginGaussianFitY);
    createParam(NDPluginGaussianFitSxString, asynParamFloat64, &NDPluginGaussianFitSx);
    createParam(NDPluginGaussianFitSyString, asynParamFloat64, &NDPluginGaussianFitSy);

    createParam(NDPluginGaussianFitOverwriteInitString, asynParamInt32, &NDPluginGaussianFitOverwriteInit);

    createParam(NDPluginGaussianFitTerminationString, asynParamInt32, &NDPluginGaussianFitTermination);
    createParam(NDPluginGaussianFitTotalTimeString, asynParamFloat64, &NDPluginGaussianFitTotalTime);

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