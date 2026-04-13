/**
 * NDPluginGaussianFit.cpp
 *
 * EPICS AreaDetector 2D Gaussian Fit Plugin
 * Main author: Luca Porzio <luca.porzio@helmholtz-berlin.de>
 *
 * Created on: 16/08/2023
 */

/* Standard Includes */
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

/* EPICS Includes */
#include <epicsMutex.h>
#include <iocsh.h>
#include <epicsExport.h>
#include <epicsString.h>

/* User Includes */
#include "NDPluginGaussianFit.h"
#include "ceres/covariance.h"

using namespace std;

static const char *driverName = "NDPluginGaussianFit";
static const double radiansToDegrees = 180.0 / M_PI;
static const double degreesToRadians = M_PI / 180.0;

static double degreesOrSentinel(double radiansValue)
{
    return radiansValue >= 0.0 ? radiansValue * radiansToDegrees : radiansValue;
}

/* Define the residual struct for fitting a 2D Gaussian */
struct RotatedGaussianResidual
{
    RotatedGaussianResidual(double x, double y, double z, double pixelNoise)
        : x_(x), y_(y), z_(z), noise_(pixelNoise) {}

    /*  Params:
        a = amplitude
        x_0 = x position of center
        s_x = std deviation on x
        y_0 = y position of center
        s_y = std deviation on y
        bkg = constant background offset
        theta = rotation angle in radians
    */
    template <typename T> bool operator()(const T *const a, const T *const x_0, const T *const s_x, const T *const y_0, const T *const s_y, const T *const bkg, const T *const theta, T *residual) const
    {
        const T dx = T(x_) - x_0[0];
        const T dy = T(y_) - y_0[0];
        const T c = ceres::cos(theta[0]);
        const T s = ceres::sin(theta[0]);
        const T dxPrime = c * dx + s * dy;
        const T dyPrime = -s * dx + c * dy;
        const T arg = (dxPrime * dxPrime) / (T(2.0) * s_x[0] * s_x[0]) + (dyPrime * dyPrime) / (T(2.0) * s_y[0] * s_y[0]);
        const T model = bkg[0] + a[0] * ceres::exp(-arg);
        residual[0] = (T(z_) - model) / T(noise_);
        return true;
    }

private:
    double x_;
    double y_;
    double z_;
    double noise_;
};

template <typename epicsType> asynStatus NDPluginGaussianFit::doComputeParamsT(NDArray *pArray, NDGaussianFit_t * pFittedParams)
{
    epicsType *pData = (epicsType *)pArray->pData;

    // Add Residual to ceres problem
    ceres::Problem problem;
    ceres::LossFunction *lossFunction = new ceres::HuberLoss(1.0);
    double pixelNoise;
    getDoubleParam(NDPluginGaussianFitPixelNoise, &pixelNoise);

    const size_t numPixels = pFittedParams->x_data.size();
    double minPixel = 0.0;
    if (numPixels > 0)
    {
        minPixel = static_cast<double>(pData[0]);
        for (size_t i = 1; i < numPixels; i++)
        {
            minPixel = std::min(minPixel, static_cast<double>(pData[i]));
        }
    }

    if (!std::isfinite(pFittedParams->background))
    {
        pFittedParams->background = minPixel;
    }

    for (size_t i = 0; i < numPixels; i++)
    {
        problem.AddResidualBlock(
            new ceres::AutoDiffCostFunction<RotatedGaussianResidual, 1, 1, 1, 1, 1, 1, 1, 1>(new RotatedGaussianResidual(
                pFittedParams->x_data[i], pFittedParams->y_data[i], static_cast<double>(pData[i]), pixelNoise)),
            lossFunction,
            &(pFittedParams->a), &(pFittedParams->x_0), &(pFittedParams->s_x), &(pFittedParams->y_0), &(pFittedParams->s_y), &(pFittedParams->background), &(pFittedParams->theta));
    }

    problem.SetParameterLowerBound(&(pFittedParams->s_x), 0, 1e-6);
    problem.SetParameterLowerBound(&(pFittedParams->s_y), 0, 1e-6);
    problem.SetParameterLowerBound(&(pFittedParams->theta), 0, -M_PI / 2.0);
    problem.SetParameterUpperBound(&(pFittedParams->theta), 0, M_PI / 2.0);

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

    double residualSqSum = 0.0;
    double imageSqSum = 0.0;
    for (size_t i = 0; i < numPixels; i++)
    {
        const double dx = pFittedParams->x_data[i] - pFittedParams->x_0;
        const double dy = pFittedParams->y_data[i] - pFittedParams->y_0;
        const double c = std::cos(pFittedParams->theta);
        const double s = std::sin(pFittedParams->theta);
        const double dxPrime = c * dx + s * dy;
        const double dyPrime = -s * dx + c * dy;
        const double arg = (dxPrime * dxPrime) / (2.0 * pFittedParams->s_x * pFittedParams->s_x) + (dyPrime * dyPrime) / (2.0 * pFittedParams->s_y * pFittedParams->s_y);
        const double model = pFittedParams->background + pFittedParams->a * std::exp(-arg);
        const double actual = static_cast<double>(pData[i]);
        const double residual = actual - model;

        residualSqSum += residual * residual;
        imageSqSum += actual * actual;
    }

    pFittedParams->nrmse = imageSqSum > 0.0 ? std::sqrt(residualSqSum / imageSqSum) : 1.0;

    auto setCovarianceErrorValues = [&](double value) {
        pFittedParams->a_err = value;
        pFittedParams->x_0_err = value;
        pFittedParams->y_0_err = value;
        pFittedParams->s_x_err = value;
        pFittedParams->s_y_err = value;
        pFittedParams->theta_err = value;
    };

    auto sigmaFromVariance = [](double variance) {
        return variance >= 0.0 ? std::sqrt(variance) : -1.0;
    };

    if (summary.IsSolutionUsable())
    {
        ceres::Covariance::Options covOptions;
        covOptions.algorithm_type = ceres::DENSE_SVD;
        // Rotated Gaussians can become numerically rank-deficient, e.g. for
        // nearly circular spots or weakly constrained angles. Let Ceres drop
        // tiny singular directions and compute a Moore-Penrose pseudoinverse.
        covOptions.null_space_rank = -1;
        ceres::Covariance covariance(covOptions);

        std::vector<std::pair<const double *, const double *> > covarianceBlocks;
        covarianceBlocks.push_back(std::make_pair(&(pFittedParams->a), &(pFittedParams->a)));
        covarianceBlocks.push_back(std::make_pair(&(pFittedParams->x_0), &(pFittedParams->x_0)));
        covarianceBlocks.push_back(std::make_pair(&(pFittedParams->y_0), &(pFittedParams->y_0)));
        covarianceBlocks.push_back(std::make_pair(&(pFittedParams->s_x), &(pFittedParams->s_x)));
        covarianceBlocks.push_back(std::make_pair(&(pFittedParams->s_y), &(pFittedParams->s_y)));
        covarianceBlocks.push_back(std::make_pair(&(pFittedParams->theta), &(pFittedParams->theta)));

        if (covariance.Compute(covarianceBlocks, &problem))
        {
            double varA = -1.0;
            double varX0 = -1.0;
            double varY0 = -1.0;
            double varSx = -1.0;
            double varSy = -1.0;
            double varTheta = -1.0;

            covariance.GetCovarianceBlock(&(pFittedParams->a), &(pFittedParams->a), &varA);
            covariance.GetCovarianceBlock(&(pFittedParams->x_0), &(pFittedParams->x_0), &varX0);
            covariance.GetCovarianceBlock(&(pFittedParams->y_0), &(pFittedParams->y_0), &varY0);
            covariance.GetCovarianceBlock(&(pFittedParams->s_x), &(pFittedParams->s_x), &varSx);
            covariance.GetCovarianceBlock(&(pFittedParams->s_y), &(pFittedParams->s_y), &varSy);
            covariance.GetCovarianceBlock(&(pFittedParams->theta), &(pFittedParams->theta), &varTheta);

            pFittedParams->a_err = sigmaFromVariance(varA);
            pFittedParams->x_0_err = sigmaFromVariance(varX0);
            pFittedParams->y_0_err = sigmaFromVariance(varY0);
            pFittedParams->s_x_err = sigmaFromVariance(varSx);
            pFittedParams->s_y_err = sigmaFromVariance(varSy);
            pFittedParams->theta_err = sigmaFromVariance(varTheta);
        }
        else
        {
            setCovarianceErrorValues(-1.0);
        }
    }
    else
    {
        setCovarianceErrorValues(-1.0);
    }

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
        if(justEnabled) asynPrint(this->pasynUserSelf, ASYN_TRACE_FLOW, "%s::%s Just Enabled. Using User Defined Initial Guesses\n", driverName, functionName);
        // Initial guesses for parameters
        getDoubleParam(NDPluginGaussianFitInitAmplitude, &pFittedParams->a);
        getDoubleParam(NDPluginGaussianFitInitX, &pFittedParams->x_0);
        getDoubleParam(NDPluginGaussianFitInitY, &pFittedParams->y_0);
        getDoubleParam(NDPluginGaussianFitInitSx, &pFittedParams->s_x);
        getDoubleParam(NDPluginGaussianFitInitSy, &pFittedParams->s_y);
        pFittedParams->background = std::numeric_limits<double>::quiet_NaN();
        pFittedParams->theta = 0.0;
        asynPrint(this->pasynUserSelf, ASYN_TRACE_FLOW, "%s::%s Overwrite Off\n", driverName, functionName);
    }
    else 
    {
        getDoubleParam(NDPluginGaussianFitAmplitude, &pFittedParams->a);
        getDoubleParam(NDPluginGaussianFitX, &pFittedParams->x_0);
        getDoubleParam(NDPluginGaussianFitY, &pFittedParams->y_0);
        getDoubleParam(NDPluginGaussianFitSx, &pFittedParams->s_x);
        getDoubleParam(NDPluginGaussianFitSy, &pFittedParams->s_y);
        getDoubleParam(NDPluginGaussianFitBackground, &pFittedParams->background);
        getDoubleParam(NDPluginGaussianFitTheta, &pFittedParams->theta);
        pFittedParams->theta *= degreesToRadians;
        asynPrint(this->pasynUserSelf, ASYN_TRACE_FLOW, "%s::%s Overwrite On\n", driverName, functionName);
    }

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
    setDoubleParam(NDPluginGaussianFitBackground, pFittedParams->background);
    setDoubleParam(NDPluginGaussianFitTheta, pFittedParams->theta * radiansToDegrees);
    setDoubleParam(NDPluginGaussianFitAmplitudeErr, pFittedParams->a_err);
    setDoubleParam(NDPluginGaussianFitXErr, pFittedParams->x_0_err);
    setDoubleParam(NDPluginGaussianFitYErr, pFittedParams->y_0_err);
    setDoubleParam(NDPluginGaussianFitSxErr, pFittedParams->s_x_err);
    setDoubleParam(NDPluginGaussianFitSyErr, pFittedParams->s_y_err);
    setDoubleParam(NDPluginGaussianFitThetaErr, degreesOrSentinel(pFittedParams->theta_err));
    setDoubleParam(NDPluginGaussianFitNRMSE, pFittedParams->nrmse);

    // Writing summary results to PVs
    setIntegerParam(NDPluginGaussianFitTermination, pFittedParams->termination_type);
    setDoubleParam(NDPluginGaussianFitTotalTime, pFittedParams->total_time_in_seconds);

    asynPrint(this->pasynUserSelf, ASYN_TRACE_FLOW, "%s::%s Fitted Params:\nAmplitude: %f\nX0: %f\nY0: %f\nSx: %f\nSy: %f\nTheta(deg): %f\nBackground: %f\nNRMSE: %f\n", driverName, functionName,
        pFittedParams->a, pFittedParams->x_0, pFittedParams->y_0, pFittedParams->s_x, pFittedParams->s_y, pFittedParams->theta * radiansToDegrees, pFittedParams->background, pFittedParams->nrmse);

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
    createParam(NDPluginGaussianFitPixelNoiseString, asynParamFloat64, &NDPluginGaussianFitPixelNoise);

    createParam(NDPluginGaussianFitAmplitudeString, asynParamFloat64, &NDPluginGaussianFitAmplitude);
    createParam(NDPluginGaussianFitXString, asynParamFloat64, &NDPluginGaussianFitX);
    createParam(NDPluginGaussianFitYString, asynParamFloat64, &NDPluginGaussianFitY);
    createParam(NDPluginGaussianFitSxString, asynParamFloat64, &NDPluginGaussianFitSx);
    createParam(NDPluginGaussianFitSyString, asynParamFloat64, &NDPluginGaussianFitSy);
    createParam(NDPluginGaussianFitBackgroundString, asynParamFloat64, &NDPluginGaussianFitBackground);
    createParam(NDPluginGaussianFitThetaString, asynParamFloat64, &NDPluginGaussianFitTheta);
    createParam(NDPluginGaussianFitAmplitudeErrString, asynParamFloat64, &NDPluginGaussianFitAmplitudeErr);
    createParam(NDPluginGaussianFitXErrString, asynParamFloat64, &NDPluginGaussianFitXErr);
    createParam(NDPluginGaussianFitYErrString, asynParamFloat64, &NDPluginGaussianFitYErr);
    createParam(NDPluginGaussianFitSxErrString, asynParamFloat64, &NDPluginGaussianFitSxErr);
    createParam(NDPluginGaussianFitSyErrString, asynParamFloat64, &NDPluginGaussianFitSyErr);
    createParam(NDPluginGaussianFitThetaErrString, asynParamFloat64, &NDPluginGaussianFitThetaErr);
    createParam(NDPluginGaussianFitNRMSEString, asynParamFloat64, &NDPluginGaussianFitNRMSE);

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
