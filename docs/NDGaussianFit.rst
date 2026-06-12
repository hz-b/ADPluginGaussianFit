.. _NDGaussianFit:

==============
NDGaussianFit
==============

**Author:** Luca Porzio (Helmholtz-Zentrum Berlin)

---

Overview
========

NDGaussianFit is an EPICS areaDetector plugin that performs 2D Gaussian fitting on NDArray data. It utilizes the `Ceres Solver <http://ceres-solver.org/>`_ library to solve the non-linear least squares optimization problem.

The plugin fits a rotated 2D Gaussian function to the intensity distribution of an input image or region of interest (ROI). This is particularly useful for beam diagnostics, centering, and analyzing beam profiles in scientific applications such as synchrotron radiation beamlines.

The fitted Gaussian function is defined as:

.. math::

    f(x, y) = A \cdot \exp\left(-\frac{(x'')^2}{2\sigma_x^2} - \frac{(y'')^2}{2\sigma_y^2}\right) + B

where:

- :math:`A` is the amplitude (peak intensity)
- :math:`(x_0, y_0)` is the centroid position
- :math:`\sigma_x` and :math:`\sigma_y` are the standard deviations along the x and y axes
- :math:`\theta` is the rotation angle of the ellipse
- :math:`B` is the constant background offset
- :math:`x''` and :math:`y''` are the coordinates after rotation:

  .. math::

      \begin{aligned}
      x'' &= \cos(\theta) \cdot (x - x_0) + \sin(\theta) \cdot (y - y_0) \\
      y'' &= -\sin(\theta) \cdot (x - x_0) + \cos(\theta) \cdot (y - y_0)
      \end{aligned}

The plugin inherits from NDPluginDriver and supports all standard plugin features (callbacks, blocking/non-blocking execution, etc.).

Installation and Configuration
==============================

To install the plugin, follow the instructions provided in the `repository README <https://github.com/ADPluginGaussianFit>`_.

The NDGaussianFit plugin is created with the NDGaussianFitConfigure command, either from C/C++ or from the EPICS IOC shell.

.. code-block:: c

    NDGaussianFitConfigure(const char *portName, const char *queueName, int blockingCallbacks, 
                           const char *NDArrayPort, int NDArrayAddr, int maxBuffers, 
                           size_t maxMemory, int priority, int stackSize)

For details on the meaning of the parameters to this function refer to the detailed documentation on the NDGaussianFitConfigure function in the NDPluginGaussianFit.cpp documentation and in the documentation for the constructor for the NDPluginGaussianFit class.

Example st.cmd Configuration
----------------------------

.. code-block:: bash

    NDGaussianFitConfigure("GAUSS1", 20, 0, "TS1", 0, 0, 0, 0, 0)
    dbLoadRecords("$(ADPLUGINGAUSSIANFIT)/db/NDPluginGaussianFit.template", "P=$(PREFIX),R=GaussianFit1:, PORT=GAUSS1,ADDR=0,TIMEOUT=1,NDARRAY_PORT=TS1")

Parameters and Records
======================

The NDGaussianFit plugin supports the following parameters in addition to those defined in the NDPluginDriver base class.

All r/w parameters have corresponding *_RBV* (readback) records. All r/o parameters are updated automatically during each processing cycle.

### Initial Guesses

| Parameter index variable | asyn interface | Access | Description | drvInfo string | EPICS record name | EPICS record type |
|--------------------------|----------------|--------|-------------|----------------|-------------------|-------------------|
| NDPluginGaussianFitInitAmplitude | asynFloat64 | r/w | Initial guess for amplitude | INITIAL_AMPLITUDE | $(P)$(R)Initial_Amplitude | ao |
| NDPluginGaussianFitInitX | asynFloat64 | r/w | Initial guess for centroid X position | INITIAL_X | $(P)$(R)Initial_X | ao |
| NDPluginGaussianFitInitY | asynFloat64 | r/w | Initial guess for centroid Y position | INITIAL_Y | $(P)$(R)Initial_Y | ao |
| NDPluginGaussianFitInitSx | asynFloat64 | r/w | Initial guess for sigma X | INITIAL_SX | $(P)$(R)Initial_SX | ao |
| NDPluginGaussianFitInitSy | asynFloat64 | r/w | Initial guess for sigma Y | INITIAL_SY | $(P)$(R)Initial_SY | ao |

### Solver Configuration

| Parameter index variable | asyn interface | Access | Description | drvInfo string | EPICS record name | EPICS record type |
|--------------------------|----------------|--------|-------------|----------------|-------------------|-------------------|
| NDPluginGaussianFitOptMaxNumIterations | asynInt32 | r/w | Maximum number of solver iterations | MAX_ITERATIONS | $(P)$(R)Max_Iterations | ao |
| NDPluginGaussianFitOptNumThreads | asynInt32 | r/w | Number of threads used by the solver | CERES_THREADS | $(P)$(R)Ceres_Threads | ao |
| NDPluginGaussianFitPixelNoise | asynFloat64 | r/w | Pixel noise value for normalization (used in threshold computation) | PIXEL_NOISE | $(P)$(R)PixelNoise | ao |

### Fitted Parameters

| Parameter index variable | asyn interface | Access | Description | drvInfo string | EPICS record name | EPICS record type |
|--------------------------|----------------|--------|-------------|----------------|-------------------|-------------------|
| NDPluginGaussianFitAmplitude | asynFloat64 | r/o | Fitted amplitude (peak intensity) | FITTED_AMPLITUDE | $(P)$(R)Amplitude_Fitted | ai |
| NDPluginGaussianFitX | asynFloat64 | r/o | Fitted centroid X position | FITTED_X | $(P)$(R)X_Fitted | ai |
| NDPluginGaussianFitY | asynFloat64 | r/o | Fitted centroid Y position | FITTED_Y | $(P)$(R)Y_Fitted | ai |
| NDPluginGaussianFitSx | asynFloat64 | r/o | Fitted sigma X (standard deviation) | FITTED_SX | $(P)$(R)SX_Fitted | ai |
| NDPluginGaussianFitSy | asynFloat64 | r/o | Fitted sigma Y (standard deviation) | FITTED_SY | $(P)$(R)SY_Fitted | ai |
| NDPluginGaussianFitBackground | asynFloat64 | r/o | Fitted constant background offset | FITTED_BACKGROUND | $(P)$(R)Background_Fitted | ai |
| NDPluginGaussianFitTheta | asynFloat64 | r/o | Fitted rotation angle in degrees | FITTED_THETA | $(P)$(R)Theta_Fitted | ai |

### Parameter Uncertainties

| Parameter index variable | asyn interface | Access | Description | drvInfo string | EPICS record name | EPICS record type |
|--------------------------|----------------|--------|-------------|----------------|-------------------|-------------------|
| NDPluginGaussianFitAmplitudeErr | asynFloat64 | r/o | Uncertainty (1-sigma) for amplitude | FITTED_AMPLITUDE_ERR | $(P)$(R)Amplitude_Err | ai |
| NDPluginGaussianFitXErr | asynFloat64 | r/o | Uncertainty (1-sigma) for centroid X | FITTED_X_ERR | $(P)$(R)X_Err | ai |
| NDPluginGaussianFitYErr | asynFloat64 | r/o | Uncertainty (1-sigma) for centroid Y | FITTED_Y_ERR | $(P)$(R)Y_Err | ai |
| NDPluginGaussianFitSxErr | asynFloat64 | r/o | Uncertainty (1-sigma) for sigma X | FITTED_SX_ERR | $(P)$(R)SX_Err | ai |
| NDPluginGaussianFitSyErr | asynFloat64 | r/o | Uncertainty (1-sigma) for sigma Y | FITTED_SY_ERR | $(P)$(R)SY_Err | ai |
| NDPluginGaussianFitThetaErr | asynFloat64 | r/o | Uncertainty (1-sigma) for rotation angle (degrees) | FITTED_THETA_ERR | $(P)$(R)Theta_Err | ai |

### Fit Quality

| Parameter index variable | asyn interface | Access | Description | drvInfo string | EPICS record name | EPICS record type |
|--------------------------|----------------|--------|-------------|----------------|-------------------|-------------------|
| NDPluginGaussianFitNRMSE | asynFloat64 | r/o | Normalized root mean square error of the fit | NRMSE | $(P)$(R)NRMSE | ai |
| NDPluginGaussianFitTermination | asynInt32 | r/o | Solver termination reason (0=CONVERGENCE, 1=NO_CONVERGENCE, 2=FAILURE, 3=USER_SUCCESS) | TERMINATION | $(P)$(R)Termination | mbbi |
| NDPluginGaussianFitTotalTime | asynFloat64 | r/o | Total solver time in seconds | TOTAL_TIME | $(P)$(R)Total_Time | ai |

### Overwrite Initial Guesses

| Parameter index variable | asyn interface | Access | Description | drvInfo string | EPICS record name | EPICS record type |
|--------------------------|----------------|--------|-------------|----------------|-------------------|-------------------|
| NDPluginGaussianFitOverwriteInit | asynInt32 | r/w | Flag to enable/disable closed loop mode (0=use user initial guesses, 1=use previous fitted values as initial conditions) | OVERWRITE | $(P)$(R)Overwrite_Initials | bo |

Termination Types
-----------------

The TERMINATION record indicates the solver termination reason:

- **0 (CONVERGENCE)**: The solver successfully converged to a solution
- **1 (NO_CONVERGENCE)**: The solver failed to converge within the maximum number of iterations
- **2 (FAILURE)**: The solver encountered an error (e.g., not enough data points)
- **3 (USER_SUCCESS)**: The solver was stopped by user callback

 fit Quality Metric
===================

The plugin provides the NRMSE (Normalized Root Mean Square Error) as a quality metric for the fit:

.. math::

    NRMSE = \sqrt{\frac{\sum_{i=1}^{N} (z_i - f(x_i, y_i))^2}{\sum_{i=1}^{N} z_i^2}}

where:
- :math:`z_i` is the measured pixel value
- :math:`f(x_i, y_i)` is the fitted Gaussian value
- :math:`N` is the total number of pixels

Lower NRMSE values indicate better fits. Typical values range from 0.01 (excellent fit) to 0.1 (acceptable fit).

Uncertainty Computation
=======================

The plugin uses the Ceres Solver covariance matrix to compute parameter uncertainties (1-sigma standard deviations). The covariance is computed using a dense SVD solver with Moore-Penrose pseudoinverse for rank-deficient problems. If the covariance computation fails or the solution is not usable, the uncertainty values are set to -1.0.

Performance Guidelines
======================

For optimal performance:

1. **Restrict to ROI**: Fitting a smaller region of interest significantly improves speed compared to fitting the full frame.

2. **Non-blocking callbacks**: If the fit time exceeds the frame rate, consider setting ``blockingCallbacks=0`` in the configuration or using decimation (via NDPluginProcess) to avoid blocking the driver.

3. **Closed loop mode**: If the input distribution doesn't change much between frames, enable closed loop mode by setting ``Overwrite_Initials=1`` (or selecting "Yes" in the OPI). This uses previous fitted parameters as initial conditions, typically reducing convergence time.

4. **Thread configuration**: Increasing the number of threads can improve performance by 10-15%. Set ``Ceres_Threads`` based on your CPU cores.

5. **Pixel noise**: Set an appropriate ``PixelNoise`` value for your detector. This affects the threshold used to exclude background pixels from the fit.

Example Performance
-------------------

On an Intel i9-12900K CPU with 1 thread, the plugin achieves approximately 10 Hz processing rate for typical 64x64 ROI.

Frequently Asked Questions
==========================

**Q: How do I choose good initial guesses?**

A: For most applications, setting the initial amplitude to the maximum pixel value, initial X/Y to the image center, and initial sigma values to 10-20% of the image size works well. The closed loop mode (Overwrite_Initials=1) can also be used to automatically refine initial guesses.

**Q: What does the rotation angle represent?**

A: The rotation angle :math:`\theta` represents the counterclockwise rotation of the Gaussian ellipse's major axis from the X-axis. It ranges from -90° to +90°.

**Q: Why are some uncertainties -1.0?**

A: Uncertainty values of -1.0 indicate that the covariance computation failed. This can happen if the fit is poorly constrained (e.g., very weak signal, insufficient data points) or if the Hessian matrix is rank-deficient.

**Q: How many pixels are needed for a reliable fit?**

A: The fit uses 7 parameters (A, x0, y0, sx, sy, theta, background). A minimum of 10 pixels above threshold is required. For reliable uncertainties, at least 50-100 pixels are recommended.

