# PhD Thesis: Multi-strange particle production in p–Pb collisions at √s<sub>NN</sub> = 8.16 TeV

- Author: RNDr. Ishaan Ahuja, PhD.
- Supervisor: doc. RNDr. Marek Bombara, PhD.
- DOI: ``https://doi.org/10.17181/cwcde-g1z94``

## File Structure and Documentation

1. GRID:

   - ``AliAnalysisTaskStrangeCascadesRun2`` : main task for code development and analysis on grid. Suited for Xi and Omega strangeness analysis.

   - ``AddTaskStrangeCascadesRun2.C`` : Adds the ``AliAnalysisTaskStrangeCascadesRun2`` task to the AliAnalysisManager. Configure before running task on Grid.

   - ``RunAnalysisTask_StrangeCascadesRun2.C`` : Configures and runs the AliRoot analysis task for analysing strange cascades. Run locally or on ALICE grid.

2. Post-processing:
   After collecting data and MC outputs from grid:

    2.1 ``FitCascades.C`` : Integrated, mult.-integrated and differential Invariant Mass projections and fitting.

    2.2 ``DrawCascades.C`` : Quick drawing and saving to images on local disk for checking Invariant Mass fitting. Optional.

    2.3 ``EfficiencyEstimation.C`` : Estimate efficiency from multiple MC files, calculate average efficiency according to weighted least squares approach (PDG). Plot efficiency comparison plots.

    2.4 ``EfficiencyCorrection.C`` : Efficiency correction on particle spectra. Generates efficiency corrected pt spectra (with statistical errors) from fit results (2.1) and average efficiency estimation (2.3).

    - ``CascadeUtils.h`` : Utilities and helper functions commonly used throughout the project.

3. Systematic Uncertainty Evaluation:

    - Multi Trial:

    3.1 Run ``AliAnalysisTaskStrangeCascadesRun2`` with the flag ``fDefOnly = false``. Should enable random cut variations and repeat the same analysis for the specified ``nVar`` (default=500).

    3.2 Repeat same post-processing framework as in (2.) using automated runner macros for result generation and file handling (e.g. ``RunVariationsFit.C``, ``RunVariationsEff.C``, etc.), preferably on a powerful PC/cluster.

    3.3 ``SysMultiTrial.C`` : Feed efficiency corrected spectra of variations and default cuts to calculate yield deviations, generate Gaussians (with possible Roger Barlow test), extract sys. uncert., and save to file/generate plots.

    - Total Systematic Uncertainty:

    3.4 ``SysSignalExtraction.C`` : Vary signal extraction parameters to calculate systematic uncertainty.

    3.5 ``TotalSystematics.C`` : Combine and calculate total systematics from various sources. Generates relative sys. uncertainty comparison, corrected spectra with stat. + sys. uncertainty, and Levy-Tsallis fits.

4. Integrated Yield and Mean p<sub>T</sub>:

    4.1 ``YieldMean.C`` : Differential calculation of yield, mean, and their statistical and systematic uncertainties by fitting and extrapolation using fit functions.

    4.2 ``RunYieldMean.C`` : Automation to perform multiple yield and mean calculations on efficiency corrected p<sub>T</sub> spectra histograms for various fit functions/models/distributions.

    4.3 ``IntYieldsMeanPtCompare.C`` : Generate final plots by calculating extrapolation systematics from fit function choice, yield vs. multiplicity plots, mean p<sub>T</sub> vs. multiplicity plots, Omega/Xi yield ratio, linear fits, comparison to published data.

    - ``ExtraUtils.C`` : Utilities and function definitions for yield extrapolation and integrated yield calculations

### Note

- Check individual files for descriptive documentation and function definitions and usage.
- For more info on running grid analysis tasks, please see [here](https://gitlab.science.upjs.sk/jf/support/-/blob/master/doc/ALICE/README.md).
