## FCMm 0.9.4 (2020-11-05)

  Fix bugs in `plot_spec_group`. Parameter `group_num` could be inputted as characters such as "%s %s".

## FCMm 0.9.3 (2020-11-02)

  1. Improve the `trapz` function. Now better to use the `trapz2` function. Much faster!
  2. Add a new usefull function namely `plot_spec_group` which helps to plot spectra lines by means of `ggplot2` with facet in different colors.

## FCMm 0.9.2 (2020-11-02)

  Not sure why is error report from travis so I just comment the travis lines in my readme page.

## FCMm 0.9.1 (2020-11-01)

  1. Parameter `stand` in function `FuzzifierDetermination`, `apply_FCM_m` or so has been deprecated now as users reported that it has some **ambiguity** about the standardized operation in FCMm. As an alternative, I set another parameter called `do.stand`. If `do.stand == TRUE` then the function will standardize the input data. The `stand` parameter is still be supported as before.
  2. A new function (Schwammle2010) to obtain optimal m value is added, namely `fDN`, but it has not been exported it. Call `FCMm:::fDN` to run it.
  3. An improved `Fuzzifierdetermination` version `Fuzzifierdetermination2` has been created which run much faster than the former one.
  4. Corresponding documents have been updated.

## FCMm 0.8.7 (2020-07-10)
  1. Updates documents in man.
  2. Cancel the export of function `read_srf_excel`.
  3. One of outputs by `Scoring_system` is changed from `opt_algorithm` to `Opt_algorithm` which is used to match with `Scoring_system_bootstrap`.
  4. Cancel the export of functions `TC2_clean` and `TC2_turbid` since their outputs could be found in `TC2_Liu20`.
  5. Add a new parameter `pos_rgb` to `apply_to_image` for definition of RGB channels.
  6. Update README.

## FCMm 0.8.6 (2020-07-01)
  1. `Assessment_via_cluster`: update the precision calculation process which now produce NA values for extreme difference between actual and predicted values.
  2. `Assessment_via_cluster`: for total precision calculation, it is unfair when outliers are included (which may result in extremely precision values). A quantile calculation for limiting the range of values to be calculated.
  3. `Assessment_via_cluster`: the function will return NA values if all prediction from j algorithm in cluster i are NA values, i.e., the algorithm to be assessed is basically failed.
  4. `Getting_Asses_result`:  add a new parameter `replace` to control the function `sample` for using replacement or non-replacement.
  5. `Scoring_system_bootstrap`: add a parameter `replace` to control the function `sample`.
  6. `Scoring_system_bootstrap`:  add an output `Results_of_scoring_system` to restore the all outputs from `Scoring_system`.
  7. `Scoring_system_bootstrap`: fix bug of determining `Opt_algorithm`. Now it is determined by maximum of score values.
  8. `Scoring_system_bootstrap`: add new item `Remove_algorithm` which is used for removing algorithms with huge errors when blending.
  9. `Scoring_system_bootstrap`: update the col plotting.
  10. `Scoring_system_bootstrap`: add Chla blending process and its corresponding plots.
  11. Add a new internal function `Chla_algorithm_blend` which is used for blending work.
  12. `apply_FCM_m`: fix the bug of using `color_palette`.
  13. `apply_to_image`: the default setting and output control of parameters is now more reasonable. 
  14. Add a new function `generate_param_ex`, which is used for generating input for `apply_to_image`when user-defined centroids come in.

## FCMm 0.8.5 (2020-06-21)

  1. Fix the bug in `Assessment_via_cluster` when metric values is extremely biased.
  2. Fix the bug in `Sampling_via_cluster` when the input `x` has noncontinuous numeric names.
  3. Fix the bug in `Scoring_system` when the `Total_score` has full NA values (for some cluster have few observation to calculate metrics.
  4. Now, the parameter `Times` of `Scoring_system_bootstrap` must be larger than 1 to avoid errors.
  5. New features! The `Scoring_system_bootstrap` now support algorithms blending by optimized candidates. 
  6. New features! `FCM.new` now support assign cluster names by sorting based on specific constrains (new parameters `sort.pos` and `sort.decreasing`).
  7. New features! Added parameter `color_palette` for spectra plotting functions, i.e., `plot_spec` and `apply_to_image`. Users could use built-in color palettes or pre-defined color codes.
  8. Documents updated.

## FCMm 0.8.4 (2020-06-13)

  1. Import two required packages into `FCMm`, i.e., `scales` and `farver`.
  2. Update documents in the man folder: `apply_to_image`, `FCM.new`, `Getting_Asses_results`, `Scoring_system`.
  3. I have changed the import way of the package `raster` since it has a conflict with `base::aggregate`. Thus, the `raster` functions are used by the way `importFrom`.
  4. Improve the function `Assessment_via_cluster`. The valid percent observation is defined by values greater than zero and `is.finite` function.
  5. Fix the bug of wrongly using `SMAPE` (actually is `CMAPE`) and so on.
  6. Update the document of function `Scoring_system`. More and more details now! Non-related comments were removed.
  7. Although the function `Scoring_system` could be used in the bootstrap way by users, in this version, I add a function `Scoring_system_bootstrap` to do that. The default sample time is set as 1000.
  8. Fix the bug in function `FCM.new`. Previously, `x.stand` is `x` when `stand = TRUE` which is unfriendly for further usage. Now, both raw and normalized data are calculated in the correct way! The centroids are saved on both raw and normalized scales.
  9. Add the function `HUE` which is the default color palette of the package `ggplot2`.
  10. Delete `rm(list=ls())` in `Builtin_centroids.Rmd`.
  11. Update documents in `Cluster_new_data.Rmd`.
  12. Update README files.


## FCMm 0.8.3 (2020-06-12)

  1. Thanks for the bug report by Xiaolan Cai, a new return `centroids` of function `FCM.new` is supported now, which indicates the cluster centroids of your training data set.
  2. Add a new `ifnotstop` for function `Assessment_via_cluster` that ZERO values in measured vector will introduce `Inf` metric results for `cal.metrics` or `cal.metrics.vector`.
  3. Add a new color palette function `scales::hue_pal` which is the default palette of package `ggplot2`.

## FCMm 0.8.2 (2020-06-11)

  0. This update is for revision from CRAN manual inspection.
  1. All examples in Rd files were unwrapped which used \dontrun{} except for apply_to_image.Rd as its executable time is > 5 sec.
  2. Functions `cal.metrics.names`, `cal.metircs.vector.names`, `Chla_algorithms_name`, and `show_sensor_names` now are alias of `cal.metrics`, `cal.metrics.vector`, `run_all_Chla_algorithms`, and `SRF_simulate` respectively.
  3. In this version, `on.exit` is used in the right position which is next to `oldoptions <- options(scipen=1000)` #247 of `Image_application.R`.
  4. The LICENSE and cph are fixed in this version.
  5. Updated the function `plot_spec`. For convenience of plotting spectra by groups. Now the `facet_wrap` are used rather than restored in a list. Therefore, its argument `HABc` is deprecated in this version.
  6. Documents are updated for grammar and spelling issues.
  7. TSM estimation algorithm `GAA_SPM` by Xiaolong Yu is supported in this version.

## FCMm 0.8.1 (2020-06-06)

  0. This is a major update for FCMm.
  1. Updated DESCRIPTION by adding fields `BugReports` and `Authors@R` for more meta information. The previous `Author` and `Maintainer` were deleted.
  2. Rd files in folder `man` were updated including `apply_FCM_m`, `apply_to_image`, `Assessment_via_cluster`, `Bloom`, `BR_Gil10`, `BR_Git11`, `C6`, `cal.metrics.names`, `cal.metrics`, `cal.metrics.vector.names`, `cal.metrics.vector`, `Chla_algorithms_name`, `FBA_Le13`, `FBA_Yang10`, `FCM.new`, `FuzzifierDetermination`, `generate_param`, `Getting_Asses_results`, `Gons08`, `level_to_variable`, `lmodel2_`, `NDCI_Mi12`, `OC4_OLCI`, `OC5_OLCI`, `OC6_OLCI`, `plot_spec`, `plot_spec_from_df`, `QAA_v5`, `RdYlBu`, `read_srf_excel`, `run_all_Chla_algorithms`, `Sampling_via_cluster`, `SCI_Shen10`, `Score_algorithm_interval`, `Score_algorithm_sort`, `Scoring_system`, `show_sensor_names`, `Spectral`, `SRF_simulate`, `TBA_Gil10`, `TBA_Git11`, `TC2`,  `TC2_clean`, `TC2_turbid`, `trapz`, and `trim_sd`.
  3. Updated README files by adding a new chunk for the assessment part.
  4. Updated NAMESPACE.
  5.  All logical variables `T` and `F` were changed by `TRUE` and `FALSE`.
  6. There were naming mistakes in `tools.R` for metrics functions `cal.metrics` and `cal.metrics.vector` which replace the `SAPE` (symmetric) series with the `CAPE` (compensated) series.
  7. Updated `run_all_Chla_algorithms` which could run without constrains for wavelength settings as algorithms with missing wavelength will return NA values alternatively.
  8. Renamed and updated vignettes for simplification `Assessment` (new), `Builtin_centrodis`,  and `Cluster_new_data`. 
  9. Codes about `options()` were removed or used by `on.exit()`.

## FCMm 0.7.5 (2020-06-02)

  # As the version 0.7.4 was submitted to CRAN, I updated `DESCRIPTION` file follow the requirement by the CRAN which also includes `LICENSE` updated by using `uesthis::use_mit_license(name="Shun Bi")`.

## FCMm 0.7.4 (2020-06-02)

  Found and fixed an Unicode `\u2010` in one Rd files of man which may result in the error check by rhub.

## FCMm 0.7.3 (2020-06-01)

  1. Refresh the Rd files in man.
  2. Export the function `trapz()`.
  3. Updated README and vignettes.

## FCMm 0.7.2 (2020-06-01)

  1. Deleted the suggested package `dplyr` as I can use `base::sample` as an alternative.
  2. The previous four vignettes are shorten to the two that are `Cluster_a_new_dataset_by_FCM` and `Usage_of_the_built_in_centroids_by_FCMm` which might provide enough help for users.
  3. Files in DOCS are updated.
  4. README updated.

## FCMm 0.7.1 (2020-05-30)

  1. Updated README files of which superlinks for vignettes are deleted as the nexting version will include other vignettes or some modifications.
  2. Renamed the filename of vignettes, also, of which superlinks are deleted. Several old version vignettes on DOC folder are also deleted.
  3. Fixed spelling errors in Rd files.
  4. In this version, the imported package - `heatmaply` - was deleted as it would include warning information when `library(FCMm)`. Anyway, the functions supported by `heatmaply`, i.e., `Spectral()` and `RdYlBu()` are imported from packages `grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral")` and  ... ,"RdYlBu").
  5. Added cran-comments.md file as I decide to release this package to CRAN.
  6. Updated man Rd files.

## FCMm 0.6.2 (2020-05-24)

  * Bugs fixed.

## FCMm 0.6.1 (2020-05-23)

  * This version updated the arrangement of codes, fixed bugs, and collected data sets.
  * Ignored files are removed from Gihub.
  * I gonna resolve the problem of `do` parameters in next versions.

## FCMm 0.5.7 (2020-05-12)

  * Bugs of `SRF_simulate` fixed.

## FCMm 0.5.6 (2020-03-26)

  * Bugs of `apply_to_image` fixed.

## FCMm 0.5.5 (2020-03-17)

  * Hey guys, due to the outbreak of coronavirus, I have stopped the update for a long time. 
  * In this version. The assessment method in `FCMm` has been improved by a new proposed scoring rule called sort-based method.  
  * I am collecting the data and references about the algorithm assessment and finally going to finish my manuscript. Hope this will help to improve the algorithm evaluation work the in water color community.  
  * Please, COVID-19 go away!

## FCMm 0.5.4 (2020-03-03)

  * Update the document of data `OLCI_TH`.  
  * Improve the effectiveness of function `apply_FCM_m` as the `matrix` may speed up the process rather than the `data.frame` (I suppose).

## FCMm 0.5.3 (2020-02-28)

  * Add newly trained `Bloom` model which was used to replace `C6` model.

## FCMm 0.5.2 (2020-02-22)

  * Bugs of `Sampling_via_cluster` were fixed.

## FCMm 0.5.1 (2020-02-21)

  * Create the new file `QAA_series.R` to add functions `TC` and `QAA_v5`. The original functions from `FCM_m_Chla_estimation` were removed accordingly.
  * Create the new file `Algorithms_assessment.R` to add functions `Assessment_via_cluster`, `Score_algorithms`, `Sampling_via_cluster`, `Getting_Asses_result`, and `Scoring_system`. The original function `Assessment_via_cluster` was removed from `Algorithms_assessment.R`.
  * See more details in help documents of these new functions.

## FCMm 0.4.17 (2020-02-20)

  * Modify the constrain of SMA linear regression such as Rsquare value and point number taken part into the calculation.
  * Next update is a first-point caning version, 0.5.X since the new assessment system (scoring system) will be load ASAP.

## FCMm 0.4.16 (2020-02-19)

  * Add several new **Error Metrics** into `tools.R`. 
  * Use `Type-II` regression, namely SMA (Standard Major Axis), to obtain the regression results, such as `Intercept` and `Slope`. The main reason for that is, inputs of regression (for validation) are obtained without fully quality control (still have uncertainties in each variable). This kind of regression, comparing the `Type-I`, cannot be used for prediction (such as using Rrs to predict Chla concentrations), but better to compare two variables on the same scale (units), useful for the error assessment.

## FCMm 0.4.15 (2020-02-17)

  * Add linear regression results (`slope` and `intercept`) to `tool.R` for hard mode assessment.

## FCMm 0.4.14 (2020-02-16)

  * Add new plot option for `Assessment_via_cluster`

## FCMm 0.4.13 (2020-02-15)

  * Result `Area` was added in the return of `apply_FCM_m`.

## FCMm 0.4.12 (2020-02-14)

  * Bugs of `TC2` fixed.

## FCMm 0.4.11 (2020-02-13)

  * Bugs of `BR_Gil10` and `NDCI_Mi12` fixed.

## FCMm 0.4.10 (2020-02-13)

  * Bugs of `FCM_m_Chla_estimation` fixed.

## FCMm 0.4.9 (2020-02-10)

  * Add `UMRPE` function to calculate unbiased error metric.

## FCMm 0.4.8 (2020-02-09)

  * Update the function `cal.metrics` which remove the double version of `.cal.mae` or `.bias`. The log10-transformation is now supported by the logical parameter `log10`. 
  * Add functions `cal.metrics.vector` and `cal.metrics.vector.names` which calculate the validation metrics via pairwise way (results are columns). 
  * Update the function `Assessment_via_cluster` by adding a fuzzy mode to calculate validation metrics via membership values.

## FCMm 0.4.7 (2020-02-08)

  * Add `QAA_v5` in file `FCM_m_Chla_estiamtion`.
  * Some bugs fixed.

## FCMm 0.4.6 (2020-02-05)

  Update the functions and references in file `FCM_m_Chla_estiamtion`.

## Notes (2020-01-21)

  * Even though the coming **phylogenetic tree method** will offer a good option to post-classify the numerous clusters, we still have to design a better allocation strategy on **cluster number**. That's to say, the nearby cluster name (for instance, `Cluster 1` and `Cluster 2`) should have the similar optical attribution. However, the current `FCM.new` or `FuzzifierDetermiantion` functions return the totally random results. Should be improved later!
  * Another **confusing** thing is the parameter `stand` in many functions (for instance, `FCM.new`). The original design of `stand` is making sure the input `x` has been standardized so that the median-process of clustering should not standardize that. Actually, the parameter could be better understood if named as `DoStand`. When `DoStand` is `True`, it means the input `x` have to be standardized in the following process. This change will make numerous source codes different so I will get to do it after other functions are gradually improved.

## FCMm 0.4.5 (2020-01-21)
  
  * Fix bugs of function `Assessment_via_cluster` to avoid produce NA values.

## FCMm 0.4.4 (2020-01-20)

 * Add new function `run_all_Chla_algorithms` to obtain the Chla concentration if required Rrs bands are get.
 * Group the function list by adding the `@family` tag using package `roxygen2`.

## FCMm 0.4.3 (2020-01-19)

 * Supply the help documents of function `SRF_simulate`.
 * Add an option of function `SRF_simulate` --- `wv_as_column` which makes the output as a dataframe with wavelength as column name. Easy to be an input of function `apply_FCM_m` or so.
 * In next version, I gonna add the Chla estimation method proposed by Liu *et al.* (2020) which has been accepted by RSE.


## FCMm 0.4.1 (2020-01-12)

 * Recently, I reread Moore's and Jackson's papers and have some new insights.
    - I should re-consider the FCM process when input data (here is Rrs spectra) is biased.
	- Mahalanobis distance should be under consideration since it eliminate the influence of variable correlation. However, it is questionable to define the covariance when facing new input data.
	- Moore *et al.* (2001) supposes the distribution of Rrs vectors belonging to one class is multivariate normal, and Rrs belongs to that class, the distance metric has a chi-sq distribution with n degrees of free (where n is the dimension of Rrs). Is this evidential in current data sets?
 
 * Phylogenetic tree is an effective method to post-classification (or post-group) of initial cluster results.
   - This mode will be added in **version 0.5.X or later** (2020-01-19).
 
 * How to assign new data to an existing clustering?
   - Still consider this problem ... If you guys have some new ideals, please let me know.


## FCMm 0.4.0 (2020-01-10)

 * The input of `apply_FCM_m()` with built-in clusters should match OLCI band settings.
   - However, the data of users may be **hyperspectral** Rrs. Need a function that convert hyperspectral Rrs into the **multispectral**.
   - Should add a function like `SRF_simulate()` to simulate the Rrs by spectral response function.
   - Meanwhile, the `SRF_list.rda` as a list with SRF data of many sensors should be added.
   
 * Need another vignette to illustrate how we can bootstrap the training set for optimizing the cluster number or fuzzifier value.
 
 * Now the result of function `apply_FCM_m()` is plotted on single picture which is hard to read when the number of spectra is large. So I gonna add a new plot function by using `facet` method from `ggplot2`.
   - A new result of `apply_FCM_m()` is now added as `p.group.facet`. You could `yourresult$p.group.facet` to obtain the spectra with facet.


## FCMm 0.3.3 (2020-01-09)

 * Add `NEWS.md` file in the top-level fold so that you could see what changes in each version of `FCMm`.
