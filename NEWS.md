## FCMm 0.8.1 (2020-06-06)

  0. This is a major update for FCMm.
  1. Updated DESCRIPTION by adding fields `BugReports` and `Authors@R` for more meta infomation. The previous `Author` and `Maintainer` were deleted.
  2. Rd files in folder `man` were updated including `apply_FCM_m`, `apply_to_image`, `Assessment_via_cluster`, `Bloom`, `BR_Gil10`, `BR_Git11`, `C6`, `cal.metrics.names`, `cal.metrics`, `cal.metrics.vector.names`, `cal.metrics.vector`, `Chla_algorithms_name`, `FBA_Le13`, `FBA_Yang10`, `FCM.new`, `FuzzifierDetermination`, `generate_param`, `Getting_Asses_results`, `Gons08`, `level_to_variable`, `lmodel2_`, `NDCI_Mi12`, `OC4_OLCI`, `OC5_OLCI`, `OC6_OLCI`, `plot_spec`, `plot_spec_from_df`, `QAA_v5`, `RdYlBu`, `read_srf_excel`, `run_all_Chla_algorithms`, `Sampling_via_cluster`, `SCI_Shen10`, `Score_algorithm_interval`, `Score_algorithm_sort`, `Scoring_system`, `show_sensor_names`, `Spectral`, `SRF_simulate`, `TBA_Gil10`, `TBA_Git11`, `TC2`,  `TC2_clean`, `TC2_turbid`, `trapz`, and `trim_sd`.
  3. Updated README files by adding a new chunk for the assessment part.
  4. Updated NAMESPACE.
  5.  All logical varaible `T` and `F` were changed by `TRUE` and `FALSE`.
  6. There were naming mistakes in `tools.R` for metrics functions `cal.metrics` and `cal.metrics.vector` which replace the `SAPE` (symmetric) series with the `CAPE` (compensated) series.
  7. Updated `run_all_CHla_algorithms` which could run without stricts for wavelength settings as algorithms with missing wavelength will return NA values alternatively.
  8. Renamed and updated vignettes for simplication `Assessment` (new), `Builtin_centrodis`,  and `Cluster_new_data`. 


## FCMm 0.7.5 (2020-06-02)

  # As the version 0.7.4 was submitted to CRAN, I updated `DESCRIPTION` file follow the requirement by the CRAN which also includes `LICENSE` updated by using `uesthis::use_mit_license(name="Shun Bi")`.

## FCMm 0.7.4 (2020-06-02)

  Found and fixed an unicode `\u2010` in one Rd files of man which may result in the error check by rhub.

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

  1. Updated README files of which superlinks for vigettes are deleted as the nexting version will include other vigettes or some modifications.
  2. Renamed the filename of vignettes, also, of which superlinks are deleted. Several old version vignettes on DOC folder are also deleted.
  3. Fixed spelling errors in Rd files.
  4. In this version, the imported packge - `heatmaply` - was deleted as it would include warning information when `library(FCMm)`. Anyway, the functions supported by `heatmaply`, i.e., `Spectral()` and `RdYlBu()` are imported from packages `grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral")` and  ... ,"RdYlBu").
  5. Added cran-comments.md file as I decide to release this packge to CRAN.
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
