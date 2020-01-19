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
