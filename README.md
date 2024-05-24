# MTinversion
python tools for seismic moment tensor inversion

This currently only compute for the data set of earthquakes in the Long Valley Caldera California at this current momemnt.

Moment tensor inversion file contains main.py which contains the moment tensor inversion code. At the moment there is an error with the final amplitude calculation which is undergoing fixing. For a semi-working look at the method check out the R script EM2_lme_iterate_gld.R. Note that main.py needs merge1.csv to work and EM2_lme_iterate_gld.R needs SOCAL_5PC_PERTURBED_ELEMGFS_updated4.nc4 to work.
