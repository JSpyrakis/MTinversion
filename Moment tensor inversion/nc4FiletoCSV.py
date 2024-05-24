import numpy as np
import pandas as pd
from statsmodels.regression.linear_model import OLS
from statsmodels.formula.api import mixedlm
from statsmodels.tools.tools import add_constant
import matplotlib.pyplot as plt
import netCDF4 as nc

# Setting path and filename
ncfname = "C:/Users/johns/Desktop/VSCode/Python/Moment tensor inversion/SOCAL_5PC_PERTURBED_ELEMGFS_updated4.nc4"
dname = "SoCal"
dname2 = "obsdata"
dname3 = "Gtensors"



# Opening the netCDF file
with nc.Dataset(ncfname, 'r') as ncin:
    elem_array = ncin.variables[dname][:]
    obs_array = ncin.variables[dname2][:]

# Convert numpy array to pandas DataFrame
def array_to_dataframe(array):
    data = array.reshape(-1, array.shape[-1])
    idx = pd.MultiIndex.from_product([range(s) for s in array.shape[:-1]], names=['time', 'element', 'channel', 'station'])
    return pd.DataFrame(data, index=idx)

elem_df = array_to_dataframe(elem_array)
obs_df = array_to_dataframe(obs_array)

# Merge slice with observed data
merge_df = pd.concat([obs_df] + [elem_df.xs(m, level='element', axis=1) for m in range(1, 7)], axis=1)

# Writing merged data to a CSV file
merge_df.to_csv('merge2.csv')

# Crude least squares fit
X = merge_df[['g1', 'g2', 'g3', 'g4', 'g5', 'g6']]
y = merge_df['value']
mm = OLS(y, X).fit()

# Shift covariates
# You can use mm.params for beta and mm.scale for sigma^2 in further calculations

# Identify shifted green's functions and fit lm model with no weights
# You can use mm.predict() for predictions

# Fit lme model with no weights
merge_df['groups'] = np.repeat(range(1, 19), 11)
mm3 = mixedlm("value ~ -1 + g1 + g2 + g3 + g4 + g5 + g6", data=merge_df, groups=merge_df['groups']).fit()

# Calculate standard errors
se = mm3.bse
print(se)

# Plotting
# You can use mm3.fittedvalues for predictions

