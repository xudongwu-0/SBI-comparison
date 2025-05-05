import rpy2.robjects as robjects
import numpy as np

from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector

robjects.r('rm(list = ls())')

# Load necessary R packages
Matrix = importr('Matrix') 

# Check the current working directory in R
print("Current working directory in R:")
print(robjects.r('getwd()'))

# Load the script "Functions.R" and the object with the Rdata
robjects.r('source("Functions.R")')
robjects.r('load("SpainPROV_CerebrovascularDiseases.Rdata")')

print("Objects in the R environment:")
print(robjects.r('ls()'))

# Print the first rows of the Data object in R
print(robjects.r('head(Data)'))

# Call the function 'samples.iCAR'
robjects.r('Obs <- samples.iCAR(Data, tau.range=c(4,400), intercept=0, k=10, s=10, l=10)')
robjects.r('str(Obs,1)')

# Access the 'Obs' object in R and create a dictionary to store the results
obs = robjects.r('Obs')
obs_python = {}

for i in range(len(obs)):
    matrices = obs[i]  # Get the list of matrices
    
    # Convert each matrix from R into a numpy object
    tau_matrices = [np.array(mat) for mat in matrices]
    
    # Store the matrices in the Python dictionary
    obs_python[i] = tau_matrices

# Print the first two matrices in obs_python[1]
print(obs_python[1][0])
print(obs_python[1][1])