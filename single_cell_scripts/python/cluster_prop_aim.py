################################################################################
## Title        : Cluster proportion - AIM+ cells
##
## Description  : Python script sourced by "Cluster Abundance 
##                Analysis - AIM+ cells" to run scanpro
##
## Author       : Giacomo S Frattari
################################################################################

# Import modules
import os
import pandas as pd
import numpy as np
from scanpro import scanpro

# Read data
# Get data directory
processed_data = os.getenv("PROCESSED_DATA")

cluster_data = pd.read_csv(processed_data + "steps/dab/clustertab.csv")

## Functions ###################################################################

### Bootstrap generator --------------------------------------------------------

def boot_generator(df, rep_col, rng):
  
    '''Takes a data frame with cluster and sample info, finds the sample with 
    lowest number of cells (n_min), returns a subsampled data frame with n_min
    cells from each sample'''
    
    # Find minimum number of cells
    n_min = df.groupby(rep_col).size().min()
    
    # Group by sample
    return df.groupby(rep_col, group_keys=False).apply(
        # Sample n_min number of cells without replacement
        lambda x: x.sample(
            n = n_min,
            replace=False,
            random_state=int(rng.integers(1e9))
        )
    ).reset_index(drop=True)

### Show clusters with no cells or low variance upon bootstrapping -------------

def show_zero_value_clusters(df):
      
      ''' For a given boot data frame, check cluster proportion values and
      variance. Report if any cluster has no cell or variance = 0, which 
      compromises the Bayes component of scanpro'''
      
      # Create df with number of cells per cluster
      prop_tab = df.groupby(
        ['sample','status','cluster'], as_index = False).size(
        ).rename(
        columns={'size':'n'})
      
      # Calculate proportion of cluster in total sample
      prop_tab['prop'] = prop_tab['n']/prop_tab.groupby(['sample','status'])['n'].transform('sum')
      
      # Flag clusters with no cells
      props = prop_tab[['sample','status','cluster','prop']]
      zero_cells = props[props.prop == 0]
      
      if len(zero_cells) > 0:
      
          print("This boot has %d conditions with no cells:" %(len(zero_cells)))
      
          print(zero_cells)
          
      # Calculate variance for each cluster across samples
      var_by_cluster = props.groupby(['cluster', 'status'])['prop'].var()
  
      # Flag clusters with low variance
      low_var = var_by_cluster[var_by_cluster == 0]
  
      if len(low_var) > 0:
          print(f"This boot has {len(low_var)} cluster(s) with variance == 0:")
          print(low_var)

## Analysis ####################################################################

### Pre-ATI --------------------------------------------------------------------

print("Analyzing pre-ATI")

# Subset cluster data frame
pre_ATI = cluster_data[cluster_data.clinical_time_point == "ART"]

# Scanpro on original data
prop_art = scanpro(
  pre_ATI, 
  clusters_col='cluster', 
  conds_col='status', 
  samples_col= 'sample', 
  transform='arcsin')

# Initiate random number generator with seed
rng = np.random.default_rng(0)

# Create empty dictionary for bootstraps and results
art_boot = {}
art_scanpro = {}

# Set i = 0
i = 0

while i < 100:
  
    print('ART Boot number %d' %(i))
    
    try:
        # Create boot
        art_boot[i] = boot_generator(pre_ATI, 'sample', rng)
        
        # Run scanpro
        art_scanpro[i] = scanpro(
          art_boot[i],
          clusters_col='cluster', 
          conds_col='status', 
          samples_col= 'sample',
          transform='arcsin',
          verbosity=0)
        
        # Advance in number of succeeded runs  
        i = i +1
          
    # Skip run if bootstrapping causes brentq error "f(a) and f(b) must have 
    # different signs" as for scanpro code (see source). Probably due to variance
    # (near) zero and Bayes conflicting
    except Exception as e:
        print(f"skip iter {i}: {e}")
        # Check which condition had zero variance
        show_zero_value_clusters(art_boot[i])    

### ATI ------------------------------------------------------------------------

print("Analyzing ATI")

# Select ATI
ATI = cluster_data[cluster_data.clinical_time_point == "ATI"]

# Select one visit pr donor
ATI = ATI[ATI["visit"].isin(["v24", "v11b"])]

# Scanpro on original data
prop_ati = scanpro(
  ATI, 
  clusters_col='cluster', 
  conds_col='status', 
  samples_col= 'sample', 
  transform='arcsin')

# Initiate random number generator with seed
rng = np.random.default_rng(0)

# Create empty dictionary for bootstraps and results
ati_boot = {}
ati_scanpro = {}

# Set i = 0
i = 0

while i < 100:
  
    print('ATI Boot number %d' %(i))
    
    try:
        # Create boot
        ati_boot[i] = boot_generator(ATI, 'sample', rng)
        
        # Run scanpro
        ati_scanpro[i] = scanpro(
          ati_boot[i],
          clusters_col='cluster', 
          conds_col='status', 
          samples_col= 'sample',
          transform='arcsin',
          verbosity=0)
        
        # Advance in number of succeeded runs   
        i = i + 1
          
    # Skip run if bootstrapping causes brentq error "f(a) and f(b) must have 
    # different signs" as for scanpro code (see source). Probably due to variance
    # (near) zero and Bayes conflicting
    except Exception as e:
        print(f"skip iter {i}: {e}")
        show_zero_value_clusters(art_boot[i])
