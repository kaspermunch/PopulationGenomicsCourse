


## Binder

    conda create -n popgen_binder -c conda-forge jupyter jupyterlab pandas numpy matplotlib ipympl nodejs seaborn  ipywidgets scikit-learn statsmodels

    conda env export --from-history -f environment.yml


[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/kaspermunch/PopulationGenomicsCourse/HEAD?filepath=Notebooks)


# How to run Jupyter notebooks

To run the notebooks you first need to install [Anaconda Python](https://www.anaconda.com/download/).

## Create a conda environment

The notebooks use various libraries that you need to install too to be able to run them. The best way to set this up is to create a conda environment with the proper libraries installed. To do that, copy/pasete this long command into your command prompt/Terminal and run it.

    conda create -c anaconda -c conda-forge -c bioconda --name popgen python=3 biopython jupyter matplotlib mpld3 numpy pandas scikit-learn scipy seaborn statsmodels pytables pyfaidx
    
Press enter once asked if you want to install a long list of libraries. 

## Running Jupyter

Once it finished there are two ways to run Jupyter:

- On Windows the easiest way is to find and open the program called Anaconca Navigator. Then under "Home" you can choose the environment "popgen" and then launch Jupyter.
- On linux and mac you can just type `conda activate popgen`. (to Deactivate type `conda deactivate`)

Jupyter opens a file browser. Navigate to the folder where you downloaded the notebook files from Blackboard and click on one of them. Then the notebook open and you can see the code.


