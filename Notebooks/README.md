
# Notebooks for students

This is a collection of notebooks that you can play with to aid your understanding of population genetic processes. This course feature is experimental. Click on the binder badge below to start a server to view and play with the notebooks.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/kaspermunch/PopulationGenomicsCourse/HEAD?filepath=Notebooks)

If you want to save your notebook, do not save it on the server. Instead, you whould select `Files > Download as` and then select `Notebook .ipynb`.

## How to run Jupyter notebooks on your own machine

To work on notebooks downloaed to your our own machine, use the terminal to navigate to the folder where you downloaded the notebooks. Then activate your environment and start the juptyer server:

```bash
conda activate popgen
jupyter notebook
```

## Set up repository for Binder (for instructors only)

Created a special environment just for running the notebooks, which is then exported to `environment.yml' in the 'binder' dir under root:

    conda create -n popgen_binder -c conda-forge jupyter jupyterlab pandas numpy matplotlib ipympl nodejs seaborn scikit-learn statsmodels

    conda env export --from-history -f ../binder/environment.yml

**NB:** If any new notebook dependencies are added they need to be added to the exported conda environment too:

```bash
conda activate popgen_binder
conda install some_needed_library
conda env export --from-history -f ../binder/environment.yml
```