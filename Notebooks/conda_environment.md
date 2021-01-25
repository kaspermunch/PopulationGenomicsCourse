Create a conda environment where you can run the notebooks:

	conda create --name pg2020 -c anaconda -c conda-forge -c bioconda python=3.7 biopython ipyparallel jupyter jupyterlab jupyter_contrib_nbextensions matplotlib mpld3 nbconvert numpy pandas scipy seaborn statsmodels scikit-bio mygene ipywidgets  scikit-allel
	
To make ipywidgets work, install NodeJS and then run:
	
	jupyter labextension install @jupyter-widgets/jupyterlab-manager