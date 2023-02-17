#!/usr/bin

echo "####################################"
echo "setup kernels in jupyterlab"
echo "####################################"

eval "$(conda shell.bash hook)"
/work/59896/NGS_aarhus/bin/R -e "IRkernel::installspec(name = 'ngs_r', displayname = 'NGS (R)')"
/work/59896/NGS_aarhus/bin/python -m ipykernel install --user --name ngs_python --display-name "NGS (Python)"

echo "####################################"
echo "add kernels settings"
echo "####################################"
cp /work/59896/Configurations/kernel_py.json ~/.local/share/jupyter/kernels/ngs_python/kernel.json 
cp /work/59896/Configurations/kernel_R.json ~/.local/share/jupyter/kernels/ngs_r/kernel.json