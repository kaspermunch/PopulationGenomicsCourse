FROM jupyter/minimal-notebook:latest

LABEL software="PopGenomicsCourses" \
      author="Samuele Soraggi, Kasper Munch" \
      version="v2023.02.01" \
      license="MIT" \
      description="Introduction to Population Genetics"


USER 0

#expose port for dash interface
EXPOSE 8850
EXPOSE 8888

RUN mkdir -p /usr/Material/ && mkdir -p ${CONDA_DIR}/envs

COPY ./Exercises /usr/Material/Exercises
COPY ./uCloud /usr/uCloud
COPY ./Scripts /usr/Material/Scripts
COPY ./Splashscreen /usr/Splashscreen
COPY start-jupyter /usr/

RUN fix-permissions "${CONDA_DIR}"
#download data
RUN fix-permissions "${CONDA_DIR}" && \
    #if [ ! -d "/usr/Material/Data" ]; then wget https://zenodo.org/record/7551294/files/Data.tar.gz?download=1 -O /usr/Material/Data.tar.gz; \
    #                                  tar -zxvf /usr/Material/Data.tar.gz -C /usr/Material/; \
    #                                  rm -f /usr/Material/Data.tar.gz; fi && \
    #wget https://zenodo.org/record/7551294/files/popgen_au_env.tar.gz?download=1 -O ${CONDA_DIR}/envs/popgen.tar.gz && \
    #tar -zxvf ${CONDA_DIR}/envs/popgen.tar.gz -C ${CONDA_DIR}/envs/ && \
    #rm -f ${CONDA_DIR}/envs/popgen.tar.gz && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}" && \
    fix-permissions "/usr/Material" && \
    fix-permissions "/usr/uCloud"

#RUN ln -s /usr/Material ./Course_Material
#RUN eval "$(mamba shell.bash hook)"

## Add JupyterLab Extensions
#RUN printf "Install JupyterLab extensions:" \
 #&& pip install --no-cache-dir "nteract-on-jupyter" \
 #&& jupyter labextension install "jupyter-threejs" \
 ## && jupyter labextension install "ipyvolume" \
 ## && jupyter lab clean -y \
 ## add support for LaTeX docs
 ## && pip install --no-cache-dir "jupyterlab-latex" \
 ## open spreadsheets such as Excel and OpenOffice
 ## && jupyter labextension install "jupyterlab-spreadsheet" \
 ## && jupyter lab clean -y \
 ## add top bar
 #&& pip install --no-cache-dir "jupyterlab-topbar" \
 #&& jupyter labextension install "jupyterlab-topbar-text" \
 #&& jupyter lab clean -y \
 ## add system monitor
 #&& pip install --no-cache-dir "jupyterlab-system-monitor" \
 ## add theme toggle bottom
 #&& jupyter labextension install "jupyterlab-theme-toggle" \
 #&& jupyter lab clean -y \
 ## add code formatter
 #&& pip install --no-cache-dir "autopep8" "yapf" "isort" "black" \
 #&& pip install --no-cache-dir "jupyterlab_code_formatter" \
 #&& jupyter lab build -y \
 #&& jupyter lab clean -y \
 #&& fix-permissions "/home/${NB_USER}" \
 ## add variableInspector
 #&& pip install --no-cache-dir "lckr-jupyterlab-variableinspector" \
 ## add nbdime
 #&& pip install --no-cache-dir "nbdime" \
 ## add Bokeh extension
 #&& pip install --no-cache-dir "jupyter_bokeh" \
 ## add Plotly extension
 #&& pip install --no-cache-dir  "plotly" \
 #&& jupyter labextension install "jupyterlab-plotly" \
 #&& fix-permissions "/home/${NB_USER}" \
## && fix-permissions "${CONDA_DIR}"

#permissions
#RUN fix-permissions "${CONDA_DIR}" && \
#    fix-permissions "/home/${NB_USER}" && \
#    fix-permissions "/usr/Material"

#create environments
#RUN mamba env create -p "${CONDA_DIR}/envs/popgen_aarhus" -f /usr/Material/Environments/environment_ucloud.yml \
#    && mamba clean --all -f -y

#install kernels
#RUN ${CONDA_DIR}/envs/popgen_au_env/bin/python -m ipykernel install --user --name="args_dashboard" --display-name "args dashboard" && \
#    ${CONDA_DIR}/envs/popgen_au_env/bin/R -e "IRkernel::installspec(user=TRUE, name = 'popgen_course', displayname = 'popgen course')" && \
#    ${CONDA_DIR}/envs/popgen_au_env/bin/python -m bash_kernel.install && \
#    fix-permissions "/home/${NB_USER}" && \
#    ln -s /usr/Material /home/${NB_USER}/work/Course_Material
    



### modify kernel files with system variables and fix library for bcftools
#RUN cp /usr/uCloud/kernelBash.json /usr/local/share/jupyter/kernels/bash/kernel.json && \
#    ln -s ${CONDA_DIR}/envs/popgen_au_env/lib/libcrypto.so.3 ${CONDA_DIR}/envs/popgen_au_env/lib/libcrypto.so.1.0.0


### install splashscreen

#RUN jupyter labextension install /usr/Splashscreen/


USER 1000

WORKDIR /home/${NB_USER}/work/Course_Material

COPY --chown="${NB_USER}":"${NB_GID}" start-jupyter "${CONDA_DIR}"/bin/
RUN pip install --no-cache-dir conda-pack && \
    chmod +x "${CONDA_DIR}"/bin/start-jupyter