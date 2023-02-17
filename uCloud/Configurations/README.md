# Population genomics

Repository for the ucloud version of the NGS summer course at Aarhus university.


# How to configure this course on Ucloud

Here the steps to configure this course on uCloud. 


## Folder structure

First of all create a project where you are administrator and there is a PIs group. All students must be in a separate students' group. Create a drive in your project and put in it a `Data` folder (with read-only rights to anyone but the administrator) that contains the necessary data for the course. 

In the same drive copy the folder `Configurations` from this repository. Then make a folder called `NGS_summer_2022`, and copy in it the folder `Notebooks` and `Scripts` from this repository. Go into the folder with `cd NGS_summer_2022` and create a symbolic link to the `Data` folder with its absolute path
```
ln -s Data /work/${DRIVE_NUMBER}/Data
```
where ${DRIVE_NUMBER} is the identifying number of your drive, that can be found in its properties (drives are not identified by name). The whole content of the `NGS_summer_2022` folder and of the `Configurations` folder must be read-only for anyone but the administrator. 


## Create the `conda` environment

Copy the `environment.yml` file in the `Environments` folder of this repository somewhere, and with conda or mamba install the environment with
```
mamba create env -f ${PATH_TO_ENVIRONMENT}/environment.yml -p /work/${DRIVE_NUMBER}/NGS_aarhus
```


## Setup the configuration file for jupyterlab

In the folder `Configurations` there is a `json` file that configures the jupyterlab sessions. You need to simply modify the two paths to the environments setup
```
"path": "/59896/Configurations/setup.sh"
```
and the drive identifier itself
```
"path": "/59896"
```
so that they match the drive in use. You can modify other options, such as the amount of resources chosen by default and the job expiration time.

## Change the setup file for the conda environment 

Open the file `setup.sh` in the folder `Configurations`. This is a script to copy the course material and link the conda environment with jupyterlab. In this file you need only to correct the drive identifier.


## Start a job with uCloud

### The first time
To use this material for the first time, you simply need to run a `jupyterlab session`. Open `jupyterlab` from the Apps menu, and then use the `browse` button to find the `json` file in the `Configurations` folder. Then click on `Add folder` and on the empty field appearing, so that you can browse and choose your personal folder (It has the name in this format: `NameSurname#Random_number`). Click on submit, and when your session starts, click on `open interface` (top-right corner). When jupyterlab opens, drag and drop from the explorer window of jupyterlab the folder `NGS_summer_2022` into the user personal folder.

### Any other time
Same procedure as above, with the difference that when you click `Add folder`, you choose the course folder you copied into your personal folder. The setup script will notice you have chosen the course folder, and will not copy it again, so you can work on your own material.

