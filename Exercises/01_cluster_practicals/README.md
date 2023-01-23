# Getting set up to use the cluster

The exercises in this course will be done using our cluster computer. It is called GenomeDK and you can [read about it here](https://genome.au.dk/). We have already recieved a username and a password for the cluster. More about that later. 

It important that you know that your account on the cluster is **temporary**. It will be deleted once the course is finished, along with any files you have on the cluster. So make sure to download any files you want to keep before the course is over. Also, your files are **not backed up**. So if you delete a file, it is gone.

The cluster is a very large collection of computers with a shared file system. Using a terminal, you can connect to the cluster from your own computer to run programs. Using other programs you can also create and edit files the same way you can on your own machine. The goal of this exercise is to make you familiar with the cluster and to get you set up to do the exercises in this course. 

## Setting up your own machine

Before we get to the cluster we need to get you properly set up on your own machine.

### Install Python

If you have not done so already, you should install a distribution of Python called *Anaconda*. Anaconda is not only an easy way of installing Python on Windows, Mac, and Linux, it also comes with the conda package management system (see below). To install Anaconda, head to [this](https://www.anaconda.com/download). When the download has completed, you must follow default installation.

### The Terminal

Most of the programs we will use in this course are command line applications. I.e. programs that are executed by writing their name and any arguments in a terminal rather than clicking on an icon and using a graphical user interface. There are many different programs that can serve as a terminal. If you have a Windows machine, you must use the *Anaconda Poweshell Prompt* (*not* the Anaconda Prompt and *not* the `CMD`). You installed Anaconda Powershell Prompt along with Anaconda Python. If you have a Mac, the terminal you will use is called *Terminal*. The Terminal application is pre-installed on Mac. So from now on, whenever we refer to the terminal, this means *Anaconda Poweshell Prompt* on Windows and *Terminal* on Mac. We will assume some familiarity with using a terminal and with executing commands on the command line. If you have not used a terminal before, or if you are a bit rusty, you should run through [this primer](https://lifehacker.com/5633909/who-needs-a-mouse-learn-to-use-the-command-line-for-almost-anything) before you go on.

<!-- 
If you are rusty or unfamiliar with the terminal, you should do this small mini-exericse:

1. Obtain the path to your home directory.
2. Create a new directory in your home directory named "ClusterPracticals". 
3. Write "This is a test" to a new file named "test_file_1" inside the directory "ClusterPracticals".
4. Change the name of the file from "test_1" to "test_file_2".
5. Make a copy of "test_file_2" named "test_file_3".
6. Delete "test_file_3".
7. Make a new directory named "test_directory_1". Delete "test_directory_1"
8. Transfer "test_file_2" to your local machine.
9. List the contents of "ClusterPracticals".
10. Make a soft link of "test_file_2" inside "ClusterPracticals".
11. See the contents of "test_file_2".
12. Create a new file named "test_file_4". Write whatever you feel like. Concatenated "test_file_2" and "test_file_4" into a new file called "test_file_5". -->

### Conda environments

You need to install packages and programs for use in your analyses and pipelines. Sometimes, however, the versions of packages you need for one project conflicts with the versions you need for other projects that you work on in parallel. Such conflicts seem like an unsolvable problem. Would it not be fantastic if you could create a small world, insulated from the rest of your Anaconda installation. Then that small world would only contain the packages you needed for a single project. If each project had its own isolated world, then there would be no such conflicts. Fortunately, there is a tool that lets you do just that, and its name is Conda. The small worlds that Conda creates are called "environments," and you can create as many as you like. You can then switch between them as you switch between your bioinformatics projects. Conda also downloads and installs the packages for you and it makes sure that the packages you install in each environment are compatible. It even makes sure that packages needed by packages (dependencies) are installed too.  By creating an enviromnet for each project, the libraries installed for each project do not interfere.

### Create a conda environment on your own machine

When you install Anaconda, conda makes a single base environment for you. It is called `base` and this is why it says "(base)" at your terminal prompt. For this course, you need a dedicated conda enviromnet for your exercises on both your local machine and on the cluster. Lets call both of them `popgen`.

The environmnet on your local machine does not need a lot of packages since it mainly serves to let you connect to the cluster. This command creates a `popgen` enviromnet and installs `slurm-jupyter` from my conda channel along with a few other packages that may come in handy:

```bash
conda create --name popgen -c kaspermunch -c conda-forge -c plotly slurm-jupyter jupyterlab pandas ipympl seaborn statsmodels popgen-dashboards
```

<!--
```bash
conda create --name popgen -c kaspermunch slurm-jupyter jupyter jupyterlab pandas numpy matplotlib ipympl nodejs seaborn scikit-learn statsmodels
```
--->

It may take `conda` a while to create the environment, so be patient. Say yes (press Enter) when asked to install packages.

**Important:** Whenever you use the terminal to do exercises, you should activate your `popgen` environment like this:

```bash
conda activate popgen
```

When the environment is active, it says `(popgen)` on the commnad prompt instead of `(base)`.

### AU VPN

To be able to connect to the cluster, you need to on the AU inernal network. You can do that by either physically being on campus, or by connecting to the AU network using VPN. To install VPN use the instructions [on this page](https://studerende.au.dk/it-support/vpn/). Before you can *use* the VPN, you need to also enable two-step verification. You can see how to do that on the same page. If you are not on phycally on campus, you need to activate your VPN before you can log in to the cluster. Your passworkd for VPN is the same as you use to log on to access Blackboard.

### Connecting to the cluster

You connect to the cluster from the terminal by executing this command (replace `username` with your cluster user name):

```bash
ssh username@login.genome.au.dk
```

When you do, you are promted for the password for your cluster username. Enter that and press enter. You are now in your home folder on the cluster. Your terminal looks the same as before but it will print:

```
  _____                                ______ _   __
 |  __ \                               |  _  \ | / /
 | |  \/ ___ _ __   ___  _ __ ___   ___| | | | |/ /
 | | __ / _ \ '_ \ / _ \| '_ ` _ \ / _ \ | | |    \
 | |_\ \  __/ | | | (_) | | | | | |  __/ |/ /| |\  \
  \____/\___|_| |_|\___/|_| |_| |_|\___|___/ \_| \_/
```

If you run the `hostname` command, you can see that you are on `fe1.genomedk.net`. Now log out of the cluster again. You do that using `exit` commannd or by presssing `Ctrl-d`. Now you are back on your own machine. Try `hostname` again and see what your own machine is called.

### Allow login without password

You will need to log in to the cluster many many times, so you should set up your `ssh` connection to the cluster so you can connect securely without typing the password every time. You do not need to know *how* this works, but if you are interested, here is roughly how:

> Firstly, you have to understand what public/private encryption keys are. A private key is a very long, random sequence of bits. A private key is kept secret and never leaves your own machine. A public key is another string of bits that is a derivative of the private key. You can generate a unique public key from the private key, but cannot get the private key from a public key: Is a one-way process. Using the public key, you can encrypt (or sign) any message, and it will only be possible to decrypt it using the private key. In other words, anyone with your public key can send you an encrypted messages that only you will be able to read. So, if the cluster has your public key saved, it can authenticate you like this: The cluster sends your machine a message that is encrypted using your public key. Your machine then decrypts the message using its private key and sends it back. If the cluster agrees it is decrypted correctly, it logs you in.

First check if you have these two authentication files on your local machine (you can do so by running `ls -a ~/.ssh` in the terminal):

```bash
~/.ssh/id_rsa
~/.ssh/id_rsa.pub
```

You most likely do not. If so, you generate a pair of authentication keys with the command below. Just press Enter when asked "Enter file in which to save the key". Do not enter a passphrase when prompted - just press enter:

```bash
ssh-keygen -t rsa
```

Now use `ssh` to create a directory `~/.ssh` on the cluster (assuming your username on the cluster is `username`). You will be prompted for your password.

```bash
ssh username@login.genome.au.dk mkdir -p .ssh
```

Finally, append the public `ssh` key on your local machine to the file `.ssh/authorized_keys` on the cluster and enter your password (replace `username` with your cluster user name):

```bash
cat ~/.ssh/id_rsa.pub | ssh username@login.genome.au.dk 'cat >> .ssh/authorized_keys'
```

From now on you can log into the cluster from your local machine without being prompted for a password.

## Setting up your home on the cluster

Now log in to the cluster

```bash
ssh username@login.genome.au.dk
```

(see, not password).

### Install Python on your cluster account

You need to install miniconda (a minimal Anaconda version) in your cluster home dir. Log in to the cluster and run this command to download the miniconda install script:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Then use this command to download and install miniconda:

```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

Follow the default installation, and say **yes** when it asks you if it should run `conda init` for you.

**NP:** Now log out out of the cluster and log back in. This is needed to make the `conda` command available to you.


### Create an environment

You also need a dedicated conda environment on the cluster. We will name this `popgen` too, but in this one we will also install all the scientific software you will use in the exercises. Run this command to create the conda envionment:

<!-- TODO: make an environment that includes bcftools that is not broken. 
Students need bcftools in the psmc exercise. For now, bcftools is installed in the software folder.
So make this work: conda create -n popgen -c bioconda -c kaspermunch bwa platypus-variant samtools beagle plink admixture gnuplot vcftools bcftools
Maybe also change exercises to only using bcftools and not samtools -->

<!--TODO: add ghostscript to the jupyter environment for use with psmc plotting: conda install -c conda-forge ghostscript -->

```bash
conda create -n popgen -c bioconda -c conda-forge -c kaspermunch bwa platypus-variant samtools beagle=4.1  plink admixture gnuplot vcftools --override-channels
```

**Important:** Whenever you log into the cluster to work on your project, you should activate your `popgen` environment like this:

```bash
conda activate popgen
```

When you environment is active it says `(popgen)` on the commnad prompt instead of `(base)`.

### Set up Jupyter

[Jupyter](https://jupyter.org/) is a notebook environment where you can easily combine text, code and plots. Using the [slurm-jupyter](https://slurm-jupyter.readthedocs.io/en/latest) tool, you can run a jupyter notebook on the cluster but see it in the browser on your own machine. That way your analysis runs on the cluster file system where your data is but the notebook interface is sent to your browser window. 

The first thing you need to do is create a separate conda environment that has jupyter installed. Do not worry about this extra environment. You will not be using it directly. We just need it to be able to run jupyter notebooks in class. 

```bash
conda create -n jupyter -c conda-forge -c bioconda -c kaspermunch slurm-jupyter popgen-dashboards jupyter jupyterlab ipyparallel pandas numpy matplotlib ipympl ipython ghostscript nodejs seaborn r-essentials rpy2 simplegeneric tzlocal r-vcfr bioconductor-biocinstaller bioconductor-snprelate r-biocmanager
```

It may take `conda` a while to create the environment, so be patient. Once created, you must activate that environemnt:

```bash
conda activate jupyter
```

and then run this this command:

```bash
config-slurm-jupyter.sh
```

That script will ask about a lot of information. You can just press enter for all of them **except** when prompted for what password you want to use: Then you must to type your cluster password.

**Now log out of the cluster**.

## Working on the the cluster

If you followed each step above, you should now be all set up. Whenever you work on your own machine, you must activate the `popgen` conda environment. When you log into the cluster, you must activate it too. The easiest way to remember, is to make sure it always says `(popgen)` on you command prompt in the terminal.

Now log in to the cluster and activate your `popgen` environemnt:

```bash
ssh usernmae@login.genome.au.dk
conda activate popgen
```
(replace `username` with your cluster user name, you probably get the drift by now)

### Your home on the cluster

When you log into the cluster, you are put in your "home folder". All users have a home folder. However, in this course you will not be using your home folder. We have made a special folder with the same name as your usename in the `populationgenomics/students` folder. You should keep everything related to the course in this folder. To get from your home folder to the this course folder you just `cd`:

```bash
cd populationgenomics/students/username
```

(replace `username` with your cluster user name)

### Running interactive commands on the cluster

When you log into the cluster you land on the "front-end" of the cluster. Think of it as the lobby of a giant hotel. If you execute the `hostname` command you will get `fe1.genomedk.net`. `fe1` is the name of the front-end machine. The "front-end" is a single machine shared by anyone who log in. So you cannot run resource intensive jobs there but quick commands are ok. Commands that finish in less than ten seconds are ok. In the exercises for this course, you will run software that takes much longer time to finish. So you need one of the computing machines on the cluster, so you can work on that instead. You ask for a computing machine by running this command:

```bash
srun --mem-per-cpu=1g --time=3:00:00 --account=populationgenomics --pty bash
```

That says that you need at most one gigabyte of memory, that you need it for at most three hours (the duration of the exercise), and that the computing expensenses should be billed to the project `populationgenomics` (which is our course). When you execute the command, your terminal will say "srun: job 40924828 queued and waiting for resources". That means that you are waiting for a machine. Once it prints "srun: job 40924828 has been allocated resources", you have been logged into a computing node. If you execute the `hostname` command, you will get something like `s05n20.genomedk.net`. `s05n20` is a computing mechine. Now you can execute any command you like without causing trouble for anyone. 

Now try to log out of the compute node by executing the `exit` command or by pressing `Ctrl-d`. If you execute the `hostname` command again you will get `fe1.genomedk.net` showing that you are back at the front-end mechine.

<!-- 

### Running commands in the terminal

When you log into the cluster you land on the "front-end" of the cluster. The "front-end" is a single machine shared by anyone who log in. You cannot run resource intensive jobs there, but quick commands are ok. Commands that finish in less than ten secons are ok. Try this command that prints "echos" the string "I can run interactive commands!" to the file `nice.txt`:

```bash
echo "I can run interactive commands!" > nice.txt
```

Use the `cat` command to show the contents of `nice.txt` in the terminal:

```bash
cat nice.txt
```

### Running interactive commands on the cluster

Say the command above was a long-running command like some population genomic analysis. Then you cannot run it on the front-end. You need to ask for one of the computing machines on the cluster so you can work on that instead. You do that by running this command:

```bash
srun --mem-per-cpu=1g --time=3:00:00 --account=populationgenomics --pty bash
```

Thay way you will use at most one gigabyte of memory, that you need at most three hours (the duration of the exercise), and that the computing expensenses should be billed to the project populationgenomics (which is our course). When you execute the command your terminal will say "srun: job 40924828 queued and waiting for resources". That means that you have asked for a machine. Once it prints "srun: job 40924828 has been allocated resources", you have been logged into a computing node. If you execute the `hostname` command you will get something like `s05n20.genomedk.net`. `s05n20` is a computing mechine. Now you can execute any command you like without causing trouble for anyone. Now try to log out of the compute node by executing the `exit` command or by pressing `Ctrl-d`. If you execute the `hostname` command again you will get `fe1.genomedk.net`. `fe1` is the front-end. -->

<!-- 
### Queueing commands on the cluster

Say the command above was a long-running command like some population genomic analysis. Then you cannot run it on the front-end. You need to submit it as a job to the cluster. When you do that, the job gets queued along with many other jobs, and as soon as the requested resources are available on the cluster, the job will start on one the many many machines. To submit a job, you must first create a file (a "batch script") that contains both the requested computer resources and the command you want to run. 

Create a file called `myscript.sh` with exactly this content:

```bash
#!/bin/bash
#SBATCH --mem=1gb
#SBATCH --time=01:00:00
#SBATCH --account=populationgenomics
#SBATCH --job-name=firstjob

echo "I can submit cluster jobs now!" > evennicer.txt
```

The first line says this is a bash script, the lines following three lines says that your job needs at most one gigabyte of memory, will run for at most one hour, that the expensenses should be billed to the project populationgenomics (which is our course). The fourth line gives the name of the name of the job. Here we have called it `firstjob`, but you should name it something sensible. 

You submit the job using the `sbatch` command: 

```bash
sbatch myscript.sh
```

Now your job is queued. Use the `mj` command to see what jobs you have queed or running. That will show something like this:

```txt
                                                                        Alloc
Job ID           Username Queue    Jobname    SessID NDS  S Elap Time   nodes
---------------- -------- -------- ---------- ------ ---  - ----------  -----
34745986         kmt      normal   firstjob       --   1  R 0-00:19:27  s03n56
```

If you want to cancel your this job before it finishes, you can use the `scancel` command:

```bash
scancel 34745986
```

Once your job finishes, it has created the file `evennicer.txt` and written "I can submit cluster jobs now!" to it. So see that you can use the `cat` command:

```bash
cat evennicer.txt
```

When you a program or script on the commandline it usually also prints some information in the terminal. When you run a job on the cluster there is no terminal to print to. Instead this is written to two files that you can read when the job finishes. In this case the fiels are called `firstjob.stdout` and `firstjob.stderr`. So see what is in them, you can use the `cat` command:

```bash
cat firstjob.stdout
```

and

```bash
cat firstjob.stderr
```

That is basically it.  -->

### How to copy files to and from the cluster

You may need to transfer files back and forth between your own machine and the cluster. To copy a file called `file` in a directory called `dir` on the cluster to the current folder on your own machine, you can use the `scp` command:

```bash
scp username@login.genome.au.dk:dir/file .
```

To copy a file called `file` in the current folder on your own machine to a folder called `dir` on the cluster, you do this:

```bash
scp ./file username@login.genome.au.dk:dir/
```

### How to run a Jupyter notebook on the cluster

Jupyter runs best in the [Chrome browser](https://www.google.com/chrome). For the best experience, install that before you go on. It does not need to be your default browser. `slurm-jupyter` will use it anyway. Now make sure you are on your own machine and that your `popgen` environment is activated. Then run this command to start a jupyter notebook on the cluster and send the display to your browser:

```bash
slurm-jupyter -u usernanme -A populationgenomics -e popgen-jupyter -m 1g -t 3h --run notebook
```

(replace `username` with your cluster user name)

Watch the terminal to see what is going on. After a while a jupyter notebook should show up in your browser window. The first time you do this, your browser may refuse to show jupyter because the connection is unsafe. In Safari you proceed to allow this. In Chrome, you can simply type the characters "thisisunsafe" while in the Chrome window:

<img src="img/thisisunsafe.png" alt="image" width="450"/>

Once ready, jupyter may ask for your cluster password. To close the jupyter notebook, press `Ctrl-c` in the terminal. Closing the browser window does **not** close down the jupyter on the cluster. You can [read this tutorial](https://www.dataquest.io/blog/jupyter-notebook-tutorial/) to learn how to use a jupyter notebook.


