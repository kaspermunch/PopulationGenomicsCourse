# Linux server for exercises
During the course there will be several exercises that use command line tools that are only available on linux machines. These programs have been installed on a linux server that you can log onto using ssh.

We will refer to [user], which is your user name. We have created an account for each of you with your first name (all lowercase) as user name and your surname (all lowercase) as password (you can use the command passwd to change the password).

## Login from linux or mac
 If you are on a linux or mac computer then ssh is usually installed by default and you can use the following command to log in to the server:

```bash
ssh -p 8922 [user]@185.45.23.197
```
To copy a file called [file] in a directory called [dir] on the server to your local machine you can use the scp command (you can also use sftp if you prefer that):

```bash
scp -P 8922 [user]@185.45.23.197:[dir]/[file] .
```
In order to copy files from your local machine to the cluster, you just need to inverse the order:

```bash
scp -P 8922 ./[file] [user]@185.45.23.197:[dir]/
```

If your file is not in the current directory, you just need to provide the entire path.

## Login from windows

There are two ways to connect from a Windows machine:

### Using Moba Xterm (which has a Graphical User Interface)

It can be downloaded here: https://mobaxterm.mobatek.net/. Click "GET MOBAXTERM NOW" and then click to get the Free version. Then choose "MobaXterm Home Edition v11.1 (Portable edition)". It is important that it is the "portable" version. You can download it to your Desktop. Double-click it to open it.

To connect to the machine, click "Session" in the top left corner. Then click SSH in the top left corner of the window that pops up. Now fill in the address (185.45.23.197), [user] and the port (8922). 

You should now see a command prompt where you can type commands. The left part should show the folder structure of the folders under your user folder. Fom here you can easily download files to your local machine.

### Using PuTTy

It can be downloaded here: https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html. To log in you can write:

```bash
plink -X -P 8922 [user]@185.45.23.197
```

And you can transfer files between the cluster and your local machine by running:

```bash
pscp -P 8922 [user]@185.45.23.197:[dir]/[file] .
```

In order to copy files from your local machine to the cluster, you just need to inverse the order:

```bash
pscp -P 8922 ./[file] [user]@185.45.23.197:[dir]/
```

If your file is not in the current directory, you just need to provide the entire path.

If you have any problems with the upper commands, there is a more interactive way for both accessing the cluster and transferring files in a windows machine.

To transfer files, you can open PuTTy and set the options like this:

![Alt text](https://user-images.githubusercontent.com/38723379/51401061-30b1cf80-1b4a-11e9-8600-c3f11228ff91.png)

To copy files from the cluster to your local machine and vice-versa in an interactive way, you will need to install WinSCP. It can be downloaded here:
https://winscp.net/eng/download.php

Then, in the program set the options like this:

![Alt text](https://user-images.githubusercontent.com/38723379/51401251-94d49380-1b4a-11e9-8f07-7c58bc7238fb.png)


## Using the command line
If you are not used to using a unix command-line interface you might find this introduction useful:
https://lifehacker.com/5633909/who-needs-a-mouse-learn-to-use-the-command-line-for-almost-anything

After going through the introduction (if necessary), try do the following mini exercises:

1) Obtain the path to your home directory.
2) Create a new directory in your home directory named "ClusterPracticals". 
3) Write "This is a test" to a new file named "test_file_1" inside the directory "ClusterPracticals".
4) Change the name of the file from "test_1" to "test_file_2".
5) Make a copy of "test_file_2" named "test_file_3".
6) Delete "test_file_3".
7) Make a new directory named "test_directory_1". Delete "test_directory_1"
8) Transfer "test_file_2" to your local machine.
9) List the contents of "ClusterPracticals".
10) Make a soft link of "test_file_2" inside "ClusterPracticals".
11) See the contents of "test_file_2".
12) Create a new file named "test_file_4". Write whatever you feel like. Concatenated "test_file_2" and "test_file_4" into a new file called "test_file_5".

If you are done and you want a bit of a challange, you can do the following:
(grep, awk, vi, bash_profile, permissions)
13) 


