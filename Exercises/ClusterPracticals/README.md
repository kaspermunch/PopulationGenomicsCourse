# Linux server for exercises
During the course there will be several exercises that use command line tools that are only available on linux machines. These programs have been installed on a linux server that you can log onto using ssh.

## Login from linux or mac
 If you are on a linux or mac computer then ssh is usually installed by default and you can use the following command to log in to the server:

```bash
ssh -p 8922 [user]@185.45.23.197
```

Where [user] is your user name. We have created an account for each of you with your first name (all lowercase) as user name and your surname (all lowercase) as password (you can use the command passwd to change the password).

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

If you are on a windows machine you need to first install PuTTy. It can be downloaded here:
https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html

To log in you can write:

```bash
plink -P 8922 [user]@185.45.23.197
```

Where [user] is your user name. We have created an account for each of you with your first name (all lowercase) as user name and your surname (all lowercase) as password (you can use the command passwd to change the password).

And you can transfer files between the cluster and your local machine by running:

```bash
pscp -P 8922 [user]@185.45.23.197:[dir]/[file] .
```

In order to copy files from your local machine to the cluster, you just need to invers the order:

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


## Libre servers with Graphical User Interface (GUI)

For windows:
https://mobaxterm.mobatek.net/

For mac:
https://cyberduck.io/


