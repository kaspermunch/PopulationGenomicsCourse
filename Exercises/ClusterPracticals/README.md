# Linux server for exercises
During the course there will be several exercises that use command line tools that are only available on linux machines. These programs have been installed on a linux server that you can log onto using ssh.

## Login from linux of mac
 If you are on a linux or mac computer then ssh is usually installed by default and you can use the following command to log in to the server:

```bash
ssh -p 8922 [user]@185.45.23.197
```

Where [user] is your user name. We have created an account for each of you with your first name (all lowercase) as user name and password = 1234[user] (you can use the command passwd to change the password).

To copy a file called [file] in a directory called [dir] on the server to your local machine you can use the scp command (you can also use sftp if you prefer that):

```bash
scp -P 8922 [user]@185.45.23.197:[dir]/[file] .
```

## Login from windows
For Windows users:
If you are on a windows machine you need to have PuTTy or Cygwin installed to use ssh. PuTTy can be downloaded here:
https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html
If you use PuTTy the ssh command is called plink (not to be confused with the PLINK program we will use to do the GWAS analyses), so to log in write:

```bash
plink -P 8922 [user]@185.45.23.197
```

And the command to copy from the server is called pscp:

```bash
pscp -P 8922 [user]@185.45.23.197:[dir]/[file] .
````

## Using the command line
If you are not used to using a unix command-line interface you might find this introduction useful:
https://lifehacker.com/5633909/who-needs-a-mouse-learn-to-use-the-command-line-for-almost-anything


## Libre servers with Graphical User Interface (GUI)

For windows:
https://mobaxterm.mobatek.net/

For mac:
https://cyberduck.io/


