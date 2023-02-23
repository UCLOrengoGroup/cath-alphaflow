
# Oracle Installation within UCL

Within UCL there is the option of accessing a database on the CS cluster. This may be either from UCL clusters such as the CS cluster or Myriad, or for local test and development.

There are two oracle client libraries that can be used in python. cath-alphaflow is currently using the old one, but upgrade to the new one is trivial and exampled in the test script. Upgrade to the new one has the future possible benefit of creating a thin client (no oracle installation needed locally) but for now a thick client (and installation of oracle client) is required so the upgrade has no immediate advantages.

## Installation

#### Oracle client
Using WSL: Instructions to setup Oracle client: https://medium.com/@arunkundgol/how-to-setup-oracle-instant-client-on-windows-subsystem-for-linux-cccee61d5b0b

General instructions: https://www.oracle.com/uk/database/technologies/instant-client.html

#### Install the python linrary

New library) https://python-oracledb.readthedocs.io/en/latest/index.html
Old library) https://oracle.github.io/python-cx_Oracle/

#### Settings to use
Note these are confidential to UCL

 ORACLE_DB_HOST="localhost"
 ORACLE_DB_PORT=1521
 ORACLE_DB_SID="ask someone in UCL"
 ORACLE_DB_USERNAME="ask someone in UCL"
 ORACLE_DB_PASSWORD="ask someone in UCL"

The settings go in the cath-alphaflow/settings.py file, replacing this function:
```
class ProductionSettings(Settings):
    pass
```
The file is in gitignore so you can safely add your settings locally.

#### Set up a tunnel so that you can access from wherever you are
We set up an ssh tunnel from localhost to the oracle database server via an intermediate server `gateway`, since the oracle server is not directly accessible from localhost The variables are `dbname` and `gatename` along with the server path. If you have direct access to the oracle server from your location then you do not need to use an ssh tunnel, but can put the oracle server's hostname directly into the connection settings for your local oracle client to use.

add this to your local .ssh/config file (replace someone with your cs cluster username)

```
Host gatename
  User someone
  HostName gatename.xx.xxx.ac.uk
```

set up a local ssh key if you don't already have one, and then copy it to gatename
```
$ ssh-keygen
$ ssh-copy-id gatename
<enter your cs cluster password>
```

you should now be able to log into gatename without entering a password (the first time you may need to enter the passphrase for you sshkey if you created one)
```
$ ssh gatename
(then exit again)
```

now try this
```
$ ssh -L 1521:dbname.xx.xx.ac.uk:1521 gatename
```
this should open a normal looking ssh session to gatename server, but also sets up the tunnel to the db in the back ground
(to close the tunnel exit the ssh gatename session)
while the special gatename session is open you should be able to talk to the oracle db on localhost port 1521


once that's all working, for convenience you can create an alias in your local .bashrc by adding this line somewhere
to avoid typing the long ssh command each time:
```
alias ssh-odb='ssh -N -L 1521:dbname.xx.xx.ac.uk:1521 gatename'
```
the -N tells ssh not to open the interactive session to gatename, but it still sets up the tunnel itself
now just use:
```
$ ssh-odb
and use CTRL-C to terminate the tunnel
```
("ssh -L 1521:dbname.xx.xxx.ac.uk:1521 gatename" means "connect localhost port 1521 to dbname.xx.xxx.ac.uk port 1521 via gatename.xx.xxx.ac.uk" )
