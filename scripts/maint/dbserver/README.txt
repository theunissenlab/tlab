Instructions for installing database and wiki
Channing Moore
12/9/08

1. Basics

The wiki comprises three, separate programs;
Apache, the web server program that handles http requests;
Mediawiki, the wiki program, which is written in PHP;
MySQL, the database that actually contains and serves the data for the wiki's pages.

To install the wiki, you need all three of these things running. For this reason these instructions deal not only with reinstalling the wiki, but the database as well.

2. Computers and Names

Currently, both the wiki and the database are served from my desktop box, roche.fet.berkeley.edu. They don't have to both be on the same box, bur right now, they are.

For ease of use, the wiki resolves on www.fet.berkeley.edu/mediawiki

and the database resolves on 

database.fet.berkeley.edu:3306

To make this happen, I had to modify the nameserver configuration files on nutmeg so that both www.fet and database.fet resolve to roche.fet. If you want to change either of these programs to run on a different computer, you'll have to change those configurations.

3. Backing Up the Database

The first thing you should do is make a backup of the current database state if possible. If the database is still working, run the script

/auto/fhome/channing/scripts/db-backup

This will back up the entire database (including my neurophysiology data) to the fileserver. Be patient, it takes 10 or 20 minutes to run. This runs as a cron job weekly, so there is always a recent backup. The files are named by date.

3. Full Installation

Now that you've backed up the database (you did that already, right?), you can proceed.

There is a script in my space on the fileserver that will do a bare-bones reinstall of the wiki and database on one box:

/auto/fhome/channing/server_install/server-install

I've commented it pretty well, so if you need details, just go read the script. It makes extensive use of the unix command 'sed' to make small changes to configuration text files, and there are several sed script files in the same directory that contain the actual code for modifications.

Be careful with this: I have no idea what would happen if you ran it on an existing install. This script does a bare-bones install on Ubuntu, fully reinstalling Apache, MySQL, PHP, and their connection libraries, and tweaking the configuration files for each.

This script restores the entire database (including neurophysiology data) from the backup file on the fileserver, and this step could take 20 minutes. Be patient. If you need to do a full reinstall, you should edit the script file and change the name of the file it restores to the most recent one.

4. Partial Reinstallation

For simpler problems you might only need part of the reinstallation script. Like I said, be careful and don't just run the whole script unless you're doing a clean install.