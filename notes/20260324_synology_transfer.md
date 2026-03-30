We used FileZilla to transfer the data from Novogene to the Gatins Lab Synology file storage. 

On the Gatins lab desktop computer, we opened FileZilla and found the local site:

```
~/Volumes/Projects/2025_MobulidPopGen_MGE/rawdata
```

and then copy/pasted the Host, Username (Account Number), Password, and Port from the Novogene email. 

Then drag and drop the smallest files first, make sure the transfer is successful. then move over the larger files


~on explorer~ I was able to check the amount of storage remaining in the /projects/gatins folder with: 

```
(base) [eppley.m@explorer-01 gatins]$ df -h /projects/gatins
Filesystem                              Size  Used Avail Use% Mounted on
vast1-mghpcc-eth.neu.edu:/work_project   36T   32T  4.1T  89% /projects
(base) [eppley.m@explorer-01 gatins]$ 
```

However when I went back into the desktop to do the transfer to discovery, the synology storage was not accessible, so I think Remy will have to log in again. 

To connect remotely to discovery, I used the following: 
```
Host:  sftp://xfer.discovery.neu.edu
Username: eppley.m
Password: my northeastern password
Port: 22
```
then press quick connect
