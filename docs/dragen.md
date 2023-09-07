## Setting up a Dragen run on MacOS

First, download the basescape binary into your home directory.

```
wget "https://launch.basespace.illumina.com/CLI/latest/amd64-osx/bs" -O $HOME/bin/bs
```

Make sure that the downloaded binary file is actually executable.

```
chmod u+x $HOME/bin/bs
```

Update basespace CLI HUB using Homebrew.

```
brew tap basespace/basespace && brew install bs-cli
```

Run the following to prove that it isuccessfully installed

```
which bs 
/usr/local/bin/bs
```

Authentificate on the command line

```
bs auth
```

If you configured it before, it may give you an error with already exists message

```
bs auth --force
```

`bs` will print a message like the one below on your terminal, simply follow the instruction in that message.

```
Please go to this URL to authenticate:  

Copy the link and paste it in the browser, you will need to login and accept the terms and conditions, then you will see the welcome message below.
https://basespace.illumina.com/oauth/device?code=s6BKU
```

After pasting the link in your browser, you will receive a message welcoming you to the basespace dashbord, something like **Welcome, YourFirstname YourLastName**.

A quick note: Basespace tracks data by project and biossamples. As such, we need to specify them before we can go on to upload data. Now, let us move back to the terminal. Just a quick note
create a bio sample and corresponding project where we shall store the data, in the command below; bio sample name is  `outbreak` and project name is `VHFnegative`.

```
bs create biosample -n outbreak -p VHFnegative
```

To upload data, we need to keep track of the projectID, which we can inspect for using the `bs` command below.

`bs list projects`

You get an output on screen looking something like this.

```
+---------------+-----------+-----------+

|     Name      |    Id     | TotalSize |

+---------------+-----------+-----------+

| CovidSeq Nick | 374662296 | 283558309 |

| SARSCoV2      | 378713340 | 1490330   |

| VHFnegative   | 387386008 | 0         |

+---------------+-----------+-----------+
```

Upload the data to the specified project ID, from the table above, we pick the project ID that corresponds to `VHFnegative`.
Note: We can choose to run make the upload within the directory where the data is kept (as in the example below). Otherwise, we need to provide a path to the data that we intend to upload.

```
bs upload dataset -p 387386008 B*-04-00*.gz
```

If all goes well, these messages will be printed on your terminal showing the progress of the data upload.

```
Creating sample: BVHF23-04-001

BVHF23-04-001_S1_L001_R1_001.fastq.gz  176.89 MiB / 176.89 MiB [==============================================] 100.00%

BVHF23-04-001_S1_L001_R2_001.fastq.gz  190.74 MiB / 190.74 MiB [==============================================] 100.00%

Creating sample: BVHF23-04-002

BVHF23-04-002_S2_L001_R2_001.fastq.gz  121.09 MiB / 121.09 MiB [==============================================] 100.00%

BVHF23-04-002_S2_L001_R1_001.fastq.gz  118.43 MiB / 118.43 MiB [==============================================] 100.00%
```

Once all the data is uploaded, get back to your browser on the basespace dashboard and initiate the analysis.
