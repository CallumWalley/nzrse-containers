---
title: "Share files with the host: BLAST, a bioinformatics example"
teaching: 10
exercises: 10
questions:
objectives:
- Learn how to mount host directories in a container
- Run a real-world bioinformatics application in a container
keypoints:
- By default Singularity mounts the host current directory, and uses it as container working directory
- Map additional host directories in the containers with the flag `-B`
---


### Access to host directories

Let's start and `cd` into the root demo directory:

```
cd $NZRSE/demos
```
{: .bash}

What directories can we access from the container?

First, let us assess what the content of the root directory `/` looks like outside *vs* inside the container, to highlight the fact that a container runs on his own filesystem:

```
ls /
```
{: .bash}

```
bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  scratch  shared  srv  sys  tmp  usr  var
```
{: .output}


Now let's look at the root directory when we're in the container
```
singularity exec library://ubuntu:18.04 ls /
```
{: .bash}

```
bin  boot  data  dev  environment  etc	home  lib  lib64  media  mnt  opt  proc  root  run  sbin  singularity  srv  sys  tmp  usr  var
```
{: .output}


> ## In which directory is the container running?
>
> For reference, let's check the host first:
>
> ```
> $ pwd
> ```
> {: .bash}
>
> ```
> /nesi/nobackup/nesi99991/${USER}/nzrse-containers/demos
> ```
> {: .output}
>
> Now inspect the container (HINT: you need to run `pwd` in the container)
>
> > ## Solution
> >
> > ```
> > $ singularity exec library://ubuntu:18.04 pwd
> > ```
> > {: .bash}
> > ```
> > /scale_wlg_nobackup/filesets/nobackup/nesi99991/${USER}/nzrse-containers/demos
> > ```
> > {: .output}
> >
> > Host and container working directories coincide! Noting that on NeSI the `/nesi/nobackup/` path is a link, this gets resolved when we start the container. Multiple symlinks in your pwd hierarchy can cause issues, if this happens you'll typically end up with your container pwd in your home dir.
> {: .solution}
{: .challenge}


> ## Can we see the content of the current directory inside the container?
>
> Not by default. Although this was the default behaviour with previous Singularity versions, Singularity v3.6 or above will neither recreate the current working directory (CWD) path inside container nor resolve any symlinks in CWD path to determine the destination path inside container. If so, what is the default path it refers to ?
>
> > ## Solution
> >
> > ```
> > $ singularity exec library://ubuntu:18.04 ls
> > ```
> > {: .bash}
> >
> > ```
> > /home/${USER}
> > ```
> > {: .output}
> >
> {: .solution}
{: .challenge}


> ## How about other directories in the host?
>
> For instance, let us inspect our training project-folder: `/nesi/project/nesi99991/`.
>
> > ## Solution
> >
> > ```
> > $ singularity exec library://ubuntu:18.04 ls /nesi/project/nesi99991/
> > ```
> > {: .bash}
> >
> > ```
> > ls: cannot access '/nesi/project/nesi99991/': No such file or directory
> > ```
> > {: .output}
> >
> > Typically, host directories external to the current directory are not visible! Though in some cases Singularity may be configured system-wide to ensure key shared storage is available by default. Assuming this path isn't visible in this case, how can we fix this? Read on...
> {: .solution}
{: .challenge}


### Bind mounting host directories

Singularity has the runtime flag `--bind`, `-B` in short, to mount host directories.

The long syntax allows to map the host dir onto a container dir with a different name/path, `-B hostdir:containerdir`.  
The short syntax just mounts the dir using the same name and path: `-B hostdir`.

Let's use the latter syntax to mount `/nesi/project/nesi99991` into the container and re-run `ls`.

```
singularity exec -B /nesi/project/nesi99991 library://ubuntu:18.04 ls /nesi/project/nesi99991/
```
{: .bash}

```
hi  there
```
{: .output}

Now we are talking!

If you need to mount multiple directories, you can either repeat the `-B` flag multiple times, or use a comma-separated list of paths, *i.e.*
```
-B dir1,dir2,dir3
```
{: .bash}

Also, if you want to keep the runtime command compact, you can equivalently specify directories to be bind mounting using an environment variable:
```
export SINGULARITY_BINBPATH="dir1,dir2,dir3"
```
{: .bash}


> ## Mounting $HOME
>
> Depending on the site configuration of Singularity,
> user home directories might or might not be mounted into containers by default.
> We do recommend avoid mounting home whenever possible,
> to avoid sharing potentially sensitive data, such as SSH keys, with the container,
> especially if exposing it to the public through a web service.
>
> If you need to share data inside the container home, you might just mount that specific file/directory, *e.g.*
> ```
> -B $HOME/.local
> ```
> {: .bash}
>
> Or, if you want a full fledged home, you might define an alternative host directory to act as your container home, as in
> ```
> -B /path/to/fake/home:$HOME
> ```
{: .callout}


### Running BLAST from a container

We'll be running a BLAST (Basic Local Alignment Search Tool) example with a container from [BioContainers](https://biocontainers.pro). BLAST is a tool bioinformaticians use to compare a sample genetic sequence to a database of known sequences; it's one of the most widely used bioinformatics tools.  
This example is adapted from the [BioContainers documentation](http://biocontainers-edu.biocontainers.pro/en/latest/running_example.html).

We're going to use the BLAST image `biocontainers/blast:v2.2.31_cv2`, which we previously downloaded in `$SIFPATH`. Let's verify it's there:

```
ls $SIFPATH/blast*
```
{: .bash}

```
/nesi/nobackup/nesi99991/bbet740/nzrse-containers/demos/sif/blast_v2.2.31_cv2.sif
```
{: .output}

> ## Run a test command
>
> To begin, let us run a simple command using the downloaded image `$SIFPATH/blast_v2.2.31_cv2.sif`,
> for instance `blastp -help`, to verify that it actually works:
>
> > ## Solution
> >
> > ```
> > $ singularity exec $SIFPATH/blast_v2.2.31_cv2.sif blastp -help
> > ```
> > {: .bash}
> >
> > ```
> > USAGE
> >   blastp [-h] [-help] [-import_search_strategy filename]
> > [..]
> >  -use_sw_tback
> >    Compute locally optimal Smith-Waterman alignments?
> > ```
> > {: .output}
> {: .solution}
{: .challenge}


Now, the `demos/03_blast` demo directory contains a human prion FASTA sequence, `P04156.fasta`, whereas another directory, `demos/03_blast_db`, contains a gzipped reference database to blast against, `zebrafish.1.protein.faa.gz`. Let us `cd` to the latter directory and uncompress the database:

```
cd $NZRSE/demos/03_blast_db
gunzip zebrafish.1.protein.faa.gz
```
{: .bash}


> ## Prepare the database
>
> We then need to prepare the zebrafish database with `makeblastdb` for the search, using the following command through a container:
>
> ```
> $ makeblastdb -in zebrafish.1.protein.faa -dbtype prot
> ```
> {: .bash}
>
> Try and run it via Singularity.
>
> > ## Solution
> >
> > ```
> > $ singularity exec -B $PWD $SIFPATH/blast_v2.2.31_cv2.sif makeblastdb -in zebrafish.1.protein.faa -dbtype prot
> > ```
> > {: .bash}
> > ```
> > Building a new DB, current time: 11/16/2019 19:14:43
> > New DB name:   /home/ubuntu/nzrse-containers/demos/03_blast_db/zebrafish.1.protein.faa
> > New DB title:  zebrafish.1.protein.faa
> > Sequence type: Protein
> > Keep Linkouts: T
> > Keep MBits: T
> > Maximum file size: 1000000000B
> > Adding sequences from FASTA; added 52951 sequences in 1.34541 seconds.
> > ```
> > {: .output}
> {: .solution}
{: .challenge}


After the container has terminated, you should see several new files in the current directory (try `ls`).  
Now let's proceed to the final alignment step using `blastp`. We need to cd into `demos/03_blast`:

```
cd ../03_blast
```
{: .bash}


> ## Run the alignment
>
> and then adapt the following command to run into the container:
>
> ```
> $ blastp -query P04156.fasta -db $NZRSE/demos/03_blast_db/zebrafish.1.protein.faa -out results.txt
> ```
> {: .bash}
>
> Note how we put the database files in a separate directory on purpose, so that you will need to bind mount its path with Singularity. Give it a go with building the syntax to run the `blastp` command.
>
> > ## Solution
> >
> > ```
> > $ singularity exec -B $PWD,$NZRSE/demos/03_blast_db $SIFPATH/blast_v2.2.31_cv2.sif blastp -query P04156.fasta -db $NZRSE/demos/03_blast_db/zebrafish.1.protein.faa -out results.txt
> > ```
> > {: .bash}
> {: .solution}
{: .challenge}

The final results are stored in `results.txt`:

```
less results.txt
```
{: .bash}

```
                                                                      Score     E
Sequences producing significant alignments:                          (Bits)  Value

  XP_017207509.1 protein piccolo isoform X2 [Danio rerio]             43.9    2e-04
  XP_017207511.1 mucin-16 isoform X4 [Danio rerio]                    43.9    2e-04
  XP_021323434.1 protein piccolo isoform X5 [Danio rerio]             43.5    3e-04
  XP_017207510.1 protein piccolo isoform X3 [Danio rerio]             43.5    3e-04
  XP_021323433.1 protein piccolo isoform X1 [Danio rerio]             43.5    3e-04
  XP_009291733.1 protein piccolo isoform X1 [Danio rerio]             43.5    3e-04
  NP_001268391.1 chromodomain-helicase-DNA-binding protein 2 [Dan...  35.8    0.072
[..]
```
{: .output}

We can see that several proteins in the zebrafish genome match those in the human prion (interesting?).
