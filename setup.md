---
layout: page
title: Setup
root: .
---

The main requirement for this workshop is a personal computer with a web browser and a command line shell program.  
These will allow you to follow the online materials and to login to a facility with the required software stack.


> ## eRNZ20 Workshop Attendee Instructions
>
> This tutorial is adapted from [https://github.com/PawseySC/sc19-containers](https://github.com/PawseySC/sc19-containers) and delivered jointly by New Zealand eScience Infrastructure (NeSI) and Pawsey Supercomputing Centre staff.
>
> The tutorial will be best experienced directly on NeSI's platforms. Though it is possible to follow along on a laptop we are not providing a standard environment for this and so will have limited ability to troubleshoot machine/environment specific issues during the tutorial. All the required software and container images have been preloaded on the NeSI platform under the `nesi99991` project.
>
> If you do not already have a NeSI account then please apply to [create a NeSI account](https://support.nesi.org.nz/hc/en-gb/articles/360000159715-Creating-a-NeSI-Account) well ahead of the tutorial.
>
> Once you have a NeSI account you will need membership in the relevant training project. Please contact NeSI Support via either the [form here](https://support.nesi.org.nz/hc/en-gb/requests/new) or [email](mailto:support@nesi.org.nz) and provide your NeSI username (subject: "container tutorial"). Updating project membership may take up to 15 minutes, so please do this ahead of the tutorial if possible.
>
> If you have a NeSI user account and access to the training project then you don't need to read further for today.
{: .callout}


### Regular users of this tutorial: read here

* **Mahuika @NeSI**: if you have access, a recent version of Singularity can be loaded with `module load singularity`.

* **BYO Device**: if you have a Linux box, you can install the required software yourself (might take a while):

  * Essential (core of the tutorial)
    - Singularity : [script]({{ page.root }}/files/install-singularity.sh) \| [docs](https://sylabs.io/guides/3.4/user-guide/installation.html)

  * Desirable (to run all Singularity examples)
    - MPICH library : [script]({{ page.root }}/files/install-mpich.sh) \| [docs](https://www.mpich.org/documentation/guides/)
    - Nvidia GPU driver (GPU card required)
    - Slurm scheduler
    - Nextflow engine : [script]({{ page.root }}/files/install-nextflow.sh) \| [docs](https://www.nextflow.io/docs/latest/getstarted.html)

  * Optional (extra applications for last two episodes)
    - Docker : [docs (unofficial)](https://www.itzgeek.com/how-tos/linux/ubuntu-how-tos/how-to-install-docker-on-ubuntu-18-04-lts-bionic-beaver.html)
    - HPCCM : [script]({{ page.root }}/files/install-hpccm.sh) \| [docs](https://github.com/NVIDIA/hpc-container-maker/blob/master/docs/getting_started.md)
    - Podman : [script]({{ page.root }}/files/install-podman.sh) \| [docs](https://podman.io/getting-started/installation)
    - Sarus : [script]({{ page.root }}/files/install-sarus.sh) \| [docs](https://sarus.readthedocs.io/en/latest/install/requirements.html)

**Notes**
* Install scripts have been tested on a Ubuntu machine through a user that can run *sudo* commands without password prompts. There's no warranty they will work in your Linux box, you should consider them as templates.
* To install Singularity on a Mac or Windows machine, follow these instructions by Sylabs on [Setting up Singularity with Vagrant](https://sylabs.io/guides/3.4/user-guide/installation.html#install-on-windows-or-mac).
