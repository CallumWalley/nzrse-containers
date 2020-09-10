---
title: "GUI enabled applications: RStudio in a container"
teaching: 10
exercises: 10
questions:
objectives:
- Learn how to run GUI interactive sessions from a container
- Learn how to setup a long running web service from a container
keypoints:
- An interactive session can essentially be executed as any other containerised application, via `singularity exec`
- Use the `%startscript` section of a def file to configure an image for long running services
- Launch/shutdown long running services in the background with `singularity instance start/stop`
---


### R and RStudio images

R is a popular language in several domains of science, mostly because of its statistical packages. In particular it is nowadays highly common in data science and bioinformatics.

The group [Rocker](https://hub.docker.com/r/rocker) has published a large number of R images we can use, including an RStudio image.  To begin, let's cd into the appropriate directory:

```
cd $NZRSE/demos/08_rstudio
```
{: .bash}


> ## Pull the container
>
> We want to use a [Tidyverse](https://www.tidyverse.org) container image (contains R, RStudio, data science packages). Can you pull the  `rocker/> tidyverse:3.6.1` from Docker Hub?
>
> > ## Solution
> >
> > ```
> > singularity pull docker://rocker/tidyverse:3.6.1
> > ```
> > {: .bash}
> {: .solution}
{: .challenge}


### Running a scripted R workflow on the shell

To begin with, we are going to run a minimalistic example taken from the workshop [Programming with R](http://swcarpentry.github.io/r-novice-inflammation/) by the Software Carpentry. In particular, their [Episode 5](http://swcarpentry.github.io/r-novice-inflammation/05-cmdline/index.html) is the source for the dataset and the R script file; the latter has been adapted for this workshop.

Let us start with running the R script through the R container; we're going to compute average values in this example:

```
singularity exec -B $PWD tidyverse_3.6.1.sif Rscript readings-density.R --mean inflammation-density.png data/inflammation-*.csv
```
{: .bash}

```
5.45
5.425
6.1
[..]
6.875
6.025
6.9
Saving 7 x 7 in image
```
{: .output}

We even got a nice plot file out of the analysis, `inflammation-density.png`. So, what if want to run our R workflow using the RStudio GUI interface?


### Run an interactive RStudio session

Beside developing R container images, Rocker has also some useful documentation on [Running Rocker R container with Singularity](https://www.rocker-project.org/use/singularity/). Let's set this up together.

Lastly, containers are read-only, but RStudio will want to be able to write configuration and temporary files in the home. Let us bind mount the current work directory as the container home.  

```

singularity exec -c -B $PWD tidyverse_3.6.1.sif rserver --www-port 8787 --www-address 0.0.0.0 --auth-none=1 --auth-validate-users=0
```
{: .bash}

Note the `-c` flag for `singularity exec`, used to avoid sharing directories such as `/tmp` with the host, and thus to better clean up the session upon exit.  
If everything is fine, no output will be printed.

In the R console, submit the analysis script we ran earlier on from the shell:

```
source("readings-density.R")
```
{: .r}

If you have a look at the bottom right panel, you can see some outputs files are generated, including `interactive.png`. Click on it, and you'll get to visualise the resulting plot!

Once you're done, click on the power icon on the top right to close the session, then go back to the shell and kill the container with `Ctrl-C`.

As a final remark, note that the setup we just described could be adapted for use from a compute node in a HPC system, too, by using the HPC scheduler.


### Setup a long running RStudio web server

The procedure we just described can be convenient for interactive sessions of relatively short duration. On the other hand, if we wanted to deploy a long running RStudio server, having to keep the terminal open isn't really handy.

Singularity has features to run containers in background. To this end we're going to explore the subcommands of `singularity instance`.

If we need an image to be run as a background instance with Singularity, this needs to be build with a special section in the def file, namely `%startscript`. Commands in this section are executed when the instance is started. If no such section is provided, by default a shell will be executed, a bit useless for our RStudio server.

Building on the experience in the past paragraph, let us design a def file for the purpose (see `tidyverse_long.def` in the demo dir:

```
Bootstrap: docker
From: rocker/tidyverse:3.6.1

%labels
  Author Pawsey Supercomputing Centre
  Version 0.0.1

%startscript
  export R_PORT=${R_PORT:-"8787"}
  export R_ADDRESS=${R_ADDRESS:-"0.0.0.0"}
  
  rserver --www-port $R_PORT --www-address $R_ADDRESS --auth-none=1 --auth-validate-users=0
```
{: .source}

Basically, we're starting from the `tidyverse` Docker image we used above, and then adding some commands under the `%startscript` header. In particular, we're adding some flexibility to the `rserver <..>` command we used above, allowing for port and address to be redefined by the user through environment variables, and at the same time by providing sensible defaults.


> ## Build an image to run a RStudio instance
>
> How would you build an image called `tidyverse_long.sif`, starting from this def file?
>
> > ## Solution
> >
> > ```
> > singularity build -r tidyverse_long.sif tidyverse_long.def
> > ```
> > {: .bash}
> {: .solution}
{: .challenge}

Once the container image is build, let's use it to start an instance via `singularity instance start`. Note how the other options are the same as in the interactive session above; the only addition is the specification of a name for the instance, `myserver` in this case, that has to follow the image name:

```
export PASSWORD=password
echo $USER && echo $PASSWORD
export R_USER=$USER && [ "$(id -u)" == "1000" ] && export R_USER=rstudio

singularity instance start -c -B $(pwd):/home/$R_USER tidyverse_long.sif myserver
```
{: .bash}

```
INFO:    instance started successfully
```
{: .output}

We can check on the running instances with

```
singularity instance list
```
{: .bash}

```
INSTANCE NAME    PID      IMAGE
myserver         18080    /home/ubuntu/ernz20-containers/demos/08_rstudio/tidyverse_long.sif
```
{: .output}

Note that we can run commands from the instance by referring to it as `instance://<INSTANCE-NAME>`, *e.g.*

```
singularity exec instance://myserver echo $USER $PASSWORD
```
{: .bash}

Once we've finished with RStudio, we can shutdown the instance with

```
singularity instance stop myserver
```
{: .bash}

```
Stopping myserver instance of /home/ubuntu/ernz20-containers/demos/08_rstudio/tidyverse_long.sif (PID=18080)
```
{: .output}
