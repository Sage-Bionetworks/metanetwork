# README

R package to build many different statistical networks on high dimensional multivariate data and run diagnostic statistics

### What is this repository for? 

* Metanetworks Package Setup
** Docker Container Setup and Package Install
* Running Metanetworks Package - Network Construction
* Running Metatnetworks Package - Module Construction
* Runing Metanetworks Package - Consensus Network Construction
* Version

### How do I get set up?

Pull the metanetworks package repo with the dockerfile and sample configuration files

```
cd <Project Working Directory>
git clone https://github.com/Sage-Bionetworks/metanetwork.git

# create support directories 
mkdir temp
mkdir out
mkdir error
```

The dockerfile to setup the proper run environment and install the metanetwork package is located within the pacakge in `metanetwork/inst/Docker/`.
To build the docker image with the proper environment run the following command from the directory you cloned the pacakge contents into.

```
docker build -t metanets metanetwork/inst/Docker/
``` 

Start a container in detached continual run mode to pass build comands into. Use the `-v` command to pass your current working directory
 into the current into the container. Any file or folder within `<Project Working Directory>` will be available within the container under
`/root/`. Conversly any file written under `/root/` inside the container will be written to `<Project Working Directory>` outside the container.
This makes checking log files and diagnosising any potential issues much easier. The container name will be `metanets` as specified by
`--name`. 

```
docker run -it -d -v "/<Project Working Directory>/:/root/" --name networks metanets
```

Basic commands to manipulate containers

```
# See running containers
docker container ls

# See all containers
docker ps -a

# Stop the container
docker stop <container name>
#eg.
docker stop metanets

# Start the container
docker start <container name>
#eg,
docker start metanets

# What if I want an interactive terminal session inside the container?
docker exec -it <container name> /bin/bash
cd ~
#eg.
docker exec -it metanets /bin/bash
cd /root/
```

### Building Individual Networks

Each network is build by passing 3 parameters into an R script called `Network_Wrapper.R`. This script is located within the metanets package within 
`metanetwork/inst/runscripts/`. Since this file was in the mounted `<Project Working Directory>` volume it is accessible inside the container with the filepath 
`/root/metanetwork/inst/runscripts/Network_Wrapper.R`. The run parameters required for the script are as follows: 
``` 
	-u, --synapse_user	Synapse User name
	-p, --synapse_pass	Synapse User Password
	-c, --config_file	Path to the complete config file
```

While the first two arguments are straight forward the third is more involved. The configuration YAML file specifies the details of your network construction. For 
detailed documentaition on this file see: `metanetwork/inst/confighelp/Network_Wrapper_Config_Documentation.Rmd`. For example configuration files see: 
`metanetwork/inst/config/network-construction/`.

Each network type is built with a single run command passed into the docker container that is already running in detached mode. An example run command from
a sample configuration file is:

```
docker exec -itd networks /bin/sh -c "export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib && 
	export PATH=/usr/lib64/openmpi/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/lib64/openmpi/bin && 
mpiexec --allow-run-as-root --mca orte_base_help_aggregate 0 --mca btl_base_warn_component_unused 0 -np 1 Rscript /root/metanetwork/inst/runscripts/Network_Wrapper.R 
-u <Synapse User Name> -p <Synapse User Password> -c /root/metanetwork/inst/config/network-construction/c3net.yml > /root/log.log"
```

The first two commands set the proper OpenMPI libraries and executable paths into the container's `LD_LIBRARY_PATH` and `PATH` environment variables.
The final command executes the network wrapper on the configuration file through implementation with OpenMPI. OpenMPI subdivides the compute environment into sub-cores
to more effciently process the network construction. NOTE: It is highly adviseable to only run one network at a time on your compute environment since OpenMPI is not
designed to run on top of itself. You can check the run progress of your network in the `log.log` file. Again because you imported your working envionment into the
container this file is available at `<Project Working Directory>/log.log` and updates live from within the running container.

If you want to run multiple networks in parallel we recomend deploying independent machines, either seperate HPC clusters or cloud based instances such as amazon EC-2
instances. This will also allow you to requisition the approriately sized machines for a given network approach. For benchmarking and expectations we have run an RNA-Seq dataset
of 750 samples and 9,740 genes on a machine consisting of 16 cores (vCPU) and 128 GiB (Memory) amazon EC-2 instance (r5.4xlarge) and the runtimes were:
```
Light Networks: c3net, mrnet, wgcna: (estimated less than 3 Hours)
Medium Networks: lassoAIC, lassoBIC, lassoCV1se, lassoCVmin, ridgeAIC, ridgeBIC, ridgeCV1se, ridgeCVmin, sparrowZ, sparrow2Z: (estimated around 10-12 Hours)
Heavy Networks: genie3, tigress: (TBD)
``` 

### Module Construction


### Who do I talk to? 

* Repo owner or admin
* Other community or team contact
