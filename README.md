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

Note that if you're going to make a consensus network the `output_name:` configuration variable is very important. All networks need to have a common string following the network type. This is required for programatic location of the induvidual networks within a given synapse Project/folder.

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
of 750 samples and 9,740 genes on a machine consisting of; 16 cores (vCPU) and 128 GiB (Memory) amazon EC-2 instance (r5.4xlarge) [Light and Medium Networks], 126 cores (vCPU) and 256 GiB (Memory) amazon EC-2 instance (c6i.32xlarge) [Heavy Networks] and the runtimes were:
```
Light Networks: c3net, mrnet, wgcna: (estimated less than 3 Hours)
Medium Networks: lassoAIC, lassoBIC, lassoCV1se, lassoCVmin, ridgeAIC, ridgeBIC, ridgeCV1se, ridgeCVmin, sparrowZ, sparrow2Z: (estimated around 10-12 Hours)
Heavy Networks: genie3, tigress: (Genie3 ~11 hours, tigress ~96 hours)
``` 
### Consensus Network Construction

Constructing a Consensus Network from the ensemble of induvidual networks is preformed with an independent configuration yaml file and a wrapper R script. A sample configuration file can be found in `inst/config/network-consensus/` and the wrapper R script is located in `inst/runscripts/Consensus_Wrapper.R `. Breifly, this wrapper script scrapes networks from a synapse Project/Folder ID. The networks are found programtically as `Network_Wrapper.R ` has named the folders containing the network as the type of network contained ie. `c3net/` contains the c3net coexpression network file. The induvidual file is identified programatically with config specified ID string(s) `pattern_id` for most networks and `run_id` for the WGCNA netwoorks (see the config read me `inst/confighelp/Consensus_Wrapper_Config_Documentation.Rmd`). These identifiers correspond to how you named the network with the `output_name` setting in the network construction yaml file. For example:

```
# If your Network(s) were saved as:
#...
output_profile:
    output_path: /root/out/
    md5_output_path: microglia_lassoCV1semd5.out
    output_name: 'microglia_lassoCV1se_Network'
    error_path: /root/error
#...

# Your Consensus configuration specification should be:
#...
pattern_id: _Network
#...
```

To run the consensus Network you need to pass the command into the runing docker container similar to before. By importing your working directory into the container you've ensured that if the configuration file is inside the container as long as it is within your working directory. Simply adapt the following command to run the consensus network:

```

docker exec -itd networks /bin/sh -c "Rscript /root/metanetwork/inst/runscripts/Consensus_Wrapper.R -u <Synapse User Name> -p <Synapse User Password> -c /root/metanetwork/inst/config/network-consensus/microglial_consensus.yml > /root/log.log"

```

Consensus network construction on the Cerebellum light and medium networks ( all networks except genie and tigress) took ~96-120 hours on a 126 cores (vCPU) and 256 GiB (Memory) amazon EC-2 instance (c6i.32xlarge).

### Module Construction

Processes the consensus network though multiple sets of module identification software to identify co-expression modules within the consensus network. Testing is still in Beta. 

```

docker exec -itd networks /bin/sh -c "Rscript /root/metanetwork/inst/runscripts/Module_Wrapper.R -u <Synapse User Name> -p <Synapse User Password> -c /root/metanetwork/inst/config/network-module/module.yml > /root/log.log"

```

### Who do I talk to? 

* Repo owner or admin
* Other community or team contact
