FROM rstudio/r-base:3.6.3-centos7

RUN yum update -y \
	&& echo "export PS1=\"\n\[$(tput bold)\]\[$(tput sgr0)\]\[\033[38;5;38m\]\u\[$(tput sgr0)\]\[$(tput sgr0)\]\[\033[38;5;15m\]@\h[\[$(tput sgr0)\]\[\033[38;5;214m\]\W\[$(tput sgr0)\]\[\033[38;5;15m\]]: \[$(tput sgr0)\]\"" >> ~/.bashrc \ && . ~/.bashrc \ 
	&& mkdir /shared \
	&& mkdir /shared/rlibs \
	&& chown root /shared/rlibs \
	&& yum install -y epel-release \
	&& yum install -y openmpi-1.10.7 \
	&& yum install -y openmpi-devel-1.10.7 \
	&& yum install -y git \
	&& yum install -y nano \
	&& yum -y install libcurl libcurl-devel \ 
	&& yum -y install openssl-devel \
	&& yum -y install libxml2-devel \
	&& yum install -y libpng-devel \
	&& yum install -y libjpeg-turbo-devel \
	&& yum -y groupinstall 'Development Tools' \
	&& yum -y install zlib-devel bzip2-devel openssl-devel ncurses-devel sqlite-devel readline-devel tk-devel \
	&& yum -y remove java \ 
	&& yum install -y java-1.8.0-openjdk \
	&& yum install -y java-1.8.0-openjdk-devel \
	&& yum install -y maven \ 
	&& yum install -y curl \
	&& yum install -y unzip \
	&& yum install -y screen \ 
	&& yum install -y centos-release-scl \
	&& yum install -y devtoolset-8-gcc* \
	&& yum install -y wget

RUN echo "local({r <- getOption('repos'); r['CRAN'] <- 'http://cran.r-project.org' ; options(repos=r)})" >> /opt/R/3.6.3/lib/R/etc/Rprofile.site

# Rmpi install
ENV LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:${LD_LIBRARY_PATH}
ENV PATH=/usr/lib64/openmpi/bin/:${PATH}

RUN R -e "install.packages('http://cran.rstudio.com/src/contrib/Archive/rjson/rjson_0.2.13.tar.gz', repos=NULL, type='source')" \
	&& R -e "options(repos = c('http://cran.r-project.org', 'http://cran.fhcrc.org', 'http://ran.synapse.org'))" \
	&& R -e "install.packages('Rmpi',configure.args=c(Rmpi='--with-Rmpi-libpath=/usr/lib64/openmpi/lib --with-Rmpi-include=/usr/include/openmpi-x86_64/ --with-Rmpi-type=OPENMPI'))"

# RandomForest
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-12.tar.gz', repos=NULL, type='source')" 

# Glmnet
RUN scl enable devtoolset-8 bash \
	&& R -e "dotR <- file.path(Sys.getenv('HOME'), '.R'); if (!file.exists(dotR)) dir.create(dotR); M <- file.path(dotR, 'Makevars'); if (!file.exists(M)) file.create(M); cat('\nCXX14FLAGS=-O3 -march=native -mtune=native -fPIC', 'CXX14=/opt/rh/devtoolset-8/root/usr/bin/g++', file = M, sep = '\n', append = TRUE)" \
	&& R -e "install.packages('glmnet')"

# Synapser
RUN R -e "install.packages('PythonEmbedInR', repos=c('http://cran.fhcrc.org', 'http://ran.synapse.org'))" \
	&& R -e "install.packages('synapser', repos=c('http://ran.synapse.org', 'http://cran.fhcrc.org'))"

# Githubr
RUN R -e "install.packages('devtools')" \
	&& R -e "devtools::install_github('brian-bot/githubr')"
# Metanets
RUN R -e "devtools::install_github('PradoVarathan/metanetwork')"
