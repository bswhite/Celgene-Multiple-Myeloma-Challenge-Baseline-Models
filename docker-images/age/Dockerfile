FROM ubuntu
RUN apt-get update
RUN apt-get -y install r-base

## Install dependencies
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("plyr")'
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('data.table')"

## Copy files over to Docker image
COPY run-mm-sc1.R /run-mm-sc1.R
COPY score_sc1.sh /score_sc1.sh
COPY model-state-metadata.Rd /model-state-metadata.Rd
