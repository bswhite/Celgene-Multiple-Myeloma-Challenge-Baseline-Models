FROM ubuntu
RUN apt-get update
RUN apt-get -y install r-base

## Install dependencies
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("plyr")'
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('data.table')"

## Copy files over to Docker image
COPY run-mm-sc3.R /run-mm-sc3.R
COPY score_sc3.sh /score_sc3.sh
COPY model-state-metadata.Rd /model-state-metadata.Rd
