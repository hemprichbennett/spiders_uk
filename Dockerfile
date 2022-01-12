FROM rocker/tidyverse

ENV TZ=Europe/Moscow \
    DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt upgrade -yq
RUN apt-get -yq install software-properties-common wget dirmngr gnupg apt-transport-https curl libxml2-dev libcurl4-openssl-dev libssl-dev cron nano git


RUN Rscript -e 'install.packages(c("vegan", "gridExtra", "forcats", "rgbif", "ggbeeswarm", "ggthemes", "picante", "dada2"))'