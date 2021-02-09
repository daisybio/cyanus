FROM rocker/shiny:4.0.0
LABEL authors="Lisi Arend, Judith Bernett, Quirin Manz"

#install system packages
RUN apt-get update && apt-get install -y libxml2-dev libssl-dev libcurl4-openssl-dev default-jdk \
   r-cran-rjava \
   && apt-get clean \
   && rm -rf /var/lib/apt/lists/
#libmariadbclient-dev libpq-dev libv8-dev liblzma-dev libbz2-dev zlib1g-dev libncurses5-dev

#copy shiny app to work dir
RUN mkdir /srv/cytof_pipeline
WORKDIR /srv/cytof_pipeline
COPY . .
RUN echo "options(shiny.maxRequestSize=1000*1024^2)" > .Rprofile

#update shiny server conf and configure it to run cytof_pipeline in single app mode
RUN sed -i 's/site_dir \/srv\/shiny-server;/app_dir \/srv\/cytof_pipeline;/g' /etc/shiny-server/shiny-server.conf


#install R packages via renv
ENV RENV_VERSION 0.12.5-2
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org')); \ 
          remotes::install_github('rstudio/renv@${RENV_VERSION}'); \
          renv::restore()"

#fix bug where app greys out on SSL connection
RUN echo 'sanitize_errors off;disable_protocols xdr-streaming xhr-streaming iframe-eventsource iframe-htmlfile;' >> /etc/shiny-server/shiny-server.conf

EXPOSE 3838