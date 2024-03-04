FROM rocker/shiny:4.2.1
LABEL authors="Lisi Arend, Judith Bernett, Quirin Manz"

#install system packages
RUN apt-get update && apt-get install -y git curl cmake libharfbuzz-dev libfribidi-dev libtiff-dev libxml2-dev libssl-dev libcurl4-openssl-dev default-jdk libgmp3-dev\
   r-cran-rjava \
   && apt-get clean \
   && rm -rf /var/lib/apt/lists/
#libmariadbclient-dev libpq-dev libv8-dev liblzma-dev libbz2-dev zlib1g-dev libncurses5-dev

#copy shiny app to work dir
RUN mkdir /srv/cyanus-dev
WORKDIR /srv/cyanus-dev

#install R packages via renv
ENV RENV_VERSION v1.0.3
ENV RENV_WATCHDOG_ENABLED = FALSE
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

COPY renv.lock.shiny renv.lock
RUN R -e "renv::restore()"

# copy shiny app
COPY ./*.R ./
COPY ./functions ./functions
COPY ./ui ./ui
COPY ./server ./server
COPY ./www ./www
COPY ./www/favicon.png ./
COPY ./data ./data

RUN echo "options(shiny.maxRequestSize=5000*1024^2)" > .Rprofile

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp(host = '0.0.0.0', port = 3838)"]
