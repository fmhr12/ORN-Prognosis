FROM rocker/shiny:4.4.2

# Install additional system dependencies if needed.
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev && \
    rm -rf /var/lib/apt/lists/*

# Install extra R packages that are not already included in the image.
RUN R -e "install.packages(c('shinythemes', 'survival', 'riskRegression', 'fastshap', 'shapviz', 'ggplot2', 'prodlim', 'gower', 'plotly'), repos='https://cran.rstudio.com/')"

# Copy the app files into the image.
COPY . /app
WORKDIR /app

# Expose port 10000 (Render's default port).
EXPOSE 10000

# Run the Shiny app.
CMD R -e "shiny::runApp('/app', host='0.0.0.0', port=as.numeric(Sys.getenv('PORT', 10000)))"
