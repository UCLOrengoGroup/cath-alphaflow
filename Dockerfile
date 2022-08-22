FROM apache/airflow:2.3.3
USER root
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
         vim wget gzip tar\
  && apt-get autoremove -yqq --purge \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*
USER airflow
