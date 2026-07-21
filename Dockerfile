# HGVS -> VCF decoder (Ruby serving image).
# Builds/serves from a mounted /data volume; the LMDB index and reference genome
# are NOT baked in — you prepare them on the host (see README) and mount them.
FROM ruby:3.2-slim-bookworm

RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential liblmdb-dev samtools gzip ca-certificates \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY Gemfile ./
RUN bundle config set --local without 'development test' && bundle install

COPY lib ./lib
COPY bin ./bin
COPY app.rb ./app.rb
COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
RUN chmod +x /usr/local/bin/docker-entrypoint.sh

ENV HOME=/tmp \
    DATA=/data \
    INDEX_DIR=/data/index \
    GENOME_FASTA=/data/GRCh38.fna \
    BUILD=GRCh38 \
    PORT=4567
EXPOSE 4567
VOLUME ["/data"]
ENTRYPOINT ["docker-entrypoint.sh"]
CMD ["serve"]
