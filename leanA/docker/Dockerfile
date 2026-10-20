FROM ubuntu:24.04
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    curl ca-certificates git bubblewrap build-essential python3 && rm -rf /var/lib/apt/lists/*
RUN useradd -m -s /bin/bash checker && mkdir -p /run/user/1000 && chown checker /run/user/1000 && chmod 700 /run/user/1000
USER checker
WORKDIR /home/checker
RUN curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -o elan-init.sh && \
    sh elan-init.sh -y --default-toolchain leanprover/lean4:v4.35.0 && rm elan-init.sh
ENV PATH=/home/checker/.elan/bin:$PATH
RUN mkdir -p /home/checker/leanA/.lake /home/checker/.cache
RUN lean --version && ls /home/checker/.elan/toolchains/*/bin
