#!/bin/bash
set -euo pipefail

docker build -t fastpli-ci - <./CI/dockerfile

# clone current commit state
rm -rf /tmp/fastpli-ci
git clone --recursive . /tmp/fastpli-ci

# generate container
if docker container inspect fastpli-cont-ci &>/dev/null; then
   docker stop fastpli-cont-ci
   docker rm fastpli-cont-ci
fi
docker create --name fastpli-cont-ci fastpli-ci

# copy repository to contaier
docker cp /tmp/fastpli-ci/. fastpli-cont-ci:/code/fastpli

# clean cp repository
rm -rf /tmp/fastpli-ci

# run
docker start -i fastpli-cont-ci
