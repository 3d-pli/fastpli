#!/bin/bash
set -euo pipefail

docker build -t fastpli-ci - <./CI/dockerfile

# clone current commit state
rm -rf /tmp/fastpli-ci
git clone . /tmp/fastpli-ci

if docker container inspect fastpli-cont-ci &>/dev/null; then

   docker stop fastpli-cont-ci
   docker rm fastpli-cont-ci
fi

# # remove fastpli container
# if [ ! "$(docker ps -q -f name=fastpli-cont-ci)" ]; then
#     if [ "$(docker ps -aq -f status=exited -f name=fastpli-cont-ci)" ]; then
#         # cleanup
#         docker rm fastpli-cont-ci
#     fi
#    #  # run your container
#    #  docker run -d --name fastpli-cont-ci fastpli-ci
# fi

# docker stop fastpli-cont-ci || true && docker rm fastpli-cont-ci || true

docker create --name fastpli-cont-ci fastpli-ci
docker cp /tmp/fastpli-ci/. fastpli-cont-ci:/code/fastpli
rm -rf /tmp/fastpli-ci
docker start -i fastpli-cont-ci
