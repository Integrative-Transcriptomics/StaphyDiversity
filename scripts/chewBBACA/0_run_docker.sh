# Create docker terminal with mounted volume
# Docker page: https://hub.docker.com/r/ummidock/chewbbaca
docker run --rm -it -v /ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis:/data ummidock/chewbbaca bash