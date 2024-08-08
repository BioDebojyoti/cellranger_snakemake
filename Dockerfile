# syntax=docker/dockerfile:1

# Use an intermediate build stage to handle dependencies
FROM --platform=$BUILDPLATFORM ubuntu:20.04 AS build

# Arguments for the build platform
ARG TARGETPLATFORM
ARG BUILDPLATFORM
RUN echo "I am running on $BUILDPLATFORM, building for $TARGETPLATFORM" > /log

# Ensure system packages are up to date and install dependencies
RUN apt-get update && apt-get install -y \
    curl \
    bzip2 \
    libz-dev \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    git \
    wget

# Download and install Cell Ranger
RUN wget -O cellranger-8.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-8.0.1.tar.gz?Expires=1722017713&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=KAFG5KZNZHPp96qn6iL1geCN3ytEOQeXJaKMlIu-X81gzv0YLUEkwrJu0gUMqcR10U6~f7QtznkLp9F9V4HHbpsZjCjPxB1Lex0PNGY8H6GAHR2a5-dYNqorhpNh068YPvgSFwa-C8Ve98Wv~jbuELdl4hLgzhUBQVDMoGxC~FRO0D0ljrGG9m4y9MDvefKIMvcKKLA~kfeHM0xnWIBKWYYhY9ObHqYMVoYmm0Op5BVtM97aah0l4XGXKFldpiDV-tqVCNEwY3IC56CrLCVAqM~s0fmNFWK0KAXafVkOVdMHDQXL~q0IbsETSE7JF2IEOcg3NbojuuIkRX2P0woFLw__"

RUN tar -xzvf cellranger-8.0.1.tar.gz && \
    mv cellranger-8.0.1 /opt/cellranger && \
    ln -s /opt/cellranger/cellranger /usr/local/bin/cellranger

# Clean up
RUN rm cellranger-8.0.1.tar.gz && apt-get clean && rm -rf /var/lib/apt/lists/*

# Final stage to create a minimal image
FROM ubuntu:20.04

# Copy Cell Ranger from the build stage
COPY --from=build /opt/cellranger /opt/cellranger
COPY --from=build /usr/local/bin/cellranger /usr/local/bin/cellranger
COPY --from=build /log /log

# Set the entrypoint to cellranger
ENTRYPOINT ["/bin/bash"]
