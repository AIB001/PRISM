# Build from the PRISM repository root:
#   docker build -f prism/generation/environments/targetdiff/Dockerfile \
#     -t local/prism-targetdiff:142f1eb .
#
# For a reproducible production build, override BASE_IMAGE with a digest-pinned
# equivalent and record the resulting image digest in generation config.
ARG BASE_IMAGE=continuumio/miniconda3:23.5.2-0
FROM ${BASE_IMAGE}

ARG TARGETDIFF_COMMIT=142f1eb7178480d435fe0b8cb95a99beb48997c7

RUN apt-get update \
    && apt-get install --yes --no-install-recommends git ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN git init /opt/targetdiff \
    && git -C /opt/targetdiff remote add origin https://github.com/guanjq/targetdiff.git \
    && git -C /opt/targetdiff fetch --depth 1 origin "${TARGETDIFF_COMMIT}" \
    && git -C /opt/targetdiff checkout --detach FETCH_HEAD \
    && test "$(git -C /opt/targetdiff rev-parse HEAD)" = "${TARGETDIFF_COMMIT}"

RUN conda env create --file /opt/targetdiff/environment.yaml \
    && conda clean --all --yes

COPY prism/generation/wrappers/targetdiff.py /opt/prism_wrappers/targetdiff.py

WORKDIR /opt/targetdiff
ENTRYPOINT ["conda", "run", "--no-capture-output", "--name", "targetdiff"]
