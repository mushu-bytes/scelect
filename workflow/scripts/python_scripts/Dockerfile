FROM continuumio/miniconda3
RUN conda install -c conda-forge scanpy python-igraph leidenalg
RUN pip install 'altair[all]'
RUN pip install scib
RUN conda install --solver=classic -c conda-forge python-annoy
RUN pip install scanorama
RUN pip install harmonypy
# env variable required for scanpy to run
ENV NUMBA_CACHE_DIR=tmp/
ENV MPLCONFIGDIR=tmp/
ENTRYPOINT ["python"]
