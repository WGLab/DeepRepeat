FROM continuumio/miniconda3
USER root

WORKDIR /app

SHELL [ "/bin/bash", "--login", "-c" ]

COPY ./DeepRepeat/bin/ ./
COPY ./DeepRepeat/environment.yml .
#RUN ls /app/

RUN conda env create -f environment.yml && conda clean --yes --all
RUN echo "conda activate py36deeprepeat" > ~/.bashrc
ENV PATH=/opt/conda/envs/py36deeprepeat/bin:$PATH

#RUN ls /opt/conda/envs/py36deeprepeat/bin
#RUN x86_64-conda_cos6-linux-gnu-g++ --version
#RUN x86_64-conda_cos6-linux-gnu-c++ --version

ENV DR_conda_base="/opt/conda/envs/py36deeprepeat/"

WORKDIR /app/scripts
RUN h5c++ -O3 -std=c++11 -o IndexF5files ComFunction.c Fast5Index.c IndexF5files.c

RUN h5c++ -O3 -std=c++11 -I ${DR_conda_base}/include -L${DR_conda_base}/lib -lhts -o genomic1FE ComFunction.c ComOption.c BamReader.c Fast5Index.c Fast5Reader.c RepeatFeatExtract.c genomic1FE.c ${DR_conda_base}/lib/libhdf5_hl_cpp.a ${DR_conda_base}/lib/libhdf5_cpp.a ${DR_conda_base}/lib/libhdf5_hl.a ${DR_conda_base}/lib/libhdf5.a -lz -ldl -lpthread

RUN chmod +x /app/scripts/IndexF5files
RUN chmod +x /app/scripts/genomic1FE
RUN chmod +x /app/DeepRepeat.py

ENV PYTHONWARNINGS="ignore"
WORKDIR /app/
ENV PATH=/app/::$PATH

#ENTRYPOINT ["python /app/DeepRepeat.py", "--help"]
#CMD ["python /app/DeepRepeat.py", "--help"]

CMD ["python", "/app/DeepRepeat.py"]


