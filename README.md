
# alphaLoop 
Accurate perturbative computation of cross-sections with Numerical Loop-Tree Duality.

## 1. Installation and requirements

Install MadGraph and alphaLoop:

```sh
wget https://launchpad.net/mg5amcnlo/2.0/2.7.x/+download/MG5aMC_3.0.2.beta.py3.tgz
cd <MG_ROOT_PATH>/PLUGIN && git clone git@bitbucket.org:vahirschi/alphaloop.git
```

As many of our example generations use an internal model to alphaLoop, it is best to soft-link it in the models directory of MadGraph:

```
cd models
ln -s ../PLUGIN/alphaloop/models/aL_sm .
```

Then

```sh
cd alphaloop
sh deploy.sh
```

to install all dependencies.

Then

```sh
cd LTD
python3 ltd_commons.py
```
to create a default hyperparameters file.

Then
```sh
cd ..
cd rust_backend
cargo build --bin ltd 
```
to compile the rust backend. Optionally, `build --release`, for a slightly faster version.

Make sure to install the python dependencies
```
pip install scipy progressbar2 networkx python-igraph wheel
```
where `pip` needs to refer to the `pip` of at least Python 3.7. On clusters you may have to run

```
python3.7 -m pip install <PACKAGES>
```


## 2. Usage

You can run one of the example cards:

```sh
python3.7 ./bin/mg5_aMC --debug --mode=alphaloop PLUGIN/alphaloop/examples/epem_a_ddxg.aL
```

and integrate from the `rust_backend` folder with:

```sh
env MG_NUMERATOR_PATH=../../MG5_aMC_v2_7_2_py3/TEST_epem_a_ddxg/ cargo run --bin ltd -- --cross_section_set ../../MG5_aMC_v2_7_2_py3/TEST_epem_a_ddxg/FORM/Rust_inputs/all_MG_supergraphs.yaml --debug 0 -c 16
```

You may have to include the following to your library path:
```sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<ALPHALOOP>/libraries/fjcore:<ALPHALOOP>/libraries/ecos:<ALPHALOOP>/libraries/scs/out:<ALPHALOOP>/libraries/Cuba-4.2
```


## 3. Working plan:

- Output a process and list all corresponding super-graphs from it.
- Then build a library that can be linked to rust_backed in order to retrieve a fast implementation of the numerator of all contributing super graphs.
- Run with the Rust backend