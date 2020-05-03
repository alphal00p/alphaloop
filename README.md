
# alphaLoop 
Accurate perturbative computation of cross-sections with Numerical Loop-Tree Duality.

## 1. Installation and requirements

`wget https://launchpad.net/mg5amcnlo/2.0/2.7.x/+download/MG5_aMC_v2.7.2.py3.tar.gz`

or

`bzr branch lp:~maddevelopers/mg5amcnlo/3.0.2.py3`

then

`cd <MG_ROOT_PATH>/PLUGIN && git clone git@bitbucket.org:vahirschi/alphaloop.git`

## 2. Usage

For now, only:

```
cd <MG_ROOT_PATH>
python3 ./bin/mg5_aMC --mode=alphaloop
alphaLoop > hello_world You
alphaLoop.Interface: Hello World, and welcome you, You!
alphaLoop > [enter_commands_here]
```

or run some example:

```
python3 ./bin/mg5_aMC PLUGIN/alphaloop/examples/epem_jjj_output.aL
```

## 3. Working plan:

Output a process and list all corresponding super-graphs from it.
Then build a library that can be linked to rust_backed in order to retrieve a fast implementation of the numerator of all contributing super graphs. As simple as that.