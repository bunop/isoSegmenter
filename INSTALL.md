
How to install isoSegmenter
=============================

isoSegmenter is written as a python module and need to be installed in a Unix distribution. No test were made in windows environment.

## Dependencies

The following programs and libraries are require to compile and run python code

* [python 2.7](https://www.python.org/downloads/)
* [libgd](http://libgd.bitbucket.org/)
* [giflib](http://sourceforge.net/projects/giflib/)

And the following python libraries are required to run isoSegmenter.py and isochoreFamilies.py:

* [gdmodule](https://github.com/Solomoriah/gdmodule)
* [Pillow](http://python-pillow.github.io/)
* [matplotlib](http://matplotlib.org/)
* [biopython](http://biopython.org/wiki/Main_Page)

To install isoSegmenter dependencies, you may use the [system package manager](https://github.com/bunop/isoSegmenter/blob/master/INSTALL.md#install-dependencies-via-package-manager) or [building packages from sources](https://github.com/bunop/isoSegmenter/blob/master/INSTALL.md#building-packages-from-sources). Choose the method more preferable for you

### Install dependencies via package manager

The most easy way to install dependencies is via your Linux distribution package manager. For example, to install all python dependencies on Debian/Ubuntu:

```bash
$ sudo apt-get install python-gd python-imaging python-matplotlib python-biopython
```

And this will also check and resolve any dependencies. Then go to [Installing isoSegmenter using GIT](https://github.com/bunop/isoSegmenter/blob/master/INSTALL.md#installing-isosegmenter-using-git)

### Building packages from sources

In case you are not a system administrator, or you want to compile the last library versions, you have to install all development libraries dependencies to compile python libraries, and then install python packages. I suggest to install libraries locally using [virtualenv](https://virtualenv.pypa.io/en/latest/). Otherwise, you can manage different python versions and installations using [pyenv](https://github.com/yyuu/pyenv). In order to compile python libraries correctly, you need to have installed the needed build libraries. To install libraries dependancies via package manager (for example in Debian/Ubuntu):

```bash
$ sudo apt-get install libgd2-xpm-dev libgif-dev
```

#### Using pyenv (optional)

[pyenv](https://github.com/yyuu/pyenv) lets you easily switch between multiple versions of Python. It could be installed as a non root user and give the opportunity to deal with different versions of python. This solution is not necessary to install isoSegmenter, however we suggest this solution if you don't want to install python dependencies in the system python environment as a root, or your system has a very old version of python. isoSegmenter was initially wrote on Python 2.6 and then tested on Python 2.7.  

#### Using virtualenv (recommeded!)

A Virtual Environment is a tool to keep the dependencies required by different projects in separate places, by creating virtual Python environments for them. [virtualenv](http://docs.python-guide.org/en/latest/dev/virtualenvs/) creates a folder which contains all the necessary executables to use the packages that a Python project would need. You can install both virtualenv and pyenv, virtualenv will create a directory in which all libraries of your pyenv python will be placed. Starting a virtualenv
environments is simple. You need to to install the `virtualenv` python package:

```bash
$ pip install virtualenv
```

[pip](http://pip.readthedocs.org/en/stable/) is the PyPA recommended tool for installing Python packages. Since Python 2.7.9 is included in python distribution; If you don't have pip installed, you could install it by [following this guide](http://pip.readthedocs.org/en/stable/installing/). You may need to install virtualenv (and pip, if needed) as a root user, if you want to install virtualenv inside system python environments. Once virtualenv is installed, we can create a python virtual environment and install packages without root privileges. Choose a directory in which the virtualenv environment will be placed, then set up the python environment:

```bash
$ virtualenv --verbose env
$ source env/bin/activate
```

The first command will create the `env` directory in which all the dependencies will be placed. Next we can start the python environment by sourcing the `activate` script that is placed inside the `env` directory, used in this example. If things gone correctly, the name of the current virtual environment will now appear on the left of the prompt (e.g. `(env)Your-Computer:your_project UserName$`) to let you know that it’s active. Now each python or pip invocation will use the python libraries and directories inside the `env` directory. You can inspect which python is currently used by typing:

```bash
$ which python
```

You will see that the current python installation is inside the `env` directory. You can follow the section [Installing dependencies using pip](https://github.com/bunop/isoSegmenter/blob/master/INSTALL.md#installing-dependencies-using-pip) and install all python libraries inside virtualenv environment. To exit from python environment, you can close the terminal or call the `deactivate` scripts:

```bash
$ deactivate
```

Once deactivated, your prompt has no more the `env` directory on the left and your python is the default python interpreter. To resume your isoSegmenter installation, you will need only to `activate` the python virtual environment with the `source env/bin/activate` command

### Installing dependencies using pip

The recommend way to install manually python packages is by using [pip](http://dubroy.com/blog/so-you-want-to-install-a-python-package/). In this case:

```bash
$ pip install Pillow
$ pip install matplotlib
$ pip install gdmodule
$ pip install biopython
```

(or use sudo if you don't have required permissions)

## Installing isoSegmenter using GIT

Once dependencies are satisfied, simply install isoSegmenter by cloning project with Git and using pip inside isoSegmenter directory:

```bash
$ git clone https://github.com/bunop/isoSegmenter.git
$ cd isoSegmenter
$ pip install .
```

To test isoSegmenter installation, simply type:

```bash
$ isoSegmenter.py --help
```

If the installation went ok, you will see the help message of isoSegmenter.py by invoking it:

```
usage: isoSegmenter.py [-h] -i INFILE [-o OUTFILE] [-g GRAPHFILE] [-b BARFILE]
                       [-w WINDOWFILE] [-v] [--windowgraph WINDOWGRAPH]
                       [--draw_legend] [--force_overwrite]
                       [--sequence_start SEQUENCE_START]
                       [--max_length MAX_LENGTH] [--draw_chname DRAW_CHNAME]
                       [--window_size WINDOW_SIZE] [--y_max Y_MAX]
                       [--y_min Y_MIN] [--isochore_min_size ISOCHORE_MIN_SIZE]

Find Isochores in sequences

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Input Fasta File (even compressed)
  -o OUTFILE, --outfile OUTFILE
                        Output isochores CSV file
  -g GRAPHFILE, --graphfile GRAPHFILE
                        Output graph filename (PNG)
  -b BARFILE, --barfile BARFILE
                        Output bar graph filename (PNG)
  -w WINDOWFILE, --windowfile WINDOWFILE
                        Output windows CSV file
  -v, --verbose         Verbosity level
  --windowgraph WINDOWGRAPH
                        Output windows Graph file
  --draw_legend         Draw legend on the right side of the image
  --force_overwrite     Force overwrite
  --sequence_start SEQUENCE_START
                        Start segmentation from this position (1-based
                        coordinates)
  --max_length MAX_LENGTH
                        Scan for isochores until for this dimension in bp
  --draw_chname DRAW_CHNAME
                        Draw chromosome name in figure
  --window_size WINDOW_SIZE
                        Set window size in bp (default: '100000')
  --y_max Y_MAX         Set max value in graph (default: '65')
  --y_min Y_MIN         Set min value in graph (default: '30')
  --isochore_min_size ISOCHORE_MIN_SIZE
                        Set how many windows an isochore need to have
                        (default: '2')

If you use isoSegmenter in your work, please cite this manuscript:

    Cozzi P, Milanesi L, Bernardi G. Segmenting the Human Genome into Isochores.
    Evolutionary Bioinformatics. 2015;11:253-261. doi:10.4137/EBO.S27693

```

## Installing isoSegmenter using docker

There are instructions on how to build and run a isoSegmenter image in the Docker subdirectory in isoSegmenter project. Please look at [README.md](https://github.com/bunop/isoSegmenter/blob/master/Docker/README.md#installing-isosegmenter-using-docker) inside `Docker` directory, or at the [isoSegmenter page on Docker Hub](https://hub.docker.com/r/bunop/isosegmenter/)
