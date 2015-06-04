
How to install isochoreFinder
=============================

## Dependencies

The following programs and libraries are require to compile and run python code

* [python 2.7](https://www.python.org/downloads/)
* [libgd](http://libgd.bitbucket.org/)
* [giflib](http://sourceforge.net/projects/giflib/)

And the following python libraries are required to run isochoreFinder.py and isochoreFamilies.py:

* [gdmodule](https://github.com/Solomoriah/gdmodule)
* [Pillow](http://python-pillow.github.io/)
* [matplotlib](http://matplotlib.org/)

To install isochoreFinder dependencies you can follow one of the proposed methods (by package, by pip, ...). Choose the method more preferable for you

### Install dependencies via package manager

The most easy way to install dependancies is via your linux distribution package manager. For example, to install all python dependencies on Debian/Ubuntu:

```bash
$ sudo apt-get install python-gd python-imaging python-matplotlib python-biopython
```

And this will also check and resolve any dependancies.

### Installing dependencies using pip

In case you are not a system administrator, or you want to compile the last library versions, you have to install all development libraries dependancies to compile python libraries, and then install python packages. I suggest to manage different python version installations using [pyenv](https://github.com/yyuu/pyenv). Otherwise, you can install libraries locally using [virtualenv](https://virtualenv.pypa.io/en/latest/). In order to compile python libraries correctly, you need to have installed the needed build libraries. To install libraries dependancies via package manager (for example in Debian/Ubuntu):

```bash
$ sudo apt-get install libgd2-xpm-dev libgif-dev
```

The recommend way to install manually python packages is by using [pip](http://dubroy.com/blog/so-you-want-to-install-a-python-package/). In this case:

```bash
$ pip install Pillow
$ pip install matplotlib
$ pip install gdmodule
$ pip install biopython
```

(or use sudo if you don't have required permissions)

## Installing isochoreFinder using GIT

Once dependencies are satisfied, simply install isochoreFinder using git:

```bash
$ git clone https://github.com/bunop/isochoreFinder.git
```

To test isochoreFinder installation, enter the isochoreFinder directory and simply type:

```bash
$ ./isochoreFinder.py --help
```

If the installation is ok, you will see the help message of isochoreFinder.py by invoking it:

```
usage: isochoreFinder.py [-h] -i INFILE [-o OUTFILE] [-g GRAPHFILE] [-w WINDOWFILE]
                    [-v] [--windowgraph WINDOWGRAPH] [--draw_legend]
                    [--force_overwrite] [--sequence_start SEQUENCE_START]
                    [--max_length MAX_LENGTH] [--draw_chname DRAW_CHNAME]
                    [--window_size WINDOW_SIZE] [--y_max Y_MAX]
                    [--y_min Y_MIN]

Find Isochores in sequences

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Input Fasta File (even compressed)
  -o OUTFILE, --outfile OUTFILE
                        Output isochores CSV file
  -g GRAPHFILE, --graphfile GRAPHFILE
                        Output graph filename (PNG)
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
```

## Installing isochoreFinder using docker

There are instructions on how to build and run a isochoreFinder image in the Docker subdirectory in isochoreFinder project. Please look at README.md inside that directory.
