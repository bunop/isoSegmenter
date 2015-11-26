
How to install isoSegmenter
=============================

## Using virtualenv (recommeded!)

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

To install isoSegmenter dependencies you can follow one of the proposed methods (by package, by pip, ...). Choose the method more preferable for you

### Install dependencies via package manager

The most easy way to install dependencies is via your linux distribution package manager. For example, to install all python dependencies on Debian/Ubuntu:

```bash
$ sudo apt-get install python-gd python-imaging python-matplotlib python-biopython
```

And this will also check and resolve any dependancies. Then go to [Installing isoSegmenter using GIT](#install-isoSegmenter)

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

## Installing isoSegmenter using GIT <a id="install-isoSegmenter"></a>

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

If you use isoSegmenter in your work, please cite this manuscripts:

    Cozzi P, Milanesi L, Bernardi G. Segmenting the Human Genome into Isochores.
    Evolutionary Bioinformatics. 2015;11:253-261. doi:10.4137/EBO.S27693
        
```

## Installing isoSegmenter using docker

There are instructions on how to build and run a isoSegmenter image in the Docker subdirectory in isoSegmenter project. Please look at README.md inside that directory.
