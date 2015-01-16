
How to install isochoreFinder
=============================

##dependencies

The following programs and libraries are require to compile and run python code

* [python 2.7](https://www.python.org/downloads/)
* [libgd](http://libgd.bitbucket.org/)
* [giflib](http://sourceforge.net/projects/giflib/)

And the following python libraries are required to run ISOfinder.py and ISOfamilies.py:

* [gdmodule](https://github.com/Solomoriah/gdmodule)
* [Pillow](http://python-pillow.github.io/)
* [matplotlib](http://matplotlib.org/)

The most easy way to install dependancies is via your linux distribution package manager. For example, to install all python dependencies on Debian/Ubuntu:

```
$ sudo apt-get install python-gd python-imaging python-matplotlib python-biopython
```

And this will also check and resolve any dependancies. In case you are not a system administrator, or you want to compile the last library versions, you have to install all development libraries dependancies to compile python libraries, and then install python packages. I suggest to manage different python version installations using [pyenv](https://github.com/yyuu/pyenv). To install libgd via package manager (for example in Debian/Ubuntu):

```
$ sudo apt-get install libgd2-xpm libgd2-xpm-dev libgif-dev
```

The recommend way to install manually python packages is by using [pip](http://dubroy.com/blog/so-you-want-to-install-a-python-package/). In this case:

```
$ pip install Pillow
$ pip install matplotlib
$ pip install gdmodule
$ pip install biopython
```

(or use sudo if you don't have required permissions)

##installing isochoreFinder using GIT

Once dependencies are satisfied, simply install isochoreFinder using git:

```
$ git clone https://github.com/bunop/isochoreFinder.git
```

To test isochoreFinder installation, enter the isochoreFinder directory and simply type:

```
$ ./ISOfinder.py --help
```

If the installation is ok, you will see the help message of ISOfinder.py:

