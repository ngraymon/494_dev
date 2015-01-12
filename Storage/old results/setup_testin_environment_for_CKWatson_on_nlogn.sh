#!/bin/bash

# location of the miniconda installer and the installation location of miniconda
INSTALL_DIRECTORY=~/CKWatson_testing_folder/python_directory
INSTALL_FILE_NAME=~/miniconda_src.sh

# make sure we haven't already installed conda
if [ -d "$INSTALL_DIRECTORY" ]; then
	echo "Conda has already been installed"
fi

# only download the installer if it isn't present already
if [ -e "$INSTALL_FILE_NAME" ]; then
	echo "The Installation file has already been downloaded\nProceeding with installation without redownloading the installer"
else
	case "$OSTYPE" in
		solaris*)
			echo "SOLARIS" ;;
		darwin*)
			curl --progress-bar --anyauth --output "$INSTALL_FILE_NAME" http://repo.continuum.io/miniconda/Miniconda-3.6.0-MacOSX-x86_64.sh ;;
		linux*)
			curl --progress-bar --anyauth --output "$INSTALL_FILE_NAME" http://repo.continuum.io/miniconda/Miniconda-3.6.0-Linux-x86_64.sh ;;
		*)	
			echo "unknown: $OSTYPE, error" ;;
	esac
fi

# install miniconda
chmod +x "$INSTALL_FILE_NAME"
"$INSTALL_FILE_NAME" -b -p "$INSTALL_DIRECTORY" 

# install binstar so we can download the packages for MMTK
"$INSTALL_DIRECTORY"/bin/conda install binstar --yes

# inform conda where the packages are located
"$INSTALL_DIRECTORY"/bin/conda config --add channels http://conda.binstar.org/ngraymon

# make an mmtk directory
"$INSTALL_DIRECTORY"/bin/conda create --yes --name marcel python=3.4 numpy=1.9 scipy=0.14 matplotlib=1.4

# install the packages
#"$INSTALL_DIRECTORY"/bin/conda install --yes --name mmpigs mmtk_pigs

# only remove the install file if we successfully created the miniconda directory
# otherwise there might be an error to address and we don't want to redownload the file multiple times
if [ -d "$INSTALL_DIRECTORY" ]; then
	rm -f "$INSTALL_FILE_NAME"
fi
