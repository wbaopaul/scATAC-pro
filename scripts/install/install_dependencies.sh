#!/bin/bash

NORMAL="\\033[0;39m"
RED="\\033[0;31m"
BLUE="\\033[0;34m"
SOFT="scATAC-pro"

die() {
    echo -e "$RED""Exit - ""$*""$NORMAL" 1>&2
    exit 1
}

function usage {
    echo -e "Usage : ./install_dependencies.sh"
    echo -e "-o"" <installation folder>"
    echo -e "-q"" <quiet>"
    echo -e "-h"" <help>"
    exit;
}


echo -e "$RED""Make sure internet connection works for your shell prompt under current user's privilege ...""$NORMAL";
echo -e "$BLUE""Starting $SOFT installation ...""$NORMAL";


############################################################
## Initialize
############################################################
quiet=0
set -- $(getopt o:qh "$@")
while [ $# -gt 0 ]
do
    case "$1" in
        (-o) install_dir=$2;shift;;
        (-q) quite=1;shift;;
        (-h) usage;;
        (--) shift;break;;
        (-*) echo "$0: error - unrecongnized option $1" 1>&2; exit 1;;
        (*) break;;
    esac
    shift
done

if [[ -z $install_dir ]]
then
    die "Error: Installation directory not defined (-o)"
fi


mkdir -p $install_dir



#############################################################################
## Search standard tools 
#############################################################################


## function to comapre version
## 0 =
## 1 >
## 2 <
vercomp () {
    if [[ $1 == $2 ]]
    then
        return 0
    fi
    local IFS=.
    local i ver1=($1) ver2=($2)
    # fill empty fields in ver1 with zeros
    for ((i=${#ver1[@]}; i<${#ver2[@]}; i++))
    do
        ver1[i]=0
    done

    for ((i=0; i<${#ver1[@]}; i++))
    do
        if [[ -z ${ver2[i]} ]]
        then
            # fill empty fields in ver2 with zeros
            ver2[i]=0
        fi
        if ((10#${ver1[i]} > 10#${ver2[i]}))
        then
            return 1
        fi
        if ((10#${ver1[i]} < 10#${ver2[i]}))
        then
            return 2
        fi
    done
    return 0
}


#check for make
which make > /dev/null 2>&1
if [ $? != "0" ]; then
    echo -e "$RED""Can not proceed without make, please install and re-run (Mac users see: http://developer.apple.com/technologies/xcode.html)""$NORMAL"
    exit 1;
fi

#check for g++
which g++ > /dev/null 2>&1
if [ $? != "0" ]; then
echo -e "$RED""Can not proceed without g++, please install and re-run""$NORMAL"
exit 1;
fi

# python
which python > /dev/null 2>&1
if [ $? != "0" ]; then
    echo -e "$RED""Can not proceed without Python, please install and re-run""$NORMAL"
    exit 1;
else
    pver=`python --version 2>&1 | cut -d" " -f2`
    vercomp $pver "3.6.0"
    if [[ $? == 2 ]]; then
        echo -e "$RED""Python v3.6.0 or higher is needed [$pver detected].""$NORMAL"
        exit 1;
    fi
fi


# R
which R > /dev/null 2>&1
if [ $? != "0" ]; then
    echo -e "$RED""Can not proceed without R, please install and re-run""$NORMAL"
    exit 1;
else
    pver=`R --version 2>&1 | head -1 | cut -d" " -f3`
    vercomp $pver "3.5.3"
    if [[ $? == 2 ]]; then
    echo -e "$RED""R v3.5.3 or higher is needed [$pver detected].""$NORMAL"
    exit 1;
    fi
fi


# perl
which perl > /dev/null 2>&1
if [ $? != "0" ]; then
    echo -e "$RED""Can not proceed without perl, please install and re-run""$NORMAL"
    exit 1;
fi



#check OS (Unix/Linux or Mac)
os=`uname`;

# get the right download program
if [ "$os" = "Darwin" ]; then
    # use curl as the download program 
    get="curl -L -o"
else
    # use wget as the download program
    get="wget --no-check-certificate -O"
fi

if [ -d ./tmp ]; then
    rm -r ./tmp
fi
mkdir ./tmp
cd ./tmp


################ Install dependencies  ###################

## By default, dependencies will be installed in the same path with scATAC-pro
PREFIX_BIN=${install_dir}

if [ ! -w $PREFIX_BIN ]; then
   quiet=1
fi

if [[ $quiet == 0 ]]; then
    #echo "Where should missing software prerequisites be installed ? [$PREFIX_BIN] "
    #read ans
    #ans=${ans:-$PREFIX_BIN}
    #PREFIX_BIN=$ans
    echo "Dependencies will be install in $PREFIX_BIN"
fi

if [ ! -d $PREFIX_BIN ]; then
    echo "Directory $PREFIX_BIN does not exist!"

    if [[ $quiet == 0 ]]; then
        echo -n "Do you want to create $PREFIX_BIN folder ? (y/n) [n] : "
        read ans
        if [ XX${ans} = XXy ]; then
            mkdir -p $PREFIX_BIN || die "Cannot create  $PREFIX_BIN folder. Maybe missing super-user (root) permissions"
        else
            die "Must specify a directory to install required softwares!"
        fi
    else
        die "Error - unable to install/check dependancies !"
    fi
fi

if [ ! -w $PREFIX_BIN ]; then
    die "Cannot write to directory $PREFIX_BIN. Maybe missing super-user (root) permissions to write there. Please run 'make configure prefix=YOUR_INSTALL_PATH' ";
fi 

##########################
##macs2
##########################
wasInstalled=0;
which macs2 > /dev/null 2>&1
if [ $? != "0" ]; then
    echo -e "$RED""macs2 not installed, trying to install it through anaconda...""$NORMAL"
    which conda > /dev/null 2>&1
    if [ $? != "0" ]; then
        echo "Install anaconda:"
        if [ "$os" = "Darwin" ]; then 
            $get tmp.sh https://repo.anaconda.com/archive/Anaconda3-2019.03-MacOSX-x86_64.sh
        else
            $get tmp.sh https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-ppc64le.sh
        fi
        bash tmp.sh
    fi
    
    conda_path=$(dirname `which conda`)
    if [ $? != '0'  ]; then 
        echo "anaconda not install successfully, please install anaconda or macs2 manually!" 
    fi
    unset PYTHONPATH 
    unset PYTHONHOME 
    conda create --name py2 python=2.7
    conda activate py2
    pip install --user --upgrade pip
    pip install --user --upgrade Numpy 
    pip install --user macs2 
    conda deactivate 
fi

macs2ver=`macs2 --version 2>&1 | cut -d " " -f 2`
if [ $? != '0' ]; then    
    which conda > /dev/null 2>&1
    if [ $? != "0" ]; then
        echo "Install anaconda:"
        if [ "$os" = "Darwin" ]; then 
            $get tmp.sh https://repo.anaconda.com/archive/Anaconda3-2019.03-MacOSX-x86_64.sh
        else
            $get tmp.sh https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-ppc64le.sh
        fi
        bash tmp.sh
    fi
    
    conda_path=$(dirname `which conda`)
    if [ $? != '0'  ]; then 
        echo "anaconda not install successfully, please install anaconda or macs2 manually!" 
    fi
    conda create --name py2 python=2.7
    source activate py2
    unset PYTHONPATH 
    pip install --user --upgrade pip
    pip install --user --upgrade Numpy 
    pip install --user macs2 
    conda deactivate 
fi

macs2ver=`macs2 --version 2>&1 | cut -d " " -f 2`
if [ $? != '0' ]; then
    echo -e "$RED"" macs2 Not installed successfully ].""NORMAL"
    exit 
fi

vercomp $macs2ver "2.1.1"
if [[ $? == 2 ]]; then
    echo -e "$RED""macs2 v2.1.1 or higher is needed [$macs2pver detected].""$NORMAL"
    exit 1;
else
    echo -e "$BLUE""macs2 appears to be installed successfully""$NORMAL"
fi


##################################################################
## samtools  
##################################################################

wasInstalled=0;
which samtools > /dev/null 2>&1
if [ $? = "0" ]; then
        
    samver=`samtools --version | grep samtools | cut -d" " -f2`
    vercomp $samver "1.9"
    if [[ $? == 2 ]]; then
        echo -e "$RED""samtools v1.9 or higher is needed [$samver detected].""$NORMAL"
        echo -e "$RED""I will try to install v1.9 now ...""$NORMAL"
    else
        echo -e "$BLUE""Samtools appears to be already installed. ""$NORMAL"
        wasInstalled=1;
    fi
fi

if [ $wasInstalled == 0 ]; then
    echo "Installing samtools ..."
    #From sources
    $get samtools-1.9.tar.bz2  http://sourceforge.net/projects/samtools/files/samtools/1.9/samtools-1.9.tar.bz2/download
    tar -xvjpf samtools-1.9.tar.bz2
    cd samtools-1.9
    make
    cd ..
    cp -r samtools-1.9 $PREFIX_BIN/
    export PATH=$PREFIX_BIN/samtools-1.9/:$PATH

    check=`samtools view -h 2>&1 | grep -i options`;
    if [ $? = "0" ]; then
        echo -e "$BLUE""samtools appears to be installed successfully""$NORMAL"
    else
        echo -e "$RED""samtools NOT installed successfully""$NORMAL"; exit 1;
    fi
fi


#########################################################################################
## bedtools 
#########################################################################################

wasInstalled=0;
which bedtools > /dev/null 2>&1
if [ $? = "0" ]; then

    bedver=`bedtools --version | cut -d" " -f2 | cut -d"-" -f1`
    vercomp $bedver "1.0"
    if [[ $? == 2 ]]; then
        echo -e "$RED""bedtools v2.27.1 or higher is needed [$bedver detected].""$NORMAL"
        echo -e "$RED""I will try to install v2.27.1 now ...""$NORMAL"
    else    
        echo -e "$BLUE""bedtools appears to be already installed. ""$NORMAL"
        wasInstalled=1;
    fi
fi

if [ $wasInstalled == 0 ]; then
    echo "Installing bedtools ..."
    #From sources
    $get bedtools-2.27.1.tar.gz https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz 
    tar -xzvf bedtools-2.27.1.tar.gz
    cd bedtools2
    make
    cd ..
    cp -r bedtools2 $PREFIX_BIN/
    export PATH=$PREFIX_BIN/bedtools2/bin:$PATH

    check=`bedtools | grep -i options`;
    if [ $? = "0" ]; then
    echo -e "$BLUE""bedtools appears to be installed successfully""$NORMAL"
    else
    echo -e "$RED""bedtools NOT installed successfully""$NORMAL"; exit 1;
    fi
fi


#########################################################################################
## deeptools 
#########################################################################################

wasInstalled=0;
which deeptools > /dev/null
if [ $? = "0" ]; then

    dver=`deeptools --version 2>&1 | cut -d" " -f2 `
    vercomp $dver "3.2.1"
    if [[ $? == 2 ]]; then
        echo -e "$RED""deeptools v3.2.1 or higher is needed [$dver detected].""$NORMAL"
        echo -e "$RED""I will try to install v3.2.1 now ...""$NORMAL"
    else
        echo -e "$BLUE""deeptools appears to be already installed. ""$NORMAL"
        wasInstalled=1;
    fi
fi

if [ $wasInstalled == 0 ]; then
    echo "Installing deeptools ..."
    ##From sources
    #$get deeptools-3.2.1.tar.gz https://github.com/deeptools/deepTools/archive/3.2.1.tar.gz 
    #tar -xzvf deeptools-3.2.1.tar.gz
    #cd deepTools-3.2.1
    #python setup.py install --prefix $PREFIX_BIN/deepTools3
    #export PATH=$PREFIX_BIN/deepTools3/bin:$PATH
    
    ## from pip
    pip install --user --upgrade deeptools

    check=`deeptools --version `;
    if [ $? = "0" ]; then
    echo -e "$BLUE""deeptools appears to be installed successfully""$NORMAL"
    else
    echo -e "$RED""deeptools NOT installed successfully""$NORMAL"; exit 1;
    fi
fi



##################################################################################
## install HINT through rgt
##################################################################################
wasInstalled=0
which rgt-hint > /dev/null

if [ $? = "0" ]; then

    hver=`rgt-hint --version 2>&1 | cut -d"-" -f3 | cut -d 'v' -f2`
    vercomp $hver "0.12.1"
    if [[ $? == 2 ]]; then
        echo -e "$RED""rgt-hint v0.12.1 or higher is needed [$dver detected].""$NORMAL"
        echo -e "$RED""I will try to install v0.12.1 now ...""$NORMAL"
    else
        echo -e "$BLUE""rgt-hint appears to be already installed. ""$NORMAL"
        wasInstalled=1;
    fi
fi


if [ $wasInstalled == 0 ]; then
    echo "Installing rgt-hint through installing RGT suite..."
    
    echo "Since it needs python2 to install RGT, I will try to use anaconda ..."
    which conda > /dev/null 2>&1
    if [ $? != "0" ]; then
        echo "I will install anaconda now ..."
        if [ "$os" = "Darwin" ]; then 
            $get tmp.sh https://repo.anaconda.com/archive/Anaconda3-2019.03-MacOSX-x86_64.sh
        else
            $get tmp.sh https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-ppc64le.sh
        fi
        bash tmp.sh
    fi
    
    if [ $? != '0'  ]; then 
        echo -e "$RED""anaconda not install successfully, thus I cannot install RGT" "$NORMAL"; 
        echo -e "$RED""If you want to run footprinting analysis, please install RGT by yourself!""$NORMAL";
    else
        conda create --name py2 python=2.7
        source activate py2
        unset PYTHONPATH 
        ## from pip
        pip install --user cython scipy numpy
        pip install --user RGT

        conda deactivate 
        check=`rgt-hint --version `;
        if [ $? = "0" ]; then
            echo -e "$BLUE""rgt-hint appears to be installed successfully""$NORMAL"
        else
            echo -e "$RED""rgt-hint NOT installed successfully""$NORMAL"; 
            echo -e "$RED""If you want to run footprinting analysis, please install RGT by yourself!""$NORMAL";
        fi
    fi
fi


echo "Installaing GEM as alternative peak calling algorithm"
$get gem.v3.4.tar.gz http://groups.csail.mit.edu/cgs/gem/download/gem.v3.4.tar.gz
tar xzvf gem.v3.4.tar.gz
cp -r gem ${PREFIX_BIN}/
GEM_PATH=${PREFIX_BIN}/gem

###################################################################################
## Install R packages 
###################################################################################
#Install R Packages
wasInstalled=0
if [ $wasInstalled == 0 ]; then
    echo "Installing missing R packages ..."
    #R CMD BATCH ../scripts/install/install_Rpackages.R install_Rpackages.Rout
    R --vanilla < ../scripts/install/install_Rpackages.R

    if [ $? == "0" ]; then
        echo -e "$BLUE""R packages appear to be installed successfully""$NORMAL"
    else
        echo -e "$RED""R packages NOT installed successfully. You may manually install some packages which I cannot installed.""$NORMAL"; exit 1;
    fi
fi


##################################################################################
## install python packages
##################################################################################
pip install --user umap-learn


##################################################################################
## Trimmer: trimgalore  and Trimmomatic
##################################################################################
echo "install Trimmomatic"
$get Trimmomatic-0.39.zip http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
mv Trimmomatic-0.39 $PREFIX_BIN/
TRIMMOMATIC_PATH=$PREFIX_BIN/Trimmomatic-0.39


wasInstalled=0;
which trim_galore > /dev/null 2>&1
if [ $? = "0" ]; then
    echo -e "$BLUE""trimgalore appears to be already installed. ""$NORMAL"
    wasInstalled=1;
fi


if [ $wasInstalled == 0 ]; then
    echo "trim_galore was not detected!"
    echo "Installing trimgalore as default trimmer ..."
    #From sources
    $get TrimGalore-0.6.3.tar.gz https://github.com/FelixKrueger/TrimGalore/archive/0.6.3.tar.gz 
    tar -xzvf TrimGalore-0.6.3.tar.gz
    cp -r TrimGalore-0.6.3 $PREFIX_BIN/
    export PATH=$PREFIX_BIN/TrimGalore-0.6.3/:$PATH
    wasInstalled=0;
    which trim_galore > /dev/null 2>&1
    if [ $? = "0" ]; then
        echo -e "$BLUE""trim_galore appears to be installed successfully""$NORMAL"
    else
        echo -e "$RED""trim_galore NOT installed successfully""$NORMAL"; exit 1;
    fi
fi

###################################################################################
## check aligner, install bwa as default aligner
###################################################################################
wasInstalled=0;
which bwa > /dev/null 2>&1
if [ $? = "0" ]; then
    echo -e "$BLUE""bwa appears to be already installed. ""$NORMAL"
fi

which bowtie2 > /dev/null 2>&1
if [ $? = "0" ]; then
    echo -e "$BLUE""bowtie2 detected. ""$NORMAL"
fi

which bowtie > /dev/null 2>&1
if [ $? = "0" ]; then
    echo -e "$BLUE""bowtie detected. ""$NORMAL"
fi

if [ $wasInstalled == 0 ]; then
    echo "Installing bwa as default aligner ..."
    #From sources
    $get bwa-0.7.17.tar.bz2 https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2/download
    tar -xvjpf bwa-0.7.17.tar.bz2
    cd bwa-0.7.17
    make
    cd ..
    mv bwa-0.7.17 $PREFIX_BIN
    export PATH=$PREFIX_BIN/bwa-0.7.17/:$PATH

    which bwa > /dev/null 2>&1
    if [ $? = "0" ]; then
        echo -e "$BLUE""bwa appears to be installed successfully""$NORMAL"
    else
        echo -e "$RED""bwa NOT installed successfully""$NORMAL"; exit 1;
    fi
fi

wasInstalled=0

## Clean up
cd ..
rm -rf ./tmp

echo -e "$RED""Dependencies checked !""$NORMAL"



#############################################################################
## Create the configure_system file 
#############################################################################

echo "Create the configure_system.txt file ..."

CUR_DIR=`pwd`
echo -e "$BLUE""Check $SOFT configuration ... ""$NORMAL"

echo "#######################################################################" > configure_system.txt
echo "## $SOFT - System settings" >> configure_system.txt
echo "#######################################################################" >> configure_system.txt

echo "#######################################################################" >> configure_system.txt
echo "## Required Software - Specified the DIRECTORY name of the executables" >> configure_system.txt
echo "## If not specified, the program will try to locate the executable" >> configure_system.txt
echo "## using the 'which' command" >> configure_system.txt
echo "#######################################################################" >> configure_system.txt


echo "INSTALL_PREFIX = "${install_dir} >> configure_system.txt


which python > /dev/null 2>&1
if [ $? = "0" ]; then
    echo "PYTHON_PATH = "`dirname $(which python)` >> configure_system.txt
else
    die "PYTHON_PATH not found. Exit."
fi

which R > /dev/null 2>&1
if [ $? = "0" ]; then
    echo "R_PATH = "`dirname $(which R)` >> configure_system.txt
else
    die "R_PATH not found. Exit." 
fi

which samtools > /dev/null 2>&1
if [ $? = "0" ]; then
    echo "SAMTOOLS_PATH = "`dirname $(which samtools)`  >> configure_system.txt
else
    die "SAMTOOLS_PATH not found. Exit." 
fi

which deeptools > /dev/null 2>&1
if [ $? = "0" ]; then
    echo "DEEPTOOLS_PATH = "`dirname $(which deeptools)` >> configure_system.txt
else
    die "DEEPTOOLS_PATH not found. Exit." 
fi


which bedtools > /dev/null 2>&1
if [ $? = "0" ]; then
    echo "BEDTOOLS_PATH = "`dirname $(which bedtools)` >> configure_system.txt
else
    die "BEDTOOLS_PATH not found. Exit." 
fi


which bwa > /dev/null 2>&1
if [ $? = "0" ]; then
    echo "BWA_PATH = "`dirname $(which bwa)` >> configure_system.txt
fi


which bowtie2 > /dev/null 2>&1
if [ $? = "0" ]; then
    echo "BOWTIE2_PATH = "`dirname $(which bowtie2)` >> configure_system.txt
fi

which bowtie > /dev/null 2>&1
if [ $? = "0" ]; then
    echo "BOWTIE_PATH = "`dirname $(which bowtie)` >> configure_system.txt
fi

which macs2 > /dev/null 2>&1
if [ $? = "0" ]; then
    echo "MACS2_PATH = "`dirname $(which macs2)` >> configure_system.txt
fi

which perl > /dev/null 2>&1
if [ $? = "0" ]; then
    echo "PERL_PATH = "`dirname $(which perl)` >> configure_system.txt
fi

which trim_galore > /dev/null 2>&1
if [ $? = "0" ]; then
    echo "TRIM_GALORE_PATH = "`dirname $(which trim_galore)` >> configure_system.txt
fi


which rgt-hint > /dev/null 2>&1
if [ $? = "0" ]; then
    echo "HINT_PATH = "`dirname $(which rgt-hint)` >> configure_system.txt
fi

echo "TRIMMOMATIC_PATH = " $TRIMMOMATIC_PATH >> configure_system.txt
echo "GEM_PATH = " $GEM_PATH >> configure_system.txt

## check rights in PREFIX folder

if [ ! -w $install_dir ]; then
    die "Cannot install scATAC-pro in $install_dir directory. Maybe missing super-user (root) permissions to write there. Please specify another directory using 'make configure prefix=YOUR_INSTALL_PATH' ";
fi 

echo -e "$RED""All Done: run 'make install' now!""$NORMAL" ;


