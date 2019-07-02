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
    echo -e "-c"" <configuration install file>"
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
set -- $(getopt c:o:qh "$@")
while [ $# -gt 0 ]
do
    case "$1" in
        (-c) conf=$2;shift;;
        (-o) install_dir=$2;shift;;
        (-q) quite=1;shift;;
        (-h) usage;;
        (--) shift;break;;
        (-*) echo "$0: error - unrecongnized option $1" 1>&2; exit 1;;
        (*) break;;
    esac;
    shift
done

if [[ -z $install_dir ]]
then
    die "Error: Installation directory not defined (-o)"
fi


if [[ ! -e $conf ]]
then
    die "Error: Configuration file not found"
fi




#############################################################################
## Read the config file
#############################################################################
while read curline_read; do
    curline=${curline_read// /}
    if [[ $curline != \#* && ! -z $curline ]]; then
        var=`echo $curline | awk -F= '{print $1}'`
        val=`echo $curline | awk -F= '{print $2}'`

        if [[ $var =~ "_PATH" ]]
        then
            if [[ ! -z $val ]]; then
            echo "export $val in PATH"
            export PATH=$val:$PATH
            fi
        else
            export $var=$val
        fi
    fi
done < $conf


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
which make > /dev/null;
if [ $? != "0" ]; then
    echo -e "$RED""Can not proceed without make, please install and re-run (Mac users see: http://developer.apple.com/technologies/xcode.html)""$NORMAL"
    exit 1;
fi

#check for g++
which g++ > /dev/null;
if [ $? != "0" ]; then
echo -e "$RED""Can not proceed without g++, please install and re-run""$NORMAL"
exit 1;
fi

# python
which python > /dev/null;
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
which R > /dev/null;
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
PREFIX_BIN=${INSTALL_PREFIX}

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
    die "Cannot write to directory $PREFIX_BIN. Maybe missing super-user (root) permissions to write there.";
fi 


################ samtools  ###################

wasInstalled=0;
which samtools > /dev/null
if [ $? = "0" ]; then
        
    samver=`samtools --version | grep samtools | cut -d" " -f2`
    vercomp $samver "1.9"
    if [[ $? == 2 ]]; then
    echo -e "$RED""samtools v1.9 or higher is needed [$samver detected].""NORMAL"
    exit 1;
    fi

    echo -e "$BLUE""Samtools appears to be already installed. ""$NORMAL"
    wasInstalled=1;
fi

if [ $wasInstalled == 0 ]; then
    echo "Installing samtools ..."
    #From sources
    $get samtools-1.9.tar.bz2  http://sourceforge.net/projects/samtools/files/samtools/1.9/samtools-1.9.tar.bz2/download
    tar -xvjpf samtools-1.9.tar.bz2
    cd samtools-1.9
    make
    cd ..
    mv samtools-1.9 $PREFIX_BIN
    export PATH=$PREFIX_BIN/samtools-1.9/:$PATH
    wasInstalled=0;
fi

if [ $wasInstalled == 0 ]; then
    check=`samtools view -h 2>&1 | grep -i options`;
    if [ $? = "0" ]; then
    echo -e "$BLUE""samtools appears to be installed successfully""$NORMAL"
    else
    echo -e "$RED""samtools NOT installed successfully""$NORMAL"; exit 1;
    fi
fi


#########################################################################################
################ bedtools 
#########################################################################################

wasInstalled=0;
which bedtools > /dev/null
if [ $? = "0" ]; then

    bedver=`bedtools --version | cut -d" " -f2 | cut -d"-" -f1`
    vercomp $bedver "1.0"
    if [[ $? == 2 ]]; then
    echo -e "$RED""bedtools v2.27.1 or higher is needed [$bedver detected].""NORMAL"
    exit 1;
    fi

    echo -e "$BLUE""bedtools appears to be already installed. ""$NORMAL"
    wasInstalled=1;
fi

if [ $wasInstalled == 0 ]; then
    echo "Installing bedtools ..."
    #From sources
    $get bedtools-2.27.1.tar.gz https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz 
    tar -xzvf bedtools-2.27.1.tar.gz
    cd bedtools-2.27
    make
    cd ..
    mv bedtools-2.27 $PREFIX_BIN
    export PATH=$PREFIX_BIN/bedtools-2.27/:$PATH
    wasInstalled=0;
fi

if [ $wasInstalled == 0 ]; then
    check=`bedtools | grep -i options`;
    if [ $? = "0" ]; then
    echo -e "$BLUE""bedtools appears to be installed successfully""$NORMAL"
    else
    echo -e "$RED""bedtools NOT installed successfully""$NORMAL"; exit 1;
    fi
fi



###################################################################################
## check aligner, if none of bwa/bowtie/bowtiew was detected, 
## install bwa as default aligner
###################################################################################
wasInstalled=0;
which bwa > /dev/null
if [ $? = "0" ]; then
    echo -e "$BLUE""bwa appears to be already installed. ""$NORMAL"
    wasInstalled=1;
fi

which bowtie2 > /dev/null
if [ $? = "0" ]; then
    echo -e "$BLUE""bowtie2 appears to be already installed. ""$NORMAL"
    wasInstalled=1;
fi

which bowtie > /dev/null
if [ $? = "0" ]; then
    echo -e "$BLUE""bowtie appears to be already installed. ""$NORMAL"
    wasInstalled=1;
fi

if [ $wasInstalled == 0 ]; then
    echo "None of bwa/bowtie2/bowtie was detected!"
    echo "Installing bwa as default aligner ..."
    #From sources
    $get bwa-0.7.17.tar.bz2 https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2/download
    tar -xvjpf bwa-0.7.17.tar.bz2
    cd bwa-0.7.17
    make
    cd ..
    mv bwa-0.7.17 $PREFIX_BIN
    export PATH=$PREFIX_BIN/bwa-0.7.17/:$PATH
    wasInstalled=0;
fi

if [ $wasInstalled == 0 ]; then
    which bwa > /dev/null
    if [ $? = "0" ]; then
        echo -e "$BLUE""bwa appears to be installed successfully""$NORMAL"
    else
        echo -e "$RED""bwa NOT installed successfully""$NORMAL"; exit 1;
    fi
fi


## Clean up
cd ..
rm -rf ./tmp

echo -e "$RED""Dependencies checked !""$NORMAL"



