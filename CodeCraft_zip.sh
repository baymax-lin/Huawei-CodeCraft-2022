#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
cd $BASEDIR

if [ ! -d CodeCraft-2022 ]
then
    echo "ERROR: $BASEDIR is not a valid directory of SDK_C++ for CodeCraft-2022."
    echo "  Please run this script in a regular directory of SDK_C++."
    exit -1
fi

rm -f CodeCraft-2022.zip
zip -r CodeCraft-2022.zip CodeCraft-2022 
