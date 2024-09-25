#!/bin/bash

mkdir ../Mixed10mResults
mkdir ../Mixed20mResults

cp MixedMain.py ../Mixed10mResults/MixedMain.py
cp MixedMain.py ../Mixed20mResults/MixedMain.py

cp Codes/* ../Mixed10mResults/Codes
cp Codes/* ../Mixed20mResults/Codes

cp Empty.sh ../Mixed10mResults/Empty.sh
cp Empty.sh ../Mixed20mResults/Empty.sh

cp Create.sh ../Mixed10mResults/Create.sh
cp Create.sh ../Mixed20mResults/Create.sh

sh ../Mixed10mResults/Create.sh
sh ../Mixed20mResults/Create.sh