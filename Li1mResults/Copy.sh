#!/bin/bash

cp LiMain.py ../Li1mResults/LiMain.py
cp LiMain.py ../Li5mResults/LiMain.py
cp LiMain.py ../Li10mResults/LiMain.py

cp Codes/* ../Li1mResults/Codes
cp Codes/* ../Li5mResults/Codes
cp Codes/* ../Li10mResults/Codes

cp Empty.sh ../Li1mResults/Empty.sh
cp Empty.sh ../Li5mResults/Empty.sh
cp Empty.sh ../Li10mResults/Empty.sh

cp Create.sh ../Li1mResults/Create.sh
cp Create.sh ../Li5mResults/Create.sh
cp Create.sh ../Li10mResults/Create.sh