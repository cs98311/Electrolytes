#!/bin/bash

cp ZnMain.py ../Zn20mResults/ZnMain.py
cp ZnMain.py ../Zn5mResults/ZnMain.py
cp ZnMain.py ../Zn10mResults/ZnMain.py

cp Codes/* ../Zn20mResults/Codes
cp Codes/* ../Zn5mResults/Codes
cp Codes/* ../Zn10mResults/Codes

cp Empty.sh ../Zn20mResults/Empty.sh
cp Empty.sh ../Zn5mResults/Empty.sh
cp Empty.sh ../Zn10mResults/Empty.sh

cp Create.sh ../Zn20mResults/Create.sh
cp Create.sh ../Zn5mResults/Create.sh
cp Create.sh ../Zn10mResults/Create.sh