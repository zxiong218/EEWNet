#!/bin/bash
python ../moni/moni_check.py check_renew_datajson data_info.json model_info.json ../plot/tmp_orgsetting.txt 6
python ../moni/netmodel.py moni model_info.json data_info.json 1
python ../plot/produce_movie.py model_info.json data_info.json tmp.mp4 600 1 
