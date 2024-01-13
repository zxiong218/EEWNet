# EEWNet: Generalized Neural Networks for Real-Time Earthquake Early Warning and Monitoring

## Install environment
The neural networks are trained and implemented with TensorFlow. You can set up a similar environment using the following commands.
- Install in the default environment:
```bash
conda env update -f=env.yml -n base
```
- Install in the "eewnet" virtual environment: 
```bash
conda env create -f env.yml
conda activate eewnet
```
## Download models and demo data
The models and demo data are available at the following links. Download two model files (\*.hdf5) and place them in the "models" folder. Download demo datasets (\*.h5) and place them in the "data" folder.

- Detection and location neural network: https://drive.google.com/file/d/1v40sHp3yYZzmIIJzBLqFRb3G1MdqFJ4P/view?usp=sharing

- Magnitude neural network: https://drive.google.com/file/d/1flzoVOUwdKoNwGAU-qQAQ9Kdo1QYL7K2/view?usp=sharing

- Demo data for Osaka (Japan) earthquakes: https://drive.google.com/file/d/1KMU7wfm-AnhcP3qQs1XqVl1BIidK9QYh/view?usp=sharing

- Demo data for Ridgecrest (US) earthquakes: https://drive.google.com/file/d/1JGhwZR_RNRdxo9nTKb_006-rWEIvLlz0/view?usp=sharing

Or the models and datasets are available at: https://doi.org/10.5061/dryad.wwpzgmsrg (Under review)
## Examples
Running the monitoring codes is simple. They will output the detection, location, and magnitude results for each time window according to your settings. The codes first load all data streams and truncating windows from your disk, and then the neural networks perform earthquake detection and parameter evaluation simultaneously. 
```bash
cd demo
python ../moni/netmodel.py moni model_info.json data_info.json 1
```
The model_info.json contains all the information for the neural networks. data_info.json contains the data information for neural networks' inputs, such as the how you set the truncating windows, stations, and the origin of the coordinates of the monitoring area. "1" represents you use the No. 1 GPU for computation.
The output files are \*_dlm.csv, \*_dlm.pkl, and "\*_dlm.mat" representing different formats of the saved results. You can use \*_dlm.pkl to make movies to visualize the results, as follows:
```bash
python ../plot/produce_movie.py model_info.json data_info.json tmp.mp4
```
You may need to install "moviepy" if have no corresponding pakages in your system. The default length of the \*.mp4 is 600 s. You may need to modify the length and the Frames Per Second (fps). The following command outputs a \*.mp4 of 150 s and 1 fps:
```bash
python ../plot/produce_movie.py model_info.json data_info.json tmp.mp4 150 1
```

## Setting files
data_info.json contains the settings for data streams.
```
{
"h5file":"../data/japan201806180758.h5",
"stations_file":"../data/station_japan.txt",
"prelist":"2018-06-18T07.58.00.000000Z.",
"tbegin":"2018-06-18T07:58:00.000000Z",
"nwin":150,
"dtwin":2,
"ampscale":0.000000001,
"E":"E",
"N":"N",
"Z":"U",
"output_3dimags_range":[0,100,1],
"data_org":{"lonlat":[135.572, 35.064],"xy":[12365.64583869963, 3898.76616],
                                      "scale":[91.21091256822669, 111.19],"rt_clockwise":0.0},
"out_nam":"test1"
}
```
"h5file" is the file location of the datasets downloaded from NIED Hi-net. "stations_file" contains the station locations and station name. "prelist" is the prefix of dataset in "h5file". For example, 2018-06-18T07.58.00.000000Z.N.TKHH.E is the key of the data from E component of station N.TKHH. "tbegin" is the starting time of the moving truncating window. "nwin" is number of time windows to truncate, "dtwin" is the intervel of the moving window. "ampscale" is used to adjust the amplitudes of the data for magnitude estimation. "E", "N", and "Z" are the component names. "output_3dimags_range" is used to set how many 3D images of the locations required to output. "out_nam" is the prefix of the output files. "data_org" is the coordinate of the monitoring area center points, and can be used for the conversion between geographical and xyz coordinates. "data_org" should be near the mean values of xy coordinates of all the monitoring stations because station distribution should well cover the monitoring area. The settings for "data_org" could be seen in "Coordinate origin of the monitoring area".

## Coordinate origin of the monitoring area
The "data_org" setting should be close to the mean of xy coordinates of the monitoring stations, and you can use the following command to test the settings:
```bash
moni_check check_json datajs_file.json modeljs_file.json 6
```
The above command can test if the distance between "data_org" and the mean of station xy coordinates is less than 6 km. You also can use "./plot/set_moni_plot.m" to visualize your station distribution and set the monitoring center points in the set_moni_plot.m file. The set_moni_plot.m file could output \*_orgsetting.txt file for your setting. You can use it to renew your data_info.json file as follows:
```bash
python ../moni/moni_check.py check_renew_datajson data_info.json model_info.json ../plot/tmp_orgsetting.txt 6
```
Your data_info.json can be updated with the new "data_org" and "stations_file": "tmp_station.txt".

## Related paper
Xiong Zhang, Miao Zhang (2023), Generalized Neural Networks for Real-Time Earthquake Early Warning, https://arxiv.org/abs/2312.15218

