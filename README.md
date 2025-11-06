# LMSC: Local Sketch Modularity-based Community Search

This project provides an implementation of two algorithms—IGA and SMA—for Local Sketch Modularity-based community search on large-scale graphs.


##  Running Algorithm

### (1) Running IGA
```
python run_o.py --algorithm IGA --q <query nodes> --l <lower bound> --h <upper bound> --t <sketch threshold> --network <network_path> 
```

### (2) Running SMA
```
python run_o.py --algorithm SMA --q <query nodes> --l <lower bound> --h <upper bound> --t <sketch threshold> --network <network_path> 
```

#### 1. single query node '1'
```
python run_o.py --algorithm SMA --q 1 --l <lower bound> --h <upper bound> --t <sketch threshold> --network <network_path> 
```

#### 2. multiple query nodes '1,2,3'
```
python run_o.py --algorithm SMA --q 1 2 3 --l <lower bound> --h <upper bound> --t <sketch threshold> --network <network_path> 
```



#### Required Arguments

```--network```: File path of the network

```--algorithm```: Specify the algorithm to use. (IGA,SMA)

```--l```: The lower bound of community size

```--h```: The upper bound of community size

```--q```: The query node(s) to be used. (e.g., 1,2,3)

```--t```: The sketch threshold value


## Datasets

Public benchmark datasets can be downloaded from the Stanford SNAP Dataset Repository:
https://snap.stanford.edu/data/index.html#communities

Datasets include:
- amazon
- dblp
- youtube
- LiveJournal


[//]: # ()
[//]: # ()
[//]: # (##### Reference)

[//]: # ([1] )
