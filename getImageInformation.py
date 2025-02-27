## in this file you can find the exposure times and get number of exp. per detector
## we can also create plots showing this 
## this is given the 2-25-25 repo path 

from lsst.daf.butler import Butler
from multiprocessing import Pool 
import matplotlib.pyplot as plt
import pickle as pkl
import numpy as np 
import os 
import argparse

# get the number of datarefs for a detector in a run 
def getNumberExposures(detectornum, runnum, rtype): # gets the Datarefs and makes trimmed hist 
    #repo_path = "/repo/ir2" # this was the repo before the camera moved to the summit 
    
    repo_path = "/sdf/group/rubin/repo/main" #the repo path has changed, allows access to EO Run 7 as of 2-5-2025 
    
    butler = Butler(repo_path, collections=['LSSTCam/photodiode','LSSTCam/raw/all'], instrument='LSSTCam')
    registry = butler.registry
    recordClasses = butler.registry.queryDimensionRecords('detector', where="instrument='LSSTCam'")
    det_raft_pairs = sorted([(rc.id, rc.full_name) for rc in recordClasses])
    sensorname = det_raft_pairs[detectornum][1]
    
    where = f"exposure.science_program='{runnum}' and exposure.observation_type ='{rtype}'"
    collections = 'LSSTCam/raw/all'
    dataId = {'detector': detectornum}
    # get datarefs 
    datarefs = list(butler.registry.queryDatasets(datasetType='raw', collections=collections, where=where, dataId=dataId))

    numberdatarefs = len(datarefs)

    print(numberdatarefs) 

    return numberdatarefs

## get the exposure times and ID for all datarefs in a run
def getTimesID(detectornum, runnum, rtype):  
    #repo_path = "/repo/ir2" # this was the repo before the camera moved to the summit 
    
    repo_path = "/sdf/group/rubin/repo/main" #the repo path has changed, allows access to EO Run 7 as of 2-5-2025 
    
    butler = Butler(repo_path, collections=['LSSTCam/photodiode','LSSTCam/raw/all'], instrument='LSSTCam')
    registry = butler.registry
    recordClasses = butler.registry.queryDimensionRecords('detector', where="instrument='LSSTCam'")
    det_raft_pairs = sorted([(rc.id, rc.full_name) for rc in recordClasses])
    sensorname = det_raft_pairs[detectornum][1]
    
    where = f"exposure.science_program='{runnum}' and exposure.observation_type ='{rtype}'"
    collections = 'LSSTCam/raw/all'
    dataId = {'detector': detectornum}
    # get datarefs 
    datarefs = list(butler.registry.queryDatasets(datasetType='raw', collections=collections, where=where, dataId=dataId))
    
    exposureTimes = [] 
    exposureIDs = []
    for dataref in datarefs:
        exp = butler.get(dataref) 
        det = exp.getDetector()
        #print(exp) # yeilds lsst.afw.image._exposure.ExposureF object 
        #print(dir(exp)) # yields all the functions this object has 
        #print(exp.getInfo()) # yields lsst.afw.image.ExposureInfo object
        #print(exp.info) # yields lsst.afw.image.ExposureInfo object
        #print(dir(exp.info)) # yields all functions of obj 
        #print(exp.info.getVisitInfo()) # yields the visit information, such as exp time, ra/dec, id, etc 
        #print(exp.info.getVisitInfo().exposureTime)
        #print(exp.info.getVisitInfo().id)
        exposuretime = exp.info.getVisitInfo().exposureTime
        exposureTimes.append(exposuretime)
        exposureid = exp.info.getVisitInfo().id
        exposureIDs.append(exposureid)

    return exposureTimes, exposureIDs

def plotExposureTimes(detectornum, runnum, rtype):
    exposureTimes, exposureIDs = getTimesID(detectornum, runnum, rtype)
    plt.hist(exposureTimes)
    plt.title(f"{runnum} det {detectornum} exposure times")
    plt.xlabel("sec")
    plt.ylabel("freq")
    plt.savefig(f"{runnum}{detectornum}expTimes")
   

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get some information about a run")
    subparsers = parser.add_subparsers(dest='command', required=True)
    
    parser_numbdatarefs = subparsers.add_parser("noData", help='tells you number datarefs ina  sensor run')
    parser_numbdatarefs.add_argument("detectornum", type=int, help='detector to examine')
    parser_numbdatarefs.add_argument("runnumber", type=str, help='run number')
    parser_numbdatarefs.add_argument("rtype", type=str, help='run type: ramp, flat, etc.') 
    parser_numbdatarefs.set_defaults(func=lambda args: getNumberExposures(args.detectornum, args.runnumber, args.rtype))
   
    parser_plotTimes = subparsers.add_parser("plotTimes", help='yields histogram of exp times in a run for a sensor')
    parser_plotTimes.add_argument("detectornum", type=int, help='detector to examine')
    parser_plotTimes .add_argument("runnumber", type=str, help='run number')
    parser_plotTimes.add_argument("rtype", type=str, help='run type: ramp, flat, etc.')
    parser_plotTimes.set_defaults(func=lambda args: plotExposureTimes(args.detectornum, args.runnumber, args.rtype)) 

    args = parser.parse_args()
    args.func(args)
