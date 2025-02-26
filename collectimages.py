## in this file we will have functions to pull any images from any run given the current USDF directory 
## saves histogram as pkl in a trimmeddict{runnum} folder 
##     repo_path = "/sdf/group/rubin/repo/main" #the repo path has changed
## allows access to EO Run 7 as of 2-26-2025 (since at least 2-5-2025) 

from lsst.daf.butler import Butler
from multiprocessing import Pool
import pickle as pkl
import numpy as np
import statistics
import os
import fast_histogram
import argparse

## get dataref, trim image, generate histogram 
def getTrimmedDatarefHistogram(detectornum, runnum, rtype):  
    #repo_path = "/repo/ir2" # this was the repo before the camera moved to the summit 
    repo_path = "/sdf/group/rubin/repo/main" #the repo path has changed, allows access to EO Run 7 (and before)
    
    butler = Butler(repo_path, collections=['LSSTCam/photodiode','LSSTCam/raw/all'], instrument='LSSTCam')
    registry = butler.registry
    recordClasses = butler.registry.queryDimensionRecords('detector', where="instrument='LSSTCam'")
    det_raft_pairs = sorted([(rc.id, rc.full_name) for rc in recordClasses])
    sensorname = det_raft_pairs[detectornum][1]
    print(sensorname)
    
    where = f"exposure.science_program='{runnum}' and exposure.observation_type ='{rtype}'"
    collections = 'LSSTCam/raw/all'
    dataId = {'detector': detectornum}
    
    # get datarefs 
    datarefs = list(butler.registry.queryDatasets(datasetType='raw', collections=collections, where=where, dataId=dataId))
    exp = butler.get(datarefs[0]) #use just the first image to get the order of amps
    det = exp.getDetector()
    amps_list = ["C00", "C01", "C02","C03","C04","C05","C06", "C07", "C10","C11","C12", "C13", "C14", "C15", "C16", "C17"]
    #ampNames = [amp.getName() for amp in det if amp.getName() in amps_list] #the order in which the channels are in the dataset 
    ## cut 13% of the time for one iteration by writing explictly the order of the amps 
    ## but run above to get this order explicitly
    ampNames = ['C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C07', 'C06', 'C05', 'C04', 'C03', 'C02', 'C01', 'C00']

    # start configuration for countsdict 
    countsdictionary = {} 
    
    trimtops = ["C00", "C01", "C02","C03","C04","C05","C06", "C07"]
    trimbottoms = ["C10","C11","C12", "C13", "C14", "C15", "C16", "C17"]
    
    for dataref in datarefs:
        exp = butler.get(dataref) 
        det = exp.getDetector()
        trimmed_ims = [exp.getMaskedImage()[amp.getRawDataBBox()].getImage() for amp in det if amp.getName() in amps_list]
        int_ims = [trimmed_im.getArray().astype(int) for trimmed_im in trimmed_ims]
        trimmeddict = {ampName: int_ims for ampName, int_ims in zip(ampNames, int_ims)}
        
        npix = 2
        for amp in trimtops:
            array = trimmeddict[amp] 
            
            if amp == 'C00': 
                side = npix
                trimmed_array = array[npix:, side:] # [top:bottom, left:right]
                trimmeddict[amp] = trimmed_array

            elif amp == 'C07':
                side = npix

                trimmed_array = array[npix:, :-side]
                trimmeddict[amp] = trimmed_array

            else: 
                trimmed_array = array[npix:, :]
                trimmeddict[amp] = trimmed_array

        for amp in trimbottoms: 
            array = trimmeddict[amp] 
            
            if amp == 'C10': 
                side = npix
                trimmed_array = array[:-npix, side:]  # [top:bottom, left:right]
                trimmeddict[amp] = trimmed_array

            elif amp == 'C17':
                side = npix
                trimmed_array = array[:-npix, :-side]
                trimmeddict[amp] = trimmed_array

            else: 
                trimmed_array = array[:-npix, :]
                trimmeddict[amp] = trimmed_array

        ## form histogram 
        
        for amp in trimmeddict: 
            array = trimmeddict[amp] 
            flattenarray = array.flatten() 
            counts = fast_histogram.histogram1d(flattenarray, bins=180000, range=[2e4, 2e5])

            if amp not in countsdictionary.keys(): 
                countsdictionary[amp] = counts  
            else: 
                fullcounts = np.add(countsdictionary[amp], counts)
                countsdictionary[amp] = fullcounts  
                #fullcounts = [x + y for x, y in zip(countsdictionary[amp], counts)] # this method is slower 
                #countsdictionary[amp] = list(map(add, countsdictionary[amp], counts )) # this method is medium 
    
    return countsdictionary 

## reduce the amount stored, since bins is a np.arange() 
def reduceBinsReported(counts):
    bins = np.arange(20000, 200000)

    firstnonzeroindex = next((i for i, x in enumerate(counts) if x != 0), None)    
    reversecounts = counts[::-1]    
    lastnonzeroindex = next((i for i, x in enumerate(reversecounts) if x!= 0), None)
    # have to subtract from the length to get the proper order index 
    finalindex = len(counts)-lastnonzeroindex
    
    firstbin = bins[firstnonzeroindex]
    lastbin = bins[finalindex] 
    binrange = [firstbin, lastbin]
    reducedcounts = counts[firstnonzeroindex: finalindex]
    return reducedcounts, binrange


## compile the histograms into a single file, write directory if it does not exist
def getHistogramsForDetector(det, runType, runNumber): 
    countsdict = getTrimmedDatarefHistogram(det, runNumber, runType)
    makeedgedict = {}

    # write edge dict for each amp in the sensor 
    for key in countsdict:
        reducedcounts, binrange  = reduceBinsReported(countsdict[key])        
        makeedgedict[key] = [reducedcounts, binrange]
    
    directory = f'/sdf/data/rubin/user/rejnicho/edgedicts/trimmeddicts{runNumber}'
    
    ## Check if a specific directory exists, if not write one 
    if os.path.isdir(directory) != True:
        os.mkdir(directory)

    filename = f'sensor{det}.pkl'
    fullfilename = os.path.join(directory, filename)
    with open(fullfilename, 'wb') as f:
        pkl.dump(makeedgedict, f)

    return det 

## call function to examine combined histograms of a specific run given runtype
def generalCollectImagesRun(runNumber, runType, startingDetector, endingDetector):
    detectors = np.arange(startingDetector, endingDetector)    
    
    passDetectorValues = [(det, runType, runNumber) for det in detectors]  # Create tuples
    print(passDetectorValues) 
    
    ## Processing for all the images is the same at this stage
    ## Trim 2 edge pixels
    ## Reduce histogram to remove 0 elements
    ## Reduce bin range to [start, end] 
    ## save this to the trimmedDict location 
    
    with Pool(4) as pool: 
        results = pool.starmap(getHistogramsForDetector, passDetectorValues)

    print("we have run the amplifeirs") 
    print(results)

if __name__ == "__main__": 
    parser = argparse.ArgumentParser(description="Generate Hist from given run of detectors specified")
    parser.add_argument("runNumber", type=str, help='run number')
    parser.add_argument("runType", type=str, help='run type: flat, ramp, etc.')
    parser.add_argument("detector1", type=int, help='first detector to run over')
    parser.add_argument("detector2", type=int, help='final detector to run over')
    
    args = parser.parse_args()
    generalCollectImagesRun(args.runNumber, args.runType, args.detector1, args.detector2) 

