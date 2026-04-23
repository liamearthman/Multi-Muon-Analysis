#
# example python code to read the NOvA cosmic ray root ntuples
#
# This uses the uproot package, https://uproot.readthedocs.io/en/latest/
#   https://masonproffitt.github.io/uproot-tutorial/03-trees/index.html
# is a really useful tutorial.
#
# You probably need to:
#
#   pip install uproot awkward numpy matplotlib
#
# in order to use this script and this data.
#
# see varlist.txt for the list of what's in there.
#
# five "runs" are included in this example.
# Each has 64 2-3 minute long 64 subruns (files), but is a contiguous
# block of data taking.


# event = each 550 microsecond run
# total_muons = the entire amount of muons or particles detected in the detecter based on the reconstruction algorithm
#



import uproot
import numpy as np
import matplotlib.pyplot as plt
import glob
from datetime import datetime
from datetime import timezone
import math

# python likes functions at the top and the real program at the bottom, so
# here are some plotting functions

# plot a histogram of the number of tracks


# plot a histogram of the length per track



# plot a histogram of the zenith angle of each track



# THE MAIN PROGRAM STARTS HERE

# open the root files.  This is one run's directory worth, 2-3 minutes of
# data in 64 subrun files.  There are more such directories available you
# could add
# root_pattern=r"C:\Users\liame\Muons\Muons1\cosmicfilter_r14604_s23.cosmic_ntuple.root"
# root_pattern=r"C:\Users\liame\Muons\ntuples\cosmicntuple_fardet_r00048628_s01_t02_S20-10-30_v1*.root" # all the root files in the data directory
def pnfs_to_xrootd(path):
    return "root://fndca1.fnal.gov:1094/" + path

root_pattern="/pnfs/nova/persistent/production/exotics_dropbox/2/**/*.root"
root_files = [pnfs_to_xrootd(f) for f in glob.glob(root_pattern, recursive=True)]
root_files.sort()
if not root_files:
    raise FileNotFoundError("No ROOT files found named ",root_pattern)



# make numpy arrays to collect things to be histogrammed

NumberOfTracks = [0]
LengthOfTracks = [0]
ThetaOfTracks = [0]
num_angular_multi_mu = [0]
multi_mu_events = []
multi_mu_times = []
event_times_dictionary = {}
double_list = []
run_set = set()

# for each subrun file, open it with uproot and see if the tree is in there
num_tracks = 0
max_tracks = 0
total_muon_events = 0
total_snapshots = 0
length = 0
same_time = 0
total_multi_muons = 0
angle_threshold = 100
for file in root_files:
    try:
        with uproot.open(file) as f:
            tree_name = "cosmicntuple/cosmicTree"
            if tree_name not in f:
                print("No tree in file ", f)
                continue

        event_list = []


        # It's here!
        tree = f[tree_name]
        branches = tree.arrays(library="np")
        # check that the variables we want are in here
        # not strictly necessary, just a safety check
        if not {"run",
                "subrun",
                "tsHigh",
                "tsLow",
                "ntrack",
                "startDirY",
                "totalLength"}.issubset(tree.keys()):
            print("can't find variable list in tree")
            continue

        # access event level variables here.  Each branch is an array
        # of size = (number of events in the file)

        # get the run and subrun number
        run_branch = branches["run"]
        subrun_branch = branches["subrun"]

        # get the two words of the time stamp of the event.
        # First is unix time_t after the 1 Jan 170 epoch.
        # second is number of nanoseconds after the second
        ts_high_branch = branches["tsHigh"]
        ts_low_branch = branches["tsLow"]
        # times of tracks within the event are in the meanTNS variable
        ts_mean_branch = branches["meanTNS"]

        # get the number of muon tracks in each event
        ntrack_tree = branches["ntrack"]
        totallength_tree = branches['totalLength']
        y_dir = branches["startDirY"]


        num_tracks += ntrack_tree.sum()
        # print(ts_mean_branch)
        # print(ntrack_tree)
        # print(totallength_tree)

        # Make a python timestamp from the first event
        # and tell us about the file just opened
        fileStart = datetime.fromtimestamp(
            ts_high_branch[0]+ts_low_branch[0]/1e9,
            timezone.utc)
        print("Run ",run_branch[0]," subrun ", subrun_branch[0],
              "at ", fileStart.strftime("%x %X.%f")," has ",
              len(ntrack_tree)," events.")

        if max(ntrack_tree) > max_tracks:
            max_tracks = max(ntrack_tree)


        # Loop over events in that file
        eventCount = 0

        # ntracks iterates over the events per ntrack_tree -> ntrack_tree stores data like [15 13 15 ... 22 12 16] - where each value represents a number of muons
        event_array = np.array(branches["event"])
        print(branches["event"])
        for ntracks in ntrack_tree:
            NumberOfTracks.append(ntracks)


            # here we can loop over the arrays of track-level info
            # this loop essentially loops over each muon in the 550microsecond window and adds
            # its angle and length to a list
            yangles = []
            azangles = []
            times = []
            for i in range(ntracks):
                # gets the time of one element of the ntrack_tree
                # current_time = branches["meanTNS"][eventCount][i]

                # builds a list of each angle per 550 microsecond event
                # ydir = branches["startDirY"][eventCount][i]
                # ang = math.degrees(math.acos(-ydir))
                # xdir = branches["startDirX"][eventCount][i]
                # zdir = branches["startDirZ"][eventCount][i]
                # xang = math.degrees(math.acos(xdir))
                # zang = math.degrees(math.acos(zdir))
                # az_angle = np.arctan2(zang, xang)

                # yangles.append(ang)
                # azangles.append(az_angle)
                mu_time = branches["meanTNS"][eventCount][i]
                times.append(mu_time)

                # LengthOfTracks.append(branches["totalLength"][eventCount][i])
                # length += branches["totalLength"][eventCount][i]
                # ThetaOfTracks.append(ang)


            sorted_pairs = sorted([(times[i], i) for i in range(ntracks)], key=lambda p: p[0])
            # print(sorted_pairs)
            clusters = []
            i=0
            while i < ntracks:
                j = i + 1
                while j < ntracks and (int(sorted_pairs[j][0]) - int(sorted_pairs[i][0])) < 100:
                    j += 1

                size = j - i



                if size >= 6:
                    NA = 'N/A'
                    run0 = int(run_branch[0])
                    subrun0 = int(subrun_branch[0])
                    art_event = int(event_array[eventCount])
                    # time0 = float(sorted_pairs[i][0]) / 1000.0

                    double_list.append(
                    (run0, subrun0, art_event, size, fileStart.strftime("%m")))
                    total_multi_muons += size

                    run_set.add((run0, subrun0, len(ntrack_tree)))




                i = j


                # for ci in clusters:
                #     num_angular_multi_mu.append(len(ci))

            event_list.append(eventCount)
            eventCount += 1
            total_muon_events += ntracks
            total_snapshots += 1
            # print(same_time)

        # print('whats good')
        # for i in(event_list):
        #     key = (run0, subrun0, eventCount, time0)
        #     # i is internal event number since we are only looping over events that are stored in this file, while our artdaq files consider
        #     # all event numbers, so we must convert between the two.
        #     if key in multi_mu_events:
        #         double_list.append([key, f"event number: {int(event_array[i])}"])








# print(f'total events: {total_events}')
# call the plotting functions
# print(f"maximum number of tracks: {max_tracks}")
# avg_tracks = total_muon_events / total_snapshots
# avg_length = length / total_muon_events
# avg_multi_mu = sum(num_angular_multi_mu) / len(num_angular_multi_mu)
# print(f"avg length: {avg_length}")
# print(f"avg tracks: {avg_tracks}")
# print(f"total multi-muon events {total_events_tracks}")
# print(f"total event muons {total_events}")
# plot_nTracks(np.array(NumberOfTracks))
# plot_lengthTracks(np.array(LengthOfTracks))
# plot_thetaTracks(np.array(ThetaOfTracks))
# double_list = []
# plot_multimu(np.array(num_angular_multi_mu))
# for i in(event_list):
#     if i in multi_mu_events:
#         # i is internal event number since we are only looping over events that are stored in this file, while our artdaq files consider
#         # all event numbers, so we must convert between the two.
#         double_list.append([f"event number: {int(event_array[i])}", event_times_dictionary[i]])

    except OSError as e:
        print(f"Skipping {file}: {e}")
        continue



print(double_list)
print(len(double_list))
print(f"{total_multi_muons / len(double_list):.3f}")

with open("muon_list_final.txt", "w") as f:
    for run, subrun, event, size, month in double_list:
        f.write(f"{run} {subrun} {event} {size} {month}\n")


with open("muon_event_size.txt", "w") as f:
    for run, subrun, size in run_set:
        f.write(f"{run} {subrun} {size}\n")
