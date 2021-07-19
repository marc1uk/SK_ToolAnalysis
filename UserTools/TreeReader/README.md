# TreeReader

A tool for reading in data from ROOT or ZEBRA files.
Access to file data is provided by two means:
1. For ROOT files an MTreeReader will be created and placed into the DataModel::Trees map.
2. For SK files (ROOT or ZEBRA), `skread` and/or `skrawread` will be called to populate fortran common blocks.
For SK ROOT files both means of data access will be available.

Typically, one entry will be loaded per Execute call (i.e per ToolChain loop).
Optionally, for SHE entries with a subsequent AFT entry, both SHE+AFT entries may be read together in one call.
The SHE data will be loaded initially, but by calling 'LoadAFT' or 'LoadSHE' the user may access the desired data.
See below for how to invoke these functions.
The reader can be configured to return only SHE+AFT pairs (skipping SHE events without an associated AFT)
or it may return all events, loading AFT data together with an SHE event only when it is available.

## Configuration
The following options are available. Defaults, when applicable, are given in parentheses.
```
verbosity 1                                    # tool verbosity (1)
readerName MyReader                            # the name to associate to this reader instance in the DataModel
inputFile /path/to/an/input/file.root          # a single input file, ROOT or ZBS
FileListName MyFileList                        # the name of a set of files prepared by the LoadFileList tool
maxEntries 10                                  # max number of entries to process before stopping the ToolChain (-1)
skFile 1                                       # whether to enable additional functionality for SK files (1)
```

When reading ROOT files the following options are also available:
```
treeName MyTree                                # the name of the tree within the file
firstEntry 10                                  # the first entry to read (0)
```

When enabling additional functionality for SK files the following options are also available:
```
LUN 10                                         # LUN to assign to the file (10)
SK_GEOMETRY 4                                  # which SK geometry this file relates to (4)
skoptn 31,30,26,25                             # options describing what to load via skread/skrawread (31)
skbadopt 23                                    # which classes of channels to mask (23)
skbadchrun 42428                               # which run to use for bad ch list for MC / calibration data
skreadMode 0                                   # which set of `skread` or `skrawread` to call (0)
skipPedestals 1                                # whether to skip pedestal and status entries (1)
readSheAftTogether 1                           # whether to read AFT data for associated SHE events together (0)
onlySheAftPairs 1                              # whether to only return SHE+AFT pairs (0)
```

When processing SK ROOT files the following additional options are also available:
```
skrootMode 0                                   # operation mode of the TreeManager (2)
outputFile /path/to/an/output/file.root        # the output file, when using the TreeManager in root2root mode
```

Notes:
* default values are given above in parentheses. If none is given, the option is mandatory.
* inputFile will take precedence over FileListName (only one is required). It supports a file path or glob.
* skFile mode will be set to 1 when reading ZBS files, based on the file extention.
* firstEntry should probably be left at 0 when in skFile mode, since event information is carried over by skread/skrawread, and event processing may fail if entries are not read sequentially from the first entry
* if maxEntries is not given or less than 0, all entries in the file will be read.
* skrootMode: 2=read, 1=write, 0=root2root copy.
* skreadMode: on each entry call... 3=both `skrawread` and `skread`, 2=`skrawread` only, 1=`skread` only, 0=`auto` - both if input file has no MC branch, only `skread` otherwise.
* if skoptn contains 25 (mask bad channels) but not 26 (get bad ch list based on current run number), then a reference run must be provided in skbadchrun. skoptn 26 cannot be used with MC data files.
* LUN will only be respected if it is not already in use. Otherwise the next free LUN will be used. Assignments start from 10.
* skipPedestals will load the next entry for which `skread` or `skrawread` did not return 3 or 4 (not pedestal or runinfo entry).
* When reading ROOT files, only enable branches you intend to use. Specify a list of input branches as follows:
```
StartInputBranchList
branchA
branchB
branchC
EndInputBranchList
```
* this will disable all branches other than `branchA`, `branchB` and `branchC`.
* for skroot files in `copy` mode, an output file will be created and entries may be copied from input to output file. Unused input branches should be disabled as above, but branches that are needed for processing but not desired in the output can be removed from the copy operation by listing only the desired output branches as follows:
```
StartOutputBranchList
branchA
branchB
EndOutputBranchList
```
* in this case branches `branchA`,`branchB` and `branchC` will be read in and accessible, but the output file will only contain branches `branchA` and `branchB`.
* outputFile is only applicable in skroot copy or write mode.
* In write mode you will need to call `skroot_set_***` and `skroot_fill_tree_` functions as required. If you need to read inputs from another file, you will need to use another TreeReader instance.
* N.B. The minimum set of active input branches for calling lf_allfit seems to be:
- SOFTWARETRG
- EVENTHEADER
- PEDESTALS
- TQLIST
- ODTQLIST
- SPACERS
- HEADER
- TQAREAL
- ATMPD
- SLE
