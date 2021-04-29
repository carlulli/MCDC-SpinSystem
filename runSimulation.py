import pandas as pd
import subprocess
# import time
# import os
# import signal


# python script to execute the main using a params txt file


def everysecond(index):
    if index % 2 == 0:
        return True

    return False


# Read in the params.txt file with prameters as index and values as column 1
paramsdf = pd.read_csv(
    './params.txt',
    sep='=',
    names=('parameters', 'value'),
    index_col=0,
    skiprows=lambda x: everysecond(x))

# paramsdf

# TRYING TO GET THE INTEGER VALUES CHANGEDD TO INT VALUES
# SOMEHOW THIS SHIT DOES WORK
# NOT THAT IMPORTANT, AS C CODE WILL DO THE SAME
doubleparams = paramsdf.loc['N':'M_number_orient', 'value'].tolist()
doubleparams.append(paramsdf.loc['start_choice', 'value'])
doubleparams

intparams = paramsdf.loc['N':'M_number_orient', 'value'].astype(int).tolist()
intparams.append(paramsdf.loc['start_choice', 'value'].astype(int))

paramsdf['value'] = paramsdf['value'].replace(doubleparams, intparams)
# paramsdf

# paramsdf['N'] = paramsdf['N'].astype(int)
# paramsdf
#
# # paramsdf.at['N', 'value'].astype(int)
# paramsdf['value'] = paramsdf['value'].replace([])
# print(paramsdf['value'].dtypes)
# print(paramsdf)

paramsvals = paramsdf['value'].astype(str).tolist()
paramsvals

args = ['./main.exe']
for i in range(len(paramsvals)):
    args.append(paramsvals[i])
print(args)


process = subprocess.run(args)
# TRYING TO FIGURE OUT HOW TO MAKE KILLING WITH STRG + C POSSIBLE
# process = subprocess.Popen(args)
# time.sleep(5)
# os.kill(process.pid, signal.SIGINT)

print(process)
