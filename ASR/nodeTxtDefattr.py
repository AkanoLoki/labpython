# Python script to convert exported ancestran node to defattr format
# Lufan Xiao
# lux011@brandeis.edu

import argparse
import os.path
import sys

# Arguments
argP = argparse.ArgumentParser(
    description='Converts ancestral node output SMP sequence posterior probability to .defattr format')
argP.add_argument('-i', '--input', type=str,
                  help='file name and path of input file', required=True)
argP.add_argument('-o', '--output', type=str,
                  help='file name and path of output file', required=True)
argP.add_argument('-v', '--verbose', action='store_true', default=False,
                  help='prints detailed conversion progress to console.')

args = argP.parse_args(sys.argv[1:])

# Load input file
try:
    if args.verbose:
        print("input file name is " + args.input)
    inF = open(args.input, 'r')
except Exception as err:
    print("[WARN] There's an error opening or reading the input file")
    if args.verbose:
        print(err)
    exit(1)

# File existence check
overwrite = 'N'
if os.path.exists(args.output):
    # File exists, overwrite without consent is dangerous
    overwrite = input(
        '[INFO] Output file "' + args.output + '" exists. Do you want to overwrite? (Y/n) ')
    while overwrite != 'Y':
        if overwrite == 'y' or overwrite == '':
            overwrite = input(
                '[INFO] Do you want to overwrite? Type Y (not y) to overwrite, type anything else to exit. (Y/n) ')
            continue
        elif overwrite != 'Y':
            if args.verbose:
                print('Output file overwriting denied. Terminating.')
            exit(0)

# Load output file
try:
    if args.verbose:
        print("output file name is " + args.output)
    outF = open(args.output, 'w')
    # write file header
    outF.writelines(['attribute: posteriors\nrecipient: residues\n'])
except Exception as err:
    print("[WARN] There's an error opening or writing the input file")
    if args.verbose:
        print(err)
    exit(1)


# read input file
alignLen = 0
alignLenLine = -1
readLen = 0
lCount = 0
Lines = inF.readlines()
for line in Lines:
    #keep track of line#
    lCount += 1

    # Dind and parse alignlen
    if line.find('alignlen:') != -1 and alignLenLine == -1:
        if args.verbose:
            print('Found alignlen value in input file at line ' + str(lCount))
        try:
            alignLen = int(line[10:])
            alignLenLine = lCount
        except:
            print('[WARN] Error parsing the alignlen value at line ' + str(lCount))
            alignLen = -1
        if args.verbose:
            print('Parsed alignment length: ' + str(alignLen))

    # Convert prob data
    elif line.find('SITE:') != -1:
        if args.verbose:
            print('Found posterior probability data at line ' + str(lCount))
        siteN = -1
        postStr = '0.000'
        try:
            siteN = int(line[6:12])
        except Exception as err:
            print("[WARN] There's an error parsing site number from the input file at line " +
                  str(lCount)+', skipping line')
            if args.verbose:
                print(err)
            continue
            # skip line if AA# cannot be parsed
        try:
            postStr = line[18:23]
        except Exception as err:
            print("[WARN] There's an error parsing posterior probability from the input file at line " +
                  str(lCount)+', default to 0.000')
            if args.verbose:
                print(err)
        if args.verbose:
            print('Posterior probability data parsed from line ' +
                  str(lCount) + ': SITE# ' + str(siteN) + ', postProb = ' + postStr)
        outF.writelines(['\t:', str(siteN), '\t', postStr, '\n'])
        readLen += 1

# completion info
print('conversion complete, please check the output file: ' + args.output)
if args.verbose:
    if alignLen == -1:
        print('alignLen data was not found in input file.')
    elif alignLen == 0:
        print('alignLen data was found in input file at line ' +
              str(alignLenLine) + ', but could not be parsed.')
    print('alignlen read from input file at line ' +
          str(alignLenLine)+': '+str(alignLen))
    print('Total sites and posterior probability parsed: ' + str(readLen))

outF.close()
exit(0)
