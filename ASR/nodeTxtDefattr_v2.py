# Python script to convert exported ancestral node to defattr format
# Lufan Xiao
# lux011@brandeis.edu

import argparse
from itertools import groupby
import os.path
import sys

# Helper method to decode fasta


def fasta_iter(rFile):
    # modified from Brent Pedersen
    # Correct Way To Parse A Fasta File In Python
    # given a fasta file. yield tuples of header, sequence
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    fList = []
    faiter = (x[1] for x in groupby(rFile, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()

        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())

        fList.append([headerStr, seq])
    return fList

# A small struct like class to store PostPs and original MSA numbering.


class SMPPosts():
    def __init__(self, site=-1, aa1='X', aa2='X', aa3='X', post1='0.000', post2='0.000', post3='0.000') -> None:
        self.site = site
        self.aa1 = aa1
        self.aa2 = aa2
        self.aa3 = aa3
        self.post1 = post1
        self.post2 = post2
        self.post3 = post3

    def __str__(self) -> str:
        return 'MSA SITE#' + str(self.site) + ', PostP: ' + self.aa1 + '=' + self.post1 + ' ' + self.aa2 + '=' + self.post2 + ' ' + self.aa3 + '=' + self.post3


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
                '[WARN] Do you really want to overwrite? Type Y (not y) to overwrite, type anything else to exit. (Y/n) ')
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
postList = []
isFasta = False
fDecoded = []
Lines = inF.readlines()
for line in Lines:
    #keep track of line#
    lCount += 1

    # alignlen
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

    # SITE and Post Probs
    elif line[:10].find('SITE:') != -1:
        if args.verbose:
            print('Found posterior probability data at line ' + str(lCount))
        siteN = -1
        post = SMPPosts()

        # parse SITE #
        try:
            siteN = int(line[6:12])
        except Exception as err:
            print("[WARN] There's an error parsing site number from the input file at line " +
                  str(lCount)+', skipping line')
            if args.verbose:
                print(err)
            continue
            # skip line if AA# cannot be parsed

        # Parse AAs and Posts
        try:
            post.site = siteN
            post.aa1 = line[12]
            post.aa2 = line[14]
            post.aa3 = line[16]
            post.post1 = line[18:23]
            post.post2 = line[25:30]
            post.post3 = line[32:37]
        except Exception as err:
            print("[WARN] There's an error parsing posterior probability from the input file at line " +
                  str(lCount)+', SITE #' + str(siteN)+'renumbered AA#'+str(len(postList)+1)+', some default values were used at this position.')
            if args.verbose:
                print(err)

        if args.verbose:
            print('Posterior probability data parsed from line ' +
                  str(lCount) + ': AA List#' + str(len(postList)) + str(post))
        postList.append(post)
        readLen += 1

    # read fastas

    elif line[0] == '>':
        fDecoded.append([line[1:].strip(), ''])
        isFasta = True
    elif isFasta:
        if str.isalpha(line.strip()):
            fDecoded[len(fDecoded) - 1][1] = fDecoded[len(fDecoded) -
                                                      1][1] + line.strip()
        else:
            isFasta = False


# List all fasta in file, ditch_gaps or sigs
fSeqs = []
fNames = []
print('Fasta formatted protein sequence in file ' + args.input + ':')
fCount = 0
for f in fDecoded:
    if f[0].find('smp_gaps') == -1 and f[0] != 'sigs':
        fCount += 1
        print('[' + str(fCount) + '] ' + f[0])
        if args.verbose:
            print(f[1])
            print('Normal FASTA sequence.')
        fNames.append(f[0])
        fSeqs.append(f[1])
fOut = -1
while fOut < 0:
    fInput = input(
        'Please input the sequence name or index (in brackets) to generate and export posterior probabilities: ')
    try:
        fOut = int(fInput) - 1
    except ValueError as error:
        fOut = fNames.index(fInput)

    if fOut < 0 or fOut >= fCount:
        print('Invalid input: sequence name or index not found.')
        fOut = -1

# generate output
if args.verbose:
    print('Generating outputs:')
outList = []
for i in range(0, len(fSeqs[fOut])):
    if fSeqs[fOut][i] == postList[i].aa1:
        outList.append('\t:' + str(i+1) + '\t' + postList[i].post1 + '\n')
    elif fSeqs[fOut][i] == postList[i].aa2:
        outList.append('\t:' + str(i+1) + '\t' + postList[i].post2 + '\n')
    elif fSeqs[fOut][i] == postList[i].aa3:
        outList.append('\t:' + str(i+1) + '\t' + postList[i].post3 + '\n')
    else:
        outList.append('\t:' + str(i+1) + '\t' + '0.000' + '\n')

    if args.verbose:
        print('Output line ' + str(i+3) +
              ' (AA# '+str(i+1) + '): ' + str(outList[i]))

outF.writelines(outList)
outF.close()

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

exit(0)
