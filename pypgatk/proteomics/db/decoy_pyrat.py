#!/software/bin/python3.2
#
# DecoyPYrat - Fast Hybrid Decoy Sequence Database Creation for Proteomic Mass Spectromtery Analyses
#
# MIT License
#
# Copyright (c) 2016 James Christopher Wright - Wellcome Trust Sanger Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#


# used to get cmd line arguments
import argparse
# used to shuffle peptides
import random
# used to rename/delete tmp file
import os

# Read command line arguments and create help documentation using argparse
parser = argparse.ArgumentParser(
    description='''Create decoy protein sequences. Each protein is reversed and the cleavage sites switched with preceding amino acid. 
		Peptides are checked for existence in target sequences if found the tool will attempt to shuffle them.
		James.Wright@sanger.ac.uk 2015''')
parser.add_argument('fasta', metavar='*.fasta|*.fa',
                    help='FASTA file of target proteins sequences for which to create decoys')
parser.add_argument('--cleavage_sites', '-c', dest='csites', default='KR',
                    help='A list of amino acids at which to cleave during digestion. Default = KR')
parser.add_argument('--anti_cleavage_sites', '-a', dest='noc', default='',
                    help='A list of amino acids at which not to cleave if following cleavage site ie. Proline. Default = none')
parser.add_argument('--cleavage_position', '-p', dest='cpos', default='c', choices=['c', 'n'],
                    help='Set cleavage to be c or n terminal of specified cleavage sites. Default = c')
parser.add_argument('--min_peptide_length', '-l', dest='minlen', default=5, type=int,
                    help='Set minimum length of peptides to compare between target and decoy. Default = 5')
parser.add_argument('--max_iterations', '-n', dest='maxit', default=100, type=int,
                    help='Set maximum number of times to shuffle a peptide to make it non-target before failing. Default=100')
parser.add_argument('--do_not_shuffle', '-x', dest='noshuf', default=False, action='store_true',
                    help='Turn OFF shuffling of decoy peptides that are in the target database. Default=false')
parser.add_argument('--do_not_switch', '-s', dest='noswitch', default=False, action='store_true',
                    help='Turn OFF switching of cleavage site with preceding amino acid. Default=false')
parser.add_argument('--decoy_prefix', '-d', dest='dprefix', default='XXX',
                    help='Set accesion prefix for decoy proteins in output. Default=XXX')
parser.add_argument('--output_fasta', '-o', dest='dout', default='decoy.fa',
                    help='Set file to write decoy proteins to. Default=decoy.fa')
parser.add_argument('--temp_file', '-t', dest='tout', default='tmp.fa',
                    help='Set temporary file to write decoys prior to shuffling. Default=tmp.fa')
parser.add_argument('--no_isobaric', '-i', dest='iso', default=False, action='store_true',
                    help='Do not make decoy peptides isobaric. Default=false')
parser.add_argument('--memory_save', '-m', dest='mem', default=False, action='store_true',
                    help='Slower but uses less memory (does not store decoy peptide list). Default=false')
args = parser.parse_args()


# Tryptic Digest - Can be modified to take 'sites' as argument and digest based on that
def digest(protein, sites, pos, no, min):
    """Return a list of cleaved peptides with minimum length in protein sequence.
        protein = sequence
        sites = string of amino acid cleavage sites
        pos = n or c for n-terminal or c-terminal cleavage
        no = amino acids following site that would prevent cleavage ie proline
        min = minimum length of peptides returned"""

    # for each possible cleavage site insert a comma with before or after depending on pos
    for s in sites:
        r = s + ','
        if pos == 'n':
            r = ',' + s
        protein = protein.replace(s, r)

    # for each possible cleavage and all none cleavage remove comma
    for s in sites:
        for n in no:
            a = s + ',' + n
            if pos == 'n':
                a = ',' + s + n
            r = s + n
            protein = protein.replace(a, r)

    # filter peptides into list by minimum size
    return list(filter(lambda x: len(x) >= min, (protein.split(','))))


def revswitch(protein, noswitch, sites):
    """Return a reversed protein sequence with cleavage residues switched with preceding residue"""
    # reverse protein sequence with a reverse splice convert to list
    revseq = list(protein[::-1])

    if noswitch == False:

        # loop sequence list
        for i, c in enumerate(revseq):
            # if value is cleavage site switch with previous amino acid
            for s in sites:
                if c == s:
                    aa = revseq[i - 1]
                    revseq[i - 1] = revseq[i]
                    revseq[i] = aa

    # return reversed with/without switched proteins as string
    return ''.join(revseq)


def shuffle(peptide):
    """shuffle peptide without moving c-terminal amino acid cleavage site"""
    # extract terminal aa
    s = peptide[-1]
    # convert peptide to list (remove K/R) and shuffle the list
    l = list(peptide[:-1])
    random.shuffle(l)
    # return new peptide
    return ''.join(l) + s


# Create empty sets to add all target and decoy peptides
upeps = set()
dpeps = set()

# Counter for number of decoy sequences
dcount = 0;

# empty protein sequence
seq = ''

# open temporary decoy FASTA file
outfa = open(args.tout, 'w')

# Open FASTA file using first cmd line argument
fasta = open(args.fasta, 'r')
# loop each line in the file
for line in fasta:
    # if this line starts with ">" then process sequence if not empty
    if line[0] == '>':
        if seq != '':

            # make sequence isobaric (check args for switch off)
            if args.iso == False:
                seq = seq.replace('I', 'L')

            # digest sequence add peptides to set
            upeps.update(digest(seq, args.csites, args.cpos, args.noc, args.minlen))

            # reverse and switch protein sequence
            decoyseq = revswitch(seq, args.noswitch, args.csites)

            # do not store decoy peptide set in reduced memory mode
            if args.mem == False:
                # update decoy peptide set
                dpeps.update(digest(decoyseq, args.csites, args.cpos, args.noc, args.minlen))

            # write decoy protein accession and sequence to file
            dcount += 1
            outfa.write('>' + args.dprefix + '_' + str(dcount) + '\n')
            outfa.write(decoyseq + '\n')

        seq = '';

    # if not accession line then append aa sequence (with no newline or white space) to seq string
    else:
        seq += line.rstrip()

# Close files
fasta.close()
outfa.close()

# Summarise the numbers of target and decoy peptides and their intersection
nonDecoys = set()
print("proteins:" + str(dcount))
print("target peptides:" + str(len(upeps)))

# Reloop decoy file in reduced memory mode to store only intersecting decoys
if args.mem:
    # open temp decoys
    with open(args.tout, "rt") as fin:
        for line in fin:
            # if line is not accession
            if line[0] != '>':
                # digest protein
                for p in digest(line.rstrip(), args.csites, args.cpos, args.noc, args.minlen):
                    # check if in target peptides if true then add to nonDecoys
                    if p in upeps:
                        nonDecoys.add(p)
    fin.close()
    print("decoy peptides: !Memory Saving Made!")
else:
    # can only report total number in normal memory mode
    print("decoy peptides:" + str(len(dpeps)))
    # find intersecting peptides
    nonDecoys = upeps.intersection(dpeps)

print("#intersection:" + str(len(nonDecoys)))

# if there are decoy peptides that are in the target peptide set
if len(nonDecoys) > 0 and args.noshuf == False:

    # create empty dictionary with bad decoys as keys
    dAlternative = dict.fromkeys(nonDecoys, '')
    noAlternative = list()

    # loop bad decoys / dictionary keys
    for dPep in dAlternative:
        i = 0;
        aPep = dPep

        # shuffle until aPep is not in target set (maximum of 10 iterations)
        while aPep in upeps and i < args.maxit:

            # increment iteration counter
            i += 1

            # shuffle peptide
            aPep = shuffle(dPep)

            # check if shuffling has an effect if not end iterations
            if (aPep == dPep):
                i = args.maxit

        # update dictionary with alternative shuffled peptide
        dAlternative[dPep] = aPep

        # warn if peptide has no suitable alternative, add to removal list
        if i == args.maxit:
            noAlternative.append(dPep)

    print(str(len(noAlternative)) + ' have no alternative peptide')
    # remove peptides with no alternative
    for p in noAlternative:
        del dAlternative[p]

    # Free up memory by clearing large sets of peptides
    upeps.clear()
    dpeps.clear()
    # open second decoy file
    with open(args.dout, "wt") as fout:
        # open original decoy file
        with open(args.tout, "rt") as fin:
            # loop each line of original decoy fasta
            for line in fin:
                # if line is not accession replace peptides in dictionary with alternatives
                if line[0] != '>':
                    # digest decoy sequence
                    for p in digest(line.rstrip(), args.csites, args.cpos, args.noc, args.minlen):
                        # store decoy peptide for final count
                        dpeps.add(p)

                        # if decoy peptide is in dictionary replace with alternative
                        if p in dAlternative:
                            line = line.replace(p, dAlternative[p])

                fout.write(line)
        fin.close()
    fout.close()

    # delete temporary file
    os.remove(args.tout)
else:
    os.rename(args.tout, args.dout)

print("final decoy peptides:" + str(len(dpeps)))