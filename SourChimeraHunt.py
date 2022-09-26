#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Juliane Hartlich
# date : Oktober 2018
# developed at Leibnz Institute DSMZ-German Collection of Migcroorganisms an Cell Cultures
# in cooperation with Technical University Braunschweig


# ---packages---------------------
# from Bio.Seq import Seq
import os, sys, argparse, sourmash, time


def command():
    Parser = argparse.ArgumentParser(
        prog="SourChimeraHunt.py",
        usage="./SourChimeraHunt.py [INPUT FASTA] [DATABASE] [-Options]",
        description="%(prog)s is used to extract chimera sequences from a multiple CCS containing FASTA file. Three output files are optained: nonchimera.fasta, multihit_chimera.fasta and misc_chimera.fasta",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    Parser.add_argument(
        "INPUT",
        nargs=1,
        help="name of input FASTA file containig multiple FASTA sequences, \nplease add directory if file is not in current working directory",
    )
    Parser.add_argument(
        "DATABASE",
        nargs=1,
        help="name of database in FASTA format used for reference based chimera detection via vsearch, \nplease add directory if file is not in current working directory",
    )
    Parser.add_argument(
        "-sdb",
        action="store",
        default="",
        dest="secondb",
        help="name of database in FASTA format used for secondary reference based chimera detection via vsearch\nplease add directory if file is not in current working directory",
    )
    Parser.add_argument(
        "-k",
        action="store",
        type=int,
        default=20,
        dest="ksize",
        help="k-mer size for hash calculation\ndefault: 20 nt",
    )
    Parser.add_argument(
        "-l",
        action="store",
        type=int,
        default=1400,
        dest="length",
        help="length cut off to filter to short sequences\ndefault: 1400 nt",
    )
    Args = Parser.parse_args()
    return Args


# getting system time and creating time stamp for log file
def time_stamp():
    Time = time.time()  # get "Unix time" in seconds
    (Yr, Mo, Da, Hr, Mi, Sc, WD, DY, ST) = time.localtime(
        Time
    )  # convert "Unix time" in time units & info
    return "%4d-%02d-%02d %02d:%02d:%02d" % (Yr, Mo, Da, Hr, Mi, Sc)


# ---main program (DEV)------------------
# parsing commandline arguments and options
ARGs = command()

# getting process ID
CWD = os.getcwd()
PID = os.getpid()

# getting name of inputfile as output label
Path = ARGs.INPUT[0]
InFileN = (Path.split("/"))[-1]
# getting name of database
DBN = ARGs.DATABASE[0]
if DBN.find("/") != -1:
    DBN = DBN.split("/")[-1]
# getting name of second database if invoked
if ARGs.secondb != "":
    SDBN = ARGs.secondb
    if SDBN.find("/") != -1:
        SDBN = SDBN.split("/")[-1]

# creating Output directory, tagged by process ID
OUTdir = "SourChimeraHunt_%d_%s" % (PID, InFileN.replace(".fasta", ""))
if not os.path.exists(OUTdir):
    os.makedirs(OUTdir)

# creating log file
LOGfile = "SourChimeraHunt_%d_log.txt" % PID
LOG = open(os.path.join(OUTdir, LOGfile), "w")
LOG.write(
    "SourChimeraHunt.py\tAlgorithm for detection and extraction of chimeric sequences\n- using sourmash and vsearch\n"
)
if ARGs.secondb != "":
    LOG.write(
        "- name of inputfile: %s\tname of used databases: %s and \n\n"
        % (InFileN, DBN, SDBN)
    )
else:
    LOG.write("- name of inputfile: %s\tname of used database: %s\n\n" % (InFileN, DBN))

File = ""
sys.stdout.write("\nReading %s\n" % InFileN)
LOG.write("%s\tReading %s\n" % (time_stamp(), InFileN))
if not os.path.isfile(Path):
    sys.stdout.write("\nFile: %s not found. Please check directory!\n" % Path)
    LOG.write(
        "%s\tFile: %s not found. Please check directory!\n" % (time_stamp(), Path)
    )
    sys.exit(0)
else:
    readfile = open(Path, "r")
    File = readfile.read()
    readfile.close()
    File = File.replace("\r", "")

# creating FASTA file containing CCS minus reverse complementary chimera
ncCCS = "%s/%s_sour_nonchimera.fasta" % (OUTdir, InFileN.replace(".fasta", ""))
# creating FASTA file comtaining reverse comlementaty chimera CCSs
mhcCCS = "%s/%s_multihit_chimera.fasta" % (OUTdir, InFileN.replace(".fasta", ""))
# creating FASTA file comtaining to short CCSs
tsCCS = "%s/%s_short_CCS.fasta" % (OUTdir, InFileN.replace(".fasta", ""))

# processing input file
Files = File.split(">")
Files = Files[1:]
sys.stdout.write("\nnumber of sequences in file: %d\n" % len(Files))
LOG.write("%s\tNumber of sequences in file: %d\n" % (time_stamp(), len(Files)))

# filtering multihit chimeras using sourmash
LOG.write("%s\tFiltering multihit chimeras using sourmash\n" % (time_stamp()))
sys.stdout.write("\nfiltering for multihit chimeras using sourmash\n")
for CurSeq in Files:
    CurSeq = CurSeq.split("\n")
    NameTag = "%s" % CurSeq[0]
    SEQ = "".join(CurSeq[1:])

    if len(SEQ) < ARGs.length:
        tsCCSwr = open(tsCCS, "a")
        tsCCSwr.write(">%s\n%s\n" % (NameTag, SEQ))
        tsCCSwr.close()
    else:
        mxhshn = len(SEQ) - ARGs.ksize + 1
        Hshs = sourmash.MinHash(n=mxhshn, ksize=ARGs.ksize)
        Hshs.add_sequence(SEQ)

        if len(Hshs.get_hashes()) < mxhshn:
            mhcCCSwr = open(mhcCCS, "a")
            mhcCCSwr.write(">%s\n%s\n" % (NameTag, SEQ))
            mhcCCSwr.close()
        elif len(Hshs.get_hashes()) == mxhshn:
            ncCCSwr = open(ncCCS, "a")
            ncCCSwr.write(">%s\n%s\n" % (NameTag, SEQ))
            ncCCSwr.close()

# futher chimera filtering with vsearch from sourmash output
if ARGs.secondb != "":
    # getting name of second database if invoked
    SDBN = ARGs.secondb
    if SDBN.find("/") != -1:
        SDBN = SDBN.split("/")[-1]

    LOG.write(
        "%s\tFiltering further chimeras using vsearch twice with reference databases: %s and %s\n"
        % (time_stamp(), DBN, SDBN)
    )
    sys.stdout.write(
        "\nfiltering further chimeras using vsearch twice with reference databases: %s and %s\n"
        % (DBN, SDBN)
    )
    LOG.write(
        "%s\tFirst vsearch run with reference database: %s \n" % (time_stamp(), DBN)
    )
    os.system(
        "vsearch --uchime_ref %s --db %s --nonchimeras %s/%s_nonchimera_db1.fasta --chimeras %s/%s_misc_chimera_db1.fasta --borderline %s/%s_unknown_db1.fasta"
        % (
            ncCCS,
            ARGs.DATABASE[0],
            OUTdir,
            InFileN.replace(".fasta", ""),
            OUTdir,
            InFileN.replace(".fasta", ""),
            OUTdir,
            InFileN.replace(".fasta", ""),
        )
    )
    os.remove(ncCCS)
    # combining nonchimera  and unknown output of vsearch
    os.system(
        "cat %s/%s_nonchimera_db1.fasta %s/%s_unknown_db1.fasta > %s/%s_nonchimera_in.fasta"
        % (
            OUTdir,
            InFileN.replace(".fasta", ""),
            OUTdir,
            InFileN.replace(".fasta", ""),
            OUTdir,
            InFileN.replace(".fasta", ""),
        )
    )
    os.remove("%s/%s_unknown_db1.fasta" % (OUTdir, InFileN.replace(".fasta", "")))
    os.remove("%s/%s_nonchimera_db1.fasta" % (OUTdir, InFileN.replace(".fasta", "")))

    # secondary analysis
    LOG.write(
        "%s\tSecond vsearch run with reference databases: %s \n" % (time_stamp(), SDBN)
    )
    os.system(
        "vsearch --uchime_ref %s/%s_nonchimera_in.fasta --db %s --nonchimeras %s/%s_nonchimera_db2.fasta --chimeras %s/%s_misc_chimera_db2.fasta --borderline %s/%s_unknown_db2.fasta"
        % (
            OUTdir,
            InFileN.replace(".fasta", ""),
            ARGs.secondb,
            OUTdir,
            InFileN.replace(".fasta", ""),
            OUTdir,
            InFileN.replace(".fasta", ""),
            OUTdir,
            InFileN.replace(".fasta", ""),
        )
    )
    os.remove("%s/%s_nonchimera_in.fasta" % (OUTdir, InFileN.replace(".fasta", "")))
    # combining nonchimera  and unknown output of vsearch
    os.system(
        "cat %s/%s_nonchimera_db2.fasta %s/%s_unknown_db2.fasta > %s/%s_nonchimera.fasta"
        % (
            OUTdir,
            InFileN.replace(".fasta", ""),
            OUTdir,
            InFileN.replace(".fasta", ""),
            OUTdir,
            InFileN.replace(".fasta", ""),
        )
    )
    os.remove("%s/%s_unknown_db2.fasta" % (OUTdir, InFileN.replace(".fasta", "")))
    os.remove("%s/%s_nonchimera_db2.fasta" % (OUTdir, InFileN.replace(".fasta", "")))

# further chimera filtering with vsearch if no second data base is invoked
elif ARGs.secondb == "":
    LOG.write(
        "%s\tFiltering for more chimeras using vsearch and reference database: %s\n"
        % (time_stamp(), DBN)
    )
    sys.stdout.write(
        "\nfiltering further chimeras using vsearch and reference database: %s\n" % DBN
    )
    os.system(
        "vsearch --uchime_ref %s --db %s --nonchimeras %s/%s_nonchimera1.fasta --chimeras %s/%s_misc_chimera.fasta --borderline %s/%s_unknown.fasta"
        % (
            ncCCS,
            ARGs.DATABASE[0],
            OUTdir,
            InFileN.replace(".fasta", ""),
            OUTdir,
            InFileN.replace(".fasta", ""),
            OUTdir,
            InFileN.replace(".fasta", ""),
        )
    )
    os.remove(ncCCS)

    # combining nonchimera  and unknown output of vsearch
    os.system(
        "cat %s/%s_nonchimera1.fasta %s/%s_unknown.fasta > %s/%s_nonchimera.fasta"
        % (
            OUTdir,
            InFileN.replace(".fasta", ""),
            OUTdir,
            InFileN.replace(".fasta", ""),
            OUTdir,
            InFileN.replace(".fasta", ""),
        )
    )
    os.remove("%s/%s_unknown.fasta" % (OUTdir, InFileN.replace(".fasta", "")))
    os.remove("%s/%s_nonchimera1.fasta" % (OUTdir, InFileN.replace(".fasta", "")))

LOG.write("%s\tProgram finished, shutting down" % time_stamp())
LOG.close()
sys.stdout.write("\nProgram finished\n")
