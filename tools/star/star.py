#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Wrapper for STAR. Assumes the STAR executable is in your PATH.

usage: star.py {genomeGenerate, align} [OPTIONS]...

STAR version:   2.4.2a
Author:         Adam Struck <strucka@ohsu.edu>
Last Modified:  06/10/2015
"""

import argparse
import logging
import subprocess


def collectArgs():
    descr = 'Spliced Transcripts Alignment to a Reference'
    parser = argparse.ArgumentParser(
        description=descr,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    subparser = parser.add_subparsers(title='runMode')
    addGenomeGenerateParser(subparser)
    addAlignParser(subparser)
    return parser


def addAlignParser(subparsers):
    parser = subparsers.add_parser(
        'align',
        description='direct STAR to run a mapping job',
        help='direct STAR to run a mapping job',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Global options
    parser.add_argument('--runMode', default='alignReads',
                        help=argparse.SUPPRESS)
    parser.add_argument('--parametersFiles', default='-',
                        help=argparse.SUPPRESS)
    parser.add_argument('--sysShell', default='-',
                        help=argparse.SUPPRESS)
    parser.add_argument('--runDirPerm', default='User_RWX',
                        help=argparse.SUPPRESS)
    parser.add_argument('--runThreadN', type=int, default=1,
                        help='Number of threads to use')

    # Genome options
    parser.add_argument('--genomeDir', required=True,
                        help='path to genome directory')
    parser.add_argument('--genomeLoad', default='NoSharedMemory',
                        help='mode of shared memory usage for the \
                        genome files')

    # Read Params
    parser.add_argument('--readFilesIn', required=True, nargs='+',
                        help='name(s) (with path) of the files \
                        containing the sequences to be mapped (e.g. \
                        RNA-seq FASTQ files). If using Illumina \
                        paired-end reads, both read1 and \
                        read2 files have to be supplied.')
    parser.add_argument('--readFilesCommand', default='-',
                        help='UncompressionCommand option, where \
                        UncompressionCommand is the un-compression \
                        command that takes the file name as input \
                        parameter, and sends the uncompressed output \
                        to stdout')
    parser.add_argument('--readMapNumber', type=int, default=-1,
                        help='number of reads to map from the \
                        beginning of the file. -1: map all reads')
    parser.add_argument('--readMatesLengthsIn', default='NotEqual',
                        choices=['Equal', 'NotEqual'],
                        help='lengths of names,sequences, qualities \
                        for both mates are the same / not the same. \
                        NotEqual is safe in all situations.')
    parser.add_argument('--clip3pNbases',
                        type=int, default=0,
                        help='number(s) of bases to clip from 3p of \
                        each mate')
    parser.add_argument('--clip5pNbases',
                        type=int, default=0,
                        help='number(s) of bases to clip from 5p of \
                        each mate')
    parser.add_argument('--clip3pAdapterSeq', default='-',
                        help='adapter sequences to clip from 3p of \
                        each mate')
    parser.add_argument('--clip3pAdapterMMp',
                        type=float, default=0.1,
                        help='max proportion of mismatches for 3p \
                        adpater clipping for each mate')
    parser.add_argument('--clip3pAfterAdapterNbases',
                        type=int, default=0,
                        help='number of bases to clip from 3p of each \
                        mate after the adapter clipping')

    # Limits
    parser.add_argument('--limitIObufferSize',
                        type=int, default=150000000,
                        help='max available buffers size (bytes) for \
                        input/output, per thread')
    parser.add_argument('--limitOutSAMoneReadBytes',
                        type=int, default=100000,
                        help='max size of the SAM record for one read.\
                        Recommended value: \
                        > (2 * (LengthMate1 + LengthMate2 + 100) * \
                        outFilterMultimapNmax)')
    parser.add_argument('--limitOutSJoneRead',
                        type=int, default=1000,
                        help='max number of junctions for one read \
                        (including all multi-mappers).')
    parser.add_argument('--limitOutSJcollapsed',
                        type=int, default=1000000,
                        help='max number of collapsed junctions')
    parser.add_argument('--limitBAMsortRAM',
                        type=int, default=0,
                        help='maximum available RAM for sorting BAM. \
                        If =0, it will be set to the genome index size. \
                        0 value can only be used with –genomeLoad \
                        NoSharedMemory option.')
    parser.add_argument('--limitSjdbInsertNsj',
                        type=int, default=1000000,
                        help='maximum number of junction to be \
                        inserted to the genome on the fly at the \
                        mapping stage, including those from annotations \
                        and those detected in the 1st step of the \
                        2-pass run.')

    # Output (general)
    parser.add_argument('--outFileNamePrefix', default='./',
                        help='outfile file name prefix')
    parser.add_argument('--outTmpDir', default='-',
                        help='path to a directory that will be used \
                        as temporary by STAR. All contents of this \
                        directory will be removed!')
    parser.add_argument('--outStd', default='Log',
                        choices=['Log', 'SAM', 'BAM_Unsorted',
                                 'BAM_SortedByCoordinate',
                                 'BAM_Quant'],
                        help='which output will be directed to stdout')
    parser.add_argument('--outReadsUnmapped', default='None',
                        choices=['None', 'Fastx'],
                        help='output of unmapped reads (besides SAM)')
    parser.add_argument('--outQSconversionAdd',
                        type=int, default=0,
                        help='add this number to the quality score')

    # Output (SAM and BAM)
    parser.add_argument('--outSAMtype', default='BAM SortedByCoordinate',
                        choices=['BAM SortedByCoordinate', 'BAM Unsorted',
                                 'SAM SortedByCoordinate', 'SAM Unsorted'],
                        help='format and how to sort output file.')
    parser.add_argument('--outSAMmode', default='Full',
                        choices=['None', 'Full', 'NoQS'],
                        help='mode of SAM output')
    parser.add_argument('--outSAMstrandField', default='None',
                        choices=['None', 'intronMotif'],
                        help='cufflinks-like strand field flag')
    parser.add_argument('--outSAMattributes', default='Standard',
                        choices=['NH', 'Standard', 'All', 'None'],
                        help='a string of desired SAM attributes, in \
                        the order desired for the output')
    parser.add_argument('--outSAMunmapped', default='None',
                        choices=['None', 'Within'],
                        help='output of unmapped reads in the SAM \
                        format')
    parser.add_argument('--outSAMorder', default='Paired',
                        choices=['Paired', 'PairedKeepInputOrder'],
                        help=' type of sorting for the SAM output')
    # TO-DO

    # BAM processing
    parser.add_argument('--bamRemoveDuplicatesType', default='-',
                        choices=['-', 'UniqueIdentical'],
                        help='if ="UniqueIdentical", mark duplicates \
                        in the BAM file, for now only works with \
                        sorted BAM feeded with inputBAMfile')
    parser.add_argument('--bamRemoveDuplicatesMate2basesN',
                        type=int, default=0,
                        help='number of bases from the 5p of mate 2 \
                        to use in collapsing (e.g. formatRAMPAGE)')

    # Output wiggle
    # TO-DO

    # Output filtering
    parser.add_argument('--outFilterType', default='Normal',
                        choices=['Normal', 'BySJout'],
                        help='reduces the number of spurious \
                        junctions')
    parser.add_argument('--outFilterMultimapScoreRange',
                        type=int, default=1,
                        help='the score range below the maximum score \
                        for multimapping alignments')
    parser.add_argument('--outFilterMultimapNmax',
                        type=int, default=10,
                        help='read alignments will be output only if \
                        the read maps fewer than this value, otherwise \
                        no alignments will be output')
    parser.add_argument('--outFilterMismatchNmax',
                        type=int, default=10,
                        help='alignment will be output only if it has \
                        fewer mismatches than this value')
    parser.add_argument('--outFilterMismatchNoverLmax',
                        type=float, default=0.3,
                        help='alignment will be output only if  its \
                        ratio of mismatches  to *mapped* length is \
                        less than this value')
    parser.add_argument('--outFilterMismatchNoverReadLmax',
                        type=float, default=1,
                        help='alignment will be output only if  its \
                        ratio of mismatches to *read* length is less \
                        than this value')
    parser.add_argument('--outFilterScoreMin', type=int, default=0,
                        help='alignment will be output only if its \
                        score is higher than this value')
    parser.add_argument('--outFilterScoreMinOverLread',
                        type=float, default=0.66,
                        help='outFilterScoreMin normalized to read \
                        length (sum of mates lengths formatpaired-end \
                        reads)')
    parser.add_argument('--outFilterMatchNmin', type=int, default=0,
                        help='alignment will be output only if the \
                        number of matched bases is higher than this \
                        value')
    parser.add_argument('--outFilterMatchNminOverLread',
                        type=float, default=0.66,
                        help='outFilterMatchNmin normalized to read \
                        length (sum of mates lengths for paired-end \
                        reads)')
    parser.add_argument('--outFilterIntronMotifs', default='None',
                        choices=['None', 'RemoveNoncanonical',
                                 'RemoveNoncanonicalUnannotated'],
                        help='filter alignment using their motifs')

    # Output filtering (splice junctions)
    parser.add_argument('--outSJfilterReads', default='All',
                        choices=['All', 'Unique'],
                        help='which reads to consider for collapsed \
                        splice junctions output')
    # TO-DO

    # Scoring
    # TO-DO

    # Alignments and seeding
    parser.add_argument('--seedSearchStartLmax',
                        type=int, default=50,
                        help='defines the search start point through \
                        the read - the read is split into pieces no \
                        longer than this value')
    # TO-DO: seed options
    parser.add_argument('--alignIntronMin',
                        type=int, default=21,
                        help='minimum intron length')
    parser.add_argument('--alignIntronMax',
                        type=int, default=0,
                        help='maximum intron length, if 0, max intron \
                        size will be determined by \
                        (2ˆwinBinNbits)*winAnchorDistNbins')
    parser.add_argument('--alignMatesGapMax',
                        type=int, default=0,
                        help='maximum genomic distance between mates, \
                        if 0, max intron gap will be determined by \
                        (2ˆwinBinNbits)*winAnchorDistNbins')
    parser.add_argument('--alignSJoverhangMin',
                        type=int, default=5,
                        help='minimum overhang for unannotated \
                        junctions')
    parser.add_argument('--alignSJDBoverhangMin',
                        type=int, default=3,
                        help='minimum overhang for annotated \
                        junctions')
    parser.add_argument('--alignSplicedMateMapLmin',
                        type=int, default=0,
                        help='minimum mapped length for a read mate \
                        that is spliced')
    parser.add_argument('--alignSplicedMateMapLminOverLmate',
                        type=float, default=0.66,
                        help='alignSplicedMateMapLmin normalized to \
                        mate length')
    parser.add_argument('--alignWindowsPerReadNmax',
                        type=int, default=10000,
                        help='max number of windows per read')
    parser.add_argument('--alignTranscriptsPerWindowNmax',
                        type=int, default=100,
                        help='max number of transcripts per window')
    parser.add_argument('--alignTranscriptsPerReadNmax',
                        type=int, default=10000,
                        help='max number of different alignments per \
                        read to consider')
    parser.add_argument('--alignEndsType', default='Local',
                        choices=['Local', 'EndToEnd',
                                 'Extend5pOfRead1'],
                        help='type of read ends alignment')
    parser.add_argument('--alignSoftClipAtReferenceEnds',
                        default='Yes', choices=['Yes', 'No'],
                        help='allow the soft-clipping of the \
                        alignments past the end of the chromosomes')

    # Windows, anchors, binning
    parser.add_argument('--winAnchorMultimapNmax',
                        type=int, default=50,
                        help='max number of loci anchors are allowed \
                        to map to')
    parser.add_argument('--winBinNbits',
                        type=int, default=16,
                        help='=log2(winBin), where winBin is the size \
                        of the bin for the windows/clustering, each \
                        window will occupy an integer number of bins.')
    parser.add_argument('--winAnchorDistNbins',
                        type=int, default=9,
                        help='max number of bins between two anchors \
                        that allows aggregation of anchors into one \
                        window')
    parser.add_argument('--winFlankNbins',
                        type=int, default=4,
                        help='=log2(winFlank), where winFlank is the \
                        size of the left and right flanking regions \
                        for each window')

    # Chimeric alignments
    parser.add_argument('--chimOutType', default='SeparateSAMold',
                        choices=['SeparateSAMold', 'WithinBAM'],
                        help='type of chimeric output')
    parser.add_argument('--chimSegmentMin',
                        type=int, default=0,
                        help='minimum length of chimeric segment \
                        length, if =0, no chimeric output')
    parser.add_argument('--chimScoreMin',
                        type=int, default=0,
                        help='minimum total (summed) score of the \
                        chimeric segements')
    parser.add_argument('--chimScoreSeparation',
                        type=int, default=10,
                        help='minimum difference (separation) between \
                        the best chimeric score and the next one')
    parser.add_argument('--chimScoreJunctionNonGTAG',
                        type=int, default=-1,
                        help='penalty for a non-GT/AG chimeric junction')
    parser.add_argument('--chimJunctionOverhangMin',
                        type=int, default=20,
                        help='minimum overhang for a chimeric junction')

    # Quantification of annotations
    parser.add_argument('--quantMode', default='-',
                        choices=['-', 'TranscriptomeSAM',
                                 'GeneCounts'],
                        help='types of quantification requested')
    parser.add_argument('--quantTranscriptomeBAMcompression',
                        type=int, default=-1,
                        help='-1 to 10 transcriptome BAM compression \
                        level, -1=default compression (6?), 0=no \
                        compression, 10=maximum compression')
    parser.add_argument('--quantTranscriptomeBan',
                        default='IndelSoftclipSingleend',
                        choices=['IndelSoftclipSingleend',
                                 'Singleend'],
                        help='prohibit various alignment type')

    # 2-pass mapping
    parser.add_argument('--twopassMode', default='None',
                        choices=['None', 'Basic'],
                        help='2-pass mapping mode')
    parser.add_argument('--twopass1readsN',
                        type=int, default=-1,
                        help='number of reads to process for the 1st \
                        step. Use very large number (or default -1) \
                        to map all reads in the first step.')

    # Splice junction database
    parser.add_argument('--sjdbGTFfile',
                        help='path to GTF file with annotations')
    parser.add_argument('--sjdbGTFchrPrefix',
                        help='prefix for chromosome names in a GTF \
                        file (e.g. ’chr’ for using ENSMEBL \
                        annotations with UCSC geneomes)')
    parser.add_argument('--sjdbGTFfeatureExon', default='exon',
                        help='feature type in GTF file to be used as \
                        exons for building transcripts')
    parser.add_argument('--sjdbGTFtagExonParentTranscript',
                        default='transcript_id',
                        help='Exons are assigned to the transcripts \
                        using parent-child relationship defined by \
                        this option')
    parser.add_argument('--sjdbGTFtagExonParentGene',
                        default='gene_id',
                        help='tag name to be used as exons’ \
                        gene-parents (Default ”gene id” works \
                        formatGTF files)')
    parser.add_argument('--sjdbFileChrStartEnd',
                        help='path to annotations formatted as a list \
                        of splice junctions coordinates in a text \
                        file')
    parser.add_argument('--sjdbOverhang', type=int, default=100,
                        help='length of the donor/acceptor sequence \
                        on each side of the junctions, \
                        ideally = (mate length - 1). \
                        If =0, splice junction database is not used')
    parser.add_argument('--sjdbInsertSave', default='Basic',
                        choices=['Basic', 'All'],
                        help='which files to save when sjdb junctions \
                        are inserted on the fly at the mapping step')


def addGenomeGenerateParser(subparsers):
    parser = subparsers.add_parser(
        'genomeGenerate',
        description='direct STAR to run a genome indices generation \
        job',
        help='direct STAR to run a genome indices generation job',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Global options
    parser.add_argument('--runMode', default='genomeGenerate',
                        help=argparse.SUPPRESS)
    parser.add_argument('--parametersFiles', default='-',
                        help=argparse.SUPPRESS)
    parser.add_argument('--sysShell', default='-',
                        help=argparse.SUPPRESS)
    parser.add_argument('--runDirPerm', default='User_RWX',
                        help=argparse.SUPPRESS)
    parser.add_argument('--runThreadN', type=int, default=1,
                        help='Number of threads to use')

    # Limits
    parser.add_argument('--limitIObufferSize',
                        type=int, default=150000000,
                        help='max available buffers size (bytes) for \
                        input/output, per thread')
    parser.add_argument('--limitGenomeGenerateRAM',
                        type=int, default=31000000000,
                        help='maximum available RAM (bytes) for genome \
                        generation')

    # Genome options
    parser.add_argument('--genomeDir', required=True,
                        help='path to genome directory')

    # Genome generation parameters
    parser.add_argument('--genomeFastaFiles', required=True,
                        help='path(s) to genome fasta files')
    parser.add_argument('--genomeChrBinNbits',
                        type=int, default=18,
                        help='=log2(chrBin), where chrBin is the size \
                        of the bins for genome storage: each \
                        chromosome will occupy an integer number of \
                        bins')
    parser.add_argument('--genomeSAindexNbases',
                        type=int, default=14,
                        help=' length (bases) of the SA pre-indexing \
                        string. Typically between 10 and 15. Longer \
                        strings will use much more memory, but allow \
                        faster searches.')
    parser.add_argument('--genomeSAsparseD',
                        type=int, default=1,
                        help='suffux array sparsity, i.e. distance \
                        between indices: use bigger numbers to \
                        decrease needed RAM at the cost of mapping \
                        speed reduction')

    # Splice junction database
    parser.add_argument('--sjdbGTFfile',
                        help='path to GTF file with annotations')
    parser.add_argument('--sjdbGTFchrPrefix', default='-',
                        help='prefix for chromosome names in a GTF \
                        file (e.g. ’chr’ for using ENSMEBL \
                        annotations with UCSC geneomes)')
    parser.add_argument('--sjdbGTFfeatureExon', default='exon',
                        help='feature type in GTF file to be used as \
                        exons for building transcripts')
    parser.add_argument('--sjdbGTFtagExonParentTranscript',
                        default='transcript_id',
                        help='Exons are assigned to the transcripts \
                        using parent-child relationship defined by \
                        this option')
    parser.add_argument('--sjdbGTFtagExonParentGene',
                        default='gene_id',
                        help='tag name to be used as exons’ \
                        gene-parents (Default ”gene id” works \
                        formatGTF files)')
    parser.add_argument('--sjdbFileChrStartEnd', default='-',
                        help='path to annotations formatted as a list \
                        of splice junctions coordinates in a text \
                        file')
    parser.add_argument('--sjdbOverhang', type=int, default=100,
                        help='length of the donor/acceptor sequence \
                        on each side of the junctions, \
                        ideally = (mate length - 1). \
                        If =0, splice junction database is not used')


def buildSTARcmd(args):
    baseCmd = 'STAR'
    formattedCmd = baseCmd
    for flag, arg, in vars(args).items():
        if arg is not None:
            if not isinstance(arg, (list, tuple)):
                formattedCmd += ' --{0} {1}'.format(flag, arg)
            else:
                formattedCmd += ' --{0}'.format(flag)
                for item in arg:
                    formattedCmd += ' {0}'.format(item)
    return formattedCmd


def execute(cmd):
    logging.info('RUNNING: %s' % (cmd))
    print 'RUNNING...\n\n', cmd
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if len(stderr):
        print stderr
    return p.returncode


def main():
    parser = collectArgs()
    args = parser.parse_args()
    cmd = buildSTARcmd(args)
    execute(cmd)


if __name__ == '__main__':
    main()
