#!/usr/bin/env python

import sys

DEPTHS = [200, 100, 50, 25, 10, 0]
EXON_BUFF = 10

def createDocDict(handle):

    print("Creating depth of coverage dictionary.")

    doc_dict = {}
    with handle as doc:
        i = 0
        next(doc)
        for line in doc:
            if i % 1000000 == 0:
                print(i)
            line = line.rstrip('\n').split('\t')
            doc_dict[line[0]] = int(line[1])
            i += 1
 
    return doc_dict


def createExonDict(handle, doc_dict_0, doc_dict_30):
    
    """
    Exons file looks like:
    1       11872   12227   OR4F5 ENSE00002234632
    OFS = '\t'
    """

    print("Creating exon dictionary.")

    exon_dict = {}
    with handle as exons:
        i = 0
        next(exons)
        for line in exons:
            if i % 100000 == 0:
                print(i)
            line = line.rstrip('\n').split('\t')

            chrom = line[0]
            start = line[1]
            stop = line[2]
            hgnc = line[3]
            ense = line[4]

            read_count_0 = 0.0
            read_count_30 = 0.0
            this_depth = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]            

            for coord in range(int(start)-EXON_BUFF+1, int(stop)+1+EXON_BUFF):
                comp_coord = line[0] + ':' + str(coord)
                if comp_coord in doc_dict_0 and comp_coord in doc_dict_30:
                    read_count_0 += doc_dict_0[comp_coord]
                    read_count_30 += doc_dict_30[comp_coord]
                    this_depth = calcDepth(doc_dict_30[comp_coord], this_depth)
            total_range = int(stop) - int(start) + 1 + (2*EXON_BUFF)

            uniq_key = ense + '_' + hgnc
            exon_dict[uniq_key] = [hgnc, read_count_0, read_count_30, total_range]
            exon_dict[uniq_key].extend(this_depth)
            exon_dict[uniq_key].extend([chrom, int(start), int(stop)])
            exon_dict[uniq_key].extend([chrom, start, stop])
            i += 1

    return exon_dict


def writeExonQC(handle, exon_dict):

    """
    exon_dict now looks like:
    {ENSE_HGNC = [HGNC, read_count_0, read_count_30, total_range, D200, D100, D50, D25, D10, D0, chrom, start, stop]}
    """

    print("Writing Exon QC.")

    handle.write("Chromosome\tStart\tStop\tGene\tExon\tAvgD\tQ30\tD200\tD100\tD50\tD25\tD10\n")

    for key in exon_dict:

        read_count_0 = exon_dict[key][1]
        read_count_30 = exon_dict[key][2]
        total_range = exon_dict[key][3]
        hgnc = exon_dict[key][0]

        if exon_dict[key][1] != 0:
            q30 = percent(read_count_30, read_count_0)
        else:
            q30 = 0.0
        avgd = "%.1f" % (read_count_30/total_range)
        d200 = percent(exon_dict[key][4], total_range)
        d100 = percent(exon_dict[key][5], total_range)
        d50 = percent(exon_dict[key][6], total_range)
        d25 = percent(exon_dict[key][7], total_range)
        d10 = percent(exon_dict[key][8], total_range)

        chrom = exon_dict[key][10]
        start = exon_dict[key][11]
        stop = exon_dict[key][12]
        
        handle.write(chrom + '\t' + str(start) + '\t' + str(stop) + '\t' + hgnc + '\t' + key.split('_')[0] + '\t' + str(avgd) + '\t' + str(q30) + '\t' + str(d200) + '\t' + str(d100) + '\t' + str(d50) + '\t' + str(d25) + '\t' + str(d10) + '\n')


def writeIntervalQC(handle, int_dict):

    """
    int_dict now looks like:
    {chrom:start-stop = [read_count_0, read_count_30, total_range, D200, D100, D50, D25, D10, D0]}
    """

    print("Writing Interval QC.")

    handle.write("Chromosome\tStart\tStop\tAvgD\tQ30\tD200\tD100\tD50\tD25\tD10\n")

    sample_total = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    for key in int_dict:
        
        read_count_0 = int_dict[key][0]
        read_count_30 = int_dict[key][1]
        total_range = int_dict[key][2]

        chrom = key.split(':')[0]
        start = key.split(':')[1].split('-')[0]
        stop = key.split(':')[1].split('-')[1].rstrip('\n')
        
        if int_dict[key][1] != 0:
            q30 = percent(read_count_30, read_count_0)
        else:
            q30 = 0.0

        avgd = "%.1f" % (read_count_30/total_range)
        d200 = percent(int_dict[key][3], total_range)
        d100 = percent(int_dict[key][4], total_range)
        d50 = percent(int_dict[key][5], total_range)
        d25 = percent(int_dict[key][6], total_range)
        d10 = percent(int_dict[key][7], total_range)

        handle.write(chrom + '\t' + start + '\t' + stop + '\t' + str(avgd) + '\t' + str(q30) + '\t' + str(d200) + '\t' + str(d100) + '\t' + str(d50) + '\t' + str(d25) + '\t' + str(d10) + '\n')

        sample_total = [x + y for x, y in zip(sample_total, int_dict[key])]
        
    if sample_total[0] != 0:
        q30 = percent(sample_total[1], sample_total[0])
    else:
        q30 = 0.0

    avgd = "%.1f" % (sample_total[1]/sample_total[2])
    d200 = percent(sample_total[3], sample_total[2])
    d100 = percent(sample_total[4], sample_total[2])
    d50 = percent(sample_total[5], sample_total[2])
    d25 = percent(sample_total[6], sample_total[2])
    d10 = percent(sample_total[7], sample_total[2])

    ### Write the final totals row, for sample level metrics.

    handle.write('TOTAL\t\t\t' + str(avgd) + '\t' + str(q30) + '\t' + str(d200) + '\t' + str(d100) + '\t' + str(d50) + '\t' + str(d25) + '\t' + str(d10) + '\n')


def percent(num1, num2):
    
    return "%.1f" % (num1*100/num2)


def calcDepth(depth, total_depths):
    
    for value in DEPTHS:
        if depth >= value:
            total_depths[DEPTHS.index(value)] += 1
    return total_depths


def createIntervalDict(handle, doc_dict_0, doc_dict_30):

    print("Creating Probe dictionary.")
    
    int_dict = {}
    with handle as intervals:
        i = 0
        for line in intervals:
            if i % 1000000 == 0:
                print(i)
            nline = line.rstrip('\n').split('-')
            read_count_0 = 0.0
            read_count_30 = 0.0
            this_depth = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            start = int(nline[0].split(':')[1])
            stop = int(nline[1])
            chrom = (nline[0].split(':')[0])

            if len(nline) == 2:
                for coord in range(start, stop+1):
                    comp_coord = chrom + ':' + str(coord)
                    if comp_coord in doc_dict_0 and comp_coord in doc_dict_30:
                        read_count_0 += doc_dict_0[comp_coord]
                        read_count_30 += doc_dict_30[comp_coord]
                        this_depth = calcDepth(doc_dict_30[comp_coord], this_depth)
                total_range = stop - start + 1
                int_dict[line.rstrip('\n')] = [read_count_0, read_count_30, total_range]
                int_dict[line.rstrip('\n')].extend(this_depth)
                i += 1

    return int_dict


def createGeneDict(ref_dict, coords, doc_dict_0, doc_dict_30):


    """
    Exons file looks like:
    1       11872   12227   OR4F5 ENSE00002234632
    OFS = '\t'
    """
    print("Creating gene dictionary.")

    gene_dict = {}

    i = 0
    for hgnc in coords:
        if i % 10000 == 0:
            print(i)

        if hgnc in ref_dict:
            chrom = ref_dict[hgnc][0]
            start = ref_dict[hgnc][1]
            stop = ref_dict[hgnc][2]


            if hgnc not in gene_dict:
                gene_dict[hgnc] = []
                this_depth = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]            
                read_count_0 = 0.0
                read_count_30 = 0.0

                total_range = 0
                for ival in coords[hgnc]:
                    for coord in range(ival[0], ival[1]):
                        comp_coord = chrom + ':' + str(coord)
                        if comp_coord in doc_dict_0 and comp_coord in doc_dict_30:
                            read_count_0 += doc_dict_0[comp_coord]
                            read_count_30 += doc_dict_30[comp_coord]
                            this_depth = calcDepth(doc_dict_30[comp_coord], this_depth)
                        total_range += 1
        
            gene_dict[hgnc] = [read_count_0, read_count_30, total_range]
            gene_dict[hgnc].extend(this_depth)
            gene_dict[hgnc].extend([chrom, start, stop])
        
        i += 1

    return gene_dict


def createRefGenes(handle):
    
    """
    Create reference gene dictionary, matching genes to coordinates.
    Ref Genes looks like:
    OR4F5	ENSG00000186092	ENST00000335137	NM_001005484	CCDS30547	1	69091	70008	+
    """

    ref_dict = {}

    with handle as genes:
        for gene in genes:
            gene = gene.rstrip('\n').split('\t')
            ref_dict[gene[0]] = [gene[5], gene[6], gene[7]]

    return ref_dict


def createCoordDict(exon_dict):

    """
    exon_dict now looks like:
    {ENSE_HGNC = [HGNC, read_count_0, read_count_30, total_range, D200, D100, D50, D25, D10, D0, chrom, start, stop]}
    """

    print("Creating coordinate dictinoary.")

    coord_dict = {}

    for exon in exon_dict:
        hgnc = exon_dict[exon][0]
        start = exon_dict[exon][11] - EXON_BUFF
        stop = exon_dict[exon][12] + EXON_BUFF

        if hgnc not in coord_dict:
            coord_dict[hgnc] = [(start,stop)]
        else:
            coord_dict[hgnc].append((start,stop))
     

    return coord_dict


def mergeCoordDict(coords):

    print("Merging coordinate dictionary intervals.")
    
    for key in coords:
        new_ivals = []
        for ival in sorted(coords[key]):
            if new_ivals == []:
                new_ivals.append(ival)
                continue
            else:
                if ival[0] > new_ivals[-1][1]+1:
                    new_ivals.append(ival)
                elif ival[1] > new_ivals[-1][1]:
                    new_ivals[-1] = (new_ivals[-1][0], ival[1])

        coords[key] = new_ivals

    return coords


def writeGeneQC(handle, gene_dict):

    """
    gene_dict now looks like:
    {HGNC = [read_count_0, read_count_30, total_range, D200, D100, D50, D25, D10, D0, chrom, start, stop]}
    """

    print("Writing gene QC.")

    handle.write("Chromosome\tStart\tStop\tGene\tAvgD\tQ30\tD200\tD100\tD50\tD25\tD10\n")

    for key in gene_dict:

        read_count_0 = gene_dict[key][0]
        read_count_30 = gene_dict[key][1]
        total_range = gene_dict[key][2]
        
        chrom = gene_dict[key][9]
        start = gene_dict[key][10]
        stop = gene_dict[key][11]

        if gene_dict[key][0] != 0:
            q30 = percent(read_count_30, read_count_0)
        else:
            q30 = 0.0
        avgd = "%.2f" % (read_count_30/total_range)
        d200 = percent(gene_dict[key][3], total_range)
        d100 = percent(gene_dict[key][4], total_range)
        d50 = percent(gene_dict[key][5], total_range)
        d25 = percent(gene_dict[key][6], total_range)
        d10 = percent(gene_dict[key][7], total_range)

        handle.write(chrom + '\t' + start + '\t' + stop + '\t' + key + '\t' + str(avgd) + '\t' + str(q30) + '\t' + str(d200) + '\t' + str(d100) + '\t' + str(d50) + '\t' + str(d25) + '\t' + str(d10) + '\n')


def main():

    all_exons = open(sys.argv[1], 'rU')
    doc_input_0 = open(sys.argv[2], 'rU')
    doc_input_30 = open(sys.argv[3], 'rU')

    doc_dict_0 = createDocDict(doc_input_0)
    doc_dict_30 = createDocDict(doc_input_30)
    
    exon_dict = createExonDict(all_exons, doc_dict_0, doc_dict_30)
    exonqc = open(sys.argv[4], 'w')
    writeExonQC(exonqc, exon_dict)
    exonqc.close()

    all_intervals = open(sys.argv[5], 'rU')
    interval_dict = createIntervalDict(all_intervals, doc_dict_0, doc_dict_30)
    intervalqc = open(sys.argv[6], 'w')
    writeIntervalQC(intervalqc, interval_dict)
    intervalqc.close()
    
    
    ref_dict = createRefGenes(open(sys.argv[8], 'rU'))
    coord_dict = createCoordDict(exon_dict)
    coord_dict = mergeCoordDict(coord_dict)
    geneqc = open(sys.argv[7], 'w')
    gene_dict = createGeneDict(ref_dict, coord_dict, doc_dict_0, doc_dict_30)
    writeGeneQC(geneqc, gene_dict)
    geneqc.close()


if __name__ == "__main__":
    main()


# def main2():

#     doc_input_0 = open(sys.argv[1], 'rU')
#     doc_input_30 = open(sys.argv[2], 'rU')
#     doc_dict_0 = createDocDict(doc_input_0)
#     doc_dict_30 = createDocDict(doc_input_30)

#     all_intervals = open(sys.argv[3], 'rU')
#     interval_dict = createIntervalDict(all_intervals, doc_dict_0, doc_dict_30)
# #    for key in interval_dict:
# #        print(key)
# #        print(interval_dict[key])

#     intervalqc = open(sys.argv[4], 'w')
#     writeIntervalQC(intervalqc, interval_dict)
#     intervalqc.close()

# main2()
