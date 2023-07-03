
import sys

#
# Short script that can be modified to convert bkgd output format
# to more standard BED format
#


# sys.stdout.write('track type=bedGraph name="B" description="Background selection statistic (B) from McVicker et al. 2009" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n')


DIR = "data/hg18/bkgd_files"

chrom_names = ["chr%d" % s for s in range(1, 23)]
chrom_names.append("chrX")

for chrom in chrom_names:
    filename = DIR + "/" + chrom + ".bkgd"
    f = open(filename, "rt")

    # use half-open BED-file coordinates (start chrom at 0)
    bed_end = 0
    
    for line in f:
        words = line.split()
        score = float(words[0])/1000
        length = int(words[1])

        bed_start = bed_end
        bed_end = bed_start + length

        ## bedgraph format:
        #sys.stdout.write("%s\t%d\t%d\t0.%d\n" %
        #                 (chrom, bed_start, bed_end, score))

        # bed format
        sys.stdout.write("%s\t%d\t%d\t.\t%.3f\n" %
                         (chrom, bed_start, bed_end, score))
        
        
    f.close()

        
        
        
        
