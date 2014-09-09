import numpy as np
import pysam
import sys
import glob
import optparse

from core import gasic
from core import tools


class CorrectAbundances(object):
    
    @staticmethod
    def similarityCorrection(samFiles, smatFile, nBootstrap, npar_boot):
        """
        Perform similarity correction step. The similarity matrix and mapping
        results must be available.
        
        Args:
        samFiles -- list of SAM files (query reads mapped to each reference).
        smatFile -- mapping information for similarity matrix with same ordering as simSamFile list.
        nBootstrap -- number of bootstrap samples, use 1 to disable bootstrapping.
        npar_boot -- number of bootstrap samples to processes in parallel.
        
        OUTPUT:
        total:             total number of reads in the dataset
        num_reads:         number of reads mapped to each genome (array)
        corr:              abundance of each genome after similarity correction
        err:               estimated standard error
        p:                 p-value for the confidence, that the true abundance is above some threshold
        """

        # find out the total number of reads for first sam file
        total = len( [1 for read in pysam.Samfile(samFiles[0], "r")] )
        sys.stderr.write("...found {} reads\n".format(total))
        
        # initialize some arrays
        #   mapping information; mapped[i,j]=1 if read j was successfully mapped to i.
        mapped = np.zeros( (len(samFiles), total) )
        
        # total number of successfully mapped reads per reference
        num_reads = np.zeros( (len(samFiles),) )
        
        # analyze the SAM files
        for n_ind,samFile in enumerate(samFiles):
            sys.stderr.write("...analyzing SAM-File {} of {}\n".format(n_ind+1, len(samFiles)))
            # samfile filename
            sf = pysam.Samfile(samFile, "r")
            
            # go through reads in samfile and check if it was successfully mapped
            mapped[n_ind,:] = np.array([int(not rd.is_unmapped) for rd in sf])
            num_reads[n_ind] = sum(mapped[n_ind,:])
            
        # run similarity correction step
        smat = np.load(smatFile)
        #p,corr,var = gasic.bootstrap(mapped, smat, nBootstrap)
        p,corr,var = gasic.bootstrap_par(mapped, smat, nBootstrap, nprocs=npar_boot)
        err = np.sqrt(var)
        return dict(total=total,num_reads=num_reads,corr=corr,err=err,p=p)
            


    @staticmethod
    def unique(names, sam_pattern):
	""" Determine the number of unique reads for every species based on the read names.
	INPUT:
	names:             array of genome names
	sam_pattern:       pattern pointing to the filenames of the SAM-files to analyze

	OUTPUT:
	unique:            number of unique reads per species.
	"""
	# one set for the names of mapped reads for each species
	mapped_read_names = [set() for nm in names]
	
	for n,nm in enumerate(names):
            # parse the samfile
            sf = pysam.Samfile(sam_pattern%nm, "r")
            for read in sf:
                # add the hash of read name to the set, if read was mapped
                if not read.is_unmapped:
                    mapped_read_names[n].add(hash(read.qname))

	unique_read_names = [set() for nm in names]
	for n in range(len(names)):
            others = set()
            for m in range(len(names)):
                if n!=m:
                    others |= mapped_read_names[m]
            unique_read_names[n] = mapped_read_names[n] - others

	return np.array([len(unq) for unq in unique_read_names])


class CorrAbundRes(object):
    """Class for results of correctAbundances() function"""

    def __init__(self, total, num_reads, corr, err, p):    
        self.total=total
        self.num_reads=num_reads
        self.corr=corr
        self.err=err
        self.p=p

    def write(self, samFiles, outFileName):
        """results into tab separated file."""
        outfh = open(outFileName,'w')
        outfh.write('\t'.join(['genome name', 'total', 'mapped reads', 'estimated reads', 'estimated error', 'p-value']))
        for n_ind,samFile in enumerate(samFiles):
            # Name, mapped reads, corrected reads, estimated error, p-value
            out = "{name}\t{total}\t{mapped}\t{corr}\t{error}\t{pval}\n"
            outfh.write(out.format(name=nm,
                                   total=total,
                                   mapped=num_reads[n_ind],
                                   corr=corr[n_ind]*total,
                                   error=err[n_ind]*total,
                                   pval=p[n_ind]))
            outfh.close()
            print "...wrote results to {}".format(outFileName)
    
        

#--- main ---#    
    # numArgs = len(args)
    # if numArgs == 1:
    #     # read the Names file
    #     names_file = args[0]
    #     names = tools.read_names(names_file)

    #     # load the similarity matrix
    #     smat = np.load(opt.smat)

    #     # start similarity correction
    #     total,num_reads,corr,err,p = similarity_correction(names, smat, opt.sam, opt.boot)

    #     # write results into tab separated file.
    #     ofile = open(opt.out,'w')
    #     ofile.write("#genome name\tmapped reads\testimated reads\testimated error\tp-value\n")
    #     for n_ind,nm in enumerate(names):
    #         # Name, mapped reads, corrected reads, estimated error, p-value
    #         out = "{name}\t{mapped}\t{corr}\t{error}\t{pval}\n"
    #         ofile.write(out.format(name=nm,mapped=num_reads[n_ind],corr=corr[n_ind]*total,error=err[n_ind]*total,pval=p[n_ind]))
    #     ofile.close()
    #     print ". wrote results to", opt.out

    # else:
    #     parser.print_help()
    #     sys.exit(1)
