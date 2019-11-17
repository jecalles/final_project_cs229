import Bio
from tqdm import tqdm
import dill as pickle

class dataparse:
    ''' A data parsing class for our ML project

    Attributes
        infile (str): path to genbank file for parsing
        outfile (str): path to write parsed data
        feature_type (list-like<str>): which genes should be extracted?
        promoter_len (int): number of bases to consider for promoter
        genome (list<Bio.SeqRecord>): list of chromosomes in organism
        genes (list<Bio.SeqFeature>): list of genes matching target annotation
        promoters (Bio.SeqRecord): promoter for given gene
        transcription_units (list<tuple(Bio.SeqFeature, Bio.SeqRecord)>):  list of gene, promoter pairs
    '''
    def __init__(self, infile=None, outfile=None, get_genes=False, get_promoters=False, parse=False, verbose=False, feature_type='CDS', promoter_len=1000):
        ''' initialize dataparser object. Optionally load a dataset
        '''
        self.infile = infile
        self.outfile = outfile
        self.feature_type = feature_type
        self.promoter_len = promoter_len
        self.genome = None
        self.genes = None
        self.promoters = None
        self.transcription_units = None

        # optionally load data
        if infile != None:
            self.load(get_genes, get_promoters, parse, infile=infile, feature_type=feature_type, inplace=True)

        # optionally write data to outfile
        if outfile != None:
            self.write(outfile)

    @classmethod
    def from_pickle(cls, pickle_file):
        ''' classmethod used to initialize dataparse object from a dill file. Uses dill package
        '''
        with open(pickle_file, 'rb') as fhandle:
            obj = pickle.load(fhandle)
        return obj

    def load(self, get_genes=False, get_promoters=False, parse=False, verbose=False, **kwargs):
        '''load genbank files for parsing. Optionally parse data while you're at it

        Args
            infile (str): path to genbank file for parsing
            get_genes (bool=False): optionally get list of genes
            get_promoters (bool=False): optionally get list of promoters
            get_pairs (bool=False): optionally get pairs of genes and promoters
            inplace (bool=False): optionally save all loaded content to instance

        Return
            genome (list<Bio.SeqRecord>): list of chromosomes in organism
            transcription_units (list<tuple(Bio.SeqFeature, Bio.SeqRecord)>):  list of gene, promoter pairs
        '''
        # parse optional arguments
        infile = kwargs.get('infile', self.infile)
        feature_type = kwargs.get('feature_type', self.feature_type)
        # get chromosomes in genbank file
        self.genome = [record for record in Bio.SeqIO.parse(infile, 'genbank')]
        # optionally parse data
        if parse:
            # perform full tilt analysis
            #   - genes
            #   - promoters
            #   - gene/promoter pairs
            self.get_transcription_units(verbose=verbose, inplace=True)
        elif get_promoters:
            # get genes and promoters as pairs
            self.get_transcription_units(verbose=verbose, inplace=True)
        elif get_genes:
            # get all genes in the genome
            self.get_genes(inplace=True)

    def get_transcription_units(self, verbose=True, inplace=False, **kwargs):
        ''' method used to get all promoter/gene pairs at once
        '''
        # parse optional arguments
        genome = kwargs.get('genome', self.genome)
        feature_type = kwargs.get('feature_type', self.feature_type)
        promoter_len = kwargs.get('promoter_len', self.promoter_len)
        promoters = []
        # loop over chromosomes in genome and extract genes
        genome_iterator = tqdm(genome, desc='dataparse.get_transcription_units: parsing chromosomes', leave=True) if verbose else genome
        for chromosome in genome_iterator:
            genes = self.get_genes(genome=chromosome, feature_type=feature_type, verbose=False)
            # extract promoter for each gene
            for gene in genes:
                promoter = self.get_promoter(gene, chromosome)
                promoters.append(promoter)
        # zip promoters and genes together
        transcription_units = zip(promoters, genes)
        # optionally update instance inplace
        if inplace:
            self.genome = genome
            self.feature_type = feature_type
            self.promoter_len = promoter_len

    def get_genes(self, inplace=False, **kwargs):
        ''' method used to get genes in a Bio.Seqrecord that match a given feature type

        Args
            genome (list<Bio.SeqRecord>): chromosomes to parse
            feature_type (list-like<str>): which genes should be extracted?
            inplace (bool=False): optionally save all loaded content to instance

        Return
            genes (list<Bio.SeqFeature>): genes that match desired type
        '''
        # parse optional arguments
        genome = kwargs.get('genome', self.genome)
        feature_type = kwargs.get('feature_type', self.feature_type)
        # parse chromosome if a single entry
        if type(genome) == Bio.SeqRecord.SeqRecord:
            genes = [feature for feature in genome.features if feature.type in feature_type]
        else:
            for chromosome in genome:
                genes = []
                genes.append(
                    [feature for feature in chromosome.features if feature.type in feature_type]
                )
        # optionally update instance inplace
        if inplace:
            self.genome = genome
            self.genes = genes
            self.feature_type = feature_type
        return genes

    def get_promoter(self, gene=None, chromosome=None, **kwargs):
        ''' method used to extract promoters from SeqRecord. Can optionally parse the entire genome in one go

        Args
            gene (Bio.SeqFeature): gene to parse
            chromosome (Bio.SeqRecord): chromosome from which gene comes
            promoter_len (int): number of bases to consider for promoter

        Return
            promoter (Bio.SeqRecord): promoter for a given gene
        '''
        # parse optional arguments
        promoter_len = kwargs.get('promoter_len', self.promoter_len)
        # get position range for promoter, modulated by strand
        if gene.strand == 1:
            # positive sense
            end = gene.location.start.position
            start = end - promoter_len
        elif gene.strand == -1:
            # negative sense
            start = gene.location.end.position
            end = start + promoter_len
        else:
            # nonsense ;)
            raise ValueError('given gene.strand is not in {1,-1}')

        # extract promoter
        promoter = chromosome[start:end]
        if gene.strand == -1:
            promoter = promoter.reverse_complement()
        return promoter

    def write(self, outfile=None):
        ''' method that writes parsed data to pickle file

        Args
            outfile (str): pickle file to write to
        '''
        # load outfile if not explicitly given
        if outfile==None:
            outfile=self.outfile
        # pickle instance using dill package
        with open(outfile, 'wb') as fhandle:
            pickle.dump(self, fhandle)
