import pandas as pd
import pybedtools
import os


def feature_intersect(inputs, output_dir, logger):   #inputs is specified as snp_df in the execution within SNP_loading4v2
	logger.write('Looking for epigenomic features...')

	snp_df = inputs
    
    #create a variable called snp_df, which is the .bed version of snps.txt, and save that snp_bed into output
	snp_bed = pd.DataFrame()
	snp_bed['chrom'] = snp_df['chrom']
	snp_bed['start'] = snp_df['locus'].sub(1)
	snp_bed['end'] = snp_df['locus']
	snp_bed.to_csv(os.path.join(output_dir, 'snps_bed.bed'), sep = '\t', index = False, header = False)
    
    #Now read back the snp_bed.bed file into a BedTool object!! so that it can be manipulated using pybedtools
	snp_bed = pybedtools.BedTool(os.path.join(output_dir, 'snps_bed.bed'))
    
    #Intersecting with H3K27ac
	logger.write('---> Looking for intersections with H3K27ac peaks...')
    
    #1. read all of the H3K27ac reference files into BedTool objects stored within dictionary keys
	dictionary = {}
	for filename in os.listdir(os.path.join(os.path.dirname(__file__),'../reference_file2/')):
		if "H3K27ac" in filename:
			f = os.path.join(os.path.join(os.path.dirname(__file__),'../reference_file2/'), filename)
			dictionary[filename] = pybedtools.BedTool(f)
	
	list = dictionary.values()
    
    #2. Intersect the 'a' file(snp_bed) with all of the 'b' files(reference file: H3K27ac broadpeak), and save the output .bed file to the user-specified args.output_dir!
	snp_bed.intersect(list, u = True).saveas(os.path.join(output_dir, 'snps_that_intersects_H3K27ac.bed'))
    
    #3. Next, read back the output file into a Data Frame!
	snps_that_intersects_H3K27ac = pd.read_csv(os.path.join(output_dir, 'snps_that_intersects_H3K27ac.bed'), sep = '\t', header = None)
    
    #4. Make a list version of snps_that_intersects_H3K27ac to create a list of SNPs as a reference for the query operation
	snps_that_intersects_H3K27ac_list = snps_that_intersects_H3K27ac[2].drop_duplicates().tolist()
    
    #5. Algorithm to query which SNPs inside snp_df are also SNPs that intersect H3K27ac, and then record the result into a new column within snp_df
	for index, row in snp_df.iterrows():
		for i in snps_that_intersects_H3K27ac_list:
			if i == row['locus']:
				snp_df.loc[index, 'H3K27ac'] = True
				break
			else:
				snp_df.loc[index, 'H3K27ac'] = False
    
    #Intersecting with H3K4me1
	logger.write('---> Looking for intersections with H3K4me1 peaks...')
	dictionary = {}
	for filename in os.listdir(os.path.join(os.path.dirname(__file__),'../reference_file2/')):
		if "H3K4me1" in filename:
			f = os.path.join(os.path.join(os.path.dirname(__file__),'../reference_file2/'), filename)
			dictionary[filename] = pybedtools.BedTool(f)
	list = dictionary.values()

	snp_bed.intersect(list, u=True).saveas(os.path.join(output_dir, 'snps_that_intersects_H3K4me1.bed'))
	snps_that_intersects_H3K4me1 = pd.read_csv(os.path.join(output_dir, 'snps_that_intersects_H3K4me1.bed'), sep = '\t', header = None)
	snps_that_intersects_H3K4me1_list = snps_that_intersects_H3K4me1[2].drop_duplicates().tolist()
	for index, row in snp_df.iterrows():
		for i in snps_that_intersects_H3K4me1_list:
			if i == row['locus']:
				snp_df.loc[index, 'H3K4me1'] = True
				break
			else:
				snp_df.loc[index, 'H3K4me1'] = False


    #Intersecting with DNase hypersensitivity sites
	logger.write('---> Looking for intersections with DNase hypersensitive sites...')
	dictionary = {}
	for filename in os.listdir(os.path.join(os.path.dirname(__file__),'../reference_file2/')):
		if "DNase" in filename:
			f = os.path.join(os.path.join(os.path.dirname(__file__),'../reference_file2/'), filename)
			dictionary[filename] = pybedtools.BedTool(f)
	list = dictionary.values()

	snp_bed.intersect(list, u=True).saveas(os.path.join(output_dir, 'snps_that_intersects_DNase.bed'))
	snps_that_intersects_DNase = pd.read_csv(os.path.join(output_dir, 'snps_that_intersects_DNase.bed'), sep = '\t', header = None)
	snps_that_intersects_DNase_list = snps_that_intersects_DNase[2].drop_duplicates().tolist()
	for index, row in snp_df.iterrows():
		for i in snps_that_intersects_DNase_list:
			if i == row['locus']:
				snp_df.loc[index, 'DNase_hypersensitive'] = True
				break
			else:
				snp_df.loc[index, 'DNase_hypersensitive'] = False
    
	return snp_df
