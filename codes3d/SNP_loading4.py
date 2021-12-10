from codes3d import *
from snps import *
import features

    
if __name__ == '__main__':
    args = parse_args() #V1
    c = CODES3D(args.config) #V2   #args.config defaults to: default=os.path.join(os.path.dirname(__file__), '../docs/codes3d.conf')
    commons_db = create_engine(c.commons_db_url, echo=False, poolclass=NullPool)
    
    if not os.path.isdir(args.output_dir): #V4   #If the specified output path is not a directory (aka it doesn't exist) then CREATE that specified directory
        os.makedirs(args.output_dir)
    logger = Logger(logfile=os.path.join(  #V18 #This initialises the "Logger" class as "logger"
        args.output_dir, 'snps_loading.log'))
    start_time = time.time()#V20
    tissues = parse_tissues(        #V16 #tissues is specified here!
        args.tissues, args.match_tissues, args.eqtl_project, commons_db)
    
    hic_df = parse_hic( #V5                   #Returns: hic_df: a dataframe columns(library, enzyme, rep_count)
        args.match_tissues,
        args.include_cell_lines,
        args.exclude_cell_lines,
        args.restriction_enzymes,
        commons_db)
        
    log_settings(args, logger) #V19 #This fn uses the Logger class to write the initial settings on the console
    
    eqtl_project = tissues['project'].tolist()[0]#V17  #This code take the rows under 'project' column -> makes it into a list -> take the first project out (aka index 0)
    
    eqtl_project_db = create_engine( #V6       
        c.eqtl_db_url.format(eqtl_project.lower()), echo=False, poolclass=NullPool) #self.eqtl_db_url = 'postgresql://{}:{}@{}:{}/{}'.format(self.user, self.password, self.host, self.port, 'eqtls_{}')
    snp_df = pd.DataFrame() #V7

    
    if args.snp_input: #V10
    
        gene_df = [] #V12         #Okay... so I guess this is a local assignment for gene_df...
        gene_info_df = None #V13 #Just need to be specified to satisfy the requirement of the get_snp() function!! under the DEFAULT setting
        
        #The function below return snp_df... but somehow with restriction enzyme info??!?!?
        snp_df = get_snp( #V14                        #Specified on line 700 as snp_df = pd.DataFrame(), so snp_df is a pandas dataframe
            args.snp_input,                           #This is the main input aka the file containing dbSNP ids
            gene_info_df,                             #If args.snp_within_gene is not specified, then gene_info_df = None as specified on line 746
            hic_df,
            args.output_dir,
            eqtl_project_db,   #For tommorow: this eqtl_project_db require parse_tissues (in a convoluted way.. somehow...), which we have analysed yesterday... now investigate this crap!!!
            c.rs_merge_arch_fp,                       #self.rs_merge_arch_fp = os.path.join(os.path.dirname(__file__), config.get("Defaults", "RS_MERGE_ARCH"))
            logger,
            args.suppress_intermediate_files)         #'--suppress-intermediate-files', action='store_true', default=False,
                                   
        if not args.suppress_intermediate_files: #V15   #defined in argparse above as -> '--suppress-intermediate-files', action='store_true', default=False
            snp_df.to_csv(os.path.join(                 #Thus if surpress intermediate file is not called, this this conditional becomes not False aka TRUE and thus executed!
                args.output_dir, 'snps.txt'), sep='\t', index=False)
        if args.features:
            features.feature_intersect(snp_df, args.output_dir, logger)
