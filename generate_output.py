def preprocess_outputs():
    ''' return preprocessing outputs'''
    outputs = expand("{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam.bai", libname = libnames, sample_label = rbps
    )+expand("{libname}/bw/COV/{sample_label}.{strand}.bw", libname = libnames, sample_label = rbps, strand = ['pos', 'neg']
    )+['QC/fastQC_basic_summary.csv',
        'QC/fastQC_passfail.csv',
        'QC/cutadapt_stat.csv',
        "QC/mapping_stats.csv",
        "QC/dup_level.csv",
        'QC/demux_read_count.txt',
        "QC/summary.csv"
        ]+expand("counts/genome/vectors/{libname}.{sample_label}.counts",
        libname = libnames, sample_label = rbps
    )+expand("QC/read_count/{libname}.{metric}.csv", libname = libnames, metric = ['region', 'genetype', 'cosine_similarity']
    )+expand("counts/repeats/tables/{repeat_type}/{experiment}.{sample_label}.tsv.gz", 
    experiment = experiments, sample_label = rbps, repeat_type = ['name', 'class', 'family']
    )+expand("counts/repeats/megatables/{repeat_type}/{libname}.tsv.gz", libname = libnames, repeat_type = ['name', 'class', 'family'])
    return outputs

def beta_binom_mixture_outputs():
    ''' generate output for beta-binomial mixture'''
    outputs = []
    if not singleplex:
        # complementary control
        outputs+=expand("beta-mixture_CC/{libname}.{sample_label}.enriched_windows.tsv",
        libname = libnames,
        sample_label = list(set(rbps)-set(config['AS_INPUT'])),
        )
        
        if config['FINEMAPPING']:
            outputs += expand("beta-mixture_CC/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.{strand}.bw",
            libname = libnames,
            sample_label = list(set(rbps)-set(config['AS_INPUT'])),
            signal_type = ['CITS', 'COV'],
            strand = ['pos', 'neg']
            )+expand("beta-mixture_CC/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html",
            libname = libnames,
            sample_label = config['RBP_TO_RUN_MOTIF'],
            signal_type = ['CITS', 'COV'])
        
        # internal control
        outputs+=expand("beta-mixture/{bg_sample_label}/{libname}.{clip_sample_label}.enriched_windows.tsv",
        libname = libnames,
        clip_sample_label = list(set(rbps)-set(config['AS_INPUT'])),
        bg_sample_label = config['AS_INPUT'] if config['AS_INPUT'] else [] #TODO: IgG secondary analysis
        )
        
        if config['FINEMAPPING']:
            outputs += expand("beta-mixture/{bg_sample_label}/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.{strand}.bw",
            libname = libnames,
            sample_label = list(set(rbps)-set(config['AS_INPUT'])),
            bg_sample_label = config['AS_INPUT'] if config['AS_INPUT'] else [],
            signal_type = ['CITS', 'COV'],
            strand = ['pos', 'neg']
            )+expand("beta-mixture/{bg_sample_label}/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html",
            libname = libnames,
            bg_sample_label = config['AS_INPUT'] if config['AS_INPUT'] else [],
            sample_label = config['RBP_TO_RUN_MOTIF'],
            signal_type = ['CITS', 'COV'])

    if external_normalization:
        outputs += expand("beta-mixture_external/{external_label}/{libname}.{clip_sample_label}.enriched_windows.tsv",
        libname = libnames, clip_sample_label = rbps,
        external_label = list(external_normalization.keys())
        )
        if config['FINEMAPPING']:
            outputs += expand("beta-mixture_external/{external_label}/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.{strand}.bw",
            libname = libnames,
            sample_label = list(set(rbps)-set(config['AS_INPUT'])),
            signal_type = ['CITS', 'COV'],
            strand = ['pos', 'neg'],
            external_label = list(external_normalization.keys())
            )+expand("beta-mixture_external/{external_label}/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html",
            libname = libnames,
            sample_label = config['RBP_TO_RUN_MOTIF'],
            external_label = list(external_normalization.keys()),
            signal_type = ['CITS', 'COV'])
        
    # complementarry control bigwigs
    if len(set(rbps)-set(config['AS_INPUT']))>2 and len(rbps)>1:
        outputs+=expand("{libname}/bw_bg/COV/{sample_label}.{strand}.bw", libname = libnames, sample_label = rbps, strand = ['pos', 'neg']
        )
    return outputs

def DMM_outputs():
    ''' generate output from Dirichlet Multinomial mixture!'''
    outputs = []
    if not singleplex:
        outputs+=expand("DMM/{libname}.{sample_label}.enriched_windows.tsv", 
            sample_label = list(set(rbps)-set(config['AS_INPUT'])), 
            libname = libnames
        )+expand("DMM_repeat/{repeat_type}/{libname}.{sample_label}.enriched_windows.tsv", 
            sample_label = list(set(rbps)-set(config['AS_INPUT'])), 
            repeat_type = ['name'],
            libname = libnames
        )+expand("DMM_repeat/{repeat_type}/{libname}.megaoutputs.tsv",
            libname = libnames,
            repeat_type = ['name']
        )+expand('mask/{libname}.genome_mask.csv',
        libname = libnames,
        )+expand('mask/{libname}.repeat_mask.csv',
        libname = libnames,
        )
        
        if config['FINEMAPPING']:
            outputs+=expand("DMM/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html", libname = libnames,
            sample_label = config['RBP_TO_RUN_MOTIF'],
            signal_type = ['CITS', 'COV']
            )+expand("DMM/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz",
            libname = libnames,
            sample_label = list(set(rbps)-set(config['AS_INPUT'])),
            signal_type = ['CITS', 'COV']
            )
        
    return outputs

def get_output(DMM, BBM):
    output = preprocess_outputs()
    
    if DMM:     
        output += DMM_outputs() 
    if BBM:
        output += beta_binom_mixture_outputs()
    return output
