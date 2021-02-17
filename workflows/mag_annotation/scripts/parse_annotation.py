#!/usr/bin/env python3

import argparse
import os
import sys
import numpy as np
import pandas as pd
import pathlib
import joblib
from sklearn.ensemble import RandomForestClassifier
from Bio import SeqIO

def parse_arguments():
    parser = argparse.ArgumentParser(description='')
    required = parser.add_argument_group('Required arguments')
    
    required.add_argument('-c',
                          '--cluster_file',
                          dest='cluster_file',
                          metavar='',
                          required=True,
                          type=str,
                          help='Path to file with all samples')
    required.add_argument('-s',
                          '--sample_table',
                          dest='sample_table',
                          metavar='',
                          required=True,
                          type=str,
                          help='Path to file with all samples')

    required.add_argument('-a',
                          '--annotation_dir',
                          dest='annotation_dir',
                          metavar='',
                          required=True,
                          type=str,
                          help='Path to Annotation directory with subdirs containing annotations of each sample') 
    required.add_argument('-o',
                          '--outdir',
                          dest='outdir',
                          metavar='',
                          required=True,
                          type=str,
                          help='Path to directory out.')

                    

    optional = parser.add_argument_group('Optional arguments')
    
    optional.add_argument('-d',
                          '--split_delimn',
                          dest='split_delimn',
                          metavar='',
                          default='_',
                          required=False,
                          type=str,
                          help='Character used to split Sample ID from contigname') 
    optional.add_argument('-m',
                          '--minimum_length',
                          dest='minimum_length',
                          metavar='',
                          default='2000',
                          required=False,
                          type=str,
                          help='Minimum contig length') 

    optional.add_argument('-g',
                          '--decontaminate',
                          dest='decontaminate',
                          action='store_true',
                          default=False,
                          help='Flag for running RF prediction of Viral Bins')
    optional.add_argument('-f',
                          '--fasta',
                          dest='fasta',
                          required=False,
                          default=None,
                          help='Same fasta used as input to VAMB')

    (args, extra_args) = parser.parse_known_args()

    return args



class bin_annotation:
    def __init__(self,bin_name,motherbin):
        self.motherbin = motherbin
        self.bin_name = bin_name
        self.ncontigs = None
        self.binsize = 0
     


### Initialising functions

def define_tasks(args):
    '''Checks the presence of files generated by the pipeline - files won't be created if already there
       User has to delete these files first.
    '''
    tasks = []

    ### Check for PVOG and MiComplete105
    for hmm in ['hmmMiComplete105','hmmVOG']:
        cleaned_hmmfile_out = os.path.join(args.outdir,hmm+'.allsamples.txt')
        if os.path.exists(cleaned_hmmfile_out):
            print('{} file already there - delete to generate again'.format(hmm))
        else:
            tasks.append(hmm)

    ### Check for DVF
    cleaned_DVF_out = os.path.join(args.outdir,'DVFpredictions.allsamples.txt')
    if os.path.exists(cleaned_DVF_out):
        print('DVF file already there - delete to generate again')
    else:
        tasks.append('DVF')

    ### Check for Prodigal
    allProdigal_out = os.path.join(args.outdir,'All.predicted_proteins.faa')
    allProdigal_out_gz =  os.path.join(args.outdir,'All.predicted_proteins.faa.gz')
    if os.path.exists(allProdigal_out) or os.path.exists(allProdigal_out_gz):
        print('Prodigal file already there - delete to generate again')
    else:
        tasks.append('Prodigal')
    

    ### Check if Bin-annotation tablee is already made 
    bintable = os.path.join(args.outdir,'vambbins_aggregated_annotation.txt')
    if os.path.exists(bintable):
        print('Bintable donee already - delete to generate again')
    else:
        tasks.append('Bintable')

    ### Run Random-Forrest Decontamination? 
    if args.decontaminate:
        tasks.append('Decontaminate')
    else:
        print('Not running Decontamination Mode (use --decontaminate flag for that.')

    return(tasks)


def parse_sample_table(args):
    '''Return list of samples'''
    samples = {}
    with open(args.sample_table,'r') as infile:
        for line in infile:
            line = line.strip().split()
            sample = line[0]
            samples[sample] = sample
    return(samples)

def read_in_clusters(args):
    cluster_file = args.cluster_file
    
    print('Reading in Clusters, this may take a while...')
    binsannos = dict()
    cls = dict()
    with open(cluster_file,'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            cluster, contig = line[0],line[1]
            sample = contig.split('_')[0]
            binid = sample + '_' + cluster
            if 'length' in contig:
                contiglength = int(contig.split('length')[1].split('_')[1])


            if binid not in binsannos:
                binsannos[binid] = bin_annotation(binid,cluster)
                binsannos[binid].ncontigs = 1
                binsannos[binid].binsize += contiglength

            else:
                binsannos[binid].ncontigs += 1
                binsannos[binid].binsize += contiglength

            if contig not in cls:
                cls[contig] = (cluster,binid)
    return cls, binsannos




############### For Blastfiles ###############
def parse_blastfiles(args,task_list,samples):
    '''Generic function parse files generated with blastn/blastp'''
    query_type_ = {'IMGVR':'contig','PLSDB':'contig','VIRALREFSEQ':'contig'}

    for db in ['IMGVR','PLSDB','VIRALREFSEQ']:
        if db not in task_list:
            continue
        cleaned_m6_out = os.path.join(args.outdir,db+'.allsamples.m6')
        
        with open(cleaned_m6_out,'w') as outfile:
            for sample in samples:
                blastfile = os.path.join(args.annotation_dir,sample,sample+'.'+db+'.m6') 
                if not os.path.isfile(blastfile):
                    print('{} does not exist for sample: {}'.format(blastfile,sample))
                    continue
                for string in blast_outline_generator(blastfile,query_type_[db],sample,args.split_delimn):
                    outfile.write(string)

def blast_outline_generator(blastfile,query_type,sample,seperator='_'):
    '''Simple generator to parse and yield blast lines (either contig or protein)
       It assumes that PLSDB hits are proteins and thus adds a contig column for identification
       While it does little to IMGVR hits
    '''
    with open(blastfile,'r') as infile:
        for line in infile:
            line = line.strip().split()

            if query_type == 'contig':
                string = '\t'.join(line) + '\n'
                yield string
            elif query_type == 'protein':
                query_split = line[0].split(seperator)
                contig = seperator.join(query_split[:-1])
                line.insert(0,contig)
                string = '\t'.join(line) + '\n'
                yield string

                
############### For HMMfiles ###############
def parse_hmmfiles(args,task_list,samples):
    '''Generic function filter and parse files generated with hmmsearch'''
    
    for hmm in ['hmmMiComplete105','hmmVOG']:
        if hmm not in task_list:
            continue

        cleaned_hmmfile_out = os.path.join(args.outdir,hmm+'.allsamples.txt')

        with open(cleaned_hmmfile_out,'w') as outfile:
            for sample in samples:
                hmmfile = os.path.join(args.annotation_dir,sample,sample+'.'+hmm+'.tbl') 
                if not os.path.isfile(hmmfile):
                    print('{} does not exist for sample: {}'.format(hmmfile,sample))
                    continue
                for string in hmm_outline_generator(hmmfile,sample,args.split_delimn):
                    outfile.write(string)

def hmm_outline_generator(hmmfile,sample,seperator='_'):
    with open(hmmfile,'r') as infile:
        target_protein = dict()
        for line in infile:
            if line[0] != '#':
                entry = line.strip().split()
                if len(entry) > 8:
                    protein, target, full_seq_eval , full_seq_score , best_dom_eval, best_dom_score =  entry[0], entry[2], entry[4], entry[5], entry[7],entry[8]
                    tmp = protein.split(seperator)
                    genome_contig = seperator.join(tmp[:-1])
                    full_seq_score = float(full_seq_score)
                    if float(full_seq_score) <= 50:
                        continue
                    if protein not in target_protein:
                        target_protein[protein] = []
                        target_protein[protein] += [(full_seq_score,genome_contig,protein,target,full_seq_eval)]
                    else:
                        target_protein[protein] += [(full_seq_score,genome_contig,protein,target,full_seq_eval)]
    for protein in target_protein:
        protein_targets = target_protein[protein]
        
        ### sort by domain score 
        protein_targets = sorted(protein_targets, key=lambda x: x[0] ,reverse=True)
        full_seq_score , contig, ORF, target, full_seq_eval = protein_targets[0]
        subentry = [contig, ORF, target, full_seq_score, full_seq_eval]
        subentry = [str(item) for item in subentry]
        string = '\t'.join(subentry) + '\n'
        yield string


############### For DVF files ############### 
def parse_DVF(args,task_list,samples):
    '''Generic function filter and parse files generated by DVF 
    '''
    if 'DVF' not in task_list:
        pass
    else:
        cleaned_DVF_out = os.path.join(args.outdir,'DVFpredictions.allsamples.txt')
        with open(cleaned_DVF_out,'w') as outfile:
            for sample in samples:
                #dvffile = os.path.join(args.annotation_dir,sample,sample + '_dvf',sample + '.flt.fa_gt' + args.minimum_length + 'bp_dvfpred.txt')
                dvfdir = os.path.join(args.annotation_dir,sample,sample + '_dvf')
                dvffiles = [f for f in os.listdir(dvfdir) if os.path.isfile(os.path.join(dvfdir, f))]
                if len(dvffiles) == 0:
                    print('DVF does not exist for sample: {}'.format(sample))
                    continue
                dvffile = os.path.join(args.annotation_dir,sample,sample + '_dvf',dvffiles[0])
                for string in DVF_generator(dvffile,args.split_delimn):
                    outfile.write(string)

def DVF_generator(dvffile,seperator):
    with open(dvffile,'r') as infile:
        for line in infile:
            if line[:4] =='name':
                continue
            line = line.strip().split('\t')
            contig,score,pval = line[0],line[2],line[3]
            sample = contig.split(seperator)[0]
            dvf_entry = [contig,score,pval,sample]
            string = '\t'.join(dvf_entry) + '\n'
            yield string

############### For Predicted Proteins ############### 

def parse_Prodigal(args,task_list,samples):
    '''Generic function to concatenate Predicted Proteins by Prodigal 
    '''
    if 'Prodigal' not in task_list:
        pass
    else:
        allProdigal_out = os.path.join(args.outdir,'All.predicted_proteins.faa')
        with open(allProdigal_out,'w') as outfile:
            for sample in samples:
                Prodigalfile = os.path.join(args.annotation_dir,sample,sample + '.predicted_proteins.faa')
                if not os.path.isfile(Prodigalfile):
                    print('{} does not exist for sample: {}'.format(Prodigalfile,sample))
                    continue
                for string in Prodigal_generator(Prodigalfile):
                    outfile.write(string)

def Prodigal_generator(Prodigalfile):
    '''There is actually not much to do but concatenate the lines of Prodigal files
       This generator can be expanded to also accomodate filters and header modifications.
    '''
    with open(Prodigalfile,'r') as infile:
        for line in infile:
            if line[0] != '#':
                yield line



############### Summarise annotation for Bins ###############

def parse_bacterial_hallmarks(args, cls):
    '''
    Count the number of Distinct Bacterial Hallmark Genes in each Cluster
    '''

    hallmark_file = os.path.join(args.outdir,'hmmMiComplete105.allsamples.txt')

    cluster_hallmarks = dict()
    with open(hallmark_file,'r') as infile:
        for line in infile:
            contig, protein, hallmark_gene, score, evalue = line.strip().split('\t')

            if contig in cls:
                cluster = cls[contig][0]

                if cluster not in cluster_hallmarks:
                    cluster_hallmarks[cluster] = set([hallmark_gene])
                else:
                    cluster_hallmarks[cluster].add(hallmark_gene)
    
    cluster_hallmark_count = dict()
    for cluster in cluster_hallmarks:
        n_distinct_hallmarks = len(cluster_hallmarks[cluster])
        cluster_hallmark_count[cluster] = n_distinct_hallmarks
    return cluster_hallmark_count

def parse_PVOG_hits(args,cls,binsannos,hmmscore = 50):
    '''
    Count the number of distinct VOG annotations scaled by ncontigs total
    '''

    vogfile = os.path.join(args.outdir,'hmmVOG.allsamples.txt')
    if not os.path.exists(vogfile):
        print('Missing VOG file !')
        sys.exit(1)
    
    
    bin_vogs = dict()
    with open(vogfile,'r') as infile:
        for line in infile:
            contig, protein, vog, score, evalue = line.strip().split('\t')
            if float(score) <= hmmscore:
                continue

            if contig in cls:
                binid = cls[contig][1]
                if binid not in bin_vogs:
                    bin_vogs[binid] = set([vog])
                else:
                    bin_vogs[binid].add(vog)
    
    ### Summarise VOGs for each bin 
    bin_vog_count = dict()
    for binid in bin_vogs:
        n_contigs = binsannos[binid].ncontigs
        vog_count_scaled = len(bin_vogs[binid])/n_contigs
        bin_vog_count[binid] = vog_count_scaled
    return bin_vog_count


def parse_IMGVR_hits(args,cls,binsannos,seqid=95,coverage = 0.5):
    '''
    Calculate How many Contigs in each bin with hits to IMGVR genomes
    '''

    imgvr_file = os.path.join(args.outdir,'IMGVR.allsamples.m6')
    if not os.path.exists(imgvr_file):
        print('Missing IMGVR file !')
        sys.exit(1)

    bin_imgvrhits = dict()
    with open(imgvr_file,'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            contig = line[0]
            sequence_identity = float(line[2])
            alignment_length = int(line[3])
            imgvr_genome_length = int(line[-1])
            if contig in cls:
                binid = cls[contig][1]

                if sequence_identity >= seqid and alignment_length/imgvr_genome_length >= coverage:
                    if binid not in bin_imgvrhits:
                        bin_imgvrhits[binid] = set([contig])
                    else:
                        bin_imgvrhits[binid].add(contig)
    
    ### Summarise strong hits for each bin
    bin_IMGVR_fraction = dict()
    for binid in bin_imgvrhits:
        n_contigs = binsannos[binid].ncontigs
        nhits = len(bin_imgvrhits[binid])
        bin_imgvr_fraction = nhits/n_contigs
        bin_IMGVR_fraction[binid] = bin_imgvr_fraction

    return bin_IMGVR_fraction

        
def parse_DVF_preds(args,cls):
    '''
    Calculate median DVF score for each bin and clustre
    '''

    dvf_file = os.path.join(args.outdir,'DVFpredictions.allsamples.txt')
    if not os.path.exists(dvf_file):
        print('Missing DVF file !')
        sys.exit(1)

    cluster_dvf_scores = dict()
    bin_dvf_scores = dict()
    with open(dvf_file,'r') as infile:
        for line in infile:
            contig, score, pvalue, sample = line.strip().split('\t')

            if contig in cls:
                score = float(score)
                binid = cls[contig][1]
                cluster = cls[contig][0]
                if cluster not in cluster_dvf_scores:
                    cluster_dvf_scores[cluster] = [score]
                else:
                    cluster_dvf_scores[cluster] = [score]

                if binid not in bin_dvf_scores:
                    bin_dvf_scores[binid] = [score]
                else:
                    bin_dvf_scores[binid] += [score]
    
    ### calculate median DVF score pr. bin 
    bin_median_dvf = dict()
    for binid in bin_dvf_scores:
        median_score = np.median(bin_dvf_scores[binid])
        bin_median_dvf[binid] = median_score

    ### Same for Cluster
    cluster_median_dvf = dict()
    for cluster in cluster_dvf_scores:
        median_score = np.median(cluster_dvf_scores[cluster])
        cluster_median_dvf[cluster] = median_score
    
    return bin_median_dvf, cluster_median_dvf


def aggregate_bin_annotation(args,cls,binsannos):


    print('Parsing Bacterial Hallmarks')
    cluster_hallmark_count = parse_bacterial_hallmarks(args, cls)

    print('Parsing VOG HMM search')
    bin_vog_count = parse_PVOG_hits(args,cls,binsannos)

    #print('Parsing IMGVR blast-hits')
    #bin_IMGVR_fraction = parse_IMGVR_hits(args,cls,binsannos,seqid=95,coverage = 0.5)

    print('Parsing DVF contig scores - this may take some timee...')
    bin_median_DVF, cluster_median_DVF = parse_DVF_preds(args,cls)

    ### Write out table with aggregated information for each Bin
    fileout = os.path.join(args.outdir,'vambbins_aggregated_annotation.txt')

    if os.path.exists(fileout):
        print(fileout,' allready exists, delete it to reproduce it!')
        sys.exit(0)

    print('Writing results to:',fileout)
    with open(fileout,'w') as out:
        header = ['binid','ncontigs','binsize','motherbin','nhallm','nVOGs','bin_DVF_score','cluster_DVF_score']

        out.write('\t'.join(header) + '\n')
        for binid in binsannos:

            lineout = []
            lineout.append(binid)

            ### add bin annotation
            ncontigs = binsannos[binid].ncontigs
            binsize = binsannos[binid].binsize

            ### Going this low...
            if binsize <= 2500:
                continue
                     
            mothercluster = binsannos[binid].motherbin
            if mothercluster not in cluster_hallmark_count:
                nhallmarks = 0
            else:
                nhallmarks = cluster_hallmark_count[mothercluster]

            if binid not in bin_vog_count:
                nVOGs = 0
            else:
                nVOGs = bin_vog_count[binid]
            
            if binid not in bin_median_DVF:
                bin_DVF_score = 0
            else:
                bin_DVF_score = bin_median_DVF[binid]
            
            if mothercluster not in cluster_median_DVF:
                cluster_DVF_score = 0
            else:
                cluster_DVF_score = cluster_median_DVF[mothercluster]

            ### 
            lineout += [ncontigs,binsize, mothercluster, nhallmarks, round(nVOGs,2), round(bin_DVF_score,2),round(cluster_DVF_score,2)]
            lineout = [str(x) for x in lineout]

            out.write('\t'.join(lineout) + '\n')



def write_concate_sequences(args,cls,_vambbins_filtered):
    '''
    Write out Contigs Concatenated as Bins 
    Write out Contigs individually
    '''
    
    viral_contigs_out = os.path.join(args.outdir,'VAMBV3.Viral_RF_predictions.contigs.fna')
    viral_bins_out = os.path.join(args.outdir,'VAMBV3.Viral_RF_predictions.bins.fna')
    binids = set(_vambbins_filtered.binid)

    print('Parsing Fasta sequences - this may take a while...')
    clusters = {}
    with open(viral_contigs_out,'w') as out:
        for record in SeqIO.parse(open(args.fasta, 'r'), 'fasta'):
            contigname = record.id 
            if contigname not in cls:
                continue
            binid = cls[contigname][1]
            if binid in binids:
                out.write('>{}\n{}\n'.format(contigname,record.seq))
                if binid not in clusters:
                    clusters[binid] = ""
                    clusters[binid] += record.seq
                else:
                    clusters[binid] += record.seq
    
    ### Write out concatenated bins
    with open(viral_bins_out,'w') as out:
        for binid in clusters:
            sequence = clusters[binid]
            out.write('>{}\n{}\n'.format(binid,sequence))




def RF_decontaminate(args):

    table_file = os.path.join(args.outdir,'vambbins_aggregated_annotation.txt')
    if not os.path.exists(table_file):
        print('Ups! Where is your VAMB bins file?:',table_file)
        sys.exit(1)

    print('Loading Model and annotation table')
    trained_model = joblib.load('dbs/RF_model.sav')
    _vambbins = pd.read_csv(table_file,sep='\t')
    _subset = _vambbins[['binsize','nhallm','nVOGs','cluster_DVF_score']]
    eval_predictions = trained_model.predict(_subset)
    _vambbins['Prediction'] = eval_predictions

    _vambbins_filtered = _vambbins[ (_vambbins.Prediction =='Viral') & (_vambbins.nVOGs >= 1) & (_vambbins.cluster_DVF_score >= 0.3) ]
    
    table_file_annotated = os.path.join(args.outdir,'vambbins_aggregated_annotation.Viral.txt')
    _vambbins_filtered.sort_values(by='binsize',ascending=False).to_csv(table_file_annotated,sep='\t',index=False)

    ### Just in case, read clusters from VAMB clusters file.
    cls, binsannos = read_in_clusters(args)


    if args.fasta is None:
        print('You need to provide the -f argument in order to write fasta files of Viral Bins and Contigs!')
        sys.exit(0)
    ### Write out contigs and Bins concatenated
    write_concate_sequences(args,cls,_vambbins_filtered)
    print('Done Writing Viral Sequences!')



    
    
    












############### Main ############### 
if __name__ == "__main__":
    args = parse_arguments()

    ### Get sample names
    samples = parse_sample_table(args)

    ### Define tasks 
    task_list = define_tasks(args)

    if len(task_list) == 0:
        print('All expected outfiles already exist in {} - We are done here'.format(args.outdir))
        sys.exit(0)
    try:
        os.mkdir(args.outdir)
    except:
        pass

    ### Run parsers
    parse_hmmfiles(args,task_list,samples)
    #parse_blastfiles(args,task_list,samples)
    parse_DVF(args,task_list,samples)
    parse_Prodigal(args,task_list,samples)
    print('Done')

    if 'Bintable' in task_list:

        ### Read in VAMB cluster information
        cls, binsannos = read_in_clusters(args)

        ### Calculate and summarise Bin-information
        aggregate_bin_annotation(args,cls,binsannos)
    
    if 'Decontaminate' in task_list:
        print('Running RF - decontaminate')
        RF_decontaminate(args)
        






         
         




