"""
"""

from __future__ import print_function
from __future__ import division

import glob
import os
import json
import gzip
import sys
import argparse

IS_PYTHON3 = sys.version[0] == "3"

VCF_DIR = "/lustre/scratch114/projects/ddd/release/20140912/final"
CHROM = "22"

consequence_counts = {"transcript_ablation": 0, "splice_donor_variant": 0,
    "splice_acceptor_variant": 0, "stop_gained": 0, "frameshift_variant": 0,
    "stop_lost": 0, "initiator_codon_variant": 0, "transcript_amplification": 0,
    "inframe_insertion": 0, "inframe_deletion": 0, "missense_variant": 0,
    "splice_region_variant": 0, "incomplete_terminal_codon_variant": 0,
    "stop_retained_variant": 0, "synonymous_variant": 0,
    "coding_sequence_variant": 0, "mature_miRNA_variant": 0,
    "5_prime_UTR_variant": 0, "3_prime_UTR_variant": 0,
    "non_coding_transcript_exon_variant": 0, "non_coding_exon_variant": 0, "intron_variant": 0,
    "NMD_transcript_variant": 0, "non_coding_transcript_variant": 0,
    "upstream_gene_variant": 0, "downstream_gene_variant": 0, "TFBS_ablation": 0,
    "TFBS_amplification": 0, "TF_binding_site_variant": 0,
    "regulatory_region_ablation": 0, "regulatory_region_amplification": 0,
    "regulatory_region_variant": 0, "feature_elongation": 0,
    "feature_truncation": 0, "intergenic_variant": 0,
    "last_base_of_exon_G": 0, "nc_transcript_variant": 0}

def get_options():
    """ get the command line options
    """
    
    parser = argparse.ArgumentParser(description="Counts singletons within \
        multiple sample VCFs.")
    parser.add_argument("--chrom", required=True, help="chromosome to investigate.")
    parser.add_argument("--last-base-sites", default="/lustre/scratch113/projects/ddd/users/jm33/last_base_sites.json",
        help="path to list of last base sites")
    parser.add_argument("--singletons", default=sys.stdout,
        help="path to send the list of singletons to.")
    parser.add_argument("--totals", default=sys.stdout,
        help="path to send the totals for each consequence class.")
    
    args = parser.parse_args()
    
    return args

def get_ddd_parents():
    """ get a dictionary of unaffected DDD parents, to their sex
    """
    
    DIR = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/"
    family_path = os.path.join(DIR, "family_relationships.txt")
    sanger_id_path = os.path.join(DIR, "person_sanger_decipher.txt")
    
    # find all the unaffected parental DDD IDs.
    parental_ids = {}
    with open(family_path, "r") as handle:
        for line in handle:
            line = line.strip().split("\t")
            person_id = line[1]
            paternal_id = line[2]
            sex = line[4]
            affected_status = line[5]
            
            if paternal_id != "0" or affected_status != "1":
                continue
            
            parental_ids[person_id] = sex
    
    # get the unaffected parental sanger IDs (since the sanger IDs are the IDs
    # used in the multisample VCFs)
    sanger_ids = {}
    with open(sanger_id_path, "r") as handle:
        for line in handle:
            line = line.strip().split("\t")
            person_id = line[0]
            decipher_id = line[1]
            sanger_id = line[2]
            
            if person_id in parental_ids:
                sanger_ids[sanger_id] = parental_ids[person_id]
    
    return sanger_ids

def get_vep_annotations(chrom):
    
    vep_dir = "/lustre/scratch113/projects/ddd/users/ddd/test_msVCF_annotate_vep/vep_annotated_vcfs/results_vcfs"
    vep_path = os.path.join(vep_dir, "{}.txt".format(chrom))
    
    vep_consequences = {}
    with open(vep_path, "r") as handle:
        exclude_header(handle)
        for line in handle:
            line = line.strip().split("\t")
            
            pos = int(line[1])
            # ref = line[3]
            # alt = line[4]
            info = line[7]
            
            consequence = info.split("CQ=")[1]
            consequence = consequence.split(";", 1)[0].split(",")
            
            vep_consequences[pos] = consequence
    
    return vep_consequences

def get_last_base_sites(path):
    """ get the list of last base sites
    """
    
    with open(path) as handle:
        last_base = json.load(handle)
    
    last_base = set([(x[0], int(x[1])) for x in last_base])
    
    return last_base

def get_header(vcf):
    """ get the full VCF header as a list of lines
    
    Args:
        vcf: handle to vcf file
    """
    
    current_pos = vcf.tell()
    vcf.seek(0)
    
    header = []
    for line in vcf:
        if not line.startswith("#"):
            break
        
        header.append(line)
    
    vcf.seek(current_pos)
    
    return(header)

def exclude_header(vcf):
    """ move the file handle to just beyond the VCF header lines
    """
    
    current_pos = vcf.tell()
    
    while vcf.readline().startswith("#"):
        current_pos = vcf.tell()
    
    vcf.seek(current_pos)

def get_sample_positions(header):
    """ make a dictionary of sample IDs (from the VCF header) vs their position
    
    Args:
        header: list of VCF header lines, the final line defines the VCF columns
    """
    
    samples = {}
    
    sample_ids = header[-1].strip().split("\t")[9:]
    sample_map = dict(zip(sample_ids, range(0, len(sample_ids))))
    
    return sample_map

def get_variant_key(variant):
    """ get the chrom, pos, and alleles for a VCF variant
    
    Args:
        variant: list of data from VCF line for the first 8 columns
    """
    
    chrom = variant[0]
    pos = int(variant[1])
    ref = variant[3]
    alts = variant[4]
    
    return (chrom, pos, ref, alts)

def get_format(format_string):
    """ figure out the format of the sample columns, map data to position
    
    Args:
        format_string: text from format column of VCF, colon separated string
    """
    
    format = format_string.split(":")
    format = dict(zip(format, range(0, len(format))))
    
    return(format)

def get_sample_genotypes(samples, format, sample_pos):
    """ get a dictionary of genotypes for specific sample IDs
    """
    
    samples = samples.strip().split("\t")
    
    genotypes = {}
    for sample_id in sample_pos:
        sample = samples[sample_pos[sample_id]]
        sample = sample.split(":")
        genotype = sample[format["GT"]]
        
        # drop missing genotypes, since they cannot contribute to the total
        # allele count
        if genotype == "./.":
            continue
        
        genotypes[sample_id] = genotype
    
    return genotypes

def reformat_chrX_genotypes(key, genotypes, ddd_parents):
    """ swap male genotypes on chrX to a hemizgous type
    """
    
    # don't alter genotypes not on the X
    if key[0] != "X":
        return genotypes
    
    # define the pseudoautosomal regions on the X chromosome
    x_par = [(60001, 2699520), (154930290, 155260560), (88456802, 92375509)]
    
    # don't alter genotypes within the pseudoautosomal regions
    if any([key[1] >= x[0] and key[1] < x[1] for x in x_par ]):
        return genotypes
    
    exclude_ids = []
    for sample_id in genotypes:
        # don't alter female parents, since they are diploid for chrX
        if ddd_parents[sample_id] != "1":
            continue
        
        geno = genotypes[sample_id].split("/")
        
        # drop genotypes for heterozygous males on chrX
        if geno[0] != geno[1]:
            exclude_ids.append(sample_id)
        
        genotypes[sample_id] = geno[0]
    
    # remove the abberrant heterozygous chrX male genotypes
    for sample_id in exclude_ids:
        del genotypes[sample_id]
    
    return genotypes

def tally_alleles(genotypes, alts):
    """ count the alleles used in the genotypes
    """
    
    # make sure we have entries for all the alt alleles, even if the count
    # ends up as zero, that way when we later match to the consequences for the alt
    # alleles, the allele counts should be consistent with the consequences
    alts = alts.split(",")
    allele_numbers = ["{}".format(x + 1) for x in range(len(alts))]
    allele_numbers.append("0")
    
    # tally each allele found in the genotypes
    counts = dict(zip(allele_numbers, [0]* len(allele_numbers)))
    for genotype in genotypes.values():
        genotype = genotype.split("/")
        
        for allele in genotype:
            counts[allele] += 1
    
    return counts

def check_singletons(key, counts, vep, is_last_base, output):
    """ checks to see if any of the alleles at a site are singletons
    
    Args:
        key: tuple of (chrom, pos, ref, alt)
        counts: dictionary of number of alleles seen in the unaffecetd DDD parents
        vep: dictionary of consequence strings for each allele, indexed by chrom position
        is_last_base: True/False for whether the variant is at the last base of an
            exon with a G reference.
    """
    
    consequences = vep[key[1]]
    
    # remove the reference allele, since we don't have a consequence for that
    del counts["0"]
    
    allele_numbers = sorted(counts)
    
    for number in allele_numbers:
        allele_count = counts[number]
        
        # don't include alleles with zero alleles in the unaffected parents,
        # otherwise they will skew the proportion of variants as singletons.
        if allele_count == 0:
            continue
        
        consequence = consequences[int(number) - 1]
        
        if is_last_base:
            consequence = "last_base_of_exon_G"
        
        consequence_counts[consequence] += 1
        
        # we only want to look at singletons
        if allele_count != 1:
            continue
        
        line = "{}\t{}\t{}\n".format(key[0], key[1], consequence)
        
        output.write(line)

def parse_vcf(chrom, ddd_parents, vep, last_base, output_path):
    """ run through a VCF, counting the alleles, looking for singletons
    """
    
    try:
        output = open(output_path, "w")
    except TypeError:
        output = output_path
    
    vcf_path = glob.glob(os.path.join(VCF_DIR, "{}:1-*.vcf.gz".format(chrom)))
    vcf_path = vcf_path[0]
    
    vcf = gzip.open(vcf_path, "r")
    if IS_PYTHON3:
        vcf = gzip.open(vcf_path, "rt")
    
    # remove the header lines
    header = get_header(vcf)
    exclude_header(vcf)
    
    # figure out where the DDD parents are in the smaple list
    index = get_sample_positions(header)
    sample_pos = dict([(x, index[x]) for x in ddd_parents if x in index])
    
    for line in vcf:
        line = line.split("\t", 9)
        variant = line[:9]
        key = get_variant_key(variant)
        
        format = get_format(variant[8])
        genotypes = get_sample_genotypes(line[9], format, sample_pos)
        genotypes = reformat_chrX_genotypes(key, genotypes, ddd_parents)
        counts = tally_alleles(genotypes, key[3])
        pos = (key[0], key[1])
        is_last_base = pos in last_base
        check_singletons(key, counts, vep, is_last_base, output)

def main():
    args = get_options()
    ddd_parents = get_ddd_parents()
    vep = get_vep_annotations(args.chrom)
    last_base = get_last_base_sites(args.last_base_sites)
    
    try:
        parse_vcf(args.chrom, ddd_parents, vep, last_base, args.singletons)
    finally:
        # and finally show how many times each consequence was seen in the DDD
        # unaffected parents, so we can determine the proportion of that each
        # consequence is seen as a singleton
        try:
            output = open(args.totals, "w")
        except TypeError:
            output = args.totals
        
        output.write("consequence\toccurrences\n")
        for consequence in sorted(consequence_counts):
            line = "{}\t{}\n".format(consequence, consequence_counts[consequence])
            output.write(line)


if __name__ == "__main__":
    main()
