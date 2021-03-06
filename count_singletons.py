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

consequence_counts = {"transcript_ablation": 0, "splice_donor_variant": 0,
    "splice_acceptor_variant": 0, "stop_gained": 0, "frameshift_variant": 0,
    "stop_lost": 0, "initiator_codon_variant": 0, "transcript_amplification": 0,
    "inframe_insertion": 0, "inframe_deletion": 0, "missense_variant": 0,
    "splice_region_variant": 0, "incomplete_terminal_codon_variant": 0,
    "stop_retained_variant": 0, "synonymous_variant": 0,
    "coding_sequence_variant": 0, "mature_miRNA_variant": 0,
    "5_prime_UTR_variant": 0, "3_prime_UTR_variant": 0,
    "non_coding_transcript_exon_variant": 0, "non_coding_exon_variant": 0,
    "intron_variant": 0, "NMD_transcript_variant": 0,
    "non_coding_transcript_variant": 0, "upstream_gene_variant": 0,
    "downstream_gene_variant": 0, "TFBS_ablation": 0,
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
    parser.add_argument("--last-base-sites",
        default="/lustre/scratch113/projects/ddd/users/jm33/last_base_sites.json",
        help="path to list of last base sites")
    parser.add_argument("--singletons", default=sys.stdout,
        help="path to send the list of singletons to.")
    parser.add_argument("--totals", default=sys.stdout,
        help="path to send the totals for each consequence class.")
    parser.add_argument("--gq-mean", type=float, default=0.0, help="minimum GQ_MEAN of \
        variants to include")
    parser.add_argument("--qual", type=float, default=0.0, help="minimum QUAL of \
        variants to include")
    
    args = parser.parse_args()
    
    return args

def get_ddd_parents():
    """ get a dictionary of unaffected DDD parents, to their sex
    
    Returns:
        dictionary of sample sanger IDs to their sex eg {"DDD_MAIN0001": "1"}
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
    """ extracts all the VEP consequence strings for variants on a chromosome
    
    We load all the VEP consequences for all the variants on a single chrom,
    known to occur in the DDD population. This way, when we actually parse the
    genotypes for the variant, we can instantly get the consequence
    
    Args:
        chrom: chromosome string e.g. "12", "X"
    
    Returns:
        dictionary of positions on the chromosome to the consequence strings
        e.g. {12000: ["missense_variant"], 120001: ["synonymous_variant"]}
    """
    
    vep_dir = "/lustre/scratch113/projects/ddd/users/ddd/test_msVCF_annotate_vep/vep_annotated_vcfs/results_vcfs"
    vep_path = os.path.join(vep_dir, "{}.txt".format(chrom))
    
    vep_consequences = {}
    with open(vep_path, "r") as handle:
        exclude_header(handle)
        for line in handle:
            line = line.strip().split("\t")
            
            pos = int(line[1])
            info = line[7]
            
            # extract the consequence field
            consequence = info.split("CQ=")[1]
            
            # tidy the consequence string, and split into a loist if there are
            # more than one consequence string for different alleles.
            consequence = consequence.split(";", 1)[0].split(",")
            
            vep_consequences[pos] = consequence
    
    return vep_consequences

def get_last_base_sites(path):
    """ get the list of last base sites
    
    Args:
        path: path to the json file containing the positions at the last base of
            an exon where the base is a G (or other bases if required).
    
    Returns:
        set of (chrom, positon) tuples, so we can quickly check whether sites in
        the VCF are within this set.
    """
    
    with open(path) as handle:
        last_base = json.load(handle)
    
    last_base = set([(x[0], int(x[1])) for x in last_base])
    
    return last_base

def get_header(vcf):
    """ get the full VCF header as a list of lines
    
    Args:
        vcf: handle to vcf file
    
    Returns:
        list of header lines
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
    
    Args:
        vcf: file handle to VCF
    """
    
    current_pos = vcf.tell()
    
    while vcf.readline().startswith("#"):
        current_pos = vcf.tell()
    
    vcf.seek(current_pos)

def get_sample_positions(header):
    """ make a dictionary of sample IDs (from the VCF header) vs their position
    
    The VCF contains many samples, and we only want to investigate a subset of
    these. We know which samples we wish to use (defined by the family
    relationships). In order to be able to get this subset in each line of the
    VCF, we need to figure out which column belongs to which sample.
    
    Args:
        header: list of VCF header lines, the final line defines the VCF columns
    
    Returns:
        dictionary of sample IDs (from the VCF header) vs their position.
    """
    
    samples = {}
    
    sample_ids = header[-1].strip().split("\t")[9:]
    sample_map = dict(zip(sample_ids, range(0, len(sample_ids))))
    
    return sample_map

def get_variant_key(variant):
    """ get the chrom, pos, and alleles for a VCF variant
    
    Args:
        variant: list of data from VCF line for the first 8 columns
    
    Returns:
        unique key for the variant as tuple (chrom, pos, ref allele, alt alleles)
    """
    
    chrom = variant[0]
    pos = int(variant[1])
    ref = variant[3]
    alts = variant[4]
    
    return (chrom, pos, ref, alts)

def get_format_index(format_string):
    """ figure out the format of the sample columns, map data to position
    
    Args:
        format_string: text from the format column of VCF, a colon separated string.
    
    Return:
        dictionary of format key to index positions
    """
    
    format = format_string.split(":")
    format = dict(zip(format, range(0, len(format))))
    
    return(format)

def parse_samples(samples, format, sample_pos):
    """ get a dictionary of genotypes for specific sample IDs
    
    We exclude samples with missing genotypes, since their data shouldn't be
    used for any subsequent analysis.
    
    Args:
        samples: list of format fields for samples e.g. ["0/1:10:.:.:.",
            "1/1:10:.:.:.", "./.:10:.:.:."].
        format: mapping of format keys to their index in the format field
        sample_pos: dictionary of sample IDs to their position in the samples
            list. This should only be a subset of the columns from the samples
            list.
    
    Returns:
        a dictionary of sample format values e.g.
        {"DDD_MAIN01": {"GT": "0/1", "DP": "10", "GQ": "20"},
         "DDD_MAIN02": {"GT": "0/0", "DP": "20", "GQ": "30"}
    """
    
    samples = samples.strip().split("\t")
    
    parsed = {}
    for sample_id in sample_pos:
        idx = sample_pos[sample_id]
        
        sample_field = samples[idx]
        sample_field = sample_field.split(":")
        
        # put all the format values for the sample into a dictionary
        data = {}
        for key in format:
            index = format[key]
            
            # sometimes the sample values format column doesn't have the same
            # number of fields as defined in the format keys column. This only
            # occurs in samples that lack genotype calls, so it's not too much
            # of a problem, as we exclude these samples anyway.
            try:
                data[key] = sample_field[index]
            except IndexError:
                data[key] = None
        
        # drop samples missing genotypes, since they cannot contribute to the
        # total allele count
        if data["GT"] == "./.":
            continue
        
        parsed[sample_id] = data
    
    return parsed

def reformat_chrX_genotypes(key, samples, ddd_parents):
    """ swap male genotypes on chrX to a hemizgous type
    
    Args:
        key: tuple of (chrom, pos, ref, alts) for the variant
        samples: dictionary of sample IDs to sample-based dictionaries
        ddd_parents: dictionary of sample IDs to sex
    
    Returns:
        samples dictionary, but with hemizygous genotypes modified if required.
    """
    
    # don't alter genotypes not on the X
    if key[0] != "X":
        return samples
    
    # define the pseudoautosomal regions on the X chromosome
    x_par = [(60001, 2699520), (154930290, 155260560), (88456802, 92375509)]
    
    # don't alter genotypes within the pseudoautosomal regions
    if any([key[1] >= x[0] and key[1] < x[1] for x in x_par ]):
        return samples
    
    exclude_ids = []
    for sample_id in samples:
        # don't alter female parents, since they are diploid for chrX
        if ddd_parents[sample_id] != "1":
            continue
        
        alleles = samples[sample_id]["GT"].split("/")
        
        # figure out if the genotype is a heterozygous for a male on chrX,
        # since these are typically errors that we want to exclude.
        if alleles[0] != alleles[1]:
            exclude_ids.append(sample_id)
        
        samples[sample_id]["GT"] = alleles[0]
    
    # remove the abberrant heterozygous chrX male genotypes
    for sample_id in exclude_ids:
        del samples[sample_id]
    
    return samples

def get_singleton_data(samples):
    """ get the data for the sample with a singleton
    
    Args:
        samples: dictionary of parsed sample values for the population
        
    Returns:
        dictionary entry for the singleton sample.
    """
    
    # find the sample ID for the sole sample without a ref genotype
    for sample_id in samples:
        geno = samples[sample_id]["GT"]
        
        if geno != "0/0" and geno != "0":
            break
    
    samples[sample_id]["id"] = sample_id
    
    return samples[sample_id]

def tally_alleles(samples, alts):
    """ count the alleles used in the genotypes
    
    Args:
        samples: dictionary of parsed sample values for the population.
    
    Returns:
        tallies of how many alleles were seen for each allele.
        e.g. ["0": 10000, "1": 100, "2": 10]
    """
    
    # make sure we have entries for all the alt alleles, even if the count
    # ends up as zero, that way when we later match to the consequences for the alt
    # alleles, the allele counts should be consistent with the consequences
    alts = alts.split(",")
    allele_numbers = ["{}".format(x + 1) for x in range(len(alts))]
    allele_numbers.append("0")
    
    # tally each allele found in the genotypes, by default give every allele a
    # zero count (note that some alleles might not occur in the unaffected
    # DDD parents).
    counts = dict(zip(allele_numbers, [0]* len(allele_numbers)))
    for sample_id in samples:
        alleles = samples[sample_id]["GT"].split("/")
        
        for allele in alleles:
            counts[allele] += 1
    
    return counts

def parse_info(info):
    """ parse the info string from the VCF INFO field.
    
    Args:
        info: VCF INFO field for a single variant.
    
    Returns:
        dictionary of key, values from the info field. Info keys that are flags
        are given True values.
    """
    
    info = info.strip().split(";")
    
    # ensure every variant has the required columns
    info_dict = {"BaseQRankSum": "NA", "GQ_MEAN": "NA", "MQ": "NA", "QD": "NA"}
    for item in info:
        if "=" in item:
            item = item.split("=")
            key = item[0]
            value = item[1]
        else:
            key = item
            value = True
        
        info_dict[key] = value
    
    return info_dict

def check_call_rate(samples):
    """ check that the variant has enough high quality, high depth genotype calls
    
    Erik Minikel said they apply the following filter to exclude variants:
    AN_ADJ > .95*max(AN_ADJ) # at least 95% of samples had a genotype call with GQ >= 20 and DP >= 10
    
    Args:
        samples: dictionary of sample information for each sample, indexed by
            sample ID.
    
    Returns:
        true/false for whether there are sufficient high quality genotype calls
    """
    
    count = 0
    for sample_id in samples:
        
        try:
            if int(samples[sample_id]["DP"]) >= 10 and int(samples[sample_id]["GQ"]) >= 20:
                count += 1
        except ValueError:
            pass
    
    return count > 0.95 * len(samples)

def check_min_alt(samples):
    """ check that the variant has at least one high quality alt allele
    
    Args:
        samples: dictionary of sample information for each sample, indexed by
            sample ID.
    
    Returns:
        true/false for whether the is at least one sample with an alt allele,
        where the genotype call is high quality.
    """
    
    for sample_id in samples:
        if samples[sample_id]["GT"] != "0/0" and int(samples[sample_id]["DP"]) >= 10 and int(samples[sample_id]["GQ"]) >= 20:
            return True
    
    return False

def check_singletons(variant, format, samples, vep, sample_pos, ddd_parents, last_base, min_gq_mean, min_qual, output):
    """ checks to see if any of the alleles at a site are singletons
    
    Args:
        key: tuple of (chrom, pos, ref, alt)
        format: dictionary of format keys, matched to their position in the sample strings
        samples: samples string from VCF line
        vep: dictionary of consequence strings for each allele, indexed by chrom position
        output: File handle for output
    """
    
    key = get_variant_key(variant)
    
    # tally the allele counts across all the genotypes
    samples = parse_samples(samples, format, sample_pos)
    samples = reformat_chrX_genotypes(key, samples, ddd_parents)
    
    if not check_call_rate(samples) or not check_min_alt(samples):
        return
    
    counts = tally_alleles(samples, key[3])
    
    pos = (key[0], key[1])
    is_last_base = pos in last_base
    
    consequences = vep[key[1]]
    
    # remove the reference allele, since we don't have a consequence for that
    total = sum(counts.values())
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
        
        info = parse_info(variant[7])
        if float(variant[5]) < min_qual:
            continue
        
        if "GQ_MEAN" in info and float(info["GQ_MEAN"]) < min_gq_mean:
            continue
        
        consequence_counts[consequence] += 1
        
        # we only want to look at singletons
        if allele_count != 1:
            continue
        
        sample = get_singleton_data(samples)
        
        line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n".format(key[0], key[1], \
            consequence, variant[5], info["BaseQRankSum"], info["GQ_MEAN"], \
            info["MQ"], info["QD"], total, sample["DP"], sample["GQ"], sample["id"])
        
        output.write(line)
        # sys.exit()

def parse_vcf(chrom, ddd_parents, vep, last_base, min_gq_mean, min_qual, output_path):
    """ run through a VCF, counting the alleles, looking for singletons
    """
    
    try:
        output = open(output_path, "w")
    except TypeError:
        output = output_path
    
    # write a file header
    output.write("chrom\tpos\tconsequence\tQUAL\tBaseQRankSum\tGQ_MEAN\tMQ\tQD\tnumber_of_alleles\tsample_DP\tsample_GQ\tsample_id\n")
    
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
        format = get_format_index(variant[8])
        
        check_singletons(variant, format, line[9], vep, sample_pos, ddd_parents, last_base, min_gq_mean, min_qual, output)

def main():
    args = get_options()
    ddd_parents = get_ddd_parents()
    vep = get_vep_annotations(args.chrom)
    last_base = get_last_base_sites(args.last_base_sites)
    
    try:
        parse_vcf(args.chrom, ddd_parents, vep, last_base, args.gq_mean, args.qual, args.singletons,)
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
