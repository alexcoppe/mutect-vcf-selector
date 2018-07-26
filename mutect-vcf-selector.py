#!/usr/bin/env python

import argparse
import re
import sets


class Variant:
    def __init__(self,vcf_line,untill_filter=False):
        splitted_vcf_line = vcf_line.split("\t")
        self.chromosome = splitted_vcf_line[0]
        self.position = splitted_vcf_line[1]
        self.ID = splitted_vcf_line[2]
        self.reference = splitted_vcf_line[3]
        self.alternative = splitted_vcf_line[4]
        self.quality = splitted_vcf_line[5]
        self.filter = splitted_vcf_line[6]
        self.info = splitted_vcf_line[7]
        if untill_filter == False:
            self.format = splitted_vcf_line[8]
            self.samples = splitted_vcf_line[9:]
        self.identifier = self.chromosome + "," + self.position + "," + self.reference + "," + self.alternative 
    def __str__(self):
        return self.chromosome + "," + self.position + "," + self.reference + "," + self.alternative + "," + self.info[:-1]



class Cosmic_Mutation(Variant):
    def __init__(self,vcf_line):
        Variant.__init__(self, vcf_line,True)
        self.is_snp,self.gene = self.__get_info(self.info)

    def __get_info(self,info_string):
        gene = None
        splitted_info = info_string.split(";")
        if "SNP" in splitted_info:
            is_snp = True
        else:
            is_snp = False
        p = re.compile('^GENE=(\S+?);')
        m = p.match(info_string)
        if m:
            gene = m.groups()[0]
        return is_snp,gene




class Clinvar_Variant(Variant):
    def __init__(self,vcf_line):
        Variant.__init__(self, vcf_line,True)
        self.info_dictionary = self.__get_info(self.info)
        self.clnsig = self.__get_info(self.info)

    def __get_info(self,info_string):
        splitted_info = info_string.split(";")
        for element in splitted_info:
            p = re.compile('CLNSIG=(\S+)')
            m = p.match(element)
            if m:
                return  m.groups()[0].split('/')
        return []



def check_clinvar(variant_line, clnsig=["Pathogenic"]):
    clinvar_variant = Clinvar_Variant(variant_line)
    clinvar_filters = clinvar_variant.filter.split(";")
    if "PASS" in clinvar_filters:
        return 1
    if sets.Set(clinvar_variant.clnsig).intersection(clnsig) :
        return 1
    return 0


def check_cosmic(variant_line):
    cosmic_mutation = Cosmic_Mutation(variant_line)
    if True in [identifier.startswith("COS") for identifier in  cosmic_mutation.ID.split(";")]:
        cosmic_filters = cosmic_mutation.filter.split(";")
        if "PASS" in cosmic_filters:
            return 1
        if not cosmic_mutation.is_snp:
            return 1
    return 0
 
 
def make_cancer_gene_census_dictionary(f):
    genes = []
    with open(f, "r") as cancer_gene_census_dictionary_file:
        for line in cancer_gene_census_dictionary_file:
            genes.append(line.split(",")[0])
    return set(genes)

 

def main():
    parser = argparse.ArgumentParser(description="Select mutect variants")

    parser.add_argument('-f', '--vcf', action='store', help="vcf file", required=True)
    parser.add_argument('-e', '--header', action='store_true', help="Show VCF's header", required=False)
    parser.add_argument('-c', '--cgc', action='store', help="The Cancer Gene Census (CGC) gene CSV file", required=True)
    args = parser.parse_args()

    if not args.clinical_significance_value:
        clinical_significance_value = ["Likely_pathogenic", "4", "Pathogenic", "5", "drug_response", "6", 
                "association", "risk_factor", "protective"]
    else:
        clinical_significance_value = args.clinical_significance_value.split(",")


    cancer_gene_census_dictionary = make_cancer_gene_census_dictionary(args.cgc)

    with open(args.vcf, "r") as vcf:
        for line in vcf:
            ok_by_clinvar,ok_by_cosmic = False,False
            if line.startswith("#"):
                if args.header == True:
                    print line[:-1]
            else:
                variant = Variant(line)
                variant_gene_name = variant.info.split("|")[3]
                if check_clinvar(line, clinical_significance_value) == 1:
                    ok_by_clinvar = True
                if check_cosmic(line) == 1:
                    ok_by_cosmic = True
                if variant.filter == "PASS":
                    print line[:-1]
                elif ok_by_cosmic == True or ok_by_clinvar == True:
                    if variant_gene_name in cancer_gene_census_dictionary:
                        print line[:-1]

if __name__ == "__main__":
    main()
