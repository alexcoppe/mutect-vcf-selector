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


def check_mutect_variant_with_cosmic(variant_line, cosmic_dictionary, mutectVersion):
    variant = Variant(variant_line)
    filters = variant.filter.split(";")
    if "PASS" in filters:
        return 1
    if 'germline_risk' in filters or mutectVersion == "mutect":
        identifier = variant.chromosome + "," + variant.position + "," + variant.reference + "," + variant.alternative
        cosmic_entry = cosmic_dictionary.get(identifier)
        if cosmic_entry:
            if not cosmic_entry.is_snp:
                return 1
        return 0

def check_mutect_variant_with_clinvar(variant_line, clinvar_dictionary, mutectVersion, clnsig=["Pathogenic"]):
    variant = Variant(variant_line)
    filters = variant.filter.split(";")
    if "PASS" in filters:
        return 1
    clinvar_entry = clinvar_dictionary.get(variant.identifier)
    if clinvar_entry == None:
        return 0
    else:
        if sets.Set(clinvar_entry.clnsig).intersection(clnsig) :
            return 1
    return 0
    

def build_variant_annotation_dictionary(file_name,db):
    f = open(file_name, "r")
    d = {}
    for line in f:
        if not line.startswith("#"):
            if db == "clinvar":
                variant = Clinvar_Variant(line)
            if db == "cosmic":
                variant = Cosmic_Mutation(line)
            d[variant.identifier] = variant
    f.close()
    return d


def check_mutect_clinvar_within_variant(variant_line, mutectVersion, clnsig=["Pathogenic"]):
    clinvar_variant = Clinvar_Variant(variant_line)
    clinvar_filters = clinvar_variant.filter.split(";")
    if "PASS" in clinvar_filters:
        return 1
    if sets.Set(clinvar_variant.clnsig).intersection(clnsig) :
        return 1
    return 0


def check_mutect_cosmic_within_variant(variant_line, mutectVersion):
    cosmic_mutation = Cosmic_Mutation(variant_line)
    if True in [identifier.startswith("COS") for identifier in  cosmic_mutation.ID.split(";")]:
        cosmic_filters = cosmic_mutation.filter.split(";")
        if "PASS" in cosmic_filters:
            return 1
        if not cosmic_mutation.is_snp:
            return 1
    return 0
 
 
 

def main():
    parser = argparse.ArgumentParser(description="Select mutect variants")

    parser.add_argument('-f', '--vcf', action='store', help="vcf file", required=True)
    parser.add_argument('-c', '--cosmic', action='store', help="VCF file with COSMIC data", required=False)
    parser.add_argument('-e', '--header', action='store_true', help="Show VCF's header", required=False)
    parser.add_argument('-l', '--clinvar', action='store', help="VCF file with ClinVar data", required=False)
    parser.add_argument('-s', '--clinical-significance-value', action='store', help="ClinVar's clinical significance value", required=False)
    args = parser.parse_args()

    if args.clinvar:
        clinvar_dictionary = build_variant_annotation_dictionary(args.clinvar, "clinvar")
    if args.cosmic:
        cosmic_dictionary = build_variant_annotation_dictionary(args.cosmic, "cosmic")
    if not args.clinical_significance_value:
        clinical_significance_value = ["Likely_pathogenic", "Pathogenic", "drug_response", 
                "association", "risk_factor", "protective"]
    else:
        clinical_significance_value = args.clinical_significance_value.split(",")

    mutectVersion = "mutect"
    with open(args.vcf, "r") as vcf:
        for line in vcf:
            ok_by_clinvar,ok_by_cosmic = False,False
            if line.startswith("#"):
                if args.header == True:
                    print line[:-1]
                if "##Mutect Version=" in line:
                    mutectVersion = "mutect2"
            else:
                if args.cosmic:
                    if check_mutect_variant_with_cosmic(line, cosmic_dictionary, mutectVersion) == 1:
                        ok_by_cosmic = True
                if args.clinvar:
                    if check_mutect_variant_with_clinvar(line, clinvar_dictionary, mutectVersion, clinical_significance_value) == 1:
                        ok_by_clinvar = True
                        
                if not args.clinvar and not args.cosmic:
                    if check_mutect_clinvar_within_variant(line, mutectVersion, clinical_significance_value) == 1:
                        ok_by_clinvar = True
                    if check_mutect_cosmic_within_variant(line, mutectVersion) == 1:
                        ok_by_cosmic = True
                if ok_by_cosmic == True or ok_by_clinvar == True:
                    print line[:-1]

if __name__ == "__main__":
    main()
