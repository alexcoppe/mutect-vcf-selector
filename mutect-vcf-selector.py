#!/usr/bin/env python

import argparse
import re


class Variant:
    def __init__(self,vcf_line,cosmic=False):
        splitted_vcf_line = vcf_line.split("\t")
        self.chromosome = splitted_vcf_line[0]
        self.position = splitted_vcf_line[1]
        self.ID = splitted_vcf_line[2]
        self.reference = splitted_vcf_line[3]
        self.alternative = splitted_vcf_line[4]
        self.quality = splitted_vcf_line[5]
        self.filter = splitted_vcf_line[6]
        self.info = splitted_vcf_line[7]
        if cosmic == False:
            self.format = splitted_vcf_line[8]
            self.samples = splitted_vcf_line[9:]



class Cosmic_Mutation(Variant):
    def __init__(self,vcf_line):
        Variant.__init__(self, vcf_line,True)
        self.identifier = self.chromosome + "," + self.position + "," + self.reference + "," + self.alternative 
        self.info_dictionary = self.__get_info(self.info)

    def __str__(self):
        return self.chromosome + "," + self.position + "," + self.reference + "," + self.alternative + "," + self.info[:-1]


    def __get_info(self,info_string):
        splitted_info = info_string.split(";")
        if "SNP" in splitted_info:
            self.is_snp = True
        else:
            self.is_snp = False
        p = re.compile('^GENE=(\S+?);')
        m = p.match(info_string)
        self.gene = m.groups()[0]



def check_mutect_variant(variant_line, cosmic_dictionary, mutectVersion):
    splitted_line = variant_line.split()
    filters = splitted_line[6].split(";")
    if "PASS" in filters:
        return 1
    if 'germline_risk' in filters or mutectVersion == "mutect":
        chromosome = splitted_line[0]
        position = splitted_line[1]
        reference = splitted_line[3]
        alternative = splitted_line[4]
        identifier = chromosome + "," + position + "," + reference + "," + alternative
        cosmic_entry = cosmic_dictionary.get(identifier)
        if cosmic_entry:
            if not cosmic_entry.is_snp:
                return 1
        return 0
    

def build_cosmic_dictionary(file_name):
    f = open(file_name, "r")
    d = {}
    for line in f:
        if not line.startswith("#"):
            cosmic_mutation = Cosmic_Mutation(line)
            d[cosmic_mutation.identifier] = cosmic_mutation
    f.close()
    return d


def main():
    parser = argparse.ArgumentParser(description="Select mutect variants")

    parser.add_argument('-f', '--vcf', action='store', help="vcf file", required=True)
    parser.add_argument('-c', '--cosmic', action='store', help="vcf file with cosmic data", required=True)
    args = parser.parse_args()

    cosmic_dictionary = build_cosmic_dictionary(args.cosmic)
    #cosmic_dictionary = []

    mutectVersion = "mutect"
    with open(args.vcf, "r") as vcf:
        for line in vcf:
            if line.startswith("#"):
                if "##Mutect Version=" in line:
                    mutectVersion = "mutect2"
            else:
                if mutectVersion == "mutect":
                    if check_mutect_variant(line, cosmic_dictionary, mutectVersion) == 1:
                        print line[:-1]
                else:
                    if check_mutect_variant(line, cosmic_dictionary, mutectVersion) == 1:
                        print line[:-1]

if __name__ == "__main__":
    main()
