import re
import pandas as pd
import argparse
import glob
import os


class breseq:
    def __init__(self, gd_dir, output_file):
        self.gd_dir = gd_dir      # directory with annotated gd files with "_annotated.gd" suffix
        self.output_file = output_file              # file for results
        self.clone_list = list()                    # list of clones
        self.mutation_table = self.create_out_table()   # table containing mutations
        self.aa_conv = self.aa_convertion_dict()        # dictionary to convert 1-letter code to 3-letter code
        self.parsed = list()                            # ids of parsed evidences
        self.ra_list = dict()                           # dictionary of ra ids for reference nucleotides

    @staticmethod
    def create_out_table():
        return pd.DataFrame(columns=[
            'Clone',
            'Reactor',
            'Time',
            'Gene_ID',
            'P/A',
            'AA_Mutation',
            'In_unevolved',
            'Annotation',
            'MIC',
            'Comment',
            '-',
            'Genome',
            'DNA_Mutation',
            'Effect'
        ])

    @staticmethod
    def aa_convertion_dict():
        return {'A': 'Ala',
                'B': 'Asx',
                'C': 'Cys',
                'D': 'Asp',
                'E': 'Glu',
                'F': 'Phe',
                'G': 'Gly',
                'H': 'His',
                'I': 'Ile',
                'K': 'Lys',
                'L': 'Leu',
                'M': 'Met',
                'N': 'Asn',
                'P': 'Pro',
                'Q': 'Gln',
                'R': 'Arg',
                'S': 'Ser',
                'T': 'Thr',
                'V': 'Val',
                'W': 'Trp',
                'X': 'Xaa',
                'Y': 'Tyr',
                'Z': 'Glx',
                '*': '*'
                }

    def read_clone_list(self):              # read file with clones
        self.clone_list = glob.glob(os.path.join(self.gd_dir, '*_annotated.gd'))

    def read_annotated_output(self):
        for clone in self.clone_list:
            out_file = open(clone, 'r')
            print('Analyzing: ' + clone)

            for line in out_file.readlines():
                line = line.rstrip()
                mutation = line.split('\t')
                if mutation[0] == 'DEL':
                    self.analyze_DEL(line, clone)
                elif mutation[0] == 'INS':
                    self.analyze_INS(line, clone)
                elif mutation[0] == 'SNP':
                    if mutation[8][:2] == 'aa':
                        self.analyze_SNP_c(line, clone)
                    else:
                        self.analyze_SNP_nc(line, clone)
                elif mutation[0] == 'MC':
                    self.analyze_MC(line, clone)
                elif mutation[0] == 'JC' and mutation[1] not in self.parsed and mutation[1]:
                    reject = re.search('reject=', line)
                    if reject is None:
                        self.analyze_JC(line, clone)
                elif mutation[0] == 'MOB':
                    self.analyze_MOB(line, clone)
                elif mutation[0] == 'RA' and mutation[1] in self.ra_list:
                    self.mutation_table.loc[self.ra_list[mutation[1]], 'DNA_Mutation'] = mutation[4] + ' ' + mutation[6] + \
                                                                         '->' + mutation[7]
                else:
                    continue

            self.ra_list = dict()       # clean RA ids for the next clone
            self.parsed = list()       # clean JC ids for the next clone

    def analyze_MOB(self, line, clone):
        temp_table = self.create_out_table()
        temp_table = temp_table.append(pd.Series(['NA'] * len(temp_table.columns), index=temp_table.columns),
                                       ignore_index=True)
        mut = line.split('\t')
        self.parsed.extend(mut[2].split(','))
        temp_table['Index'] = clone + '_' + mut[2] + '_' + mut[4]
        temp_table = temp_table.set_index('Index')

        gid = re.search('gene_name=([^\t]+)', line)
        aa = re.search('gene_position=([^\t]+)', line)
        anno = re.search('gene_product=([^\t]+)', line)
        eff = re.search('mutation_category=([^\t]+)', line)
        pstart = re.search('position_start=([^\t]+)', line)
        pend = re.search('position_end=([^\t]+)', line)

        temp_table['Clone'] = clone
        temp_table['Reactor'] = ''
        temp_table['Time'] = ''
        temp_table['Gene_ID'] = gid.group(1)
        temp_table['P/A'] = ''
        if aa is None:
            temp_table['AA_Mutation'] = ''
        else:
            temp_table['AA_Mutation'] = aa.group(1)
        temp_table['In_unevolved'] = ''
        temp_table['Annotation'] = anno.group(1)
        temp_table['MIC'] = ''
        temp_table['Comment'] = ''  # write frequency here
        temp_table['-'] = 'MOB'  # write mutation type here
        temp_table['Genome'] = mut[3]
        temp_table['DNA_Mutation'] = mut[5] + ' ins ' + pstart.group(1) + '-' + pend.group(1)
        temp_table['Effect'] = eff.group(1)
        self.mutation_table = self.mutation_table.append(temp_table)

    def analyze_INS(self, line, clone):
        temp_table = self.create_out_table()
        temp_table = temp_table.append(pd.Series(['NA']*len(temp_table.columns), index=temp_table.columns),
                                       ignore_index=True)

        mut = line.split('\t')
        temp_table['Index'] = clone + '_' + mut[2] + '_' + mut[4]
        temp_table = temp_table.set_index('Index')

        self.parsed.extend(mut[2].split(','))
        gid = re.search('gene_name=([^\t]+)', line)
        aa = re.search('gene_position=([^\t]+)', line)
        anno = re.search('gene_product=([^\t]+)', line)
        eff = re.search('mutation_category=([^\t]+)', line)

        temp_table['Clone'] = clone
        temp_table['Reactor'] = ''
        temp_table['Time'] = ''
        temp_table['Gene_ID'] = gid.group(1)
        temp_table['P/A'] = ''
        temp_table['AA_Mutation'] = aa.group(1)
        temp_table['In_unevolved'] = ''
        temp_table['Annotation'] = anno.group(1)
        temp_table['MIC'] = ''
        temp_table['Comment'] = ''        # write frequency here
        temp_table['-'] = 'INS'                      # write mutation type here
        temp_table['Genome'] = mut[3]
        temp_table['DNA_Mutation'] = mut[4] + ' ins ' + mut[5]
        temp_table['Effect'] = eff.group(1)
        self.mutation_table = self.mutation_table.append(temp_table)

    def analyze_DEL(self, line, clone):
        temp_table = self.create_out_table()
        temp_table = temp_table.append(pd.Series(['NA'] * len(temp_table.columns), index=temp_table.columns),
                                       ignore_index=True)
        mut = line.split('\t')
        temp_table['Index'] = clone + '_' + mut[2] + '_' + mut[4]
        temp_table = temp_table.set_index('Index')

        self.parsed.extend(mut[2].split(','))
        gid = re.search('gene_name=([^\t]+)', line)
        aa = re.search('gene_position=([^\t]+)', line)
        anno = re.search('gene_product=([^\t]+)', line)
        eff = re.search('mutation_category=([^\t]+)', line)

        temp_table['Clone'] = clone
        temp_table['Reactor'] = ''
        temp_table['Time'] = ''
        temp_table['Gene_ID'] = gid.group(1)
        temp_table['P/A'] = ''
        if aa is None:
            temp_table['AA_Mutation'] = ''
        else:
            temp_table['AA_Mutation'] = aa.group(1)
        temp_table['In_unevolved'] = ''
        temp_table['Annotation'] = anno.group(1)
        temp_table['MIC'] = ''
        temp_table['Comment'] = ''        # write frequency here
        temp_table['-'] = 'DEL'                      # write mutation type here
        temp_table['Genome'] = mut[3]
        temp_table['DNA_Mutation'] = mut[4] + ' del ' + mut[5] + 'nt'
        temp_table['Effect'] = eff.group(1)
        self.mutation_table = self.mutation_table.append(temp_table)

    def analyze_SNP_nc(self, line, clone):
        temp_table = self.create_out_table()
        temp_table = temp_table.append(pd.Series(['NA'] * len(temp_table.columns), index=temp_table.columns),
                                       ignore_index=True)
        mut = line.split('\t')
        temp_table['Index'] = clone + '_' + mut[2] + '_' + mut[4]
        temp_table = temp_table.set_index('Index')

        gid = re.search('gene_name=([^\t]+)', line)
        aa = re.search('gene_position=([^\t]+)', line)
        anno = re.search('gene_product=([^\t]+)', line)
        eff = re.search('snp_type=([^\t]+)', line)

        temp_table['Clone'] = clone
        temp_table['Reactor'] = ''
        temp_table['Time'] = ''
        temp_table['Gene_ID'] = gid.group(1)
        temp_table['P/A'] = ''
        temp_table['AA_Mutation'] = aa.group(1)
        temp_table['In_unevolved'] = ''
        temp_table['Annotation'] = anno.group(1)
        temp_table['MIC'] = ''
        temp_table['Comment'] = ''
        temp_table['-'] = 'SNP'                      # write mutation type here
        temp_table['Genome'] = mut[3]
        temp_table['DNA_Mutation'] = ''
        if eff is None:
            temp_table['Effect'] = ''
        else:
            temp_table['Effect'] = eff.group(1)
        temp_table['Index'] = clone + '_' + mut[2] + '_' + mut[4]           # set index to find the line later
        temp_table = temp_table.set_index('Index')
        self.ra_list[mut[2]] = clone + '_' + mut[2] + '_' + mut[4]          # remember to look in RA line to find dna mutation

        self.mutation_table = self.mutation_table.append(temp_table)

    def analyze_SNP_c(self, line, clone):
        temp_table = self.create_out_table()
        temp_table = temp_table.append(pd.Series(['NA'] * len(temp_table.columns), index=temp_table.columns),
                                       ignore_index=True)
        mut = line.split('\t')
        temp_table['Index'] = clone + '_' + mut[2] + '_' + mut[4]
        temp_table = temp_table.set_index('Index')

        gid = re.search('gene_name=([^\t]+)', line)
        aa_r = re.search('aa_ref_seq=([^\t]+)', line)
        aa_a = re.search('aa_new_seq=([^\t]+)', line)
        aa_pos = re.search('aa_position=([^\t]+)', line)
        anno = re.search('gene_product=([^\t]+)', line)
        eff = re.search('snp_type=([^\t]+)', line)
        codon_r = re.search('codon_ref_seq=([^\t]+)', line)
        codon_a = re.search('codon_new_seq=([^\t]+)', line)
        codon_p = re.search('codon_position=(\d+)', line)

        temp_table['Clone'] = clone
        temp_table['Reactor'] = ''
        temp_table['Time'] = ''
        temp_table['Gene_ID'] = gid.group(1)
        temp_table['P/A'] = ''
        temp_table['AA_Mutation'] = self.aa_conv[aa_r.group(1)] + aa_pos.group(1) + self.aa_conv[aa_a.group(1)]
        temp_table['In_unevolved'] = ''
        temp_table['Annotation'] = anno.group(1)
        temp_table['MIC'] = ''
        temp_table['Comment'] = ''
        temp_table['-'] = 'SNP'                      # write mutation type here
        temp_table['Genome'] = mut[3]
        temp_table['DNA_Mutation'] = ''
        temp_table['Effect'] = eff.group(1)
        self.ra_list[mut[2]] = clone + '_' + mut[2] + '_' + mut[4]  # remember to look in RA line to find dna mutation
        self.mutation_table = self.mutation_table.append(temp_table)

    def analyze_MC(self, line, clone):
        temp_table = self.create_out_table()
        temp_table = temp_table.append(pd.Series(['NA'] * len(temp_table.columns), index=temp_table.columns),
                                       ignore_index=True)
        mut = line.split('\t')
        temp_table['Index'] = clone + '_' + mut[2] + '_' + mut[4]
        temp_table = temp_table.set_index('Index')

        gid = re.search('gene_name=([^\t]+)', line)

        temp_table['Clone'] = clone
        temp_table['Reactor'] = ''
        temp_table['Time'] = ''
        temp_table['Gene_ID'] = gid.group(1)
        temp_table['P/A'] = ''
        temp_table['AA_Mutation'] = ''
        temp_table['In_unevolved'] = ''
        temp_table['Annotation'] = ''
        temp_table['MIC'] = ''
        temp_table['Comment'] = ''
        temp_table['-'] = 'MC'
        temp_table['Genome'] = mut[3]
        temp_table['DNA_Mutation'] = mut[4] + '-' + mut[5]
        temp_table['Effect'] = 'region deletion'
        self.mutation_table = self.mutation_table.append(temp_table)

    def analyze_JC(self, line, clone):
        temp_table = self.create_out_table()
        temp_table = temp_table.append(pd.Series(['NA'] * len(temp_table.columns), index=temp_table.columns),
                                       ignore_index=True)
        mut = line.split('\t')
        temp_table['Index'] = clone + '_' + mut[2] + '_' + mut[4]
        temp_table = temp_table.set_index('Index')

        freq = re.search('polymorphism_frequency=([^\t]+)', line)

        temp_table['Clone'] = clone
        temp_table['Reactor'] = ''
        temp_table['Time'] = ''
        temp_table['Gene_ID'] = ''
        temp_table['P/A'] = ''
        temp_table['AA_Mutation'] = ''
        temp_table['In_unevolved'] = ''
        temp_table['Annotation'] = ''
        temp_table['MIC'] = ''
        temp_table['Comment'] = freq.group(1)  # write frequency here
        temp_table['-'] = 'JC'
        temp_table['Genome'] = mut[3] + '<>' + mut[6]
        temp_table['DNA_Mutation'] = mut[4] + '<>' + mut[7]
        temp_table['Effect'] = 'junction'
        self.mutation_table = self.mutation_table.append(temp_table)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert breseq data to table format')
    parser.add_argument('-g', '--gd_dir', help='Directory with GD annotated files')
    parser.add_argument('-o', '--output', help='Output file')
    args = parser.parse_args()

    gd_dir = args.gd_dir
    out_file = args.output

    x = breseq(gd_dir, out_file)
    x.read_clone_list()
    x.read_annotated_output()
    x.mutation_table.to_csv(out_file, sep='\t', index=False)
