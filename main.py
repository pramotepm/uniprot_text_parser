from Bio import SwissProt
import json
import re
import optparse
import sys
from cc import Comments


def print_all_properties(rec):
    tmp_dict = vars(rec)
    for k in tmp_dict.keys():
        print "%s : %s\n" % (k, tmp_dict[k])


# Find the pattern that we want using RegEx
def aggregrate_reference_id(regex_pattern, _list):
    match_result = [re.search(regex_pattern, member) for member in _list]
    return [m.group(0) for m in match_result if m is not None]


# In this function, there are many database where this protein was referenced. But, we just
# concerned the reference that contained in PDB, RefSeq and Ensembl.
def parse_cross_reference(cr):
    o = dict()
    o['PDB'] = list()
    o['RefSeq'] = list()
    o['Ensembl'] = list()
    for ref in cr:
        if 'PDB' in ref:
            o['PDB'].append(ref[1])
        elif 'RefSeq' in ref:
            # In RefSeq, we concerned the name of proein that begins with NP_xxxxxx and NP_xxxxxx
            o['RefSeq'].extend(aggregrate_reference_id('N[M|P]_\d{5,7}\.\d', ref[1:]))
        elif 'Ensembl' in ref:
            # In Ensembl, we concerned the name of proein that begins with ENST_xxxxxxxxxxx,
            # ENSP_xxxxxxxxxxx and ENSG_xxxxxxxxxxx
            o['Ensembl'].extend(aggregrate_reference_id('ENS[T|P|G]\d{10,12}', ref[1:]))
    return o


def usage():
    p = optparse.OptionParser()
    p.add_option('-i', '--input', default=False, metavar="FILE",
                 help="read input from FILE of Swiss-Prot Database text format (.dat)")
    p.add_option('-o', '--output', default=False, metavar="FILE", help="write output to FILE in JSON format")
    p.add_option('-t', '--test', default=False, metavar="INT", help="test running by print n line(s) to standard output, n is positive"
                                                     "integer")
    return p.parse_args()


def read_record(input_file):
    target_specie = 'Homo sapiens (Human)'
    with open(input_file, 'r') as handle_in:
        for record in SwissProt.parse(handle_in):
            if target_specie in record.organism:
                yield record


def main():
    # program's option
    options, arguments = usage()

    input_file = options.input
    output_file = options.output
    n_test_print = options.test if type(options.test) is bool else int(options.test)
    print n_test_print

    with (open(output_file, 'w') if output_file else sys.stdout) as handle_out:
        for record in read_record(input_file):
            # put the parsed value from each record to dictionary,
            # and prepare to parse again in JSON object format
            out_dict = dict()
            out_dict['name'] = record.entry_name

            # for logging, we print name
            print 'Parse: ' + record.entry_name

            out_dict['sequence'] = record.sequence
            out_dict['cross_reference'] = parse_cross_reference(record.cross_references)
            out_dict['length'] = record.sequence_length

            # In the CC section, there are a lot of information,
            # so we treat they as a class of "Comments" class
            cc = Comments(record.comments)
            out_dict['comments'] = dict()
            out_dict['comments']['function'] = cc.get_topic_function()
            out_dict['comments']['alternative_products'] = cc.get_topic_alternative_products()

            out_dict['accessions'] = record.accessions
            out_dict['proved'] = record.data_class

            # parse into JSON format
            json_record = json.dumps(out_dict, sort_keys=True, separators=(',', ':'))

            # write out to file
            if output_file is not None:
                handle_out.write(json_record + '\n')
            else:
                print (json_record + '\n')

            if n_test_print:
                if n_test_print-1 > 0:
                    n_test_print -= 1
                else:
                    break

def main_test():
    options, arguments = usage()
    input_file = options.input
    for record in read_record(input_file):
        # do something to process each record
        print_all_properties(record)
        break

if __name__ == '__main__':
    main()
    # main_test()