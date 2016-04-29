from Bio import SwissProt
import json
import re
import optparse

def print_all_properties(rec):
    tmp_dict = vars(rec)
    for k in tmp_dict.keys():
        print "%s : %s" % (k, tmp_dict[k])

# Find the pattern that we want using RegEx
def aggregrate_reference_id(regex_pattern, _list):
    match_result = [re.search(regex_pattern, member) for member in _list]
    return [m.group(0) for m in match_result if m is not None]


# In this section, there are also many database's names where this protein was referenced. But this time, we just
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
            # In Ensembl, we concerned the name of proein that begins with ENST_xxxxxxxxxxx, ENSP_xxxxxxxxxxx and ENSG_xxxxxxxxxxx
            o['Ensembl'].extend(aggregrate_reference_id('ENS[T|P|G]\d{10,12}', ref[1:]))
    return o


# In CC section, there are explaination topics; "FUNCTION" is one of them. Therefore, after we could analyze with
# text mining technique, these topics, including "FUNCTION" may be mined.
def parse_comments(cc):
    return [comment.replace('FUNCTION: ', '').replace('\n', '') for comment in cc if comment.startswith('FUNCTION:')]


def help_func():
    p = optparse.OptionParser()
    p.add_option('-i', '--input', default=True, metavar="FILE", help="read input from FILE of Swiss-Prot Database text format (.dat)")
    p.add_option('-o', '--output', default=True, metavar="FILE", help="write output to FILE in JSON format")
    return p.parse_args()

def main():

    # program's option
    options, arguments = help_func()

    ### field for database ###
    seq = ''
    cross_ref = ''
    seq_len = ''
    cc = ''
    accessions = ''
    record_name = ''
    record_class = ''

    source_file = options.input
    output_file = options.output
    target_specie = 'Homo sapiens (Human)'

    with open(source_file, 'rU') as handle_in, open(output_file, 'w') as handle_out:
        for record in SwissProt.parse(handle_in):
            if target_specie in record.organism:
                record_name = record.entry_name
                seq = record.sequence
                cross_ref = parse_cross_reference(record.cross_references)
                seq_len = record.sequence_length
                cc = parse_comments(record.comments)
                accessions = record.accessions
                record_class = record.data_class
                # convert to python's dictionary > JSON object
                out_dict = dict()
                out_dict['sequence'] = seq
                out_dict['cross_reference'] = cross_ref
                out_dict['length'] = seq_len
                out_dict['comment'] = cc
                out_dict['accessions'] = accessions
                out_dict['name'] = record_name
                out_dict['proved'] = record_class
                json_record = json.dumps(out_dict, sort_keys=True, separators=(',', ':'))
                handle_out.write(json_record + '\n')
if __name__ == '__main__':
    main()