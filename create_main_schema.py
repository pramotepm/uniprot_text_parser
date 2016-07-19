from Bio import SwissProt
import json
import re
import optparse
import sys
from collections import deque
from cc import Comments
from cc import ParsingException
from myhash import get_digest_md5
from config import Config


def print_all_properties(rec):
    tmp_dict = vars(rec)
    for k in tmp_dict.keys():
        print "%s : %s\n" % (k, tmp_dict[k])
    print "================"


# Find the pattern that we want using RegEx
def aggregate_reference_id(regex_pattern, _list):
    tmp = []
    match_result = [re.search(regex_pattern, member) for member in _list]
    for m in match_result:
        tmp.append(m.group(1))
        if m.group(3) is not None:
            tmp.append(m.group(3))
    return tmp


def _foo(con1, con2, con3, in_list):
    if len(in_list) > 4:
        err_mesg = "Cross database references are out of condition"
        raise ParsingException(err_mesg)
    d = dict()
    in_list = deque(in_list)
    while len(in_list) != 0:
        t = in_list.popleft()
        if t.startswith(con1):
            d['transcript_id'] = t
        elif t.startswith(con2):
            d['protein_id'] = t
        elif t.startswith(con3):
            d['gene_id'] = t
        else:
            d['isoform_id'] = t
    if 'transcript_id' not in d and 'protein_id' not in d:
        return None
    else:
        return d


# In this function, there are many database where this protein was referenced. But, we just
# concerned the reference that contained in PDB, RefSeq and Ensembl.
def parse_cross_reference(cr):
    o = dict()
    o['pdb'] = list()
    o['refseq'] = list()
    o['ensembl'] = list()
    o['go'] = list()
    for ref in cr:
        if 'PDB' in ref:
            o['pdb'].append(ref[1])
        elif 'Ensembl' in ref:
            # In Ensembl, we concerned the name of proein that begins with ENST_xxxxxxxxxxx,
            # ENSP_xxxxxxxxxxx and ENSG_xxxxxxxxxxx
            ensembl_list = _foo('ENST', 'ENSP', 'ENSG',
                                aggregate_reference_id('(ENS[T|P|G]\d+)(. \[(.*)\])?', ref[1:]))
            o['ensembl'].append(ensembl_list)
            # print o['Ensembl']
        elif 'RefSeq' in ref:
            # In RefSeq, we concerned the name of proein that begins with NP_xxxxxx and NP_xxxxxx
            refseq_list = _foo('NM', 'NP', 'NG',
                               aggregate_reference_id('([A-Z][A-Z]_\d+\.\d)(. \[(.*)\])?', ref[1:]))
            if refseq_list is not None:
                o['refseq'].append(refseq_list)
        elif 'GO' in ref:
            a, b = ref[2].split(':', 1)
            o['go'].append({'id':ref[1], 'domain':a.lower(), 'def':b, 'flag':ref[3]})
    return o


def usage():
    p = optparse.OptionParser()
    p.add_option('-t', '--test', default=False, metavar="INT", help="test running by print n line(s) to standard "
                                                                    "output, n is positive integer", type="int")
    p.add_option('--all-property', action="store_const", default=False, const=True, dest="verbose",
                 help="use default biopython parser to print all information of record(s)")
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

    config = Config()
    input_file = config.get_sprot_path()
    n_test_print = options.test if type(options.test) is bool else options.test
    db = config.get_connection_collection()

    for record in read_record(input_file):
        if options.verbose:
            print_all_properties(record)
        else:
            # put the parsed value from each record to dictionary,
            # and prepare to parse again in JSON object format
            out_dict = dict()
            out_dict['uniprot_id'] = record.entry_name

            # for logging, print name
            # print 'Parsing from SwissProt: ' + record.entry_name

            # Add Ensembl, PDB, RefSeq field to DB
            out_dict.update(parse_cross_reference(record.cross_references))

            # In the CC section, there are a lot of information,
            # so we treat them as a class of "Comments" class
            cc = Comments(record.comments)
            out_dict['function'] = cc.get_topic_function()
            out_dict.update(cc.get_topic_alternative_products())
            out_dict['isoform_product'] = cc.get_topic_alternative_product_isoform_product()
            # In case of there is no isoform, use accession[0] as its isoform ID
            if out_dict['isoform_product'][0]['isoform_id'] == '':
                out_dict['isoform_product'][0]['isoform_id'] = record.accessions[0] + "-1"
            out_dict['isoform_product'][0]['seq'] = record.sequence
            out_dict['isoform_product'][0]['length'] = record.sequence_length
            # Set the first isoform as reference isoform (REF)
            out_dict['isoform_product'][0]['is_ref'] = True
            out_dict['isoform_product'][0]['digest'] = get_digest_md5(record.sequence)

            out_dict['uniprot_accession'] = record.accessions
            out_dict['prove'] = record.data_class

            # parse into JSON format
            if n_test_print == False:
                json_record = json.dumps(out_dict, sort_keys=True, separators=(',', ':'))

            # insert output to db

            db.insert_one(json.loads(json_record))

        # how many record will printed
        if n_test_print:
            if n_test_print-1 > 0:
                n_test_print -= 1
            else:
                break
    # print "Run complete"


if __name__ == '__main__':
    main()
