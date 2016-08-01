import re
import sys
from myhash import get_digest_md5
from config import Config
from pymongo import ASCENDING

def read_f(path):
    with open(path, 'rU') as f_in:
        for line in f_in:
            yield line.strip().rstrip()


def chk_seq_length_correct(seq, length):
    if len(seq) != length:
        sys.exit()


def extract_info(lines):
    read_seq_state = False
    seq_frag = []
    is_match_NP = False
    for line in lines:
        if line.startswith('LOCUS'):
            NP_regex = re.search('(NP_\d+)', line)
            length_regex = re.search('(\d+) aa', line)
            if NP_regex is not None:
                NP_rec = NP_regex.group(1)
                length = int(length_regex.group(1))
                is_match_NP = True
            else:
                is_match_NP = False
        elif line.startswith('DBSOURCE'):
            try:
                NM_rec = re.search('(NM_\d+.\d+)', line).group(1)
            except:
                pass
        elif line.startswith('VERSION'):
            try:
                NP_ver_rec = re.search('(NP_\d+.\d+)', line).group(1)
            except:
                pass
        elif line.startswith('//'):
            seq = ''.join(seq_frag).upper()
            if is_match_NP:
                chk_seq_length_correct(seq, length)
                yield NP_rec, NP_ver_rec, NM_rec, seq
            # set to default value
            read_seq_state = False
            seq_frag = []
        elif read_seq_state:
            seq_frag.extend(map(lambda x: x.group(1), list(re.finditer('([a-z]+)', line))))
        elif line.startswith('ORIGIN'):
            read_seq_state = True


def compare_to_db(NP, NP_rec, NM_rec, seq):
    cursor = db.find({'isoform_product.digest': get_digest_md5(seq)})
    for doc in cursor:
        uniprot_id = doc['uniprot_id']
        isoform_id = find_isoform_id_for_seq(doc, seq)
        chk_insert = not is_exists(doc, isoform_id, NP_rec, NM_rec) and (add_to_matched_NP_NM(doc, uniprot_id, isoform_id, NP_rec, NM_rec) or add_new_NP_NM(doc, uniprot_id, isoform_id, NP_rec, NM_rec))
        if chk_insert:
            print 'Add RefSeq %s -> SwissProt %s(%s)' % (NP, uniprot_id, isoform_id)
        # else:
            # print 'Collision', uniprot_id, isoform_id, NP_rec, NM_rec, seq


def is_exists(doc, isoform_id, NP_rec, NM_rec):
    for db_isoform_id, db_NP_rec, db_NM_rec in [ (i['isoform_id'], i['protein_id'], i['transcript_id']) for i in doc['refseq'] if 'protein_id' in i and 'transcript_id' in i and 'isoform_id' in i ]:
        if isoform_id == db_isoform_id and NP_rec == db_NP_rec and NM_rec == db_NM_rec:
            return True
    return False


def add_to_matched_NP_NM(doc, uniprot_id, isoform_id, NP_rec, NM_rec):
    if 'refseq' in doc:
        for NP_rec_db, NM_rec_db in [(i['protein_id'], i['transcript_id']) for i in doc['refseq'] if 'isoform_id' not in i and 'protein_id' in i and 'transcript_id' in i]:
            if NM_rec == NM_rec_db and NP_rec == NP_rec_db:
                return db.update_one({'uniprot_id': uniprot_id, 'refseq.protein_id': NP_rec, 'refseq.transcript_id': NM_rec}, {'$set': {'refseq.$.isoform_id': isoform_id}}).modified_count == 1
    return False


def add_new_NP_NM(doc, uniprot_id, isoform_id, NP_rec, NM_rec):
    if 'refseq' in doc and isoform_id not in [ i['isoform_id'] for i in doc['refseq'] if 'isoform' in i ]:
        return db.update_one({'uniprot_id':uniprot_id}, {'$push':{'refseq':{'transcript_id':NM_rec, 'protein_id':NP_rec, 'isoform_id':isoform_id}}}).modified_count == 1
    return False


def find_isoform_id_for_seq(doc, seq):
    for r in doc['isoform_product']:
        chk_note = (r['note'] != 'Not described' and r['note'] != 'External') if 'note' in r else True
        if chk_note:
            # print doc['uniprot_id']
            if r['seq'] == seq:
                return r['isoform_id']
    return False


def main():
    config = Config()
    global db
    db = config.get_connection_collection()
    db.create_index([("isoform_product.digest", ASCENDING)], name="temp1")
    db.create_index([("uniprot_id", ASCENDING), ("refseq.protein_id", ASCENDING), ("refseq.transcript_id", ASCENDING)], name="temp2")
    for NP, NP_ver, NM, seq in extract_info(read_f(config.get_refseq_protein_path())):
        compare_to_db(NP, NP_ver, NM, seq)
    db.drop_index("temp1")
    db.drop_index("temp2")

if __name__ == "__main__":
    main()
