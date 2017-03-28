from config import Config
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from myhash import get_digest_md5
import os


def main():
    config = Config()
    db = config.get_connection_collection()
    for doc in db.find({}, {'uniprot_id':1, 'isoform_product':1}):
        for doc2 in doc['isoform_product']:
            if 'seq' in doc2 and (True if 'note' not in doc2 else doc2['note'] != 'External'):
                if doc2['digest'] != get_digest_md5(doc2['seq']):
                    exit(1)
                seq = Seq.Seq(doc2['seq'])
                is_ref = 'is_ref' in doc2
                isoform_id = '%s.%s' % (doc['uniprot_id'], doc2['isoform_id'].split('-')[1])
                seq_record = SeqRecord(seq, isoform_id, '', '')
                with open(os.path.join(config.get_fasta_dir_path(), ''.join([isoform_id, '.ref' if is_ref else '.var'])), 'w') as f_out:
                    SeqIO.write(seq_record, f_out, 'fasta')
                f_out.close()

if __name__ == '__main__':
    main()
