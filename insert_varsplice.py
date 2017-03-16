from Bio import SeqIO
from myhash import get_digest_md5
from config import Config
from pymongo import ASCENDING

def find_human_protein(record):
    foo, isoform_id, uniprot_id = str(record.id).split('|')
    if uniprot_id.endswith("_HUMAN"):
        return (isoform_id, uniprot_id)
    else:
        return None


def main():
    config = Config()
    db = config.get_connection_collection()
    db.create_index([("isoform_product.isoform_id", ASCENDING), ("uniprot_id", ASCENDING)], name="temp1", unique=True)
    with open(config.get_varsplic_path(), 'rU') as file_input:
        for record in SeqIO.parse(file_input, "fasta"):
            ids = find_human_protein(record)
            if ids is not None:
                isoform_id, uniprot_id = ids
                seq = str(record.seq)
                length = len(seq)
                # print uniprot_id   # print for logging
                digest = get_digest_md5(seq)
                db.update_one({"uniprot_id":uniprot_id, "isoform_product.isoform_id": isoform_id}, {"$set": {"isoform_product.$.length": length, "isoform_product.$.seq": seq, "isoform_product.$.digest": digest}})
    db.drop_index("temp1")

if __name__ == '__main__':
    main()
