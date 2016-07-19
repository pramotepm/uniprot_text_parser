from config import Config
from myhash import get_digest_md5

def main():
    config = Config()
    db = config.get_connection_collection()
    cursor = db.find({"isoform_product.note":"External", "isoform_product.seq":None})
    for document in cursor:
        uniprot_id = document['uniprot_id']
        isoform_ids = [i['isoform_id'] for i in document['isoform_product'] if i['note'] == "External"]
        for isoform_id in isoform_ids:
            seq_find_query = "{isoform_product.isoform_id:\"%s\"}, {isoform_product:{$elemMatch:{isoform_id:\"%s\"}}}" % (
            isoform_id, isoform_id)
            cursor2 = db.find({"isoform_product.isoform_id": isoform_id}, {"isoform_product": {"$elemMatch": {"isoform_id": isoform_id}}})
            foo = [document2['isoform_product'][0]['seq'] for document2 in cursor2 if 'seq' in document2['isoform_product'][0] and document2['isoform_product'][0]['note'] != 'External']
            if len(foo) != 1:
                print "Ambiguous External Sequence ! : %s %s" % (uniprot_id, isoform_id)
                exit(1)
            else:
                print "Correct > %s %s" % (uniprot_id, isoform_id)
                seq = foo[0]
                length = len(seq)
                digest = get_digest_md5(seq)
                db.update_one({"uniprot_id": uniprot_id, "isoform_product.isoform_id": isoform_id}, {"$set": {"isoform_product.$.length": length, "isoform_product.$.seq": seq, "isoform_product.$.digest": digest}})

if __name__ == '__main__':
    main()