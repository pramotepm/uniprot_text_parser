import hashlib

def get_digest_md5(seq):
    m = hashlib.md5()
    m.update(seq)
    return m.hexdigest()
