import hashlib

def str_to_md5(in_str):
    md5 = hashlib.md5(in_str.encode("utf-8")).hexdigest()
    return md5
