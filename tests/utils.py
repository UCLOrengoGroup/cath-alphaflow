def write_uniprot_ids(fh, uniprot_ids):
    fh.write("uniprot_id\n")
    for uniprot_id in uniprot_ids:
        fh.write(uniprot_id + "\n")
