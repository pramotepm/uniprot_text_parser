import create_main_schema
import insert_varsplice
import match_refseq_entry
import correct_external
import generate_fasta
import fire


def create_main_table():
    print "PARSE SWISSPROT DB"
    create_main_schema.main()
    print "\tcomplete"

    print "INSERTING VARIANT REFERENCE FROM SWISSPROT VARSPLICE"
    insert_varsplice.main()
    print "\tcomplete"

    print "CORRECT THE EXTERNAL SEQUENCE"
    correct_external.main()
    print "\tcomplete"

    print "MAP SWISSPROT TO REFSEQ"
    match_refseq_entry.main()
    print "\tcomplete"


def gen_fasta():
    print "Generating FASTA files..."
    generate_fasta.main()
    print "\tcomplete"


if __name__ == '__main__':
  fire.Fire({
      'create-main-table': create_main_table,
      'generate-fasta': gen_fasta,
  })