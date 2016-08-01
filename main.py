import create_main_schema
import insert_varsplice
import match_refseq_entry
import correct_external

# print "PARSING FROM SWISSPROT DATABASE"
# create_main_schema.main()
# print "\tcomplete"

# print "INSERTING VARIANT REFERENCE FROM SWISSPROT VARSPLICE"
# insert_varsplice.main()
# print "\tcomplete"

# print "CORRECT THE EXTERNAL SEQUENCE"
# correct_external.main()
# print "\tcomplete"

print "DISCOVERY NEW SEQUENCE IN REFSEQ PROTEIN"
match_refseq_entry.main()
print "\tcomplete"
