import ConfigParser
import sys
import os
from pymongo import MongoClient


class Config:
    __conf = None

    def __init__(self):
        self.__parse()

    def __parse(self):
        self.__conf = ConfigParser.ConfigParser()
        self.__conf.read(os.path.join(os.path.dirname(sys.argv[0]), "config.ini"))\

    def __connect(self):
        uri = "mongodb://%s:%s/%s" % (self.__conf.get('database', 'host'), self.__conf.get('database', 'port'), self.__conf.get('database', 'db'))
        return MongoClient(uri)

    def get_connection_collection(self):
        client = self.__connect()
        print client.get_default_database()
        return client.get_default_database().get_collection(self.__conf.get('database', 'collection'))

    def get_connection_db(self):
        return self.__connect()

    def get_sprot_path(self):
        return self.__conf.get('uniprot', 'sprot_path')

    def get_varsplice_path(self):
        return self.__conf.get('uniprot', 'varsplic_path')

    def get_refseq_protein_path(self):
        return self.__conf.get('refseq', 'protein_seq_path')
