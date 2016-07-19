import re


class ParsingException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

# In CC section, they are descriptions of each topic; "FUNCTION" is one of topics.
# RegEx python playground: http://pythex.org/
class Comments:
    def __init__(self, raw_text):
        self.raw_text = raw_text

    def get_raw_text(self):
        return self.raw_text

    # In CC section, we use the variable 'topic' to collect its TEXT STRING
    def __get_content(self, cc, topic):
        topic += ":"
        return [comment.replace(topic + ' ', '').replace('\n', '') for comment in cc if comment.startswith(topic)]

    # Get a topic "FUNCTION" from CC section
    def get_topic_function(self):
        t = self.__get_content(self.get_raw_text(), 'FUNCTION')
        return t

    # Get a topic "ALTERNATIVE PRODUCTS" from CC section
    def get_topic_alternative_products(self):
        d = dict()
        t = self.__get_content(self.get_raw_text(), 'ALTERNATIVE PRODUCTS')
        if len(t) > 1:
            err_mesg = "In CC section, the \'ALTERNATIVE PRODUCTS\' topic was found more than 1"
            raise ParsingException(err_mesg)
        # Put into JSON format
        elif len(t) == 1:
            # Number of alternative splicing isoform from topic "ALTERNATIVE PRODUCTS"
            n_as_isoform = int(re.search('Named isoforms=(\d+);', t[0]).group(1))
            d['n_isoform'] = n_as_isoform
            # Experiment type :
            #     1. Alternative isoform (get AS by natural method?)
            #     2. Alternative initiation (get AS by chemical substance)
            #     etc..
            exp_type = re.search('Event=(.*?);', t[0]).group(1).split(', ')
            d['as_event'] = exp_type
        elif len(t) == 0:
            d['n_isoform'] = 1
        return d

    def get_topic_alternative_product_isoform_product(self):
        l = list()
        t = self.__get_content(self.get_raw_text(), 'ALTERNATIVE PRODUCTS')
        if len(t) > 1:
            err_mesg = "In CC section, the \'ALTERNATIVE PRODUCTS\' topic was found more than 1"
            raise ParsingException(err_mesg)
        # Put into JSON format
        elif len(t) == 1:
            # IsoId=A, B, C ; are all of the protein's alternative splicing isoform from past to present, the first
            # order of the list of IsoId (A) is the current isoform ID
            for matched_string in re.finditer('IsoId=(.*?); +Sequence=(.*?);', t[0]):
                # print matched_string
                history, isoform_note = matched_string.groups()
                isoform_id_history = history.split(',')
                isoform_id_current = isoform_id_history[0]
                d = dict()
                d['isoform_id'] = isoform_id_current
                d['note'] = isoform_note
                if len(isoform_id_history) > 1:
                    d['obsolete_isoform_id'] = []
                for isoform_id in isoform_id_history[1:]:
                    d['obsolete_isoform_id'].append(isoform_id)
                l.append(d)
        elif len(t) == 0:
            d = dict()
            d['isoform_id'] = ''
            l.append(d)
        # print l
        return l

if __name__ == '__main__':
    pass
