import re


class ParsingException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class Comments:
    def __init__(self, raw_text):
        self.raw_text = raw_text

    def get_raw_text(self):
        return self.raw_text

    # In CC section, they are descriptions of each topic; "FUNCTION" is one of them.
    # Therefore, after we could analyze with text mining technique, these topics, including "FUNCTION" may be included.
    def __get_content(self, cc, topic):
        topic += ":"
        return [comment.replace(topic + ' ', '').replace('\n', '') for comment in cc if comment.startswith(topic)]

    def get_topic_function(self):
        t = self.__get_content(self.raw_text, 'FUNCTION')
        return t

    def get_topic_alternative_products(self):
        d = dict()
        t = self.__get_content(self.raw_text, 'ALTERNATIVE PRODUCTS')
        if len(t) > 1:
            err_mesg = "In CC section, the \'ALTERNATIVE PRODUCTS\' topic was found more than 1"
            raise ParsingException(err_mesg)
        # Number of alternative splicing isoform
        n_as_isoform = int(re.search('Named isoforms=(\d)', t[0]).group(1)) if len(t) != 0 else 0
        d['nASisoform'] = n_as_isoform
        return d


if __name__ == '__main__':
    pass
