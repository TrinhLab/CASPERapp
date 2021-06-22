import logging

"""File stores all global variables, mainly used for file information/directory monitoring"""

# setup logger in global space
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler('CASPER.log', mode='w')
fh_formatter = logging.Formatter('%(asctime)s %(levelname)s %(lineno)d:%(filename)s(%(process)d) - %(message)s')
fh.setFormatter(fh_formatter)
fh.setLevel(logging.DEBUG)
logger.addHandler(fh)


OPERATING_SYSTEM_ID = ""

CSPR_DB = ""

appdir = ""